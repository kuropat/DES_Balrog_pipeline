def _add_noise_to_obs(obs, noise_image, noise):
    """
    parameters
    ----------
    obs: Observation, ObsList, MultiBandObsList
        The obs
    noise_image: ndarray
        Array of noise to add
        if obs is a MultiBandObsList, noise_image can also be a list
    noise: float
        sigma of noise, for modifying the weight
        if obs is a MultiBandObsList, noise can also be a list
    """
    if isinstance(obs, MultiBandObsList):
        new_mb_obs = MultiBandObsList()
        for i,obslist in enumerate(obs):

            if isinstance(noise_image,list):
                use_noise_image = noise_image[i]
                use_noise = noise[i]
            else:
                use_noise_image = noise_image
                use_noise = noise

            new_obslist=_add_noise_to_obs(obslist,use_noise_image,use_noise)
            new_mb_obs.append(new_obslist)

        return new_mb_obs

    elif isinstance(obs, ObsList):
        new_obslist = ObsList()
        for i, tobs in enumerate(obs):
            new_obs = _add_noise_to_obs(tobs, noise_image, noise)
            new_obslist.append(new_obs)

        return new_obslist

    elif isinstance(obs, Observation):

        new_weight = obs.weight.copy()
        w=numpy.where(new_weight > 0)
        new_weight[w] = 1.0/(1.0/new_weight[w] + noise**2)

        new_obs = Observation(obs.image + noise_image,
                              weight=new_weight,
                              jacobian=obs.jacobian.copy(),
                              psf=deepcopy(obs.psf) )

        return new_obs

    else:
        raise ValueError("obs should be an Observation,ObsList,MultiBandObsList")

class MetacalDetrendNGMixFitter(MetacalNGMixBootFitter):
    def _fit_galaxy(self, model, coadd, guess=None,**kw):

        # this runs the psfs, max fitter, and metacal
        super(MetacalDetrendNGMixFitter,self)._fit_galaxy(
            model,
            coadd,
            guess=guess,
            **kw
        )

        boot=self.boot
        obs_dict_orig=boot.get_metacal_max_result()['obs_dict']

        mb_obs_list=boot.mb_obs_list
        if len(mb_obs_list) > 1 or len(mb_obs_list[0]) > 1:
            raise NotImplementedError("fix to work with multiple obs/bands")

        obs = mb_obs_list[0][0]
        im=obs.image
        wt=obs.weight

        noise_image1 = numpy.random.normal(loc=0.0,
                                           scale=1.0,
                                           size=im.shape)
        w=numpy.where(wt > 0)
        base_noise = numpy.sqrt( numpy.median(1.0/wt[w]) )

        #
        # now add extra noise before and after metacal
        #

        Rnoise_types=['1p','1m','2p','2m']
        #Rnoise_types=None
        new_results=[]

        detrend_noises = self['target_noise']*numpy.array(self['detrend_factors'])
        print("    doing detrend noise")

        for i, dtnoise in enumerate(detrend_noises):
            extra_noise = numpy.sqrt(dtnoise**2 - base_noise**2)

            #print("    doing detrend noise: %.3f "
            #      "extra noise: %.3f" % (dtnoise,extra_noise))


            # same noise image, just scaled
            noise_image = noise_image1*extra_noise

            #
            # adding noise *before* metacal
            # new psf observations are generated, psf models are refit currently
            #
            mb_obs_before = _add_noise_to_obs(mb_obs_list, noise_image, extra_noise)
            mcal_obs_before = ngmix.metacal.get_all_metacal(
                mb_obs_before,
                types=Rnoise_types,
                **self['metacal_pars']
            )
            self._do_metacal(model, boot, metacal_obs=mcal_obs_before)
            res_before = boot.get_metacal_max_result()

            #
            # adding noise *after* metacal
            # psf models get copied over
            #
            mcal_obs_after = {}
            for key,tobs in obs_dict_orig.iteritems():
                if key in Rnoise_types:
                    new_obs = _add_noise_to_obs(tobs, noise_image, extra_noise)
                    mcal_obs_after[key] = new_obs

            self._do_metacal(model, boot, metacal_obs=mcal_obs_after)
            res_after = boot.get_metacal_max_result()

            Rnoise = res_before['mcal_R']    - res_after['mcal_R']
            #Rnoise_psf = res_before['mcal_Rpsf']    - res_after['mcal_Rpsf']

            new_res={
                'mcal_Rnoise':Rnoise,
                #'mcal_Rnoise_psf':Rnoise_psf,
            }
            new_results.append(new_res)

        res=self.gal_fitter.get_result()
        res['mcal_dt_results'] = new_results

    def _copy_galaxy_result(self, model, coadd):
        """
        copy parameters specific to this class
        """
        super(MetacalDetrendNGMixFitter,self)._copy_galaxy_result(model, coadd)

        dindex=0
        res=self.gal_fitter.get_result()

        if res['flags'] == 0:
            tmodel = '%s_mcal' % model
            n = self._get_namer(tmodel, coadd)

            d=self.data

            for idt,dtres in enumerate(res['mcal_dt_results']):

                f=n('dt_Rnoise')
                d[f][dindex,idt,:,:] = dtres['mcal_Rnoise']

                #f=n('dt_Rnoise_psf')
                #d[f][dindex,idt,:,:] = dtres['mcal_Rnoise_psf']


    def _get_metacal_dtype(self, coadd):
        dt=super(MetacalDetrendNGMixFitter,self)._get_metacal_dtype(coadd)

        ndetrend = len(self['detrend_factors'])

        for model in self._get_all_models(coadd):
            front='%s_mcal' % (model)
            n=Namer(front)
            dt += [
                (n('dt_Rnoise'),'f8',(ndetrend,2,2)),
                #(n('dt_Rnoise_psf'),'f8',(ndetrend,2)),
            ]

        # not adding these to mcal_flist
        #self.mcal_flist = [d[0] for d in dt]

        return dt

class PostcalNGMixBootFitter(MetacalNGMixBootFitter):
    def _fit_galaxy(self, model, coadd, guess=None,**kwargs):
        mb_obs_list = self.boot.mb_obs_list

        super(MetacalNGMixBootFitter,self)._fit_galaxy(model,
                                                       coadd,
                                                       guess=guess,
                                                       **kwargs)
        postcal_res=self._do_postcal(model)

        res=self.gal_fitter.get_result()
        res.update(postcal_res)

    def _get_postcal_obsdict(self, mb_obs_list):
        import galsim

        obs = mb_obs_list[0][0]
        psf_obs = obs.psf

        im = obs.image.copy()
        psf_im = psf_obs.image.copy()
        gs_im = galsim.Image(im, scale=1.0)
        gs_psf = galsim.Image(psf_im, scale=1.0)

        i_im = galsim.InterpolatedImage(gs_im)
        i_psf = galsim.InterpolatedImage(gs_psf)

        step=self['postcal_pars']['step']

        shts = [('1p',Shape( step, 0.0)),
                ('1m',Shape(-step, 0.0)),
                ('2p',Shape( 0.0,  step)),
                ('2m',Shape( 0.0, -step))]

        odict={}
        for t in shts:
            name = t[0]
            shear = t[1]

            s_i_im = i_im.shear(g1=shear.g1, g2=shear.g2)
            s_i_psf = i_psf.shear(g1=shear.g1, g2=shear.g2)

            s_im = s_i_im.drawImage(ny=im.shape[0],
                                    nx=im.shape[1],
                                    scale=1.0,
                                    method='no_pixel')
            s_psf_im = s_i_psf.drawImage(ny=psf_im.shape[0],
                                         nx=psf_im.shape[1],
                                         scale=1.0,
                                         method='no_pixel')

            spsf_obs = Observation(
                s_psf_im.array,
                weight=psf_obs.weight.copy(),
                jacobian=psf_obs.jacobian.copy()
            )
            sobs = Observation(
                s_im.array,
                weight=obs.weight.copy(),
                jacobian=obs.jacobian.copy(),
                psf=spsf_obs
            )

            odict[name] = sobs

        return odict

    def _do_postcal(self,
                    model,
                    boot=None):
        """
        the basic fitter for this class
        """

        print("    doing postcal")

        if boot is None:
            boot=self.boot

        odict = self._get_postcal_obsdict(boot.mb_obs_list)

        Tguess = boot.mb_obs_list[0][0].psf.meta['Tguess']

        psf_pars = {}
        for k,v in self['psf_pars'].iteritems():
            if k != 'model' and k != 'ntry':
                psf_pars.update({k:v})

        fits={}
        for key in odict:
            obs = odict[key]

            tboot=self._get_bootstrapper(model, obs)
            
            tboot.fit_psfs(self['psf_pars']['model'],
                           Tguess,
                           ntry=self['psf_pars']['ntry'],
                           fit_pars=psf_pars)


            self._fit_max(model, boot=tboot)

            tboot.set_round_s2n()
            res=tboot.get_max_fitter().get_result()
            rres=tboot.get_round_result()
            res['s2n_r'] = rres['s2n_r']
            res['T_r'] = rres['T_r']

            fits[key] = res

        res = self._extract_postcal_responses(fits)
        return res

    def _extract_postcal_responses(self, fits):
        """
        pars pars_cov gpsf, s2n_r, T_r, psf_T_r required

        expect the shape to be in pars[2] and pars[3]
        """
        step = self['postcal_pars']['step']

        res1p = fits['1p']
        res1m = fits['1m']
        res2p = fits['2p']
        res2m = fits['2m']

        pars_mean = (res1p['pars']+
                     res1m['pars']+
                     res2p['pars']+
                     res2m['pars'])/4.0

        pars_cov_mean = (res1p['pars_cov']+
                         res1m['pars_cov']+
                         res2p['pars_cov']+
                         res2m['pars_cov'])/4.0

        pars_mean[2] = 0.5*(fits['1p']['pars'][2] + fits['1m']['pars'][2])
        pars_mean[3] = 0.5*(fits['2p']['pars'][3] + fits['2m']['pars'][3])

        s2n_r_mean = (res1p['s2n_r']
                      + res1m['s2n_r']
                      + res2p['s2n_r']
                      + res2m['s2n_r'])/4.0

        if self['verbose']:
            print_pars(pars_mean, front='    parsmean:   ')

        R=numpy.zeros( (2,2) ) 
        Rpsf=numpy.zeros(2)

        fac = 1.0/(2.0*step)

        R[0,0] = (fits['1p']['pars'][2]-fits['1m']['pars'][2])*fac
        R[0,1] = (fits['1p']['pars'][3]-fits['1m']['pars'][3])*fac
        R[1,0] = (fits['2p']['pars'][2]-fits['2m']['pars'][2])*fac
        R[1,1] = (fits['2p']['pars'][3]-fits['2m']['pars'][3])*fac

        #Rpsf[0] = (pars['1p_psf'][2]-pars['1m_psf'][2])*fac
        #Rpsf[1] = (pars['2p_psf'][3]-pars['2m_psf'][3])*fac


        #gpsf_name = 'pcal_%spsf' % shape_type
        #raw_gpsf_name = '%spsf' % shape_type
        res = {
            'pcal_pars':pars_mean,
            'pcal_pars_cov':pars_cov_mean,
            'pcal_g':pars_mean[2:2+2],
            'pcal_g_cov':pars_cov_mean[2:2+2, 2:2+2],
            'pcal_R':R,
            #'pcal_Rpsf':Rpsf,
            #'pcal_gpsf':fits['gpsf'],
            'pcal_s2n_r':s2n_r_mean,
        }
        return res




    def _get_fit_data_dtype(self,coadd):
        dt=super(MetacalNGMixBootFitter,self)._get_fit_data_dtype(coadd)

        dt_pcal = self._get_postcal_dtype(coadd)
        dt += dt_pcal
        return dt

    def _copy_galaxy_result(self, model, coadd):
        super(MetacalNGMixBootFitter,self)._copy_galaxy_result(model,coadd)

        dindex=0
        res=self.gal_fitter.get_result()

        tmodel='%s_pcal' % model
        n = self._get_namer(tmodel, coadd)

        if res['flags'] == 0:

            for f in self.pcal_flist:
                front='%s_' % (model)
                mf = f.replace(front,'')
                #print("copying %s -> %s" % (f,mf), res[mf])
                self.data[f][dindex] = res[mf]

    def _get_postcal_dtype(self, coadd):

        nband=self['nband']
        np=5+nband

        dt=[]
        for model in self._get_all_models(coadd):
            front='%s_pcal' % (model)
            n=Namer(front)
            dt += [
                (n('pars'),'f8',np),
                (n('pars_cov'),'f8',(np,np)),
                (n('g'),'f8',2),
                (n('g_cov'),'f8', (2,2) ),
                #(n('c'),'f8',2),
                (n('s2n_r'),'f8'),
                (n('R'),'f8',(2,2)),
                #(n('Rpsf'),'f8',2),
                #(n('gpsf'),'f8',2),
            ]

        self.pcal_flist = [d[0] for d in dt]

        return dt


    def _make_struct(self,coadd):
        data=super(MetacalNGMixBootFitter,self)._make_struct(coadd)

        models=self._get_all_models(coadd)
        for model in models:
            for f in self.pcal_flist:
                data[f] = DEFVAL

        return data


class PostcalNGMixSimFitter(PostcalNGMixBootFitter):
    def _do_postcal(self, model):

        mb_obs_list = boot.mb_obs_list
        if len(mb_obs_list) > 1 or len(mb_obs_list[0]) > 1:
            raise NotImplementedError("only a single obs for now")


        res=super(PostcalNGMixSimFitter,self)._do_postcal(model)


        print("    Calculating Rnoise")

        fitter = self.boot.get_max_fitter()
        gmix_list = []
        for band in xrange(self['nband']):
            gm = fitter.get_gmix(band=band)
            gmix_list.append(gm)


        # for noise added *before* metacal steps
        mobs_before = ngmix.simobs.simulate_obs(
            gmix_list,
            mb_obs_list,
            add_noise=True,
            convolve_psf=True
        )
        # for noise added *after* metacal steps
        mobs_after = ngmix.simobs.simulate_obs(
            gmix_list,
            mb_obs_list,
            add_noise=False,
            convolve_psf=True
        )


        boot_model_before = self._get_bootstrapper(model,mobs_before)
        boot_model_after = self._get_bootstrapper(model,mobs_after)

        mcal_obs_after = ngmix.metacal.get_all_metacal(
            mobs_after[0][0],
            self['metacal_pars']['step'],
        )

        # now add noise after creating the metacal observations
        # using the same noise image!

        noise = mobs_before[0][0].noise_image
        for key in mcal_obs_after:
            obs=mcal_obs_after[key]
            obs.image = obs.image + noise

        res_before=self._do_postcal(
            model,
            boot=boot_model_before
        )
        res_after=self._do_postcal(
            model,
            boot=boot_model_after,
            metacal_obs=mcal_obs_after
        )

        res_before = boot_model_before.get_metacal_max_result()
        res_after = boot_model_after.get_metacal_max_result()

        gnoise = res_before['mcal_g'] - res_after['mcal_g']
        Rnoise = res_before['mcal_R'] - res_after['mcal_R']
        Rpsf_noise = res_before['mcal_Rpsf'] - res_after['mcal_Rpsf']

        res = self.boot.get_max_fitter().get_result()
        res['mcal_gnoise'] = gnoise
        res['mcal_Rnoise'] = Rnoise
        res['mcal_Rpsf_noise'] = Rpsf_noise





    def _get_postcal_dtype(self, coadd):
        dt=super(PostcalNGMixSimFitter,self)._get_postcal_dtype(coadd)

        nband=self['nband']
        np=5+nband

        for model in self._get_all_models(coadd):
            front='%s_pcal' % (model)
            n=Namer(front)
            dt += [
                (n('Rnoise'),'f8',(2,2)),
            ]

        self.pcal_flist += [d[0] for d in dt]

        return dt




class MetacalSubnNGMixBootFitter(MetacalNGMixBootFitter):
    def _get_subn_metacal_obs(self):
        """
        subtract a correlated noise image sheared by
        negative of the shear applied to the real obs
        """
        from ngmix.metacal import get_all_metacal

        print("    subtracting sheared correlated noise")
        mb_obs_list = self.boot.mb_obs_list

        # simulated noise for these observations

        step = self['metacal_pars']['step']
        # currentlly only works for single obs
        mcal_obs = get_all_metacal(mb_obs_list, step)

        nrand = self.get('subn_nrand',1)

        for irand in xrange(nrand):
            tnoise_mb_obs = ngmix.simobs.simulate_obs(None,
                                                      mb_obs_list)
            tmcal_noise_obs = get_all_metacal(tnoise_mb_obs, step)

            if irand==0:
                mcal_noise_obs = tmcal_noise_obs 
            else:
                #print("    adding extra realization:",irand)
                for key in mcal_noise_obs:

                    # these are MultiBandObsLists
                    mobs = mcal_obs[key]
                    nmobs = mcal_noise_obs[key]
                    tnmobs = tmcal_noise_obs[key] 

                    assert len(mobs)==len(nmobs)

                    # loop over bands
                    for band in xrange(len(mobs)):
                        # mobs[band] etc are ObsLists

                        # append copies of the original, for adding
                        # the noise images to later
                        nn = len(tnmobs[band])

                        mobs[band].extend( deepcopy(mobs[band][0:nn] ) )

                        # append more observations from the ObsList
                        nmobs[band].extend( tnmobs[band] )

                        assert len(mobs[band])==len(nmobs[band])


        # add the noise sheared by negative of the shear
        # applied to observation

        for ipairs in [('1p','1m'),
                       ('1m','1p'),
                       ('2p','2m'),
                       ('2m','2p'),
                       ('1p_psf','1m_psf'),
                       ('1m_psf','1p_psf'),
                       ('2p_psf','2m_psf'),
                       ('2m_psf','2p_psf')]:


            mk=ipairs[0]
            nk=ipairs[1]

            imbobs = mcal_obs[mk]
            nmbobs = mcal_noise_obs[nk]

            for imb in xrange(len(imbobs)):
                iolist=imbobs[imb]
                nolist=nmbobs[imb]

                for iobs in xrange(len(iolist)):

                    obs  = iolist[iobs]
                    nobs = nolist[iobs]

                    im  = obs.image
                    nim = nobs.image

                    obs.image = im + nim
                    obs.weight = 0.5*obs.weight

        return mcal_obs

    def _fit_galaxy(self, model, coadd, guess=None,**kwargs):
        mb_obs_list = self.boot.mb_obs_list

        super(MetacalNGMixBootFitter,self)._fit_galaxy(model,
                                                       coadd,
                                                       guess=guess,
                                                       **kwargs)
        metacal_obs = self._get_subn_metacal_obs()
        self._do_metacal(model, self.boot, metacal_obs=metacal_obs)

        metacal_res = self.boot.get_metacal_max_result()

        res=self.gal_fitter.get_result()
        res.update(metacal_res)


class MetacalSimnNGMixBootFitter(MetacalNGMixBootFitter):

    def _fit_galaxy(self, model, coadd, guess=None,**kwargs):
        mb_obs_list = self.boot.mb_obs_list
        if len(mb_obs_list) > 1 or len(mb_obs_list[0]) > 1:
            raise NotImplementedError("only a single obs for now")


        super(MetacalSimnNGMixBootFitter,self)._fit_galaxy(
            model,
            coadd,
            guess=guess,
            **kwargs
        )

        print("    Calculating Rnoise")

        fitter = self.boot.get_max_fitter()
        gmix_list = []
        for band in xrange(self['nband']):
            gm = fitter.get_gmix(band=band)
            gmix_list.append(gm)

        # for noise added *before* metacal steps
        mobs_before = ngmix.simobs.simulate_obs(
            gmix_list,
            mb_obs_list,
            add_noise=True,
            convolve_psf=True
        )
        # for noise added *after* metacal steps
        mobs_after = ngmix.simobs.simulate_obs(
            gmix_list,
            mb_obs_list,
            add_noise=False,
            convolve_psf=True
        )


        boot_model_before = self._get_bootstrapper(model,mobs_before)
        boot_model_after = self._get_bootstrapper(model,mobs_after)

        mcal_obs_after = ngmix.metacal.get_all_metacal(
            mobs_after[0][0],
            self['metacal_pars']['step'],
        )

        # now add noise after creating the metacal observations
        # using the same noise image!

        noise = mobs_before[0][0].noise_image
        for key in mcal_obs_after:
            obs=mcal_obs_after[key]
            obs.image = obs.image + noise

        self._do_metacal(
            model,
            boot_model_before,
        )
        self._do_metacal(
            model,
            boot_model_after,
            metacal_obs=mcal_obs_after
        )

        res_before = boot_model_before.get_metacal_max_result()
        res_after = boot_model_after.get_metacal_max_result()

        gnoise = res_before['mcal_g'] - res_after['mcal_g']
        Rnoise = res_before['mcal_R'] - res_after['mcal_R']
        Rpsf_noise = res_before['mcal_Rpsf'] - res_after['mcal_Rpsf']

        res = self.boot.get_max_fitter().get_result()
        res['mcal_gnoise'] = gnoise
        res['mcal_Rnoise'] = Rnoise
        res['mcal_Rpsf_noise'] = Rpsf_noise

    def _get_metacal_dtype(self, coadd):
        dt = super(MetacalSimnNGMixBootFitter,self)._get_metacal_dtype(coadd)

        dt_extra=[]
        for model in self._get_all_models(coadd):
            tmodel='%s_pcal' % model
            n = self._get_namer(tmodel, coadd)
            dt_extra += [
                (n('gnoise'), 'f8', 2),
                (n('Rnoise'), 'f8', (2,2)),
                (n('Rpsf_noise'), 'f8', 2),
            ]

        self.mcal_flist += [d[0] for d in dt_extra]

        dt += dt_extra

        return dt

class MetacalAddnNGMixBootFitter(MetacalSimnNGMixBootFitter):
    def _get_noisier_mobs(self):

        fac=self['simnoise']['noise_fac']
        ofacsq=1.0/(1.0 + fac**2)
        #ofacsq=1.0

        mobs=self.mb_obs_list

        new_mobs=MultiBandObsList()
        for olist in mobs:

            new_olist=ObsList()
            for obs in olist:
                noise_image = ngmix.simobs.get_noise_image(obs.weight)

                noise_image *= fac
                new_weight = obs.weight.copy()
                new_weight *= ofacsq

                new_image = obs.image.copy()
                new_image += noise_image

                new_obs = Observation(
                    new_image,
                    weight=new_weight,
                    jacobian=obs.jacobian.copy(),
                    psf=obs.psf
                )

                new_obs.noise_image=noise_image

                new_olist.append(new_obs)
            new_mobs.append(new_olist)

        return new_mobs

    def _fit_galaxy(self, model, coadd, guess=None,**kwargs):
        # note jumping over parent here
        super(MetacalSimnNGMixBootFitter,self)._fit_galaxy(model,
                                                           coadd,
                                                           guess=guess,
                                                           **kwargs)

        if len(self.mb_obs_list) > 1 or len(self.mb_obs_list[0]) > 1:
            raise NotImplementedError("only a single obs for now")

        print("    Calculating Rnoise by adding noise to image")

        mobs_before = self._get_noisier_mobs()

        boot_model_before = self._get_bootstrapper(model,mobs_before)
        boot_model_after = self._get_bootstrapper(model,self.mb_obs_list)

        mcal_obs_after = ngmix.metacal.get_all_metacal(
            self.mb_obs_list[0][0],
            self['metacal_pars']['step'],
        )

        # now add noise after creating the metacal observations
        # using the same noise image!

        noise_image = mobs_before[0][0].noise_image
        new_weight = mobs_before[0][0].weight
        for key in mcal_obs_after:
            obs=mcal_obs_after[key]
            obs.image = obs.image + noise_image
            obs.weight = new_weight.copy()

        self._do_metacal(model, boot_model_before)
        self._do_metacal(model, boot_model_after,
                         metacal_obs=mcal_obs_after)

        res_before = boot_model_before.get_metacal_max_result()
        res_after = boot_model_after.get_metacal_max_result()

        fac2inv=1.0/self['simnoise']['noise_fac']**2
        #print("        multiplying by:",fac2inv)
        Rnoise = (res_before['mcal_R'] - res_after['mcal_R'])*fac2inv
        Rpsf_noise = (res_before['mcal_Rpsf'] - res_after['mcal_Rpsf'])*fac2inv

        print("    Rnoise[0,0]: %g" % Rnoise[0,0])

        res = self.boot.get_max_fitter().get_result()
        res['mcal_Rnoise'] = Rnoise
        res['mcal_Rpsf_noise'] = Rpsf_noise


class MetacalRegaussBootFitter(MaxNGMixBootFitter):

    def _guess_and_run_boot(self,model,mb_obs_list,coadd,**kw):

        flags=0

        boot=get_bootstrapper(mb_obs_list, find_cen=False, **self)
        self.boot=boot

        res={}
        ppars=self['psf_pars']
        rpars=self['regauss_pars']

        try:
            boot.fit_metacal_regauss(
                ppars['Tguess'],
                rpars['Tguess'],
                psf_ntry=ppars['ntry'],
                ntry=rpars['ntry'],
                metacal_pars=self['metacal_pars'],
            )
            mres=boot.get_metacal_regauss_result()
            res.update(mres)

            # yikes, monkey patching
            # TODO figure out how to fix this
            boot._result = res

            self._copy_galaxy_result()
            self._print_galaxy_result()

        except BootPSFFailure as err:
            print("    psf fitting failed")
            flags = PSF_FIT_FAILURE
        except (BootGalFailure,GMixRangeError) as err:
            print("    galaxy fitting failed")
            flags = GAL_FIT_FAILURE

        except BootPSFFailure:
            print("    psf fitting failed")
            flags = PSF_FIT_FAILURE

        res['flags'] = flags
        return flags, boot

    def _print_galaxy_result(self):
        print("implement print")

    def _copy_galaxy_result(self):
        dindex=0
        res=self.boot.get_metacal_regauss_result()

        data=self.data

        if res['flags'] == 0:
            for f in ['pars','pars_cov','e','e_cov',
                      'c','R', 'Rpsf','epsf']:

                mf = 'mcal_%s' % f
                data[f][dindex] = res[mf]

    def _get_fit_data_dtype(self, coadd):

        np=6

        nband=self['nband']
        bshape=(nband,)

        dt = [
            ('nimage_use','i4',bshape),
            ('pars','f8',np),
            ('pars_cov','f8',(np,np)),
            ('e','f8',2),
            ('e_cov','f8', (2,2) ),
            ('c','f8',2),
            ('R','f8',(2,2)),
            ('Rpsf','f8',2),
            ('epsf','f8',2),
        ]

        return dt

    def _make_struct(self,coadd):

        dt=self._get_fit_data_dtype(coadd)
        num=1
        data=numpy.zeros(num, dtype=dt)

        for n in data.dtype.names:
            if n != 'flags':

                data[n] = DEFVAL

        return data
