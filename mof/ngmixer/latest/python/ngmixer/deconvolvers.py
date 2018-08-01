from __future__ import print_function
import time
import numpy
from numpy import array
import os
import scipy.stats

# local imports
from .defaults import DEFVAL, PDEFVAL, NO_ATTEMPT, \
    PSF_FIT_FAILURE, GAL_FIT_FAILURE, \
    LOW_PSF_FLUX, PSF_FLUX_FIT_FAILURE, \
    NBR_HAS_NO_PSF_FIT, METACAL_FAILURE
from .fitting import BaseFitter
from .util import Namer, print_pars
from .render_ngmix_nbrs import RenderNGmixNbrs

# ngmix imports
import ngmix
from ngmix import Observation, ObsList, MultiBandObsList, GMixRangeError
from ngmix.gexceptions import BootPSFFailure, BootGalFailure
from ngmix.gmix import GMix
from ngmix.jacobian import Jacobian

from pprint import pprint

from .bootfit import NGMixBootFitter


class Deconvolver(NGMixBootFitter):
    def __init__(self,*args,**kw):
        super(Deconvolver,self).__init__(*args,**kw)
        self['normalize_psf'] = self.get('normalize_psf',True)
        self['deconv_pars']['trim_images'] = self['deconv_pars'].get('trim_images',False)

    def _setup(self):
        """
        ngmix-specific setups for all ngmix fitters
        """

        # whether to use the coadd_ prefix for names when no me fit was done
        self['use_coadd_prefix'] = self.get('use_coadd_prefix',True)

    def __call__(self,mb_obs_list,coadd=False, **kw):
        """
        fit the obs list
        """

        if 'target_noise' in self:
            self._add_extra_sim_noise(mb_obs_list)

        # only fit stuff that is not flagged
        new_mb_obs_list = self._get_good_mb_obs_list(mb_obs_list)
        self.new_mb_obs_list = new_mb_obs_list
        self.mb_obs_list = mb_obs_list

        mb_obs_list.update_meta_data({'fit_data':self._make_struct(coadd)})
        self.data = mb_obs_list.meta['fit_data']

        if self['make_plots']:
            self.plot_dir = './%d-plots' % new_mb_obs_list.meta['id']
            if not os.path.exists(self.plot_dir):
                os.makedirs(self.plot_dir)

        '''
        cres=self._find_center(mb_obs_list)
        if cres['flags'] != 0:
            print("could not find center")
            return GAL_FIT_FAILURE
        '''

        dpars=self['deconv_pars']

        flags=0
        res=None
        try:
            # this sets self.boot, which is only used for psf stuff
            self._fit_psfs(new_mb_obs_list)
            flags |= self._fit_psf_flux(coadd)

            try:
                # to get a center
                if dpars['trim_images']:
                    self._fit_gauss()

                res = self._do_deconv(new_mb_obs_list)
                self._copy_galaxy_result(res, coadd)

                if res['flags'] != 0:
                    print("    deconv failed with flags:",res['flags'])
                    flags |= GAL_FIT_FAILURE

            except BootGalFailure as err:
                print("    galaxy fitting failed: %s" % str(err))
                flags = GAL_FIT_FAILURE

        except BootPSFFailure as err:
            print("    psf fitting failed")
            flags = PSF_FIT_FAILURE

        # fill the epoch data
        self._fill_epoch_data(mb_obs_list,new_mb_obs_list)

        # fill in PSF stats in data rows
        if (flags & PSF_FIT_FAILURE) == 0:
            self._do_psf_stats(mb_obs_list,coadd)

        if res is not None:
            self._fill_nimage_used(res, coadd)

        return flags

    def _do_deconv(self, mb_obs_list):
        import deconv

        mfitter=self.boot.get_max_fitter()
        mres=mfitter.get_result()

        cen=mres['pars'][0:0+2].copy()

        if self['deconv_pars']['trim_images']:
            mb_obs_list=self._trim_images(mb_obs_list, cen)

        meas=deconv.measure.calcmom_ksigma_obs(
            mb_obs_list,
            self['sigma_weight'],  # arcsec
            dk=self['dk'],         # 1/arcsec or None
        )

        return meas

    def _trim_images(self, mbo_input, censky):
        """
        cen in sky coords, relative to jacobian
        center
        """

        mbo = MultiBandObsList()
        for obslist in mbo_input:

            new_obslist=ObsList()
            for obs in obslist:

                j=obs.jacobian
                scale=j.get_scale()
                cenpix=array( j.get_cen() )

                new_cen = cenpix + censky/scale
                #print("cen0:",cenpix,"newcen:",new_cen)

                new_im=_trim_image(obs.image, new_cen)
                new_wt=_trim_image(obs.weight, new_cen)
                new_j = j.copy()
                new_j.set_cen(row=new_cen[0], col=new_cen[1])

                newobs = Observation(
                    new_im,
                    weight=new_wt,
                    jacobian=new_j,
                    psf=obs.psf,
                )

                new_obslist.append( newobs )
            mbo.append( new_obslist )

        return mbo

    '''
    def _find_center(self, mb_obs_list):
        scale=mb_obs_list[0][0].jacobian.get_scale()
        Tguess=4.0 * scale**2

        max_pars=self['max_pars']
        fit_pars=max_pars['fit_pars']

        guesser=ngmix.guessers.TFluxGuesser(
            Tguess,
            Fguess,
            scaling='linear'
        )
        runner=ngmix.bootstrap.MaxRunner(
            mb_obs_list,
            'gauss',
            Tguess,
            fit_pars,
        )
        runner.go(ntry=max_pars['ntry'])
        res=runner.fitter.get_result()
        return res
    '''


    def _scale_image(self, imin):
        im=imin.copy()

        maxval=im.max()
        #im = numpy.log10(im.clip(min=1.0e-4*maxval))
        im = numpy.log10(im.clip(min=1.0e-3*maxval))
        im -= im.min()
        im /= im.max()

        return im


    def _fill_nimage_used(self,res,coadd):
        nim_used = numpy.zeros(self['nband'],dtype='i4')
        for band in xrange(self['nband']):
            nim_used[band] = res['nimage_use_band'][band]

        n=self._get_namer('', coadd)

        self.data[n('nimage_use')][0,:] = nim_used

    def _fit_psfs(self, mb_obs_list):
        """
        fit the psf model to every observation's psf image
        """

        print('    fitting the PSFs')
        boot=ngmix.bootstrap.Bootstrapper(mb_obs_list)

        psf_pars=self['psf_pars']
        fit_pars=psf_pars['fit_pars']

        boot.fit_psfs(psf_pars['model'],
                      None,
                      Tguess_key='Tguess',
                      ntry=psf_pars['ntry'],
                      fit_pars=fit_pars,
                      norm_key='psf_norm')

        # check for no obs in a band if PSF fit fails
        for band,obs_list in enumerate(boot.mb_obs_list):
            if len(obs_list) == 0:
                raise BootPSFFailure("psf fitting failed - band %d has no obs" % band)
            #for obs in obs_list:
            #    sigma=numpy.sqrt( obs.psf.gmix.get_T()/2.0 )
            #    sigmapix = sigma/obs.psf.jacobian.get_scale()
            #    print("sigma pixels:",sigmapix)

        self.boot=boot

    def _fit_gauss(self):
        """
        to get a center
        """

        prior=self['model_pars']['gauss']['prior']

        max_pars=self['max_pars']

        try:
            self.boot.fit_max(
                'gauss',
                max_pars,
                ntry=max_pars['ntry'],
                prior=prior,
            )
        except GMixRangeError:
            raise BootGalFailure("failure fitting gauss")

        res=self.boot.get_max_fitter().get_result()
        if res['flags'] != 0:
            raise BootGalFailure("failure fitting gauss")


    def _copy_galaxy_result(self, meas, coadd):
        """
        Copy from the result dict to the output array
        """

        dindex=0
        data=self.data

        res=meas.get_result()

        data['dc_flags'][dindex] = res['flags']
        data['dc_orflags'][dindex] = res['orflags']
        data['dc_flags_band'][dindex] = res['flags_band']
        data['dc_orflags_band'][dindex] = res['orflags_band']

        data['dk'][dindex] = res['dk']

        data['T'][dindex] = res['T']
        data['e'][dindex] = res['e']
        data['wflux'][dindex] = res['wflux']
        data['wflux_band'][dindex] = res['wflux']

    def _print_galaxy_result(self):
        res=self.gal_fitter.get_result()
        if 'pars' in res:
            print_pars(res['pars'],    front='    gal_pars: ')
        if 'pars_err' in res:
            print_pars(res['pars_err'],front='    gal_perr: ')

        mres=self.boot.get_max_fitter().get_result()
        if 's2n_w' in mres:
            rres=self.boot.get_round_result()            
            tup=(mres['s2n_w'],
                 rres['s2n_r'],
                 mres['chi2per'],
                 mres['chi2per']*mres['dof'],
                 mres['dof'],
                 scipy.stats.chi2.sf(mres['chi2per']*mres['dof'],mres['dof']))
            print("    s2n: %.1f s2n_r: %.1f chi2per: %.3f (chi2: %.3g dof: %.3g pval: %0.3g)" % tup)

    def get_num_pars_psf(self):
        npdict = {'em1':6,
                  'em2':6*2,
                  'em3':6*3,
                  'coellip2':6*2, # this is the "full pars" count, not that used in fit
                  'coellip3':6*3,
                  'turb':6*3,
                  'gauss':6}
        model = self['psf_pars']['model'].lower()
        assert model in npdict,"psf model %s not allowed in NGMixBootFitter" % model
        return npdict[model]

    def get_fit_data_dtype(self,me,coadd):
        dt = []
        if me:
            dt += self._get_fit_data_dtype(False)
        if coadd:
            dt += self._get_fit_data_dtype(True)
        return dt

    def _get_fit_data_dtype(self,coadd):
        dt=[]

        nband=self['nband']
        bshape=(nband,)
        simple_npars=5+nband

        n=self._get_namer('psf', coadd)

        dt += [(n('flags'),   'i4',bshape),
               (n('flux'),    'f8',bshape),
               (n('flux_err'),'f8',bshape),
               (n('flux_s2n'),'f8',bshape)]

        n=self._get_namer('', coadd)

        dt += [(n('nimage_use'),'i4',bshape)]

        dt += [(n('mask_frac'),'f8'),
               (n('psfrec_T'),'f8'),
               (n('psfrec_g'),'f8', 2)]

        if nband==1:
            fcov_shape=(nband,)
        else:
            fcov_shape=(nband,nband)

        dt += [
            ('dc_flags','i4'),
            ('dc_orflags','i4'),
            ('dc_flags_band','i4',bshape),
            ('dc_orflags_band','i4',bshape),

            ('dk','f8'),

            ('T','f8'),
            ('e','f8',2),
            ('wflux','f8',bshape),
            ('wflux_band','f8',bshape),
        ]

        return dt

    def _make_struct(self,coadd):
        """
        make the output structure
        """
        num = 1
        dt = self._get_fit_data_dtype(coadd)
        data = numpy.zeros(num,dtype=dt)

        n=self._get_namer('psf', coadd)

        data[n('flags')] = NO_ATTEMPT
        data[n('flux')] = DEFVAL
        data[n('flux_err')] = DEFVAL
        data[n('flux_s2n')] = DEFVAL

        n=self._get_namer('', coadd)

        data[n('mask_frac')] = DEFVAL
        data[n('psfrec_T')] = DEFVAL
        data[n('psfrec_g')] = DEFVAL

        # the deconvolve parameters
        data['dc_flags']        = NO_ATTEMPT
        data['dc_orflags']      = NO_ATTEMPT
        data['dc_flags_band']   = NO_ATTEMPT
        data['dc_orflags_band'] = NO_ATTEMPT
        data['dk']              = DEFVAL
        data['T']               = DEFVAL
        data['e']               = DEFVAL
        data['wflux']           = DEFVAL
        data['wflux_band']      = DEFVAL

        return data

    def get_default_fit_data(self,me,coadd):
        dt = self.get_fit_data_dtype(me,coadd)
        d = numpy.zeros(1,dtype=dt)
        if me:
            dme = self._make_struct(False)
            for tag in dme.dtype.names:
                d[tag] = dme[tag]

        if coadd:
            dcoadd = self._make_struct(True)
            for tag in dcoadd.dtype.names:
                d[tag] = dcoadd[tag]

        return d

    def get_default_epoch_fit_data(self):
        d = self._make_epoch_struct()
        return d


class MetacalDeconvolver(Deconvolver):

    def __init__(self,*args,**kw):
        super(MetacalDeconvolver,self).__init__(*args,**kw)
        self._types=['noshear','1p','1m','2p','2m']

        self._set_measurer_class()

    def _make_shears(self):
        step=self['metacal_pars'].get('step',0.01)

        shears={}
        for t in self._types:
            if t == 'noshear':
                shears[t]=None
            else:
                if t=='1p':
                    shears[t] = ngmix.Shape( step, 0.0)
                elif t=='1m':
                    shears[t] = ngmix.Shape(-step, 0.0)
                elif t=='2p':
                    shears[t] = ngmix.Shape(0.0,   step)
                elif t=='2m':
                    shears[t] = ngmix.Shape(0.0,  -step)

        return shears

    def _do_deconv(self, mb_obs_list):



        if self['deconv_pars']['trim_images']:
            mfitter=self.boot.get_max_fitter()
            mres=mfitter.get_result()
            cen=mres['pars'][0:0+2].copy()
            mb_obs_list=self._trim_images(mb_obs_list, cen)

        res = self._do_metacal_deconv(mb_obs_list)

        return res

    def _set_measurer_class(self):
        import deconv
        dpars=self['deconv_pars']
        weight_type=dpars.get('weight_type','ksigma')
        if weight_type=='ksigma':
            self._measure_class=deconv.measure.ObsKSigmaMoments
        elif weight_type=='ksigmac':
            self._measure_class=deconv.measure.ObsKSigmaMomentsC
        elif weight_type=='gauss':
            self._measure_class=deconv.measure.ObsGaussMoments
        else:
            raise NotImplementedError("bad weight type: '%s'" % weight_type)

    def _do_metacal_deconv(self, mb_obs_list):
        import deconv

        doplots=False

        types=['noshear','1p','1m','2p','2m']

        shears = self._make_shears()


        dpars={}
        dpars.update(self['deconv_pars'])

        if 'sigma_weight_factor' in dpars:
            sigma_weight=self._get_sigma_weight(mb_obs_list)
        else:
            sigma_weight=dpars.pop('sigma_weight')

        print("sigma weight:",sigma_weight)

        moments = self._measure_class(
            mb_obs_list,
            sigma_weight,
            **dpars
        )


        try:

            flags=0
            res={}

            for type in shears:
                moments.go(shear=shears[type],doplots=doplots)
                res[type]=moments.get_result()

                flags |= res[type]['flags']

            if doplots:
                if 'q'==raw_input('hit a key: '):
                    stop

        except deconv.DeconvRangeError as err:
            raise BootGalFailure(str(err))

        res['flags'] = flags
        return res


    def _get_sigma_weight(self, mbo):
        """
        we need to use the same weight function for each
        epoch/band.  We could use the largest psf size,
        or smallest, or mean...
        """

        n=0
        ssum=0.0
        for obslist in mbo:
            for obs in obslist:
                ssum += numpy.sqrt(obs.psf.gmix.get_T()/2.0)
                n+=1

        scale = obs.jacobian.get_scale()

        sigma = ssum/n
        sigma_weight = sigma*self['deconv_pars']['sigma_weight_factor']
        """
        print(
            "psf sigma:",sigma,
            "pixels:",sigma/scale,
            "sigma_weight:",sigma_weight,
            "pixels:",sigma_weight/scale,
        )
        """
        return sigma_weight

    def _get_metacal_namer(self, type='noshear'):
        if type=='noshear':
            back=None
        else:
            back=type
        return Namer(back=back)

    def _copy_galaxy_result(self, allres, coadd):
        """
        Copy from the result dict to the output array
        """

        dindex=0
        data=self.data

        for type in self._types:
            res=allres[type]

            n=self._get_metacal_namer(type=type)

            data[n('dc_flags')][dindex] = res['flags']
            data[n('dc_orflags')][dindex] = res['orflags']
            data[n('dc_flags_band')][dindex] = res['flags_band']
            data[n('dc_orflags_band')][dindex] = res['orflags_band']

            data[n('T')][dindex] = res['T']
            data[n('e')][dindex] = res['e']
            data[n('wflux')][dindex] = res['wflux']
            data[n('wflux_band')][dindex] = res['wflux']

            for nn in ['T_err','e_cov','flux_s2n','s2n_w']:
                if nn in res:
                    data[n(nn)][dindex] = res[nn]

            if 'flux_s2n' in res:
                print("    s2n_w: %.3g flux_s2n: %.3g" % (res['s2n_w'],res['flux_s2n']))

    def _make_struct(self,coadd):
        """
        make the output structure
        """
        num = 1
        dt = self._get_fit_data_dtype(coadd)
        data = numpy.zeros(num,dtype=dt)

        n=self._get_namer('psf', coadd)

        data[n('flags')] = NO_ATTEMPT
        for nn in ['flux','flux_err']:
            data[n(nn)] = DEFVAL

        n=self._get_namer('', coadd)

        for nn in ['mask_frac','psfrec_T','psfrec_g']:
            data[n(nn)] = DEFVAL

        # the deconvolve parameters
        for type in self._types:
            n=self._get_metacal_namer(type=type)

            for nn in ['dc_flags','dc_orflags','dc_flags_band','dc_orflags_band']: 
                data[n(nn)] = NO_ATTEMPT

            for nn in ['T','e','wflux','wflux_band','s2n_w','flux_s2n']:
                data[n(nn)] = DEFVAL

            for nn in ['T_err','e_cov']:
                data[n(nn)] = PDEFVAL

        return data


    def _get_fit_data_dtype(self,coadd):
        dt=[]

        nband=self['nband']
        bshape=(nband,)
        simple_npars=5+nband

        n=self._get_namer('psf', coadd)

        dt += [(n('flags'),   'i4',bshape),
               (n('flux'),    'f8',bshape),
               (n('flux_err'),'f8',bshape),
               (n('flux_s2n'),'f8',bshape)]

        n=self._get_namer('', coadd)

        dt += [(n('nimage_use'),'i4',bshape)]

        dt += [(n('mask_frac'),'f8'),
               (n('psfrec_T'),'f8'),
               (n('psfrec_g'),'f8', 2)]

        if nband==1:
            fcov_shape=(nband,)
        else:
            fcov_shape=(nband,nband)

        for type in self._types:

            if type=='noshear':
                back=None
            else:
                back=type

            n=Namer(back=back)
            dt += [
                (n('dc_flags'),'i4'),
                (n('dc_orflags'),'i4'),
                (n('dc_flags_band'),'i4',bshape),
                (n('dc_orflags_band'),'i4',bshape),

                (n('T'),'f8'),
                (n('e'),'f8',2),

                (n('T_err'),'f8'),
                (n('e_cov'),'f8', (2,2) ),
                (n('flux_s2n'),'f8'),
                (n('s2n_w'),'f8'),

                (n('wflux'),'f8',bshape),
                (n('wflux_band'),'f8',bshape),
            ]

        return dt

    def _fill_nimage_used(self,res,coadd):
        nim_used = numpy.zeros(self['nband'],dtype='i4')
        for band in xrange(self['nband']):
            nim_used[band] = res['noshear']['nimage_use_band'][band]

        n=self._get_namer('', coadd)

        self.data[n('nimage_use')][0,:] = nim_used


class MetacalDeconvolverPSFBase(MetacalDeconvolver):

    def _do_metacal_deconv(self, mb_obs_list):
        import deconv

        doplots=False

        types=['noshear','1p','1m','2p','2m']

        shears = self._make_shears()

        dpars = self['deconv_pars']

        moments = deconv.measure.KSigmaMomentsPSFBase(
            mb_obs_list,
            **dpars
        )


        try:

            flags=0
            res={}

            for type in shears:
                moments.go(shear=shears[type],doplots=doplots)
                res[type]=moments.get_result()

                flags |= res[type]['flags']

            if doplots:
                if 'q'==raw_input('hit a key: '):
                    stop

        except deconv.DeconvRangeError as err:
            raise BootGalFailure(str(err))

        res['flags'] = flags
        return res


def _trim_image(im, cen):

    drow_low = cen[0]
    drow_high = im.shape[0] - cen[0] - 1

    dcol_low = cen[1]
    dcol_high = im.shape[1] - cen[1] - 1

    drow = min(drow_low,drow_high)
    dcol = min(dcol_low,dcol_high)

    frlow=cen[0]-drow
    frhigh=cen[0]+drow
    fclow=cen[1]-dcol
    fchigh=cen[1]+dcol


    rlow=int(frlow)
    rhigh=int(frhigh)
    clow=int(fclow)
    chigh=int(fchigh)

    return im[
        rlow:rhigh+1,
        clow:chigh+1,
    ]
