from __future__ import print_function
import time
import numpy
import os
import scipy.stats

# local imports
from .defaults import DEFVAL, NO_ATTEMPT, \
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

def get_bootstrapper(obs, type='boot', **keys):
    from ngmix.bootstrap import Bootstrapper
    from ngmix.bootstrap import CompositeBootstrapper
    from ngmix.bootstrap import MaxMetacalBootstrapper

    if type=='boot':
        boot=Bootstrapper(obs, **keys)
    elif type=='composite':
        boot=CompositeBootstrapper(obs, **keys)
    elif type=='metacal':
        boot=MaxMetacalBootstrapper(obs, **keys)
    elif type=='metacal-analytic':
        boot=ngmix.bootstrap.MetacalAnalyticPSFBootstrapper(obs, **keys)
    else:
        raise ValueError("bad bootstrapper type: '%s'" % type)

    return boot

class NGMixBootFitter(BaseFitter):
    """
    Use an ngmix bootstrapper
    """
    def __init__(self,*args,**kw):
        super(NGMixBootFitter,self).__init__(*args,**kw)

    def _setup(self):
        """
        ngmix-specific setups for all ngmix fitters
        """
        # LM doesn't calculate a very good covariance matrix
        self['replace_cov'] = self.get('replace_cov',False)

        # in the fitters use log(flux) and log(T)
        self['use_logpars'] = self.get('use_logpars',False)

        # which models to fit
        self['fit_models'] = self.get('fit_models',list(self['model_pars'].keys()))

        # allow pre-selection based on psf flux
        self['min_psf_s2n'] = self.get('min_psf_s2n',-numpy.inf)

        # find the center and reset jacobians before doing model fits
        self['pre_find_center'] = self.get('pre_find_center',False)

        # do we normalize the psf to unity when doing the PSF mags?
        self['normalize_psf'] = self.get('normalize_psf',True)

        # how to mask unmodeled nbrs
        self['unmodeled_nbrs_masking_type'] = self.get('model_nrbs_unmodmask','nbrs-seg')

        # whether to use the coadd_ prefix for names when no me fit was done
        self['use_coadd_prefix'] = self.get('use_coadd_prefix',True)

        # add in estimate for intrinsic profile variance 
        self['intr_prof_var_fac'] = self.get('intr_prof_var_fac',0.0)

        # this actually means calculating things like
        # lensfit sens from samples or B&A PQR
        self['do_shear'] = self.get('do_shear',False)

    def _get_namer(self, model, coadd):
        if coadd and (self['fit_me_galaxy'] or self['use_coadd_prefix']):
            n = Namer('coadd_%s' % model)
        else:
            n = Namer(model)

        return n

    def get_models_for_checking(self):
        models = [modl for modl in self['fit_models']]
        pars = [modl+'_max_pars' for modl in self['fit_models']]
        covs = [modl+'_max_pars_cov' for modl in self['fit_models']]

        coadd_models = ['coadd_'+modl for modl in models]
        coadd_pars = ['coadd_'+modl for modl in pars]
        coadd_covs = ['coadd_'+modl for modl in covs]

        return models,pars,covs,coadd_models,coadd_pars,coadd_covs,5+self['nband']

    def _get_good_mb_obs_list(self,mb_obs_list):
        new_mb_obs_list = MultiBandObsList()
        for obs_list in mb_obs_list:
            new_obs_list = ObsList()
            for obs in obs_list:
                if obs.meta['flags'] == 0:
                    new_obs_list.append(obs)
            new_mb_obs_list.append(new_obs_list)
        new_mb_obs_list.update_meta_data(mb_obs_list.meta)
        new_mb_obs_list.update_meta_data({'old_mb_obs_list':mb_obs_list})
        return new_mb_obs_list

    def _add_extra_sim_noise(self, mb_obs_list):
        from numpy.random import normal

        target_noise=self['target_noise']
        target_var = target_noise**2
        target_ivar = 1.0/target_var

        print("    target noise:",target_noise)

        for obslist in mb_obs_list:
            for obs in obslist:

                obs.image_orig=obs.image.copy()
                obs.weight_orig=obs.weight.copy()

                wt=obs.weight
                w=numpy.where(wt > 0)
                if w[0].size > 0:

                    im = obs.image

                    extra_var_values = numpy.zeros(im.shape)
                    orig_var  = numpy.zeros(im.shape)
                    orig_var[w] = 1.0/wt[w]

                    wbad=numpy.where(orig_var > target_var)
                    if wbad[0].size > 0:
                        raise ValueError("attempting to set target var "
                                         "smaller than existing var")

                    # is zero where weight <= 0
                    extra_var_values[w] = target_var - orig_var[w]
                    extra_noise_values = numpy.sqrt(extra_var_values)

                    noise_image = normal(loc=0.0, scale=1.0, size=im.shape)
                    noise_image *= extra_noise_values
                    im += noise_image
                
                    wt[w] = target_ivar

                    obs.image = im
                    obs.weight = wt


    def __call__(self,mb_obs_list,coadd=False,nbrs_fit_data=None,nbrs_meta_data=None):
        """
        fit the obs list
        """

        if 'target_noise' in self:
            self._add_extra_sim_noise(mb_obs_list)

        # only fit stuff that is not flagged
        new_mb_obs_list = self._get_good_mb_obs_list(mb_obs_list)
        self.new_mb_obs_list = new_mb_obs_list
        self.mb_obs_list = mb_obs_list
        
        # FIXME - removed if, this might have sid effects?
        #if 'fit_data' not in mb_obs_list.meta:
        mb_obs_list.update_meta_data({'fit_data':self._make_struct(coadd)})
        self.data = mb_obs_list.meta['fit_data']

        if self['make_plots']:
            self.plot_dir = './%d-plots' % new_mb_obs_list.meta['id']
            if not os.path.exists(self.plot_dir):
                os.makedirs(self.plot_dir)

        fit_flags = 0
        for model in self['fit_models']:
            print('    fitting: %s' % model)

            if self['model_nbrs'] and nbrs_fit_data is not None:
                # put back nbrs info if needed
                if nbrs_meta_data is not None:
                    self._restore_nbrs_meta_data(mb_obs_list,nbrs_meta_data,coadd=coadd)
                    
                # render nbrs
                self._render_nbrs(model,new_mb_obs_list,coadd,nbrs_fit_data)

            model_flags, boot = self._guess_and_run_boot(model,
                                                         new_mb_obs_list,
                                                         coadd,
                                                         nbrs_fit_data=nbrs_fit_data)

            fit_flags |= model_flags

            # fill the epoch data
            self._fill_epoch_data(mb_obs_list,boot.mb_obs_list)

            # fill in PSF stats in data rows
            if (model_flags & PSF_FIT_FAILURE) == 0:
                self._do_psf_stats(mb_obs_list,coadd)
            else:
                break

        self._fill_nimage_used(mb_obs_list,boot.mb_obs_list,coadd)

        if self['model_nbrs']:
            self._fill_nbrs_data(mb_obs_list)

        return fit_flags

    def _guess_and_run_boot(self,model,new_mb_obs_list,coadd,nbrs_fit_data=None):
        n=self._get_namer(model, coadd)

        guess = None
        guess_errs = None
        guess_TdbyTe = 1.0

        # should always be True for mof runs. But when just subtracting
        # neighbors one ususally does not want to guess from the nbrs, since
        # mof will generally use all bands
        guess_from_nbrs=self.get('guess_from_nbrs',True)
        if nbrs_fit_data is not None and guess_from_nbrs:

            ind = new_mb_obs_list.meta['cen_ind']

            if (nbrs_fit_data[n('flags')][ind] == 0
                and nbrs_fit_data['flags'][ind] == 0):

                guess = nbrs_fit_data[n('pars')][ind]

                # lots of pain to get good guesses...
                # the ngmix ParsGuesser does this
                #    for pars 0 through 3 inclusive - uniform between -width to +width
                #    for pars 4 through the end - guess = pars*(1+width*uniform(low=-1,high=1))
                # thus for pars 4 through the end, I divide the error by the pars so that guess is
                #  between 1-frac_err to 1+frac_err where frac_err = err/pars
                # I also scale the errors by scale
                scale = 0.5

                # get the errors (cov in this case)
                guess_errs = numpy.diag(nbrs_fit_data[n('max_pars_cov')][ind]).copy()

                #if less than zero, set to zero
                w, = numpy.where(guess_errs < 0.0)
                if w.size > 0:
                    guess_errs[w[:]] = 0.0

                # take sqrt
                guess_errs = numpy.sqrt(guess_errs)*scale

                # get pars to scale by
                # don't divide by zero! - if zero set to 0.1 (default val in ngmix)
                w, = numpy.where(guess == 0.0)
                guess_scale = guess.copy()
                if w.size > 0:
                    guess_scale[w] = 0.1
                w = numpy.arange(4,guess.size,1)

                # final equation - need sqrt then apply scale and then divide by pars
                guess_errs[w[:]] = guess_errs[w]/numpy.abs(guess_scale[w])

                # don't guess to wide for the shear
                if guess_errs[2] > 0.1:
                    guess_errs[2] = 0.1

                if guess_errs[3] > 0.1:
                    guess_errs[3] = 0.1

                print_pars(guess,front='    guess pars:  ')
                print_pars(guess_errs,front='    guess errs:  ')

                if model == 'cm':
                    guess_TdbyTe = nbrs_fit_data[n('TdByTe')][ind]

        if model == 'cm':
            return self._run_boot(model,new_mb_obs_list,coadd,
                                  guess_TdbyTe=guess_TdbyTe,
                                  guess=guess,
                                  guess_widths=guess_errs)
        else:
            return self._run_boot(model,new_mb_obs_list,coadd,
                                  guess=guess,
                                  guess_widths=guess_errs)

    def _fill_nbrs_data(self,mb_obs_list):
        nd = len(mb_obs_list.meta['nbrs_ids'])

        # if no nbrs, return
        if nd == 0:
            return

        for band,obs_list in enumerate(mb_obs_list):
            for obs in obs_list:
                # if we use this obs, grab nbrs
                if obs.meta['flags'] == 0:
                    od = self.get_default_nbrs_data(nd+1)

                    # do nbrs
                    for i in xrange(nd):
                        od['nbr_id'][i] = mb_obs_list.meta['nbrs_ids'][i]
                        od['nbr_flags'][i] = obs.meta['nbrs_flags'][i]

                        if od['nbr_flags'][i] == 0 and obs.meta['nbrs_psfs'][i].has_gmix():
                            od['nbr_jac_row0'][i] = obs.meta['nbrs_jacs'][i].get_cen()[0]
                            od['nbr_jac_col0'][i] = obs.meta['nbrs_jacs'][i].get_cen()[1]

                            od['nbr_jac_dudrow'][i] = obs.meta['nbrs_jacs'][i].get_dudrow()
                            od['nbr_jac_dudcol'][i] = obs.meta['nbrs_jacs'][i].get_dudcol()

                            od['nbr_jac_dvdrow'][i] = obs.meta['nbrs_jacs'][i].get_dvdrow()
                            od['nbr_jac_dvdcol'][i] = obs.meta['nbrs_jacs'][i].get_dvdcol()

                            od['nbr_psf_fit_pars'][i,:] = obs.meta['nbrs_psfs'][i].get_gmix().get_full_pars()
                        else:
                            od['nbr_flags'][i] |= NBR_HAS_NO_PSF_FIT

                    # add in cen as "nbr" of self
                    i = nd
                    od['nbr_id'][i] = mb_obs_list.meta['id']

                    jac = obs.get_jacobian()
                    od['nbr_jac_row0'][i] = jac.get_cen()[0]
                    od['nbr_jac_col0'][i] = jac.get_cen()[1]
                    od['nbr_jac_dudrow'][i] = jac.get_dudrow()
                    od['nbr_jac_dudcol'][i] = jac.get_dudcol()
                    od['nbr_jac_dvdrow'][i] = jac.get_dvdrow()
                    od['nbr_jac_dvdcol'][i] = jac.get_dvdcol()

                    flags=0
                    if obs.has_psf_gmix():
                        od['nbr_psf_fit_pars'][i,:] = obs.get_psf_gmix().get_full_pars()
                    else:
                        flags |= NBR_HAS_NO_PSF_FIT

                    od['nbr_flags'][i] = flags

                    # add to metadata
                    obs.update_meta_data({'nbrs_data':od})

    def _restore_nbrs_meta_data(self,mb_obs_list,nbrs_meta_data,coadd=False):
        """
        restore the nbrs info - reverse of _fill_nbrs_data function
        """
        nd = len(mb_obs_list.meta['nbrs_ids'])
        
        # if no nbrs, return
        if nd == 0:
            return
        
        cen_id = mb_obs_list.meta['id']

        for band,obs_list in enumerate(mb_obs_list):
            band_num = obs_list.meta['band_num']
            
            for obs in obs_list:
                cutout_index = obs.meta['cutout_index']
                
                # if we use this obs, grab nbrs
                if obs.meta['flags'] == 0:
                    # do central
                    q, = numpy.where((nbrs_meta_data['id'] == cen_id) &
                                     (nbrs_meta_data['nbr_id'] == cen_id) &
                                     (nbrs_meta_data['band_num'] == band_num) &
                                     (nbrs_meta_data['cutout_index'] == cutout_index))                    
                    
                    if len(q) != 1:
                        raise ValueError('cen not found in nbrs_meta_data during restore!'\
                                             ' - cen_id = %d, band = %d, cutout_index = %d' \
                                             % (cen_id,band_num,cutout_index))
                    
                    cen_flags = nbrs_meta_data['nbr_flags'][q[0]]
                    if cen_flags == 0:
                        cen_jac = Jacobian(row=nbrs_meta_data['nbr_jac_row0'][q[0]],
                                       col=nbrs_meta_data['nbr_jac_col0'][q[0]],
                                       dudrow=nbrs_meta_data['nbr_jac_dudrow'][q[0]],
                                       dudcol=nbrs_meta_data['nbr_jac_dudcol'][q[0]],
                                       dvdrow=nbrs_meta_data['nbr_jac_dvdrow'][q[0]],
                                       dvdcol=nbrs_meta_data['nbr_jac_dvdcol'][q[0]])
                        cen_psf_gmix = GMix(pars=nbrs_meta_data['nbr_psf_fit_pars'][q[0],:])
                        cen_psf_obs = Observation(numpy.zeros((1,1)),gmix=cen_psf_gmix)
                    else:
                        cen_psf_obs = Observation(numpy.zeros((1,1)))
                        cen_jac = None
                        
                    obs.update_meta_data({'cen_flags':cen_flags,
                                          'cen_jac':cen_jac,
                                          'cen_psf':cen_psf_obs})
                        
                    # do nbrs
                    nbrs_flags = []
                    nbrs_psfs = []
                    nbrs_jacs = []
                    for i in xrange(nd):
                        # get
                        nbrs_id = mb_obs_list.meta['nbrs_ids'][i]
                        
                        # find the nbr
                        q, = numpy.where((nbrs_meta_data['id'] == cen_id) &
                                         (nbrs_meta_data['nbr_id'] == nbrs_id) &
                                         (nbrs_meta_data['band_num'] == band_num) &
                                         (nbrs_meta_data['cutout_index'] == cutout_index))
                        
                        if len(q) != 1:
                            raise ValueError('more than one nbr or no nbr found in nbrs_meta_data during restore!'\
                                                 ' - cen_id = %d, nbrs_id = %d, band = %d, cutout_index = %d' \
                                                 % (cen_id,nbrs_id,band_num,cutout_index))
                        
                        nbrs_flags.append(nbrs_meta_data['nbr_flags'][q[0]])
                        if nbrs_meta_data['nbr_flags'][q[0]] == 0:
                            jac = Jacobian(row=nbrs_meta_data['nbr_jac_row0'][q[0]],
                                           col=nbrs_meta_data['nbr_jac_col0'][q[0]],
                                           dudrow=nbrs_meta_data['nbr_jac_dudrow'][q[0]],
                                           dudcol=nbrs_meta_data['nbr_jac_dudcol'][q[0]],
                                           dvdrow=nbrs_meta_data['nbr_jac_dvdrow'][q[0]],
                                           dvdcol=nbrs_meta_data['nbr_jac_dvdcol'][q[0]])
                            nbrs_jacs.append(jac)
                            psf_gmix = GMix(pars=nbrs_meta_data['nbr_psf_fit_pars'][q[0],:])
                            psf_obs = Observation(numpy.zeros((1,1)),gmix=psf_gmix)
                            nbrs_psfs.append(psf_obs)
                        else:
                            nbrs_psfs.append(Observation(numpy.zeros((1,1))))
                            nbrs_jacs.append(None)
                            
                    obs.update_meta_data({'nbrs_jacs':nbrs_jacs,
                                          'nbrs_psfs':nbrs_psfs,
                                          'nbrs_flags':nbrs_flags})
                        
    def _render_nbrs(self,model,mb_obs_list,coadd,nbrs_fit_data):
        """
        render nbrs
        """

        print('    rendering nbrs')

        if len(mb_obs_list.meta['nbrs_inds']) == 0:
            return

        n=self._get_namer(model, coadd)

        pars_name = 'max_pars'
        if n(pars_name) not in nbrs_fit_data.dtype.names:
            pars_name = 'pars'

        pars_tag = n(pars_name)
        assert pars_tag in nbrs_fit_data.dtype.names

        if 'max' in pars_tag:
            fit_flags_tag = n('max_flags')
        else:
            fit_flags_tag = n('flags')

        if model == 'cm':
            fracdev_tag = n('fracdev')
            TdByTe_tag = n('TdByTe')
        else:
            fracdev_tag=None
            TdByTe_tag=None

        cen_ind = mb_obs_list.meta['cen_ind']
        for band,obs_list in enumerate(mb_obs_list):
            for obs in obs_list:
                if obs.meta['flags'] != 0:
                    continue
                
                # get properties of cen 
                if 'cen_flags' in obs.meta and obs.meta['cen_flags'] == 0:
                    cen_jac = obs.meta['cen_jac']
                    if obs.meta['cen_psf'].has_gmix():
                        cen_psf_gmix = obs.meta['cen_psf'].get_gmix()
                    else:
                        cen_psf_gmix = None
                else:
                    if obs.has_psf_gmix():
                        cen_psf_gmix = obs.get_psf_gmix()
                    else:
                        cen_psf_gmix = None                    
                    cen_jac = obs.get_jacobian()
                cen_seg = obs.seg
                
                # get props of nbrs
                nbrs_inds = mb_obs_list.meta['nbrs_inds']
                nbrs_flags = obs.meta['nbrs_flags']
                nbrs_psf_gmixes = []
                for psf_obs in obs.meta['nbrs_psfs']:
                    if psf_obs is not None and psf_obs.has_gmix():
                        nbrs_psf_gmixes.append(psf_obs.get_gmix())
                    else:
                        nbrs_psf_gmixes.append(None)
                nbrs_jacs = obs.meta['nbrs_jacs']
                
                # call the nbrs code
                cenim, nbrs_imgs, nbrs_masks = RenderNGmixNbrs._render_nbrs(model, band, 
                                                                            obs.image.shape,
                                                                            cen_ind, cen_psf_gmix, cen_jac, cen_seg,
                                                                            nbrs_inds, nbrs_flags,
                                                                            nbrs_jacs, nbrs_psf_gmixes,
                                                                            pars_tag, fit_flags_tag, nbrs_fit_data,
                                                                            unmodeled_nbrs_masking_type=self['unmodeled_nbrs_masking_type'],
                                                                            verbose=True,
                                                                            fracdev_tag=fracdev_tag,TdByTe_tag=TdByTe_tag)
                
                # do central
                if cenim is not None:
                    sub_nbrs_from_cenim = False
                else:                
                    cenim = obs.image_orig.copy()
                    sub_nbrs_from_cenim = True

                if self['intr_prof_var_fac'] > 0.0:
                    varim = numpy.zeros_like(cenim)  
                    if not sub_nbrs_from_cenim:
                        varim += self['intr_prof_var_fac']*cenim*cenim

                # now do nbrs
                nbrsim = numpy.zeros_like(cenim)
                for curr_nbrsim in nbrs_imgs:
                    if curr_nbrsim is not None:
                        nbrsim += curr_nbrsim
                        
                        if self['intr_prof_var_fac'] > 0.0:
                            varim += self['intr_prof_var_fac']*curr_nbrsim*curr_nbrsim

                masked_pix = numpy.zeros_like(cenim)
                masked_pix[:,:] = 1.0
                for msk in nbrs_masks:
                    if msk is not None:
                        masked_pix *= msk
                            
                # get total image and adjust central if needed
                if sub_nbrs_from_cenim:
                    cenim -= nbrsim                    
                totim = cenim + nbrsim
                
                if self['model_nbrs_method'] == 'subtract':
                    obs.image = obs.image_orig - nbrsim
                elif self['model_nbrs_method'] == 'frac':
                    frac = numpy.zeros_like(totim)
                    frac[:,:] = 1.0
                    msk = totim > 0.0
                    frac[msk] = cenim[msk]/totim[msk]
                    obs.image = obs.image_orig*frac
                else:
                    assert False,'nbrs model method %s not implemented!' % self['model_nbrs_method']
                    
                # mask unmodeled nbrs
                new_weight = obs.weight_orig.copy()
                new_weight *= masked_pix

                if self['intr_prof_var_fac'] > 0.0:
                    qnz = numpy.where(new_weight != 0.0)
                    if qnz[0].size > 0:
                        new_weight[qnz] = 1.0/(1.0/new_weight[qnz] + varim[qnz])

                obs.weight = new_weight
                    
                if self['make_plots']:
                    self._plot_nbrs_model(band,model,obs,nbrsim,cenim,coadd)

    def _plot_nbrs_model(self,band,model,obs,nbrsim,cenim,coadd):
        """
        plot nbrs model
        """
        if coadd:
            ptype='coadd'
        else:
            ptype='mb'

        obj_id = obs.meta['id']
        ptype = '%s-%s-%s' % (ptype,model,'max')
        title='%d %s' % (obj_id,ptype)

        def plot_seg(seg):
            seg_new = seg.copy()
            seg_new = seg_new.astype(float)
            uvals = numpy.unique(seg)
            if len(uvals) > 1:
                mval = 1.0*(len(uvals)-1)
                ind = 1.0
                for uval in uvals:
                    if uval > 0:
                        qx,qy = numpy.where(seg == uval)
                        seg_new[qx[:],qy[:]] = ind/mval
                        ind += 1
            else:
                seg_new[:,:] = 1.0

            return seg_new

        icut_cen = obs.meta['icut']
        
        totim = nbrsim + cenim
        
        import images
        import biggles
        width = 1920
        height = 1200
        biggles.configure('screen','width', width)
        biggles.configure('screen','height', height)
        tab = biggles.Table(2,3)
        tab.title = title

        tab[0,0] = images.view(obs.image_orig,title='original image',show=False,nonlinear=0.075)
        tab[0,1] = images.view(nbrsim,title='models of nbrs',show=False,nonlinear=0.075)

        tab[0,2] = images.view(plot_seg(obs.seg),title='seg map',show=False)

        tab[1,0] = images.view(obs.image,title='corrected image',show=False,nonlinear=0.075)
        msk = totim != 0
        frac = numpy.zeros(totim.shape)
        frac[msk] = cenim[msk]/totim[msk]
        tab[1,1] = images.view(frac,title='fraction of flux due to central',show=False)
        tab[1,2] = images.view(obs.weight,title='weight map',show=False)

        try:
            if icut_cen > 0:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-icut%d.png' % (ptype,band,icut_cen))
            else:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-coadd.png' % (ptype,band))
            print("        making plot %s" % fname)
            tab.write_img(1920,1200,fname)
        except:
            print("        caught error plotting nbrs")
            pass

        # for testing
        if False:
            import fitsio
            if icut_cen > 0:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-icut%d.fits' % (ptype,band,icut_cen))
            else:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-coadd.fits' % (ptype,band))
            print("        making file %s" % fname)
            fitsio.write(fname,obs.image,clobber=True)

            if icut_cen > 0:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-icut%d-nbrs.fits' % (ptype,band,icut_cen))
            else:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-coadd-nbrs.fits' % (ptype,band))
            print("        making file %s" % fname)
            fitsio.write(fname,nbrsim,clobber=True)

            if icut_cen > 0:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-icut%d-cen.fits' % (ptype,band,icut_cen))
            else:
                fname = os.path.join(self.plot_dir,'%s-nbrs-band%d-coadd-cen.fits' % (ptype,band))
            print("        making file %s" % fname)
            fitsio.write(fname,cenim,clobber=True)

    def _fill_epoch_data(self,mb_obs_list,new_mb_obs_list):
        print('    filling PSF data')
        for band,obs_list in enumerate(mb_obs_list):
            for obs in obs_list:
                used = False
                res = None

                if obs.meta['flags'] != 0:
                    obs.update_meta_data({'fit_flags':obs.meta['flags']})
                    ed = self._make_epoch_struct()
                    ed['psf_fit_flags'] = obs.meta['flags']
                    obs.update_meta_data({'fit_data':ed})
                    continue

                if obs.meta['flags'] == 0 and obs.has_psf():
                    psf_obs = obs.get_psf()

                    ed = self._make_epoch_struct()
                    ed['npix'] = obs.image.size
                    ed['wsum'] = obs.weight.sum()
                    ed['wmax'] = obs.weight.max()
                    ed['psf_counts'] = psf_obs.image.sum()

                    if 'fitter' in psf_obs.meta:
                        res = obs.get_psf().meta['fitter'].get_result()
                        ed['psf_fit_flags'] = res['flags']

                    if psf_obs.has_gmix():
                        used = True
                        psf_gmix = psf_obs.get_gmix()
                        g1,g2,T = psf_gmix.get_g1g2T()
                        pars = psf_gmix.get_full_pars()

                        ed['psf_fit_g'][0,0] = g1
                        ed['psf_fit_g'][0,1] = g2
                        ed['psf_fit_T'] = T
                        ed['psf_fit_pars'] = pars

                    if obs in new_mb_obs_list[band] and used:
                        obs.update_meta_data({'fit_flags':0})
                    else:
                        obs.update_meta_data({'fit_flags':PSF_FIT_FAILURE})

                    obs.update_meta_data({'fit_data':ed})

    def _do_psf_stats(self,mb_obs_list,coadd):
        print('    doing PSF stats')

        n=self._get_namer('', coadd)

        Tsum = 0.0
        g1sum = 0.0
        g2sum = 0.0
        wsum = 0.0
        wrelsum = 0.0
        npix = 0.0
        did_one = False
        did_one_max = False
        for band,obs_list in enumerate(mb_obs_list):
            for obs in obs_list:

                fdata=obs.meta.get('fit_data',None)
                if (obs.meta['flags'] == 0
                        and fdata is not None
                        and fdata['psf_fit_flags'][0] == 0):

                    assert obs.meta['fit_flags'] == 0
                    assert obs.get_psf().has_gmix()

                    if fdata['wmax'][0] > 0.0:
                        did_one_max = True
                        npix += fdata['npix'][0]
                        wrelsum += fdata['wsum'][0]/fdata['wmax'][0]

                    did_one = True

                    wsum += fdata['wsum'][0]
                    Tsum += fdata['wsum'][0]*fdata['psf_fit_T'][0]

                    g1sum += fdata['wsum'][0]*fdata['psf_fit_g'][0,0]
                    g2sum += fdata['wsum'][0]*fdata['psf_fit_g'][0,1]

        if did_one_max:
            self.data[n('mask_frac')][0] = 1.0 - wrelsum/npix

        if did_one:
            self.data[n('psfrec_g')][0,0] = g1sum/wsum
            self.data[n('psfrec_g')][0,1] = g2sum/wsum
            self.data[n('psfrec_T')][0] = Tsum/wsum

    def _fill_nimage_used(self,mb_obs_list,new_mb_obs_list,coadd):
        nim = numpy.zeros(self['nband'],dtype='i4')
        nim_used = numpy.zeros_like(nim)
        for band in xrange(self['nband']):
            nim[band] = len(mb_obs_list[band])
            nim_used[band] = len(new_mb_obs_list[band])

        n=self._get_namer('', coadd)

        self.data[n('nimage_use')][0,:] = nim_used

    def _get_bootstrapper(self, model, mb_obs_list):
        """
        get the bootstrapper for fitting psf through galaxy
        """
        
        find_cen=self.get('pre_find_center',False)
        if model == 'cm':
            fracdev_prior=self['model_pars']['cm']['fracdev_prior']
            boot=get_bootstrapper(mb_obs_list,
                                  type='composite',
                                  fracdev_prior=fracdev_prior,
                                  find_cen=find_cen,
                                  **self)
        else:
            boot=get_bootstrapper(mb_obs_list, find_cen=find_cen, **self)

        return boot

    def _run_boot(self, model, mb_obs_list, coadd, guess=None, **kwargs):
        """
        run a boot strapper
        """

        flags=0
        boot=self._get_bootstrapper(model,mb_obs_list)
        self.boot=boot

        n=self._get_namer('psf', coadd)

        try:
            self._fit_psfs(coadd)


            flags |= self._fit_psf_flux(coadd)

            if flags == 0:
                dindex = 0
                s2n = self.data[n('flux_s2n')]
                max_s2n = numpy.nanmax(s2n)
                if max_s2n < self['min_psf_s2n']:
                    flags |= LOW_PSF_FLUX

            if flags == 0:
                try:
                    self._fit_galaxy(model,coadd,guess=guess,**kwargs)
                    self._copy_galaxy_result(model,coadd)
                    self._print_galaxy_result()
                except (BootGalFailure,GMixRangeError) as err:
                    print("    galaxy fitting failed: %s" % err)
                    flags = GAL_FIT_FAILURE

        except BootPSFFailure:
            print("    psf fitting failed")
            flags = PSF_FIT_FAILURE

        return flags, boot

    def _fit_psf_flux(self,coadd):

        print("    fitting psf flux")
        self.boot.fit_gal_psf_flux(normalize_psf=self['normalize_psf'])

        res=self.boot.get_psf_flux_result()

        n=self._get_namer('psf', coadd)

        flagsall=0
        for band in xrange(self['nband']):
            flags=res['flags'][band]
            flagsall |= flags

            flux=res['psf_flux'][band]
            flux_err=res['psf_flux_err'][band]

            self.data[n('flags')][0,band] = flags
            self.data[n('flux')][0,band] = flux
            self.data[n('flux_err')][0,band] = flux_err

            if flux_err > 0:
                s2n = flux/flux_err
                self.data[n('flux_s2n')][0,band] = s2n
                tup=(band,flux,flux_err,s2n)
                print("        psf flux(%s): %g +/- %g s2n: %g" % tup)
            else:
                print("        psf flux(%s): %g +/- %g" % (band,flux,flux_err))

        if flagsall != 0:
            # we only propagate this bit to the main 'flags' field
            flagsall = PSF_FLUX_FIT_FAILURE

        return flagsall

    def _fit_psfs(self,coadd, boot=None):
        """
        fit the psf model to every observation's psf image
        """

        print('    fitting the PSFs')

        if boot is None:
            boot=self.boot

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

        if (self['make_plots']
            and (('made_psf_plots' not in self.mb_obs_list.meta) or
                 ('made_psf_plots' in self.mb_obs_list.meta and
                  self.mb_obs_list.meta['made_psf_plots'] == False)) ):

            self.mb_obs_list.update_meta_data({'made_psf_plots':True})
            for band,obs_list in enumerate(boot.mb_obs_list):
                for obs in obs_list:
                    psf_obs = obs.get_psf()
                    if psf_obs.has_gmix():
                        self._do_psf_plot(psf_obs,obs.meta['id'],band,obs.meta['band_id'],coadd)

    def _do_psf_plot(self,obs,obs_id,band,band_id,coadd):
        """
        make residual plots for psf
        """
        import images

        title='%d band: %s' % (obs_id, band)
        if coadd:
            title='%s coadd' % title
        else:
            title='%s %d' % (title,band_id)

        im=obs.image

        gmix = obs.get_gmix()
        model_im=gmix.make_image(im.shape, jacobian=obs.jacobian)
        modflux=model_im.sum()
        if modflux <= 0:
            log.log("        psf model flux too low: %f" % modflux)
            return

        model_im *= ( im.sum()/modflux )

        ims = 1e3/numpy.max(im)
        plt=images.compare_images(im*ims, model_im*ims,
                                  label1='psf', label2='model',
                                  show=False, nonlinear=0.075)
        plt.title=title

        if coadd:
            fname=os.path.join(self.plot_dir,'%d-psf-resid-band%d-coadd.png' % (obs_id,band))
        else:
            fname=os.path.join(self.plot_dir,'%d-psf-resid-band%d-%d.png' % (obs_id,band,band_id))

        print("        making plot %s" % fname)
        plt.write_img(1920,1200,fname)

    def _fit_galaxy(self, model, coadd, guess=None, **kwargs):
        """
        over-ride for different fitters
        """
        raise RuntimeError("over-ride me")

    def _plot_resids(self, obj_id, fitter, model, coadd, fitter_type):
        """
        make plots
        """
        import images

        if coadd:
            ptype='coadd'
        else:
            ptype='mb'

        ptype = '%s-%s-%s' % (ptype,model,fitter_type)
        title='%d %s' % (obj_id,ptype)
        try:
            res_plots = None
            if fitter_type != 'isample':
                res_plots=fitter.plot_residuals(title=title)
            if res_plots is not None:
                for band, band_plots in enumerate(res_plots):
                    for icut, plt in enumerate(band_plots):
                        fname=os.path.join(self.plot_dir,'%d-%s-resid-band%d-im%d.png' % (obj_id,ptype,band,icut))
                        print("        making plot %s" % fname)
                        plt.write_img(1920,1200,fname)

        except GMixRangeError as err:
            print("        caught error plotting resid: %s" % str(err))

    def _plot_trials(self, obj_id, fitter, model, coadd, fitter_type, wgts):
        """
        make plots
        """
        if coadd:
            ptype='coadd'
        else:
            ptype='mb'

        ptype = '%s-%s-%s' % (ptype,model,fitter_type)
        title='%d %s' % (obj_id,ptype)

        try:
            if fitter_type == 'isample':
                pdict=fitter.make_plots(title=title)
            else:
                pdict=fitter.make_plots(title=title,weights=wgts)
        
            pdict['trials'].aspect_ratio=1.5
            pdict['wtrials'].aspect_ratio=1.5
            
            trials_png=os.path.join(self.plot_dir,'%d-%s-trials.png' % (obj_id,ptype))
            wtrials_png=os.path.join(self.plot_dir,'%d-%s-wtrials.png' % (obj_id,ptype))
            
            print("        making plot %s" % trials_png)
            pdict['trials'].write_img(1200,1200,trials_png)
            
            print("        making plot %s" % wtrials_png)
            pdict['wtrials'].write_img(1200,1200,wtrials_png)
        except:
            print("        caught error plotting trials")

        try:
            from .util import plot_autocorr
            trials=fitter.get_trials()
            plt=plot_autocorr(trials)
            plt.title=title
            fname=os.path.join(self.plot_dir,'%d-%s-autocorr.png' % (obj_id,ptype))
            print("        making plot %s" % fname)
            plt.write_img(1000,1000,fname)
        except:
            print("        caught error plotting autocorr")

    def _plot_images(self, obj_id, model, coadd):
        import images
        imlist = []
        titles = []
        for band,obs_list in enumerate(self.new_mb_obs_list):
            for obs in obs_list:
                imlist.append(obs.image*obs.weight)
                titles.append('band: %d %d' % (band,obs.meta['band_id']))
        if coadd:
            coadd_png=os.path.join(self.plot_dir,'%d-coadd-images.png' % (obj_id))
        else:
            coadd_png=os.path.join(self.plot_dir,'%d-mb-images.png' % (obj_id))
        plt=images.view_mosaic(imlist, titles=titles, show=False)
        print("        making plot %s" % coadd_png)
        plt.write_img(1200,1200,coadd_png)

    def _copy_galaxy_result(self, model, coadd):
        """
        Copy from the result dict to the output array
        """

        dindex=0

        res=self.gal_fitter.get_result()
        mres=self.boot.get_max_fitter().get_result()

        rres=self.boot.get_round_result()

        n=self._get_namer(model, coadd)

        data=self.data

        data[n('flags')][dindex] = res['flags']

        fname, Tname = self._get_lnames()

        if res['flags'] == 0:
            pars=res['pars']
            pars_cov=res['pars_cov']

            flux=pars[5:]
            flux_cov=pars_cov[5:, 5:]

            data[n('max_flags')][dindex] = mres['flags']
            data[n('max_pars')][dindex,:] = mres['pars']
            data[n('max_pars_cov')][dindex,:,:] = mres['pars_cov']

            data[n('pars')][dindex,:] = pars
            data[n('pars_cov')][dindex,:,:] = pars_cov

            data[n('g')][dindex,:] = res['g']
            data[n('g_cov')][dindex,:,:] = res['g_cov']

            data[n('flags_r')][dindex]  = rres['flags']
            data[n('s2n_r')][dindex]    = rres['s2n_r']
            data[n(Tname+'_r')][dindex] = rres['pars'][4]
            data[n('psf_T_r')][dindex]  = rres['psf_T_r']

            for sn in ['s2n_w','chi2per','dof']:
                data[n(sn)][dindex] = res[sn]

            if self['do_shear']:
                if 'g_sens' in res:
                    data[n('g_sens')][dindex,:] = res['g_sens']

                if 'R' in res:
                    data[n('P')][dindex] = res['P']
                    data[n('Q')][dindex,:] = res['Q']
                    data[n('R')][dindex,:,:] = res['R']

            for f in ['fracdev',
                      'fracdev_noclip',
                      'fracdev_err',
                      'TdByTe',
                      'TdByTe_noclip']:
                if f in res:
                    data[n(f)][dindex] = res[f]

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

    def get_epoch_fit_data_dtype(self):
        npars = self.get_num_pars_psf()
        dt=[('npix','i4'),
            ('wsum','f8'),
            ('wmax','f8'),
            ('psf_fit_flags','i4'),
            ('psf_counts','f8'),
            ('psf_fit_g','f8',2),
            ('psf_fit_T','f8'),
            ('psf_fit_pars','f8',npars)]

        return dt

    def _make_epoch_struct(self,num=1):
        dt = self.get_epoch_fit_data_dtype()

        epoch_data = numpy.zeros(num, dtype=dt)
        epoch_data['npix'] = DEFVAL
        epoch_data['wsum'] = DEFVAL
        epoch_data['wmax'] = DEFVAL
        epoch_data['psf_fit_flags'] = NO_ATTEMPT
        epoch_data['psf_counts'] = DEFVAL
        epoch_data['psf_fit_g'] = DEFVAL
        epoch_data['psf_fit_T'] = DEFVAL
        epoch_data['psf_fit_pars'] = DEFVAL

        return epoch_data

    def get_fit_data_dtype(self,me,coadd):
        dt = []
        if me:
            dt += self._get_fit_data_dtype(False)
        if coadd:
            dt += self._get_fit_data_dtype(True)
        return dt

    def _get_lnames(self):
        if self['use_logpars']:
            fname='log_flux'
            Tname='log_T'
        else:
            fname='flux'
            Tname='T'

        return fname, Tname

    def _get_all_models(self,coadd):
        """
        get all model names, includeing the coadd_ ones
        """
        self['fit_models'] = self.get('fit_models',list(self['model_pars'].keys()))

        models=[]
        #if coadd:
        if coadd and (self['fit_me_galaxy'] or self['use_coadd_prefix']):
            models = models + ['coadd_%s' % model for model in self['fit_models']]
        else:
            models = models + self['fit_models']

        return models

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

        fname, Tname=self._get_lnames()

        models=self._get_all_models(coadd)
        for model in models:
            n=Namer(model)
            np=simple_npars

            dt+=[(n('flags'),'i4'),
                 (n('pars'),'f8',np),
                 (n('pars_cov'),'f8',(np,np)),
                 (n('g'),'f8',2),
                 (n('g_cov'),'f8',(2,2)),

                 (n('max_flags'),'i4'),
                 (n('max_pars'),'f8',np),
                 (n('max_pars_cov'),'f8',(np,np)),

                 (n('s2n_w'),'f8'),
                 (n('chi2per'),'f8'),
                 (n('dof'),'f8'),

                 (n('flags_r'),'i4'),
                 (n('s2n_r'),'f8'),
                 (n(Tname+'_r'),'f8'),
                 (n('psf_T_r'),'f8')]

            if self['do_shear']:
                dt += [(n('g_sens'), 'f8', 2)]

            if 'cm' in model:
                dt += [(n('fracdev'),'f4'),
                       (n('fracdev_noclip'),'f4'),                       
                       (n('fracdev_err'),'f4'),
                       (n('TdByTe'),'f4'),
                       (n('TdByTe_noclip'),'f4')]

        return dt

    def _make_struct(self,coadd):
        """
        make the output structure
        """
        num = 1
        dt = self._get_fit_data_dtype(coadd=coadd)
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

        fname, Tname=self._get_lnames()

        models=self._get_all_models(coadd)
        for model in models:
            n=Namer(model)

            data[n('flags')] = NO_ATTEMPT

            data[n('pars')] = DEFVAL
            data[n('pars_cov')] = DEFVAL

            data[n('g')] = DEFVAL
            data[n('g_cov')] = DEFVAL

            data[n('s2n_w')] = DEFVAL
            data[n('chi2per')] = DEFVAL

            data[n('max_flags')] = NO_ATTEMPT
            data[n('max_pars')] = DEFVAL
            data[n('max_pars_cov')] = DEFVAL

            data[n('flags_r')] = NO_ATTEMPT
            data[n('s2n_r')] = DEFVAL
            data[n(Tname+'_r')] = DEFVAL
            data[n('psf_T_r')] = DEFVAL

            if self['do_shear']:
                data[n('g_sens')] = DEFVAL

            if 'cm' in model:
                data[n('fracdev')] = DEFVAL
                data[n('fracdev_noclip')] = DEFVAL
                data[n('fracdev_err')] = DEFVAL
                data[n('TdByTe')] = DEFVAL
                data[n('TdByTe_noclip')] = DEFVAL

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

    def get_nbrs_data_dtype(self):
        npars = self.get_num_pars_psf()
        dt = [('nbr_id','i8'),
              ('nbr_flags','i4'),
              ('nbr_jac_row0','f8'),
              ('nbr_jac_col0','f8'),
              ('nbr_jac_dudrow','f8'),
              ('nbr_jac_dudcol','f8'),
              ('nbr_jac_dvdrow','f8'),
              ('nbr_jac_dvdcol','f8'),
              ('nbr_psf_fit_pars','f8',(npars,)),
              ]
        
        return dt
        
    def get_default_nbrs_data(self,n=1):
        dt = self.get_nbrs_data_dtype()
        d = numpy.zeros(n,dtype=dt)

        for tag in d.dtype.names:
            d[tag] = DEFVAL
        d['nbr_flags'] = 0

        return d

class MaxNGMixBootFitter(NGMixBootFitter):
    def _fit_max(self, model, guess=None, boot=None, **kwargs):
        """
        do a maximum likelihood fit

        note prior applied during
        """

        if boot is None:
            boot=self.boot

        max_pars=self['max_pars']
        prior=self['model_pars'][model]['prior']

        guess_widths = kwargs.get('guess_widths',None)

        # now with prior
        if model == 'cm' and 'guess_TdbyTe' in kwargs:
            boot.fit_max(model,
                         max_pars,
                         prior=prior,
                         ntry=max_pars['ntry'],
                         guess=guess,
                         guess_TdbyTe=kwargs['guess_TdbyTe'],
                         guess_widths=guess_widths)
        else:
            boot.fit_max(model,
                         max_pars,
                         prior=prior,
                         ntry=max_pars['ntry'],
                         guess=guess,
                         guess_widths=guess_widths)


        if self['replace_cov']:
            print("        replacing cov")
            cov_pars=self['cov_pars']
            boot.try_replace_cov(cov_pars)

    def _fit_galaxy(self, model, coadd, guess=None, **kwargs):
        self._fit_max(model,guess=guess,**kwargs)

        rpars=self['round_pars']
        self.boot.set_round_s2n(fitter_type=rpars['fitter_type'])

        self.gal_fitter=self.boot.get_max_fitter()

        if self['make_plots']:
            self._plot_resids(self.new_mb_obs_list.meta['id'],
                              self.boot.get_max_fitter(),
                              model, coadd, 'max')
            self._plot_images(self.new_mb_obs_list.meta['id'], model, coadd)

class ISampNGMixBootFitter(MaxNGMixBootFitter):
    def _setup(self):
        super(ISampNGMixBootFitter,self)._setup()

        # verbose True so we can see isamp output
        self['verbose'] = True

    def _fit_galaxy(self, model, coadd, guess=None, **kwargs):
        self._fit_max(model,guess=guess)
        self._do_isample(model)
        self._add_shear_info(model)

        self.gal_fitter=self.boot.get_isampler()

        if self['make_plots']:
            self._plot_resids(self.new_mb_obs_list.meta['id'],
                              self.boot.get_max_fitter(),
                              model,
                              coadd,
                              'max')

            self._plot_images(self.new_mb_obs_list.meta['id'], model, coadd)

            self._plot_trials(self.new_mb_obs_list.meta['id'],
                              self.boot.get_isampler(),
                              model,
                              coadd,
                              'isample',
                              self.boot.get_isampler().get_iweights())

    def _do_isample(self, model):
        """
        run isample on the bootstrapper
        """
        ipars=self['isample_pars']
        prior=self['model_pars'][model]['prior']
        self.boot.isample(ipars, prior=prior)

        rpars=self['round_pars']
        self.boot.set_round_s2n(fitter_type=rpars['fitter_type'])

    def _add_shear_info(self, model):
        """
        add shear information based on the gal_fitter
        """

        boot=self.boot
        max_fitter=boot.get_max_fitter()
        sampler=boot.get_isampler()

        # this is the full prior
        prior=self['model_pars'][model]['prior']
        g_prior=prior.g_prior

        iweights = sampler.get_iweights()
        samples = sampler.get_samples()
        g_vals=samples[:,2:2+2]

        res=sampler.get_result()

        # keep for later if we want to make plots
        self.weights=iweights

        # we are going to mutate the result dict owned by the sampler
        stats = max_fitter.get_fit_stats(res['pars'])
        res.update(stats)

        ls=ngmix.lensfit.LensfitSensitivity(g_vals,
                                            g_prior,
                                            weights=iweights,
                                            remove_prior=True)
        g_sens = ls.get_g_sens()
        g_mean = ls.get_g_mean()

        res['g_sens'] = g_sens
        res['nuse'] = ls.get_nuse()

    def _copy_galaxy_result(self, model, coadd):
        super(ISampNGMixBootFitter,self)._copy_galaxy_result(model,coadd)

        dindex=0
        res=self.gal_fitter.get_result()

        n = self._get_namer(model, coadd)

        if res['flags'] == 0:
            for f in ['efficiency','neff']:
                self.data[n(f)][dindex] = res[f]

    def _get_fit_data_dtype(self,coadd):
        dt=super(ISampNGMixBootFitter,self)._get_fit_data_dtype(coadd)

        for model in self._get_all_models(coadd):
            n=Namer(model)
            dt += [(n('efficiency'),'f4'),
                   (n('neff'),'f4')]

        return dt

    def _make_struct(self,coadd):
        d = super(ISampNGMixBootFitter,self)._make_struct(coadd)

        for model in self._get_all_models(coadd):
            n=Namer(model)
            d[n('efficiency')] = DEFVAL
            d[n('neff')] = DEFVAL

        return d

class MetacalNGMixBootFitter(MaxNGMixBootFitter):
    def _setup(self):
        from .priors import set_priors

        super(MetacalNGMixBootFitter,self)._setup()

        self['nrand'] = self.get('nrand',1)
        if self['nrand'] is None:
            self['nrand']=1

        self['use_original_weight'] = self.get('use_original_weight',False)

        defpars={'step':0.01}
        self['metacal_pars'] = self.get('metacal_pars',defpars)
        self['metacal_pars']['nband'] = self['nband']
        print("metacal pars:")
        pprint(self['metacal_pars'])

        types=self['metacal_pars'].get('types',None)
        if types is None:
            self['metacal_pars']['types']=ngmix.metacal.METACAL_TYPES

        if 'model_pars' in self['metacal_pars']:
            print("setting separate metacal prior")
            set_priors(self['metacal_pars'])

    def _get_metacal_bootstrapper(self, model, mb_obs_list):
        """
        get the bootstrapper for fitting psf through galaxy
        """

        if 'analytic_psf' in self['metacal_pars']:
            type='metacal-analytic'
        else:
            type='metacal'

        find_cen=self.get('pre_find_center',False)
        boot=get_bootstrapper(
            mb_obs_list,
            find_cen=find_cen,
            type=type,
            **self
        )

        return boot


    def _fit_galaxy(self, model, coadd, guess=None,**kwargs):
        mb_obs_list = self.boot.mb_obs_list

        super(MetacalNGMixBootFitter,self)._fit_galaxy(model,
                                                       coadd,
                                                       guess=guess,
                                                       **kwargs)


        self.mcal_boot=self._do_metacal(model, self.boot)

        metacal_res = self.mcal_boot.get_metacal_result()

        res=self.gal_fitter.get_result()
        res.update(metacal_res)

    def _get_metacal_stuff(self, model):
        metacal_pars=self['metacal_pars']
        model_pars=metacal_pars.get('model_pars',None)
        if model_pars is not None:
            #print("getting model and priors from metacal pars")
            model=list(model_pars.keys())[0]
            prior=model_pars[model]['prior']
        else:
            #print("using main prior")
            prior=self['model_pars'][model]['prior']

        psf_pars=metacal_pars.get('psf_pars',None)
        if psf_pars is None:
            #print("using main psf model")
            psf_pars=self['psf_pars']
        else:
            #print("getting psf pars from metacal")
            pass

        psf_fit_pars = psf_pars.get('fit_pars',None)
        return metacal_pars, model, prior, psf_pars, psf_fit_pars

    def _do_metacal(self,
                    model,
                    boot,
                    metacal_obs=None):
        """
        the basic fitter for this class
        """

        if self['replace_bad_pixels']:
            print("        replacing masked pixels")
            boot.replace_masked_pixels(inplace=True)
        else:
            print("        not replacing masked pixels")

        metacal_pars, model, prior, psf_pars, psf_fit_pars = \
                self._get_metacal_stuff(model)
        max_pars=self['max_pars']
        Tguess=boot.mb_obs_list[0][0].psf.gmix.get_T()

        # new bootstrapper for metacal
        mcal_boot=self._get_metacal_bootstrapper(
            model,
            boot.mb_obs_list,
        )

        try:

            mcal_boot.fit_metacal(
                psf_pars['model'],
                model,
                max_pars,
                Tguess,
                psf_fit_pars=psf_fit_pars,
                prior=prior,
                ntry=max_pars['ntry'],
                metacal_pars=self['metacal_pars'],
            )

        except BootPSFFailure as err:
            # the _run_boot code catches this one
            raise BootGalFailure(str(err))

        return mcal_boot


    def _get_fit_data_dtype(self,coadd):
        dt=super(MetacalNGMixBootFitter,self)._get_fit_data_dtype(coadd)

        dt_mcal = self._get_metacal_dtype(coadd)
        dt += dt_mcal
        return dt

    def _copy_galaxy_result(self, model, coadd):
        super(MetacalNGMixBootFitter,self)._copy_galaxy_result(model,coadd)

        dindex=0
        d=self.data
        res=self.gal_fitter.get_result()

        n=self._get_namer(model, coadd)

        mcal_flags=res['mcal_flags']
        d['mcal_flags'][dindex] = mcal_flags

        if res['mcal_flags'] != 0:
            print("    metacal flags set:",mcal_flags)
            d[n('flags')][dindex] |= METACAL_FAILURE

        for type in self['metacal_pars']['types']:
            tres=res[type]
            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            d['mcal_pars%s' % back][dindex] = tres['pars']
            d['mcal_g%s' % back][dindex] = tres['g']
            d['mcal_g_cov%s' % back][dindex] = tres['g_cov']

            # no error for T_r
            T     = tres['pars'][4]
            T_err = numpy.sqrt(tres['pars_cov'][4,4])

            d['mcal_T%s' % back][dindex] = T
            d['mcal_T_err%s' % back][dindex] = T_err

            d['mcal_T_r%s' % back][dindex] = tres['T_r']

            d['mcal_s2n_r%s' % back][dindex] = tres['s2n_r']

            if type=='noshear':
                for p in ['pars_cov','gpsf','Tpsf']:
                    name='mcal_%s' % p
                    d[name][dindex] = tres[p]

    def _get_metacal_dtype(self, coadd):

        nband=self['nband']
        simple_npars=5+nband
        np=simple_npars

        dt=[('mcal_flags','i8')]

        models=self._get_all_models(coadd)
        if len(models) > 1:
            raise RuntimeError("for metacal, only fit one model and "
                               "either coadd or me")

        for type in self['metacal_pars']['types']:

            if type=='noshear':
                back=''
            else:
                back='_%s' % type

            dt += [
                ('mcal_g%s' % back,'f8',2),
                ('mcal_g_cov%s' % back,'f8',(2,2)),  # might be used for weights
                ('mcal_pars%s' % back,'f8',np),
            ]

            if type=='noshear':
                dt += [
                    ('mcal_pars_cov','f8',(np,np)),
                    ('mcal_gpsf','f8',2),
                    ('mcal_Tpsf','f8'),
                ]

            dt += [
                ('mcal_T%s' % back,'f8'),
                ('mcal_T_err%s' % back,'f8'),
                ('mcal_T_r%s' % back,'f8'),
                ('mcal_s2n_r%s' % back,'f8'),
            ]

        self.mcal_flist = [d[0] for d in dt]

        return dt


    def _make_struct(self,coadd):
        data=super(MetacalNGMixBootFitter,self)._make_struct(coadd)

        models=self._get_all_models(coadd)
        for model in models:
            for f in self.mcal_flist:
                data[f] = DEFVAL

        return data





