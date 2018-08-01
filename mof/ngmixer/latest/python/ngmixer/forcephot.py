from __future__ import print_function
import time
import numpy
import os
import scipy.stats
import fitsio
from .ngmixing import NGMixer

# local imports
from .defaults import (
    DEFVAL, PDEFVAL, NO_ATTEMPT,
    PSF_FIT_FAILURE, GAL_FIT_FAILURE,
    LOW_PSF_FLUX, PSF_FLUX_FIT_FAILURE,
    NBR_HAS_NO_PSF_FIT, METACAL_FAILURE,
    FORCEPHOT_BAD_MODEL,
    FORCEPHOT_FAILURE,
)
from .fitting import BaseFitter
from .bootfit import NGMixBootFitter
from .util import Namer, print_pars

# ngmix imports
import ngmix
from ngmix import Observation, ObsList, MultiBandObsList, GMixRangeError
from ngmix.gexceptions import BootPSFFailure, BootGalFailure
from ngmix.gmix import GMix
from ngmix.jacobian import Jacobian

from pprint import pprint

class ForcedPhotometryNGMixer(NGMixer):
    def __init__(self, *args, **kw):

        self.models_file=kw.get('models_file',None)
        if self.models_file is None:
            raise ValueError("send models_file= for forced photometry")

        super(ForcedPhotometryNGMixer,self).__init__(*args, **kw)


    def go(self):
        """
        also setup models etc.
        """
        self._check_forced_photometry()
        self._set_fp_models(self.models_file)
        super(ForcedPhotometryNGMixer,self).go()

    def _check_forced_photometry(self):
        fp=self.get('forced_photometry',None)
        if fp is None:
            raise ValueError("forced_photometry must be in config")

        if 'model' not in fp:
            raise ValueError("model must be in forced_photometry config")
        if 'min_T' not in fp:
            raise ValueError("min_T must be in forced_photometry config")
        
        if not self['center_psf']:
            raise ValueError("center_psf must be True for forced photometry")

    def fit_all_obs_lists(self,
                          coadd_mb_obs_list,
                          mb_obs_list,
                          nbrs_fit_data=None,
                          nbrs_meta_data=None,
                          make_epoch_data=True,
                          make_nbrs_data=True):
        """
        add the model information to the observation lists.  Use the meta
        data for this
        """

        
        gm, flags = self._make_model()

        if flags != 0:
            return flags

        coadd_mb_obs_list.update_meta_data({'forcephot_model':gm.copy()})
        mb_obs_list.update_meta_data({'forcephot_model':gm.copy()})

        flags = super(ForcedPhotometryNGMixer,self).fit_all_obs_lists(
            coadd_mb_obs_list,
            mb_obs_list,
            nbrs_fit_data=None,
            nbrs_meta_data=None,
            make_epoch_data=True,
            make_nbrs_data=True,
        )
        return flags

    def _make_model(self):
        """
        make the gaussian mixture model
        """
        model_data, flags = self._extract_current_model_data()
        if flags != 0:
            return None, flags

        model=self['forced_photometry']['model']
        pars = model_data['%s_pars' % model][0:6].copy()
        pars[-1] = 1.0

        if model=='cm':
            gm = ngmix.gmix.GMixCM(
                model_data['cm_fracdev'],
                model_data['cm_TdByTe'],
                pars,
            )
        else:
            gm = ngmix.gmix.GMixModel(
                model_data['%s_pars' % model],
                pars,
            )

        return gm, 0

    def _extract_current_model_data(self):
        """
        extract the current model paramters
        """
        index = self.curr_fofindex

        data = self.forcephot_models[index]

        if data['flags'] != 0:
            return None, FORCEPHOT_BAD_MODEL
        else:
            return data, 0

    def _set_fp_models(self, models_file):
        """
        load models for use in forced photometry
        """
        if models_file is None:
            raise ValueError("for forced photometry you must send a "
                             "models file")

        model_name = self['forced_photometry']['model']
        columns = [
            'id',
            'flags',
            '%s_pars' % model_name,
        ]
        if model_name == 'cm':
            columns += ['cm_fracdev', 'cm_TdByTe']

        if self.fof_range is not None:
            beg = self.fof_range[0]
            end = self.fof_range[1]+1
        else:
            beg=0
            end = self.imageio.get_num_fofs()

        with fitsio.FITS(models_file) as fits:
            if end > fits['model_fits'].get_nrows():
                raise ValueError("requested fof range is larger than "
                                 "models file: [%d,%d]" % (beg,end-1))
            self.forcephot_models = fits['model_fits'][columns][beg:end]

    def _get_nbrs_dtype(self):
        return None




#class ForcedPhotometryFitter(BaseFitter):
class ForcedPhotometryFitter(NGMixBootFitter):
    #def __init__(self, conf):
    #    super(NGMixBootFitter,self).__init__(conf)

    #def _setup(self):
    #    """
    #    put setups, particularly checking for some stuff that is optional
    #    only steups that would be used for *any* fitter go here
    #    """
    #    pass
    def _set_models(self):
        pass
        #self['fit_models'] = self.get('fit_models',list(self['model_pars'].keys()))

    def __call__(self,
                 mb_obs_list,
                 coadd=False,
                 **kw):
        """
        do forced photometry
        """

        # only fit stuff that is not flagged
        new_mb_obs_list = self._get_good_mb_obs_list(mb_obs_list)
        self.new_mb_obs_list = new_mb_obs_list
        self.mb_obs_list = mb_obs_list
        
        mb_obs_list.update_meta_data({'fit_data':self._make_struct(coadd)})
        self.data = mb_obs_list.meta['fit_data']

        # do the forced photometry for each band
        psf_flags, fit_flags = self._do_fits(
            new_mb_obs_list,
            coadd,
        )

        # fill the epoch data
        self._fill_epoch_data(mb_obs_list,new_mb_obs_list)

        self._fill_nimage_used(mb_obs_list,new_mb_obs_list,coadd)
        self._calc_mask_frac(mb_obs_list, coadd)

        return fit_flags

    def _set_priors(self):
        """
        do nothing
        """
        pass

    def _do_fits(self, mb_obs_list, coadd):
        """
        use psf as a template, measure flux (linear)
        """

        psf_flags = self._do_psf_fits(mb_obs_list, coadd)
        #psf_flags=0

        nband = len(mb_obs_list)

        flags=[]
        flux = numpy.zeros(nband) - 9999.0
        flux_err = numpy.zeros(nband) + 9999.0

        n=self._get_namer('fp', coadd)

        gm_model = mb_obs_list.meta['forcephot_model']

        flagsall=0
        for band,obs_list in enumerate(mb_obs_list):
            res = self._fit_one_band(obs_list, gm_model)
            flagsall |= res['flags']

            self.data[n('flags')][0,band] = res['flags']

            if res['flags'] == 0:

                self.data[n('flux')][0,band] = res['flux']
                self.data[n('flux_err')][0,band] = res['flux_err']
                if res['flux_err'] > 0:
                    s2n = res['flux']/res['flux_err']
                    self.data[n('flux_s2n')][0,band] = s2n
     
            else:
                print("failed to fit template flux for band",band)

        if flagsall != 0:
            flagsall = FORCEPHOT_FAILURE 

        return psf_flags, flagsall

    def _do_psf_fits(self, mb_obs_list, coadd):
        # use a bootstrapper just for the psf fits
        self.boot = self._get_bootstrapper("gauss", mb_obs_list)
        psf_flags = 0
        try:
            self._fit_psfs(coadd)
        except BootPSFFailure as err:
            print("    psf fitting failed: %s" % str(err))
            psf_flags = PSF_FIT_FAILURE

        return psf_flags

    def _fit_one_band(self, obs_list, gm_model):
        """
        fit one of the bands
        """
        if len(obs_list)==0:
            raise BootPSFFailure("no epochs for band %d" % i)

        try:
            fitter = self._get_template_fitter(obs_list, gm_model)
            fitter.go()
            res=fitter.get_result()
        except GMixRangeError as err:
            print(err)
            res = {'flags':FORCEPHOT_FAILURE}

        return res


    def _get_template_fitter(self, obs_list, gm_model):

        T = gm_model.get_T()
        if T < self['forced_photometry']['min_T']:
            print("    T too small, doing pure psf flux")
            model = None
        else:
            model = gm_model.make_galsim_object()

        return self._get_image_psf_template_fitter(obs_list, model)

    def _get_image_psf_template_fitter(self, obs_list, model):
        if not obs_list[0].has_psf():
            raise RuntimeError("you need to have psfs for image-psf fits")

        fitter=ngmix.galsimfit.GalsimTemplateFluxFitter(
            obs_list,
            model=model,
            normalize_psf=self['normalize_psf'],
        )

        return fitter



    def _make_struct(self,coadd):
        """
        make the output structure
        """
        num = 1
        dt = self._get_fit_data_dtype(coadd)
        data = numpy.zeros(num,dtype=dt)

        n=self._get_namer('fp', coadd)

        data[n('flags')] = NO_ATTEMPT
        data[n('flux')] = DEFVAL
        data[n('flux_err')] = DEFVAL
        data[n('flux_s2n')] = DEFVAL

        n=self._get_namer('', coadd)

        data[n('mask_frac')] = DEFVAL
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

        # fp stands for ForcedPhotometry
        n=self._get_namer('fp', coadd)

        dt += [
            #('flags','i4'),
            ('nimage_use','i4',bshape),
            ('mask_frac','f8'),
            (n('flags'),   'i4',bshape),
            (n('flux'),    'f8',bshape),
            (n('flux_err'),'f8',bshape),
            (n('flux_s2n'),'f8',bshape),
        ]

        return dt

    def get_default_nbrs_data(self):
        return None


    # copied from bootfit
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


    def _get_namer(self, model, coadd):
        if coadd and (self['fit_me_galaxy'] or self['use_coadd_prefix']):
            n = Namer('coadd_%s' % model)
        else:
            n = Namer(model)

        return n


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


    def get_default_epoch_fit_data(self):
        d = self._make_epoch_struct()
        return d

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

    def _fill_nimage_used(self,mb_obs_list,new_mb_obs_list,coadd):
        nim = numpy.zeros(self['nband'],dtype='i4')
        nim_used = numpy.zeros_like(nim)
        for band in xrange(self['nband']):
            nim[band] = len(mb_obs_list[band])
            nim_used[band] = len(new_mb_obs_list[band])

        n=self._get_namer('', coadd)

        self.data[n('nimage_use')][0,:] = nim_used

    def _calc_mask_frac(self,mb_obs_list,coadd):
        print('    doing PSF stats')

        n=self._get_namer('', coadd)

        wrelsum = 0.0
        npix = 0.0
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

        if did_one_max:
            self.data[n('mask_frac')][0] = 1.0 - wrelsum/npix
