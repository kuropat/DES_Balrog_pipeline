"""
code for fitting
"""

import numpy
from .util import UtterFailure

def get_fitter_class(ftype):
    """
    return the fitter class for this obs
    """
    from .__init__ import FITTERS

    cftype = ftype.lower()
    assert cftype in FITTERS,'could not find fitter class %s' % cftype

    return FITTERS[cftype]

class BaseFitter(dict):
    """
    abstract base class for fitting
    """
    def __init__(self,conf):
        self.update(conf)

        self._setup()

    def _setup(self):
        """
        put setups, particularly checking for some stuff that is optional
        only steups that would be used for *any* fitter go here
        """
        pass

    def get_models_for_checking(self):
        """
        return a list of model pars for checking for convergence in mof
        also return a second list of the cov estimate for each parameter

            model_name,pars_tag,cov_tag,coadd_model_name,coadd_pars_tag,coadd_cov_tag,npars = fitter.get_models_for_checking()

        """
        raise NotImplementedError("get_models_for_checking method of BaseFitter must be defined in subclass.")

    def get_fit_data_dtype(self,me,coadd):
        """
        returns a numpy dtype for the galaxy fit data as a list

        if me == True, then return with me columns
        if coadd == True, then return with coadd columns
        if both true, return with both sets of columns

        For example
            return [('x','f8'),('x','f8')]
        """
        raise NotImplementedError("get_fit_data_dtype method of BaseFitter must be defined in subclass.")

    def get_default_fit_data(self,me,coadd):
        """
        returns the default values for a line in the fit data table

        me and coadd behave as in get_fit_data_dtype
        """
        raise NotImplementedError("get_default_fit_data method of BaseFitter must be defined in subclass.")

    def get_epoch_fit_data_dtype(self):
        """
        returns a numpy dtype for the galaxy per epoch fit data as a list
        For example
            return [('x','f8'),('y','f8')]
        """
        raise NotImplementedError("get_epoch_fit_data_dtype method of BaseFitter must be defined in subclass.")

    def get_default_epoch_fit_data(self):
        """
        returns a default line in the per galaxy epoch fit data
        """
        raise NotImplementedError("get_default_epoch_fit_data method of BaseFitter must be defined in subclass.")

    def get_nbrs_data_dtype(self):
        """
        returns a numpy dtype for the galaxy per epoch nbrs data as a list
        For example
            return [('x','f8'),('y','f8')]
        """
        raise NotImplementedError("get_nbrs_data_dtype method of BaseFitter must be defined in subclass.")

    def get_default_nbrs_data(self):
        """
        returns a default line in the per galaxy epoch nbrs data
        """
        raise NotImplementedError("get_default_nbrs_data method of BaseFitter must be defined in subclass.")

    def __call__(self,mb_obs_list,coadd=False,make_epoch_data=True,nbrs_fit_data=None,make_plots=False):
        """
        do fit of single obs list

        method be used like this

            flags = Fitter(mb_obs_list,coadd=False)

        It should also flag any indv. obs that was not used in the fit with a
        non-zero value in the 'fit_flags' field of the meta data dict.

        The fit data should be returned as a numpy array in the mb_obs_list
        meta data dict in 'fit_data'.

        The per epoch fit data should be returned as a numpy array in
        'fit_data' in the indv. obs's meta data dicts.

        These numpy arrays should have the same dtype as returned by the
        methods get_fit_data_dtype and get_epoch_fit_data_dtype respectively.

        If something goes really wrong, the method can rise the UtterFailure exception above.

        The coadd keyword will be set to true if this obs is a coadd.

        If make_epoch_data is False, then epoch_fit_data can be set to None in
        the meta data dict. Even if this value is not None if make_epoch_data
        is False, this value will be ignored be the calling routine.

        If make_plots is set, fitter should make some plots.

        Nbrs Modeling
        -------------
        if nbrs_fit_data is not None, then all of the nbrs for this obejct should be modeled.

        Each mb_obs_list in mb_obs_lists has four meta data fields that are used to do this

            nbrs_inds: index into nbrs_fit_data
            nbrs_psfs: psf obs for each nbr
            nbrs_jacs: jaobians for each nbr
            nbrs_flags: only render nbrs with flags == 0
            cen_ind: index into nbrs_fit_data

        The idea is that the rows in nbrs_fit_data correspond to the objects in mb_obs_lists in order.

        These fields are lists that can be iterated over in parallel like this

            for ind,psf,jac,flags in zip(nbrs_inds,nbrs_psfs,nbrs_jacs,nbrs_flags):
                # render and model the nbr for a given obs

        The original image and weight map in each obs are here

            obs.image_orig
            obs.weight_orig

        Any meta data (like PSF models or WCS jacobians) for the nbrs should be set in the nbrs_data field of
        the meta data of each obs in the mb_obs_list. This field is a numpy array of dtype returned by 
        'get_nbrs_data_dtype'.

        """
        raise NotImplementedError("__call__ method of BaseFitter must be defined in subclass.")
