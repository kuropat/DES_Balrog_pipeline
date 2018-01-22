#!/usr/bin/env python
from __future__ import print_function
import os
import ngmix
import fitsio
import numpy

def set_priors(conf):
    """
    Sets priors for each model.

    Currently only separable priors can be set.
    """

    from ngmix.joint_prior import PriorSimpleSep
    from ngmix.priors import ZDisk2D

    if 'model_pars' not in conf:
        print("no model_pars set, not setting priors")
        return

    g_prior_flat=ZDisk2D(1.0)
    model_pars=conf['model_pars']

    assert 'nband' in conf,'# of bands nband must be in config dict conf when setting priors'

    # set comps
    for model,params in model_pars.iteritems():

        if params is None:
            print("    no priors for model '%s'" % model)
            continue 

        print("loading prior for: %s" % model)

        params['cen_prior'] = get_cen_prior(params['cen'])
        params['g_prior'] = get_g_prior(params['g'])
        params['T_prior'] = get_T_prior(params['T'])
        params['counts_prior'] = get_counts_prior(params['counts'],
                                                  conf['nband'])

        if model == 'cm':
            params['fracdev_prior'] = get_fracdev_prior(params['fracdev'])

        print("    full")
        prior = PriorSimpleSep(params['cen_prior'],
                               params['g_prior'],
                               params['T_prior'],
                               params['counts_prior'])

        # for the exploration, for which we do not apply g prior during
        print("    gflat")
        gflat_prior = PriorSimpleSep(params['cen_prior'],
                                     g_prior_flat,
                                     params['T_prior'],
                                     params['counts_prior'])

        params['prior'] = prior
        params['gflat_prior'] = gflat_prior

def get_T_prior(params):

    typ=params['type']

    if typ == 'flat':
        pars=params['pars']
        prior = ngmix.priors.FlatPrior(*pars)

    elif typ=='TwoSidedErf':
        pars=params['pars']
        prior = ngmix.priors.TwoSidedErf(*pars)

    elif typ =='lognormal':

        mean=params['mean']
        sigma=params['sigma']
        prior = ngmix.priors.LogNormal(mean, sigma)

    elif typ=="cosmos_exp":
        prior = ngmix.priors.TPriorCosmosExp()

    elif typ=="cosmos_dev":
        prior = ngmix.priors.TPriorCosmosDev()

    elif typ=="gmixnd":
        prior = load_gmixnd(params)
    else:
        raise ValueError("bad T prior type: %s" % typ)

    return prior

def get_counts_prior(params, nband):
    typ=params['type']

    if typ=="gmixnd":
        if nband > 1:
            raise NotImplementedError("make gmixnd for counts "
                                      "work for N bands")
        prior_list = [load_gmixnd(params)]
    else:

        if typ == 'flat':
            pclass = ngmix.priors.FlatPrior
        elif typ=='TwoSidedErf':
            pclass = ngmix.priors.TwoSidedErf
        else:
            raise ValueError("bad counts prior type: %s" % typ)

        if params['repeat']:
            pars=params['pars']
            # we assume this is one that will be repeated
            prior_list = [pclass(*pars)]*nband
        else:
            pars=params['pars']
            # assume this is a list of lists
            prior_list=[]
            for tpars in pars:
                cp = pclass(*tpars)
                prior_list.append(cp)

            mess=("counts prior must be length "
                  "%d, got %d" % (nband,len(prior_list)) )
            assert len(prior_list) == nband,mess


    return prior_list

def get_fracdev_prior(params):
    if 'file' in params:
        fname=os.path.expanduser( params['file'] )
        fname=os.path.expandvars( fname )

        print("reading fracdev_prior:",fname)
        data = fitsio.read(fname)

        weights=data['weights']
        means=data['means']
        covars=data['covars']

    elif 'means' in params:
        means = numpy.array(params['means'])
        weights = numpy.array(params['weights'])
        covars= numpy.array(params['covars'])

    else:
        raise ValueError("set either the weights,means,covars or a "
                         "file name for the fracdev prior")

    if len(means.shape) == 1:
        means = means.reshape( (means.size, 1) )
    if len(covars.shape) == 1:
        covars = covars.reshape( (covars.size, 1, 1) )

    prior = ngmix.gmix.GMixND(weights,
                              means,
                              covars)
    return prior

def get_g_prior(params):
    typ=params['type']

    if typ=='cosmos-sersic':
        g_prior = ngmix.priors.make_gprior_cosmos_sersic(type='erf')
    elif typ=='cosmos-exp':
        g_prior = ngmix.priors.make_gprior_cosmos_exp()
    elif typ=='cosmos-dev':
        g_prior = ngmix.priors.make_gprior_cosmos_dev()
    elif typ =='ba':
        sigma=params['sigma']
        g_prior = ngmix.priors.GPriorBA(sigma)
    elif typ=='flat':
        g_prior=ngmix.priors.ZDisk2D(1.0)
    else:
        raise ValueError("implement gprior '%s'" % typ)

    return g_prior

def get_cen_prior(params):
    width=params['width']
    prior=ngmix.priors.CenPrior(0.0, 0.0, width, width)
    return prior

def load_gmixnd(spec):

    fname = os.path.expandvars(spec['file'])

    pdf=ngmix.gmix.GMixND(file=fname)

    if 'cov_factor' in spec:
        print("    using cov factor:",spec['cov_factor'])
        pdf.covars *= spec['cov_factor']

    return pdf
