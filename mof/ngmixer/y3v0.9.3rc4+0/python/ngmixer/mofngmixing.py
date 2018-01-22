#!/usr/bin/env python
from __future__ import print_function
import numpy
import time
import pprint
import os
import fitsio
import copy

# local imports
from . import imageio
from . import fitting
from . import files
from .ngmixing import NGMixer
from .defaults import DEFVAL,_CHECKPOINTS_DEFAULT_MINUTES
from .defaults import NO_ATTEMPT,NO_CUTOUTS,BOX_SIZE_TOO_BIG,IMAGE_FLAGS
from .defaults import MOF_SKIPPED_IN_CONV_CHECK, \
    MOF_NOT_CONVERGED, \
    MOF_NBR_NOT_CONVERGED, \
    MOF_NBR_BAD_FIT, \
    MOF_NBR_NOT_FIT, \
    MOF_NBR_SKIPPED_IN_CONV_CHECK, \
    MOF_FOFMEM_NOT_CONVERGED, \
    MOF_FOFMEM_BAD_FIT, \
    MOF_FOFMEM_NOT_FIT, \
    MOF_FOFMEM_SKIPPED_IN_CONV_CHECK, \
    MOF_NOT_CONVERGED, \
    MOF_SKIPPED_IN_CONV_CHECK

from .util import UtterFailure,Namer,print_pars
from .util import print_with_verbosity

class MOFNGMixer(NGMixer):
    def _set_defaults(self):
        super(MOFNGMixer,self)._set_defaults()
        self['mof']['write_convergence_data'] = self['mof'].get('write_convergence_data',False)

    def _get_models_to_check(self):
        me_models_to_check,me_pars_models_to_check,me_cov_models_to_check, \
            coadd_models_to_check,coadd_pars_models_to_check,coadd_cov_models_to_check, \
            npars = self.fitter.get_models_for_checking()

        models_to_check = []
        pars_models_to_check = []
        cov_models_to_check = []

        if self['fit_me_galaxy']:
            models_to_check.extend(me_models_to_check)
            pars_models_to_check.extend(me_pars_models_to_check)
            cov_models_to_check.extend(me_cov_models_to_check)

        if self['fit_coadd_galaxy']:
            models_to_check.extend(coadd_models_to_check)
            pars_models_to_check.extend(coadd_pars_models_to_check)
            cov_models_to_check.extend(coadd_cov_models_to_check)

        return models_to_check,pars_models_to_check,cov_models_to_check,npars

    def _check_convergence(self,foflen,itr,coadd_mb_obs_lists,mb_obs_lists):
        """
        check convergence of fits
        """

        models_to_check,pars_models_to_check,cov_models_to_check,npars = self._get_models_to_check()

        maxabs = numpy.zeros(npars,dtype='f8')
        maxabs[:] = -numpy.inf
        maxfrac = numpy.zeros(npars,dtype='f8')
        maxfrac[:] = -numpy.inf
        maxerr = numpy.zeros(npars,dtype='f8')
        maxerr[:] = -numpy.inf

        for fofind in xrange(foflen):
            print_with_verbosity('    fof obj: %ld' % (fofind+1),verbosity=1)

            for model,pars_model,model_cov in zip(models_to_check,pars_models_to_check,cov_models_to_check):
                if pars_model not in self.curr_data.dtype.names:
                    continue
                
                n = Namer(model)
                
                self.curr_data[n('mof_flags')][fofind] = 0                           
                self.curr_data[n('mof_num_itr')][fofind] = itr+1
                
                if self.curr_data['flags'][fofind] or self.prev_data['flags'][fofind]:
                    print('    skipping fof obj %s in convergence check' % (fofind+1))
                    self.curr_data[n('mof_flags')][fofind] = MOF_SKIPPED_IN_CONV_CHECK
                    continue

                old = self.prev_data[pars_model][fofind]
                new = self.curr_data[pars_model][fofind]
                absdiff = numpy.abs(new-old)
                absfracdiff = numpy.abs(new/old-1.0)
                abserr = numpy.abs((old-new)/numpy.sqrt(numpy.diag(self.curr_data[model_cov][fofind])))
                
                self.curr_data[n('mof_abs_diff')][fofind] = absdiff
                self.curr_data[n('mof_frac_diff')][fofind] = absfracdiff
                self.curr_data[n('mof_err_diff')][fofind] = abserr
                if numpy.all((absdiff <= self['mof']['maxabs_conv'])      | \
                             (absfracdiff <= self['mof']['maxfrac_conv']) | \
                             (abserr <= self['mof']['maxerr_conv'])):
                    self.curr_data[n('mof_flags')][fofind] = 0
                else:
                    self.curr_data[n('mof_flags')][fofind] = MOF_NOT_CONVERGED


                for i in xrange(npars):
                    if absdiff[i] > maxabs[i]:
                        maxabs[i] = copy.copy(absdiff[i])
                    if absfracdiff[i] > maxfrac[i]:
                        maxfrac[i] = copy.copy(absfracdiff[i])
                    if abserr[i] > maxerr[i]:
                        maxerr[i] = copy.copy(abserr[i])

                print_with_verbosity('        %s:' % model,verbosity=1)
                print_pars(old,        front='            old      ',verbosity=1)
                print_pars(new,        front='            new      ',verbosity=1)
                print_pars(absdiff,    front='            abs diff ',verbosity=1)
                print_pars(absfracdiff,front='            frac diff',verbosity=1)
                print_pars(abserr,     front='            err diff ',verbosity=1)


        fmt = "%8.3g "*len(maxabs)
        print("    max abs diff : "+fmt % tuple(maxabs))
        print("    max frac diff: "+fmt % tuple(maxfrac))
        print("    max err diff : "+fmt % tuple(maxerr))

        self.maxabs = maxabs
        self.maxfrac = maxfrac
        self.maxerr = maxerr                

        self._set_nbr_mof_flags(foflen,coadd_mb_obs_lists,mb_obs_lists)

        if numpy.all((maxabs <= self['mof']['maxabs_conv'])   | \
                     (maxfrac <= self['mof']['maxfrac_conv']) | \
                     (maxerr <= self['mof']['maxerr_conv'])):
            return True
        else:
            return False

    def _set_nbr_mof_flags(self,foflen,coadd_mb_obs_lists,mb_obs_lists):
        models_to_check,pars_models_to_check,cov_models_to_check,npars = self._get_models_to_check()

        for model,pars_model,model_cov in zip(models_to_check,pars_models_to_check,cov_models_to_check):
            if pars_model not in self.curr_data.dtype.names:
                continue

            n = Namer(model)
            
            any_not_conv = 0
            any_bad_fit = 0
            any_not_fit = 0
            any_skip_conv = 0

            for cen_ind in xrange(foflen):

                if mb_obs_lists[cen_ind].meta['obj_flags'] != 0:
                    any_not_fit = 1
                    
                if self.curr_data[n('mof_flags')][cen_ind]&MOF_NOT_CONVERGED != 0:
                    any_not_conv = 1
                    
                if self.curr_data['flags'][cen_ind]:
                    any_bad_fit = 1
                    
                if self.curr_data[n('mof_flags')][cen_ind]&MOF_SKIPPED_IN_CONV_CHECK != 0:
                    any_skip_conv = 1
                
                nbr_not_conv = 0
                nbr_bad_fit = 0
                nbr_not_fit = 0
                nbr_skip_conv = 0
                
                nbrs_inds = mb_obs_lists[cen_ind].meta['nbrs_inds']
                for nbrs_ind in nbrs_inds:
                    
                    if mb_obs_lists[nbrs_ind].meta['obj_flags'] != 0:
                        nbr_not_fit = 1
                    
                    if self.curr_data[n('mof_flags')][nbrs_ind]&MOF_NOT_CONVERGED != 0:
                        nbr_not_conv = 1
                                            
                    if self.curr_data['flags'][nbrs_ind]:
                        nbr_bad_fit = 1
                        
                    if self.curr_data[n('mof_flags')][nbrs_ind]&MOF_SKIPPED_IN_CONV_CHECK != 0:
                        nbr_skip_conv = 1
                    
                # set flag for nbrs lists
                if nbr_not_conv:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_NBR_NOT_CONVERGED
                    
                if nbr_bad_fit:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_NBR_BAD_FIT
                    
                if nbr_not_fit:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_NBR_NOT_FIT
                    
                if nbr_skip_conv:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_NBR_SKIPPED_IN_CONV_CHECK


            # set flags for the whole fof
            for cen_ind in xrange(foflen):
                if any_not_conv:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_FOFMEM_NOT_CONVERGED
                    
                if any_bad_fit:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_FOFMEM_BAD_FIT
                
                if any_not_fit:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_FOFMEM_NOT_FIT
                
                if any_skip_conv:
                    self.curr_data[n('mof_flags')][cen_ind] |= MOF_FOFMEM_SKIPPED_IN_CONV_CHECK
                
    def _set_default_data_for_fofind(self,fofind):
        for tag in self.default_data.dtype.names:
            self.curr_data[tag][fofind] = self.default_data[tag]
        self.curr_data['fofind'][fofind] = fofind

    def do_fits(self):
        """
        Fit all objects in our list
        """

        self.done = False

        print('doing fits')

        t0=time.time()
        num = 0
        numfof = 0
        numtot = self.imageio.get_num_fofs()

        print('fof index: %d:%d' % (self.curr_fofindex+1-self.start_fofindex,numtot))
        for coadd_mb_obs_lists,mb_obs_lists in self.imageio:
            numfof += 1

            foflen = len(mb_obs_lists)
            print('    num in fof: %d' % foflen)

            # get data to fill
            self.curr_data = self._make_struct(num=foflen)
            for i in xrange(foflen):
                self._set_default_data_for_fofind(i)
                
            #####################################################################
            # fit the fof once with no nbrs
            # sort by stamp size
            # set weight to uberseg if more than one thing in fof
            for coadd_mb_obs_list,mb_obs_list in zip(coadd_mb_obs_lists,mb_obs_lists):
                for obs_list in mb_obs_list:
                    for obs in obs_list:
                        if obs.meta['flags'] == 0:
                            if foflen > 1:
                                obs.weight = getattr(obs,'weight_us',obs.weight)
                            else:
                                obs.weight = getattr(obs,'weight_raw',obs.weight)
                            obs.weight_orig = obs.weight.copy()
                for obs_list in coadd_mb_obs_list:
                    for obs in obs_list:
                        if obs.meta['flags'] == 0:
                            if foflen > 1:
                                obs.weight = getattr(obs,'weight_us',obs.weight)
                            else:
                                obs.weight = getattr(obs,'weight_raw',obs.weight)
                            obs.weight_orig = obs.weight.copy()

            bs = []
            for coadd_mb_obs_list,mb_obs_list in zip(coadd_mb_obs_lists,mb_obs_lists):
                box_size = self._get_box_size(mb_obs_list)
                if box_size < 0:
                    box_size = self._get_box_size(coadd_mb_obs_list)
                bs.append(box_size)
            bs = numpy.array(bs)
            q = numpy.argsort(bs)
            q = q[::-1] # sort to fit biggest to smallest
            for i in q:
                self.curr_data_index = i
                coadd_mb_obs_list = coadd_mb_obs_lists[i]
                mb_obs_list = mb_obs_lists[i]
                if foflen > 1:
                    print('  fof obj: %d:%d' % (self.curr_data_index+1,foflen))
                print('    id: %d' % mb_obs_list.meta['id'])

                num += 1
                ti = time.time()
                self.fit_obj(coadd_mb_obs_list,mb_obs_list,nbrs_fit_data=None,make_epoch_data=True)
                ti = time.time()-ti
                print('    time: %f' % ti)


            #####################################################################
            # now fit again with nbrs if needed
            if foflen > 1:

                if self['mof']['write_convergence_data']:
                    self._write_convergence_data(mb_obs_lists,self.curr_data, \
                                                 self['mof']['convergence_model'],init=True)

                converged = False
                for itr in xrange(self['mof']['max_itr']):
                    print('itr %d - fof index %d:%d ' % (itr+1,\
                                                         self.curr_fofindex+1-self.start_fofindex,\
                                                         numtot))

                    # switch back to non-uberseg weights
                    if itr >= self['mof']['min_useg_itr']:
                        for coadd_mb_obs_list,mb_obs_list in zip(coadd_mb_obs_lists,mb_obs_lists):
                            for obs_list in mb_obs_list:
                                for obs in obs_list:
                                    if obs.meta['flags'] == 0:
                                        obs.weight = getattr(obs,'weight_raw',obs.weight)
                                        obs.weight_orig = obs.weight.copy()
                            for obs_list in coadd_mb_obs_list:
                                for obs in obs_list:
                                    if obs.meta['flags'] == 0:
                                        obs.weight = getattr(obs,'weight_raw',obs.weight)
                                        obs.weight_orig = obs.weight.copy()

                    # data
                    self.prev_data = self.curr_data.copy()
                    
                    # fitting
                    for i in numpy.random.choice(foflen,size=foflen,replace=False):
                        self.curr_data_index = i
                        
                        coadd_mb_obs_list = coadd_mb_obs_lists[i]
                        mb_obs_list = mb_obs_lists[i]
                        print('  fof obj: %d:%d - itr %d' % (self.curr_data_index+1,foflen,itr+1))
                        print('    id: %d' % mb_obs_list.meta['id'])

                        num += 1
                        ti = time.time()
                        self.fit_obj(coadd_mb_obs_list,mb_obs_list,
                                     nbrs_fit_data=self.curr_data,
                                     make_epoch_data=False,
                                     make_nbrs_data=True if itr == 0 else False)
                        ti = time.time()-ti
                        print('    time: %f' % ti)

                    if self['mof']['write_convergence_data']:
                        self._write_convergence_data(mb_obs_lists,self.curr_data, \
                                                     self['mof']['convergence_model'],init=False)

                    print('  convergence itr %d:' % (itr+1))
                    if self._check_convergence(foflen,itr,coadd_mb_obs_lists,mb_obs_lists) and itr >= self['mof']['min_itr']:
                        converged = True
                        break

                print('  convergence fof index: %d' % (self.curr_fofindex+1-self.start_fofindex))
                print('    converged: %s' % str(converged))
                print('    num itr: %d' % (itr+1))
            else:
                # one object in fof, so set mof flags
                models_to_check,pars_models_to_check,cov_models_to_check,npars = self._get_models_to_check()
                for model in models_to_check:
                    n = Namer(model)
                    self.curr_data[n('mof_flags')] = 0

            # append data and incr.
            self.data.extend(list(self.curr_data))
            self.curr_fofindex += 1

            tm=time.time()-t0
            self._try_checkpoint(tm)

            if self.curr_fofindex-self.start_fofindex < numtot:
                print('fof index: %d:%d' % (self.curr_fofindex+1-self.start_fofindex,numtot))

        tm=time.time()-t0
        print("time: %f" % tm)
        if num > 0:
            print("time per fit: %f" % (tm/num))
        if numfof > 0:
            print("time per fof: %f" % (tm/numfof))

        self.done = True

    def _write_convergence_data(self,mb_obs_lists,curr_data,model,init=False):
        for i in xrange(len(mb_obs_lists)):
            iter_fname = 'iter_pars_%d.dat' % (mb_obs_lists[i].meta['id'])
            if init:
                os.system('rm -f %s' % iter_fname)
                os.system('touch %s' % iter_fname)
            with open(iter_fname,'a+') as fp:
                if init:
                    fp.write('#row col g1 g2 T')
                    for j in xrange(self['nband']):
                        fp.write(' b%d' % j)
                    fp.write(' row_err col_err g1_err g2_err T_err')
                    for j in xrange(self['nband']):
                        fp.write(' b%d_err' % j)
                    fp.write('\n')

                if self['fit_coadd_galaxy']:
                    n = Namer('coadd_%s' % model)
                else:
                    n = Namer(model)
                npars = len(curr_data[n('max_pars')][i])

                for j in xrange(npars):
                    fp.write('%0.20g ' % curr_data[n('max_pars')][i,j])
                for j in xrange(npars):
                    fp.write('%0.20g ' % numpy.sqrt(curr_data[n('max_pars_cov')][i,j,j]))
                fp.write('\n')

    def _get_dtype(self):
        dt = super(MOFNGMixer,self)._get_dtype()
        models_to_check,pars_models_to_check,cov_models_to_check,npars = self._get_models_to_check()
        
        for model in models_to_check:
            n = Namer(model)
            dt += [(n('mof_flags'),'i4'),
                   (n('mof_num_itr'),'i4'),
                   (n('mof_abs_diff'),'f8',npars),
                   (n('mof_frac_diff'),'f8',npars),
                   (n('mof_err_diff'),'f8',npars)]
            
        dt += [('fofind','i8')]
            
        return dt

    def _make_struct(self,num=1):
        """
        make an output structure
        """
        data = super(MOFNGMixer,self)._make_struct(num=num)
        models_to_check,pars_models_to_check,cov_models_to_check,npars = self._get_models_to_check()
        for model in models_to_check:
            n = Namer(model)
            data[n('mof_flags')] = NO_ATTEMPT
            data[n('mof_num_itr')] = DEFVAL
            data[n('mof_abs_diff')] = DEFVAL
            data[n('mof_frac_diff')] = DEFVAL
            data[n('mof_err_diff')] = DEFVAL

        data['fofind'] = DEFVAL

        return data
