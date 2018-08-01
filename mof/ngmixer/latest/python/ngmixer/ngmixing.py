#!/usr/bin/env python
from __future__ import print_function
import numpy
import time
import pprint
import os
import fitsio

# local imports
from . import imageio
from . import fitting
from . import files
from .defaults import object_blacklist
from .defaults import DEFVAL,_CHECKPOINTS_DEFAULT_MINUTES
from .defaults import NO_ATTEMPT,NO_CUTOUTS,BOX_SIZE_TOO_BIG,IMAGE_FLAGS,BAD_OBJ,UTTER_FAILURE,OBJECT_IN_BLACKLIST

from .util import UtterFailure, seed_numpy


class NGMixer(dict):
    def __init__(self,
                 config_file,
                 data_files,
                 work_dir='.',
                 output_file=None,
                 fof_range=None,
                 fof_file=None,
                 mof_file=None,
                 models_file=None,
                 extra_data={},
                 random_seed=None,
                 init_only=False,
                 profile=False,
                 make_plots=False,
                 verbosity=0,
                 config=None):

        # parameters
        if config is not None:
            self.update(config)
        else:
            self.update(files.read_yaml(config_file))
        self['config_file'] = config_file
        self['work_dir'] = work_dir
        self['make_plots'] = self.get('make_plots',make_plots)
        self['fit_coadd_galaxy'] = self.get('fit_coadd_galaxy',False)
        self['fit_me_galaxy'] = self.get('fit_me_galaxy',True)
        self['max_box_size']=self.get('max_box_size',2048)
        self['verbosity'] = verbosity

        self.profile = profile
        self.extra_data = extra_data
        
        # random numbers
        seed_numpy(random_seed)

        self._set_defaults()
        self.fof_range=fof_range
        self._set_imageio(data_files, fof_range, fof_file, mof_file, extra_data)


        self._set_priors()
        self._set_fitter_and_data()

        pprint.pprint(self)

        # checkpointing and outputs
        self.output_file = output_file

        # this may be over-ridden if we restarted from a checkpoint
        self.start_fofindex = 0
        self._setup_checkpoints()
        self._setup_output_data()


        # run the code
        if not init_only:
            if self.profile:
                self.go_profile()
            else:
                self.go()

    def go(self):
        self.do_fits()
        self.write_data()
        if self.done:
            self.cleanup_checkpoint()

    def go_profile(self):
        import cProfile
        import pstats

        print("doing profile")

        cProfile.runctx('self.go()',
                        globals(),locals(),
                        'profile_stats')
        p = pstats.Stats('profile_stats')
        p.sort_stats('time').print_stats()


    def _set_defaults(self):
        if self['verbosity'] > 0:
            # over-ride anything in the config
            self['verbose'] = True
        else:
            self['verbose'] = self.get('verbose',False)

    def _set_imageio(self, data_files, fof_range, fof_file, mof_file, extra_data):
        """
        determine and instantiate the imageio class.  Set some file
        related attributes
        """
        # read the data
        io_type = self['imageio']['io_type']
        imageio_class = imageio.get_imageio_class(io_type)
        self.imageio = imageio_class(self,
                                     data_files,
                                     fof_range=fof_range,
                                     fof_file=fof_file,
                                     mof_file=mof_file,
                                     extra_data=extra_data)
        self.curr_fofindex = 0

        self['nband'] = self.imageio.get_num_bands()
        self.iband = range(self['nband'])

    def _set_fitter_and_data(self):
        # get the fitter
        fitter_class = fitting.get_fitter_class(self['fitter_type'])
        self.fitter = fitter_class(self)

        def_data = self.fitter.get_default_fit_data(self['fit_me_galaxy'],
                                                    self['fit_coadd_galaxy'])

        def_edata = self.fitter.get_default_epoch_fit_data()

        def_ndata = self.fitter.get_default_nbrs_data()

        self.default_data = def_data
        self.default_epoch_data = def_edata
        self.default_nbrs_data = def_ndata

    def _set_priors(self):
        """
        Set priors on the parameters we will fit
        """
        from .priors import set_priors
        set_priors(self)

    def get_data(self):
        return numpy.array(self.data,dtype=self.data_dtype)

    def get_epoch_data(self):
        return numpy.array(self.epoch_data,dtype=self.epoch_data_dtype)

    def get_nbrs_data(self):
        return numpy.array(self.nbrs_data,dtype=self.nbrs_data_dtype)

    def get_file_meta_data(self):
        return self.imageio.get_file_meta_data()

    def _extract_nbrs_data(self,coadd_mb_obs_lists,mb_obs_lists):
        if 'mof_fit_data' not in self.extra_data:
            raise ValueError('MOF fit data must be given to extract nbrs_fit_data!')

        cids = []
        for mb_obs_list in mb_obs_lists:
            cids.append(mb_obs_list.meta['id'])
        
        nbrs_fit_data = []
        for cid in cids:
            q, = numpy.where(self.extra_data['mof_fit_data']['id'] == cid)
            if len(q) != 1:
                raise ValueError('MOF data for object %d not found!' % cid)
            nbrs_fit_data.append(self.extra_data['mof_fit_data'][q[0]])
        nbrs_fit_data = numpy.array(nbrs_fit_data,dtype=self.extra_data['mof_fit_data'].dtype.descr)

        return nbrs_fit_data

    def do_fits(self):
        """
        Fit all objects in our list
        """

        self.done = False

        print('doing fits')

        t0=time.time()
        num = 0
        numtot = self.imageio.get_num_fofs()

        print('fof index: %d:%d' % (self.curr_fofindex+1-self.start_fofindex,numtot))
        for coadd_mb_obs_lists,mb_obs_lists in self.imageio:
            foflen = len(mb_obs_lists)

            # get data to fill
            self.curr_data = self._make_struct(num=foflen)
            for tag in self.default_data.dtype.names:
                self.curr_data[tag][:] = self.default_data[tag]
            self.curr_data_index = 0

            if 'mof_fit_data' in self.extra_data:
                nbrs_fit_data = self._extract_nbrs_data(coadd_mb_obs_lists,mb_obs_lists)
                nbrs_meta_data = self.extra_data['mof_nbrs_data']
            else:
                nbrs_fit_data = None
                nbrs_meta_data = None

            # fit the fof
            for coadd_mb_obs_list,mb_obs_list in zip(coadd_mb_obs_lists,mb_obs_lists):
                if foflen > 1:
                    print('fof obj: %d:%d' % (num,foflen))
                print('    id: %d' % mb_obs_list.meta['id'])

                num += 1
                ti = time.time()
                self.fit_obj(
                    coadd_mb_obs_list,
                    mb_obs_list,
                    nbrs_fit_data=nbrs_fit_data,
                    nbrs_meta_data=nbrs_meta_data,
                )
                ti = time.time()-ti
                print('    time: %f' % ti)

                self.curr_data_index += 1

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
            print("time per: %f" % (tm/num))

        self.done = True

    def _check_basic_things(self, coadd_mb_obs_list, mb_obs_list):

        if mb_obs_list.meta['id'] in object_blacklist:
            print("    skipping bad object:",mb_obs_list.meta['id'])
            return OBJECT_IN_BLACKLIST
        # get the box size
        for obsl in [coadd_mb_obs_list, mb_obs_list]:
            box_size = self._get_box_size(obsl)
            if box_size > 0:
                break
        self.curr_data['box_size'][self.curr_data_index] = box_size
        print('    box_size: %d' % self.curr_data['box_size'][self.curr_data_index])

        # check flags
        flags = 0

        if self['fit_coadd_galaxy']:
            if coadd_mb_obs_list.meta['obj_flags'] != 0:
                flags |= BAD_OBJ
                print('    skipping bad object')
            flags |= self._obj_check(coadd_mb_obs_list)

        if self['fit_me_galaxy']:
            if mb_obs_list.meta['obj_flags'] != 0:
                flags |= BAD_OBJ
                print('    skipping bad object')
            flags |= self._obj_check(mb_obs_list)

        return flags

    def _get_box_size(self, mb_obs_list):
        box_size = DEFVAL
        for band,obs_list in enumerate(mb_obs_list):
            for obs in obs_list:
                if obs.meta['flags'] == 0:
                    box_size = obs.image.shape[0]
                    break
        return box_size

    def _obj_check(self, mb_obs_list):
        """
        Check box sizes, number of cutouts
        Require good in all bands
        """
        for band,obs_list in enumerate(mb_obs_list):
            flags=self._obj_check_one(band,obs_list)
            if flags != 0:
                break
        return flags

    def _obj_check_one(self, band, obs_list):
        """
        Check box sizes, number of cutouts, flags on images
        """
        flags=0

        ncutout = len(obs_list)

        if ncutout == 0:
            print('    no cutouts')
            flags |= NO_CUTOUTS
            return flags

        num_use = 0
        for obs in obs_list:
            if obs.meta['flags'] == 0:
                num_use += 1

        if num_use < ncutout:
            print("    for band %d removed %d/%d images due to flags"
                     % (band, ncutout-num_use, ncutout))

        if num_use == 0:
            flags |= IMAGE_FLAGS
            return flags

        box_size = self.curr_data['box_size'][self.curr_data_index]
        if box_size > self['max_box_size']:
            print('    box size too big: %d' % box_size)
            flags |= BOX_SIZE_TOO_BIG

        return flags

    def fit_obj(self,coadd_mb_obs_list,mb_obs_list,
                nbrs_fit_data=None,nbrs_meta_data=None,
                make_epoch_data=True,make_nbrs_data=False):
        """
        fit a single object
        """

        t0 = time.time()

        #check flags
        flags = self._check_basic_things(coadd_mb_obs_list,mb_obs_list)

        if flags == 0:
            fit_flags = self.fit_all_obs_lists(coadd_mb_obs_list,mb_obs_list,
                                               nbrs_fit_data=nbrs_fit_data,
                                               nbrs_meta_data=nbrs_meta_data,
                                               make_epoch_data=make_epoch_data,
                                               make_nbrs_data=make_nbrs_data)
            flags |= fit_flags

        # add in data
        self.curr_data['flags'][self.curr_data_index] = flags
        self.curr_data['time_last_fit'][self.curr_data_index] = time.time()-t0
        self.curr_data['obj_flags'][self.curr_data_index] = mb_obs_list.meta['obj_flags']

        # fill in from mb_obs_meta
        for tag in mb_obs_list.meta['meta_data'].dtype.names:
            self.curr_data[tag][self.curr_data_index] = mb_obs_list.meta['meta_data'][tag][0]

    def _fill_epoch_data(self,mb_obs_list):
        # fill in epoch data
        for band,obs_list in enumerate(mb_obs_list):
            for obs in obs_list:
                if 'fit_data' in obs.meta and obs.meta['fit_data'] is not None \
                   and 'meta_data' in obs.meta and obs.meta['flags'] == 0:

                    ed = self._make_epoch_struct()

                    for tag in self.default_epoch_data.dtype.names:
                        ed[tag] = self.default_epoch_data[tag]

                    for tag in obs.meta['fit_data'].dtype.names:
                        ed[tag] = obs.meta['fit_data'][tag][0]

                    for tag in obs.meta['meta_data'].dtype.names:
                        ed[tag] = obs.meta['meta_data'][tag][0]

                    self.epoch_data.extend(list(ed))

    def _fill_nbrs_data(self,mb_obs_list):
        # fill in nbrs data
        for band,obs_list in enumerate(mb_obs_list):
            for obs in obs_list:
                if 'nbrs_data' in obs.meta and obs.meta['nbrs_data'] is not None \
                        and 'fit_data' in obs.meta and obs.meta['fit_data'] is not None \
                        and 'meta_data' in obs.meta and obs.meta['flags'] == 0:

                    ed = self._make_nbrs_struct(len(obs.meta['nbrs_data']))

                    for tag in self.default_nbrs_data.dtype.names:
                        ed[tag] = self.default_nbrs_data[tag][0]

                    for tag in obs.meta['nbrs_data'].dtype.names:
                        ed[tag] = obs.meta['nbrs_data'][tag]

                    for i in xrange(len(obs.meta['nbrs_data'])):
                        for tag in obs.meta['meta_data'].dtype.names:
                            ed[tag][i] = obs.meta['meta_data'][tag][0]

                    self.nbrs_data.extend(list(ed))

    def fit_all_obs_lists(self,coadd_mb_obs_list,mb_obs_list,
                          nbrs_fit_data=None,nbrs_meta_data=None,
                          make_epoch_data=True,make_nbrs_data=True):
        """
        fit all obs lists
        """

        fit_flags = None

        if self['fit_me_galaxy']:
            nfit=sum([len(ol) for ol in mb_obs_list])
            print('    fitting',nfit,'me galaxy')
            try:
                me_fit_flags = self.fitter(mb_obs_list,coadd=False,
                                           nbrs_fit_data=nbrs_fit_data,
                                           nbrs_meta_data=nbrs_meta_data)

                # fill in epoch data
                if make_epoch_data:
                    self._fill_epoch_data(mb_obs_list)

                if make_nbrs_data and self['model_nbrs']:
                    self._fill_nbrs_data(mb_obs_list)

                # fill in fit data
                for tag in mb_obs_list.meta['fit_data'].dtype.names:
                    self.curr_data[tag][self.curr_data_index] = mb_obs_list.meta['fit_data'][tag][0]

            except UtterFailure as err:
                print("    me fit got utter failure error: %s" % str(err))
                me_fit_flags = UTTER_FAILURE

            if fit_flags is None:
                fit_flags = 0
            fit_flags |= me_fit_flags

        if self['fit_coadd_galaxy']:
            print('    fitting coadd galaxy')
            try:
                coadd_fit_flags = self.fitter(coadd_mb_obs_list,coadd=True,
                                              nbrs_fit_data=nbrs_fit_data,
                                              nbrs_meta_data=nbrs_meta_data)

                # fill in epoch data
                if make_epoch_data:
                    self._fill_epoch_data(coadd_mb_obs_list)

                if make_nbrs_data and self['model_nbrs']:
                    self._fill_nbrs_data(coadd_mb_obs_list)

                # fill in fit data
                for tag in coadd_mb_obs_list.meta['fit_data'].dtype.names:
                    self.curr_data[tag][self.curr_data_index] = \
                            coadd_mb_obs_list.meta['fit_data'][tag][0]

            except UtterFailure as err:
                print("    coadd fit got utter failure error: %s" % str(err))
                coadd_fit_flags = UTTER_FAILURE

            if fit_flags is None:
                fit_flags = 0
            fit_flags |= coadd_fit_flags

        if fit_flags is None:
            fit_flags = NO_ATTEMPT

        return fit_flags

    def _setup_output_data(self):
        """
        set the output data structures, if not already set based
        on the input checkpoint data
        """
        if self.checkpoint_data is None:
            self.data = []
            self.data_dtype = self._get_dtype()
            self.epoch_data = []
            self.epoch_data_dtype = self._get_epoch_dtype()
            self.nbrs_data = []
            self.nbrs_data_dtype = self._get_nbrs_dtype()

    def _get_epoch_dtype(self):
        """
        makes epoch dtype
        """
        dt = self.imageio.get_epoch_meta_data_dtype()
        dt += self.fitter.get_epoch_fit_data_dtype()
        return dt

    def _make_epoch_struct(self,ncutout=1):
        """
        returns ncutout epoch structs to be filled
        """
        dt = self._get_epoch_dtype()
        epoch_data = numpy.zeros(ncutout, dtype=dt)
        return epoch_data

    def _get_nbrs_dtype(self):
        """
        make nbrs dtype
        """
        dt = self.imageio.get_epoch_meta_data_dtype()
        dt += self.fitter.get_nbrs_data_dtype()
        return dt

    def _make_nbrs_struct(self,ncutout=1):
        """
        returns ncutout nbrs structs to be filled
        """
        dt = self._get_nbrs_dtype()
        nbrs_data = numpy.zeros(ncutout, dtype=dt)
        return nbrs_data

    def _get_dtype(self):
        dt = self.imageio.get_meta_data_dtype()
        dt += [('flags','i4'),
               ('time_last_fit','f8'),
               ('box_size','i2'),
               ('obj_flags','i4')]
        dt += self.fitter.get_fit_data_dtype(self['fit_me_galaxy'],self['fit_coadd_galaxy'])
        return dt

    def _make_struct(self,num=1):
        """
        make an output structure
        """
        dt = self._get_dtype()
        data=numpy.zeros(num, dtype=dt)
        data['flags'] = NO_ATTEMPT
        data['time_last_fit'] = DEFVAL
        data['box_size'] = DEFVAL
        data['obj_flags'] = NO_ATTEMPT
        return data

    def _setup_checkpoints(self):
        """
        Set up the checkpoint times in minutes and data

        self.checkpoint_data and self.checkpoint_file
        """
        self.checkpoints = self.get('checkpoints',_CHECKPOINTS_DEFAULT_MINUTES)
        self.n_checkpoint    = len(self.checkpoints)
        self.checkpointed    = [0]*self.n_checkpoint

        self._set_checkpoint_data()

        if self.checkpoint_file is not None:
            self.do_checkpoint=True
        else:
            self.do_checkpoint=False

    def _set_checkpoint_data(self):
        """
        See if checkpoint data was sent
        """
        import fitsio
        import cPickle

        self.checkpoint_data = None

        if self.output_file is not None:
            if self.output_file[-5:] == '.fits':
                self.checkpoint_file = self.output_file.replace('.fits','-checkpoint.fits')
            else:
                self.checkpoint_file = self.output_file.replace('.fit','-checkpoint.fits')
            if os.path.exists(self.checkpoint_file):
                self.checkpoint_data={}
                print('reading checkpoint data: %s' % self.checkpoint_file)
                with fitsio.FITS(self.checkpoint_file) as fobj:
                    self.checkpoint_data['data'] = fobj['model_fits'][:]

                    if 'epoch_data' in fobj:
                        self.checkpoint_data['epoch_data']=fobj['epoch_data'][:]

                    if 'nbrs_data' in fobj:
                        self.checkpoint_data['nbrs_data']=fobj['nbrs_data'][:]

                    if 'checkpoint_data' in fobj:
                        self.checkpoint_data['checkpoint_data'] = fobj['checkpoint_data'][:]

        if self.checkpoint_data is not None:
            # data
            self.data = self.checkpoint_data['data']
            #fitsio.fitslib.array_to_native(self.data, inplace=True)
            self.data = self.data.byteswap().newbyteorder()
            self.data_dtype = self._get_dtype()

            # for nband==1 the written array drops the arrayness
            if self['nband']==1:
                self.data.dtype = self.data_dtype

            self.data = list(self.data)

            # epoch data
            if 'epoch_data' in self.checkpoint_data:
                self.epoch_data = self.checkpoint_data['epoch_data']
                self.epoch_data = self.epoch_data.byteswap().newbyteorder()
                self.epoch_data_dtype = self._get_epoch_dtype()
                self.epoch_data = list(self.epoch_data)
            else:
                self.epoch_data_dtype = self._get_epoch_dtype()
                self.epoch_data = []

            # nbrs data
            if 'nbrs_data' in self.checkpoint_data:
                self.nbrs_data = self.checkpoint_data['nbrs_data']
                self.nbrs_data = self.nbrs_data.byteswap().newbyteorder()
                self.nbrs_data_dtype = self._get_nbrs_dtype()
                self.nbrs_data = list(self.nbrs_data)
            else:
                self.nbrs_data_dtype = self._get_nbrs_dtype()
                self.nbrs_data = []

            # checkpoint data
            rs = cPickle.loads(self.checkpoint_data['checkpoint_data']['random_state'][0])
            numpy.random.set_state(rs)
            self.curr_fofindex = self.checkpoint_data['checkpoint_data']['curr_fofindex'][0]
            self.imageio.set_fof_start(self.curr_fofindex)
            self.start_fofindex = self.checkpoint_data['checkpoint_data']['curr_fofindex'][0]

    def _try_checkpoint(self, tm):
        """
        Checkpoint at certain intervals.
        Potentially modified self.checkpointed
        """

        should_checkpoint, icheck = self._should_checkpoint(tm)

        if should_checkpoint:
            self._write_checkpoint(tm)
            self.checkpointed[icheck]=1

    def _should_checkpoint(self, tm):
        """
        Should we write a checkpoint file?
        """

        should_checkpoint=False
        icheck=-1

        if self.do_checkpoint:
            tm_minutes=tm/60

            for i in xrange(self.n_checkpoint):

                checkpoint=self.checkpoints[i]
                checkpointed=self.checkpointed[i]

                if tm_minutes > checkpoint and not checkpointed:
                    should_checkpoint=True
                    icheck=i

        return should_checkpoint, icheck

    def _write_checkpoint(self, tm):
        """
        Write out the current data structure to a temporary
        checkpoint file.
        """
        import fitsio
        from .files import StagedOutFile
        import cPickle

        print('checkpointing at %f minutes' % (tm/60))
        print(self.checkpoint_file)

        # make checkpoint data
        cd = numpy.zeros(1,dtype=[('curr_fofindex','i8'),('random_state','|S16384')])
        cd['curr_fofindex'][0] = self.curr_fofindex
        cd['random_state'][0] = cPickle.dumps(numpy.random.get_state())

        with StagedOutFile(self.checkpoint_file, tmpdir=self['work_dir']) as sf:
            with fitsio.FITS(sf.path,'rw',clobber=True) as fobj:

                fobj.write(numpy.array(self.data,dtype=self.data_dtype), extname="model_fits")

                if len(self.epoch_data) > 0:
                    fobj.write(numpy.array(self.epoch_data,dtype=self.epoch_data_dtype), extname="epoch_data")

                if len(self.nbrs_data) > 0:
                    fobj.write(numpy.array(self.nbrs_data,dtype=self.nbrs_data_dtype), extname="nbrs_data")

                fobj.write(cd, extname="checkpoint_data")

    def cleanup_checkpoint(self):
        """
        if we get this far, we have succeeded in writing the data. We can remove
        the checkpoint file
        """
        if os.path.exists(self.checkpoint_file):
            print('removing checkpoint file: %s' % self.checkpoint_file)
            os.remove(self.checkpoint_file)

    def write_data(self):
        """
        write the actual data.  clobber existing
        """
        if self.output_file is not None:
            self.meta = self.get_file_meta_data()

            from .githash import hash as ngmixer_hash
            try:
                from ngmix.githash import hash as ngmix_hash
            except:
                ngmix_hash = ' '
            self.githashes = numpy.zeros(1,dtype=[('ngmixer_githash','S%d' % len(ngmixer_hash)),
                                                  ('ngmix_githash','S%d' % len(ngmix_hash))])
            self.githashes['ngmixer_githash'][:] = ngmixer_hash
            self.githashes['ngmix_githash'][:] = ngmix_hash
                        
            from .files import StagedOutFile
            work_dir = self['work_dir']
            with StagedOutFile(self.output_file, tmpdir=work_dir) as sf:
                print('writing: %s' % sf.path)
                with fitsio.FITS(sf.path,'rw',clobber=True) as fobj:
                    fobj.write(self.get_data(),extname="model_fits")

                    if self.epoch_data is not None:
                        fobj.write(self.get_epoch_data(),extname="epoch_data")

                    if self.nbrs_data is not None and len(self.nbrs_data) > 0:
                        fobj.write(self.get_nbrs_data(),extname="nbrs_data")

                    if self.meta is not None:
                        fobj.write(self.meta,extname="meta_data")

                    if self.githashes is not None:
                        fobj.write(self.githashes,extname="githashes")


