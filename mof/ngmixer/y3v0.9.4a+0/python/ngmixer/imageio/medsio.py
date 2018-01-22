#!/usr/bin/env python
from __future__ import print_function
import os
import numpy
import copy
import fitsio

# meds and ngmix imports
import meds
from ngmix import Jacobian
from ngmix import Observation, ObsList, MultiBandObsList

# local imports
from .imageio import ImageIO
from .extractor_corrector import MEDSExtractorCorrector
from ..defaults import DEFVAL,IMAGE_FLAGS
from .. import nbrsfofs

class MEDSImageIO(ImageIO):
    """
    Class for MEDS image I/O.

    The following methods must be defined in a subclass (see the docstrings for the
    precise function signature).

        MEDSImageIO._load_psf_data() - load the psf files if needed

        MEDSImageIO._get_psf_image(band,mindex,icut) - get PSF image for band, index
            in MEDS file (mindex) and cutout index (icut)


    """

    def __init__(self,*args,**kwargs):
        super(MEDSImageIO, self).__init__(*args,**kwargs)
        self.conf = args[0]
        self._set_defaults()
        if self.conf['reject_outliers']:
            print("will reject outliers")

        # deal with extra data
        # leaving this just in case we need it
        if self.extra_data is None:
            self.extra_data = {}

        meds_files = args[1]
        if not isinstance(meds_files, list):
            self.meds_files = [meds_files]
        else:
            self.meds_files =  meds_files

        # do sub range files if needed
        self._setup_work_files()

        # load meds files and image flags array
        self._load_meds_files()

        # set extra config
        self.iband = range(len(self.meds_list))
        self.conf['nband'] = len(self.meds_list)

        self.conf['max_cutouts'] = self.conf.get('max_cutouts',None)
        self.conf['psfs_in_file'] = self.conf.get('psfs_in_file',False)

        # indexing of fofs
        self._set_and_check_index_lookups()

        # psfs
        if not self.conf['psfs_in_file']:
            self._load_psf_data()

        # make sure if we are doing nbrs we have the info we need
        if self.conf['model_nbrs']:
            assert 'nbrs' in self.extra_data,"You must supply a nbrs file to model nbrs!"

    def _set_defaults(self):
        self.conf['min_weight'] = self.conf.get('min_weight',-numpy.inf)
        self.conf['reject_outliers'] = self.conf.get('reject_outliers',True) # from cutouts
        self.conf['model_nbrs'] = self.conf.get('model_nbrs',False)
        self.conf['ignore_zero_images'] = self.conf.get('ignore_zero_images',False)


        # max fraction of image masked in bitmask or that has zero weight
        self.conf['max_bmask_frac'] = self.conf.get('max_bmask_frac',1.0)
        self.conf['max_zero_weight_frac'] = self.conf.get('max_zero_weight_frac',1.0)
        self.conf['symmetrize_weight'] = self.conf.get('symmetrize_weight',False)


        # check this region around the center
        self.conf['central_bmask_radius'] = \
                self.conf.get('central_bmask_radius',None)

    def _load_psf_data(self):
        pass

    def _get_psf_image(self, band, mindex, icut):
        """
        Get an image representing the psf

        should return

        image, center of psf, estimate of width in pixels, filename

        The filename is optional and is just for debugging purposes.
        """

        if not self.conf['psfs_in_file']:
            raise NotImplementedError("only use base class method when "
                                      "psfs are in the MEDS file")

        meds=self.meds_list[band]
        psfim = meds.get_psf(mindex, icut)

        cat = meds.get_cat()
        names=cat.dtype.names
        if 'psf_sigma' in names:
            sigma_pix = cat['psf_sigma'][mindex, icut]
        else:
            sigma_pix = 2.5

        if 'psf_cutout_row' in names:
            cen = (
                cat['psf_cutout_row'][mindex, icut],
                cat['psf_cutout_col'][mindex, icut],
            )
        else:
            cen = (numpy.array(psfim.shape)-1.0)/2.0

        return psfim, cen, sigma_pix, 'none'

    def _get_sub_fname(self,fname):
        rng_string = '%s-%s' % (self.fof_range[0], self.fof_range[1])
        bname = os.path.basename(fname)
        bname = bname.replace('.fits.fz','').replace('.fits','')
        bname = '%s-%s.fits' % (bname, rng_string)
        newf = os.path.expandvars(os.path.join(self.conf['work_dir'], bname))
        return newf

    def _get_sub(self):
        """
        Local files will get cleaned up
        """
        extracted=[]


        corrmeds=self.conf.get('correct_meds',None)

        if corrmeds is not None:
            if self.mof_file is None:
                raise ValueError(
                    "you must send mof outputs to make "
                    "corrected MEDS files"
                )

            min_weight=self.conf.get('min_weight',0.0)

            band_not_in_name=self.conf.get("band_not_in_name",False)

            for iband,meds_file in enumerate(self.meds_files):
                print(meds_file)
                newf = self._get_sub_fname(meds_file)


                # we must assume the bands line up with the MOF run
                if band_not_in_name:
                    band=iband
                else:
                    band=None

                ex=MEDSExtractorCorrector(
                    self.mof_file,
                    meds_file,
                    self.fof_range[0],
                    self.fof_range[1],
                    newf,
                    replace_bad=corrmeds['replace_bad'],
                    reject_outliers=corrmeds['reject_outliers'],
                    min_weight=min_weight,
                    cleanup=True,
                    verbose=False,
                    make_plots=self.conf['make_plots'],
                    band=band,
                )
                extracted.append(ex)
            extracted.append(None)

        else:
            if self.fof_file is None:
                for meds_file in self.meds_files:
                    print(meds_file)
                    newf = self._get_sub_fname(meds_file)
                    ex=meds.MEDSExtractor(
                        meds_file,
                        self.fof_range[0],
                        self.fof_range[1],
                        newf,
                        copy_all=True,
                        cleanup=True,
                    )
                    extracted.append(ex)
                extracted.append(None)
            else:
                # do the fofs first
                print(self.fof_file)
                newf = self._get_sub_fname(self.fof_file)
                fofex = nbrsfofs.NbrsFoFExtractor(
                    self.fof_file,
                    self.fof_range[0],
                    self.fof_range[1],
                    newf,
                    cleanup=True,
                )

                # now do the meds
                for meds_file in self.meds_files:
                    print(meds_file)
                    newf = self._get_sub_fname(meds_file)
                    ex=meds.MEDSNumberExtractor(
                        meds_file,
                        fofex.numbers,
                        newf,
                        cleanup=True,
                    )
                    extracted.append(ex)
                extracted.append(fofex)

        return extracted

    def _setup_work_files(self):
        """
        Set up local, possibly sub-range meds files
        """
        self.meds_files_full = self.meds_files
        self.fof_file_full = self.fof_file
        self.extracted=None

        if 'correct_meds' in self.conf and self.fof_range is None:
            with meds.MEDS(self.meds_files_full[0]) as m:
                self.fof_range=[0,m.size-1]

        if self.fof_range is not None:
            extracted=self._get_sub()
            meds_files=[ex.sub_file for ex in extracted if ex is not None]
            if extracted[-1] is not None:
                self.fof_file = extracted[-1].sub_file
                meds_files = meds_files[:-1]
            self.meds_files = meds_files
            self.extracted = extracted

    def _set_and_check_index_lookups(self):
        """
        Deal with common indexing issues in one place

        Indexing notes:

        Each MEDS file consists of a list of objects to be fit, which are indexed by

           self.mindex = 0-offset from start of file

        The code runs by processing groups of objects, which we call FoFs. Note however
        that these groupings can be arbitrary. The FoFs are specified by a numpy array
        (read from a FITS table) which has two columns

            fofid - the ID of the fof group, 0-offset
            number - the ID number of the object in the coadd detection tile seg map

        We require that the FoF file number column matches the MEDS file, line-by-line.

        We do however build some lookup tables to enable easy translation. They are

            self.fofid2mindex = this is a dict keyed on the fofid - it returns the set of mindexes
                that give the members of the FoF group
            self.number2mindex = lookup table for converting numbers to mindexes
        These are both dictionaries which use python's hashing of keys. There may be a performance
        issue here, in terms of building the dicts, for large numbers of objects.
        """

        # warn the user
        print('making fof indexes')

        if self.fof_file is not None:
            read_fofs = True
            self.fof_data = fitsio.read(self.fof_file)
        else:
            read_fofs = False
            nobj = len(self.meds_list[0]['number'])
            self.fof_data = numpy.zeros(nobj,dtype=[('fofid','i8'),('number','i8')])
            self.fof_data['fofid'][:] = numpy.arange(nobj)
            self.fof_data['number'][:] = self.meds_list[0]['number'][:]

        # first, we error check
        for band,meds in enumerate(self.meds_list):
            msg = "FoF number is not the same as MEDS number for band %d!" % band
            assert numpy.array_equal(meds['number'],self.fof_data['number']),msg


        #set some useful stuff here
        self.fofids = numpy.unique(self.fof_data['fofid'])
        self.fofids = numpy.sort(self.fofids)
        self.num_fofs = len(self.fofids)

        if read_fofs:
            #build the fof hash
            self.fofid2mindex = {}
            for fofid in self.fofids:
                q, = numpy.where(self.fof_data['fofid'] == fofid)
                assert len(q) > 0, 'Found zero length FoF! fofid = %ld' % fofid
                assert numpy.array_equal(self.fof_data['number'][q],self.meds_list[0]['number'][q])
                self.fofid2mindex[fofid] = q.copy()
        else:
            #use a list of lists since is much faster to make
            self.fofid2mindex = []
            for fofid in self.fofids:
                self.fofid2mindex.append([fofid])
                assert self.fofid2mindex[fofid][0] == fofid

    def _flag_objects(coadd_mb_obs_lists,me_mb_obs_lists,mindexes):
        qnz, = numpy.where(self.extra_data['obj_flags']['flags'] != 0)
        for mindex,coadd_mb_obs_list,me_mb_obs_list in zip(mindexes,coadd_mb_obs_lists,me_mb_obs_lists):
            q, = numpy.where(self.extra_data['obj_flags']['id'][qnz] == me_mb_obs_list.meta['id'])
            if len(q) > 0:
                assert len(q) == 1
                assert me_mb_obs_list.meta['id'] ==  self.extra_data['obj_flags']['id'][qnz[q[0]]]
                coadd_mb_obs_list.meta['obj_flags'] |= self.extra_data['obj_flags']['flags'][qnz[q[0]]]
                me_mb_obs_list.meta['obj_flags'] |= self.extra_data['obj_flags']['flags'][qnz[q[0]]]

    def _add_nbrs_info(self,coadd_mb_obs_lists,me_mb_obs_lists,mindexes):
        """
        adds nbr info to obs lists
        """

        # save orig images and weights
        for mb_obs_list in coadd_mb_obs_lists:
            for obs_list in mb_obs_list:
                for obs in obs_list:
                    obs.image_orig = obs.image.copy()
                    obs.weight_orig = obs.weight.copy()

        for mb_obs_list in me_mb_obs_lists:
            for obs_list in mb_obs_list:
                for obs in obs_list:
                    obs.image_orig = obs.image.copy()
                    obs.weight_orig = obs.weight.copy()

        # do indexes
        for cen,mindex in enumerate(mindexes):
            nbrs_inds = []
            nbrs_ids = []

            # if len is 1, then only a single galaxy in the FoF and do nothing
            if len(mindexes) > 1:
                q, = numpy.where(self.extra_data['nbrs']['number'] == self.meds_list[0]['number'][mindex])
                for ind in q:
                    if self.extra_data['nbrs']['nbr_number'][ind] != -1:
                        qi, = numpy.where(self.meds_list[0]['number'][mindexes] == self.extra_data['nbrs']['nbr_number'][ind])
                        assert len(qi) == 1
                        nbrs_inds.append(qi[0])
                        nbrs_ids.append(self.meds_list[0]['id'][mindexes[qi[0]]])
                        assert coadd_mb_obs_lists[nbrs_inds[-1]].meta['id'] == nbrs_ids[-1]
                        assert me_mb_obs_lists[nbrs_inds[-1]].meta['id'] == nbrs_ids[-1]

                assert cen not in nbrs_inds,'weird error where cen_ind is in nbrs_ind!'

            coadd_mb_obs_lists[cen].update_meta_data({'nbrs_inds':nbrs_inds,'nbrs_ids':nbrs_ids,'cen_ind':cen})
            me_mb_obs_lists[cen].update_meta_data({'nbrs_inds':nbrs_inds,'nbrs_ids':nbrs_ids,'cen_ind':cen})

        # now do psfs and jacobians
        self._add_nbrs_psfs_and_jacs(coadd_mb_obs_lists,mindexes)
        self._add_nbrs_psfs_and_jacs(me_mb_obs_lists,mindexes)

    def _add_nbrs_psfs_and_jacs(self,mb_obs_lists,mindexes):
        # for each object
        for cen,mindex in enumerate(mindexes):
            # for each band per object
            for band,obs_list in enumerate(mb_obs_lists[cen]):
                # for each obs per band per object
                for obs in obs_list:
                    if obs.meta['flags'] == 0:
                        # for each nbr per obs per band per object
                        nbrs_psfs = []
                        nbrs_flags = []
                        nbrs_jacs = []
                        for ind in mb_obs_lists[cen].meta['nbrs_inds']:
                            psf_obs,jac = self._get_nbr_psf_obs_and_jac(band,cen,mindex,obs,ind,mindexes[ind],mb_obs_lists[ind])
                            nbrs_psfs.append(psf_obs)
                            nbrs_jacs.append(jac)

                            if psf_obs is None or jac is None or mb_obs_lists[ind].meta['obj_flags'] != 0:
                                nbrs_flags.append(1)
                            else:
                                nbrs_flags.append(0)

                        obs.update_meta_data({'nbrs_psfs':nbrs_psfs,'nbrs_flags':nbrs_flags,'nbrs_jacs':nbrs_jacs})

    def _get_nbr_psf_obs_and_jac(self,band,cen_ind,cen_mindex,cen_obs,nbr_ind,nbr_mindex,nbrs_obs_list):
        assert nbrs_obs_list.meta['id'] ==  self.meds_list[band]['id'][nbr_mindex]
        assert cen_obs.meta['id'] ==  self.meds_list[band]['id'][cen_mindex]

        cen_file_id = cen_obs.meta['meta_data']['file_id'][0]
        nbr_obs = None
        for obs in nbrs_obs_list[band]:
            if obs.meta['meta_data']['file_id'][0] == cen_file_id and self.meds_list[band]['id'][nbr_mindex] == obs.meta['id']:
                nbr_obs = obs

        if nbr_obs is not None:            
            # cause the object to be flagged above
            if nbr_obs.meta['flags'] != 0:
                # for debug
                if False:
                    assert False, "nbr obs has flags != 0 when cen does not! band = %d, cen id = %d, nbr_id = %d, file_id = %d" % \
                        (band,self.meds_list[band]['id'][cen_mindex],self.meds_list[band]['id'][nbr_mindex],cen_file_id)
                return None,None
            
            nbr_psf_obs = nbr_obs.get_psf()
            nbr_icut = nbr_obs.meta['icut']
            assert self.meds_list[band]['file_id'][nbr_mindex,nbr_icut] == cen_file_id

            # construct jacobian
            # get fiducial location of object in postage stamp
            row = self.meds_list[band]['orig_row'][nbr_mindex,nbr_icut] - cen_obs.meta['orig_start_row']
            col = self.meds_list[band]['orig_col'][nbr_mindex,nbr_icut] - cen_obs.meta['orig_start_col']
            nbr_jac = Jacobian(row=row,
                               col=col,
                               dudrow=self.meds_list[band]['dudrow'][nbr_mindex,nbr_icut],
                               dudcol=self.meds_list[band]['dudcol'][nbr_mindex,nbr_icut],
                               dvdrow=self.meds_list[band]['dvdrow'][nbr_mindex,nbr_icut],
                               dvdcol=self.meds_list[band]['dvdcol'][nbr_mindex,nbr_icut])
            # FIXME - the code below is wrong...I think - commented out for now
            #pixscale = jacob.get_scale()
            #row += pars_obj[0]/pixscale
            #col += pars_obj[1]/pixscale
            #nbr_jac.set_cen(row=row,col=col)

            return nbr_psf_obs,nbr_jac
        else:
            return self._get_offchip_nbr_psf_obs_and_jac(band,cen_ind,cen_mindex,cen_obs,nbr_ind,nbr_mindex,nbrs_obs_list)

    def _get_offchip_nbr_psf_obs_and_jac(self,band,cen_ind,cen_mindex,cen_obs,nbr_ind,nbr_mindex,nbrs_obs_list):
        assert False,'        FIXME: off-chip nbr %d for cen %d' % (nbr_ind+1,cen_ind+1)
        return None,None

    def get_num_fofs(self):
        return copy.copy(self.num_fofs - self.fof_start)

    def get_num_bands(self):
        """"
        returns number of bands for galaxy images
        """
        return copy.copy(self.conf['nband'])

    def get_meta_data_dtype(self):
        dt=[('id','i8'),
            ('number','i4'),
            ('ra','f8'),
            ('dec','f8'),
            ('nimage_tot','i4',(self.conf['nband'],)),
            ('fofid','i8')]
        return dt

    def _get_meta_row(self,num=1):
        # build the meta data
        dt = self.get_meta_data_dtype()
        meta_row = numpy.zeros(num,dtype=dt)
        for tag in meta_row.dtype.names:
            meta_row[tag][:] = DEFVAL
        return meta_row

    def get_epoch_meta_data_dtype(self):
        dt=[('id','i8'),     # could be coadd_objects_id
            ('number','i4'), # 1-n as in sextractor
            ('band_num','i2'),
            ('cutout_index','i4'), # this is the index in meds
            ('orig_row','f8'),
            ('orig_col','f8'),
            ('file_id','i4'),
            ('pixel_scale','f8')]   # id in meds file
        return dt

    def _get_epoch_meta_row(self,num=1):
        # build the meta data
        dt = self.get_epoch_meta_data_dtype()
        meta_row = numpy.zeros(num,dtype=dt)
        for tag in meta_row.dtype.names:
            meta_row[tag][:] = DEFVAL
        return meta_row


    def get_file_meta_data(self):
        meds_meta_list = self.meds_meta_list
        dt = meds_meta_list[0].dtype.descr

        if 'config_file' in self.conf:
            tmp,config_file = os.path.split(self.conf['config_file'])
            clen=len(config_file)
            dt += [('ngmixer_config','S%d' % clen)]

        flen=max([len(mf) for mf in self.meds_files_full] )
        dt += [('meds_file','S%d' % flen)]

        nband=len(self.meds_files_full)
        meta=numpy.zeros(nband, dtype=dt)

        for band in xrange(nband):
            meds_file = self.meds_files_full[band]
            meds_meta=meds_meta_list[band]
            mnames=meta.dtype.names
            for name in meds_meta.dtype.names:
                if name in mnames:
                    meta[name][band] = meds_meta[name][0]

            if 'config_file' in self.conf:
                meta['ngmixer_config'][band] = config_file
            meta['meds_file'][band] = meds_file

        return meta

    def __next__(self):
        if self.fofindex >= self.num_fofs:
            raise StopIteration
        else:
            fofid = self.fofids[self.fofindex]
            mindexes = self.fofid2mindex[fofid]
            coadd_mb_obs_lists = []
            me_mb_obs_lists = []
            for mindex in mindexes:
                print('  getting obj w/ id %d' % self.meds_list[0]['id'][mindex])
                
                c,me = self._get_multi_band_observations(mindex)
                
                # add fof ids here
                if self.fof_file is not None:
                    c.meta['meta_data']['fofid'][:] = fofid
                    me.meta['meta_data']['fofid'][:] = fofid
                
                coadd_mb_obs_lists.append(c)
                me_mb_obs_lists.append(me)

            if 'obj_flags' in self.extra_data:
                self._flag_objects(coadd_mb_obs_lists,me_mb_obs_lists,mindexes)

            if self.conf['model_nbrs']:
                self._add_nbrs_info(coadd_mb_obs_lists,me_mb_obs_lists,mindexes)

            self.fofindex += 1
            return coadd_mb_obs_lists,me_mb_obs_lists

    next = __next__

    def _get_multi_band_observations(self, mindex):
        """
        Get an ObsList object for the Coadd observations
        Get a MultiBandObsList object for the SE observations.
        """

        coadd_mb_obs_list=MultiBandObsList()
        mb_obs_list=MultiBandObsList()

        for band in self.iband:
            cobs_list, obs_list = self._get_band_observations(band, mindex)
            coadd_mb_obs_list.append(cobs_list)
            mb_obs_list.append(obs_list)

        meta_row = self._get_meta_row()
        meta_row['id'][0] = self.meds_list[0]['id'][mindex]
        meta_row['number'][0] = self.meds_list[0]['number'][mindex]
        meta_row['ra'][0] = self.meds_list[0]['ra'][mindex]
        meta_row['dec'][0] = self.meds_list[0]['dec'][mindex]

        # to account for max_cutouts limit, we count the actual number
        #meta_row['nimage_tot'][0,:] = numpy.array([self.meds_list[b]['ncutout'][mindex]-1 for b in xrange(self.conf['nband'])],dtype='i4')
        meta_row['nimage_tot'][0,:] = numpy.array([len(mb_obs_list[b]) for b in xrange(self.conf['nband'])],dtype='i4')

        meta = {'meta_data':meta_row,'meds_index':mindex,'id':self.meds_list[0]['id'][mindex],'obj_flags':0}

        coadd_mb_obs_list.update_meta_data(meta)
        mb_obs_list.update_meta_data(meta)

        return coadd_mb_obs_list, mb_obs_list

    def _reject_outliers(self, obs_list):
        for attr in ['weight','weight_raw','weight_us','weight_orig']:
            imlist=[]
            wtlist = []
            for obs in obs_list:
                if obs.meta['flags'] == 0 and hasattr(obs,attr) and getattr(obs,attr) is not None:
                    imlist.append(obs.image)
                    wtlist.append(getattr(obs,attr))

            # weight map is modified
            if len(wtlist) > 0:
                nreject=meds.reject_outliers(imlist,wtlist)
                if nreject > 0:
                    print('    rejected pixels using %s: %d' % (attr,nreject))

    def _get_image_flags(self, band, mindex):
        meds=self.meds_list[band]
        ncutout=meds['ncutout'][mindex]
        return numpy.zeros(ncutout, dtype='i8')
        #return numpy.zeros(self.conf['nband'])

    def _should_use_obs(self, band, mindex, icut):
        max_cutouts=self.conf['max_cutouts']
        if (max_cutouts is not None and icut > max_cutouts):
            return False
        else:
            return True

    def _get_band_observations(self, band, mindex):
        """
        Get an ObsList for the coadd observations in each band
        If psf fitting fails, the ObsList will be zero length
        note we have already checked that we have a coadd and a single epoch
        without flags
        """

        meds=self.meds_list[band]
        ncutout=meds['ncutout'][mindex]

        image_flags=self._get_image_flags(band, mindex)

        coadd_obs_list = ObsList()
        obs_list       = ObsList()

        fake=numpy.zeros((0,0))
        for icut in xrange(ncutout):

            flags=0
            if not self._should_use_obs(band, mindex, icut):
                obs = Observation(fake)
                flags = IMAGE_FLAGS
            else:
                iflags = image_flags[icut]
                if iflags != 0:
                    flags = IMAGE_FLAGS
                    obs = Observation(fake)
                else:
                    obs = self._get_band_observation(band, mindex, icut)
                    if obs is None:
                        flags = IMAGE_FLAGS
                        obs = Observation(fake)

            # fill the meta data
            self._fill_obs_meta_data(obs,band,mindex,icut)

            # set flags
            meta = {'flags':flags}
            obs.update_meta_data(meta)

            if icut==0:
                coadd_obs_list.append(obs)
            else:
                obs_list.append(obs)

        if self.conf['reject_outliers'] and len(obs_list) > 0:
            self._reject_outliers(obs_list)

        obs_list.update_meta_data({'band_num':band})
        coadd_obs_list.update_meta_data({'band_num':band})

        return coadd_obs_list, obs_list

    def _get_meds_orig_filename(self, meds, mindex, icut):
        """
        Get the original filename
        """
        return ''

    def _get_band_observation(self, band, mindex, icut):
        """
        Get an Observation for a single band.
        """

        meds=self.meds_list[band]

        bmask,skip = self._get_meds_bmask(meds, mindex, icut)
        if skip:
            return None


        wt,wt_us,wt_raw,seg,skip = self._get_meds_weight(meds, mindex, icut)
        if skip:
            return None

        im = self._get_meds_image(meds, mindex, icut)

        jacob = self._get_jacobian(meds, mindex, icut)


        if self.conf['ignore_zero_images'] and 0.0==im.sum():
            print("    image all zero, skipping")
            return None

        psf_obs = self._get_psf_observation(band, mindex, icut, jacob)

        obs=Observation(im,
                        weight=wt,
                        bmask=bmask,
                        jacobian=jacob,
                        psf=psf_obs)
        if wt_us is not None:
            obs.weight_us = wt_us
        else:
            obs.weight_us = None

        obs.weight_raw = wt_raw
        obs.seg = seg

        fname = self._get_meds_orig_filename(meds, mindex, icut)
        obs.filename=fname

        if 'trim_image' in self.conf:
            self._trim_obs(obs)

        return obs

    def _get_trimmed_startstop(self, dim, cen, new_dim):
        start=int(cen-new_dim/2.0+0.5)
        end=int(cen+new_dim/2.0+0.5)

        if start < 0:
            start=0
        if end > (dim-1):
            end=dim-1

        return start,end

    def _trim_obs(self, obs):
        dims=obs.image.shape
        new_dims=self.conf['trim_image']['dims']
        print("trimming obs:",new_dims)

        j=obs.jacobian
        cen=j.get_cen()

        rowstart, rowend=self._get_trimmed_startstop(
            dims[0],
            cen[0],
            new_dims[0],
        )
        colstart, colend=self._get_trimmed_startstop(
            dims[1],
            cen[1],
            new_dims[1],
        )

        obs.image  = obs.image[rowstart:rowend, colstart:colend]
        obs.weight = obs.weight[rowstart:rowend, colstart:colend]
        obs.weight_raw = obs.weight_raw[rowstart:rowend, colstart:colend]
        if obs.weight_us is not None:
            obs.weight_us = obs.weight_us[rowstart:rowend, colstart:colend]

        j.set_cen(
            row=cen[0]-rowstart,
            col=cen[1]-colstart,
        )

        print("new dims:",obs.image.shape)



    def _fill_obs_meta_data(self,obs, band, mindex, icut):
        """
        fill meta data to be included in output files
        """

        meds=self.meds_list[band]

        meta_row = self._get_epoch_meta_row()
        meta_row['id'][0] = meds['id'][mindex]
        meta_row['number'][0] = meds['number'][mindex]
        meta_row['band_num'][0] = band
        meta_row['cutout_index'][0] = icut
        meta_row['orig_row'][0] = meds['orig_row'][mindex,icut]
        meta_row['orig_col'][0] = meds['orig_col'][mindex,icut]
        file_id  = meds['file_id'][mindex,icut].astype('i4')
        meta_row['file_id'][0]  = file_id
        meta_row['pixel_scale'][0] = obs.get_jacobian().get_scale()
        meta={'icut':icut,
              'cutout_index':icut,
              'orig_start_row':meds['orig_start_row'][mindex, icut],
              'orig_start_col':meds['orig_start_col'][mindex, icut],
              'meta_data':meta_row,
              'id':meds['id'][mindex],
              'band_id':icut}
        obs.update_meta_data(meta)

    def _get_meds_image(self, meds, mindex, icut):
        """
        Get an image cutout from the input MEDS file
        """
        try:
            self.imname= os.path.basename(meds.get_source_path(mindex,icut))
        except IndexError:
            self.imname=''

        return meds.get_cutout(mindex, icut)

    def _badfrac_too_high(self, icut, nbad, shape, maxfrac, type):
        ntot=shape[0]*shape[1]
        frac = float(nbad)/ntot

        if maxfrac=='one-over-side':
            maxfrac=1.0/min(shape)

        if frac > maxfrac:
            print("    skipping cutout",icut,"due to high ",type,"frac:",frac)
            return True
        else:
            return False


    def _get_meds_bmask(self, meds, mindex, icut):
        """
        Get an image cutout from the input MEDS file
        """
        maxfrac=self.conf['max_bmask_frac']
        skip=False

        if 'bmask_cutouts' in meds._fits:
            bmask=meds.get_cutout(mindex, icut, type='bmask')
            bmask=numpy.array(bmask, dtype='i4', copy=False)

            if self.conf['symmetrize_bmask']:
                if bmask.shape[0] == bmask.shape[1]:
                    borig=bmask.copy()
                    rotmask=numpy.rot90(bmask)
                    bmask |= rotmask
                else:
                    raise RuntimeError("cannot symmetrize non-square bmask")

            w=numpy.where(bmask != 0)

            notok=self._badfrac_too_high(
                icut, w[0].size, bmask.shape, maxfrac, 'bmask'
            )
            if notok:
                skip=True
                return None,skip

            if 'bmask_skip_flags' in self.conf:
                w=numpy.where( (bmask & self.conf['bmask_skip_flags']) != 0)
                if w[0].size > 0:
                    print("    skipping cutout",icut,
                          "because mask bits set:",
                          self.conf['bmask_skip_flags'])
                    skip=True
                    return None,skip

            if 'central_bmask_radius' in self.conf:
                rad=self.conf['central_bmask_radius']
                if rad is not None:
                    row0 = meds['cutout_row'][mindex,icut]
                    col0 = meds['cutout_col'][mindex,icut]

                    row_start = _clip_pixel(row0-rad, bmask.shape[0])
                    row_end   = _clip_pixel(row0+rad, bmask.shape[0])
                    col_start = _clip_pixel(col0-rad, bmask.shape[1])
                    col_end   = _clip_pixel(col0+rad, bmask.shape[1])

                    bmask_sub = bmask[row_start:row_end,
                                      col_start:col_end]
                    wcen=numpy.where(bmask_sub != 0)
                    if wcen[0].size > 0:
                        print("    skipping cutout",icut,"due center masked")
                        skip=True
                        return None,skip

        else:
            bmask=None

        return bmask, skip

    def _clip_weight(self,wt):
        wt = wt.astype('f8', copy=False)

        w = numpy.where(wt < self.conf['min_weight'])
        if w[0].size > 0:
            wt[w] = 0.0

            if self.conf['symmetrize_weight']:
                #print("    symmetrizing weight")
                wt_rot = numpy.rot90(wt)
                w_rot = numpy.where(wt_rot < self.conf['min_weight'])
                wt[w_rot] = 0.0

        wt=wt.clip(min=0.0)
        return wt

    def _get_meds_weight(self, meds, mindex, icut):
        """
        Get a weight map from the input MEDS file
        """
        maxfrac=self.conf['max_zero_weight_frac']
        skip=False

        wt_raw = meds.get_cutout(mindex, icut, type='weight')
        if self.conf['region'] == 'mof':
            wt=wt_raw.copy()
            wt_us = meds.get_cweight_cutout_nearest(mindex, icut)
        elif self.conf['region'] == "cweight-nearest":
            wt = meds.get_cweight_cutout_nearest(mindex, icut)
            wt_us = None
        elif self.conf['region'] == 'seg_and_sky':
            wt=meds.get_cweight_cutout(mindex, icut)
            wt_us = None
        elif self.conf['region'] == 'weight':
            wt=wt_raw.copy()
            wt_us = None
        else:
            raise ValueError("no support for region type %s" % self.conf['region'])

        wt = self._clip_weight(wt)
        wt_raw = self._clip_weight(wt_raw)
        if wt_us is not None:
            wt_us = self._clip_weight(wt_us)

        try:
            seg = meds.interpolate_coadd_seg(mindex, icut)
        except:
            seg = meds.get_cutout(mindex, icut, type='seg')


        '''
        if self.conf['symmetrize_weight']:
            raise RuntimeError("this is bogus!  Need to zero the map not add")
            wt     = wt     + numpy.rot90(wt)
            wt_raw = wt_raw + numpy.rot90(wt_raw)

            if wt_us is not None:
                wt_us  = wt_us  + numpy.rot90(wt_us)
        '''

        # check raw weight map for zero pixels
        wzero=numpy.where(wt_raw == 0.0)

        notok=self._badfrac_too_high(
            icut, wzero[0].size, wt_raw.shape, maxfrac, 'zero weight'
        )
        if notok:
            skip=True

        return wt,wt_us,wt_raw,seg, skip

    def _get_jacobian(self, meds, mindex, icut):
        """
        Get a Jacobian object for the requested object
        """
        jdict = meds.get_jacobian(mindex, icut)
        jacob = Jacobian(row=jdict['row0'],
                         col=jdict['col0'],
                         dudrow=jdict['dudrow'],
                         dudcol=jdict['dudcol'],
                         dvdrow=jdict['dvdrow'],
                         dvdcol=jdict['dvdcol'])
        return jacob

    def _get_psf_observation(self, band, mindex, icut, image_jacobian):
        """
        Get an Observation representing the PSF and the "sigma"
        from the psf model
        """
        im, cen, sigma_pix, fname = self._get_psf_image(band, mindex, icut)

        psf_jacobian = image_jacobian.copy()
        psf_jacobian.set_cen(row=cen[0], col=cen[1])

        psf_obs = Observation(im, jacobian=psf_jacobian)
        psf_obs.filename=fname

        # convert to sky coords
        sigma_sky = sigma_pix*psf_jacobian.get_scale()

        psf_obs.update_meta_data({'sigma_sky':sigma_sky})
        psf_obs.update_meta_data({'Tguess':sigma_sky*sigma_sky})
        psf_obs.update_meta_data({'psf_norm':im.sum()})

        return psf_obs

    def _load_meds_files(self):
        """
        Load all listed meds files
        """
        self.meds_list=[]
        self.meds_meta_list=[]

        for i,funexp in enumerate(self.meds_files):
            f = os.path.expandvars(funexp)
            print('band %d meds: %s' % (i,f))
            medsi=meds.MEDS(f)
            medsi_meta=medsi.get_meta()

            if i==0:
                nobj_tot=medsi.size
            else:
                nobj=medsi.size
                if nobj != nobj_tot:
                    raise ValueError("mismatch in meds "
                                     "sizes: %d/%d" % (nobj_tot,nobj))
            self.meds_list.append(medsi)
            self.meds_meta_list.append(medsi_meta)

        self.nobj_tot = self.meds_list[0].size

def _clip_pixel(pixel, npix):
    pixel=int(pixel)
    if pixel < 0:
        pixel=0
    if pixel > (npix-1):
        pixel = (npix-1)
    return pixel


