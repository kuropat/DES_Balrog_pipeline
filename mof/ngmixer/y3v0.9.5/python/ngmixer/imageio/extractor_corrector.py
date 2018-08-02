from __future__ import print_function
import numpy
import os
import time
import fitsio
import meds

from ..render_ngmix_nbrs import RenderNGmixNbrs, RenderNGmixNbrsFromFile

# these are higher than anything in Y1 DES masks
NBRS_MASKED = 2**29
CEN_MODEL_MISSING = 2**30

class MEDSExtractorCorrector(meds.MEDSExtractor):
    """
    parameters
    ----------
    mof_file: string
        Result from an mof run on this tile.  It is assumed the
        model was 'cm', cmodel
    meds_file: string
        The MEDS file from which to extract a subset
    start: integer
        First index to extract
    end: intege
        Last index to extract
    sub_file: string
        Where to write the result

    reject_outliers: bool
        Set the weight to zero for pixels that are outliers.  Happens
        before bad pixel replacement

    replace_bad: bool
        If True, replace pixels the value for pixels with zero weight using the
        central model.  If bad_flags is also sent, also perform the replacement
        for pixels with those bits set.  If the central model did not converge
        (bad fit) then set the flag CEN_MODEL_MISSING in the bmask and force
        the weight to zero.  Default False.

    bad_flags: int or sequence
        If replace_bad is True, and bad_flags is sent, replacement also occurs
        for pixels with those flags set

    min_weight: float
        Min weight to consider "bad".  If the compression preserves
        zeros, this can be left at zero.  Default 0.0

    model: string
        Model used in the MOF, default 'cm'.  We will soon add
        this strig to the MOF outputs so we don't need the
        keyword.

    cleanup: bool
        Set to True to clean up the file when this object is
        cleaned up.
    verbose: bool
        Be verbose
    """
    def __init__(self,
                 mof_file,
                 meds_file,
                 start,
                 end,
                 sub_file,
                 reject_outliers=False,
                 replace_bad=False,
                 bad_flags=None,
                 min_weight=0.0,
                 # these are the bands in the mof
                 band_names = ['g','r','i','z'],
                 band=None,
                 model='cm',
                 cleanup=False,
                 make_plots=False,
                 copy_all=False,
                 verbose=False):

        self.cleanup=cleanup
        self.mof_file=mof_file

        self.band_names=band_names
        self.model=model

        self.reject_outliers=reject_outliers
        self.replace_bad=replace_bad

        self._set_bad_flags(bad_flags)

        self.min_weight=min_weight
        self.make_plots=make_plots
        self.verbose=verbose

        self._set_band(meds_file, band=band)

        self._setup_plotting(meds_file)

        # this will run extraction
        super(MEDSExtractorCorrector,self).__init__(
            meds_file,
            start,
            end,
            sub_file,
            cleanup=cleanup,
            copy_all=copy_all,
        )

    def _set_bad_flags(self, bad_flags):
        if bad_flags is not None:
            try:
                bflags=sum(bad_flags)
            except:
                bflags=bad_flags
        else:
            bflags=None

        self.bad_flags=bflags

    def _setup_plotting(self, meds_file):
        if self.make_plots:
            self.pdir='plots'
            if not os.path.exists(self.pdir):
                os.makedirs(self.pdir)

            bname=os.path.basename(meds_file)
            front=bname.replace('.fits.fz','').replace('.fits','')
            self.front=front

    def _reject_outliers(self, imlist, wtlist):
        """
        See meds.reject_outliers for the algorithm
        """
        nreject=meds.reject_outliers(imlist, wtlist)
        if nreject > 0:
            print("    rejected",nreject,"outliers")

    def _do_render_nbrs(self, mfile, mindex, icut):
        """
        render the neighbor and central galaxy

        Deal specially with the case of no neighbors, when
        there is not central fit in the neighbors structure
        """
        try:
            seg = mfile.interpolate_coadd_seg(mindex, icut)
        except:
            seg = mfile.get_cutout(mindex, icut, type='seg')

        ncutout=mfile['ncutout'][mindex]
        coadd_object_id=mfile['id'][mindex]

        print("  cutout %d/%d" % (icut+1,ncutout))

        res = self.renderer.render_nbrs(
            coadd_object_id,
            icut,
            seg,
            self.model,
            self.band,
            total=True,
            verbose=self.verbose,
        )
        if res is not None:
            cen_img, nbrs_img, nbrs_mask, nbrs_ids, pixel_scale = res

            if cen_img is None:
                print("    bad central fit")

        else:
            if self.verbose:
                print('    no nbrs, rendering central')
            nbrs_img=None
            nbrs_mask=None
            pixel_scale=None

            res = self.renderer.render_central(
                coadd_object_id,
                mfile,
                mindex,
                icut,
                self.model,
                self.band,
                seg.shape,
            )
            if res is None:
                cen_img,pixel_scale=None,None
            else:
                cen_img, pixel_scale=res

        return cen_img, nbrs_img, nbrs_mask, pixel_scale

    def _sub_nbrs(self, img, weight, nbrs_img, nbrs_mask, pixel_scale):
        """
        Subtract off the neighbor image and modify
        the weight map if any neighbor could not be rendered
        """
        # don't subtract in regions where the weight
        # is zero; this is currently a bug workaround because
        # in the MEDS maker we are including cutouts that
        # cross the edge
        pw=numpy.where(weight > self.min_weight)
        if pw[0].size > 0:
            # subtract neighbors if they exist
            img[pw] = img[pw] - nbrs_img[pw]*pixel_scale**2

        # if the nbr fit failed, nbrs_mask will be set to
        # zero in the seg map of the neighbor
        weight *= nbrs_mask

    def _mask_nbrs_mask(self, icut, nbrs_mask, bmask):
        """
        set the NBRS_MASKED flag in the bitmask where
        the neighbors could not be rendered
        """
        if nbrs_mask is not None:
            w=numpy.where(nbrs_mask != 1)
            if w[0].size > 0:
                print("     modifying",w[0].size,
                      "bmask pixels in cutout",icut,"for nbrs_mask")
                bmask[w] = bmask[w] | NBRS_MASKED

    def _replace_bad_pixels(self,
                            icut,
                            img,
                            weight,
                            bmask,
                            cen_img,
                            pixel_scale):
        """
        set masked or zero weight pixels to the value of the central.
        For codes that use the weight, such as max like, this makes
        no difference, but it may be important for codes that
        take moments or use FFTs

        Note if NBRS_MASKED  is set, the weight map is already zerod
        """

        imravel = img.ravel()
        wtravel = weight.ravel()
        bmravel = bmask.ravel()

        weight_logic = (wtravel <= self.min_weight)
        wbad_wt,  = numpy.where(weight_logic==True)
        wgood_wt, = numpy.where(weight_logic==False)

        logic = weight_logic
        if self.bad_flags is not None:
            bmask_logic  = (bmravel & self.bad_flags) != 0
            logic = logic & bmask_logic

        wbad,=numpy.where(logic)

        if wbad.size > 0 and wbad.size != imravel.size:

            if cen_img is None:

                print("        bad central fit for",icut,
                      "could not replace bad pixels. forcing "
                      "weight to zero")
                bmravel[wbad] = bmravel[wbad] | CEN_MODEL_MISSING

                if self.bad_flags is not None:
                    # make sure we don't use these pixels, at least in
                    # maximum likelihood codes
                    wtravel[wbad] = 0.0

            else:

                # weight map for adding noise; use median where the
                # weight is zero
                weight_for_noise = wtravel.copy()

                if wbad_wt.size > 0 and wgood_wt.size > 0:

                    medwt = numpy.median(weight_for_noise[wgood_wt])
                    print("        replacing",wbad_wt.size,
                          "weight in replaced pixels with",medwt)
                    weight_for_noise[wbad_wt] = medwt

                # replace zero weight and bad pixels with model
                scaled_cen = cen_img.ravel()*pixel_scale**2
                print("        setting",wbad.size,
                      "bad pixels in cutout",icut,
                      "to central model")
                imravel[wbad] = scaled_cen[wbad]

                # add appropriate noise as well
                print("        adding noise to",wbad.size,"replaced")
                err = numpy.sqrt( 1.0/weight_for_noise )
                noise = numpy.random.normal(
                    loc=0.0,
                    scale=1.0,
                    size=err.size,
                )
                noise *= err
                imravel[wbad] += noise[wbad]


    def _extract(self):
        # first do the ordinary extraction
        super(MEDSExtractorCorrector,self)._extract()

        # self.sub_file should now have been closed
        time.sleep(2)
        self._set_renderer()

        with fitsio.FITS(self.sub_file,mode='rw') as fits:
            mfile = RefMeds(fits)
            cat = mfile.get_cat()

            nobj=cat.size

            # loop through objects and correct each
            for mindex in xrange(cat['id'].size):

                coadd_object_id=cat['id'][mindex]
                ncutout=cat['ncutout'][mindex]
                box_size=cat['box_size'][mindex]
                start_row=cat['start_row'][mindex]

                print("%d/%d  %d" % (mindex+1, nobj, coadd_object_id))
                #if ncutout > 1 and box_size > 0:
                if ncutout > 0 and box_size > 0:

                    #imlist = mfile.get_cutout_list(mindex, type='image')[1:]
                    #wtlist = mfile.get_cutout_list(mindex, type='weight')[1:]
                    #bmlist = mfile.get_cutout_list(mindex, type='bmask')[1:]
                    imlist = mfile.get_cutout_list(mindex, type='image')
                    wtlist = mfile.get_cutout_list(mindex, type='weight')
                    bmlist = mfile.get_cutout_list(mindex, type='bmask')

                    if self.reject_outliers:
                        self._reject_outliers(imlist, wtlist)

                    # do coadd too
                    # if we ran MOF without fitting the coadd, we will get
                    # a warning "central not in epochs data"
                    for icut in xrange(ncutout):

                        #img_orig = imlist[icut-1]
                        #weight   = wtlist[icut-1]
                        #bmask    = bmlist[icut-1]

                        img_orig = imlist[icut]
                        weight   = wtlist[icut]
                        bmask    = bmlist[icut]

                        img=img_orig.copy()

                        cen_img, nbrs_img, nbrs_mask, pixel_scale= \
                                self._do_render_nbrs(mfile, mindex, icut)

                        if nbrs_img is not None:
                            self._sub_nbrs(
                                img, weight, nbrs_img, nbrs_mask, pixel_scale,
                            )

                        # zeros the weight and sets flag for neighbors that were
                        # not rendered
                        self._mask_nbrs_mask(icut, nbrs_mask, bmask)

                        if self.replace_bad:
                            self._replace_bad_pixels(
                                icut,
                                img,
                                weight,
                                bmask,
                                cen_img,
                                pixel_scale,
                            )

                        # now overwrite pixels on disk
                        fits['image_cutouts'].write(img.ravel(),
                                                    start=start_row[icut])
                        fits['weight_cutouts'].write(weight.ravel(),
                                                     start=start_row[icut])
                        fits['bmask_cutouts'].write(bmask.ravel(),
                                                    start=start_row[icut])

                        if self.make_plots:
                            self._doplots(
                                coadd_object_id,
                                mindex, icut, img_orig, img, bmask, weight,
                            )


                else:
                    # we always want to see this
                    print("    not writing ncutout:",ncutout,"box_size:",box_size)


    def _extract_old(self):
        # first do the ordinary extraction
        super(MEDSExtractorCorrector,self)._extract()

        # self.sub_file should now have been closed
        time.sleep(2)
        self._set_renderer()

        with fitsio.FITS(self.sub_file,mode='rw') as fits:
            mfile = RefMeds(fits)
            cat = mfile.get_cat()

            nobj=cat.size

            # loop through objects and correct each
            for mindex in xrange(cat['id'].size):

                coadd_object_id=cat['id'][mindex]
                ncutout=cat['ncutout'][mindex]
                box_size=cat['box_size'][mindex]
                start_row=cat['start_row'][mindex]

                print("%d/%d  %d" % (mindex+1, nobj, coadd_object_id))
                if ncutout > 1 and box_size > 0:
                    for icut in xrange(1,ncutout):

                        try:
                            seg = mfile.interpolate_coadd_seg(mindex, icut)
                        except:
                            seg = mfile.get_cutout(mindex, icut, type='seg')

                        print("  cutout %d/%d" % (icut+1,ncutout))

                        res = self.renderer.render_nbrs(
                            coadd_object_id,
                            icut,
                            seg,
                            self.model,
                            self.band,
                            total=True,
                            verbose=self.verbose,
                        )
                        if res is not None:
                            cen_img, nbrs_img, nbrs_mask, nbrs_ids, pixel_scale = res

                            if cen_img is None:
                                print("    bad central fit")

                        else:
                            if self.verbose:
                                print('    no nbrs, rendering central')
                            nbrs_img=None
                            nbrs_mask=None
                            nbrs_ids=None
                            pixel_scale=None

                            res = self.renderer.render_central(
                                coadd_object_id,
                                mfile,
                                mindex,
                                icut,
                                self.model,
                                self.band,
                                seg.shape,
                            )
                            if res is None:
                                cen_img,pixel_scale=None,None
                            else:
                                cen_img, pixel_scale=res

                        img_orig = mfile.get_cutout(mindex, icut, type='image')
                        weight   = mfile.get_cutout(mindex, icut, type='weight')
                        bmask    = mfile.get_cutout(mindex, icut, type='bmask')

                        img=img_orig.copy()
                        imravel = img.ravel()
                        wtravel = weight.ravel()
                        bmravel = bmask.ravel()

                        if nbrs_img is not None:

                            # don't subtract in regions where the weight
                            # is zero; this is currently a bug workaround because
                            # in the MEDS maker we are including cutouts that
                            # cross the edge
                            pw=numpy.where(weight > self.min_weight)
                            if pw[0].size > 0:
                                # subtract neighbors if they exist
                                img[pw] = img[pw] - nbrs_img[pw]*pixel_scale**2

                            # if the nbr fit failed, nbrs_mask will be set to
                            # zero in the seg map of the neighbor
                            weight *= nbrs_mask

                        # mark zero weight pixels in the mask; this is to
                        # deal with forgetting to do so in the meds code
                        # when we set the weight to zero
                        # we should do this in the meds reader
                        weight_logic = (weight <= self.min_weight)
                        wbad_wt=numpy.where(weight_logic)
                        if wbad_wt[0].size > 0:
                            bmask[wbad_wt] = bmask[wbad_wt] | ZERO_WEIGHT

                        if nbrs_mask is not None:
                            w=numpy.where(nbrs_mask != 1)
                            if w[0].size > 0:
                                print("     modifying",w[0].size,
                                      "bmask pixels in cutout",icut,"for nbrs_mask")
                                bmask[w] = bmask[w] | NBRS_MASKED

                        # set masked or zero weight pixels to the value of the central.
                        # For codes that use the weight, such as max like, this makes
                        # no difference, but it may be important for codes that 
                        # take moments or use FFTs
                        # note zero weight has been marked in the bmask
                        if self.replace_bad:

                            # working with unravelled images for efficiency
                            bmask_logic  = (bmravel != 0)
                            wbad,=numpy.where(bmask_logic)

                            if wbad.size > 0:
                                if cen_img is None:
                                    print("        bad central fit for",icut,
                                          "could not replace bad pixels")
                                    bmravel[wbad] = bmravel[wbad] | CEN_MODEL_MISSING

                                else:

                                    wgood_wt=numpy.where(weight_logic == False)
                                    if wgood_wt[0].size > 0 and wbad_wt[0].size > 0:
                                        # fix the weight map.  We are going to add noise
                                        # as well as signal in the bad pixels.  The problem
                                        # pixels will still be marked in the bad pixel mask
                                        print("        replacing",wbad_wt[0].size,
                                              "weight in replaced pixels")
                                        medwt = numpy.median(weight[wgood_wt])
                                        weight[wbad_wt] = medwt

                                    scaled_cen = cen_img.ravel()*pixel_scale**2
                                    print("        setting",wbad.size,
                                          "bad pixels in cutout",icut,
                                          "to central model")
                                    imravel[wbad] = scaled_cen[wbad]


                                    print("        adding noise to",wbad.size,"replaced")
                                    err = numpy.sqrt( 1.0/wtravel[wbad] )
                                    rand = numpy.random.normal(
                                        loc=0.0,
                                        scale=1.0,
                                        size=err.size,
                                    )
                                    rand *= err
                                    imravel[wbad] = imravel[wbad] + rand


                        # now overwrite pixels on disk
                        fits['image_cutouts'].write(imravel,  start=start_row[icut])
                        fits['weight_cutouts'].write(wtravel, start=start_row[icut])
                        fits['bmask_cutouts'].write(bmravel,  start=start_row[icut])

                        if self.make_plots:
                            self._doplots(
                                coadd_object_id,
                                mindex, icut, img_orig, img, bmask, weight,
                            )


                else:
                    # we always want to see this
                    print("    not writing ncutout:",ncutout,"box_size:",box_size)

    def _load_ngmix_data(self):
        #self.fit_data = fitsio.read(self.mof_file)
        #self.nbrs_data = fitsio.read(self.mof_file,ext='nbrs_data')
        #self.epoch_data = fitsio.read(self.mof_file,ext='epoch_data')
        pass

    def _set_renderer(self):
        self._load_ngmix_data()

        # build the renderer, set options
        conf = {'unmodeled_nbrs_masking_type':'nbrs-seg'}
        self.renderer = RenderNGmixNbrsFromFile(
            self.mof_file,
            **conf
        )
        """
        self.renderer = RenderNGmixNbrs(
            self.fit_data,
            self.nbrs_data,
            epoch_data=self.epoch_data,
            **conf
        )
        """


    def _set_band(self, meds_file, band=None):
        if band is not None:
            self.band=band
            return

        # get the band for the file
        band = -1
        for band_name in self.band_names:
            #btest = '-%s-' % band_name
            btest = '_%s_' % band_name
            if btest in meds_file:
                band = self.band_names.index(band_name)
                break
        if band == -1:
            raise ValueError("Could not find band for file '%s'!" % meds_file)

        self.band = band

    def _doplots(self, id, mindex, icut, img_orig, img, bmask, weight):
        import images
        imdiff = img_orig-img

        lbmask = numpy.log10( 1.0 + bmask )
        imlist=[img_orig, img, imdiff, lbmask, weight]
        imlist=[t - t.min() for t in imlist]

        titles=['image','imfix','image-imfix','bmask','weight']

        title='%d-%06d-%02d' % (id, mindex, icut)
        pname='%s-%d-%06d-%02d.png' % (self.front,id, mindex, icut)
        pname=os.path.join(self.pdir, pname)
        tab=images.view_mosaic(
            imlist,
            titles=titles,
            show=False,
        )
        tab.title=title

        ctab=images.multiview(img,show=False)
        tab[1,2] = ctab[0,1]

        print("    writing:",pname)
        tab.write_img(800,800,pname)

class RefMeds(meds.MEDS):
    """
    meds file object that accepts a ref to an open fitsio FITS object

    does not close the file if used in python context
    """
    def __init__(self,fits,filename=None):
        self._filename = filename
        self._fits = fits
        self._cat = self._fits["object_data"][:]
        self._image_info = self._fits["image_info"][:]
        self._meta = self._fits["metadata"][:]

    def close(self):
        pass


