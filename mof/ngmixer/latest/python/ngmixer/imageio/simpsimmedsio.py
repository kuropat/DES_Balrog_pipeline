"""
simple MEDS simulation file
"""
from __future__ import print_function
import os
import numpy
import copy
import fitsio

# meds and ngmix imports
import meds

# local imports
from .medsio import MEDSImageIO

PSF_IND_FIELD='ind_psf'
PSF_IM_FIELD='psf_im'

class SimpSimMEDSImageIO(MEDSImageIO):
    def _set_defaults(self):
        super(SimpSimMEDSImageIO, self)._set_defaults()
        self.conf['psf_ind_field'] = self.conf.get('psf_ind_field',PSF_IND_FIELD)
        self.conf['psf_im_field'] = self.conf.get('psf_im_field',PSF_IM_FIELD)
        self.conf['psfs_in_file'] = self.conf.get('psfs_in_file',False)

        if self.conf['psfs_in_file']:
            self.conf['psfs_are_psfex'] = self.conf.get('psfs_are_psfex',False)

        self.conf['center_psf'] = self.conf.get('center_psf',False)

    def _load_psf_data(self):
        if not self.conf['psfs_in_file']:
            if 'psf_file' in self.extra_data:
                self.psf_file = self.extra_data['psf_file']
            else:
                pth,bname = os.path.split(self.meds_files_full[0])
                bname = bname.replace('meds','psf')
                self.psf_file = os.path.join(pth,bname)
            print('psf file: %s' % self.psf_file)
            self.psf_data = fitsio.read(self.psf_file)

        elif self.conf['psfs_are_psfex']:
            # load all the PSFEx objects from file extensions
            last_ext = self.meds_list[0]._fits[-1]
            print(last_ext)
            extname = last_ext.get_extname()

            self._psf_ext_front=extname[0:12]
            self._psf_ext_end='psfcat'

            self.psfex_lists = self._get_psfex_lists()

    def _get_psfex_lists(self):
        """
        Load psfex objects for each of the SE images
        include the coadd so we get  the index right
        """
        print('loading psfex')

        psfex_lists=[]
        for band in self.iband:
            meds=self.meds_list[band]

            psfex_list = self._get_psfex_objects(meds)
            psfex_lists.append( psfex_list )

        return psfex_lists

    def _get_psfex_objects(self, meds):
        """
        Load psfex objects for all images, including coadd
        """
        import psfex

        psfex_list=[]

        info=meds.get_image_info()
        nimage=info.size
        for i in xrange(nimage):

            impath=info['image_path'][i]
            psf_ext = self._psfex_ext_from_image_path(info, impath)

            pex = psfex.PSFEx(
                meds._fits,
                ext=psf_ext,
            )

            psfex_list.append(pex)

        return psfex_list

    def _psfex_ext_from_image_path(self, info, impath):
        ccdname=os.path.basename( info['image_path'][35].strip() ).split('.')[0]

        return '%s_%s_%s' % (self._psf_ext_front, ccdname, self._psf_ext_end)

    def _get_psfex_image(self, band, mindex, icut):
        """
        Get an image representing the psf
        """

        meds=self.meds_list[band]
        file_id=meds['file_id'][mindex,icut]

        pex=self.psfex_lists[band][file_id]

        row=meds['orig_row'][mindex,icut]
        col=meds['orig_col'][mindex,icut]

        if self.conf['center_psf']:
            use_row,use_col=round(row),round(col)
        else:
            use_row,use_col=row,col

        im=pex.get_rec(use_row,use_col)
        cen=pex.get_center(use_row,use_col)

        im=im.astype('f8', copy=False)

        sigma_pix=pex.get_sigma()

        if 'trim_psf' in self.conf and icut > 0:
            im,cen=self._trim_psf(im, cen)

        return im, cen, sigma_pix, pex['filename']

    def _get_psf_image(self, band, mindex, icut):
        """
        Get an image representing the psf
        """
        
        if self.conf['psfs_in_file']:
            if self.conf['psfs_are_psfex']:
                return self._get_psfex_image(band, mindex, icut)
            else:
                im = self.meds_list[band].get_psf(mindex,icut)
                pfile = self.meds_files[band]
        else:
            meds=self.meds_list[band]
            psf_ind_field=self.conf['psf_ind_field']
            
            ind_psf = meds[psf_ind_field][mindex,icut]
            
            psf_im_field=self.conf['psf_im_field']
            im = self.psf_data[psf_im_field][ind_psf].copy()
            
            pfile = self.psf_file

        im /= im.sum()
        cen = ( numpy.array(im.shape) - 1.0)/2.0
        sigma_pix = 2.5
            
        return im, cen, sigma_pix, pfile
