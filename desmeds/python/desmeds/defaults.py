"""
default configuration

These get overridden by user configuration on constructing a DESMEDSMaker
object

"""

__version__= "0.9.2.2"

# sextractor names
default_config = {
    # we read the coadd catalog from detband
    'refband':'i',
    
    # types of cutout images to make
    'cutout_types': ['image','weight','seg','bmask'],

    # row,col names in the coadd catalog
    'row_name':'y_image',
    'col_name':'x_image',

    # for fpacking the file
    'fpack_dims': [10240,1],

    # we will put everything onto this magnitude scale
    'magzp_ref':30.0,

    # scamp uses 1-offset positions
    'position_offset':1.0,

    #
    # box size defaults
    #

    # number of object "sigma" to draw the box
    'sigma_fac':5.0,

    # 2**N or 3*2**N for fast FFTS
    'allowed_box_sizes': [
        2,3,4,6,8,12,16,24,32,48,
        64,96,128,192,256,
        384,512,768,1024,1536,
        2048,3072,4096,6144
    ],
    'min_box_size': 32,
    'max_box_size': 256,

    # astrometry gets refined during coaddition
    'use_astro_refine': True,

    # extensions for single-epoch images
    'se_image_ext': 1,
    'se_weight_ext': 3,

    'se_bmask_ext': 2,

    'se_bkg_ext': 1,
    'se_seg_ext': 1,

    # coadd images
    'coadd_image_ext': 1,
    'coadd_weight_ext': 2,
    # currentlly no bitmask for coadds!
    'coadd_bmask_ext': -1,

    # coadds are already background subtracted
    'coadd_bkg_ext': -1,

    'coadd_seg_ext': 1,

}
