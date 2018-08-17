#!/usr/bin/env python

import fitsio
import numpy 
import time
import argparse
import logging
import os,sys

import despyfits
from despyfits import maskbits
from despyfits import DESImage
from despyfits import compressionhdu as chdu
from despyastro import astrometry
from despyastro import zipper_interp as zipp
from despyastro import CCD_corners 

def build_parser():

    desc = """
    Creates a co-added MEF fits file from a set of single-plane fits
    containing the SCI/MSK/WGT planes. Interpolates the SCI plane
    using information in the 'custom'-weight mask, and also creates the MSK
    plane to be used by SExtractor for IMAFLAG_ISO.
    Felipe Menanteau (NCSA)
    """

    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--sci_file", default=None,required=True,
                        help="Input SCI fits file")
    parser.add_argument("--msk_file", default=None,required=False,
                        help="Input WGT fits file")
    parser.add_argument("--wgt_file", default=None,required=True,
                        help="Input WGT fits file")
    parser.add_argument("-o","--outname", default=None, required=True,
                        help="Name of output FITS file.")
    parser.add_argument("--clobber", action='store_true', default=False,
                        help="Clobber output MEF fits file")
    parser.add_argument("--add_noise", action='store_true', default=False,
                        help="Add Poisson Noise to the zipper")
    parser.add_argument("--xblock", default=1, type=int,
                        help="Block size of zipper in x-direction")
    parser.add_argument("--yblock", default=1, type=int,
                        help="Block size of zipper in y-direction")
    parser.add_argument("--ydilate", default=0, type=int,
                        help="Dilate pixels in the y-axis")
    parser.add_argument("--maxcols",dest="DEFAULT_MAXCOLS", default=100, type=int,
                        help="Widest feature to interpolate.  Default=None means no limit.")
    parser.add_argument("--mincols",dest="DEFAULT_MINCOLS", default=1, type=int,
                        help="Narrowest feature to interpolate.")
    parser.add_argument("--interp_image", action='store', choices=['WGT', 'MSK'], default='MSK',
                        help="Image to use that define pixels to interpolate over (MSK or WGT)")
    parser.add_argument("--region_file", default=None, type=str, required=False,
                        help="Write ds9 region file with interpolated region")

    # Keep zeros in SCI (yes/no)
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--keep_sci_zeros',    dest='keep_sci_zeros', action='store_true', default=True,
                       help="Keep zeros in SCI frame")
    group.add_argument('--no-keep_sci_zeros', dest='keep_sci_zeros', action='store_false')

    # Header options for DESDM Framework
    parser.add_argument("--band", default=None, type=str, required=False,
                        help="Add (optional) BAND to SCI header if not present")
    parser.add_argument("--magzero", default=None, type=float, required=False,
                        help="Add (optional) MAGZERO to SCI header")
    parser.add_argument("--tilename", default=None, type=str, required=False,
                        help="Add (optional) TILENAME to SCI header")
    parser.add_argument("--tileid", default=None, type=int, required=False,
                        help="Add (optional) TILE_ID to SCI header")
    return parser

def cmdline():
    
    parser = build_parser()
    args = parser.parse_args()

    # Sanity ckecks for inputs...
    # Make sure that outname is defined
    if not args.outname:
        raise ValueError('Undefined outname as %s' % args.outname)
    
    # Output file exits
    if os.path.isfile(args.outname) and args.clobber is False:
        raise ValueError("Output file exists, try --clobber option, no files created")
    return args

def create_logger(level=logging.NOTSET):

    logging.basicConfig(level=level,
                        format='[%(asctime)s] [%(levelname)s] %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger('Merge')
    return logger

def merge(**kwargs):

    """
    Merge the coadded science (SCI), mask (MSK) and weight (WGT) files
    into a single MEF file using fitsio, and perform zipper
    interpolation along columns on the SCI plane using the 'custom'
    MSK weight plane created by SWarp.

    Felipe Menanteau
    """

    sci_file    = kwargs.get('sci_file')
    msk_file    = kwargs.get('msk_file',None)
    wgt_file    = kwargs.get('wgt_file')
    logger      = kwargs.get('logger',None)
    outname     = kwargs.get('outname',False)
    interp_mask = kwargs.get('interp_mask',1)
    xblock      = kwargs.get('xblock',0)
    ydilate     = kwargs.get('ydilate',0)
    BADPIX_INTERP = kwargs.get('BADPIX_INTERP',None)
    # Header options (all optional)
    BAND        = kwargs.get('band',None)
    MAGZERO     = kwargs.get('magzero',None)
    TILENAME    = kwargs.get('tilename',None)
    TILEID      = kwargs.get('tileid',None)
    interp_image  = kwargs.get('interp_image',None)
    keep_sci_zeros = kwargs.get('keep_sci_zeros',True) 
    # The rest (and most) will be passed as kwargs to the zipper interpolation routine

    if not logger:
        logger = create_logger(level=logging.NOTSET)
        kwargs['logger'] = logger # Add logger to kwargs to pass on
        
    logger.info("Reading in %s" % sci_file)
    SCI,sci_hdr = fitsio.read(sci_file, ext=0, header=True)
    logger.info("Reading in %s" % wgt_file)
    WGT,wgt_hdr = fitsio.read(wgt_file, ext=0, header=True)
    
    if msk_file:
        logger.info("Reading in %s" % msk_file)
        MSK,msk_hdr = fitsio.read(msk_file, ext=0, header=True)
        # Make the MSK for MSK-like weight
        MSK = numpy.where(MSK == 0,1,0)
    else: # Create mask from WGT plane
        MSK = numpy.copy(WGT)
        MSK = numpy.where(MSK == 0,1,0)
        msk_hdr = wgt_hdr

    # Make sure that we do not interpolate over zeroes on the SCI frame
    if keep_sci_zeros:
        logger.info("Preserving zeros in SCI frame")
        MSK  = numpy.where(SCI == 0,0,MSK)
    else:
        logger.info("Ignoring zeros in SCI frame")
        
    # Define the mask we'll use for interpolation
    if interp_image == 'MSK':
        MSK_interp = MSK
    elif interp_image == 'WGT':
        MSK_interp = numpy.copy(WGT)
        MSK_interp = numpy.where(MSK_interp == 0,1,0)
        MSK_interp = numpy.where(SCI == 0,0,MSK_interp)
    else:
        exit("ERROR: wrong interp_image was passed")

    # Perform column interpolation -- Axis=2 -- only xblock > 0
    if xblock > 0:
        SCI,MSK_interp = zipp.zipper_interp(SCI,MSK_interp,interp_mask,axis=2,**kwargs)
        # Interpolate the WGT plane if we don't have a msk_file
        if not msk_file:
            WGT,MSK_interp = zipp.zipper_interp(WGT,MSK_interp,interp_mask,axis=2,**kwargs)
    else:
        logger.info("Skipping interpolation of SCI plane -- xblock=0")

    # In case the MSK_interp was dilated, we want to pass that info the MSK plane
    MSK = numpy.where(MSK_interp == 1,1,MSK)

    # Update compression settings
    logger.info("Updating compression settings")
    sci_hdr = DESImage.update_hdr_compression(sci_hdr,'SCI')
    msk_hdr = DESImage.update_hdr_compression(msk_hdr,'MSK')
    wgt_hdr = DESImage.update_hdr_compression(wgt_hdr,'WGT')

    # Use special compression when COMBINE_TYPE is 'CHI-MEAN' 
    if wgt_hdr.get('COMBINET') == 'CHI-MEAN': 
        logger.info("Overriding compression settings for WGT with 'CHI-MEAN' COMBINET will use GZIP")
        wgt_hdr['FZALGOR'] = 'GZIP_1'
        wgt_hdr['FZQVALUE'] = 0
        wgt_hdr['FZQMETHD'] = 'NO_DITHER'

    # Add corners, centers and extend
    logger.info("Updating Image Corners")
    sci_hdr = CCD_corners.update_DESDM_corners(sci_hdr,get_extent=True, verb=False)
    msk_hdr = CCD_corners.update_DESDM_corners(msk_hdr,get_extent=True, verb=False)

    # Add BAND if present
    if BAND:
        record={'name':'BAND', 'value':BAND, 'comment':'Short name for filter'}
        sci_hdr.add_record(record)

    # Add MAGZERO if present
    if MAGZERO:
        record={'name':'MAGZERO', 'value':MAGZERO, 'comment':'Mag Zero-point in magnitudes/s'}
        sci_hdr.add_record(record)

    # Add TILENAME if present
    if TILENAME:
        record={'name':'TILENAME', 'value':TILENAME, 'comment':'DES Tilename'}
        sci_hdr.add_record(record)

    # Add TILEID if present
    if TILEID:
        record={'name':'TILEID', 'value':TILEID, 'comment':'Tile ID for DES Tilename'}
        sci_hdr.add_record(record)

    # Insert EUPS PIPEPROD and PIPEVER to SCI HDU")
    if DESImage.pipekeys_write:
        logger.info("Inserting EUPS PIPEPROD and PIPEVER to SCI HDU")
        sci_hdr = DESImage.insert_eupspipe(sci_hdr)

    # Add to image history
    sci_hdr['HISTORY'] = time.asctime(time.localtime()) + \
                         ' column_interp over mask 0x{:04X}'.format(interp_mask)

    # Write it out now
    logger.info("Writing %s" % outname)
    ofits = fitsio.FITS(outname,'rw',clobber=True)
    ofits.write(SCI,extname='SCI',header=sci_hdr)
    ofits.write(MSK,extname='MSK',header=msk_hdr)
    ofits.write(WGT,extname='WGT',header=wgt_hdr)
    ofits.close()

    return

