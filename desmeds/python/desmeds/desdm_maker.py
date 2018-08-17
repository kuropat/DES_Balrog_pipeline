from __future__ import print_function
import os
from os.path import basename
import numpy
from numpy import zeros, sqrt, log, vstack, array
import subprocess
import shutil
import yaml

import fitsio
import esutil as eu

import meds
from meds.util import \
    make_wcs_positions, \
    get_meds_input_struct, \
   get_image_info_struct

from . import blacklists
from . import util


from . import files
from .defaults import default_config

from .files import \
        TempFile, \
        StagedInFile, \
        StagedOutFile

# desdb is not needed in all scenarios

try:
    import desdb
except ImportError:
    pass


fwhm_fac = 2*sqrt(2*log(2))

from .maker import DESMEDSMaker

class DESMEDSMakerDESDM(DESMEDSMaker):
    """
    This is the class for use by DESDM.  For this version,
    all inputs are explicit rather than relying on database
    queries

    No "stubby" meds file is created, because DESDM does
    not allow pipelines

    parameters
    ----------
    medconf: string
        path to a meds config file.  see docs for DESMEDSMaker
    fileconf: string
        path to a yaml file configuration

        Required fields in the yaml file:
            band: band in string form
            coadd_image_url: string
            coadd_seg_url: string
            coadd_magzp: float
            nwgint_flist: string
                path to the nwgint file list
            seg_flist: string
                path to the seg file list
            bkg_flist: string
                path to the bkg file list
            meds_url: string
                path to the output meds file
    """
    def __init__(self,
                 medsconf,
                 fileconf,
                 tmpdir=None):

        self.medsconf=medsconf
        self.fileconf=fileconf
        self.tmpdir=tmpdir

        self._load_config(medsconf)
        self._load_file_config(fileconf)

        self._set_extra_config('none', self.file_dict['band'])

        # not relevant for this version
        self.DESDATA = 'rootless'

    def go(self):
        """
        make the MEDS file
        """

        self._load_coadd_info()
        self._read_coadd_cat()
        self._build_image_data()
#        print(" build_meta_data \n")
        self._build_meta_data()
#        print("before build_object_data \n")
        self._build_object_data()

        self._write_meds_file() # does second pass to write data

    def _get_image_id_len(self, srclist):
        """
        for y3 using string ids
        """
        image_id_len=len(self.cf['image_id'])

        slen = len(self._get_portable_url(self.cf,'image_url'))
        for s in srclist:
            tlen = len(s['id'])
            if tlen > image_id_len:
                image_id_len = tlen

        return image_id_len

    def _load_coadd_info(self):
        """
        Mock up the results of querying the database for Coadd
        info
        """
        print('getting coadd info and source list')

        fd=self.file_dict
        cf={}

        iid = self._get_filename_as_id(fd['coadd_image_url'])

        cf['image_url'] = fd['coadd_image_url']
        cf['seg_url']   = fd['coadd_seg_url']
        cf['image_id']  = iid

        # probably from from header MAGZERO
        cf['magzp']     = fd['coadd_magzp']

        cf['srclist'] = self._load_srclist()

        # In this case, we can use refband==input band, since
        # not using a db query or anything
        self.cf=cf
        self.cf_refband=cf

    def _read_coadd_cat(self):
        """
        read the DESDM coadd catalog, sorting by the number field (which
        should already be the case)
        """

        fname=self.file_dict['coadd_cat_url']

        print('reading coadd cat:',fname)
        self.coadd_cat = fitsio.read(fname, lower=True)

        # sort just in case, not needed ever AFIK
        q = numpy.argsort(self.coadd_cat['number'])
        self.coadd_cat = self.coadd_cat[q]

    def _get_srclist(self):
        """
        mock up the interface for the Coadd class
        """
        return self.cf['srclist']

    def _load_srclist(self):
        """
        get all the necessary information for each source image
        """
        # this is a list of dicts
        srclist=self._load_nwgint_info()
        nepoch = len(srclist)

        # now add in the other file types
        bkg_info=self._read_generic_flist('bkg_flist')
        seg_info=self._read_generic_flist('seg_flist')

        if len(bkg_info) != nepoch:
            raise ValueError("bkg list has %d elements, nwgint "
                             "list has %d elements" % (len(bkg_info),nepoch))
        if len(seg_info) != nepoch:
            raise ValueError("seg list has %d elements, nwgint "
                             "list has %d elements" % (len(seg_info),nepoch))

        for i,src in enumerate(srclist):
            src['red_bkg'] = bkg_info[i]
            src['red_seg'] = seg_info[i]

        return srclist

    def _read_generic_flist(self, key):
        """
        read a list of file paths, one per line
        """
        fname=self.file_dict[key]
        print("reading:",key)

        flist=[]
        with open(fname) as fobj:
            for line in fobj:
                line=line.strip()
                if line=='':
                    continue

                flist.append(line)
        return flist

    def _extract_nwgint_line(self, line):
        """
        the nwgint (red image) lines are 
            path magzp
        """
        line=line.strip()
        if line=='':
            return None,None

        ls=line.split()
        if len(ls) != 2:
            raise ValueError("got %d elements for line in "
                             "nwgint list: '%s'" % line)

        path=ls[0]
        magzp=float(ls[1])

        return path, magzp


    def _load_nwgint_info(self):
        """
        Load all meta information needed from the
        ngmwint files
        """
        fname=self.file_dict['nwgint_flist']
        print("reading nwgint list and loading headers:",fname)

        red_info=[]
        with open(fname) as fobj:
            for line in fobj:

                path, magzp = self._extract_nwgint_line(line)
                if path==None:
                    continue

                sid = self._get_filename_as_id(path)

                # now mock up the structure of the Coadd.srclist

                wcs_hdr = fitsio.read_header(path, ext=self['se_image_ext'])
                wcs_header = util.fitsio_header_to_dict(wcs_hdr)

                s={
                    'id':sid,
                    'flags':0,  # assume no problems!
                    'red_image':path,
                    'magzp':magzp,
                    'wcs_header':wcs_header,
                }

                red_info.append(s)

        return red_info

    def _get_coadd_objects_ids(self):
        """
        mock up the query to the database
        """

        dt=[
            ('object_number','i4'),
            ('coadd_objects_id','i8')
        ]

        nobj=self.coadd_cat.size

        iddata=numpy.zeros(nobj, dtype=dt)

        idmap = fitsio.read(
            self.file_dict['coadd_object_map'],
            lower=True,
        )

        s=numpy.argsort(idmap['object_number'])

        iddata['object_number']    = idmap['object_number'][s]
        iddata['coadd_objects_id'] = idmap['id'][s] + 1

        return iddata

    def _get_portable_url(self, file_dict, name):
        """
        We don't have DESDATA defined when DESDM is running
        the code, so just return the path
        """

        return file_dict[name]

    def _load_config(self, medsconf):
        """
        load the default config, then load the input config

        note desdm names things randomly so we can't assert
        any consistency between the internal config name
        and the file name
        """

        self.update(default_config)

        if isinstance(medsconf,dict):
            conf=medsconf
        else:
            with open(medsconf) as fobj:
                conf=yaml.load( fobj )

        util.check_for_required_config(conf, ['medsconf'])
        self.update(conf)


    def _load_file_config(self, fileconf):
        """
        load the yaml file config
        """
        with open(fileconf) as fobj:
            self.file_dict=yaml.load( fobj )

    def _write_meds_file(self):
        """
        write the data using the MEDSMaker
        """

        maker=meds.MEDSMaker(
            self.obj_data,
            self.image_info,
            config=self,
            meta_data=self.meta_data,
        )

        fname=self.file_dict['meds_url']
        print("writing MEDS file:",fname)

        # this will do nothing if tmpdir is None; sf.path will
        # in fact equal fname and no move is performed

        with StagedOutFile(fname,tmpdir=self.tmpdir) as sf:
            if sf.path[-8:] == '.fits.fz':
                local_fitsname = sf.path.replace('.fits.fz','.fits')

                with TempFile(local_fitsname) as tfile:
                    maker.write(tfile.path)

                    # this will fpack to the proper path, which
                    # will then be staged out if tmpdir is not None
                    # if the name is wrong, the staging will fail and
                    # an exception raised
                    self._fpack_file(tfile.path)

            else:
                maker.write(sf.path)

    def _fpack_file(self, fname):
        cmd='fpack %s' % fname
        print("fpacking with command: '%s'" % cmd)
        subprocess.check_call(cmd,shell=True)

class Preparator(dict):
    """
    class to prepare inputs for the DESDM version
    of the MEDS maker

    This is not used by DESDM, but is useful for testing
    outside of DESDM
    N.K. We nedd to copy all but coadd and create part of the file config yaml file
    The second part of it as well as object map will be created after coadd and sextractor runs
    on the nullweighted files
    
    TODO: 
        - write psf map file
        - write part of file config
    """
    def __init__(self, medsconf, tilename, band):
        from .coaddinfo import Coadd
        from .coaddsrc import CoaddSrc

        if isinstance(medsconf,dict):
            conf=medsconf
        else:
            conf=files.read_meds_config(medsconf)
        self.update(conf)

        self['tilename']=tilename
        self['band']=band

        csrc=CoaddSrc(
            self['medsconf'],
            self['tilename'],
            self['band'],
            campaign=self['campaign'],
        )

        self.coadd=Coadd(
            self['medsconf'],
            self['tilename'],
            self['band'],
            campaign=self['campaign'],
            sources=csrc,
        )
        self['nullwt_dir']=files.get_nullwt_dir(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )


    def go(self):
        """
        download the data and make the null weight images
        """
        print("downloading all data")
        info = self.coadd.download()
        print(info)
#        self._make_objmap(info)
        self._copy_psfs(info)
        self._make_nullwt(info)

        fileconf=self._write_file_config(info)

        self._write_nullwt_flist(info['src_info'], fileconf)
        self._write_seg_flist(info['src_info'], fileconf)
        self._write_bkg_flist(info['src_info'], fileconf)

    def clean(self):
        """
        remove all sources and nullwt files
        """
        self.clean_sources()
        self.clean_nullwt()

    def clean_sources(self):
        """
        remove the downloaded source files
        """
        self.coadd.clean()

    def clean_nullwt(self):
        """
        remove all the generated nullwt files
        """
        print("removing nullwt images:",self['nullwt_dir'])
        shutil.rmtree(self['nullwt_dir'])



    def _make_objmap(self, info):
        fname=files.get_desdm_objmap(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )
        print("File name:",fname)
        if not os.path.exists(fname):

            dir=os.path.dirname(fname)
            if not os.path.exists(dir):
                print("making directory:",dir)
                os.makedirs(dir)

            sources = self.coadd.get_sources()
#            objmap = sources.cache.get_objmap(info)
            objmap = self.coadd.get_objmap(info)
#            print(objmap)
            print("writing objmap:",fname)
            fitsio.write(fname, objmap, extname='OBJECTS',clobber=True)

    def _write_nullwt_flist(self, src_info, fileconf):
        fname=fileconf['nwgint_flist']
        print("writing:",fname)
        with open(fname, 'w') as fobj:
            for sinfo in src_info:
                fobj.write("%s %.16g\n" % (sinfo['nullwt_path'], sinfo['magzp'] ))

    def _write_seg_flist(self, src_info, fileconf):
        fname=fileconf['seg_flist']
        print("writing:",fname)
        with open(fname, 'w') as fobj:
            for sinfo in src_info:
                fobj.write("%s\n" % sinfo['seg_path'])

    def _write_bkg_flist(self, src_info, fileconf):
        fname=fileconf['bkg_flist']
        print("writing:",fname)
        with open(fname, 'w') as fobj:
            for sinfo in src_info:
                fobj.write("%s\n" % sinfo['bkg_path'])

    def _write_file_config(self, info):
        fname=files.get_desdm_file_config(
            self['medsconf'],
            self['tilename'],
            self['band'],
        ) 
        files.makedir_fromfile(fname)
        nullwt_flist=files.get_desdm_nullwt_flist(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )
        seg_flist=files.get_desdm_seg_flist(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )
        bkg_flist=files.get_desdm_bkg_flist(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )
#        objmap=files.get_desdm_objmap(
#            self['medsconf'],
#            self['tilename'],
#            self['band'],
#        )


        meds_file=files.get_meds_file(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )
        output={
            'band':self['band'],
            'tilename':self['tilename'],
            'coadd_image_url':info['image_path'],
            'coadd_cat_url':info['cat_path'],
            'coadd_seg_url':info['seg_path'],
            'coadd_magzp':info['magzp'],
            'nwgint_flist':nullwt_flist,
            'seg_flist':seg_flist,
            'bkg_flist':bkg_flist,
            'meds_url':meds_file,
        }
#            'coadd_object_map':objmap,
        print ("output")
        print (output)
        print("writing file config:",fname)
        with open(fname,'w') as fobj:
            for key,value in output.iteritems():
                if key=="coadd_magzp":
                    value = '%.16g' % value

                fobj.write("%s: %s\n" % (key, value))

        return output

    def _make_nullwt(self, info):

        src_info=info['src_info']
        self._add_nullwt_paths(src_info)

        dir=self['nullwt_dir']
        if not os.path.exists(dir):
            print("making directory:",dir)
            os.makedirs(dir)

        print("making nullweight images")
        for sinfo in src_info:
            if os.path.exists(sinfo['nullwt_path']):
                continue

            cmd = _NULLWT_TEMPLATE % sinfo

            subprocess.check_call(cmd,shell=True)


    def _add_nullwt_paths(self, src_info):
        for sinfo in src_info:

            sinfo['nullwt_config'] = files.get_nwgint_config(self['campaign'])

            sinfo['nullwt_path'] = files.get_nullwt_file(
                self['medsconf'],
                self['tilename'],
                self['band'],
                sinfo['image_path'],
            )

    def _copy_psfs(self, info):

        medsdir=files.get_meds_base()

        psfmap_file=files.get_psfmap_file(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )

        psf_dir=files.get_psf_dir(self['medsconf'], self['tilename'])
        if not os.path.exists(psf_dir):
            print("making directory:",psf_dir)
            os.makedirs(psf_dir)

        print("writing psfmap:",psfmap_file)
        with open(psfmap_file,'w') as psfmap_fobj:
            print("copying psf files")

            psfs = self._get_psf_list(info)
            for psf_file in psfs:
                bname=os.path.basename(psf_file)
                ofile = os.path.join(psf_dir, bname)

                fs=bname.split('_')
                if 'DES' in fs[0]:
                    # this is the coadd psf, so fake it
                    expnum = -9999
                    ccdnum = -9999

                else:
                    # single epoch psf
                    expnum = fs[0][1:]
                    ccdnum = fs[2][1:]

#                ofile_medsdir=ofile.replace(medsdir, '$MEDS_DIR')
#                ofile_medsdir=ofile.replace(medsdir, self['MEDS_DIR'])
                ofile_medsdir = ofile
                psfmap_fobj.write("%s %s %s\n" % (expnum, ccdnum, ofile_medsdir))

                if os.path.exists(ofile):
                    continue

                print("copying: %s -> %s" % (psf_file, ofile))
                shutil.copy(psf_file, ofile)

    def _get_psf_list(self, info):
        psfs = []
        psfs.append(info['psf_path'])

        for sinfo in info['src_info']:
            psfs.append(sinfo['psf_path'])

        return psfs


_NULLWT_TEMPLATE=r"""
coadd_nwgint                  \
   -i "%(image_path)s"        \
   -o "%(nullwt_path)s"       \
   --headfile "%(head_path)s" \
   --max_cols 50              \
   -v                         \
   --interp_mask TRAIL,BPM    \
   --invalid_mask EDGE        \
   --null_mask BPM,BADAMP,EDGEBLEED,EDGE,CRAY,SSXTALK,STREAK,TRAIL \
   --me_wgt_keepmask STAR     \
   --block_size 5             \
   --tilename %(tilename)s    \
   --hdupcfg "%(nullwt_config)s"
"""


