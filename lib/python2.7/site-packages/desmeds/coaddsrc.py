from __future__ import print_function
import os
import shutil
import tempfile
import numpy
import fitsio

from . import files
from .coaddinfo import Coadd

class CoaddSrc(Coadd):
    """
    class to work with coadd sources (se images, etc.)
    """
    def __init__(self, *args, **kw):
        super(CoaddSrc,self).__init__(*args, **kw)
        self._set_finalcut_campaign()

    def get_info(self):
        """
        get info for the specified tilename and band
        """

        if hasattr(self,'_info_list'):
            info_list=self._info_list
        else:
            info_list = self._do_query()

            # add full path info
            self._add_full_paths(info_list)

            self._info_list=info_list

        return info_list

    def _do_query(self):
        """
        get info for the specified tilename and band
        """

        query = _QUERY_COADD_SRC_BYTILE % self

        print(query)
        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute(query)

        info_list=[]

        for row in curs:
            tile,path,fname,comp,band,pai,magzp = row
            info = {
                'tilename':tile,
                'filename':fname,
                'compression':comp,
                'path':path,
                'band':band,
                'pfw_attempt_id':pai,
                'magzp': magzp,
            }

            info_list.append(info)

        return info_list


    def _add_full_paths(self, info_list):
        """
        seg maps have .fz for finalcut
        """


        for info in info_list:

            dirdict=self._get_all_dirs(info)

            info['image_path'] = os.path.join(
                dirdict['image']['local_dir'],
                info['filename']+info['compression'],
            )

            info['bkg_path'] = os.path.join(
                dirdict['bkg']['local_dir'],
                info['filename'].replace('immasked.fits','bkg.fits')+info['compression'],
            )

            info['seg_path'] = os.path.join(
                dirdict['seg']['local_dir'],
                info['filename'].replace('immasked.fits','segmap.fits')+info['compression'],
            )

            info['psf_path'] = os.path.join(
                dirdict['psf']['local_dir'],
                info['filename'].replace('immasked.fits','psfexcat.psf')
            )


    def _get_all_dirs(self, info):
        dirs={}

        path=info['path']
        dirs['image'] = self._get_dirs(path)
        dirs['seg']   = self._get_dirs(path, type='seg')
        dirs['bkg']   = self._get_dirs(path, type='bkg')
        dirs['psf']   = self._get_dirs(path, type='psf')
        return dirs

    def _extract_alt_dir(self, path, type):
        """
        extract the catalog path from an image path, e.g.

        OPS/finalcut/Y2A1v3/20161124-r2747/D00596130/p01/red/immask/

        would yield

        OPS/finalcut/Y2A1v3/20161124-r2747/D00596130/p01/red/bkg/
        OPS/finalcut/Y2A1v3/20161124-r2747/D00596130/p01/seg

        """

        ps = path.split('/')

        assert ps[-1]=='immask'

        if type=='bkg':
            ps[-1] = type
        elif type in ['seg','psf']:
            ps = ps[0:-1]
            assert ps[-1]=='red'
            ps[-1] = type

        return '/'.join(ps)

    def _set_finalcut_campaign(self):
        if self['campaign']=='Y3A1_COADD':
            self['finalcut_campaign']='Y3A1_FINALCUT'
        elif self['campaign']=='Y3A2_COADD':
            self['finalcut_campaign']='Y3A1_FINALCUT'
        else:
            raise ValueError("determine finalcut campaign "
                             "for '%s'" % self['campaign'])


    def download(self, *args):
        raise NotImplementedError("use Coadd to download")
    def remove(self, *args):
        raise NotImplementedError("use Coadd to remove")


'''
class CoaddSrcCache(CoaddCache):
    """
    cache to hold path info for the sources of all
    coadds in the given campaign
    """
    def __init__(self, campaign='Y3A1_COADD'):
        self.campaign=campaign.upper()
        self._set_finalcut_campaign()

    def get_info(self, tilename, band):
        """
        get info for the specified tilename and band
        """

        query = _QUERY_COADD_SRC_BYTILE.format(
            campaign=self.campaign,
            finalcut_campaign=self.finalcut_campaign,
            tilename=tilename,
            band=band,
        )

        print(query)
        conn=self.get_conn()
        curs = conn.cursor()
        curs.execute(query)

        info_list=[]

        for row in curs:
            tile,path,fname,comp,band,pai,magzp = row
            info = {
                'tilename':tile,
                'filename':fname,
                'compression':comp,
                'path':path,
                'band':band,
                'pfw_attempt_id':pai,
                'magzp': magzp,
            }

            info_list.append(info)

        return info_list


    def get_info_old(self, tilename, band):
        """
        get info for the specified tilename and band
        """
        cache=self.get_data()

        key = make_cache_key(tilename, band)
        w,=numpy.where(cache['key']==key)

        if w.size == 0:
            raise ValueError("tilename '%s' and band '%s' "
                             "not found" % (tilename,band))

        entries=[]

        for i in w:
            c=cache[i]

            tilename=c['tilename'].strip()
            path=c['path'].strip()
            filename=c['filename'].strip()
            band=c['band'].strip()
            comp=c['compression'].strip()

            entry = {
                'tilename':tilename,
                'filename':filename,
                'compression':comp,
                'path':path,
                'band':band,
                'pfw_attempt_id':c['pfw_attempt_id'],

                # need to add this to the cache
                'magzp': 30.0,
            }

            entries.append(entry)


        return entries

    def get_data(self):
        """
        get the full cache data
        """
        if not hasattr(self,'_cache'):
            self.load_cache()
        return self._cache

    def load_cache(self):
        """
        load the cache into memory
        """

        fname=self.get_filename()
        zp_fname=self.get_zp_filename()
        if not os.path.join(fname):
            self.make_cache()

        print("loading cache:",fname)
        self._cache=fitsio.read(fname)
        self._cache['filename'] = \
                numpy.char.rstrip(self._cache['filename'])
        print("loading zp cache:",fname)
        self._zp_cache = fitsio.read(zp_filename)
        self._zp_cache['imagename'] = \
                numpy.char.rstrip(self._zp_cache['imagename'])

    def get_filename(self):
        """
        path to the cache
        """
        return files.get_coadd_src_cache_file(self.campaign)

    def get_zp_filename(self):
        """
        path to the cache
        """
        return files.get_zp_cache_file(self.campaign)


    def _get_query(self):
        query = _QUERY_COADD_SRC.format(
            campaign=self.campaign,
            finalcut_campaign=self.finalcut_campaign,
        )
        return query

    def make_zp_cache(self):
        """
        cache all the relevant information for this campaign
        """

        fname=self.get_zp_filename()

        print("writing to:",fname)
        query = _ZP_QUERY
        curs = self._doquery(query)

        dt=self._get_zp_dtype()

        info=numpy.fromiter(curs, dtype=dt)

        print("writing to:",fname)
        fitsio.write(fname, info, clobber=True)


    def _get_dtype(self):
        dt = super(CoaddSrcCache,self)._get_dtype()
        dt += [('magzp','f8')]
        return dt

    def _get_zp_dtype(self):
        return [
            ('imagename','S40'),
            ('magzp','f8'),
        ]

    def _set_finalcut_campaign(self):
        if self.campaign=='Y3A1_COADD':
            self.finalcut_campaign='Y3A1_FINALCUT'
        else:
            raise ValueError("determine finalcut campaign "
                             "for '%s'" % self.campaig)

'''
#select imagename, mag_zero from ZEROPOINT where IMAGENAME='D00504555_z_c41_r2378p01_immasked.fits' and source='FGCM' and version='v2.0';

#_QUERY_COADD_SRC="""
#select
#    i.tilename || '-' || j.band as key,
#    i.tilename,
#    fai.path,
#    j.filename as filename,
#    fai.compression,
#    j.band as band,
#    i.pfw_attempt_id,
#    z.mag_zero as magzp
#from
#    image i,
#    image j,
#    proctag tme,
#    proctag tse,
#    file_archive_info fai,
#    zeropoint z
#where
#    tme.tag='{campaign}'
#    and tme.pfw_attempt_id=i.pfw_attempt_id
#    and i.filetype='coadd_nwgint'
#    -- and i.tilename='DES0215-0458'
#    and i.expnum=j.expnum
#    and i.ccdnum=j.ccdnum
#    and j.filetype='red_immask'
#    and j.pfw_attempt_id=tse.pfw_attempt_id
#    and tse.tag='{finalcut_campaign}'
#    and fai.filename=j.filename
#    -- and z.imagename = j.filename
#    -- and z.source='FGCM'
#    -- and z.version='v2.0'
#    -- and rownum < 1000
#"""
_QUERY_COADD_SRC="""
select
    i.tilename || '-' || j.band as key,
    i.tilename,
    fai.path,
    j.filename as filename,
    fai.compression,
    j.band as band,
    i.pfw_attempt_id,
    z.mag_zero as magzp
from
    y3a2_image i,
    y3a2_image j,
    y3a2_proctag tme,
    y3a2_proctag tse,
    y3a2_file_archive_info fai,
    y3a2_zeropoint z
where
    tme.tag='{campaign}'
    and tme.pfw_attempt_id=i.pfw_attempt_id
    and i.filetype='coadd_nwgint'
    and i.tilename='DES0215-0458'
    and i.expnum=j.expnum
    and i.ccdnum=j.ccdnum
    and j.filetype='red_immask'
    and j.pfw_attempt_id=tse.pfw_attempt_id
    and tse.tag='{finalcut_campaign}'
    and fai.filename=j.filename
    and z.imagename = j.filename
    and z.source='FGCM'
    and z.version='v2.0'
    -- and rownum < 1000
"""
#_QUERY_COADD_SRC_BYTILE="""
#select
#    i.tilename,
#    fai.path,
#    j.filename as filename,
#    fai.compression,
#    j.band as band,
#    i.pfw_attempt_id,
#    z.mag_zero as magzp
#from
#    image i,
#    image j,
#    proctag tme,
#    proctag tse,
#    file_archive_info fai,
#    zeropoint z
#where
#    tme.tag='%(campaign)s'
#    and tme.pfw_attempt_id=i.pfw_attempt_id
#    and i.filetype='coadd_nwgint'
#    and i.tilename='%(tilename)s'
#    and i.expnum=j.expnum
#    and i.ccdnum=j.ccdnum
#    and j.filetype='red_immask'
#    and j.pfw_attempt_id=tse.pfw_attempt_id
#    and j.band='%(band)s'
#    and tse.tag='%(finalcut_campaign)s'
#    and fai.filename=j.filename
#    and z.imagename = j.filename
#    and z.source='FGCM'
#    and z.version='v2.0'
#order by
#    filename
#"""
_QUERY_COADD_SRC_BYTILE="""
select
    i.tilename,
    fai.path,
    j.filename as filename,
    fai.compression,
    j.band as band,
    i.pfw_attempt_id,
    z.mag_zero as magzp
from
    y3a2_image i,
    y3a2_image j,
    y3a2_proctag tme,
    y3a2_proctag tse,
    y3a2_file_archive_info fai,
    y3a2_zeropoint z
where
    tme.tag='%(campaign)s'
    and tme.pfw_attempt_id=i.pfw_attempt_id
    and i.filetype='coadd_nwgint'
    and i.tilename='%(tilename)s'
    and i.expnum=j.expnum
    and i.ccdnum=j.ccdnum
    and j.filetype='red_immask'
    and j.pfw_attempt_id=tse.pfw_attempt_id
    and j.band='%(band)s'
    and tse.tag='%(finalcut_campaign)s'
    and fai.filename=j.filename
    and z.imagename = j.filename
    and z.source='FGCM'
    and z.version='v2.0'
order by
    filename
"""

#_ZP_QUERY="""
#select
#    imagename,
#    mag_zero as magzp
#from
#    zeropoint
#where
#    source='FGCM'
#    and version='v2.0'
#"""
_ZP_QUERY="""
select
     imagename,
     mag_zero as magzp
from
     y3a2_zeropoint
where
     source='FGCM'
     and version='v2.0'
"""


_QUERY_COADD_SRC_old2="""
select
    i.tilename || '-' || j.band as key,
    i.tilename,
    fai.path,
    j.filename as filename,
    fai.compression,
    j.band as band,
    i.pfw_attempt_id
from
    image i,
    image j,
    proctag tme,
    proctag tse,
    file_archive_info fai 
where
    tme.tag='{campaign}'
    and tme.pfw_attempt_id=i.pfw_attempt_id
    and i.filetype='coadd_nwgint'
    -- and i.tilename='DES0215-0458'
    and i.expnum=j.expnum
    and i.ccdnum=j.ccdnum
    and j.filetype='red_immask'
    and j.pfw_attempt_id=tse.pfw_attempt_id
    and tse.tag='{finalcut_campaign}'
    and fai.filename=j.filename
    --and rownum < 1000
"""



