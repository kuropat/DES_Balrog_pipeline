from __future__ import print_function
import os
import shutil
import numpy
import fitsio
import tempfile
import subprocess
import easyaccess
from . import files

class Coadd(dict):
    """
    information for coadds.  Can use the download() method to copy
    to the local disk heirchy
    """
    def __init__(self, medsconf, tilename, band, campaign='Y3A2_COADD', src=None, sources=None):
        dbname = 'dessci'
        desfile = os.getenv("DES_SERVICES")
        if not desfile: desfile = os.path.join(os.getenv("HOME"), ".desservices.ini")
#
        self.desconfig = easyaccess.config_ea.get_desconfig(desfile, dbname)
        self['user'] = self.desconfig.get('db-' + dbname, 'user')
        self['password'] = self.desconfig.get('db-' + dbname, 'passwd')

        
        self['medsconf'] = medsconf
        self['tilename'] = tilename
        self['band'] = band

        self['source_dir']=files.get_source_dir(
            self['medsconf'],
            self['tilename'],
            self['band'],
        )

        self['campaign']=campaign.upper()
        self.sources=sources

    def get_info(self):
        """
        get info for the tilename and band

        if sources were sent to the constructor, add source info as well
        """

        if hasattr(self, '_info'):
            info = self._info
        else:
            info = self._do_query()

            # add full path info
            self._add_full_paths(info)

            sources=self.get_sources()
            if sources is not None:
                self._add_src_info(info)

            self._info=info

        return info

    def download(self):
        """
        download sources for a single tile and band
        """
        if not os.path.exists(self['source_dir']):
            print("making source dir:",self['source_dir'])
            os.makedirs(self['source_dir'])

        info=self.get_info()

        self['flist_file']=self._write_download_flist(info)

        if 'DESREMOTE_RSYNC_USER' in os.environ:
            self['userstring'] = os.environ['DESREMOTE_RSYNC_USER']+'@'
        else:
            self['userstring'] = ''

        cmd=_DOWNLOAD_CMD % self
        print("command :",cmd)

        try:
            subprocess.check_call(cmd,shell=True)
        finally:
            files.try_remove(self['flist_file'])

        return info
        
    def clean(self):
        """
        remove downloaded files for the specified tile and band
        """

        print("removing sources:",self['source_dir'])
        shutil.rmtree(self['source_dir'])

    def _get_objmap_query(self, info):
        #return _OBJECT_MAP_QUERY
        filename=os.path.basename(info['cat_path'])
        #filename=os.path.basename(info['filename'])
        return _OBJECT_MAP_QUERY % filename

    def _get_objmap_dtype(self):
        return [ 
            ('object_number','i4'),
            ('id','i8'),
        ]


    def get_sources(self):
        """
        get the source list
        """
        return self.sources

    def _do_query(self):
        """
        get info for the specified tilename and band
        """
        
        query = _QUERY_COADD_TEMPLATE_BYTILE % self

        print(query)
        conn=self.get_conn()
        curs = conn.cursor()
        curs.execute(query)

        c=curs.fetchall()
#        print(str(c), sep=' ', end='\n',)
        tile,path,fname,comp,band,pai = c[0]

        entry = {
            'tilename':tile,
            'filename':fname,
            'compression':comp,
            'path':path,
            'band':band,
            'pfw_attempt_id':pai,

            # need to add this to the cache?  should always
            # be the same...
            'magzp': 30.0,
        }

        return entry



    def _add_full_paths(self, info):
        """
        seg maps don't have .fz extension for coadd
        """
        dirdict=self._get_all_dirs(info)
        info['image_path'] = os.path.join(
            dirdict['image']['local_dir'],
            info['filename']+info['compression'],
        )
        info['cat_path'] = os.path.join(
            dirdict['cat']['local_dir'],
            info['filename'].replace('.fits','_cat.fits'),
        )
        info['seg_path'] = os.path.join(
            dirdict['seg']['local_dir'],
            info['filename'].replace('.fits','_segmap.fits'),
        )
        info['psf_path'] = os.path.join(
            dirdict['psf']['local_dir'],
            info['filename'].replace('.fits','_psfcat.psf'),
        )


    def _get_download_flist(self, info, no_prefix=False):
        """
        get list of files for this tile

        parameters
        ----------
        info: dict
            The info dict for this tile/band, possibly including
            the src_info

        no_prefix: bool 
            If True, the {source_dir} is removed from the front
        """
        source_dir=self['source_dir']

        if source_dir[-1] != '/':
            source_dir += '/'

        types=self._get_download_types()
        stypes=self._get_source_download_types()

        flist=[]
        for type in types:
            tname='%s_path' % type

            fname = info[tname]

            if no_prefix:
                fname = fname.replace(source_dir, '')

            flist.append(fname)

        if 'src_info' in info:
            for sinfo in info['src_info']:
                for type in stypes:
                    tname='%s_path' % type

                    fname = sinfo[tname]

                    if no_prefix:
                        fname = fname.replace(source_dir, '')

                    flist.append(fname)


        return flist


    def _write_download_flist(self, info):

        flist_file=self._get_tempfile()
        flist=self._get_download_flist(info, no_prefix=True)
#        print ( flist, sep=' ', end='\n',)

        print("writing file list to:",flist_file)
        with open(flist_file,'w') as fobj:
            for fname in flist:
                fobj.write(fname)
                fobj.write('\n')
                
        return flist_file

    def _get_tempfile(self):
        return tempfile.mktemp(
            prefix='coadd-flist-',
            suffix='.dat',
        )


    def _get_download_types(self):
        return ['image','cat','seg','psf']

    def _get_source_download_types(self):
        return ['image','bkg','seg','psf','head']


    def _add_src_info(self, info):
        """
        get path info for the input single-epoch sources
        """

        sources=self.get_sources()
        src_info = self.sources.get_info()

        self._add_head_full_paths(info, src_info)

        info['src_info'] = src_info

    def _add_head_full_paths(self, info, src_info):
        dirdict=self._get_all_dirs(info)

        # this is a full path
        auxdir=dirdict['aux']['local_dir']

        head_front=info['filename'][0:-7]

        for src in src_info:
            fname=src['filename']

            fid = fname[0:15]

            head_fname = '%s_%s_scamp.ohead' % (head_front, fid)

            src['head_path'] = os.path.join(
                auxdir,
                head_fname,
            )

    def get_conn(self):
        if not hasattr(self, '_conn'):
            self._make_conn()

        return self._conn

    def _make_conn(self):
        sources = self.get_sources()
        if sources is not None:
            # share connection with the sources
            conn=sources.get_conn()
        else:
            import easyaccess as ea
            conn=ea.connect(section='dessci')

        self._conn=conn

    def _get_all_dirs(self, info):
        dirs={}

        path=info['path']
        dirs['image'] = self._get_dirs(path)
        dirs['cat'] = self._get_dirs(path, type='cat')
        dirs['aux'] = self._get_dirs(path, type='aux')
        dirs['seg'] = self._get_dirs(path, type='seg')
        dirs['psf'] = self._get_dirs(path, type='psf')
        return dirs

    def _get_dirs(self, path, type=None):
        #local_dir = '$DESDATA/%s' % path
        local_dir = '%s/%s' % (self['source_dir'], path)
        remote_dir = '$DESREMOTE_RSYNC/%s' % path

        local_dir=os.path.expandvars(local_dir)
        remote_dir=os.path.expandvars(remote_dir)

        if type is not None:
            local_dir=self._extract_alt_dir(local_dir,type)
            remote_dir=self._extract_alt_dir(remote_dir,type)

        return {
            'local_dir':local_dir,
            'remote_dir':remote_dir,
        }

    def _extract_alt_dir(self, path, type):
        """
        extract the catalog path from an image path, e.g.

        OPS/multiepoch/Y3A1/r2577/DES0215-0458/p01/coadd/

        would yield

        OPS/multiepoch/Y3A1/r2577/DES0215-0458/p01/{type}/
        """

        ps = path.split('/')

        assert ps[-1]=='coadd'

        ps[-1] = type
        return '/'.join(ps)
    


'''
class CoaddCache(object):
    """
    cache to hold path info for the coadds
    """
    def __init__(self, campaign='Y3A1_COADD'):
        self.campaign=campaign.upper()

    def get_data(self):
        """
        get the full cache data
        """
        if not hasattr(self,'_cache'):
            self.load_cache()

        return self._cache

    def get_info(self, tilename, band):
        """
        get info for the specified tilename and band
        """
        
        query = _QUERY_COADD_TEMPLATE_BYTILE.format(
            campaign=self.campaign,
            tilename=tilename,
            band=band,
        )

        print(query)
        conn=self.get_conn()
        curs = conn.cursor()
        curs.execute(query)

        c=curs.fetchall()

        tile,path,fname,comp,band,pai = c[0]

        entry = {
            'tilename':tile,
            'filename':fname,
            'compression':comp,
            'path':path,
            'band':band,
            'pfw_attempt_id':pai,

            # need to add this to the cache?  should always
            # be the same...
            'magzp': 30.0,
        }

        return entry

    def get_info_old(self, tilename, band='i'):
        """
        get info for the specified tilename and band
        """

        cache = self.get_data()

        key = make_cache_key(tilename, band)
        w,=numpy.where(cache['key'] == key)
        if w.size == 0:
            raise ValueError("%s not found in cache" % key)

        c=cache[w[0]]

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

            # need to add this to the cache?  should always
            # be the same...
            'magzp': 30.0,
        }

        return entry


    def get_objmap(self, info):
        """
        get the mapping between OBJECT_NUMBER and ID
        """
        query=self._get_objmap_query(info)
        print(query)

        conn = self.get_conn()
        curs = conn.cursor()
        curs.execute(query)

        dtype=self._get_objmap_dtype()
        return numpy.fromiter(curs,dtype=dtype)

    def _get_objmap_query(self, info):
        #return _OBJECT_MAP_QUERY
        filename=os.path.basename(info['cat_path'])
        #filename=os.path.basename(info['filename'])
        return _OBJECT_MAP_QUERY % filename

    def load_cache(self):
        """
        load the cache into memory
        """

        fname=self.get_filename()
        if not os.path.join(fname):
            self.make_cache()

        print("loading cache:",fname)
        self._cache = fitsio.read(fname)

    def make_cache(self):
        """
        cache all the relevant information for this campaign
        """

        fname=self.get_filename()

        print("writing to:",fname)
        query = self._get_query()
        curs = self._doquery(query)

        dt=self._get_dtype()

        info=numpy.fromiter(curs, dtype=dt)


        print("writing to:",fname)
        fitsio.write(fname, info, clobber=True)

    def _get_dtype(self):
        return [ 
            ('key','S14'),
            ('tilename','S12'),
            ('path','S65'),
            ('filename','S40'),
            ('compression','S3'),
            ('band','S1'),
            ('pfw_attempt_id','i8'),

        ]

    def _get_objmap_dtype(self):
        return [ 
            ('object_number','i4'),
            ('id','i8'),
        ]


    def get_filename(self):
        """
        path to the cache
        """
        return files.get_coadd_cache_file(self.campaign)

    def _doquery(self, query):

        print(query)
        conn=self.get_conn()
        curs = conn.cursor()
        curs.execute(query)

        return curs

    def _get_query(self):
        query = _QUERY_COADD_TEMPLATE.format(
            campaign=self.campaign,
        )
        return query

    def get_conn(self):
        if not hasattr(self, '_conn'):
            self._make_conn()

        return self._conn

    def _make_conn(self):
        import easyaccess as ea
        self._conn=ea.connect(section='dessci')


def make_cache_key(tilename, band):
    return '%s-%s' % (tilename, band)

'''


#_QUERY_COADD_TEMPLATE="""
#select
#    m.tilename || '-' || m.band as key,
#    m.tilename as tilename,
#    fai.path as path,
#    fai.filename as filename,
#    fai.compression as compression,
#    m.band as band,
#    m.pfw_attempt_id as pfw_attempt_id

#from
#    prod.proctag t,
#    prod.coadd m,
#    prod.file_archive_info fai
#where
#    t.tag='{campaign}'
#    and t.pfw_attempt_id=m.pfw_attempt_id
#    and m.filetype='coadd'
#    and fai.filename=m.filename
#    and fai.archive_name='desar2home'\n"""

_QUERY_COADD_TEMPLATE="""

select
     m.tilename || '-' || m.band as key,
     m.tilename as tilename,
     fai.path as path,
     fai.filename as filename,
     fai.compression as compression,
     m.band as band,
     m.pfw_attempt_id as pfw_attempt_id

from
     y3a2_proctag t,
     y3a2_coadd m,
     y3a2_file_archive_info fai
where
     t.tag='{campaign}'
     and t.pfw_attempt_id=m.pfw_attempt_id
     and m.filetype='coadd'
     and fai.filename=m.filename
     and fai.archive_name='desar2home' """

#_QUERY_COADD_TEMPLATE_BYTILE="""
#select
#    m.tilename as tilename,
#    fai.path as path,
#    fai.filename as filename,
#    fai.compression as compression,
#    m.band as band,
#    m.pfw_attempt_id as pfw_attempt_id

#from
#    prod.proctag t,
#    prod.coadd m,
#    prod.file_archive_info fai
#where
#    t.tag='%(campaign)s'
#    and t.pfw_attempt_id=m.pfw_attempt_id
#    and m.tilename='%(tilename)s'
#    and m.band='%(band)s'
#    and m.filetype='coadd'
#    and fai.filename=m.filename
#    and fai.archive_name='desar2home'\n"""

_QUERY_COADD_TEMPLATE_BYTILE="""
select
    m.tilename as tilename,
    fai.path as path,
    fai.filename as filename,
    fai.compression as compression,
    m.band as band,
    m.pfw_attempt_id as pfw_attempt_id
from
    y3a2_proctag t,
    y3a2_coadd m,
    y3a2_file_archive_info fai
where
    t.tag='%(campaign)s'
    and t.pfw_attempt_id=m.pfw_attempt_id
    and m.tilename='%(tilename)s'
    and m.band='%(band)s'
    and m.filetype='coadd'
    and fai.filename=m.filename
    and fai.archive_name='desar2home' """

_DOWNLOAD_CMD = r"""
    wget \
        --http-user=%(user)s --http-password=%(password)s --no-check-certificate \
        -i %(flist_file)s -P %(source_dir)s -x -nH --cut-dirs=2 \
        --base https://desar2.cosmology.illinois.edu/DESFiles/desarchive/ \
"""

_OBJECT_MAP_QUERY = """
select
    object_number,
    id
from
    -- coadd_object
    -- prod.COADD_OBJECT_SAVE
    prod.COADD_OBJECT
where
    filename='%s'
order by
    object_number
"""

#
# not used
#
_QUERY_TEMPLATE="""
select
    fai.path as path,
    fai.filename as filename,
    fai.compression as compression,
    m.band as band,
    m.pfw_attempt_id as pfw_attempt_id

from
    prod.proctag t,
    prod.coadd m,
    prod.file_archive_info fai
where
    t.tag='{campaign}'
    and t.pfw_attempt_id=m.pfw_attempt_id
    and m.tilename='{tilename}'
    and m.filetype='coadd'
    and fai.filename=m.filename
    and fai.archive_name='desar2home' and rownum <= 10\n"""
