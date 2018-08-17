"""
deprecated
"""
from __future__ import print_function
import sys
import os
import numpy
import fitsio

from . import files

try:
    import desdb
except ImportError:
    pass


_wq_template="""
command: |
    source ~/.bashrc
    %(cmd)s

job_name: %(job_name)s

mode: bynode
"""

_script_template="""#!/bin/bash
# we write stdout from make-meds-input to the meds_input file; log messages
# go to stderr
#
# make-cutouts only writes to stderr
#
# the exit status from those programs (the last to run) is written to the
# meds_status file and becomes the exit status of this program.
#
# there should be no stdout from this script, only stderr
#
# in production this script is called by a pbs submit script and the
# stderr is written to a log file (meds_log)

#nsetup_ess
module unload shapelets && module load shapelets/{vers}

magzp_ref={magzp_ref}
coadd_file="{coadd_image}"
coadd_image_id="{coadd_image_id}"
coaddseg_file="{coadd_seg}"
coadd_seg_hdu={coadd_seg_hdu}
coaddcat_file="{coadd_cat}"
coadd_objects_id_file="{coadd_objects_id_file}"
coadd_srclist="{meds_srclist}"
cutout_file="{meds_file}"
meds_input="{meds_input}"
medsconf={medsconf}
meds_status="{meds_status}"

min_boxsize={min_boxsize}
max_boxsize={max_boxsize}
use_alt_wcs={use_alt_wcs}
se_wcs_hdu={se_wcs_hdu}

# temporary directory and filenames
cutout_file_bname=$(basename $cutout_file)

db=$(echo $cutout_file_bname | sed 's/\.fits\.fz//')
tmpdir="$TMPDIR/meds-tmp-$db-$RANDOM"
mkdir -p "$tmpdir"

tmpfile_fz="$tmpdir/$cutout_file_bname"
tmpfile_uc=$(echo $tmpfile_fz | sed 's/\.fz//')

echo "tmp file fz: $tmpfile_fz"
echo "tmp file uc: $tmpfile_uc"


# make final output directory
output_dir=$(dirname "$cutout_file")
mkdir -p "$output_dir"

# remove existing file
if [[ -e "$cutout_file" ]]; then
    rm -v "$cutout_file"
fi

echo "creating meds input $meds_input"
make-meds-input "$coaddcat_file" $min_boxsize $max_boxsize $coadd_objects_id_file > "$meds_input"

exit_status=$?
if [[ $exit_status == "0" ]]; then

    # keywords after cutout_file= are just extra medatdata to be written into
    # the metadata table; they are not used by the code

    cmd="
    make-cutouts
        magzp_ref=$magzp_ref
        coadd_file=$coadd_file
        coadd_image_id=$coadd_image_id
        coaddseg_file=$coaddseg_file
        coadd_seg_hdu=$coadd_seg_hdu
        cat_file=$meds_input
        coadd_srclist=$coadd_srclist
        cutout_file=$tmpfile_uc
        medsconf=$medsconf
        coaddcat_file=$coaddcat_file
        min_boxsize=$min_boxsize
        max_boxsize=$max_boxsize
        use_alt_wcs=$use_alt_wcs
        se_wcs_hdu=$se_wcs_hdu
    "

    $cmd

    exit_status=$?
    if [[ $exit_status == "0" ]]; then
        size=$(du -h "$tmpfile_uc" | awk '{{print $1}}')
        echo "size of uncompressed file: $size"
        echo "fpacking to $tmpfile_fz"
        fpack -t {fpack_dim0},{fpack_dim1} "$tmpfile_uc" 

        ok="yes"
        exit_status=$?
        if [[ $exit_status == "0" ]]; then
            if [[ ! -e "$tmpfile_fz" ]]; then
                ok="no"
                echo "compressed cutout file not made $tmpfile_fz"
                echo "copying uncompressed file $tmpfile_uc"
                exit_status=1
            fi
        fi

        if [[ $ok == "yes" ]]; then
            cp -v "$tmpfile_fz" "$output_dir/"
        else
            cp -v "$tmpfile_uc" "$output_dir/"
        fi
    fi
fi

echo "cleaning up tmpdir: $tmpdir"
rm -rf $tmpdir

mess="writing status $exit_status to meds_status:
    $meds_status"
echo $mess

echo "exit_status: $exit_status" > "$meds_status"

echo "time: $SECONDS"
exit $exit_status
"""



def release_is_sva1(release):
    if isinstance(release,basestring):
        return 'sva1' in release.lower()
    else:
        for r in release:
            if 'sva1' in r.lower():
                return True

    return False

def get_magzp_offset(conf):
    """
    For SVA1 an artificial offset was added to the zeropoints
    in the zeropoint table.
    """

    if release_is_sva1(conf['release']):
        offset=0.057
        print('using magzp offset:',offset)
    else:
        offset=0.0

    return offset

class Generator(object):
    def __init__(self, medsconf, check=False, version='work'):

        self.medsconf=medsconf
        self.conf=files.read_meds_config(medsconf)
        self.conn=desdb.Connection()

        self.check=check
        self.version=version
        self.magzp_offset = get_magzp_offset(self.conf)

        self.df=desdb.files.DESFiles()

    def load_coadd(self, coadd_run, band):
        """
        load all the relevant info for the specified coadd and band
        """
        self.coadd_run=coadd_run
        self.band=band

        print("loading coadd info and srclist")
        self.cf=desdb.files.Coadd(coadd_run=coadd_run,
                                  band=band,
                                  conn=self.conn)

        self.cf.load(srclist=True)

        self.set_srclist()

    def write_all(self):
        """
        write all part
        """

        self.write_stats()

        if len(self.srclist) > 0:
            self.write_srclist()
            self.write_wq()
            self.write_idfile()
            self.write_script()
        else:
            print("not writing srclist/wq/idfile/script "
                  "because srclist is empty")

    def write_stats(self):
        """
        write stats about the inputs
        """
        import yaml
        stats_path = files.get_meds_stats_file(self.medsconf,
                                               self.coadd_run,
                                               self.cf['tilename'],
                                               self.band)

        make_dirs(stats_path)
        stats={'nsource':len(self.srclist)}
        with open(stats_path,'w') as fobj:
            yaml.dump(stats, stream=fobj)

    def write_srclist(self):
        """
        write the source list file
        """

        srclist_path = files.get_meds_srclist_file(self.medsconf,
                                                   self.coadd_run,
                                                   self.cf['tilename'],
                                                   self.band)

        make_dirs(srclist_path)

        print('writing:',srclist_path)
        nmissing=0
        with open(srclist_path,'w') as fobj:
            for r in self.srclist:

                if self.check:
                    nmissing += self.do_check_inputs(r)

                magzp = r['magzp'] + self.magzp_offset

                if self.include_wcs:
                    line='%d %d %s %s %s %s %s\n'
                    line = line % (r['id'],r['flags'],r['red_image'],r['red_bkg'],r['red_seg'],r['magzp'],r['wcs_file'])
                else:
                    line='%d %d %s %s %s %s\n'
                    line = line % (r['id'],r['flags'],r['red_image'],r['red_bkg'],r['red_seg'],r['magzp'])
                fobj.write(line)
        if self.check:
            print('nmissing: ',nmissing)

    def write_wq(self):
        """
        write the wq script
        """

        script_file = files.get_meds_script_file(self.medsconf,
                                                 self.cf['tilename'],
                                                 self.band)

        wq_file = files.get_meds_wq_file(self.medsconf,
                                         self.cf['tilename'],
                                         self.band)

        job_name='%s-%s' % (self.cf['tilename'],self.band)

        cmd="bash %s" % script_file

        text=_wq_template % {'job_name':job_name, 'cmd':cmd}

        make_dirs(wq_file)
        print('writing wq script:',wq_file)

        with open(wq_file,'w') as fobj:
            fobj.write(text)

    def get_coadd_cat_file(self, band):
        coadd_cat_file=self.df.url(type='coadd_cat',
                                   coadd_run=self.coadd_run,
                                   tilename=self.cf['tilename'],
                                   band=band)
        return coadd_cat_file

    def write_idfile(self):
        """
        generate the file with the coadd_objects_ids
        """
        detband=self.conf['detband']


        coadd_cat_file=self.get_coadd_cat_file(detband)

        # 1-column ascii file holding the coadd_objects_id
        cid_file = files.get_meds_coadd_objects_id_file(self.medsconf,
                                                        self.coadd_run,
                                                        self.cf['tilename'],
                                                        self.band)

        coadd_info = self.get_coadd_object_info()

        print("reading:",coadd_cat_file)
        coadd_cat=fitsio.read(coadd_cat_file, lower=True)

        print("verifying")
        verify_coadd_ids(coadd_info, coadd_cat)

        make_dirs(cid_file)

        print("writing:",cid_file)
        with open(cid_file,'w') as fobj:
            for i in xrange(coadd_info.size):
                fobj.write("%s\n" % coadd_info['coadd_objects_id'][i])


    def write_script(self):
        """
        write the script to create the meds file
        """

        df=self.df
        coadd_run=self.coadd_run
        band=self.band
        conf=self.conf
        cf=self.cf
        medsconf=self.medsconf

        if self.include_wcs:
            use_alt_wcs=1
            se_wcs_hdu=1
        else:
            use_alt_wcs=0
            se_wcs_hdu=2

        detband=conf['detband']

        if 'magzp_ref' in conf:
            magzp_ref = conf['magzp_ref']
            print('using fixed magzp_ref',magzp_ref)
        else:
            magzp_ref = desdb.files.get_release_magzp_ref(conf['release'], band)

        coadd_image=df.url(type='coadd_image',
                           coadd_run=coadd_run,
                           tilename=cf['tilename'],
                           band=band)

        coadd_seg=df.url(type='coadd_seg',
                         coadd_run=coadd_run,
                         tilename=cf['tilename'],
                         band=band)

        coadd_seg, is_fz =check_fz(coadd_seg)
        if is_fz:
            coadd_seg_hdu=2
        else:
            coadd_seg_hdu=1

        coadd_cat=df.url(type='coadd_cat',
                         coadd_run=coadd_run,
                         tilename=cf['tilename'],
                         band=detband)

        meds_srclist = files.get_meds_srclist_file(self.medsconf,
                                                   self.coadd_run,
                                                   self.cf['tilename'],
                                                   self.band)
        meds_input = files.get_meds_input_file(self.medsconf,
                                                 self.coadd_run,
                                                 self.cf['tilename'],
                                                 self.band)

        cid_file = files.get_meds_coadd_objects_id_file(self.medsconf,
                                                        self.coadd_run,
                                                        self.cf['tilename'],
                                                        self.band)
        meds_file = files.get_meds_file(self.medsconf,
                                        self.coadd_run,
                                        self.cf['tilename'],
                                        self.band)
        meds_status = files.get_meds_status_file(self.medsconf,
                                                 self.coadd_run,
                                                 self.cf['tilename'],
                                                 self.band)
        script_file = files.get_meds_script_file(self.medsconf,
                                                 self.cf['tilename'],
                                                 self.band)


        text=_script_template.format(medsconf=medsconf,
                                     magzp_ref=magzp_ref,
                                     vers=self.version,
                                     coadd_image=coadd_image,
                                     coadd_image_id=cf['image_id'],
                                     coadd_seg=coadd_seg,
                                     coadd_seg_hdu=coadd_seg_hdu,
                                     coadd_cat=coadd_cat,
                                     meds_input=meds_input,
                                     coadd_objects_id_file=cid_file,
                                     meds_srclist=meds_srclist,
                                     min_boxsize=conf['min_boxsize'],
                                     max_boxsize=conf['max_boxsize'],
                                     meds_file=meds_file,
                                     meds_status=meds_status,
                                     fpack_dim0=conf['fpack_dims'][0],
                                     fpack_dim1=conf['fpack_dims'][1],
                                     use_alt_wcs=use_alt_wcs,
                                     se_wcs_hdu=se_wcs_hdu)


        make_dirs(script_file, meds_file)

        print(script_file)
        with open(script_file,'w') as fobj:
            fobj.write(text)




    def get_coadd_object_info(self):
        query="""
    select
        object_number, coadd_objects_id
    from
        coadd_objects
    where
        imageid_%(band)s = %(coadd_id)d 
    order by
        object_number
        """ % {'coadd_id':self.cf['image_id'],
               'band':self.band}

        print("getting coadd_object_ids info for:",self.cf['image_id'])

        res = self.conn.quick(query, array=True)
        return res


    def set_srclist(self):
        """
        set the srclist, checking possibly for redone astrometry.
        also check against blacklist
        """
        srclist=self.cf.srclist

        for s in srclist:
            s['flags']=0

        add_bigind(srclist)

        self.include_wcs=False
        if 'astro_rerun_file' in self.conf:
            self.include_wcs=True
            self.wcs_type='astro_rerun'
            srclist=match_to_astro_rerun(srclist,
                                         self.conf,
                                         self.cf['tilename'])
        elif 'use_astro_refine' in self.conf:
            if self.conf['use_astro_refine']:
                self.include_wcs=True
                self.wcs_type='astro_refine'

                # set wcs_file in the srclist (the fits version of .head)
                # convert head files to fits files if needed
                set_astro_refine(self.cf['coadd_run'],srclist)

        add_blacklist_flags(srclist)

        self.srclist=srclist


    def do_check_inputs(self, r):
        nmissing=0
        for ftype in ['red_image','red_bkg','red_seg']:
            if not os.path.exists(r[ftype]):
                print("missing %s %s %s: %s" % (r['run'],r['expname'],ftype,r[ftype]))
                nmissing+=1

        return nmissing

def set_astro_refine(coadd_run, srclist):
    df=desdb.files.DESFiles()

    outdir=df.dir('astro_refine_fits',coadd_run=coadd_run)
    if not os.path.exists(outdir):
        try:
            print("making dir:",outdir)
            os.makedirs(outdir)
        except:
            pass

    for s in srclist:
        head_file=s['astro_refine']
        fits_file=df.url('astro_refine_fits',
                         coadd_run=coadd_run,
                         expname=s['expname'],
                         ccd=s['ccd'])
        s['wcs_file'] = fits_file

        if not os.path.exists(fits_file):

            print("reading:",head_file)
            hdata = fitsio.read_scamp_head(head_file)

            print("writing:",fits_file)
            fitsio.write(fits_file, None, header=hdata, clobber=True)


def match_to_astro_rerun(srclist, conf, tilename):
    """
    So the ASTROM_FLAG has the following bits set.  Good ones have ASTROM_FLAG == 0.  The flags mean:
    2^0: no Scamp solution (probably at boundary
    2^1: too few scamp stars (rare)
    2^2: no stats available: scamp crashed.  A few images at the boundary with this problem.
    2^3: too few (<50) good matched stars.  Just a problem for boundary images that aren't in gold anyway.
    2^4: internal use
    2^5: bad MAG_PSF_OFFSET.  > 4sigma outlier for the band in question.  Probably bad GCM solution; use with caution!
    2^6: bad RA offset: > 4sigma outlier for RA offset.  Probably bad astrometry solution.  use with caution!
    2^7: bad Dec offset
    2^8: bad RA rms (>30 mas for riz, 40mas for g)
    2^9: bad Dec rms

    You really want to filter out images with flags 2^0-2^3.  (And some of
    these are missing from my new header files!) You probably want to filter
    out images with the higher flags, but I leave this up to you.  These are
    probably bad images anyway.  And many of them are at the boundary.
    """
    import esutil as eu
    from esutil.numpy_util import ahelp
    from pprint import pprint

    fname=conf['astro_rerun_file']
    print("reading:",fname)
    t=fitsio.read(fname,lower=True)

    remove=2**0 + 2**1 + 2**2 + 2**3

    #w,=numpy.where( (t['astrom_flag'] & remove) == 0 )
    #if w.size == 0:
    #    raise RuntimeError("no good entries")

    #print("found %d/%d good entries" % (w.size, t.size))

    #bigind = make_bigind(t['expnum'][w], t['ccdnum'][w])
    bigind = make_bigind(t['expnum'], t['ccdnum'])

    mdict={}
    for i,ind in enumerate(bigind):
        #mdict[ind] = t[w[i]]
        mdict[ind] = t[i]

    new_srclist=[]
    for s in srclist:
        sbigind=s['bigind']
        mdata=mdict.get(sbigind,None)
        if mdata is not None:
            wcs_file=get_wcs_file_old(s)
            flags=mdata['astrom_flag']
            if (flags & remove) != 0:
                print("skipping bad:",s['expname'],s['ccd'])
            else:
                s['flags'] |= flags
                s['wcs_file'] = wcs_file
                new_srclist.append(s)
        else:
            print("error: not found:",sbigind)
            pprint(s)
            print("continuing anyway")
            #raise ValueError("image was not found!")

            #pprint(s)
            #ahelp(mdata)

    print("kept %d/%d from srclist" % (len(new_srclist),len(srclist)))
    return new_srclist


def read_blacklist(fname):
    import esutil as eu
    dt=[('expnum','i8'),
        ('ccd','i8')]

    #print("reading:",fname)
    with eu.recfile.Recfile(fname, 'r', delim=' ', dtype=dt) as fobj:
        data=fobj.read()
    return data

def read_blacklist_as_dict(fname):
    data=read_blacklist(fname)

    bigind=make_bigind(data['expnum'], data['ccd'])
    d={}
    for i in xrange(data.size):
        d[bigind[i]] = data[i]

    return d

def get_exp_blacklists():
    subdirs=['EXTRA','blacklists']
    dir=desdb.files.get_dir_generic(subdirs)

    # Eli's flags go to 2**9

    print("reading blackists")
    ldict={}
    fname=os.path.join(dir, 'ghost-scatter-sv-uniq.txt')
    ldict['ghost-sv'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**10}

    fname=os.path.join(dir, 'ghost-scatter-y1-uniq.txt')
    ldict['ghost-y1'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**11}

    fname=os.path.join(dir, 'noise-y1-uniq.txt')
    ldict['noise-y1'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**12}

    fname=os.path.join(dir, 'streak-sv-uniq.txt')
    ldict['streak-sv'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**13}

    fname=os.path.join(dir, 'streak-y1-uniq.txt')
    ldict['streak-y1'] = {'blacklist':read_blacklist_as_dict(fname),
                         'flag': 2**14}

    return ldict

def get_wcs_file_old(sdict):
    subdirs=['EXTRA',
             'red',
             sdict['run'],
             'astrorerun',
             sdict['expname']]

    ccdstr='%02d' % sdict['ccd']
    fileparts=[sdict['expname'], ccdstr, 'header']

    path = desdb.files.get_path_generic(subdirs, fileparts, ext='fits')
    return path

def add_blacklist_flags(srclist):
    """
    bigind and flags must be present already
    """
    blacklists = get_exp_blacklists()
    for s in srclist:
        for bname,blist in blacklists.iteritems():
            if s['bigind'] in blist['blacklist']:
                print("    found in blacklist:",bname)
                s['flags'] |= blist['flag']


def add_bigind(srclist):
    for s in srclist:
        expname=s['expname']
        expnum=int( expname.split('_')[1] )
        s['bigind'] = make_bigind(expnum, s['ccd'])


def make_bigind(expnum, ccdnum):
    return expnum + ccdnum*10**7


def verify_coadd_ids(coadd_info, coadd_cat):

    w,=numpy.where(coadd_cat['number'] != coadd_info['object_number'])
    if w.size > 0:
        raise ValueError("number fields don't "
                         "match %d/%d" % (w.size,cat.size))


def make_dirs(*args):

    for f in args:
        d=os.path.dirname(f)
        if not os.path.exists(d):
            print('making dir:',d)
            try:
                os.makedirs(d)
            except:
                pass

def check_fz(name):
    if not os.path.exists(name):
        is_fz=False
        name=name.replace('.fits.fz','.fits')
        if not os.path.exists(name):
            raise ValueError("fits file missing: %s{.fz}" % name)
    else:
        is_fz=True

    return name, is_fz


