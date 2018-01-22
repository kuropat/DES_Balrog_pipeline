#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import meds
import fitsio
import numpy as np
import glob
from ..files import read_yaml
from .concat_io import get_concat_class
from . import concat
import time

class BaseNGMegaMixer(dict):
    def __init__(self,conf,extra_cmds='',seed=None):
        self.update(conf)
        self.ngmix_conf = read_yaml(os.path.expandvars(self._get_ngmix_config()))

        # set this
        self['make_corrected_meds'] = self.ngmix_conf.get('make_corrected_meds',False)

        # do not set this to subtract neighbors
        self['model_nbrs'] = self.ngmix_conf.get('model_nbrs',False)

        if self['model_nbrs'] and self['make_corrected_meds']:
            raise ValueError("do not set both 'make_corrected_meds' and 'model_nbrs'")

        self['extra_cmds'] = extra_cmds
        self.rng = np.random.RandomState(seed=seed)

    def get_mof_file(self,full_coadd_tile,DESDATA,mof_version):
        coadd_run=full_coadd_tile.split('/')[-1]
        coadd_tile = full_coadd_tile.split('_')[-1]
        moff = os.path.join(DESDATA,
                            'EXTRA',
                            'meds',
                            self['meds_version'],
                            'mof-data',
                            mof_version,
                            coadd_run,
                            '%s-meds-%s-mof-%s.fits' % (coadd_tile,self['meds_version'],mof_version))

        return moff

    def get_files(self,full_coadd_tile):
        """
        get all paths to files
        """
        files = {}
        files['full_coadd_tile'] = full_coadd_tile
        coadd_tile = full_coadd_tile.split('_')[-1]
        files['coadd_tile'] = coadd_tile

        # desdata
        DESDATA = os.environ.get('DESDATA')
        files['DESDATA'] = DESDATA

        # meds files
        files['meds_files'] = []
        for band in self['bands']:
            medsf = os.path.join(DESDATA,
                                 'meds',
                                 self['meds_version'],
                                 files['full_coadd_tile'].split('/')[-1],
                                 '%s-%s-meds-%s.fits.fz' % (coadd_tile,band,self['meds_version']))
            assert os.path.exists(medsf),"MEDS file %s for band %s does not exist!" % (medsf,band)
            files['meds_files'].append(medsf)

        # now look for nbrs
        nbrsf = os.path.join(DESDATA,
                             'EXTRA',
                             'meds',
                             self['meds_version'],
                             'nbrs-data',
                             self['nbrs_version'],
                             files['full_coadd_tile'].split('/')[-1],
                             '%s-meds-%s-nbrslist-%s.fits' % (coadd_tile,self['meds_version'],self['nbrs_version']))
        files['nbrs_file'] = nbrsf
        files['nbrs_log'] = nbrsf.replace('.fits','.log')

        # now look for mof
        if 'mof_version' in self:
            moff = self.get_mof_file(full_coadd_tile,DESDATA,self['mof_version'])
        else:
            moff = None
        files['mof_file'] = moff

        # do the fofs
        foff = os.path.join(DESDATA,
                            'EXTRA',
                            'meds',
                            self['meds_version'],
                            'nbrs-data',
                            self['nbrs_version'],
                            files['full_coadd_tile'].split('/')[-1],
                            '%s-meds-%s-nbrsfofs-%s.fits' % (coadd_tile,self['meds_version'],self['nbrs_version']))
        files['fof_file'] = foff

        # finally look for flags
        flagsf = os.path.join(DESDATA,
                              'EXTRA',
                              'meds',
                              self['meds_version'],
                              'obj-flags-data',
                              self['obj_flags_version'],
                              files['full_coadd_tile'].split('/')[-1],
                              '%s-meds-%s-flags-%s.fits' % (coadd_tile,self['meds_version'],self['obj_flags_version']))
        files['obj_flags'] = flagsf

        # get output files
        self._set_output_files(files)

        # ngmix config
        files['ngmix_config'] = self._get_ngmix_config() 

        # nbrs config
        files['nbrs_config'] = self._get_nbrs_config()

        return files

    def _set_output_files(self,files):
        # main output dir
        if 'NGMIXER_OUTPUT_DIR' in os.environ:
            odir = os.path.expandvars(os.environ['NGMIXER_OUTPUT_DIR'])
        elif 'output_dir' in self:
            odir = self['output_dir']
        else:
            odir = '.'

        files['main_output_dir'] = os.path.join(odir,
                                                self['run'],
                                                files['full_coadd_tile'])

        files['work_output_dir'] = os.path.join(files['main_output_dir'],
                                                'work')

        files['run_output_dir'] = os.path.join(odir,self['run'])

    def _get_ngmix_config(self):
        if 'NGMIXER_CONFIG_DIR' in os.environ:
            ngmix_conf = os.path.join(os.environ['NGMIXER_CONFIG_DIR'],
                                      'ngmix_config',
                                      'ngmix-'+self['run']+'.yaml')
            ngmix_conf = os.path.expandvars(ngmix_conf)
            ngmix_conf = ngmix_conf.replace(os.environ['NGMIXER_CONFIG_DIR'],'${NGMIXER_CONFIG_DIR}')
        elif 'ngmix_config' in self:
            ngmix_conf = self['ngmix_config']
        else:
            ngmix_conf = 'ngmix-'+self['run']+'.yaml'

        return ngmix_conf

    def _get_nbrs_config(self):
        assert 'nbrs_config' in self or 'nbrs_version' in self,"nbrs_config or nbrs_version must be given to build nbrs files!"
        if 'nbrs_version' in self and 'NGMIXER_CONFIG_DIR' in os.environ:
            conf = os.path.join(os.environ['NGMIXER_CONFIG_DIR'],
                                'nbrs_config',
                                'nbrs-'+self['nbrs_version']+'.yaml')
            conf = os.path.expandvars(conf)
            conf = conf.replace(os.environ['NGMIXER_CONFIG_DIR'],'${NGMIXER_CONFIG_DIR}')
        elif 'nbrs_config' in self:
            conf = self['ngmix_config']
        else:
            conf = 'nbrs-'+self['nbrs_version']+'.yaml'

        return conf

    def setup_nbrs_coadd_tile(self,coadd_tile):
        print("setting up nbrs for tile '%s'" % coadd_tile)
        files = self.get_files(coadd_tile)
        self.make_output_dirs(files,[])
        os.system('cp %s %s' % (files['nbrs_config'],os.path.join(files['work_output_dir'],'.')))
        self.write_nbrs_script(files)
        self.write_nbrs_job_script(files)

        return files

    def get_nbrs_script_file(self, files):
        """
        script to run the neighbors generation code
        """
        return os.path.join(files['work_output_dir'],'runnbrs.sh')

    def write_nbrs_script(self,files):
        fmt="""#!/bin/bash
config={nbrs_config}
obase={base_name}
lfile=$obase".log"
meds_file={tile}

ngmixer-meds-make-nbrs-data $config $meds_file &> $lfile
"""
        args = {}
        args['nbrs_config'] = files['nbrs_config']
        args['base_name'] = '%s-nbrs' % files['coadd_tile']
        args['tile'] = files['meds_files'][0]
        args['cmd'] = 'ngmixer-meds-make-nbrs-data'

        scr = fmt.format(**args)

        scr_name = self.get_nbrs_script_file(files)
        with open(scr_name,'w') as fp:
            fp.write(scr)

        os.system('chmod 755 %s' % scr_name)

    def run_nbrs_coadd_tile(self,coadd_tile):
        print("running nbrs for tile '%s'" % coadd_tile)
        files = self.get_files(coadd_tile)
        self.run_nbrs(files)

    def get_chunk_output_dir(self,files,chunk,rng):
        return os.path.join(files['work_output_dir'], \
                            'chunk%05d_%d_%d' % (chunk,rng[0],rng[1]))

    def get_chunk_output_basename(self,files,chunk,rng):
        return '%s-%s-%d-%d' % (files['coadd_tile'],self['run'],rng[0],rng[1])

    def get_fof_ranges(self,files):
        if self['model_nbrs']:
            fofs = fitsio.read(files['fof_file'])
            num_fofs = len(np.unique(fofs['fofid']))
        else:
            m = meds.MEDS(files['meds_files'][0])
            num_fofs = m.size
            m.close()

        nchunks = num_fofs/self['num_fofs_per_chunk']
        if nchunks*self['num_fofs_per_chunk'] < num_fofs:
            nchunks += 1

        fof_ranges = []
        for chunk in xrange(nchunks):
            sr = chunk*self['num_fofs_per_chunk']
            sp = sr + self['num_fofs_per_chunk'] - 1
            if sp >= num_fofs:
                sp = num_fofs-1
            fof_ranges.append([sr,sp])

        return fof_ranges

    def make_output_dirs(self,files,fof_ranges):
        """
        make output dirs
        """
        odir = files['main_output_dir']
        wdir = files['work_output_dir']
        for dr in [odir,wdir]:
            if not os.path.exists(dr):
                os.makedirs(dr)

        os.system('cp %s %s' % (self['run_config'],os.path.join(wdir,'.')))

        for chunk,rng in enumerate(fof_ranges):
            dr = self.get_chunk_output_dir(files,chunk,rng)
            if not os.path.exists(dr):
                os.mkdir(dr)

        if len(self['extra_cmds']) > 0:
            os.system('cp %s %s' % (self['extra_cmds'],os.path.join(wdir,'.')))

    def get_tmp_dir(self):
        return '`mktemp -d /tmp/XXXXXXXXXX`'

    def make_scripts(self,files,fof_ranges):
        os.system('cp %s %s' % (files['ngmix_config'],os.path.join(files['work_output_dir'],'.')))

        for i,rng in enumerate(fof_ranges):
            self.write_script(files,i,rng)
            self.write_job_script(files,i,rng)

    def get_chunk_script_file(self, files, chunk, rng):
        """
        path to the script to run a chunk
        """
        return os.path.join(self.get_chunk_output_dir(files,chunk,rng),'runchunk.sh')

    def get_chunk_output_file(self, files, chunk, rng):
        """
        path to the script to run a chunk
        """
        dr = self.get_chunk_output_dir(files,chunk,rng)
        base = self.get_chunk_output_basename(files,self['run'],rng)
        fname = os.path.join(dr,base+'.fits')

        return fname

    def get_chunk_log_file(self, files, chunk, rng):
        """
        path to the script to run a chunk
        """
        outfile=self.get_chunk_output_file(files, chunk, rng)
        logfile=outfile.replace('.fits','.log')
        assert outfile!=logfile,"output file is not .fits"
        return logfile

    def get_script_template(self):
        template = r"""#!/bin/bash
chunk={chunk}
config={ngmix_config}
obase={base_name}
ofile=$obase".fits"
lfile=$obase".log"
meds="{meds_files}"
tmpdir={tmpcmd}


# call with -u to avoid buffering
if [ ! -f "$ofile" ]
then
    python -u `which {cmd}` \
        --fof-range={start},{stop} \
        --work-dir=$tmpdir \
        {fof_opt} \
        {nbrs_opt} \
        {flags_opt} \
        {seed_opt} \
        {mof_opt} \
        $config $ofile $meds &> $tmpdir/$lfile

    mv -v $tmpdir/$lfile .
fi
"""
        return template

    def write_script(self,files,chunk,rng):

        fmt=self.get_script_template()
        args = {}
        args['output_dir'] = self.get_chunk_output_dir(files,chunk,rng)
        args['chunk'] = chunk
        args['start'] = rng[0]
        args['stop'] = rng[1]
        if 'NGMIXER_CONFIG_DIR' not in os.environ:
            args['ngmix_config'] = os.path.join('..',files['ngmix_config'])
        else:
            args['ngmix_config'] = files['ngmix_config']
        args['meds_files'] = ' '.join([medsf.replace(files['DESDATA'],'${DESDATA}') for medsf in files['meds_files']])
        args['base_name'] = self.get_chunk_output_basename(files,self['run'],rng)
        args['tmpcmd'] = self.get_tmp_dir()
        args['cmd'] = 'ngmixit'


        args['fof_opt'] = ''
        args['nbrs_opt'] = ''
        args['flags_opt'] = ''

        if self['model_nbrs']:
            if os.path.exists(files['fof_file']):
                args['fof_opt'] = '--fof-file=%s'% files['fof_file'].replace(files['DESDATA'],'${DESDATA}')

            if os.path.exists(files['nbrs_file']):
                args['nbrs_opt'] = '--nbrs-file=%s'% files['nbrs_file'].replace(files['DESDATA'],'${DESDATA}')

            if os.path.exists(files['obj_flags']):
                args['flags_opt'] = '--obj-flags=%s'% files['obj_flags'].replace(files['DESDATA'],'${DESDATA}')

        if 'mof_version' in self:
            # for example will be used when 'make_corrected_meds' is set in the
            # ngmix config
            if os.path.exists(files['mof_file']):
                args['mof_opt'] = "--mof-file=%s"% files['mof_file'].replace(files['DESDATA'],'${DESDATA}')
        else:
            args['mof_opt'] = ''

        if 'seed' not in self:
            seed = self.rng.randint(low=1,high=1000000000)
            args['seed_opt'] = '--seed=%d' % seed
        else:
            args['seed_opt'] = ''

        scr = fmt.format(**args)

        scr_name = self.get_chunk_script_file(files,chunk,rng)
        with open(scr_name,'w') as fp:
            fp.write(scr)

        os.system('chmod 755 %s' % scr_name)

    def get_files_fof_ranges(self,coadd_tile):
        files = self.get_files(coadd_tile)

        # error check
        if self['model_nbrs']:
            assert os.path.exists(files['fof_file']),"fof file %s must be made to model nbrs!" % files['fof_file']
            assert os.path.exists(files['nbrs_file']),"nbrs file %s must be made to model nbrs!" % files['nbrs_file']

        if 'mof_version' in self:
            assert os.path.exists(files['mof_file']),"MOF file %s must be made to subtract nbrs on the fly!" % files['mof_file']

        fof_ranges = self.get_fof_ranges(files)

        return files,fof_ranges

    def setup_coadd_tile(self,coadd_tile):
        print("setting up tile '%s'" % coadd_tile)
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)

        self.make_output_dirs(files,fof_ranges)
        self.make_scripts(files,fof_ranges)

        return files, fof_ranges

    def run_coadd_tile(self,coadd_tile):
        print("running tile '%s'" % coadd_tile)
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)

        for chunk,rng in enumerate(fof_ranges):
            fname=self.get_chunk_output_file(files,chunk,rng)
            if not os.path.exists(fname):
                self.run_chunk(files,chunk,rng)

    def clean_coadd_tile(self,coadd_tile):
        print("cleaning tile '%s'" % coadd_tile)
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)

        for chunk,rng in enumerate(fof_ranges):
            dr = self.get_chunk_output_dir(files,chunk,rng)
            base = self.get_chunk_output_basename(files,self['run'],rng)

            fname = os.path.join(dr,base+'-checkpoint.fits')
            if os.path.exists(fname):
                os.remove(fname)

            fname = os.path.join(dr,base+'.fits')
            if os.path.exists(fname):
                os.remove(fname)

    def rerun_coadd_tile(self,coadd_tile):
        print("re-running tile '%s'" % coadd_tile)

        # first clean
        self.clean_coadd_tile(coadd_tile)

        # then run
        self.run_coadd_tile(coadd_tile)

    def get_collated_file(self, coadd_tile, blind=True):
        files = self.get_files(coadd_tile)
        collated_file = concat._get_collated_file(
            files['main_output_dir'],
            files['coadd_tile'],
            self['run'],
            blind=blind,
        )

        return collated_file, files

    def collate_coadd_tile(self,coadd_tile,verify=False,blind=True,clobber=True,skip_errors=False):


        collated_file, files=self.get_collated_file(coadd_tile, blind=blind)

        if not clobber and os.path.exists(collated_file):
            print('file',collated_file,'already exists, skipping')
            return

        # fof ranges is the slow part
        print("collating tile '%s'" % coadd_tile)
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)

        clist = []
        for chunk,rng in enumerate(fof_ranges):
            dr = self.get_chunk_output_dir(files,chunk,rng)
            base = self.get_chunk_output_basename(files,self['run'],rng)
            fname = os.path.join(dr,base+'.fits')
            clist.append(fname)

        tc = get_concat_class(self['concat_type'])
        tc = tc(self['run'],
                files['ngmix_config'],
                clist,
                files['main_output_dir'],
                files['coadd_tile'],
                bands=self['bands'],
                blind=blind,
                clobber=clobber,
                skip_errors=skip_errors)

        if verify:
            tc.verify()
        else:
            tc.concat()

    def archive_coadd_tile(self,coadd_tile,compress=True):
        print("archiving tile '%s'" % coadd_tile)

        # remove outputs
        self.clean_coadd_tile(coadd_tile)

        # tar (and maybe compress) the rest
        #files,fof_ranges = self.get_files_fof_ranges(coadd_tile)
        files = self.get_files(coadd_tile)
        if compress:
            tar_cmd = 'tar -czvf'
            tail = '.tar.gz'
        else:
            tar_cmd = 'tar -cvf'
            tail = '.tar'
        work_dir = files['work_output_dir']
        ofile = '%s%s' % (work_dir,tail)
        if os.path.exists(ofile):
            os.remove(ofile)
        cmd = '%s %s %s' % (tar_cmd,ofile,work_dir)
        os.system(cmd)

        # remove untarred work dir
        os.system('rm -rf %s' % work_dir)

    def link_coadd_tile(self,coadd_tile,verify=False,blind=True,clobber=True,skip_errors=False):
        print("linking tile '%s'" % coadd_tile)

        collated_file, files=self.get_collated_file(coadd_tile, blind=blind)

        bname = os.path.basename(collated_file)

        odir = '/'.join(files['main_output_dir'].split('/')[:-1])
        odir = os.path.join(odir,'output')
        if not os.path.exists(odir):
            os.makedirs(odir)

        cwd = os.getcwd()

        cfile = collated_file.split('/')
        cfile = os.path.join('..',cfile[-2],cfile[-1])

        try:
            os.chdir(odir)
            #os.system('rm -f %s && ln -s %s %s' % (bname,cfile,bname))
            os.system('rm -f %s' % bname)
            if os.path.exists(cfile):
                print("linking to:",os.path.join(odir,bname))
                os.system('ln -s %s %s' % (cfile,bname))
            else:
                print("missing")
            os.chdir(cwd)
        except:
            print("failed to link tile '%s'" % coadd_tile)

    def install_mof_outputs(self,coadd_tile,clobber=True):
        print("installing MOF outputs for tile '%s'" % coadd_tile)

        collated_file, files=self.get_collated_file(coadd_tile, blind=blind)

        '''
        # startup concat to get output name
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)

        clist = []
        for chunk,rng in enumerate(fof_ranges):
            dr = self.get_chunk_output_dir(files,chunk,rng)
            base = self.get_chunk_output_basename(files,self['run'],rng)
            fname = os.path.join(dr,base+'.fits')
            clist.append(fname)

        tc = get_concat_class(self['concat_type'])
        tc = tc(self['run'],
                files['ngmix_config'],
                clist,
                files['main_output_dir'],
                files['coadd_tile'],
                bands=self['bands'],
                blind=False,
                clobber=False,
                skip_errors=False)

        # now copy to proper spot in DESDATA
        collated_file = tc.collated_file
        '''

        moff = self.get_mof_file(coadd_tile,files['DESDATA'],self['run'])
        odir = os.path.split(moff)[0]
        if not os.path.exists(odir):
            os.makedirs(odir)

        if os.path.exists(moff):
            if clobber:
                try:
                    os.remove(moff)
                except:
                    pass
            else:
                raise IOError("MOF output file '%s' already exists in DESDATA!" % moff)

        print("    moving '%s' -> '%s'"%(collated_file,moff))
        os.system('cp %s %s' % (collated_file,moff))

    def write_job_script(self,files,i,rng):
        """
        method that writes a script to run the runchunk.sh script

        The script must run the extra cmds in the file in self['extra_cmds'],
        if this file exists, and then run 'runchunk.sh'.

        The script should assume it is in the same working dir as runchunk.sh.

        See the example below.

        """
        raise NotImplementedError("write_job_script method of BaseNGMegaMixer must be defined in subclass.")

    def run_chunk(self,files,chunk,rng):
        """
        This method must make some sort of system call to actually submit a single chunk to a queue
        or to run it on the local node.

        See the example below.
        """
        raise NotImplementedError("run_chunk method of BaseNGMegaMixer must be defined in subclass.")

    def write_nbrs_job_script(self,files):
        """
        method that writes a script to run the runnbrs.sh script

        The script must run the extra cmds in the file in self['extra_cmds'],
        if this file exists, and then run 'runnbrs.sh'.

        The script should assume it is in the same working dir as runnbrs.sh.

        See the example below.
        """
        raise NotImplementedError

    def run_nbrs(self,files):
        """
        This method must make some sort of system call to actually submit the nbrs job to a queue
        or to run it on the local node.

        See the example below.
        """
        raise NotImplementedError("run_nbrs method of BaseNGMegaMixer must be defined in subclass.")

class NGMegaMixer(BaseNGMegaMixer):
    def write_job_script(self,files,i,rng):
        fname = os.path.join(self.get_chunk_output_dir(files,i,rng),'job.sh')

        if len(self['extra_cmds']) > 0:
            with open(self['extra_cmds'],'r') as f:
                ec = f.readlines()
            ec = '\n'.join([e.strip() for e in ec])
        else:
            ec = ''

        with open(fname,'w') as fp:
            fp.write("""#!/bin/bash
{extracmds}
./runchunk.sh

""".format(extracmds=ec))

        os.system('chmod 755 %s' % fname)

    def run_chunk(self,files,chunk,rng):
        dr = self.get_chunk_output_dir(files,chunk,rng)
        os.system('cd %s && ./job.sh && cd -' % dr)

    def write_nbrs_job_script(self,files):
        fname = os.path.join(files['work_output_dir'],'jobnbrs.sh')

        if len(self['extra_cmds']) > 0:
            with open(self['extra_cmds'],'r') as f:
                ec = f.readlines()
            ec = '\n'.join([e.strip() for e in ec])
        else:
            ec = ''

        with open(fname,'w') as fp:
            fp.write("""#!/bin/bash
{extracmds}
./runnbrs.sh

""".format(extracmds=ec))

        os.system('chmod 755 %s' % fname)

    def run_nbrs(self,files):
        dr = files['work_output_dir']
        os.system('cd %s && ./jobnbrs.sh && cd -' % dr)
