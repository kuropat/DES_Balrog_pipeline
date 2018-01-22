#!/usr/bin/env python
from __future__ import print_function
import os
import sys
import meds
import fitsio
import numpy as np
import glob
from ..files import read_yaml

from .megamixer import NGMegaMixer

class BNLCondorMegaMixer(NGMegaMixer):
    def __init__(self,*args,**kwargs):
        super(BNLCondorMegaMixer,self).__init__(*args,**kwargs)
        self['queue'] = self.get('queue','medium')


    def setup_coadd_tile(self,coadd_tile):
        """
        Set up the condor script
        """

        # this makes the individual scripts
        files,fof_ranges=super(BNLCondorMegaMixer,self).setup_coadd_tile(coadd_tile)

        # above sets self['work_output_dir']
        condor_script=self.get_condor_file(files)
        #condor_script_all=self.get_condor_file(files,doall=True)

        _try_remove_file(condor_script)
        _try_remove_file(condor_script+'.submitted')
        #_try_remove_file(condor_script_all)
        #_try_remove_file(condor_script_all+'.submitted')

        self.write_master_script(files)

        # this holds everything and is not automatically used.
        # when using the megamixer run command, another file
        # is written just for those outputs that don't exist
        #self.write_condor(files,fof_ranges,doall=True)

        # this one just holds the jobs for which the output file
        # was not found.
        # note when using megamixit run this will be over-written
        # just in case some more files were completed or removed
        # after running setup
        self.write_condor(files,fof_ranges,doall=False)

    def run_nbrs_coadd_tile(self,coadd_tile):
        """
        run the neighbors creation
        """
        self.run_coadd_tile(coadd_tile, nbrs=True)

    def run_coadd_tile(self,coadd_tile, nbrs=False):
        """
        Write the condor file and submit it, only doing those
        for which the output file doesn't exist
        """

        if nbrs:
            files = self.get_files(coadd_tile)
            fname,nwrite=self.write_nbrs_condor(files)
        else:
            files,fof_ranges = self.get_files_fof_ranges(coadd_tile)
            fname,nwrite=self.write_condor(files,fof_ranges)

        if nwrite == 0:
            print("    no unfinished jobs left to run")
            return
        fname=os.path.basename(fname)

        dr = files['work_output_dir']
        cmd='cd %s && condor_submit %s && cd -' % (dr,fname)
        print("condor command: '%s'" % cmd)
        os.system(cmd)


    def setup_nbrs_coadd_tile(self,coadd_tile):
        """
        Set up the condor script
        """

        # this makes the individual scripts
        files=super(BNLCondorMegaMixer,self).setup_nbrs_coadd_tile(coadd_tile)

        # above sets self['work_output_dir']
        condor_script=self.get_nbrs_condor_file(files)
        condor_script_all=self.get_nbrs_condor_file(files,doall=True)

        _try_remove_file(condor_script)
        _try_remove_file(condor_script+'.submitted')
        _try_remove_file(condor_script_all)
        _try_remove_file(condor_script_all+'.submitted')

        # always written for posterity
        self.write_nbrs_condor(files,doall=True)

        # written only if the output file doesn't exist
        self.write_nbrs_condor(files,doall=False)


    def get_tmp_dir(self):
        return '${TMPDIR}'

    #def get_chunk_output_dir(self,files,chunk,rng):
    #    return os.path.join(files['work_output_dir'],
    #                        'chunk%06d' % (chunk+1))

    def get_chunk_jobname(self, files, rng):
        return '%s-%s-%05d-%05d' % (self['run'],files['coadd_tile'], rng[0], rng[1])

    def get_chunk_output_dir(self,files,chunk,rng):
        return os.path.join(files['work_output_dir'], \
                            'chunk%03d-%05d-%05d' % (chunk,rng[0],rng[1]))

    def get_chunk_output_basename(self,files,chunk,rng):
        return '%s-%s-%05d-%05d' % (files['coadd_tile'],self['run'],rng[0],rng[1])

    def get_script_template(self):
        """
        for condor we need to use the directory made by the system
        So we ignore tmpdir
        """
        template = r"""#!/bin/bash

output_dir={output_dir}
chunk={chunk}
config={ngmix_config}
obase={base_name}
output_file="$output_dir/$obase.fits"
meds="{meds_files}"

# ignoring tmpcmd
tmpdir=$TMPDIR

pushd $tmpdir

{cmd} \
    --fof-range={start},{stop} \
    --work-dir=$tmpdir \
    {fof_opt} \
    {nbrs_opt} \
    {flags_opt} \
    {seed_opt} \
    {mof_opt} \
    $config $output_file $meds

popd
        """
        return template


    def write_master_script(self,files):
        """
        this is a master script to run the individual scripts,
        move to the temporary directory made by condor, and
        move the log file to the output
        """
        path=self.get_master_script_file(files)
        text=self.get_master_script_text()

        print("    writing:",path)
        with open(path,'w') as fobj:
            fobj.write(text)
        os.system('chmod 755 %s' % path)

    def get_master_script_file(self, files):
        """
        location of the master script
        """
        path=os.path.join(files['work_output_dir'], 'master.sh')
        return path

    def get_nbrs_master_script_file(self, files):

        """
        location of the master script
        """
        path=os.path.join(files['work_output_dir'], 'master.sh')
        return path


    def get_master_script_text(self):
        """
        for condor we need to use the directory made by the system
        We cd there, write the log file and move to the final dir
        """
        template = r"""#!/bin/bash

if [[ -n $_CONDOR_SCRATCH_DIR ]]; then
    # the condor system creates a scratch directory for us,
    # and cleans up afterward
    tmpdir=$_CONDOR_SCRATCH_DIR
    export TMPDIR=$tmpdir
else
    # otherwise use the TMPDIR
    tmpdir=$TMPDIR
    mkdir -p $tmpdir
fi


# this will be a path to a script
cmd=$1
logfile=$2

# move the log file into that directory
odir=$(dirname $logfile)
tmplog=$(basename $logfile)

mkdir -p $odir

# any temp files written to cwd, plus log, go into the
# temporary directory
pushd $tmpdir

echo $cmd
$cmd &> $tmplog

mv -v "$tmplog" "$odir/"

popd
        \n"""
        return template



    def get_condor_file(self, files, doall=False):
        """
        need to have run get_files* to set the dir below

        will hold a line for every job
        """
        if doall:
            name='submit_all.condor'
        else:
            name='submit.condor'

        path=os.path.join(files['work_output_dir'], name)
        return path

    def get_nbrs_condor_file(self, files, doall=False):
        """
        need to have run get_files* to set the dir below

        will hold a line for every job
        """
        if doall:
            name='submit_all_nbrs.condor'
        else:
            name='submit_nbrs.condor'

        path=os.path.join(files['work_output_dir'], name)
        return path



    def get_condor_head(self,script):
        text="""
Universe        = vanilla

Notification    = Never

# Run this exe with these args
Executable      = {script}

#Image_Size      = 1000000
Image_Size       =  500000

GetEnv = True

kill_sig        = SIGINT

#requirements = (cpu_experiment == "star") || (cpu_experiment == "phenix")
#requirements = (cpu_experiment == "star")

+Experiment     = "astro"\n\n"""
        return text.format(script=script)

    def get_condor_job_lines(self, jobname, script_file, log_file):
        text="""
+job_name = "%(jobname)s"
Arguments = %(script_file)s %(log_file)s
Queue\n"""
        text=text % dict(
            jobname=jobname,
            script_file=script_file,
            log_file=log_file,
        )
        return text

    def get_nbrs_condor_job_lines(self, jobname):
        text="""
+job_name = "%(jobname)s"
Arguments = ""
Queue\n"""
        text=text % dict(
            jobname=jobname,
        )
        return text


    def write_condor(self,files,fof_ranges, doall=False):
        """
        write the condor submit file
        """

        fname=self.get_condor_file(files,doall=doall)
        master=self.get_master_script_file(files)

        nwrite=0
        with open(fname,'w') as fobj:
            head=self.get_condor_head(master)
            fobj.write(head)

            for chunk,rng in enumerate(fof_ranges):

                output_file=self.get_chunk_output_file(files,chunk, rng)

                if doall:
                    dowrite=True
                elif not os.path.exists(output_file):
                    dowrite=True
                else:
                    dowrite=False

                if dowrite:
                    nwrite+=1
                    jobname=self.get_chunk_jobname(files, rng)
                    script_file=self.get_chunk_script_file(files, chunk, rng)
                    log_file=self.get_chunk_log_file(files, chunk, rng)

                    lines=self.get_condor_job_lines(jobname, script_file, log_file)

                    fobj.write(lines)

        if doall:
            print("        total of",nwrite,"jobs",fname)
        else:
            print("        wrote",nwrite,"jobs",fname)
            if nwrite==0:
                for i in xrange(10):
                    _try_remove_file(fname)
                    _try_remove_file(fname+'.submitted')

        return fname, nwrite

    def write_nbrs_script(self,files):
        fmt=r"""#!/bin/bash

if [[ -n $_CONDOR_SCRATCH_DIR ]]; then
    # the condor system creates a scratch directory for us,
    # and cleans up afterward
    tmpdir=$_CONDOR_SCRATCH_DIR
    export TMPDIR=$tmpdir
else
    # otherwise use the TMPDIR
    tmpdir=$TMPDIR
    mkdir -p $tmpdir
fi

config=%(nbrs_config)s
logfile=%(logfile)s
meds_file=%(tile)s

# make sure output dir exists
odir=$(dirname $logfile)
mkdir -p $odir

# move the log file into that directory
tmplog=$(basename $logfile)

# any temp files written to cwd, plus log, go into the
# temporary directory
pushd $tmpdir

ngmixer-meds-make-nbrs-data ${config} ${meds_file} &> ${tmplog}

mv -v ${tmplog} ${logfile}
popd
"""
        args = {}
        args['nbrs_config'] = files['nbrs_config']
        args['logfile'] = files['nbrs_log']
        args['tile'] = files['meds_files'][0]

        scr = fmt % args

        script_file = self.get_nbrs_script_file(files)
        with open(script_file,'w') as fobj:
            fobj.write(scr)

        os.system('chmod 755 %s' % script_file)

    def write_nbrs_condor(self, files, doall=False):
        """
        write the condor submit file
        """

        nwrite=0
        output_file=files['nbrs_file']
        condor_file=self.get_nbrs_condor_file(files,doall=doall)
        script_file = self.get_nbrs_script_file(files)
        head=self.get_condor_head(script_file)

        if doall:
            dowrite=True
        elif not os.path.exists(output_file):
            dowrite=True
        else:
            dowrite=False

        if not dowrite:
            print("    output already exists")
            return condor_file, nwrite

        with open(condor_file,'w') as fobj:
            fobj.write(head)

            jobname=files['coadd_tile']
            lines=self.get_nbrs_condor_job_lines(jobname)
            fobj.write(lines)

        nwrite=1

        if doall:
            print("        total of 1 job",condor_file)
        else:
            print("        wrote 1 job",condor_file)

        return condor_file, nwrite

    def write_job_script(self,files,i,rng):
        """
        for condor we don't have a separate job script
        """
        return

    def write_nbrs_job_script(self, files):
        """
        for condor we don't have a separate job script
        """
        return


    def write_nbrs_job_script(self,files):
        """
        no individual job script in condor
        """
        return

    def run_nbrs(self,files):
        raise NotImplementedError("set up run nbrs")


def _try_remove_file(fname):
    try:
        os.remove(fname)
    except:
        pass


