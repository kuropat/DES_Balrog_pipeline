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

# It is unlikely the user has cake installed
try:
    import cake
except:
    pass

class CakeNGMegaMixer(NGMegaMixer):
    def setup_coadd_tile(self,coadd_tile):
        print("setting up tile '%s'" % coadd_tile)
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)
        
        self._cake_chunks = []
        self._cake_chunk_ids = []
        
        self.make_output_dirs(files,fof_ranges)
        self.make_scripts(files,fof_ranges)

        # add to cake
        cake_file = os.path.join(files['work_output_dir'],'cake_chunks.db')
        try:
            os.remove(cake_file)
        except:
            pass
        with cake.taskdb.SQLiteTaskDB(name=cake_file) as taskdb:
            for cmd,id in zip(self._cake_chunks,self._cake_chunk_ids):
                taskdb.add(cmd,id=id)

        # add to global cake
        cake_file = os.path.join(files['run_output_dir'],'cake_chunks.db')
        with cake.taskdb.SQLiteTaskDB(name=cake_file) as taskdb:
            for cmd,id in zip(self._cake_chunks,self._cake_chunk_ids):
                try:
                    taskdb.add(cmd,id=id)
                except:
                    pass

    def write_job_script(self,files,i,rng):
        fname = os.path.join(self.get_chunk_output_dir(files,i,rng),'job.sh')
        jobname = self.get_chunk_output_basename(files,self['run'],rng)

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
        
        dr = self.get_chunk_output_dir(files,i,rng)
        self._cake_chunks.append('cd %s && ./job.sh && cd -' % dr)
        self._cake_chunk_ids.append(jobname)
        
    def run_chunk(self,files,chunk,rng):
        print("To run the chunks, build a job using cake (i.e., cake run /path/to/cake_chunks.db)")

    def write_nbrs_job_script(self,files):
        fname = os.path.join(files['work_output_dir'],'jobnbrs.sh')
        jobname = '%s-nbrs' % files['coadd_tile']

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

        # add to global cake
        cake_file = os.path.join(files['run_output_dir'],'cake_nbrs.db')
        dr = files['work_output_dir']
        with cake.taskdb.SQLiteTaskDB(name=cake_file) as taskdb:
            taskdb.add('cd %s && ./jobnbrs.sh && cd -' % dr,id=jobname)

    def run_nbrs(self,files):
        print("To run the chunks, build a job using cake (i.e., cake run /path/to/cake_nbrs.db)")
