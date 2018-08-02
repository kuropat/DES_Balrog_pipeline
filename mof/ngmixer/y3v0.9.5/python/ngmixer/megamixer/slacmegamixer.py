#!/usr/bin/env python
import os
import sys
import meds
import fitsio
import numpy as np
import glob
from ..files import read_yaml

from .megamixer import NGMegaMixer

class SLACNGMegaMixer(NGMegaMixer):
    def __init__(self,*args,**kwargs):
        super(SLACNGMegaMixer,self).__init__(*args,**kwargs)
        self['queue'] = self.get('queue','medium')

    def get_tmp_dir(self):
        return '${TMPDIR}'

    def write_job_script(self,files,chunk,rng):
        fname = os.path.join(self.get_chunk_output_dir(files,chunk,rng),'job.sh')
        if os.path.exists(fname):
            os.remove(fname)
        subf=fname+'.submitted'
        if os.path.exists(subf):
            os.remove(subf)

        output_file=self.get_chunk_output_file(files,chunk, rng)

        if not os.path.exists(output_file):

            jobname = self.get_chunk_output_basename(files,self['run'],rng)

            if len(self['extra_cmds']) > 0:
                with open(self['extra_cmds'],'r') as f:
                    ec = f.readlines()
                ec = '\n'.join([e.strip() for e in ec])
            else:
                ec = ''

            with open(fname,'w') as fp:
                fp.write("""#!/bin/bash
#BSUB -J {jobname}
#BSUB -oo ./{jobname}.oe
#BSUB -R "linux64 && rhel60 && scratch > 6"
#BSUB -n 1
#BSUB -We 24:00
#BSUB -W 48:00

{extracmds}
export TMPDIR=/scratch/$USER/$LSB_JOBID-$LSB_JOBINDEX
mkdir -p $TMPDIR
./runchunk.sh
rm -rf $TMPDIR

""".format(extracmds=ec,jobname=jobname))

            os.system('chmod 755 %s' % fname)

    def run_chunk(self,files,chunk,rng):
        dr = self.get_chunk_output_dir(files,chunk,rng)
        os.system('cd %s && bsub -q %s < job.sh && cd -' % (dr,self['queue']))

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
#BSUB -J {jobname}
#BSUB -oo ./{jobname}.oe
#BSUB -R "linux64 && rhel60 && scratch > 6"
#BSUB -n 1
#BSUB -We 24:00
#BSUB -W 48:00

{extracmds}
export TMPDIR=/scratch/$USER/$LSB_JOBID
mkdir -p $TMPDIR
./runnbrs.sh
rm -rf $TMPDIR

""".format(extracmds=ec,jobname=jobname))

        os.system('chmod 755 %s' % fname)

    def run_nbrs(self,files):
        dr = files['work_output_dir']
        os.system('cd %s && bsub -q %s < jobnbrs.sh && cd -' % (dr,self['queue']))

class SLACArrayNGMegaMixer(SLACNGMegaMixer):
    def __init__(self,*args,**kwargs):
        super(SLACNGMegaMixer,self).__init__(*args,**kwargs)
        self['queue'] = self.get('queue','medium')

    def get_chunk_output_dir(self,files,chunk,rng):
        return os.path.join(files['work_output_dir'],
                            'chunk%d' % (chunk+1))

    def write_job_script(self,files,i,rng):
        super(SLACArrayNGMegaMixer,self).write_job_script(files,i,rng)

        array_name = os.path.join(files['work_output_dir'],'array_jobs.dat')
        dr = os.path.join('.',*self.get_chunk_output_dir(files,i,rng).split('/')[-1:])
        with open(array_name,'a') as fp:
            fp.write('cd %s && ./job.sh && cd -\n' % dr)

    def run_chunk(self,files,chunk,rng):
        self._run_chunks.append(chunk)

    def setup_coadd_tile(self,coadd_tile):
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)
        try:
            os.remove(os.path.join(files['work_output_dir'],'array_jobs.dat'))
        except:
            pass

        super(SLACArrayNGMegaMixer,self).setup_coadd_tile(coadd_tile)

        raname = os.path.join(files['work_output_dir'],'runarray.py')
        with open(raname,'w') as fp:
            fp.write("#!/usr/bin/env python\n"
                     "import os\n"
                     "import sys\n"
                     "\n"
                     "chunk = int(sys.argv[1])-1\n"
                     "\n"
                     "loc = 0\n"
                     "for line in open('array_jobs.dat','r'):\n"
                     "    cmd = line.strip()\n"
                     "    if loc == chunk:\n"
                     "        break\n"
                     "    loc += 1\n"
                     "\n"
                     "os.system(cmd)\n"
                     "\n")

        os.system('chmod 755 %s' % raname)

    def run_coadd_tile(self,coadd_tile):
        self._run_chunks = []
        super(SLACArrayNGMegaMixer,self).run_coadd_tile(coadd_tile)
        files,fof_ranges = self.get_files_fof_ranges(coadd_tile)

        if len(self['extra_cmds']) > 0:
            with open(self['extra_cmds'],'r') as f:
                ec = f.readlines()
            ec = '\n'.join([e.strip() for e in ec])
        else:
            ec = ''

        raname = os.path.join(files['work_output_dir'],'jobarray.sh')
        jobname = files['coadd_tile']

        if len(self._run_chunks) > 0:
            rngs = []
            curr_rng = [self._run_chunks[0],self._run_chunks[0]]
            for i in xrange(1,len(self._run_chunks)):
                if self._run_chunks[i] == self._run_chunks[i-1]+1:
                    curr_rng[1] = self._run_chunks[i]
                else:
                    curr_rng[1] = self._run_chunks[i-1]
                    rngs.append(curr_rng)
                    curr_rng = [self._run_chunks[i],self._run_chunks[i]]

            rngs.append(curr_rng)

            istr = '['
            for rng in rngs:
                if rng[0] == rng[1]:
                    istr += '%d,' % (rng[0]+1)
                else:
                    istr += '%d-%d,' % (rng[0]+1,rng[1]+1)
            istr = istr[:-1]
            istr += ']'
            jobnamearr = jobname+istr

            with open(raname,'w') as fp:
                fp.write("#!/bin/bash\n"
                         "#BSUB -J {jobnamearr}\n"
                         "#BSUB -oo ./chunk%I/{jobname}.%I.oe\n"
                         "#BSUB -R \"linux64 && rhel60 && scratch > 6\"\n"
                         "#BSUB -We 24:00\n"
                         "#BSUB -W 48:00\n"
                         "\n"
                         "{extracmds}\n"
                         "./runarray.py $LSB_JOBINDEX\n"
                         "\n"
                         .format(extracmds=ec,
                                 jobname=jobname,
                                 jobnamearr=jobnamearr,
                                 odir='array_job_output'))

            os.system('chmod 755 %s' % raname)

            os.system('cd %s && bsub -q %s < jobarray.sh && cd -' % (files['work_output_dir'],self['queue']))
