#!/usr/bin/env python
""" The program to run sof with megamixer on meds files produced by BalrogBase program 
 By: N. Kuropatkin  08/08/2018
 
 """
import os
import sys
import string
import shutil
import math
import getopt
import glob
import subprocess
import yaml
import balrogutils.ngmixit_tools as ngmt 
from multiprocessing import Pool


" The routine to run chunk scripts in multiprocessing pool "

def runChunkScript(inpar):
    chunkbase = inpar[0]
    chunkdir = inpar[1]
    chunkext = chunkdir.split('chunk')[1]
    workdir = chunkbase + '/' + chunkdir
    tilename = inpar[2]
    medsconf = inpar[3]
#    os.chdir(workdir)
    scriptname = workdir + '/' + tilename+'-'+medsconf+'-sof-'+chunkext+'.sh'
    res = ''
    try:
        res=subprocess.check_output(scriptname,shell=True)  
    except subprocess.CalledProcessError as e:
        print "error %s"% e
        print res
    print res

class BalrogSofMegamixer():
    
    def __init__(self, confile, tilename):
        '''
        Constructor
        '''
        self.confile = confile
        self.workdir = os.getcwd()
        self.tilename = tilename
        self.tilesF = 'tiles.yaml'

        self.medsdir = os.getenv('MEDS_DIR')
        self.medsdata = os.getenv('MEDS_DATA')
        self.bbase = os.getenv('BALROG_BASE')
        self.medsconf = os.getenv('medsconf')
        self.outputbase = os.getenv('NGMIXER_OUTPUT_DIR')
        self.myenv = dict(os.environ).copy()
#        self.desdata = os.getenv('DESDATA')
        self.simData = self.medsdir + '/'+self.medsconf+'/balrog_images/'
#        self.BU = BalrogUtils(confile,tilename)
        print "medsconf = %s \n" % self.medsconf
        print "meds_dir=%s \n" % self.medsdir
        
        
    " The method to run megamixer that will create chunk files " 
    def runMegamixer(self,config,cmd,tiles):
        """   Possible commands:
        setup - setup jobs, this command works
        setup-nbrs - setup jobs for neighbors/fof finding
        collate - combine job outputs into a single file, this works
        verify - verify that all job outputs are present and OK, does not work
        clean - clean all outputs from a run, does not work
        archive - run after collate to delete intermediate files and tar logs, does not work
        """
        
        command = "megamix-meds --skip-errors --noblind %s %s %s" % (config,cmd,tiles)
        print command
#        envir = dict(os.environ).copy()
        try:
            #subprocess.check_output(command,shell=True,env=envir)
            subprocess.check_output(command,shell=True)
        except subprocess.CalledProcessError as e:
            print "error %s"% e 
            return 1
        return 0
    " Creates tiles info file for Megamixer "
    def makeTiles(self):
        " creates tiles.yaml file to be used by megamixer "
        
        tiledic = {"medsconf": self.medsconf ,"tile_ids":[str(self.tilename),]}
              
        with open(self.tilesF, 'w') as conf_f:
            yaml.dump(tiledic, conf_f)   
    
    def getChunkSize(self,ncpu):
        workdir = self.medsdir + '/' + self.medsconf+'/'+self.tilename
        allfiles = os.listdir(str(workdir))
        listofmeds = []
        for line in allfiles:
            if line.find('meds') > 0:
                medsF = line
                listofmeds.append(workdir+'/'+medsF)
        nfit = ngmt.find_number_meds(listofmeds[0])
        chunksize = int(nfit/(ncpu - 0.5))

        return chunksize
      
    def makeNewConf(self,outconf,chunksize):
        outconf = open(outconf,'w')
        for line in open(self.confile):
            if line.find('chunksize') > 0:
                tokens = line.split(':')
                line = tokens[0]+': '+str(chunksize)
                outconf.write(line)
            else:
                outconf.write(line)
        outconf.close()
            
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 4:
        "Usage: BalrogSofMegamixer.py  <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -n <number of CPUs to use> "
        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:t:n:",["confile=","tilename","ncpu"])
    except getopt.GetoptError:
        print "Usage: BalrogSofMegamixer.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -n <number of CPUs to use> "
        sys.exit(2)
    c_flag = 0
    t_flag = 0
    n_flag = 0
    ncpu = 1
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage:  BalrogSofMegamixer.py <required inputs>"
            print "  Required inputs:"
            print "  -c <confile> - configuration file"
            print " -t <tile> - tile name"
            print " -n <number of CPUs to use> "
            sys.exit(2)
           
        elif opt in ("-c","--confile"):
            c_flag = 1
            confile = arg 
        elif opt in ("-t","--tilename"):
            t_flag = 1
            tilename = arg
        elif opt in ("-n","--ncpu"):
            n_flag = 1
            ncpu = int(arg)
    sumF = c_flag + t_flag + n_flag
    if sumF != 3:
        print "Usage: BalrogSofMegamixer.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -n <number of CPUs to use> "
        sys.exit(-2)
    print " Start with megamixer  \n"
    """  """
    balP = BalrogSofMegamixer(confile,tilename)
    balP.makeTiles()
    tilesF = balP.tilesF
#    ncpu = 16
    chunksize = balP.getChunkSize(ncpu)
    confname = confile.split('/')[-1]
    print "cofig file name = %s \n" % confname
    newconf = str(os.getenv('TMPDIR'))+'/' +confname
    balP.makeNewConf(newconf,chunksize)
    print "New cofig file=%s \n" % newconf
    " create chunks "
    status = balP.runMegamixer(newconf,'setup',tilesF)
    print "Megamixer status= % d \n" % status
    " at this point we have chunk scripts ready "
    chankdirs = balP.medsdir + '/'+balP.medsconf+'-sof/' + tilename
    chunklist = os.listdir(chankdirs)
    pars = [(chankdirs,chunk,tilename,balP.medsconf) for chunk in chunklist]
    print pars
#    ncpu = 16
    pool = Pool(processes=ncpu)
    it = pool.imap(runChunkScript, pars)
    for i in  it:
        continue
    pool.close()
#    pool.join()
    " Now collate chunks "
    status = balP.runMegamixer(newconf,'collate',tilesF)
    print "collate status=%d \n" % status
    
    sys.exit(status)
