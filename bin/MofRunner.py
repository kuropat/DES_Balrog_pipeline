#!/usr/bin/env python
""" The program to organize all steps necessary to produce
 mof files for Balrog project
 
 By N. Kuropatkin  12/21/2017
 """
import os
import sys
#import math
#import string
#import shutil
import getopt
import subprocess
import yaml
#import fitsio
#from despyastro import wcsutil

#
#import time
#import timeit
#from time import sleep
#import urllib2, ssl
#import config_ea as config_mod
#import cx_Oracle
#from Cython.Compiler.Naming import self_cname
#from pox import shutils
#import numpy
import glob
from multiprocessing import Pool
#import pathos
#from pathos.multiprocessing import ProcessingPool as Pool
#from Cython.Compiler.PyrexTypes import BuiltinObjectType
from fileinput import filename
from astropy.modeling.tests.test_projections import pars
#try:
#    from termcolor import colored
#except:
#    def colored(line, color):
#        return line
    
" The method to create one mof chunk to be run by pool "
def makeChunk(inpar):
    " create mof chunk  "
    args=inpar[0]
    cnum = inpar[1]
    tilename = args['tilename']
    mofdir = args['mofdir']
    medslist = args['medslist'] #mofdir+'/meds_list.txt'
    medsdir = args['medsdir'] #/data/des51.a/data/kuropat/meds_test/y3v02/DES0347-5540'
    psfmapF = args['psfmap'] #medsdir+'/DES0347-5540_all_psfmap.dat'
    mofconf =args['mofconf']
    command = ['run_ngmixit','--nranges', '42',  '--wrange', '1' ,'--tilename', tilename]
    command += ['--meds_list',medslist,'--bands','g,r,i,z','--fof-file', mofdir+'/'+tilename+'_fofslist.fits']
    command += ['--psf-map', psfmapF]
    command += ['--nbrs-file',mofdir+'/'+tilename+'_nbrslist.fits']
    command += [mofconf,mofdir+'/'+tilename+'_mof-chunk-01.fits']
    command[4] = str(cnum)
    command[18] = mofdir+'/'+tilename+'_mof-chunk-'+('%02d'%cnum)+'.fits'       
    try:
        subprocess.check_output(command)
    except:
        print "failed to run run_ngmixit \n"
         


class MofRunner():
    
    def __init__(self, confile, tilename, mof_conf):
        '''
        Constructor
        '''
        self.confile = confile
        self.workdir = os.getcwd()
        self.tilename = tilename
        self.mofconfile = mof_conf
        self.conf = self.read_config(confile)
        self.bands = self.conf['bands']
        self.bendstr = self.makeBandStr()
        self.det = self.conf['det']
        self.medsconf = self.conf['medsconf']
        print self.medsconf
        self.medsdir = os.getenv('MEDS_DATA')
        self.tiledir = self.medsdir+'/'+self.medsconf+'/'+self.tilename
        self.bbase = os.getenv('BALROG_BASE')
        self.desdata = os.getenv('DESDATA')
        
        print "meds_dir=%s \n" % self.medsdir
        self.mofdir = self.tiledir +'/mof'
        if not os.path.exists(self.mofdir):
            os.makedirs(self.mofdir)
        print "mof_dir=%s \n" % self.mofdir
        self.curdir = os.getcwd()

    def makeBandStr(self):
        bandstr = ''
        for band in self.bands:
            bandstr+=band+','
        bandstr = bandstr[:-1]
        return bandstr

    def read_config(self,confile):
        """
        read the  config file

        parameters
        ----------
        dbname: string
            Name of DESDM databese to query tile info as 'desoper'
        """


        print("reading:",confile)
        with open(confile) as parf:
            data=yaml.load(parf)


        return data
    
    def make_meds_list(self):
        self.medslF = self.mofdir+'/meds_list.txt'
        outF = open(self.medslF,'w')
        listF = glob.glob(self.tiledir+"/*meds*.fits.fz")
        for line in listF:
            filename = line.split('/')[-1]
            bandN = filename.split('_')[1]
            outF.write('%s %s \n' % (self.tiledir+'/'+filename,bandN))
        return self.medslF
     
    def make_psf_map(self):
        self.psfmapF = self.tiledir +'/'+self.tilename+'_all_psfmap.dat'
        if os.path.exists(self.psfmapF):
            os.remove(self.psfmapF)
        listF = glob.glob(self.tiledir+"/*psfmap*.dat")
        outF = open(self.psfmapF,'w')
        for fileN in listF:
            for line in open(fileN,'r'):
                outF.write(line)
        outF.close()

    def makeChunkList(self):
        chunklistF = self.mofdir +'/'+self.tilename+'_mof-chunk.list'
        if os.path.exists(chunklistF):
            os.remove(chunklistF)
        listF = glob.glob(self.mofdir+"/*mof-chunk*.fits")
        outF = open(chunklistF,'w')
        for fileN in listF:
                outF.write('%s \n' % fileN)
        outF.close()
        return chunklistF
    
#    def makeChunk(self,cnum):
#        " create mof chunk  "
#        command = ['run_ngmixit','--nranges', '42',  '--wrange', '1' ,'--tilename', self.tilename]
#        command += ['--meds_list',self.medslist,'--bands','g,r,i,z','--fof-file', self.mofdir+'/'+self.tilename+'_fofslist.fits']
#        command += ['--psf-map', self.psfmapF]
#        command += ['--nbrs-file',self.mofdir+'/'+self.tilename+'_nbrslist.fits']
#        command += ['mof-config/run-Y3A1-v4-mof.yaml',self.mofdir+'/'+self.tilename+'_mof-chunk-01.fits']
#        command[4] = str(cnum)
#        command[18] = self.mofdir+'/'+self.tilename+'_mof-chunk-'+('%02d'%cnum)+'.fits'       

    def collate_chunks(self):
        self.chunklistF = self.makeChunkList()
        command=['megamix-meds-collate-desdm','--noblind', 'mof-config/run-Y3A1-v4-mof.yaml']
        command += [self.chunklistF,self.mofdir+'/'+self.tilename+'_mof.fits']
        try:
            subprocess.check_output(command)
        except:
            print "failed to collate chunks \n"
    
    def getMofPars(self):
        self.psfmapF = self.tiledir +'/'+self.tilename+'_all_psfmap.dat'
        self.medslF = self.mofdir+'/meds_list.txt'
        pars ={}
        pars['medsdir'] = self.medsdir
        pars['mofdir'] = self.mofdir
        pars['tilename'] = self.tilename
        pars['psfmap'] = self.psfmapF
        pars['medslist'] = self.medslF
        pars['mofconf'] = self.mofconfile 
        return pars
        
    def make_nbrs_data(self):
        " Firsc create mads list "
        self.medslist = self.make_meds_list()
        " Now create psfmap for all bands "
        self.make_psf_map()
        "  Second run run_ngmixer-meds-make-nbrs-data for all bands "
        command = ['run_ngmixer-meds-make-nbrs-data','./'+self.mofconfile,'--nbrs-file',self.mofdir+'/'+self.tilename+'_nbrslist.fits']
        command +=['--fof-file',self.mofdir+'/'+self.tilename+'_fofslist.fits','--meds_list',self.medslist]
        print command
        try:
            subprocess.check_output(command)
        except:
            print "failed to run run_ngmixer-meds-make-nbrs-data \n"


                
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 4:
        "Usage: MofRunner.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <mof conf file>"

        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:t:m:",["confile=","tilename","mofconfile"])
    except getopt.GetoptError:
        print "Usage:MofRunner.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <mof conf file>"
        sys.exit(2)
    c_flag = 0
    t_flag = 0
    m_flag = 0
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage: MofRunner.py<required inputs>"
            print "  Required inputs:"
            print "  -c <confile> - configuration file"
            print " -t <tile> - tile name"
            print " -m <mof conf file>"
            sys.exit(2)
           
        elif opt in ("-c","--confile"):
            c_flag = 1
            confile = arg 
        elif opt in ("-t","--tilename"):
            t_flag = 1
            tilename = arg
        elif opt in ("-m","--mofconfile"):
            m_flag = 1
            mofconf = arg
    sumF = c_flag + t_flag + m_flag
    if sumF != 3:
        print "Usage:  MakeSWarpCoadd.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <mof conf file>"
        sys.exit(-2)
    print " Start with MakeSWarpCoadd \n"
    mofR = MofRunner(confile,tilename,mofconf)
    mofR.make_nbrs_data()
    args= mofR.getMofPars()
    pars = [(args, chunks) for chunks in range(1,43) ]
    print pars
#
    pool = Pool(processes=16)
    pool.map(makeChunk, pars) 
           
    pool.close()
    pool.join()
    mofR.collate_chunks()
