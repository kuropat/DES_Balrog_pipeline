#!/usr/bin/env python
""" The program to run test Balrog simulation
 
 By N. Kuropatkin  12/06/2017
 """
import os
import sys
import math
import string
import shutil
import getopt
import subprocess
import yaml
import fitsio
import pickle
from despyastro import wcsutil

#
import time
import timeit
from time import sleep
import urllib2, ssl
import easyaccess
import cx_Oracle
from pox import shutils
import numpy
import glob
#
from multiprocessing import Pool

from ImageMath import ops
from pyparsing import line



try:
    from termcolor import colored
except:
    def colored(line, color):
        return line
from BalrogPipelineSof import BalrogPipelineSof

if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 4:
        "Usage: RunSim.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <sof conf file>"

        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:t:m:",["confile=","tilename","mofconfile"])
    except getopt.GetoptError:
        print "Usage:RunSim.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <sof conf file>"
        sys.exit(2)
    c_flag = 0
    t_flag = 0
    m_flag = 0
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage:RunSim.py<required inputs>"
            print "  Required inputs:"
            print "  -c <confile> - configuration file"
            print " -t <tile> - tile name"
            print " -m <sof conf file>"
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
        print "Usage:  RunSim.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <sof conf file>"
        sys.exit(-2)
    print " Start with Meds  \n"
    """ We start with the main line - create mads and sof for original tile """
    balP = BalrogPipelineSof(confile,tilename,mofconf)
#    balP.prepMeds()
    datadir = balP.tiledir
    basedir = str(balP.medsdir)+'/' + str(balP.medsconf)
    saveS = True
    print "save seeds ",saveS
    print "\n"
    

    ncpu = len(balP.bands)#
#
    tilelistF = 'TileList.csv'
    tilelistS = './TileList.csv'
    outf=open(tilelistF,'w')
    outf.write(tilename+'\n')
    outf.close()
    confFile = './Balrog-GalSim/config/bal_config.yaml'
    geom_file = './Balrog-GalSim/inputs/Y3A2_COADDTILE_GEOM.fits'
    config_dir = './Balrog-GalSim/config/'
    psf_dir = datadir
    tile_dir = basedir
    config_file = 'bal_config.yaml'
   
    output_dir =  basedir

    command = "./Balrog-GalSim/balrog/balrog_injection.py -l %s -g %s -t %s -c %s  -o %s -v %s " % (tilelistS,geom_file,tile_dir,config_dir,output_dir,config_file) 

    print command
    res=''
    try:
        retval = subprocess.call(command.split(),stderr=subprocess.STDOUT)
#        res = subprocess.check_output(command,shell=True)  
    except subprocess.CalledProcessError as e:
        print "error %s"% e
        print retval
    print retval