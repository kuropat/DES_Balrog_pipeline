#!/usr/bin/env python
""" The program to run Sof on the base files
 meds
 
 By N. Kuropatkin  05/03/2018
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



import numpy
from numpy.random.mtrand import seed
import glob
#
from multiprocessing import Pool
from ngmixer import ngmixit_tools 
try:
    from termcolor import colored
except:
    def colored(line, color):
        return line



""" create a list of random numbers to be used as seeds """
def makeSeedList(nchunks):
    arrayV = numpy.random.rand(nchunks)
#    print arrayV
    seedlist=[]
    for c in range(0,nchunks):
       seedlist.append(int(arrayV[c]*10000))
    return seedlist

" The method to create one sof chunk to be run by pool "
def makeChunk(inpar):
    " create sof chunk  "
    args=inpar[0]
    cnum = inpar[1]
    tilename = args['tilename']
    mofdir = args['mofdir']
    datadir = args['datadir']
    medslist = args['medslist'] #mofdir+'/meds_list.txt'
    medsdir = args['medsdir'] #/data/des71.a/data/kuropat/meds_test/y3v02/DES0347-5540'
    psfmapF = args['psfmap'] #medsdir+'/DES0347-5540_all_psfmap.dat'
    mofconf =args['mofconf']
    nchunks = args['nchunks']
    seedlist = args['seedlist']
    
    listofmeds = []
    for line in open(medslist):
        medsF = line.split(' ')[0]
        listofmeds.append(medsF)
    nfit = ngmixit_tools.find_number_meds(listofmeds[0])

    # Get the job bracketing
    j1,j2 = ngmixit_tools.getrange(int(cnum),nfit,nchunks)
    print "# Found %s nfit" % nfit
    print "# will run chunk %s, between %s-%s jobs" % (cnum,j1,j2)
#    print seedlist
    seedN = seedlist[cnum-1]
    print seedN
    print '\n'
    outfile =  mofdir+'/'+tilename+'_sof-chunk-'+'%02d'%cnum +'.fits'
    command=['ngmixit', '--seed' ,'%d' %seedN, '--fof-range', '%d,%d'%(j1,j2)]
    command += ['--nbrs-file',mofdir+'/'+tilename+'_nbrslist.fits']
    command+=['--psf-map',psfmapF,mofconf,outfile]
    for line in listofmeds:
        medsN = line
        command+=[medsN]

    print command
    print '\n'
    try:
        subprocess.check_output(command)
    except:
        print "failed to run ngmixit \n"
        
            
#" The method to create one sof chunk to be run by pool "
#def makeChunk(inpar):
#    " create sof chunk  "
#    args=inpar[0]
#    cnum = inpar[1]
#    tilename = args['tilename']
#    mofdir = args['mofdir']
#    datadir = args['datadir']
#    medslist = args['medslist'] #mofdir+'/meds_list.txt'
#    medsdir = args['medsdir'] #/data/des71.a/data/kuropat/meds_test/y3v02/DES0347-5540'
#    psfmapF = args['psfmap'] #medsdir+'/DES0347-5540_all_psfmap.dat'
#    mofconf =args['mofconf']
#    seedlist = args['seedlist']
#    nchunks = args['nchunks']
#    print seedlist
#    seedN = seedlist[cnum-1]
#    print seedN
#    command = ['run_ngmixit','--nranges', '%d'%int(nchunks),  '--wrange', '1' ,'--tilename', tilename]
#    command += ['--meds_list',medslist,'--bands','g,r,i,z']
#    command += ['--psf-map', psfmapF]
#    command += ['--nbrs-file',mofdir+'/'+tilename+'_nbrslist.fits']
#    command += [mofconf,mofdir+'/'+tilename+'_sof-chunk-01.fits','--seed','%d' %seedN]
#    command[4] = str(cnum)
#    wrange=int(cnum)
#    nranges=nchunks 
#    command[16] = mofdir+'/'+tilename+'_sof-chunk-'+('%02d'%cnum)+'.fits' 
#
#    print command
#    try:
#        res=subprocess.check_output(command)
#    except:
#        print "failed to run run_ngmixit \n"
#    print res   

class BalrogSof():
    
    def __init__(self, confile, tilename, mof_conf):
        '''
        Constructor
        '''
        self.confile = confile
        self.workdir = os.getcwd()
        self.tilename = tilename
        self.conf = self.read_config(confile)
        self.dbname = self.conf['dbname']
        self.bands = self.conf['bands']
        self.det = self.conf['det']
        self.medsconf = self.conf['medsconf']
        print self.medsconf
        self.medsdir = os.getenv('MEDS_DATA')
        self.bbase = os.getenv('BALROG_BASE')
        self.desdata = os.getenv('DESDATA')
        self.simData = self.medsdir + '/'+self.medsconf+'/balrog_images/'
        print "simData = %s \n" % self.simData 
        print "meds_dir=%s \n" % self.medsdir
        self.mofconfile = mof_conf
        self.pixscale=0.2636
        self.magbase = 30.

        self.curdir = os.getcwd()
        self.autocommit = True
        self.quiet = False
        print " Make coadds for bands: \n"
        print self.bands
        print " detection image : \n"
        print self.det

#

        self.tiledir = self.medsdir+'/'+self.medsconf+'/'+self.tilename
        self.mofdir = self.tiledir +'/sof'
        if not os.path.exists(self.mofdir):
            os.makedirs(self.mofdir)
        self.realDir = ''
        self.curdir = os.getcwd()
#
        self.outpath = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename+'/coadd/')
        self.logpath = shutil.abspath(self.outpath)+'/LOGS/'

    def getRealizations(self):
        self.realizationslist  = os.listdir(self.simData)
        return self.realizationslist
            
    def setArgs(self,**args):
        self.inargs = args
        
    def ConfigMap(self):
        dict1 = {}
        dict1['MEDS_DATA'] = self.medsdir
        dict1['BALROG_BASE'] = self.bbase
        dict1['bands'] = self.bands
        dict1['det'] = self.det
        dict1['TILENAME'] = self.tilename
        dict1['TILEINFO'] = self.tileinfo
        dict1['TILEDIR'] = self.tiledir
        dict1['MOFDIR'] = self.mofdir
        
        return dict1     

        
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
    

                    
    " read catalog to find its length and create a fake objmap"
    def make_fake_objmap(self,catfile,fname):
        fitsio.read(catfile)
        coadd_cat = fitsio.read(catfile, lower=True)
        print "Creating objmap %s \n" % fname
        # sort just in case, not needed ever AFIK
        q = numpy.argsort(coadd_cat['number'])
        coadd_cat = coadd_cat[q]
        nobj = len(coadd_cat) 
        print " Number of objects in the cat = %d \n" % nobj 
        # now write some data
        if os.path.exists(fname):
            os.remove(fname)
        fits = fitsio.FITS(fname,'rw')

        data = numpy.zeros(nobj, dtype=[('object_number','i4'),('id','i8')])
       
        data['object_number'] = [num for num in range(nobj)]
        data['id'] = [long(num) for num in range(nobj)]

        print("writing objmap:",fname)
        fitsio.write(fname, data, extname='OBJECTS',clobber=True)  

    def make_meds_list(self,datadir):
        mofdir = datadir + '/sof'
        medsLF = mofdir+'/meds_list.txt'
        outF = open(medsLF,'w')
        listF = glob.glob(datadir+"/*meds*.fits.fz")
        medsdic = {}
        for line in listF:
            filename = line.split('/')[-1]
            bandN = filename.split('_')[1]
            medsdic[bandN] = datadir + '/' + filename            
        for band in self.bands:
            bandN = str(band)
            outF.write("%s %s \n" % (medsdic[bandN],bandN)) 
        outF.close()    
        return medsLF
    
    def makeChunkList(self,datadir):
        mofdir = datadir + '/sof'
        chunklistF = mofdir +'/'+self.tilename+'_sof-chunk.list'
        if os.path.exists(chunklistF):
            os.remove(chunklistF)
        listF = glob.glob(mofdir+"/*sof-chunk*.fits")
        outF = open(chunklistF,'w')
        for fileN in listF:
                outF.write('%s \n' % fileN)
        outF.close()
        return chunklistF
     

    def collate_chunks(self,datadir):
        chunklistF = self.makeChunkList(datadir)
        mofdir = datadir + '/sof'
        command=['megamix-meds-collate-desdm','--noblind', 'sof-config/run-Y3A1-v4-sof.yaml']
        command += [chunklistF,mofdir+'/'+self.tilename+'_sof.fits']
        try:
            subprocess.check_output(command)
        except:
            print "failed to collate chunks \n"
            

    
    def getMofPars(self,datadir):
        mofdir = datadir + '/sof'
        psfmapF = datadir +'/'+self.tilename+'_all_psfmap.dat'
        medslF = mofdir+'/meds_list.txt'

        pars ={}
        pars['medsdir'] = self.medsdir
        pars['mofdir'] = mofdir
        pars['datadir'] = datadir
        pars['tilename'] = self.tilename
        pars['psfmap'] = psfmapF
        pars['medslist'] = medslF
        pars['mofconf'] = self.mofconfile 
        return pars
        
    def make_nbrs_data(self,datadir):
        " First create mads list "
        mofdir = datadir+'/sof'
        if not os.path.exists(mofdir):
            os.makedirs(mofdir)
        medslist = self.make_meds_list(datadir)
        
        " Now create psfmap for all bands "
        self.make_psf_map(datadir)
        "  Second run run_ngmixer-meds-make-nbrs-data for all bands "
        command = ['run_ngmixer-meds-make-nbrs-data','./'+self.mofconfile,'--nbrs-file',mofdir+'/'+self.tilename+'_nbrslist.fits']
        command +=['--fof-file',mofdir+'/'+self.tilename+'_fofslist.fits','--meds_list',medslist]
        print command
        try:
            subprocess.check_output(command)
        except:
            print "failed to run run_ngmixer-meds-make-nbrs-data \n"
            
    def make_psf_map(self,datadir):
        psfmapF = datadir +'/'+self.tilename+'_all_psfmap.dat'
        if os.path.exists(psfmapF):
            os.remove(psfmapF)
        listF = glob.glob(datadir+"/*psfmap*.dat")
        outF = open(psfmapF,'w')
        for fileN in listF:
            for line in open(fileN,'r'):
                outF.write(line)
        outF.close()
        return psfmapF
        
    
                
    " This is the sequence of command composing the pipeline "    
    def prepMeds(self):
        "  First run desmeds-prep-tile "
        
        for band in self.bands:
            command = [self.bbase+'/bin/desmeds-prep-tile',self.medsconf,self.tilename,band]
            print command
            try:
                subprocess.check_output(command)
            except:
                print "failed to copy files for band %s \n" % band
   
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
        if not os.path.exists(self.logpath):
            os.makedirs(self.logpath)
            
    def makeIngection(self):
        " At this point we need to include the image injection commands "
        "                                                               "
        "                                                               "     
        " prepare information for the coadd tile "
        
    def getCoaddPars(self,datadir):
        outpath = datadir+'/coadd/'
        logpath = outpath +'/LOGS/'
        outpathL = outpath+'/lists/'
        
        pard = {}
        pard['ra_cent'] = self.tileinfo['RA_CENT']
        pard['dec_cent'] = self.tileinfo['DEC_CENT']
        pard['naxis1'] = self.tileinfo['NAXIS1']
        pard['naxis2'] = self.tileinfo['NAXIS2']
        self.pixscale = self.tileinfo['PIXELSCALE']
        pard['pixscale'] = self.pixscale
        
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        if not os.path.exists(logpath):
            os.makedirs(logpath)

        if not os.path.exists(outpathL):
            os.makedirs(outpathL)
        pard['datadir'] = datadir
        pard['outpath'] = outpath
        pard['logpath'] = logpath
        pard['listpath'] = outpathL
        pard['tilename'] = self.tilename
        pard['medsdir'] = self.medsdir
        pard['medsconf'] = self.medsconf
        pard['magbase'] = self.magbase
        pard['confile'] = self.confile
        return pard
    
    def run(self):
        self.prepMeds()
        byband = {}
        "  Run injection here "
        " Prepare data for injected images "
        self.prepInData()
        " Create coadd images for each band and each revision "
        for realV in self.realizationslist:
            outpathR = os.path.join(self.medsdir,self.medsconf+'/balrog_images/'+str(realV)+'/'+self.tilename+'/coadd')
            realD =  os.path.join(self.medsdir,self.medsconf+'/balrog_images/'+str(realV)+'/'+self.tilename)
            print " outpathR = %s \n" % outpathR
            logpathR = shutil.abspath(outpathR)+'/LOGS/'
            listpathR = shutil.abspath(outpathR)+'/lists/'
            for band in self.bands:
               self.makeCoadd(band,outpathR,self.tiledir)
            " Now make detection image "
            self.makeDetectionImage(outpathR)
            " Create s-extractor catalogs "
            for band in self.bands:
                self.makeCatalog(band,outpathR)
            " Now create fake objectMaps "
            self.makeObjMaps(realD)
            " Now crete meds for all bands "
            for band in self.bands:
                self.makeMeds(band,realD)
                
    def CorrectPath(self,psfList,basedir):
        outf = open("tempList",'w')
        for line in open(psfList,'r'):
            tokens = line.split(' ')
            outline = tokens[0]
            subtok = tokens[1].split(self.medsconf)
            outline +=' '+self.medsdir + '/'+self.medsconf+'/'+subtok[1]
            outf.write(outline)
        outf.close()
        shutil.move("tempList",psfList)           
    " Prepare data for meds run on simulated images "
    def prepInData(self):
        " make list of psfmaps in base "
        basedir = self.medsdir + '/'+self.medsconf+'/' + self.tilename
        allPSFmapFiles = glob.glob(basedir+"/*_psfmap*.dat")
        for psfLfile in allPSFmapFiles:
            self.CorrectPath(psfLfile,self.medsdir + '/'+self.medsconf+'/')
        " make list of realizations "
        realizationslist = os.listdir(self.simData)
        " each realization contains tilename subrirectory "
        for realD in realizationslist:
            " new meds  dir is simData+/realization+ tilename"
            self.realDir = os.path.abspath(self.simData+'/'+str(realD)+'/'+self.tilename)
            " correct psf reference in psfmaps "
            allPSFmapFiles = glob.glob(self.realDir+"/*_psfmap*.dat")
            for psfLfile in allPSFmapFiles:
                self.CorrectPath(psfLfile,self.realDir+'/')            

    def YamlCorrect(self,band,realDir):
        listDir = realDir+'/lists/'
        fileName = self.tilename+'_'+str(band)+'_fileconf-'+self.medsconf+'.yaml'
        fileList = listDir+fileName
#        print "fileList=%s \n" % fileList

        with open(fileList,'r') as fconf:
            conf=yaml.load( fconf)
            fconf.close()
        comp = 'nwgint_flist'
        line = conf['nwgint_flist']
        nwgFlist = line.split('/')[-1]
        line= realDir + '/lists/'+ nwgFlist
        conf['nwgint_flist'] = line
        line = conf['meds_url']
        medsF = line.split('/')[-1]
        line = realDir +'/'+medsF
        conf['meds_url'] = line
        with open(fileList, 'w') as conf_f:
            yaml.dump(conf, conf_f) 
# 
                

        
                
    def changeList(self,listF,srcDir,realDir):
        "  change path in each list "
        outFile = realDir+'/lists/' + listF.split('/')[-1]
        outF = open(outFile,'w')
        for line in open(listF,'r'):
            line.strip()
            newline = realDir+'/nullwt-'+line.split('nullwt-')[1]
            outF.write(newline)
        outF.close()

                
    def makeObjMaps(self,outpath):  
        print "makeObjMaps outpath=%s \n" % outpath     
        for band in self.bands:
#            fname = os.path.join(outpath,'/lists/'+self.tilename+'_'+band+'_objmap-'+self.medsconf+'.fits')
            fname = shutil.abspath(outpath) + '/lists/'+self.tilename+'_'+band+'_objmap-'+self.medsconf+'.fits'
            print "Obj map file %s \n" % fname
            catname = outpath+'/coadd/'+self.tilename+'_'+band+'_cat.fits'
            print "makeObjMaps catname=%s n" % catname
            if os.path.exists(fname):
                os.remove(fname)
            self.make_fake_objmap(catname, fname)
 
    def makeMeds(self,band,outpath):
        fname = outpath +'/lists/'+self.tilename+'_'+band+'_fileconf-'+self.medsconf+'.yaml'
        with open(fname,'r') as fconf:
            conf=yaml.load( fconf)
        print "file conf \n"
        print conf
        objmap_url = outpath +'/lists/'+self.tilename+'_'+band+'_objmap-'+self.medsconf+'.fits'
        seg_url = outpath+'/coadd/'+self.tilename+'_'+band+'_segmap.fits'
        cimage_url = outpath+'/coadd/'+self.tilename+'_'+band+'.fits'
        cat_url = outpath+'/coadd/'+self.tilename+'_'+band+'_cat.fits'
        conf['coadd_object_map'] = objmap_url
        conf['coadd_image_url'] = cimage_url
        conf['coadd_cat_url'] = cat_url
        conf['coadd_seg_url'] = seg_url
        with open(fname, 'w') as conf_f:
            yaml.dump(conf, conf_f) 
        command = [self.bbase+'/bin/desmeds-make-meds-desdm',self.confile,fname]
        print command
        try:
            subprocess.check_output(command)
        except:
            print "on meds for band %s \n" % band
              
    def makeDetectionImage(self,outpathR):   
        print " Images are created start with detection image \n"
        "Now we have coadd images let's make detection one "
        restemp = outpathR+'/'+self.tilename
        ra_cent = self.tileinfo['RA_CENT']
        dec_cent = self.tileinfo['DEC_CENT']
        naxis1 = self.tileinfo['NAXIS1']
        naxis2 = self.tileinfo['NAXIS2']
        self.pixscale = self.tileinfo['PIXELSCALE']

        (im_list,wt_list,msk_list) = self.create_det_list(restemp)

        restemp =  outpathR+'/'+self.tilename
        self.MakeDetImage(im_list,wt_list,msk_list,ra_cent,dec_cent,naxis1,naxis2,restemp) 
        restemp = outpathR +'/'+self.tilename+'_det'

        self.coadd_assemble(restemp)
    
    " clean chunk file in the sof directory "    
    def cleanChunks(self,datadir):
        chunklist = glob.glob(datadir+'/sof/*chunk*.fits')
        for chunkF in chunklist:
            os.remove(chunkF)   

 
    def makeCatalog(self,band,outpath):
        restemp = outpath+'/'+self.tilename
        logfile = outpath+'/LOGS/'+self.tilename    
        self.SEXcaller(restemp,band,logfile)
        
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 4:
        "Usage: BalrogSof.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <sof conf file>"

        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:t:m:",["confile=","tilename","mofconfile"])
    except getopt.GetoptError:
        print "Usage:BalrogSof.py <required inputs>"
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
            print "Usage:BalrogSof.py<required inputs>"
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
        print "Usage:  BalrogSof.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <sof conf file>"
        sys.exit(-2)

    """ We start with the main line - create mads and sof for original tile """
    balP = BalrogSof(confile,tilename,mofconf)
    datadir = balP.tiledir
    balP.prepInData() 
    saveS = True
    print "save seeds ",saveS
    print "\n"
    

    nchunks = 24
    seedlist = makeSeedList(nchunks)
#    if saveS:
#        with open('seedL1', 'wb') as fp:
#           pickle.dump(seedlist, fp)
#    else:
#        with open ('seedL1_base', 'rb') as fp:
#           seedlist = pickle.load(fp)
#    print " Seeds for base \n"
#    print seedlist
#    args =balP.getMofPars(datadir)
#    args['mofdir'] = datadir+'/sof'
#    args['datadir'] = datadir
#    args['seedlist'] = seedlist
#    args['nchunks'] = nchunks
    


    """ ---------------------------------------------------------"""

    reallist = balP.getRealizations()
    basedir = str(balP.medsdir)+'/' + str(balP.medsconf)
    for real in reallist:
        datadir = basedir+'/balrog_images/' + str(real)+'/'+balP.tilename
        

        
        seedlist = makeSeedList(nchunks)
        if saveS:
            with open('seedL2', 'wb') as fp:
                pickle.dump(seedlist, fp)
        else:
            with open ('seedL2_base', 'rb') as fp:
                seedlist = pickle.load(fp)
        print " Seeds for realization 0 \n"
        print seedlist
        args =balP.getMofPars(datadir)
        args['mofdir'] = datadir+'/sof'
        args['datadir'] = datadir
        args['seedlist'] = seedlist
        args['nchunks'] = nchunks
#
    
        balP.make_nbrs_data(datadir)

        pars = [(args, chunks) for chunks in range(1,nchunks+1) ]
#
#
        pool = Pool(processes=nchunks)
        pool.map(makeChunk, pars) 
           
        pool.close()
        pool.join()

        balP.collate_chunks(datadir)
        balP.cleanChunks(datadir)
