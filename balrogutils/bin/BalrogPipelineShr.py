#!/usr/bin/env python
""" The program to run mcal (shear )  processing on moffiles
  It suppose to be an extra step in mof production
  It does not work on sof.
 
 By N. Kuropatkin  16/04/2018
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
from Cython.Compiler.Naming import self_cname
from pox import shutils
import numpy
import glob
#
from multiprocessing import Pool
from Cython.Compiler.PyrexTypes import BuiltinObjectType
from IcoImagePlugin import IcoImageFile
from ImageMath import ops
from pyparsing import line
from scipy.weave.catalog import os_dependent_catalog_name
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
    
" Method to run mads maker ----------- "        
def makeMeds(inpar):
    args = inpar[0]
    band = inpar[1]
    ra_cent = args['ra_cent']
    dec_cent = args['dec_cent']
    naxis1 = args['naxis1']
    naxis2 = args['naxis2']
    pixscale = args['pixscale']
   
    outpath = args['outpath']
    logpath = args['logpath'] 
    medsdir = args['medsdir']
    medsconf = args['medsconf']
    tilename = args['tilename']
    confile = args['confile']
    datadir = args['datadir']
    restemp = outpath+'/'+tilename
    logfile = logpath+'/'+tilename    
    fname = datadir+'/lists/'+tilename+'_'+band+'_fileconf-'+medsconf+'.yaml'
    with open(fname,'r') as fconf:
        conf=yaml.load( fconf)
        print "file conf \n"
        print conf
        objmap_url = datadir+'/lists/'+tilename+'_'+band+'_objmap-'+medsconf+'.fits'
        seg_url = outpath+'/'+tilename+'_'+band+'_segmap.fits'
        cimage_url = outpath+'/'+tilename+'_'+band+'.fits'
        cat_url = outpath+'/'+tilename+'_'+band+'_cat.fits'
        conf['coadd_object_map'] = objmap_url
        conf['coadd_image_url'] = cimage_url
        conf['coadd_cat_url'] = cat_url
        conf['coadd_seg_url'] = seg_url
        with open(fname, 'w') as conf_f:
            yaml.dump(conf, conf_f) 
        command = ['desmeds-make-meds-desdm',confile,fname]
        print command
        try:
            subprocess.check_output(command)
        except:
            print "on meds for band %s \n" % band
            
            
" The method to create one shr chunk to be run by pool "

def makeChunk(inpar):
    " create shr chunk  "
    args=inpar[0]
    cnum = inpar[1]
    tilename = args['tilename']
    mofdir = args['mofdir']
    sheardir = args['sheardir']
    datadir = args['datadir']
    nchunks = args['nchunks']
    psfmapF = args['psfmap'] #medsdir+'/DES0347-5540_all_psfmap.dat'
    mofconf =args['mofconf']
    pardir = args['pardir']
    moffile = tilename+'_'+pardir+'.fits'
    medslist = datadir+'/'+pardir+'/meds_list.txt'
    seedlist = args['seedlist']

    listofmeds = []

    for line in open(medslist):
        medsF = line.split(' ')[0]
        listofmeds.append(medsF)
    print listofmeds
    nfit = ngmixit_tools.find_number_meds(listofmeds[0])

    # Get the job bracketing
    j1,j2 = ngmixit_tools.getrange(int(cnum),nfit,nchunks)
    print "# Found %s nfit" % nfit
    print "# will run chunk %s, between %s-%s jobs" % (nchunks,j1,j2)
    seedN = seedlist[cnum-1]
    mfile = datadir+'/'+pardir+'/'+moffile
    outfile =  sheardir+'/'+tilename+'_shr-chunk-'+'%02d'%cnum +'.fits'
    print seedN
    print '\n'

    command=['ngmixit', '--seed' ,'%d' %seedN, '--fof-range', '%d,%d'%(j1,j2)]
    command+=['--mof-file', mfile,'--psf-map',psfmapF,mofconf,outfile]
    for medsName in listofmeds:
        command+=[medsName]
    print command
    print '\n'
    try:
        subprocess.check_output(command)
    except:
        print "failed to run ngmixit \n"
       

class BalrogPipelineShr():
    
    def __init__(self, confile, tilename, mof_conf,pardir):
        '''
        Constructor
        '''
        urlbase = "https://desar2.cosmology.illinois.edu/DESFiles/"
        self.confile = confile
        self.workdir = os.getcwd()
        self.tilename = tilename
        self.pardir = pardir
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

        self.tiledir = self.medsdir+'/'+self.medsconf+'/'+self.tilename
        self.mofdir = self.tiledir +'/shr'
        if not os.path.exists(self.mofdir):
            os.makedirs(self.mofdir)
        self.realDir = ''
        self.curdir = os.getcwd()


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
        dict1['MOFFILE'] = self.moffile
        
        return dict1     

    """ establish connection to the database """
    def connectDB(self):        
        self.user = self.desconfig.get('db-' + self.dbname, 'user')
        self.dbhost = self.desconfig.get('db-' + self.dbname, 'server')
        self.port = self.desconfig.get('db-' + self.dbname, 'port')
        self.password = self.desconfig.get('db-' + self.dbname, 'passwd')
        kwargs = {'host': self.dbhost, 'port': self.port, 'service_name': self.dbname}
        self.dsn = cx_Oracle.makedsn(**kwargs)
        if not self.quiet: print('Connecting to DB ** %s ** ...' % self.dbname)
        connected = False
        
        for tries in range(3):
            try:
                self.con = cx_Oracle.connect(self.user, self.password, dsn=self.dsn)
                if self.autocommit: self.con.autocommit = True
                connected = True
                break
            except Exception as e:
                lasterr = str(e).strip()
                print(colored("Error when trying to connect to database: %s" % lasterr, "red"))
                print("\n   Retrying...\n")
                sleep( 8 )
        if not connected:
            print('\n ** Could not successfully connect to DB. Try again later. Aborting. ** \n')
            os._exit(0)
        self.cur = self.con.cursor()
        
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
    
    def get_tile_infor(self):
        coadd_tile_table = 'COADDTILE_GEOM'
        if self.dbname == 'desoper':
            coadd_tile_table='COADDTILE_GEOM'
        elif self.dbname == 'dessci':
            coadd_tile_table='Y3A2_COADDTILE_GEOM'
        query = QUERY_GEOM.format(tilename=self.tilename,coadd_tile_table=coadd_tile_table)
        print query
        self.cur.execute(query)
        desc = [d[0] for d in self.cur.description]
        # cols description
        line = self.cur.fetchone()
#        self.cur.close()
        # Make a dictionary/header for all the columns from COADDTILE table
        tileinfo = dict(zip(desc, line))
        return tileinfo

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

    def makeChunkList(self,datadir):
        mofdir = datadir + '/shr'
        chunklistF = mofdir +'/'+self.tilename+'_shr-chunk.list'
        if os.path.exists(chunklistF):
            os.remove(chunklistF)
        listF = glob.glob(mofdir+"/*shr-chunk*.fits")
        outF = open(chunklistF,'w')
        for fileN in listF:
                outF.write('%s \n' % fileN)
        outF.close()
        return chunklistF
     

    def collate_chunks(self,datadir):
        chunklistF = self.makeChunkList(datadir)
        mofdir = datadir + '/shr'
        command=['megamix-meds-collate-desdm','--noblind', 'shr-config/run-Y3A1-v001-mcal.yaml']
        command += [chunklistF,mofdir+'/'+self.tilename+'_shr.fits']
        try:
            subprocess.check_output(command)
        except:
            print "failed to collate chunks \n"
            

    
    def getMofPars(self,datadir):
        mofdir = datadir + '/shr'
        psfmapF = datadir +'/'+self.tilename+'_all_psfmap.dat'
        medslF = datadir+ '/'+ self.pardir +'/meds_list.txt'

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
        mofdir = datadir+'/shr'
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
            
            
    " Prepare data for meds run on simulated images "
    def prepInData(self):
        " make list of realizations "
        realizationslist = os.listdir(self.simData)
        " each realization contains tilename subrirectory "
        for realD in realizationslist:
            " new meds dir is simData+/realisatiuon+ tilename"
            self.realDir = os.path.abspath(self.simData+'/'+str(realD)+'/'+self.tilename)
            print " realDir=%s \n" % self.realDir
            " create subdirectories in the realization directory "
            for band in self.bands:
                nullw = self.realDir + '/nullwt-'+str(band)
#                print "nullweights = %s \n" % nullw
                if not os.path.exists(nullw):
                    os.makedirs(nullw)
            coaddir = self.realDir + '/coadd/'
            if not os.path.exists(coaddir):
                os.makedirs(coaddir)
                os.makedirs(coaddir+'/lists')
                os.makedirs(coaddir+'/LOGS')
            listdir = self.realDir + '/lists'
            if not os.path.exists(listdir):
                os.makedirs(listdir)
            " now copy common files "
            psffiles = glob.glob(self.tiledir+'/*.dat')
            for datF in psffiles:
                dstF = datF.split('/')[-1] 
                shutil.copyfile(datF, self.realDir+'/'+dstF)
            " Copy list files "
            listFiles = glob.glob(self.tiledir+'/lists/*')
#            newLists = []
            for listF in listFiles:
                destF = listF.split('/')[-1]
#                newLists.append(destF)
                if listF.find('nullwt') > 0:
                    self.changeList(listF,self.tiledir,self.realDir)
                else:
                    shutil.copyfile(listF, self.realDir+'/lists/'+destF)
 
            " now nullweight images are prepared in injection step "
            for band in self.bands:
                " We just need to copy and rename injected files in nullweight "
                srcDir = self.realDir + '/'+str(band) + '/'
                dstDir = self.realDir + '/nullwt-'+str(band) +'/'
                flist = glob.glob(srcDir+'*')
                for fileN in flist:
                    realF = fileN.split('/'+band+'/')[1]
                    dstFN = realF.split('balrog')[0]+'immasked_nullwt.fits'
                    shutil.copyfile(fileN, dstDir+dstFN)
                " Finally correct file config yaml "
                self.YamlCorrect(band,self.realDir)
            
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
    " The method to create one shr chunk  "
    def makeChunk(self,inpar):
        print " create shear chunk  \n"
        args=inpar[0]
        cnum = inpar[1]
        print " input Cnum=%d\n" %cnum
        tilename = args['tilename']
        mofdir = args['mofdir']
        sheardir = args['sheardir']
        datadir = args['datadir']
        nchunks = args['nchunks']
        psfmapF = args['psfmap'] #medsdir+'/DES0347-5540_all_psfmap.dat'
        mofconf =args['mofconf']
        pardir = args['pardir']
        moffile = tilename+'_'+pardir+'.fits'
        medslist = datadir+'/'+pardir+'/meds_list.txt'
        seedlist = args['seedlist']

        listofmeds = []

        for line in open(medslist):
            medsF = line.split(' ')[0]
            listofmeds.append(medsF)
        print listofmeds
        nfit = ngmixit_tools.find_number_meds(listofmeds[0])

    # Get the job bracketing
        j1,j2 = ngmixit_tools.getrange(int(cnum),nfit,nchunks)
        print "# Found %s nfit \n" % nfit
        print "# will run chunk %s, between %s-%s jobs \n" % (cnum,j1,j2)
        seedN = seedlist[cnum-1]
        mfile = datadir+'/'+pardir+'/'+moffile
        outfile =  sheardir+'/'+tilename+'_shr-chunk-'+'%02d'%int(cnum) +'.fits'
        print seedN
        print '\n'

        command=['ngmixit', '--seed' ,'%d' %seedN, '--fof-range', '%d,%d'%(j1,j2)]
        command+=['--mof-file', mfile,'--psf-map',psfmapF,mofconf,outfile]
        for medsName in listofmeds:
            command+=[medsName]
        print " created command \n"
        print command
        print '\n'
        res = ''
        try:
            res=subprocess.check_output(command)
        except:
            print "failed to run ngmixit \n"
            print res          
        print res
        
                
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
              

    
    " clean chunk file in the shr directory "    
    def cleanChunks(self,datadir):
        chunklist = glob.glob(datadir+'/shr/*chunk*.fits')
        for chunkF in chunklist:
            os.remove(chunkF)   


        
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 5:
        "Usage: BalrogPipelineShr.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <shr conf file>"
        print " -p <parent dir mof or sof >"

        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:t:m:p:",["confile=","tilename","mofconfile","parentdir"])
    except getopt.GetoptError:
        print "Usage:BalrogPipelineShr.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <shr conf file>"
        print " -p <parent dir mof or sof >"
        
        sys.exit(2)
    c_flag = 0
    t_flag = 0
    m_flag = 0
    p_flag = 0
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage:BalrogPipelineShr.py<required inputs>"
            print "  Required inputs:"
            print "  -c <confile> - configuration file"
            print " -t <tile> - tile name"
            print " -m <shr conf file>"
            print " -p <parent dir mof or sof >"
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
        elif opt in ("-p","--pardir"):
            p_flag = 1
            pardir = arg    
    sumF = c_flag + t_flag + m_flag + p_flag
    if sumF != 4:
        print "Usage:  BalrogPipelineShr.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <shr conf file>"
        print " -p <parent dir mof or sof >"
        sys.exit(-2)

    """ We start with the main line - create shear for original tile """
    balP = BalrogPipelineShr(confile,tilename,mofconf,pardir)
#    balP.prepMeds()
    datadir = balP.tiledir
     
    saveS = True
    print "save seeds ",saveS
    print "\n"

    nchunks = 16
    seedlist = makeSeedList(nchunks)
    if saveS:
        with open('seedL1', 'wb') as fp:
           pickle.dump(seedlist, fp)
    else:
        with open ('seedL1_base', 'rb') as fp:
           seedlist = pickle.load(fp)
    print " Seeds for base \n"
    print seedlist
    args =balP.getMofPars(datadir)
    args['mofdir'] = datadir+'/' + pardir
    args['sheardir'] = datadir+'/shr'
    args['datadir'] = datadir
    args['seedlist'] = seedlist
    args['nchunks'] = nchunks
    args['pardir'] = pardir
    if not os.path.exists(datadir + '/shr'):
        os.makedirs(datadir + '/shr')
#    balP.make_nbrs_data(datadir)

    pars = [(args, chunks) for chunks in range(1,nchunks+1) ]
#    print pars
#
    pool = Pool(processes=nchunks)

    pool.map(makeChunk, pars) 
      
    pool.close()
    pool.join()
#    for chunk in range(0,nchunks):
#        balP.makeChunk(pars[chunk])
    print " Finish producing chunks \n"
    balP.collate_chunks(datadir)
    balP.cleanChunks(datadir)


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
        args['mofdir'] = datadir+'/' + pardir
        args['sheardir'] = datadir+'/shr'
        args['datadir'] = datadir
        args['seedlist'] = seedlist
        args['nchunks'] = nchunks
        args['pardir'] = pardir

        if not os.path.exists(datadir + '/shr'):
            os.makedirs(datadir + '/shr')
        pars = [(args, chunks) for chunks in range(1,nchunks+1) ]
#        print pars
#
        pool = Pool(processes=nchunks)
        pool.map(makeChunk, pars) 
#        for chunk in range(0,nchunks):
#            balP.makeChunk(pars[chunk])  
        pool.close()
        pool.join()
        print "Finish with realization %d chunks \n" % int(real)
        balP.collate_chunks(datadir)
        balP.cleanChunks(datadir)
