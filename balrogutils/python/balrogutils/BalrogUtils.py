#!/usr/bin/env python
""" The helper class to work with BAlrog Pipelines.
 By: N. Kuropatkin  08/06/2018
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
import numpy
import re

from time import sleep
#import urllib2, ssl

import easyaccess
import cx_Oracle
import glob
#
#from multiprocessing import Pool
#from numpy.random.mtrand import seed
#from ngmixer import ngmixit_tools 


try:
    from termcolor import colored
except:
    def colored(line, color):
        return line


# The query template used to get the geometry of the tile
QUERY_GEOM = """
    SELECT ID, PIXELSCALE, NAXIS1, NAXIS2,
        RA_CENT, DEC_CENT,
        RA_SIZE,DEC_SIZE,
        RAC1, RAC2, RAC3, RAC4,
        DECC1, DECC2, DECC3, DECC4,
        RACMIN,RACMAX,DECCMIN,DECCMAX,
        CROSSRA0
    FROM {coadd_tile_table}
    WHERE tilename='{tilename}'
"""



class BalrogUtils():
    
    def __init__(self, confile, tilename):
        '''
        Constructor
        '''
        urlbase = "https://desar2.cosmology.illinois.edu/DESFiles/"
        self.confile = confile
        self.workdir = os.getcwd()
        self.tilename = tilename
        self.conf = self.read_config(confile)
        self.dbname = self.conf['dbname']
        self.bands = self.conf['bands']
        self.det = self.conf['det']
        self.medsconf = self.conf['medsconf'] # 'y3v02'
        print self.medsconf
        self.medsdir = os.getenv('MEDS_DATA')
        self.bbase = os.getenv('BALROG_BASE')
        self.desdata = os.getenv('DESDATA')
        self.simData = self.medsdir + '/'+self.medsconf+'/balrog_images/'
        print "simData = %s \n" % self.simData 
        print "meds_dir=%s \n" % self.medsdir
        self.mofconfile = ''
        self.pixscale=0.2636
        self.magbase = 30.

        self.curdir = os.getcwd()
        self.autocommit = True
        self.quiet = False
        print " Make coadds for bands: \n"
        print self.bands
        print " detection image : \n"
        print self.det
        desfile = os.getenv("DES_SERVICES")
        if not desfile: desfile = os.path.join(os.getenv("HOME"), ".desservices.ini")
#
        self.desconfig = easyaccess.config_ea.get_desconfig(desfile, self.dbname)
        self.connectDB()
        self.tileinfo = self.get_tile_infor()
        self.tiledir = self.medsdir+'/'+self.medsconf+'/'+self.tilename
        self.mofdir = self.tiledir +'/mof'
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


    def create_det_list(self,restemp):
        im_list = ''
        wt_list = ''
        msk_list = ''
        first = True
        for band in self.det:
            im_file = restemp+'_'+band+'_sci.fits'
            msk_file = restemp+'_'+band+'_msk.fits'
            wt_file = restemp+'_'+band+'_wgt.fits'
            if first:
                im_list +=str(im_file)
                msk_list += str(msk_file)
                wt_list += str(wt_file)
                first = False
            else:
                im_list += ','+im_file
                msk_list += ','+str(msk_file)
                wt_list += ','+wt_file 
        return (im_list,wt_list,msk_list)

    def CombineFits(self,resFile,simFile,origFile):
        template = resFile.split('/')[-1]
        restemp = template.split('.fits')[0]
        imName = restemp +'_sci.fits'
        weightName = restemp +'_wgt.fits'
        maskName = restemp + '_msk.fits'
        outName = restemp + '.fits'
        " read sci image "
        sciIm,imHdr = fitsio.read(simFile, ext =0, header=True)
        mskIm,mskHdr = fitsio.read(origFile, ext=1, header=True)
        wgtIm,wgtHdr = fitsio.read(origFile, ext=2, header=True)
        i1 = imHdr['NAXIS1']
        i2 = imHdr['NAXIS2']
        " Now write all 3 extensions in the output file "
#        print "resFile=%s \n" % resFile
        if os.path.exists(resFile):
            os.remove(resFile)
        fits = fitsio.FITS(resFile,'rw')
        
        fits.write( sciIm,header=imHdr)

        fits.write(mskIm,header=mskHdr)

        fits.write(wgtIm,header=wgtHdr)
        " and as fits clear extname we will put it back into headers "
        fits[0].write_key('EXTNAME', 'SCI', comment="Extension name")
        fits[1].write_key('EXTNAME', 'MSK', comment="Extension name")
        fits[2].write_key('EXTNAME', 'WGT', comment="Extension name")
        " clean original files "
        if os.path.exists(imName):
            os.remove(imName)
        if os.path.exists(weightName):
            os.remove(weightName)
        if os.path.exists(maskName):
            os.remove(maskName)

    def coadd_assembleF(self,restemp):
        imName = restemp +'_sci.fits'
        weightName = restemp +'_wgt.fits'
        maskName = restemp + '_msk.fits'
        outName = restemp + '.fits'
        commandN = ['coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
        commandN +=[ '--msk_file', '%s' % maskName,  '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
        commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% self.tilename]
        commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
        if string.find(restemp,'det') >0:
            commandN = ['coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
            commandN +=[ '--msk_file', '%s' % maskName, '--band','det', '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
            commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% self.tilename]
            commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
        print commandN

        try:
            subprocess.check_output(commandN)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
            
    """ Use swarp to make detection image """        
    def MakeDetImage(self,ilist,wlist,msklist,ra_cent,dec_cent,naxis1,naxis2,restemp):
        imName = restemp+'_det_sci.fits'
        weightName = restemp+'_det_wgt.fits'
        mskName = restemp+'_det_msk.fits'
        tmpImName =  restemp+'_det_tmp_sci.fits'
#
        command = ['swarp',"%s"%ilist]
        command +=['-NTHREADS','1','-c',os.path.join('./','etc','Y3A1_v1_swarp.config')]
        command +=["-PIXELSCALE_TYPE","MANUAL","-PIXEL_SCALE","%f"%self.pixscale]
        command +=["-CENTER_TYPE","MANUAL","-CENTER","%s,%s"%(ra_cent,dec_cent)]
        command +=['-IMAGE_SIZE',"%d,%d"%(naxis1,naxis2)]
        command +=["-RESAMPLE","N","-BLANK_BADPIXELS","Y"]
        command +=["-COMBINE_TYPE","CHI-MEAN"] 
        command +=["-IMAGEOUT_NAME",imName]
        command +=["-WEIGHTOUT_NAME",weightName]
        command +=["-WEIGHT_IMAGE","%s"%wlist]
        command +=["-COPY_KEYWORDS","BUNIT,TILENAME,TILEID"]
        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
        " Now compose detection mask "    
        command = ['swarp',"%s"%ilist]
        command +=['-NTHREADS','1','-c',os.path.join('./','etc','Y3A1_v1_swarp.config')]       
        command +=["-PIXEL_SCALE","%f"%self.pixscale]
        command +=["-CENTER","%f,%f"%(ra_cent,dec_cent)]
        command +=['-IMAGE_SIZE',"%d,%d"%(naxis1,naxis2)]
        command +=["-RESAMPLE","N","-BLANK_BADPIXELS","Y"]
        command +=["-COMBINE_TYPE","CHI-MEAN"] 
        command +=["-IMAGEOUT_NAME",tmpImName]
        command +=["-WEIGHTOUT_NAME",mskName]
        command +=["-WEIGHT_IMAGE","%s"%msklist]
        command +=["-COPY_KEYWORDS","BUNIT,TILENAME,TILEID"] 
        print command
        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
        if os.path.exists(tmpImName):
            os.remove(tmpImName)  
            
    " Clean coadd files we do not need any more " 
    def cleanCoadds(self,coaddir):
        scifiles = glob.glob(coaddir+'/*_sci.fits')
        for sciF in scifiles:
            os.remove(sciF)
        mskfiles = glob.glob(coaddir+'/*_msk.fits')
        for mskF in mskfiles:
            os.remove(mskF)
        wgtfiles = glob.glob(coaddir+'/*_wgt.fits')
        for wgtF in wgtfiles:
            os.remove(wgtF)
                    
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
        mofdir = datadir + '/mof'
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
         
    def getMofPars(self,datadir):
        mofdir = datadir + '/mof'
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
        mofdir = datadir+'/mof'
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
            command = ['desmeds-prep-tile',self.medsconf,self.tilename,band]
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

    " Prepare data for meds run on simulated images "
    def prepInData(self):
        " make list of realizations "
        realizationslist = os.listdir(self.simData)
        " each realization contains tilename subrirectory "
        for realD in realizationslist:
            " new meds dir is simData+/realization+/medsconf+/ tilename"
            self.realDir = os.path.abspath(self.simData+'/'+str(realD)+'/'+self.medsconf+'/'+self.tilename)
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
 
            " now   nullweight images are prepared in injection step"
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
#        shutil.move(tempFile, fileList)

                
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

        self.coadd_assembleF(restemp)
    
    " clean chunk file in the mof directory "    
    def cleanChunks(self,datadir):
        chunklist = glob.glob(datadir+'/mof/*chunk*.fits')
        for chunkF in chunklist:
            os.remove(chunkF)

                  
if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    curdir = os.getcwd()
    confile = curdir + '/desmesd-config/meds-y3v02.yaml'
    tilename = 'DES0239+0126'
    balP = BalrogUtils(confile,tilename)
    balP.prepMeds()