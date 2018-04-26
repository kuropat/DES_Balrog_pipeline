#!/usr/bin/env python
""" The program to organize all steps necessary to produce
 meds and sof files for Balrog project
 
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
import config_ea as config_mod
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

""" The method to make coadded images for sextractor and meds"""
def makeCoadd(inpar):
    args = inpar[0]
    band = inpar[1]
    ra_cent = args['ra_cent']
    dec_cent = args['dec_cent']
    naxis1 = args['naxis1']
    naxis2 = args['naxis2']
    pixscale = args['pixscale']
   
    outpath = args['outpath']
    logpath = args['logpath']
    tilename = args['tilename']
    print "Make coadd band=%s \n" % band
    restemp = outpath+'/'+tilename+'_'+band
    (images,weights,masks,fluxes) = create_swarp_lists(args,band)
    SWARPcaller(images,weights,fluxes,masks,ra_cent,dec_cent,naxis1,naxis2,pixscale,restemp)
#        byband[band] = (images,weights,fluxes)
    coadd_assemble(tilename,restemp)     
    return

" -----------------  SWARP ----------------------- "
def SWARPcaller(ilist,wlist,msklist,flist,ra_cent,dec_cent,naxis1,naxis2,pixscale,restemp):
    imName = restemp +'_sci.fits'
    weightName = restemp +'_wgt.fits'
    maskName = restemp +'_msk.fits'
    tmpImName = restemp +'_tmp_sci.fits'
    print "imName=%s weightName=%s \n" % (imName,weightName)
    tokens= ilist.split(',')
    tokens1 = wlist.split(',')
    tokens2 = flist.split(',')
    print "input images %d weights %d scales %d \n" % (len(tokens),len(tokens1),len(tokens2))
    print "images %s \n" % ilist
    print "weights %s \n" % wlist
    print "scales %s \n" % flist        
    command = ['swarp',"@%s"%ilist]
    command +=['-NTHREADS','1','-c',os.path.join('./','etc','Y3A1_v1_swarp.config')]
        
    command +=["-PIXEL_SCALE","%f" % pixscale]
    command +=["-CENTER","%f,%f"%(ra_cent,dec_cent)]
    command +=['-IMAGE_SIZE',"%d,%d"%(naxis1,naxis2)]
    command +=["-BLANK_BADPIXELS","Y"]
    command +=["-BLANK_BADPIXELS","Y"]
    command +=["-DELETE_TMPFILES","Y"]
    command +=["-COMBINE","Y","-COMBINE_TYPE","WEIGHTED"] 
    command +=["-IMAGEOUT_NAME",imName]
    command +=["-WEIGHTOUT_NAME",weightName]
    command +=["-FSCALE_DEFAULT","@%s"%flist]
    command +=["-WEIGHT_IMAGE","@%s"%wlist]
    command +=["-COPY_KEYWORDS","BUNIT,TILENAME,TILEID"]

    try:
        subprocess.check_output(command)  
    except subprocess.CalledProcessError as e:
        print "error %s"% e
    " Now make a mask image "        
    command = ['swarp',"@%s"%ilist]
    command +=['-NTHREADS','1','-c',os.path.join('./','etc','Y3A1_v1_swarp.config')]       
    command +=["-PIXEL_SCALE","%f"%pixscale]
    command +=["-CENTER","%f,%f"%(ra_cent,dec_cent)]
    command +=['-IMAGE_SIZE',"%d,%d"%(naxis1,naxis2)]
    command +=["-BLANK_BADPIXELS","Y"]
    command +=["-BLANK_BADPIXELS","Y"]
    command +=["-COMBINE_TYPE","WEIGHTED"]
    command +=["-IMAGEOUT_NAME",tmpImName]
    command +=["-WEIGHTOUT_NAME",maskName]
    command +=["-FSCALE_DEFAULT","@%s"%flist]
    command +=["-WEIGHT_IMAGE","@%s"%msklist]
    command +=["-COPY_KEYWORDS","BUNIT,TILENAME,TILEID"]

    try:
        subprocess.check_output(command)  
    except subprocess.CalledProcessError as e:
        print "error %s"% e
    if os.path.exists(tmpImName):
        os.remove(tmpImName)

" ---------- coadd assemble like in DESDM ---------- "
def coadd_assemble(tilename,restemp):
    imName = restemp +'_sci.fits'
    weightName = restemp +'_wgt.fits'
    maskName = restemp + '_msk.fits'
    outName = restemp + '.fits'
    commandN = ['./bin/coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
    commandN +=[ '--msk_file', '%s' % maskName,  '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
    commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% tilename]
    commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
    if string.find(restemp,'det') >0:
        commandN = ['./bin/coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
        commandN +=[ '--msk_file', '%s' % maskName, '--band','det', '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
        commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% tilename]
        commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
    print commandN

    try:
        subprocess.check_output(commandN)  
    except subprocess.CalledProcessError as e:
        print "error %s"% e
        
""" -------- create list of files for SWARP 
need to be modified after injection will be implemented  ------------ """
def create_swarp_lists(args,band):

    print " make swarp list band=%s \n" % band
    medsdir = args['medsdir']
    medsconf = args['medsconf']
    tilename = args['tilename']
    magbase = args['magbase']
    outpath = args['outpath']
    datadir = args['datadir']
    outpathL = args['listpath']

    flistP = str(datadir)+'/lists/'
    filename = tilename+'_'+band+'_nullwt-flist-'+medsconf+'.dat'
    flistF = flistP+filename

    imext = "[0]"
    mskext = "[1]"
    wext = "[2]"
    list_imF = outpathL+'/'+tilename+'_'+band+'_im.list'
    list_wtF = outpathL+'/'+tilename+'_'+band+'_wt.list'
    list_mskF = outpathL+'/'+tilename+'_'+band+'_msk.list'      
    list_fscF = outpathL+'/'+tilename+'_'+band+'_fsc.list'
    fim = open(list_imF,'w')
    fwt = open(list_wtF,'w')
    scf = open(list_fscF,'w')
    mskf = open(list_mskF,'w')
    first = True
    for line in open(flistF,'r'):
        fname = line.split()[0]
        zerp = float(line.split()[1])
        if first:
            flxscale      = 10.0**(0.4*(magbase - zerp))
            fim.write(str(fname+imext)+'\n')
            fwt.write(str(fname+wext)+'\n')
            mskf.write(str(fname+mskext)+'\n')
            scf.write(str(flxscale)+'\n')
            first = False
        else:
            flxscale      = 10.0**(0.4*(magbase - zerp))
            fim.write(str(fname+imext)+'\n')
            fwt.write(str(fname+wext)+'\n')
            mskf.write(str(fname+mskext)+'\n')
            scf.write(str(flxscale)+'\n')
    fim.close()
    fwt.close()
    mskf.close()
    scf.close()
    return (list_imF,list_wtF,list_fscF,list_mskF)        

" ---------------- run s-extractor to create objects catalog --------- "        
def makeCatalog(inpar):
#    BalP, band = args
#    BalP.makeCatalog(band)
    args = inpar[0]
    band = inpar[1]
    ra_cent = args['ra_cent']
    dec_cent = args['dec_cent']
    naxis1 = args['naxis1']
    naxis2 = args['naxis2']
    pixscale = args['pixscale']
   
    outpath = args['outpath']
    logpath = args['logpath']
    tilename = args['tilename']
    restemp = outpath+'/'+tilename
    logfile = logpath+'/'+tilename    
    SEXcaller(restemp,band,logfile)
    
" -------------SEX ------------------- "    
def SEXcaller(restemp,band,logtemp):
        logFile = logtemp+'_'+band+'_sextractor.log'

        flog = open(logFile,'w')
        image = restemp+'_det.fits[0]'+','+restemp+'_'+band+'.fits[0]'
        flagim = restemp+'_'+band+'.fits[1]'
        weight = restemp+'_det.fits[2]'+','+restemp+'_'+band+'.fits[2]'
        catName = restemp+'_'+band+'_cat.fits'
        command=['sex',"%s"%image,'-c',os.path.join('./','etc','Y3A2_v2_sex.config')]
        command+=['-WEIGHT_IMAGE',"%s"%weight]
        command+=['-CATALOG_NAME',"%s"%catName]
        command+=['-MAG_ZEROPOINT', '30', '-DEBLEND_MINCONT', '0.001', '-DETECT_THRESH', '1.1',  '-ANALYSIS_THRESH', '1.1']
        command+=['-CHECKIMAGE_TYPE', 'SEGMENTATION', '-CHECKIMAGE_NAME', restemp+'_'+band+'_segmap.fits']
        command+=['-FLAG_IMAGE', "%s"%flagim]
#
        command+=['-WEIGHT_TYPE','MAP_WEIGHT']
        command+=['-PARAMETERS_NAME', os.path.join('./','etc','balrog_sex.param')]
        command+=['-FILTER_NAME', os.path.join('./','etc','gauss_3.0_7x7.conv')]
        
           
        print "command= %s \n" % command
        try:

            output,error = subprocess.Popen(command,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

        except subprocess.CalledProcessError as e:
            print "error %s"% e
        flog.write(output)
        flog.write(error)    


        flog.close()   

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
    seedlist = args['seedlist']

#    print seedlist
    seedN = seedlist[cnum-1]
    print seedN
    command = ['run_ngmixit','--nranges', '16',  '--wrange', '1' ,'--tilename', tilename]
    command += ['--meds_list',medslist,'--bands','g,r,i,z']
    command += ['--psf-map', psfmapF]
    command += ['--nbrs-file',mofdir+'/'+tilename+'_nbrslist.fits']
    command += [mofconf,mofdir+'/'+tilename+'_sof-chunk-01.fits','--seed','%d' %seedN]
    command[4] = str(cnum)
    wrange=int(cnum)
    nranges=16 
    command[16] = mofdir+'/'+tilename+'_sof-chunk-'+('%02d'%cnum)+'.fits' 
#
    print command
    res=''
    try:
        res=subprocess.check_output(command)
    except:
        print "failed to run run_ngmixit \n"
    print res   

class BalrogPipelineSof():
    
    def __init__(self, confile, tilename, mof_conf):
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
        desfile = os.getenv("DES_SERVICES")
        if not desfile: desfile = os.path.join(os.getenv("HOME"), ".desservices.ini")
#
        self.desconfig = config_mod.get_desconfig(desfile, self.dbname)
        self.connectDB()
        self.tileinfo = self.get_tile_infor()
        self.tiledir = self.medsdir+'/'+self.medsconf+'/'+self.tilename
        self.mofdir = self.tiledir +'/sof'
        if not os.path.exists(self.mofdir):
            os.makedirs(self.mofdir)
        self.realDir = ''
        self.curdir = os.getcwd()
#        self.realisationslist  = os.listdir(self.simData)
        self.outpath = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename+'/coadd/')
        self.logpath = shutil.abspath(self.outpath)+'/LOGS/'

    def getRealisations(self):
        self.realisationslist  = os.listdir(self.simData)
        return self.realisationslist
            
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

    def create_swarp_lists(self,band,outpath,listpath):
#        outpath = './coadd/'+self.tilename+'/lists/'
        
        outpathL = listpath
        if not os.path.exists(outpathL):
            os.makedirs(outpathL)
        flistP = outpath.split('coadd')[0] +'lists/'
        filename = self.tilename+'_'+band+'_nullwt-flist-'+self.medsconf+'.dat'
        flistF = flistP+filename
        print "Flist file = %s \n" % flistF
        imext = "[0]"
        mskext = "[1]"
        wext = "[2]"
        list_imF = outpathL+self.tilename+'_'+band+'_im.list'
        list_wtF = outpathL+self.tilename+'_'+band+'_wt.list'
        list_mskF = outpathL+self.tilename+'_'+band+'_msk.list'      
        list_fscF = outpathL+self.tilename+'_'+band+'_fsc.list'
        fim = open(list_imF,'w')
        fwt = open(list_wtF,'w')
        scf = open(list_fscF,'w')
        mskf = open(list_mskF,'w')
        first = True
        for line in open(flistF,'r'):
            fname = line.split()[0]
            zerp = float(line.split()[1])
            if first:
                flxscale      = 10.0**(0.4*(self.magbase - zerp))
                fim.write(str(fname+imext)+'\n')
                fwt.write(str(fname+wext)+'\n')
                mskf.write(str(fname+mskext)+'\n')
                scf.write(str(flxscale)+'\n')
                first = False
            else:
                flxscale      = 10.0**(0.4*(self.magbase - zerp))
                fim.write(str(fname+imext)+'\n')
                fwt.write(str(fname+wext)+'\n')
                mskf.write(str(fname+mskext)+'\n')
                scf.write(str(flxscale)+'\n')
        fim.close()
        fwt.close()
        mskf.close()
        scf.close()
        return (list_imF,list_wtF,list_fscF,list_mskF)

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
     
    def makeCoadd(self,band,outpath,tilepath):
        logpath = outpath+'/LOGS/'
        listpath = outpath+'/lists/'
        ra_cent = self.tileinfo['RA_CENT']
        dec_cent = self.tileinfo['DEC_CENT']
        naxis1 = self.tileinfo['NAXIS1']
        naxis2 = self.tileinfo['NAXIS2']
        self.pixscale = self.tileinfo['PIXELSCALE']
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        if not os.path.exists(logpath):
            os.makedirs(logpath)
        " Create coadd images for given band "

        restemp = outpath+'/'+self.tilename+'_'+band
        (images,weights,masks,fluxes) = self.create_swarp_lists(band,outpath,listpath)
        self.SWARPcaller(images,weights,fluxes,masks,ra_cent,dec_cent,naxis1,naxis2,restemp)
#        byband[band] = (images,weights,fluxes)
        self.coadd_assemble(restemp)     
        
    """ call swarp to create a stamp if more than one ccd """
    def SWARPcaller(self,ilist,wlist,msklist,flist,ra_cent,dec_cent,naxis1,naxis2,restemp):
        imName = restemp +'_sci.fits'
        weightName = restemp +'_wgt.fits'
        maskName = restemp +'_msk.fits'
        tmpImName = restemp +'_tmp_sci.fits'
        print "imName=%s weightName=%s \n" % (imName,weightName)
        tokens= ilist.split(',')
        tokens1 = wlist.split(',')
        tokens2 = flist.split(',')
        print "input images %d weights %d scales %d \n" % (len(tokens),len(tokens1),len(tokens2))
        print "images %s \n" % ilist
        print "weights %s \n" % wlist
        print "scales %s \n" % flist        
        command = ['swarp',"@%s"%ilist]
        command +=['-NTHREADS','1','-c',os.path.join('./','etc','Y3A1_v1_swarp.config')]
        
        command +=["-PIXEL_SCALE","%f"%self.pixscale]
        command +=["-CENTER","%f,%f"%(ra_cent,dec_cent)]
        command +=['-IMAGE_SIZE',"%d,%d"%(naxis1,naxis2)]
        command +=["-BLANK_BADPIXELS","Y"]
        command +=["-BLANK_BADPIXELS","Y"]
        command +=["-DELETE_TMPFILES","Y"]
        command +=["-COMBINE","Y","-COMBINE_TYPE","WEIGHTED"] 
        command +=["-IMAGEOUT_NAME",imName]
        command +=["-WEIGHTOUT_NAME",weightName]
        command +=["-FSCALE_DEFAULT","@%s"%flist]
        command +=["-WEIGHT_IMAGE","@%s"%wlist]
        command +=["-COPY_KEYWORDS","BUNIT,TILENAME,TILEID"]

        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
        " Now make a mask image "        
        command = ['swarp',"@%s"%ilist]
        command +=['-NTHREADS','1','-c',os.path.join('./','etc','Y3A1_v1_swarp.config')]       
        command +=["-PIXEL_SCALE","%f"%self.pixscale]
        command +=["-CENTER","%f,%f"%(ra_cent,dec_cent)]
        command +=['-IMAGE_SIZE',"%d,%d"%(naxis1,naxis2)]
        command +=["-BLANK_BADPIXELS","Y"]
        command +=["-BLANK_BADPIXELS","Y"]
        command +=["-COMBINE_TYPE","WEIGHTED"]
        command +=["-IMAGEOUT_NAME",tmpImName]
        command +=["-WEIGHTOUT_NAME",maskName]
        command +=["-FSCALE_DEFAULT","@%s"%flist]
        command +=["-WEIGHT_IMAGE","@%s"%msklist]
        command +=["-COPY_KEYWORDS","BUNIT,TILENAME,TILEID"]

        try:
            subprocess.check_output(command)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e
        if os.path.exists(tmpImName):
            os.remove(tmpImName)
            
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

    def coadd_assemble(self,restemp):
        imName = restemp +'_sci.fits'
        weightName = restemp +'_wgt.fits'
        maskName = restemp + '_msk.fits'
        outName = restemp + '.fits'
        commandN = ['./bin/coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
        commandN +=[ '--msk_file', '%s' % maskName,  '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
        commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% self.tilename]
        commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
        if string.find(restemp,'det') >0:
            commandN = ['./bin/coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
            commandN +=[ '--msk_file', '%s' % maskName, '--band','det', '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
            commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% self.tilename]
            commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
        print commandN

        try:
            subprocess.check_output(commandN)  
        except subprocess.CalledProcessError as e:
            print "error %s"% e


    def SEXcaller(self,restemp,band,logtemp):
        logFile = logtemp+'_'+band+'_sextractor.log'

        flog = open(logFile,'w')
        image = restemp+'_det.fits[0]'+','+restemp+'_'+band+'.fits[0]'
#        flagim = restemp+'_det.fits[1]'+','+restemp+'_'+band+'.fits[1]'
        flagim = restemp+'_'+band+'.fits[1]'
        weight = restemp+'_det.fits[2]'+','+restemp+'_'+band+'.fits[2]'
        catName = restemp+'_'+band+'_cat.fits'
        command=['sex',"%s"%image,'-c',os.path.join('./','etc','Y3A2_v2_sex.config')]
        command+=['-WEIGHT_IMAGE',"%s"%weight]
        command+=['-CATALOG_NAME',"%s"%catName]
        command+=['-MAG_ZEROPOINT', '30', '-DEBLEND_MINCONT', '0.001', '-DETECT_THRESH', '1.1',  '-ANALYSIS_THRESH', '1.1']
        command+=['-CHECKIMAGE_TYPE', 'SEGMENTATION', '-CHECKIMAGE_NAME', restemp+'_'+band+'_segmap.fits']
        command+=['-FLAG_IMAGE', "%s"%flagim]
#
        command+=['-WEIGHT_TYPE','MAP_WEIGHT']
        command+=['-PARAMETERS_NAME', os.path.join('./','etc','balrog_sex.param')]
        command+=['-FILTER_NAME', os.path.join('./','etc','gauss_3.0_7x7.conv')]
        
           
        print "command= %s \n" % command
        try:

            output,error = subprocess.Popen(command,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

        except subprocess.CalledProcessError as e:
            print "error %s"% e
        flog.write(output)
        flog.write(error)    


        flog.close()
        
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
        mofdir = datadir + '/sof'
        medslF = mofdir+'/meds_list.txt'
        outF = open(medslF,'w')
        listF = glob.glob(datadir+"/*meds*.fits.fz")
        for line in listF:
            filename = line.split('/')[-1]
            bandN = filename.split('_')[1]
            outF.write('%s %s \n' % (self.tiledir+'/'+filename,bandN))
        return medslF
     
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
            
    def make_meds_list(self,datadir):
        mofdir=datadir+'/sof'
        medslF = mofdir+'/meds_list.txt'
        outF = open(medslF,'w')
        listF = glob.glob(datadir+"/*meds*.fits.fz")
        for line in listF:
            filename = line.split('/')[-1]
            bandN = filename.split('_')[1]
            outF.write('%s %s \n' % (datadir+'/'+filename,bandN))
        return medslF    
    
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
        for realV in self.realisationslist:
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
        " make list of realisations "
        realisationslist = os.listdir(self.simData)
        " each realisation contains tilename subrirectory "
        for realD in realisationslist:
            " new meds dir is simData+/realisatiuon+ tilename"
            self.realDir = os.path.abspath(self.simData+'/'+str(realD)+'/'+self.tilename)
            print " realDir=%s \n" % self.realDir
            " create subdirectories in the realisation directory "
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
                    realF = fileN.split(band+'/')[1]
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
        "Usage: BalrogPipelineSof.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -m <sof conf file>"

        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:t:m:",["confile=","tilename","mofconfile"])
    except getopt.GetoptError:
        print "Usage:BalrogPipelineSof.py <required inputs>"
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
            print "Usage:BalrogPipelineSof.py<required inputs>"
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
        print "Usage:  BalrogPipelineSof.py <required inputs>"
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
     
    saveS = True
    print "save seeds ",saveS
    print "\n"
    

    ncpu = len(balP.bands)#
#
    args= balP.getCoaddPars(datadir)
    pars = [(args, band) for band in balP.bands ]
#    print pars
    

#    pool = Pool(processes=ncpu)
#    pool.map(makeCoadd, pars) 
#    pool.close()
#    pool.join()
    coaddir = datadir+'/coadd'
#    print 'det image data dir = %s \n' % coaddir
#    balP.makeDetectionImage(coaddir)
#
#    balP.cleanCoadds(coaddir)
#    
#    pool = Pool(processes=ncpu)
#    pool.map(makeCatalog,pars)
#    
#    balP.makeObjMaps(datadir)
#    pool.close()
#    pool.join()
#    pool = Pool(processes=ncpu)
#    pool.map(makeMeds,pars)
#    
#    pool.close()
#    pool.join()
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
    args['mofdir'] = datadir+'/sof'
    args['datadir'] = datadir
    args['seedlist'] = seedlist
    
#    balP.make_nbrs_data(datadir)

    pars = [(args, chunks) for chunks in range(1,17) ]
#    print pars
#
#    pool = Pool(processes=16)

#    pool.map(makeChunk, pars) 
      
#    pool.close()
#    pool.join()

#    balP.collate_chunks(datadir)
#    balP.cleanChunks(datadir)

    #
    """ That's the place where star and galaxy injection should be inserted if we
    will run them in per band mode. After injection if new source filew will be created
    we need to make new file lists for coadd """
    " Now do thins for each realisation "
    datadir = balP.tiledir
    basedir = str(balP.medsdir)+'/' + str(balP.medsconf)
#    tilelistF = 'TileList.csv'
#    tilelistS = './TileList.csv'
#    outf=open(tilelistF,'w')
#    outf.write(tilename+'\n')
#    outf.close()
#    confFile = './Balrog-GalSim/config/bal_config.yaml'
#    geom_file = './Balrog-GalSim/inputs/Y3A2_COADDTILE_GEOM.fits'
#    config_dir = './Balrog-GalSim/config/'
#    psf_dir = datadir
#    tile_dir = basedir
#    config_file = 'bal_config.yaml'
   
#    output_dir =  basedir
#    command = "./Balrog-GalSim/balrog/balrog_injection.py -l %s -g %s -t %s -c %s   -o %s -v %s " % (tilelistS,geom_file,tile_dir,config_dir,output_dir,config_file) 
    
#    print command
#    res=''
#    try:
#        retval = subprocess.call(command.split(),stderr=subprocess.STDOUT)
#        
#    except subprocess.CalledProcessError as e:
#        print "error %s"% e
#        print retval
#    print retval
    """ ---------------------------------------------------------"""
    balP.prepInData()
    reallist = balP.getRealisations()
    basedir = str(balP.medsdir)+'/' + str(balP.medsconf)
    for real in reallist:
        datadir = basedir+'/balrog_images/' + str(real)+'/'+balP.tilename
        
        ncpu = len(balP.bands)
#
        args= balP.getCoaddPars(datadir)
        pars = [(args, band) for band in balP.bands ]
 
        pool = Pool(processes=ncpu)
        pool.map(makeCoadd, pars) 
        pool.close()
        pool.join()
#       
        coaddir = datadir+'/coadd'
        balP.makeDetectionImage(coaddir)
#
        balP.cleanCoadds(coaddir)
#    
        pool = Pool(processes=ncpu)
        pool.map(makeCatalog,pars)
#    
        balP.makeObjMaps(datadir)
        pool.close()
        pool.join()
        pool = Pool(processes=ncpu)
        pool.map(makeMeds,pars)
#    
        pool.close()
        pool.join()
        
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
#    balP.run()  
    
        balP.make_nbrs_data(datadir)

        pars = [(args, chunks) for chunks in range(1,17) ]
#        print pars
#
        pool = Pool(processes=16)
        pool.map(makeChunk, pars) 
           
        pool.close()
        pool.join()

        balP.collate_chunks(datadir)
        balP.cleanChunks(datadir)
