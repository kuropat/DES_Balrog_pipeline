#!/usr/bin/env python
""" The program to organize all steps necessary to produce
 meds files for Balrog project
 
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
#from pathos.multiprocessing import ProcessingPool as Pool
from multiprocessing import Pool
from Cython.Compiler.PyrexTypes import BuiltinObjectType
from IcoImagePlugin import IcoImageFile
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
        
def create_swarp_lists(args,band):

    print " make swarp list band=%s \n" % band
    medsdir = args['medsdir']
    medsconf = args['medsconf']
    tilename = args['tilename']
    magbase = args['magbase']
    outpathL = os.path.join(medsdir,medsconf+'/'+tilename+'/coadd/lists/')
#    if not os.path.exists(outpathL):
#        os.makedirs(outpathL)
    flistP = os.path.join(medsdir,medsconf+'/'+tilename+'/lists/')
    filename = tilename+'_'+band+'_nullwt-flist-'+medsconf+'.dat'
    flistF = flistP+filename

    imext = "[0]"
    mskext = "[1]"
    wext = "[2]"
    list_imF = outpathL+tilename+'_'+band+'_im.list'
    list_wtF = outpathL+tilename+'_'+band+'_wt.list'
    list_mskF = outpathL+tilename+'_'+band+'_msk.list'      
    list_fscF = outpathL+tilename+'_'+band+'_fsc.list'
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
    restemp = outpath+'/'+tilename
    logfile = logpath+'/'+tilename    
    fname = os.path.join(medsdir,medsconf+'/'+tilename+'/lists/')+tilename+'_'+band+'_fileconf-'+medsconf+'.yaml'
    with open(fname,'r') as fconf:
        conf=yaml.load( fconf)
        print "file conf \n"
        print conf
        objmap_url = os.path.join(medsdir,medsconf+'/'+tilename+'/lists/')+tilename+'_'+band+'_objmap-'+medsconf+'.fits'
        seg_url = outpath+tilename+'_'+band+'_segmap.fits'
        cimage_url = outpath+tilename+'_'+band+'.fits'
        cat_url = outpath+tilename+'_'+band+'_cat.fits'
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
         

class BalrogPipeline():
    
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
        self.mofdir = self.tiledir +'/mof'
        if not os.path.exists(self.mofdir):
            os.makedirs(self.mofdir)
        print "mof_dir=%s \n" % self.mofdir
        self.curdir = os.getcwd()
        self.outpath = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename+'/coadd/')
        self.logpath = shutil.abspath(self.outpath)+'/LOGS/'

            
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

    def create_swarp_lists(self,band):
#        outpath = './coadd/'+self.tilename+'/lists/'
        outpathL = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename+'/coadd/lists/')
        if not os.path.exists(outpathL):
            os.makedirs(outpathL)
        flistP = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename+'/lists/')
        filename = self.tilename+'_'+band+'_nullwt-flist-'+self.medsconf+'.dat'
        flistF = flistP+filename

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
     
    def makeCoadd(self,band):

        ra_cent = self.tileinfo['RA_CENT']
        dec_cent = self.tileinfo['DEC_CENT']
        naxis1 = self.tileinfo['NAXIS1']
        naxis2 = self.tileinfo['NAXIS2']
        self.pixscale = self.tileinfo['PIXELSCALE']
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
        if not os.path.exists(self.logpath):
            os.makedirs(self.logpath)
        " Create coadd images for given band "

        restemp = self.outpath+'/'+self.tilename+'_'+band
        (images,weights,masks,fluxes) = self.create_swarp_lists(band)
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
            
    def CombineFits(self,restemp):
        imName = restemp +'_sci.fits'
        weightName = restemp +'_wgt.fits'
        maskName = restemp + '_msk.fits'
        outName = restemp + '.fits'
        " read sci image "
        sciIm,imHdr = fitsio.read(imName, ext =0, header=True)
        wgtIm,wgtHdr = fitsio.read(weightName, ext=0, header=True)
        mskIm,mskHdr = fitsio.read(maskName, ext=0, header=True)
        i1 = imHdr['NAXIS1']
        i2 = imHdr['NAXIS2']
        " Now write all 3 extensions in the output file "
        fits = fitsio.FITS(outName,'rw')
        
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
            
                    
    " read catalog to find its length and create a fake objmap"
    def make_fake_objmap(self,catfile,fname):
        fitsio.read(catfile)
        self.coadd_cat = fitsio.read(catfile, lower=True)

        # sort just in case, not needed ever AFIK
        q = numpy.argsort(self.coadd_cat['number'])
        self.coadd_cat = self.coadd_cat[q]
        nobj = len(self.coadd_cat) 
        print " Number of objects in the cat = %d \n" % nobj 
        # now write some data
        fits = fitsio.FITS(fname,'rw')

        data = numpy.zeros(nobj, dtype=[('object_number','i4'),('id','i8')])
       
        data['object_number'] = [num for num in range(nobj)]
        data['id'] = [long(num) for num in range(nobj)]

        print("writing objmap:",fname)
        fitsio.write(fname, data, extname='OBJECTS',clobber=True)  

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
     

    def collate_chunks(self):
        self.chunklistF = self.makeChunkList()
        command=['megamix-meds-collate-desdm','--noblind', 'mof-config/run-Y3A1-v4-mof.yaml']
        command += [self.chunklistF,self.mofdir+'/'+self.tilename+'_mof.fits']
        try:
            subprocess.check_output(command)
        except:
            print "failed to collate chunks \n"
            
    def make_meds_list(self):
        self.medslF = self.mofdir+'/meds_list.txt'
        outF = open(self.medslF,'w')
        listF = glob.glob(self.tiledir+"/*meds*.fits.fz")
        for line in listF:
            filename = line.split('/')[-1]
            bandN = filename.split('_')[1]
            outF.write('%s %s \n' % (self.tiledir+'/'+filename,bandN))
        return self.medslF    
    
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
    
    def collate_chunks(self):
        self.chunklistF = self.makeChunkList()
        command=['megamix-meds-collate-desdm','--noblind', 'mof-config/run-Y3A1-v4-mof.yaml']
        command += [self.chunklistF,self.mofdir+'/'+self.tilename+'_mof.fits']
        try:
            subprocess.check_output(command)
        except:
            print "failed to collate chunks \n"
            
                
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
        
    def getCoaddPars(self):
        pard = {}
        pard['ra_cent'] = self.tileinfo['RA_CENT']
        pard['dec_cent'] = self.tileinfo['DEC_CENT']
        pard['naxis1'] = self.tileinfo['NAXIS1']
        pard['naxis2'] = self.tileinfo['NAXIS2']
        self.pixscale = self.tileinfo['PIXELSCALE']
        pard['pixscale'] = self.pixscale
        
        if not os.path.exists(self.outpath):
            os.makedirs(self.outpath)
        if not os.path.exists(self.logpath):
            os.makedirs(self.logpath)
        outpathL = self.outpath+'/lists'

        print "oupathL=%s \n" % outpathL
        if not os.path.exists(outpathL):
            os.makedirs(outpathL)
        
        pard['outpath'] = self.outpath
        pard['logpath'] = self.logpath
        pard['tilename'] = self.tilename
        pard['medsdir'] = self.medsdir
        pard['medsconf'] = self.medsconf
        pard['magbase'] = self.magbase
        pard['confile'] = self.confile
        return pard
    
    def run(self):
#        self.prepMeds()
#        byband = {}
        
        " Create coadd images for each band "
        for band in self.bands:
            self.makeCoadd(band)
        " Now make detection image "
        self.makeDetectionImage()
        " Create s-extractor catalogs "
        for band in self.bands:
            self.makeCatalog(band)
        " Now create fake objectMaps "
        self.makeObjMaps()
        " Now crete meds for all bands "
        for band in self.bands:
            self.makeMeds(band)
            
    def makeObjMaps(self):        
        for band in self.bands:
            fname = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename)+'/lists/'+self.tilename+'_'+band+'_objmap-'+self.medsconf+'.fits'
            catname = self.outpath+self.tilename+'_'+band+'_cat.fits'
            self.make_fake_objmap(catname, fname)
 
    def makeMeds(self,band):
        fname = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename+'/lists/')+self.tilename+'_'+band+'_fileconf-'+self.medsconf+'.yaml'
        with open(fname,'r') as fconf:
            conf=yaml.load( fconf)
        print "file conf \n"
        print conf
        objmap_url = os.path.join(self.medsdir,self.medsconf+'/'+self.tilename+'/lists/')+self.tilename+'_'+band+'_objmap-'+self.medsconf+'.fits'
        seg_url = self.outpath+self.tilename+'_'+band+'_segmap.fits'
        cimage_url = self.outpath+self.tilename+'_'+band+'.fits'
        cat_url = self.outpath+self.tilename+'_'+band+'_cat.fits'
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
              
    def makeDetectionImage(self):   
        print " Images are created start with detection image \n"
        "Now we have coadd images let's make detection one "
        restemp = self.outpath+'/'+self.tilename
        ra_cent = self.tileinfo['RA_CENT']
        dec_cent = self.tileinfo['DEC_CENT']
        naxis1 = self.tileinfo['NAXIS1']
        naxis2 = self.tileinfo['NAXIS2']
        self.pixscale = self.tileinfo['PIXELSCALE']

        (im_list,wt_list,msk_list) = self.create_det_list(restemp)

        restemp = self.outpath+'/'+self.tilename
        self.MakeDetImage(im_list,wt_list,msk_list,ra_cent,dec_cent,naxis1,naxis2,restemp) 
        restemp = self.outpath+'/'+self.tilename+'_det'

        self.coadd_assemble(restemp)
        
        

 
    def makeCatalog(self,band):
        restemp = self.outpath+'/'+self.tilename
        logfile = self.logpath+'/'+self.tilename    
        self.SEXcaller(restemp,band,logfile)
        
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
    print " Start with Meds  \n"
    
    balP = BalrogPipeline(confile,tilename,mofconf)
    balP.prepMeds()
#    balP.run()

    ncpu = len(balP.bands)
    pool = Pool(ncpu)
    args= balP.getCoaddPars()
    pars = [(args, band) for band in balP.bands ]
    print pars
#
    pool = Pool(processes=ncpu)
    pool.map(makeCoadd, pars) 
    
    
    balP.makeDetectionImage()
    
    pool.map(makeCatalog,pars)
    
    balP.makeObjMaps()
    
    pool.map(makeMeds,pars)
    
    pool.close()
    pool.join()
        
    args =balP.getMofPars()
    
#    bal.run()  
    
    balP.make_nbrs_data()

    pars = [(args, chunks) for chunks in range(1,43) ]
    print pars
#
    pool = Pool(processes=16)
    pool.map(makeChunk, pars) 
           
    pool.close()
    pool.join()
    balP.collate_chunks()
    