#!/usr/bin/env python
""" The program to run base operations  for Balrog grid production.
This includes following stepe:
1. prepare single epoch data copying them from NCSA
2. Create coadd images and catalogs
3. create mads files for original files.
4. run injection with GalSim
5. create coadd images and catalogs for injected images
6. Create MEDS files for injected images.
 
 By: N. Kuropatkin  04/30/2018
 reorganized by: N. Kuropatkin 08/06/2018
 """
import os
import sys
import string
import getopt
import subprocess
import yaml
import numpy
#
from multiprocessing import Pool

from balrogutils.BalrogUtils import BalrogUtils

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
    (images,weights,fluxes,masks) = create_swarp_lists(args,band)
    SWARPcaller(images,weights,masks,fluxes,ra_cent,dec_cent,naxis1,naxis2,pixscale,restemp)
#
    coadd_assembleF(tilename,restemp)     
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
def coadd_assembleF(tilename,restemp):
    imName = restemp +'_sci.fits'
    weightName = restemp +'_wgt.fits'
    maskName = restemp + '_msk.fits'
    outName = restemp + '.fits'
    commandN = ['coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
    commandN +=[ '--msk_file', '%s' % maskName,  '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
    commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% tilename]
    commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
    if string.find(restemp,'det') >0:
        commandN = ['coadd_assemble', '--sci_file', '%s'% imName,  '--wgt_file', '%s'% weightName ]
        commandN +=[ '--msk_file', '%s' % maskName, '--band','det', '--outname', '%s'% outName, '--xblock', '10',  '--yblock', '3' ]
        commandN +=[ '--maxcols', '100',  '--mincols', '1',  '--no-keep_sci_zeros',  '--magzero', '30',  '--tilename', '%s'% tilename]
        commandN +=[  '--interp_image', 'MSK',  '--ydilate', '3']
    print commandN

    try:
        subprocess.check_output(commandN)  
    except subprocess.CalledProcessError as e:
        print "error %s"% e
        
""" create a list of random numbers to be used as seeds """
def makeSeedList(nchunks):
    arrayV = numpy.random.rand(nchunks)
    seedlist=[]
    for c in range(0,nchunks):
        seedlist.append(int(arrayV[c]*10000))
        
    return seedlist
        
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
            
            


if __name__ == "__main__":
    print sys.argv
    nbpar = len(sys.argv)
    if nbpar < 6:
        "Usage: BalrogBase.py  <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -g <galsimconf> - GalSim configuration file "
        print " -n <number of CPUs to use in GalSim "
        print " -m <mode> - positional code 1 -prep only; 2 - coadd and catalog \n"
        print "    4 - meds for base; 8 - injection; 16 - coadd and catalog for injected; \n"
        print "    32 - meds for injected; 63 - all together; \n"

        sys.exit(-2)

    try:
        opts, args = getopt.getopt(sys.argv[1:],"hc:t:g:n:m:",["confile","tilename","galsimconf","ncpu","mode"])
    except getopt.GetoptError:
        print "Usage: BalrogBase.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -g <galsimconf> - GalSim configuration file "
        print " -n <number of CPUs to use in GalSim "
        print " -m <mode> - positional code 1 -prep only; 2 - coadd and catalog \n"
        print "    4 - meds for base; 8 - injection; 16 - coadd and catalog for injected; \n"
        print "    32 - meds for injected; 63 - all together; \n"
        sys.exit(2)
    c_flag = 0
    t_flag = 0
    g_flag = 0
    m_flag = 0
    n_flag = 0
    ncpu = 4
    for opt,arg in opts:
        print "%s %s"%(opt,arg)
        if opt == "-h":
            print "Usage:  BalrogBase.py <required inputs>"
            print "  Required inputs:"
            print "  -c <confile> - configuration file"
            print " -t <tile> - tile name"
            print " -g <galsimconf> - GalSim configuration file "
            print " -n <number of CPUs to use in GalSim "
            print " -m <mode> - positional code 1 -prep only; 2 - coadd and catalog \n"
            print "    4 - meds for base; 8 - injection; 16 - coadd and catalog for injected; \n"
            print "    32 - meds for injected; 63 - all together; \n"
            sys.exit(2)
           
        elif opt in ("-c","--confile"):
            c_flag = 1
            confile = arg 
        elif opt in ("-t","--tilename"):
            t_flag = 1
            tilename = arg
        elif opt in ("-g","--galsimconf"):
            g_flag = 1
            gconf = arg
        elif opt in ("-m","--mode"):
            m_flag = 1
            mode = int(arg)
        elif opt in ("-n","--ncpu"):
            n_flag = 1
            ncpu = int(arg)
    sumF = c_flag + t_flag + g_flag + m_flag + n_flag
    if sumF != 5:
        print "Usage: BalrogBase.py <required inputs>"
        print "  Required inputs:"
        print "  -c <confile> - configuration file"
        print " -t <tile> - tile name"
        print " -g <galsimconf> - GalSim configuration file "
        print " -n <number of CPUs to use in GalSim "
        print " -m <mode> - positional code 1 -prep only; 2 - coadd and catalog \n"
        print "    4 - meds for base; 8 - injection; 16 - coadd and catalog for injected; \n"
        print "    32 - meds for injected; 63 - all together; \n"
        sys.exit(-2)
    print " Start with Meds  \n"
    """ We start with the main line - create mads  for original tile """
    workdir = os.path.realpath(os.curdir)
    balP = BalrogUtils(confile,tilename)
    " Check if we need prep stage "
    if (mode & 1) == 1:
        print "Start with PrepMeds \n"
        balP.prepMeds()
    datadir = balP.tiledir 
    print " Data dir=%s \n" % datadir
    

    ncpub = len(balP.bands)
#
    args= balP.getCoaddPars(datadir)


    pars = [(args, band) for band in balP.bands ]
    print pars
    " Check if we need coadd "
    if (mode & 2) == 2:
        pool = Pool(processes=ncpub)
        pool.map(makeCoadd, pars) 
        pool.close()
        pool.join()
        coaddir = datadir+'/coadd'
        print 'det image data dir = %s \n' % coaddir
        balP.makeDetectionImage(coaddir)
#
        balP.cleanCoadds(coaddir)
#    
        pool = Pool(processes=ncpub)
        pool.map(makeCatalog,pars)
        pool.close()
        pool.join()
        balP.makeObjMaps(datadir)
        balP.make_psf_map(datadir)
    " Check if we need base meds "
    if (mode & 4 ) == 4:
        print "Start with base meds \n"
        pool = Pool(processes=ncpub)
        pool.map(makeMeds,pars)
#    
        pool.close()
        pool.join()

    realdir = str(balP.medsdir)+'/' + str(balP.medsconf) +'/balrog_images/' 
    if not os.path.exists(realdir):
        os.makedirs(realdir)


    if ( mode & 8) == 8:
        datadir = balP.tiledir
        basedir = str(balP.medsdir)+'/' + str(balP.medsconf)
        tilelistF = 'TileList.csv'
        tilelistS = './TileList.csv'
        outf=open(tilelistF,'w')
        outf.write(tilename+'\n')
        outf.close()
        confFile =  gconf
        geom_file = './inputs/Y3A2_COADDTILE_GEOM.fits'
        config_dir = './Balrog-GalSim/config/'
        psf_dir = datadir+'/psfs'
        tile_dir = basedir
        config_file = 'bal_config.yaml'
   
        output_dir =  basedir
#        command = "python ./Balrog-GalSim/balrog/balrog_injection.py %s -l %s -g %s -t %s -c %s -o %s -n %d -v 1 " % (config_file,tilelistS,geom_file,tile_dir,config_dir,output_dir,ncpu) 
        command = "./balrogutils/bin/RunInject.sh %s %s %s %d " % (basedir,tilename,confFile,ncpu) 
        print command
        print '\n'
        retval=''
        try:
            retval = subprocess.call(command.split(),stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print e.output
            print 'Error running command: \n' 
            print e.cmd  
            print ' \n see above shell error \n'
            print 'Return code: ' + str(e.returncode)
            sys.exit(e.returncode)

        print retval
        print '\n'
        if mode == 8:
            sys.exit(retval)
    "-------------------------------------------------------"
    balP.prepInData()
    reallist = balP.getRealizations()
    basedir = str(balP.medsdir)+'/' + str(balP.medsconf)
    print "basedir=%s \n" % basedir
    for real in reallist:
        datadir = basedir+'/balrog_images/' + str(real)+'/'+ str(balP.medsconf)+'/'+balP.tilename
       
        print "realization=%s \n" % real

        if not os.path.exists(datadir):
            os.makedirs(datadir)
        " Now check if we need coadds for injected images "
        if ( mode & 16) == 16:
            ncpu = len(balP.bands)
            args= balP.getCoaddPars(datadir)
            pars = [(args, band) for band in balP.bands ]
#
            pool = Pool(processes=ncpub)
            pool.map(makeCoadd, pars) 
            pool.close()
            pool.join()
#       
            coaddir = datadir+'/coadd'
            print "Make detection Image \n"
            balP.makeDetectionImage(coaddir)
#
            balP.cleanCoadds(coaddir)
#    
            print " Start with MakeCatalog \n"
            pool = Pool(processes=ncpub)
            pool.map(makeCatalog,pars)
#    
            balP.makeObjMaps(datadir)
            pool.close()
            pool.join()
            " Now check if we need MEDS for injected images "
        if ( mode & 32) == 32:
            ncpu = len(balP.bands)
            args= balP.getCoaddPars(datadir)
            pars = [(args, band) for band in balP.bands ]
            print "run makeMeds with %d cpus \n" % ncpub
            pool = Pool(processes=ncpub)
            pool.map(makeMeds,pars)
#    
            pool.close()
            pool.join()
