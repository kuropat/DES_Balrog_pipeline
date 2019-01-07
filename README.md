# DES_Balrog_pipeline
## Based on Erin Sheldon meds and mof and megamixer packages
#### Modified and repacked version. Version 2.0
#### Updated version to work with latest Balrog-GalSim, python 2.7.15+1 and
ngmix v1.2. The code was tested on both SL7 and SL6 hosts.  

###### code specific to des Balrog production is located in __balrogutils__ directory.

The directory contains `bin` and `python` subdirectories.

Executable modules are:

#### BalrogBase.py

 The code provides a baseline for Balrog production 

Implemented steps:

###### Block 1: 

   Query DESDM database for a tile configuration.    
   Run `meds_prep` to copy all necessary files from DESDM to `$MEDS_DATA` directory

###### Block 2: 

   Create coadded images and catalogs from files downloaded by Block 1

###### Block 3: 

   Run desmeds to create meds files on input images

###### Block 4: 

   Generate images with injected objects in
   `<meds_data>/<medsconf>/balrog_images/<realisation>/<medsconf>/<tilename>` 
   prepare configurations, filelists and images in each realization subdirectory   

###### Block 5: 

Create coadded images and catalogs for ech realization.

###### Block 6: 

Run desmeds to create meds files for each realization

###### Input parameters: 

*
          -c <confile>  configuration file like `desmeds-config/meds-y3v02.yaml`
          -t <tile>  tile name like `DES0239+0216`
	  -g <galsimconfig> - name of the configuration file for GalSim, that
             should be in Balrog-GalSim/config/ directory. 
          -n <number of CPUs to use in GalSim> like 16
          -m <mode>  positional code 1 -prep only; 2 - coadd and catalog 
                      4 - meds for base; 8 - generate injected images
                     16 - coadd and catalog for injected images
                     32 - meds for injected images; 63 - all blocks
*

#### RunBase.sh  
An example how to create environment and run the BalrogBase.py

###### input parameters:

*
       <meds base> directory where data will be
       <tilename>
       <galsin config file> - GalSim configuration yaml file
       <mode> positional code:
             1 prep only 
             2 coadd and catalog 
             4 meds for base
             8 injection 
            16 coadd and catalog for injected images 
            32  meds for injected 
            63  all together
       <ncpu>  number of CPUs to be used by GalSim, by defauld the progrm will use
                            ncpu equal the number of bands
*

The program requires that Spencer Everett
[Balrog-GalSim package](https://github.com/sweverett/Balrog-GalSim.git)
 be installed in the base directory, and inputs/ subdirectory with star and
 galaxy catalogs should be present.

To create coadded images and catalogs the program is using configuration
files found in `etc` subdirectory

To run GalSim simulation the program requares the `config_template_COSMOS.yaml`
file be copied to `Balrog-GalSim/config/` subdirectory. And example of the
file is included in the repository.

In general the base directory structure should look like following:
```      
drwxr-xr-x 7 kuropat sdss     155 Aug 10 13:37 Balrog-GalSim
drwxr-xr-x 4 kuropat sdss      41 Aug  6 14:08 balrogutils
drwxr-xr-x 3 kuropat sdss      24 Aug  3 10:49 DESDATA
drwxr-xr-x 4 kuropat sdss      85 Aug 10 17:18 desmeds
drwxr-xr-x 2 kuropat sdss      36 Aug  3 10:49 desmeds-config
drwxr-xr-x 2 kuropat sdss    4096 Aug  3 10:49 etc
drwxr-xr-x 3 kuropat sdss    4096 Aug  3 10:50 inputs
drwxr-xr-x 2 kuropat sdss      70 Aug  8 13:13 mof-config
drwxr-xr-x 2 kuropat sdss     115 Aug  3 10:50 shr-config
drwxr-xr-x 2 kuropat sdss     138 Aug  3 14:25 sof-config
```

#### BalrogSofMegamixer.py
 
The program to create sof files for given data set.
Usage: BalrogSofMegamixer.py  `<required inputs>`
###### Required inputs:

*
        -c <confile>  configuration file like `sof-config/run-y3v02-sof.yaml`
        -t <tile>  tile name like `DES0239+0126`
        -n <number of CPUs to use> like 16
*

Beside of input parameters the program requires following environment 
variables be defined:
*        
        MEDS_DIR - a base directory where input data are in the structure
         `<medsversion>/<tilename>` for example `${MEDS_DIR}/y3v02/DES0239+0126/`
        MEDS_DATA equal to MEDS_DIR  
        BALROG_BASE  the base directory from which the program is running
        medsconf  meds version like y3v02
        NGMIXER_OUTPUT_DIR  is directory where output of the ngmixer will bestored
*

Results will be put in `${NGMIXER_OUTPUT_DIR}/y3v02-sof/` directory


#### RunSof.sh
 
An example script how to create environment and run BalrogSofMegamixer.py
###### input parameters:

* 
         <meds base> directory where data will be. It is place where 
                               `<medsconf>/<tilename>` subdirectories are.
         <tilename> name of the tile
         <mof configfuration file>  like `sof-config/run-y3v02-sof.yaml`
         <ncpu> number of CPUs to use
*
#### BalrogMofMegamixer.py 

A program to create mof using Sheldon's megamixer.
Input parameters are the same as in BalrogSofMegamixer.py

#### RunMof.sh 

An example script how to create environment and run BalrogMofMegamixer.py
###### input parameters:

*
       <meds base> directory where data will be. It is place where `<medsconf>/<tilename>` subdirectories are;
       <tilename>
       <mof configfuration file>  like `mof-config/run-y3v02-mof.yaml`
       <ncpu> number of CPUs to use
*

#### RunBfd.sh
An example script how to run BFD program.
###### input parameters:
*
       <meds base> directory where data will be. It is place where `<medsconf>/<tilename>` subdirectories are;
       <tilename>
*

#### RunShr.sh
An example script how to run metacal program.
###### input parameters:
*
     <meds base> directory where data will be. It is place where `<medsconf>/<tilename>` subdirectories are; 
     <tilename>
     <mcal configfuration file>  like `shr-config/run-y3v02-mcal.yaml` 
     <ncpu> number of CPUs to use
*
#### BalrogShrMegamixer.py 
The program to tun megamixer metacal program on MOF files  produced by BalrogMofMegamixer.py program 
###### input parameters:
*
        -c <confile> - configuration file
        -t <tile> - tile name
        -n <number of CPUs to use> 
*
