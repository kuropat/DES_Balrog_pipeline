# DES_Balrog_pipeline
# Based on Erin Sheldon meds and mof packages
#  Those packages are modified to work with Balrog and are included with
# this pipeline

code specific to des Balrog production

# The code provides a baseline for Balrog production
Implemented steps:
1.Query DESDM database for a tile configuration.
2. Run meds_prep to copy all necessary files from DESDM to $MEDS_DATA directory
3. Create coadded images and catalogs
4. Run meds to create meds files
5. Run mof to create mof file.
6. Generate images with injected objects in <meds_data>/<medsconf>/balrog_images/<realisation>/<tilename> (unimplemented )
7. prepare configurations, filelists and images in each realisation subdirectory to ru the rest of the pipeline.
8. Create coadded images and catalogs for given combination of <realisation>/<tilename>
9. Run meds to create meds files
10. Run mof to create mof file. ( runs complete mof at present time, will chane in future )

# Environment
  See the setup_balrog.sh file to create appropriate environment.
  As the program is using DESDM database one need to have .desservices.ini file
  properly configured in his home directory.

  See the meds configuration file in the meds-config subdirectory.
  One can change a set of bands for which the meds and mof files are created by
  modifying the config file.

The program runs in multiprocessor environment. Presently we have reserved 16
CPUs. It is hardcoded right now and will be a parameter in future.

To run the program:
BalrogPipeline.py -c desmeds-config/Pipeline.py -t <tilename> -m mof-config/run-Y3A1-v4-mof.yaml

