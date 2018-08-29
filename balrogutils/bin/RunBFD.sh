#!/bin/bash
#
# The script will prepare data and run rundesmultiband.py  prorgam
#
# Warning: By sure that psf files are in the directory cpecified in the script.
#          It is not true when working with realizations. So copy it there before running the program
#

echo $1
echo $HOME
HOSTNAME=`/bin/hostname`
DATE=`/bin/date +%H%M%S`
export medsconf="y3v02"
medsbase=$1  # medsbase point to directory where ${medsconf}/${tilename} structure is
tilename=$2


workdir=`pwd`
echo "workdir=" $workdir

export workdir=`pwd`


echo "Start balrog bfd at " $DATE > bfd_${tilename}.log

#
source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
DESTCACHE="persistent"
setup pyyaml 3.11+2
setup astropy 2.0.4+0
setup galsim 1.5.1tmv
setup ngmixer y3v0.9.5
setup meds esdevel
setup pyyaml 3.11+2
setup scipy 0.14.0+10
setup fitsio 1.0.0rc1+0
setup esutil esdevel
setup pixmappy 1.0
setup bfd 2.0
setup bfdmeds 1.0
#
export BALROG_BASE=`pwd`
export PYTHONPATH=${BALROG_BASE}/balrogutils/python/:${PYTHONPATH}
export PATH=${BALROG_BASE}/balrogutils/bin:$PATH
#
export PYTHONPATH=${BALROG_BASE}/desmeds/python/:$PYTHONPATH
export PATH=${BALROG_BASE}/desmeds/bin/:$PATH
export DESMEDS_CONFIG_DIR=${BALROG_BASE}/desmeds-config/
#
cd ../
mkdir -p ${medsbase}/

cd ${medsbase}
MEDS_DATA=`pwd`
export realdir=${MEDS_DATA}/${medsconf}/${tilename}
mkdir -p ${realdir}


cd ${BALROG_BASE}
echo "BALROG_BASE="${BALROG_BASE}
echo "MEDS_DATA="${MEDS_DATA}
# 

ls -l >>${BALROG_BASE}/bfd_${tilename}.log 2>&1 

#
#
export DESDATA=${BALROG_BASE}/DESDATA



#

cd ${BALROG_BASE}

#
echo "Start with rundesmultiband.py "  >> ${BALROG_BASE}/bfd_${tilename}.log 2>&1

echo "python $BFD_DIR/bin/rundesmultiband.py --targets --infile ${tilename} --outfile ${tilename}_BFDgriz.fits --bands g r i z   --wt_n 4 --wt_sigma 1.0 \ "
echo "--band_weights 0 0.55 0.3 0.15 --bad_pixel_threshold 0.05 --save_bands_separately yes \ "
echo "--pad_factor 2 --psfdir ${realdir}/psfs --dir ${realdir}/ \ "
echo "--num_proc 16  --moffile ${MEDS_DATA}/${medsconf}-mof/output/${tilename}-${medsconf}-mof.fits \ "
echo " -l ${BALROG_BASE}/${tilename}_multiband.log "
 
python $BFD_DIR/bin/rundesmultiband.py --targets --infile ${tilename} \
--outfile ${tilename}_BFDgriz.fits --bands g r i z  \
 --wt_n 4 --wt_sigma 1.0 --band_weights \
0 0.55 0.3 0.15 --bad_pixel_threshold 0.05 --save_bands_separately yes \
--pad_factor 2 --psfdir ${realdir}/psfs --dir ${realdir}/ \
--num_proc 16  --moffile ${MEDS_DATA}/${medsconf}-mof/output/${tilename}-${medsconf}-mof.fits \
-l ${BALROG_BASE}/${tilename}_multiband.log >> ${BALROG_BASE}/bfd_${tilename}.log 2>&1
#
retcod=$?
echo "rundesmultiband.py ret code=" $retcod >> ${BALROG_BASE}/bfd_${tilename}.log 2>&1
if [ $retcod != 0 ]; then
    echo " rundesmultiband.py failed "
 


fi
exit $retcod