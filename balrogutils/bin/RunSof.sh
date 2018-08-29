#!/bin/bash
#
# CALL PARAMETERS - <meds base> derectory where data will be;
#                   <tilename>;
#                   <sof configfuration file>;
#                   <ncpu> number of CPUs to use                  
#
echo $1
echo $HOME
HOSTNAME=`/bin/hostname`
DATE=`/bin/date +%H%M%S`
export medsconf="y3v02"
medsbase=$1  # medsbase point to directory where ${medsconf}/${tilename} structure is
tilename=$2
conf=$3
ncpu=$4
source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
setup mofpsfex 0.3.2
setup despyastro 0.3.9+2
setup fitsverify
setup esutil 0.6.2rc1+1
setup pixcorrect 0.5.3+12
setup despydb
setup ngmixer y3v0.9.4b
setup meds 0.9.3+0
setup cfitsio 3.370+0
setup easyaccess
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
export MEDS_DATA=`pwd`
mkdir -p ${medsconf}/${tilename}

cd ${BALROG_BASE}

echo "MEDS_DATA="${MEDS_DATA}
# We have to redifine MEDS_DIR variable as Erin expect it be 
# the place where we put all results
export MEDS_DIR=${MEDS_DATA}
echo "MEDS_DIR="${MEDS_DIR}
mkdir -p ${MEDS_DIR}/temp
export TMPDIR=${MEDS_DIR}/temp
export NGMIXER_OUTPUT_DIR=${MEDS_DIR}
#
export DESDATA=${BALROG_BASE}/DESDATA


cd ${BALROG_BASE}
ls -l  >> sof_${tilename}.log 2>&1
#
#
# Now execute megamix-test-tiles.sh
#
BalrogSofMegamixer.py -c ${conf} -t ${tilename} -n ${ncpu}  >> sof_${tilename}.log 2>&1
echo "Exit status: "$?
exit $?
