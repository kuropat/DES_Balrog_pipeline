#!/bin/bash
# The scrip should be run from balrog-base directory.
# It expects definite structure be present:
#    Balrog-GalSim/config directory and config file in it.
#    ./Balrog-GalSim/inputs/Y3A2_COADDTILE_GEOM.fits should be preloaded
#    $tile_dir/$tilename/psfs/ directory should contain pfs files,
#     tis is created by running RunBase.sh script in modes befor injection
#     (1+2+4) = 7, See RunBase.sh modes.
#
#     Also file ./TileList.csv will be prepared by RunBase.sh before this script,
#     otherwise it should be created by user.
#
# CALL PARAMETERS - <meds base> derectory where data will be;
#                   <tilename>
#                   <galsim cofig> GalSim configuration file
#                   <ncpu> - number of CPUs to be used by GalSim
#
echo $1
export HOME=`pwd`
echo $HOME
HOSTNAME=`/bin/hostname`
DATE=`/bin/date +%H%M%S`
medsbase=$1
tilename=$2
gconf=$3
ncpu=$4
unset PYTHONPATH
LD_LIBRARY_PATH_SAVE=$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/cvmfs/des.opensciencegrid.org/2015_Q2/eeups/SL6/eups/packages/Linux64/python/2.7.9+1/lib:$LD_LIBRARY_PATH
source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
export IFDH_CP_MAXRETRIES=2
DESTCACHE="persistent"
setup ngmix v1.2 
setup python 2.7.15+1
setup extralibs

workdir=`pwd`

#
export BALROG_BASE=${workdir}
echo "BALROG_BASE="${BALROG_BASE}

cd ${BALROG_BASE}
echo "BALROG_BASE="${BALROG_BASE}
export PYTHONPATH=${BALROG_BASE}/Balrog-GalSim/balrog:${PYTHONPATH}
#
export DESMEDS_CONFIG_DIR=${BALROG_BASE}/desmeds-config/
if [ ! -f ${BALROG_BASE}/base_${tilename}.log ];
	then
		${BALROG_BASE}/touch base_${tilename}.log
fi
echo "Start with Injection " $DATE >> ${BALROG_BASE}/base_${tilename}.log
#
cd ../
mkdir -p ${medsbase}/
#
cd ${BALROG_BASE}

export medsconf="y3v02"

ls -l  >> base_${tilename}.log 2>&1
#
#
export PYTHONPATH=${BALROG_BASE}/Balrog-GalSim/balrog:${PYTHONPATH}
cd ${BALROG_BASE}
echo " Now in "`pwd`
tile_dir=$medsbase
config_dir="./Balrog-GalSim/config/"
config_file=$gconf
geom_file="./Balrog-GalSim/inputs/Y3A2_COADDTILE_GEOM.fits"
output_dir=$tile_dir
psf_dir=$tile_dir"/"$tilename"/psfs"
echo "psf_dir="$psf_dir
#
# the .TileList.csv will be created before the program will be called.
#
tile_list="./TileList.csv"
# Now injection
#
echo "Start injection "
echo "python ./Balrog-GalSim/balrog/balrog_injection.py $config_file -l $tile_list -g $geom_file -t $tile_dir -c $config_dir -p $psf_dir -o $output_dir -n $ncpu -v 1"
python ./Balrog-GalSim/balrog/balrog_injection.py $config_file -l $tile_list -g $geom_file -t $tile_dir -c $config_dir -p $psf_dir -o $output_dir -n $ncpu -v 1  >> base_${tilename}.log
LD_LIBRARY_PATH=$LD_LIBRARY_PATH_SAVE
