#!/bin/bash
# The scrip should be run from balrog-base directory.
# CALL PARAMETERS - <meds base> derectory where data will be;
#                   <tilename>
#                   <galsim cofig> GalSim configuration file
#                   <ncpu> - number of CPUs to be used by GalSim
#
echo $1
export $HOME=`pwd`
echo $HOME
HOSTNAME=`/bin/hostname`
DATE=`/bin/date +%H%M%S`
medsbase=$1
tilename=$2
gconf=$3
ncpu=$4
source /cvmfs/des.opensciencegrid.org/eeups/startupcachejob21i.sh
export IFDH_CP_MAXRETRIES=2
DESTCACHE="persistent"
setup python 2.7.9+1
setup numpy 1.15.2
setup fitsio 0.9.12rc1
setup  -j ngmix v0.9.5
setup  galsim 1.5.1tmv # Brian's updated builds
setup numpy 1.15.2
setup python 2.7.9+1

#setup extralibs

workdir=`pwd`

#
export BALROG_BASE=${workdir}
echo "BALROG_BASE="${BALROG_BASE}

cd ${BALROG_BASE}
echo "BALROG_BASE="${BALROG_BASE}
export PYTHONPATH=${BALROG_BASE}/Balrog-GalSim/balrog:${PYTHONPATH}
#
export DESMEDS_CONFIG_DIR=${BALROG_BASE}/desmeds-config/

echo "Start with Injection " $DATE > ${BALROG_BASE}/base_${tilename}.log
#
cd ../
mkdir -p ${medsbase}/
#
cd ${BALROG_BASE}

export medsconf="y3v02"
#
# Now edit bal_config.yaml to make absolute path for input data
#
#cd Balrog-GalSim/config
#sed "s(BALROG-BASE(${BALROG_BASE}(g" config_template_COSMOS.yaml > bal_config.yaml
#sed "s(BALROG-BASE(${BALROG_BASE}(g" config_template.yaml > bal_config.yaml
#cat bal_config.yaml >> ${BALROG_BASE}/base_${tilename}.log 2>&1
#cd ${BALROG_BASE}
if [ ! -f  base_${tilename}.log ];
	then
		touch base_${tilename}.log
fi
ls -l  >& base_${tilename}.log 2>&1
#
#
export PYTHONPATH=${BALROG_BASE}/Balrog-GalSim/balrog:${PYTHONPATH}
cd ${BALROG_BASE}
echo " Now in "`pwd`
tile_dir=$medsbase
config_dir="./Balrog-GalSim/config/"
#config_file="bal_config.yaml"
config_file=$gconf
geom_file="./inputs/Y3A2_COADDTILE_GEOM.fits"
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
python ./Balrog-GalSim/balrog/balrog_injection.py $config_file -l $tile_list -g $geom_file -t $tile_dir -c $config_dir -p $psf_dir -o $output_dir -n $ncpu -v 1  >& base_${tilename}.log 

