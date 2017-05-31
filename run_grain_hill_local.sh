#! /usr/bin/env bash
#
# Version of run script for grain hill model that runs it locally, in serial
# mode.
#

echo "run_grain_hill_local.sh here"

SOURCEDIR=$1

echo "  copying files..."
cp $SOURCEDIR/grain_hill_as_class.py $PWD
cp $SOURCEDIR/cts_model.py $PWD
cp $SOURCEDIR/lattice_grain.py $PWD

runcmd="python $SOURCEDIR/grain_hill_dakota_friendly_driver.py"

echo "  starting python script:"
echo $runcmd

$runcmd

echo "run_grain_hill_local.sh DONE."

# 
# # Create a working directory on the compute node. Copy the contents of
# # the original PBS working directory to it.
# working=/state/partition1/$PBS_JOBNAME-$PBS_JOBID #${TMPDIR}
# if [ ! -d $working ]; then
#     mkdir $working
# fi
# trap "rm -rf $working" EXIT
# cd $working && cp $PBS_O_WORKDIR/* .
# 
# # echo "--> Running on nodes: " `uniq $PBS_NODEFILE`
# # echo "--> Number of available cpus: " $ncpu
# # echo "--> Number of available nodes: " $nnodes
# echo "--> Run command: " $runcmd
# echo "--> Working directory: " $working
# echo "init workdir: " $PBS_O_WORKDIR
# $runcmd
# 
# # Copy the completed run to scratch storage.
# #cp -R $working/* $PBS_O_WORKDIR
# cp -R $working /scratch/gtucker/GrainHill/ParamStudy5x5_Mar2017
# 
