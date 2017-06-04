#! /usr/bin/env bash
# A `qsub` submission script that runs GrainHill 
#
## Send email when the job is aborted, started, or stopped
#PBS -m abe

## Send email here
#PBS -M gtucker@colorado.edu

echo "run_grain_hill_beach.sh here"

HOMEDIR=$1
WORKDIR=$2
SCRATCHDIR=$3

echo $HOMEDIR
echo $WORKDIR
echo $SCRATCHDIR

source $HOMEDIR/.bashrc

#echo "  copying files..."
#cp $SOURCEDIR/grain_hill_as_class.py $PWD
#cp $SOURCEDIR/cts_model.py $PWD
#cp $SOURCEDIR/lattice_grain.py $PWD
#cp $DRIVERDIR/grain_hill_dakota_friendly_driver.py $PWD

runcmd="python grain_hill_dakota_friendly_driver.py"

# Create a working directory on the compute node. Copy the contents of
# the original PBS working directory to it.
working=$WORKDIR/$PBS_JOBNAME-$PBS_JOBID
if [ ! -d $working ]; then
    mkdir $working
fi
trap "rm -rf $working" EXIT
cd $working && cp $PBS_O_WORKDIR/* .

# echo "--> Running on nodes: " `uniq $PBS_NODEFILE`
# echo "--> Number of available cpus: " $ncpu
# echo "--> Number of available nodes: " $nnodes
echo "--> Run command: " $runcmd
echo "--> Working directory: " $working
echo "init workdir: " $PBS_O_WORKDIR
$runcmd

# Copy the completed run to scratch storage.
#cp -R $working $SCRATCHDIR
