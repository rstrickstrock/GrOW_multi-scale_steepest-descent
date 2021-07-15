#!/bin/bash -l

# ---------------------------
# our name 
#PBS -N Job_$T
#
#              
#PBS -q default
#PBS -l nodes=1:ppn=1
#PBS -l mem=5G
#PBS -l walltime=72:00:00
#PBS -j oe
#PBS -r n
#PBS -A FH
#PBS -S /bin/bash
#PBS -V


# give the values of some environment variables

echo "HOSTNAME $HOSTNAME is the host name of the node running the job"
echo "JOB_ID $JOBID"
echo "JOB_NAME $JOBNAME"
echo "PE $PE is the parallel environment"
echo "PE_HOSTFILE $PE_HOSTFILE is the hostfile"
echo "QUEUE $QUEUE is the queue running this job"
echo "MPIRUN: $MPIRUN"
echo "BINDIR: $BINDIR"

NNODES=1
NCORES=$NODES

#module load $MPIRUN
#module load $BINDIR

export DO_PARALLEL="mpirun -np $NCORES "

python $sim_prog $OUTPATH/$CONFIGFILENAME

echo ENDE
date
