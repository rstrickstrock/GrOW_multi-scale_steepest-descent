#!/bin/bash -l

WALLTIME=01:00:00

for WD in $@ ; do
    export OUTPATH=$WD

    let r=r+1
    str="\$NAME$r"
    export NAME=`eval echo $str`

    if test -e $BCROOT/err$NAME; then
       rm -f $BCROOT/err$NAME
    fi

    if test -e $BCROOT/out$NAME; then
       rm -f $BCROOT/out$NAME
    fi
  
    if ! test -e $OUTPATH/terminated.txt; then
       if [[ $BATCH == "y" ]]; then
          $qsub -q $qsub_queue -N $OUTPATH $QSUBOPTIONS -V -e $BCROOT/err$NAME -o $BCROOT/out$NAME -l nodes=1:ppn=$NODES -l vmem=120GB -l walltime=$WALLTIME $BATCHSCRIPT
          echo $qsub -q $qsub_queue -N $OUTPATH $QSUBOPTIONS -V -e $BCROOT/err$NAME -o $BCROOT/out$NAME -l nodes=1:ppn=$NODES -l vmem=120GB -l walltime=$WALLTIME $BATCHSCRIPT
       else
          $sim_prog $OUTPATH/$CONFIGFILENAME
       fi
    fi
    sleep 1
done

