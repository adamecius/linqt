#This script is intended to Run several realizations of my Software

#!/bin/bash

#COMPILE DE SOFTWARE FOR THE COMMAND
ROOT=~/Dropbox/GraphenePerioBC
bash $ROOT/scripts/FullCompile.sh
#RUN ACTUAL COMMAND

echo "Bash version ${BASH_VERSION}..."
R=5;
PCID=`uname -n`;
#COMMAND=$ROOT/bin/PureGrapheneConductivity
COMMAND=$ROOT/bin/HaldaneHoneyComb
JOBID=`echo $1 | sed 's/[^0-9]*//g'`;
SEED=$(($RANDOM*$RANDOM))
for i in {1..10..1}
  do
	GPU=1
	${COMMAND}$PCID 256 256 2048 0.1 0.5 0.4 0 $R 0.90 1 $PCID-${i}GPU${GPU}-PBS${JOBID} $GPU $SEED
	${COMMAND}$PCID 256 256 2048 0.1 0.1 0.0 0 $R 0.90 1 $PCID-${i}GPU${GPU}-PBS${JOBID} $GPU $SEED
	${COMMAND}$PCID 256 256 2048 0.9 0.5 0.4 0 $R 0.90 1 $PCID-${i}GPU${GPU}-PBS${JOBID} $GPU $SEED
	${COMMAND}$PCID 256 256 2048 0.9 0.1 0.0 0 $R 0.90 1 $PCID-${i}GPU${GPU}-PBS${JOBID} $GPU $SEED
#	${COMMAND}$PCID 512 512 4096 0.1 0.1 0.0 0 $R 0.90 1 $PCID-${i}GPU${GPU}-PBS${JOBID} $GPU $SEED

 done


