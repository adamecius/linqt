#!/bin/bash
#COMPILE DE SOFTWARE FOR THE COMMAND
source /opt/intel/bin/compilervars.sh intel64
bash /home/jgarcia/Dropbox/GraphenePerioBC/scripts/FullCompile.sh DosThesis
#bash $ROOT/scripts/FullCompileHaldane.sh
#RUN ACTUAL COMMAND
echo "Bash version 4.2.37(1)-release..."
PCID=`uname -n`;
COMMAND=/home/jgarcia/Dropbox/GraphenePerioBC/bin/DosThesis
JOBID=`echo $1 | sed 's/[^0-9]*//g'`;
SEED=$RANDOM
#for M in {50..1000..50}
#for M in {1050..2000..50}
#for M in {2050..3000..50}
for M in {3050..4000..50}
#for M in {4050..5000..50}
	do
        GPU=0
	i=0;
			 #Nx  Ny  M   W R   MachineName, GPU, seed
	${COMMAND}$PCID 500 500 $M 0.0 100 $PCID-${i}GPU${GPU}-PBS${JOBID} $GPU $SEED

 	done
