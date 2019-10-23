#!/bin/bash
#COMPILE DE SOFTWARE FOR THE COMMAND
source /opt/intel/bin/compilervars.sh intel64
bash /home/jgarcia/Dropbox/GraphenePerioBC/scripts/FullCompile.sh GrapheneMagneticField
#bash $ROOT/scripts/FullCompileHaldane.sh
#RUN ACTUAL COMMAND
echo "Bash version 4.3.11(1)-release..."
PCID=`uname -n`;
COMMAND=/home/jgarcia/Dropbox/GraphenePerioBC/bin/GrapheneMagneticField
JOBID=`echo $1 | sed 's/[^0-9]*//g'`;
SEED=$RANDOM

#for M in 50 61 75 92 113 139 171 210 257 316 387 475 583 716 878 1078
for M in 3000
 
	do
        GPU=0
			 #Nx  Ny  M   W R   MachineName, GPU, seed
	${COMMAND}$PCID 20 20 $M 0.0 1 10 Magnetic $GPU $SEED

 	done
