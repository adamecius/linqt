#!/bin/bash
#COMPILE DE SOFTWARE FOR THE COMMAND
source /opt/intel/bin/compilervars.sh intel64
bash /home/jgarcia/Dropbox/GraphenePerioBC/scripts/FullCompile.sh DosThesis
#bash $ROOT/scripts/FullCompileHaldane.sh
#RUN ACTUAL COMMAND
echo "Bash version 4.3.11(1)-release..."
PCID=`uname -n`;
COMMAND=/home/jgarcia/Dropbox/GraphenePerioBC/bin/DosThesis
JOBID=`echo $1 | sed 's/[^0-9]*//g'`;
SEED=$RANDOM
for R in 1000
        do
        GPU=1
                         #Nx  Ny  M   W R   MachineName, GPU, seed
        ${COMMAND}$PCID 1000 1000 2000 0.0 $R RandomConvergence0 0 $SEED

        done








