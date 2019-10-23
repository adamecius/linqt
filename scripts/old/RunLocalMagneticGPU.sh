#This is script is intented to be use together with PBS Torque in order to run programs in a distributed way.

#!/bin/bash

#Verify Bash Version
echo "Bash version ${BASH_VERSION}..."
#Compile the Softeware
echo "Compile the Software"
ROOT=~/Dropbox/GraphenePerioBC
bash $ROOT/scripts/FullCompileMagnetic.sh
#Simulation Parameter 
#RUN ACTUAL COMMAND
R=5;
PCID=`uname -n`;
JOBID=`echo $1 | sed 's/[^0-9]*//g'`;
GPU="$""GPU";
i="i";
COMMAND="$ROOT/bin/PureGrapheneConductivity"
#Type Simulation
echo "attemping to run simulation: "
SEED=$((529049+$(date +%Y%m%d)))
#20670078
for i in {1..5..1}
  do
        GPU=1
	$COMMAND$PCID 128 1024 6144 0.1 1 $R $PCID-${i}GPU${GPU}-PBS${JOBID} $GPU $SEED
 done
	


