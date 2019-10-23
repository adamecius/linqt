#! /bin/bash SIMULATION PARAMETER
GPU=0
Nx=100
Ny=$Nx
M=700
W=1.5;
RUNNAME=AndersonNx${Nx}Ny${Ny}W${W}
RUNDIR=$(pwd)
R=3;
S=100;
ID=$RANDOM
if [ $GPU == 0 ]; then
	TGPU=""
else
 if [ $GPU == 1 ]; then 
	TGPU="2"
else
	echo "There is no computer with more than 3 gpu in our cluster"
	exit
 fi
fi


#LOCAL SCRIPT
cat > $RUNDIR/scripts/local/RunLocal$ID.sh << !
#!/bin/bash
#COMPILE DE SOFTWARE FOR THE COMMAND
source /opt/intel/bin/compilervars.sh intel64
bash $RUNDIR/scripts/FullCompile.sh AndersonGraphene
#bash \$ROOT/scripts/FullCompileHaldane.sh
#RUN ACTUAL COMMAND
echo "Bash version ${BASH_VERSION}..."
PCID=\`uname -n\`;
COMMAND=$RUNDIR/bin/AndersonGraphene
JOBID=\`echo \$1 | sed 's/[^0-9]*//g'\`;
SEED=\$RANDOM


for i in {1..$S..1}
        do
        GPU=$GPU
                         #Nx  Ny  M   W  R   MachineName, GPU, seed
        \${COMMAND}\$PCID $Nx $Ny $M $W $R \$PCID-\${i}GPU\${GPU}-PBS\${JOBID} $GPU \$SEED

        done
!
chmod +x $RUNDIR/scripts/local/RunLocal$ID.sh

#TORQUE SCRIPT
cat > $RUNDIR/scripts/TorqueRun.sh << !
#PBS -l nodes=1:ppn=1:GPU$TGPU
##PBS -l nodes=1:ppn=1:main
#PBS -N $RUNNAME
#PBS -m abe
#PBS -M jgarcia@if.ufrj.br

# Bookeeping
echo    "Executando no no  : \$HOSTNAME"
echo -n "Data              : "
date
echo    "Job ID            : \$PBS_JOBID"
echo -n "Diretorio         : "
cd $RUNDIR
pwd
$RUNDIR/scripts/local/RunLocal$ID.sh \$PBS_JOBID
rm $RUNDIR/scripts/local/RunLocal$ID.sh 
!

qsub scripts/TorqueRun.sh




