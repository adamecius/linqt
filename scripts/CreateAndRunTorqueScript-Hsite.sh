
#! /bin/bash

#SIMULATION PARAMETER
GPU=0
Nx=100
Ny=$Nx
M=10000;
nflux=1;
p=0.05;
V0=0.685;
R=10;
S=1000
W=0.1;

RUNNAME=MagneticNx${Nx}Ny${Ny}M${M}W${W}V0${V0}p$p
RUNDIR=/home/jgarcia/Dropbox/GraphenePerioBC/
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
cat > scripts/local/RunLocal$ID.sh << !
#!/bin/bash
#COMPILE DE SOFTWARE FOR THE COMMAND
ROOT=/home/jgarcia/Dropbox/GraphenePerioBC
bash \$ROOT/scripts/FullCompileHsite.sh
#bash \$ROOT/scripts/FullCompileHaldane.sh
#RUN ACTUAL COMMAND
echo "Bash version ${BASH_VERSION}..."
PCID=\`uname -n\`;
COMMAND=\$ROOT/bin/PureGrapheneConductivity
#COMMAND=\$ROOT/bin/HaldaneHoneyComb
JOBID=\`echo \$1 | sed 's/[^0-9]*//g'\`;
SEED=\$((\$RANDOM*\$RANDOM))
SEED=24455532
for i in {1..$S..1}
	do
        GPU=$GPU
			 #Nx  Ny  M   W    Index   R   MachineName, GPU, seed
	\${COMMAND}\$PCID $Nx $Ny $M $W  $V0 $p  $R \$PCID-\${i}GPU\${GPU}-PBS\${JOBID} $GPU \$SEED

 	done
!
chmod +x scripts/local/RunLocal$ID.sh

#TORQUE SCRIPT
cat > scripts/TorqueRun.sh << !
#PBS -l nodes=1:ppn=1:GPU$TGPU
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
/home/jgarcia/Dropbox/GraphenePerioBC/scripts/local/RunLocal$ID.sh \$PBS_JOBID
rm /home/jgarcia/Dropbox/GraphenePerioBC/scripts/local/RunLocal$ID.sh 
!

qsub scripts/TorqueRun.sh




