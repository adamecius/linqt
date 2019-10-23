#! /bin/bash SIMULATION PARAMETER
GPU=0
Nx=1000
Ny=$Nx
M=2500;
W=0.0;
RUNNAME=RandomConvNx${Nx}Ny${Ny}M${M}
RUNDIR=$(pwd)
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
bash $RUNDIR/scripts/FullCompile.sh DosThesis
#bash \$ROOT/scripts/FullCompileHaldane.sh
#RUN ACTUAL COMMAND
echo "Bash version ${BASH_VERSION}..."
PCID=\`uname -n\`;
COMMAND=$RUNDIR/bin/DosThesis
JOBID=\`echo \$1 | sed 's/[^0-9]*//g'\`;
SEED=\$RANDOM
#for R in 1 2 4 5 6 8 10 13 16 20 25 32 40 50 63 79 100 126 158 200 251 316 398 501 631 794 1000
for R in 631 794 1000
	do
        GPU=$GPU
			 #Nx  Ny  M   W R   MachineName, GPU, seed
	\${COMMAND}\$PCID $Nx $Ny $M $W \$R RandomConvergence0 $GPU \$SEED

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




