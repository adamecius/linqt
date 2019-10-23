#! /bin/bash
ROOT=$(pwd)
#SIMULATION PARAMETER
GPU=1;
Nx=128;
Ny=128;
M=512;
W=0.1;
U=0.0;
tU=0.0;
ISO=0.4;
p=0.1;
Lx=100;
Ly=100;
R=1;
RUNNAME=GrapheneLatticeISOaNx${Nx}Ny${Ny}M${M}W${W}ISO${ISO}p${p}Lx${Lx}Ly${Ly}
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
echo "Bash version ${BASH_VERSION}..."
PCID=\`uname -n\`;
SRCID=HaldaneMap
bash $ROOT/scripts/FullCompile.sh \$SRCID
COMMAND=$ROOT/bin/\$SRCID
JOBID=\`echo \$1 | sed 's/[^0-9]*//g'\`;
SEED=\$RANDOM
GPU=$GPU
		 #Nx  Ny  M  W  U  tU  ISO  p  R  MachineName, GPU, seed
\${COMMAND}\$PCID $Nx $Ny $M $W $ISO $p $R $Lx $Ly  \${PCID}-GPU\${GPU}-PBS\${JOBID} $GPU \$SEED
!
chmod +x scripts/local/RunLocal$ID.sh

#TORQUE SCRIPT
cat > scripts/TorqueRun.sh << !
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
cd $(pwd)
pwd
scripts/local/RunLocal$ID.sh \$PBS_JOBID
rm scripts/local/RunLocal$ID.sh 
!

qsub scripts/TorqueRun.sh

