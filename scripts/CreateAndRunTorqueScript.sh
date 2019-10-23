#! /bin/bash
#SIMULATION PARAMETER
GPU=0
Nx=100
Ny=$Nx
M=300;
nflux=1;
W=0.2;
RUNNAME=MagneticFieldNx${Nx}Ny${Ny}M${M}W${W}flux$nflux
RUNDIR=/home/jgarcia/Dropbox/GraphenePerioBC/
R=1;
S=1;
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
bash \$ROOT/scripts/FullCompile.sh PureGrapheneConductivity
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
	\${COMMAND}\$PCID $Nx $Ny $M $W  $nflux  $R \$PCID-\${i}GPU\${GPU}-PBS\${JOBID} $GPU \$SEED

 	done
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
cd $RUNDIR
pwd
/home/jgarcia/Dropbox/GraphenePerioBC/scripts/local/RunLocal$ID.sh \$PBS_JOBID
rm /home/jgarcia/Dropbox/GraphenePerioBC/scripts/local/RunLocal$ID.sh 
!

qsub scripts/TorqueRun.sh




