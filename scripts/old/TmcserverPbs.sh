#PBS -l nodes=1:ppn=1:GPU
##PBS -N HaldaneNxNy256M2048
#PBS -N MagneticM6144
#PBS -m abe
#PBS -M jgarcia@if.ufrj.br
#
#=========================================================================================
# JobSkel.sh - Esqueleto para Jobs no PBS/Torque
#
# 2006 (C) D.C & P.H
#
# Instituto de Fisica da UFRJ
#
#=========================================================================================

# CUDA
cd $PBS_O_WORKDIR

# Bookeeping
echo    "Executando no no  : $HOSTNAME"
echo -n "Data              : "
date
echo    "Job ID            : $PBS_JOBID"
echo -n "Diretorio         : "
pwd

########INICIALIZACION##########
#~/Dropbox/GraphenePerioBC/scripts/RunLocalHaldane.sh $PBS_JOBID
#~/Dropbox/GraphenePerioBC/scripts/RunLocalHaldaneGPU.sh $PBS_JOBID
#~/Dropbox/GraphenePerioBC/scripts/RunLocalMagnetic.sh $PBS_JOBID
~/Dropbox/GraphenePerioBC/scripts/RunLocalMagneticGPU.sh $PBS_JOBID
