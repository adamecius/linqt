#!/bin/bash


#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=16
#SBATCH --ntasks=768
#SBATCH --tasks-per-node=48
#SBATCH --cpus-per-task=1      

# You can choose the parallel environment through modules
module purge
module load impi intel mkl gsl fftw


echo "The simulation started in:"
date
SECONDS=0


srun $myExec 

echo "The simulation ended in:"
date
duration=$SECONDS
echo "And last:"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds
elapsed."


