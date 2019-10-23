#PBS -l nodes=1:ppn=1:GPU
##PBS -l nodes=1:ppn=1:main
#PBS -N MomentConvNx20Ny20W0.0
#PBS -m abe
#PBS -M jgarcia@if.ufrj.br

# Bookeeping
echo    "Executando no no  : $HOSTNAME"
echo -n "Data              : "
date
echo    "Job ID            : $PBS_JOBID"
echo -n "Diretorio         : "
cd /home/jgarcia/Dropbox/GraphenePerioBC
pwd
/home/jgarcia/Dropbox/GraphenePerioBC/scripts/local/RunLocal23907.sh $PBS_JOBID
rm /home/jgarcia/Dropbox/GraphenePerioBC/scripts/local/RunLocal23907.sh 
