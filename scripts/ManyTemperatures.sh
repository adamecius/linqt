for M in 500 600 700 800 1000 1100 1200 1300 1400 1500 1600 1700 1800 2000
#for M in 1800 2000
do
python PythonScripts/ConductivityTemperature.py Average/GrapheneMagneticConductivityXYNx200Ny200M${M}* Average/GrapheneMagneticConductivityXYNx200Ny200M${M}*
done
