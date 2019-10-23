

for N in 400
do	
 	for M in 400 800 1600 3200 5000
 	do
 		Nx=$N
 		Ny=$N
# 		python PythonScripts/Mean.py "data/GrapheneMagneticConductivityXXNx${Nx}Ny${Ny}M${M}U0.1B1578R10*"
# 		python PythonScripts/Mean.py "data/GrapheneMagneticConductivityXYNx${Nx}Ny${Ny}M${M}U0.1B1578R10*"
		python PythonScripts/ConductivityTemperature.py Average/GrapheneMagneticConductivityXXNx${Nx}Ny${Ny}M${M}U0.1B1578R10* Average/GrapheneMagneticConductivityXYNx${Nx}Ny${Ny}M${M}U0.1B1578R10*
 	done
#	python PythonScripts/PolynomialConvergence.py "Temperature/GrapheneMagneticConductivityXXNx${Nx}Ny${Ny}M*"
	python PythonScripts/PolynomialConvergence.py "Temperature/GrapheneMagneticConductivityXYNx${Nx}Ny${Ny}M*"	
done
