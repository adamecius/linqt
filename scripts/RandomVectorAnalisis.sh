for N in  400 
do
 	for M in 3200  
	do
		for S in 1 5 10 15 20 25 	
#		for S in 1 4 8 12 16 20 24 28	
#		for S in  1 2 4 8 16 27
#		for S in  1 3 6 9 12 15 18 21 24 27 
	 	do
 			Nx=$N
 			Ny=$N
 			python PythonScripts/SelectedMean.py "data/GrapheneMagneticConductivityXXNx${Nx}Ny${Ny}M${M}U0.1B1578R10*" $S
 			python PythonScripts/SelectedMean.py "data/GrapheneMagneticConductivityXYNx${Nx}Ny${Ny}M${M}U0.1B1578R10*" $S
			python PythonScripts/ConductivityTemperature.py Average/GrapheneMagneticConductivityXXNx${Nx}Ny${Ny}M${M}U0.1B1578R10*S${S}*.dat Average/GrapheneMagneticConductivityXYNx${Nx}Ny${Ny}M${M}U0.1B1578R10*S${S}*.dat
	 	done
		python PythonScripts/RandomVectorConvergenceVer2.py "Temperature/GrapheneMagneticConductivityXYNx${Nx}Ny${Ny}M${M}U0.1B1578R10S*T0.0001*"
	done
done
