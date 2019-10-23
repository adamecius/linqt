for N in 100 300 500 700 1000
do	

 	for M in 10 12 15 19 23 29 36 44 55 68 84 104 129 160 198 245 303 375 464 575 711 880 1090 1349 1669 2066 2557 3165 3917 4848 6000 
 	do
 		Nx=$N
 		Ny=$N
		python PythonScripts/MeanNotnan.py "data/GrapheneDOSNx${Nx}Ny${Ny}M${M}U0.0R100Realization*"
 	done
done
