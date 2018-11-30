MYFFTW_LIB=/home/ICN2/jgarcia/.local/fftw/lib
MYFFTW_INC=/home/ICN2/jgarcia/.local/fftw/include

#mpic++ -O3 -I../../shared/include/ -Iinclude src/band-structure.cpp -o band-structure  -llapacke -llapack -lblas 
#mpic++ -O3 -I../../shared/include/ -Iinclude src/density-of-states.cpp -o density-of-states  -llapacke -llapack -lblas 
##mpic++ -O3 -DUSE_DENSE -I../../shared/include/ -Iinclude src/compute-nonEqOp-ED.cpp -o compute-nonEqOp-ED  -llapacke -llapack -lblas 
#mpic++ -funroll-loops -march=native -mfpmath=sse -O3 -DUSE_DENSE -I../../shared/include/ -Iinclude src/compute-nonEqOp-ED-random.cpp -o compute-nonEqOp-ED-random
#mpiicc -O3 src/mpi-fourier.cpp  -I../../shared/include/ -I$MYFFTW_INC -L$MYFFTW_LIB -lgsl -lgslcblas -lfftw3_mpi -lfftw3 -lm -o test
#mpiicc -O3 src/rand_generator_mpi.cpp  -I../../../shared/include/ -I$MYFFTW_INC -L$MYFFTW_LIB -lgsl -lgslcblas -lfftw3_mpi -lfftw3 -lm -o rand_generator_mpi

OLAG="-O3 -funroll-loops -march=native -mfpmath=sse"

icpc -c -Iinclude/ src/parser/parser_util.cpp -o lib/parser_util.o
icpc -c $OLAG -Iinclude/ src/kpm/kpm.cpp -o lib/kpm.o
icpc -c $OLAG -Iinclude/ src/uck_tb_operator/add_k_dep.cpp -o lib/add_k_dep.o
icpc -c $OLAG -Iinclude/ src/uck_tb_operator/read_operator.cpp -o lib/read_operator.o
icpc -c $OLAG -Iinclude/ src/uck_tb_operator/read_lattice.cpp -o lib/read_lattice.o
icpc -c $OLAG -Iinclude/ src/uck_tb_operator/mat_op.cpp -o lib/mat_op.o

mpiicc -c $OLAG -Iinclude/ -I$MYFFTW_INC  src/lattice_fftw3/MPI_init.cpp  -o lib/MPI_init.o 
mpiicc -c $OLAG -Iinclude/ -I$MYFFTW_INC -L$MYFFTW_LIB src/lattice_fftw3/ffw.cpp  -o lib/ffw.o 


mpiicc -O3  -Iinclude/ -I$MYFFTW_INC -L$MYFFTW_LIB src/compute-nonEqOp-KPM+FFT.cpp  -lgsl -lgslcblas -o compute-nonEqOp-KPM+FFT lib/parser_util.o lib/kpm.o lib/add_k_dep.o lib/read_operator.o lib/read_lattice.o lib/mat_op.o lib/MPI_init.o lib/ffw.o -lfftw3_mpi -lfftw3 -lm
