MYFFTW_LIB=/home/ICN2/jgarcia/.local/fftw/lib
MYFFTW_INC=/home/ICN2/jgarcia/.local/fftw/include

OLAG="-O3 -funroll-loops -march=native -mfpmath=sse"

#icpc -c -Iinclude/ src/parser/parser_util.cpp -o lib/parser_util.o
#icpc -c $OLAG -Iinclude/ src/kpm/kpm.cpp -o lib/kpm.o
#icpc -c $OLAG -Iinclude/ src/uck_tb_operator/add_k_dep.cpp -o lib/add_k_dep.o
#icpc -c $OLAG -Iinclude/ src/uck_tb_operator/read_operator.cpp -o lib/read_operator.o
#icpc -c $OLAG -Iinclude/ src/uck_tb_operator/write_operator.cpp -o lib/write_operator.o
#icpc -c $OLAG -Iinclude/ src/uck_tb_operator/read_lattice.cpp -o lib/read_lattice.o
#icpc -c $OLAG -Iinclude/ src/uck_tb_operator/mat_op.cpp -o lib/mat_op.o

#mpiicc -c $OLAG -Iinclude/ -I$MYFFTW_INC  src/lattice_fftw3/MPI_init.cpp  -o lib/MPI_init.o 
#mpiicc -c $OLAG -Iinclude/ -I$MYFFTW_INC -L$MYFFTW_LIB src/lattice_fftw3/ffw.cpp  -o lib/ffw.o 



icpc -O3 -Iinclude/ src/proj-band-structure.cpp -o proj-band-structure lib/parser_util.o lib/add_k_dep.o lib/read_operator.o lib/read_lattice.o -llapack -mkl=sequential
#icpc -O3 -Iinclude/ src/band-structure.cpp -o band-structure lib/parser_util.o lib/add_k_dep.o lib/read_operator.o lib/read_lattice.o -llapack -mkl=sequential

#mpiicc -O3  -Iinclude/ -I$MYFFTW_INC -L$MYFFTW_LIB src/compute-nonEqOp-KPM+FFT.cpp  -lgsl -lgslcblas -o compute-nonEqOp-KPM+FFT lib/parser_util.o lib/kpm.o lib/add_k_dep.o lib/read_operator.o lib/read_lattice.o lib/mat_op.o lib/MPI_init.o lib/ffw.o -lfftw3_mpi -lfftw3 -lm
