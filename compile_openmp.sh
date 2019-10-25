CC=icpc
CFLAG=" -DMKL_ILP64 -I${MKLROOT}/include"
OFLAG=" -fopenmp -O3 -march=core-avx2 -fma -ftz -fomit-frame-pointer "
INC="-Iinclude"
SRC="src"
LIB=" -lconfig++ -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -lgslcblas -lm -ldl"


$CC $INC -c $SRC/mkl_sparse_matrix.cpp -o $SRC/mkl_sparse_matrix.o
$CC $INC -c $SRC/sparse_matrix.cpp -o $SRC/sparse_matrix.o
$CC $INC -c $SRC/linear_algebra.cpp -o $SRC/linear_algebra.o
$CC $INC -c $SRC/chebyshev_solver.cpp -o $SRC/chebyshev_solver.o

$CC $INC -c $SRC/main.cpp -o $SRC/main.o
$CC -o main $SRC/main.o $SRC/chebyshev_solver.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o -mkl=sequential
#$LIB
#$CC $INC $CFLAG $OFLAG $SRC/kuboBastinFromChebmom.cpp -o ~/.local/bin/kuboBastinFromChebmom $LIB
#$CC $INC $CFLAG $OFLAG $SRC/kuboBastinKernelChebmom.cpp -o ~/.local/bin/kuboBastinKernelChebmom $LIB
##$CC $INC $CFLAG $OFLAG $SRC/kuboBastinFromChebmom.cpp -o kuboBastinFromChebmom $LIB




#$CC $INC $CFLAG $OFLAG $SRC/neqvstimesFromKPMmom.cpp -o neqvstimesFromKPMmom $LIB
#$CC $INC $CFLAG $OFLAG $SRC/condFromKPMmom.cpp -o ~/.local/bin/condFromKPMmom $LIB
#$CC $INC $CFLAG $OFLAG $SRC/running_condFromKPMmom.cpp -o ~/.local/bin/running_condFromKPMmom $LIB
#$CC $INC $CFLAG $OFLAG $SRC/NonEqOp_FL_FromKPMmom.cpp -o ~/.local/bin/NonEqOp_FL_FromKPMmom $LIB
#$CC $INC $CFLAG $OFLAG $SRC/running_NonEqOp_FL_FromKPMmom.cpp -o ~/.local/bin/running_NonEqOp_FL_FromKPMmom $LIB
#$CC $INC $CFLAG $OFLAG $SRC/kuboBastinFromKPMmom.cpp -o ~/.local/bin/kuboBastinFromKPMmom $LIB
#$CC $INC $CFLAG $OFLAG $SRC/compute-kpmMom-nonEqOp.cpp -o ~/.local/bin/compute-kpmMom-nonEqOp $LIB
