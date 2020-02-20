INC="-DNDEBUG  -I../../include"
SRC="../../src"
BUILD_DIR="../.."

#ICC WITH MKL
CC=icpc
CFLAG="-fopenmp -std=c++11 "
OFLAG="-O3 "
CLINK=" -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl"


#G++ WITH MKL
#CC=g++
#CFLAG="-std=c++11  "
#OFLAG="-O3 -mavx2 -m3dnow -fomit-frame-pointer" # -mfma -mavx2 -m3dnow -fomit-frame-pointer
#CLINK="-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl"

#Compile linear algebra library
$CC $CFLAG $INC $OFLAG -c $SRC/mkl_linear_algebra.cpp -o $SRC/linear_algebra.o

#Compile the sparse matrix library
$CC $CFLAG $INC $OFLAG -c $SRC/sparse_matrix.cpp -o $SRC/sparse_matrix.o
$CC $CFLAG $INC $OFLAG -c $SRC/mkl_sparse_matrix.cpp -o $SRC/mkl_sparse_matrix.o

#Compile chebyshev moments
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_moments.cpp -o $SRC/chebyshev_moments.o
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_moments1D.cpp -o $SRC/chebyshev_moments1D.o
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_moments2D.cpp -o $SRC/chebyshev_moments2D.o
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_vectors.cpp -o $SRC/chebyshev_vector.o

#Compile the solvers
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_solver.cpp -o $SRC/chebyshev_solver.o

#compute main
$CC $CFLAG $INC -c inline_compute-kpm-nonEqOp.cpp  -o $SRC/inline_compute-kpm-nonEqOp.o
$CC $CFLAG $INC -c inline_compute-kpm-local_nonEqOp.cpp  -o $SRC/inline_compute-kpm-local_nonEqOp.o


$CC $CFLAG $OFLAG -o inline_compute-kpm-nonEqOp       $SRC/inline_compute-kpm-nonEqOp.o $SRC/chebyshev_solver.o $SRC/chebyshev_moments.o  $SRC/chebyshev_moments2D.o $SRC/chebyshev_vector.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o $CLINK
$CC $CFLAG $OFLAG -o inline_compute-kpm-local_nonEqOp $SRC/inline_compute-kpm-local_nonEqOp.o $SRC/chebyshev_solver.o $SRC/chebyshev_moments.o  $SRC/chebyshev_moments2D.o $SRC/chebyshev_vector.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o $CLINK


rm -f ../*.o *.o
