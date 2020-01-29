CC=icpc
CFLAG=" -qopenmp -std=c++17 -mkl=parallel "
OFLAG="-O3 -march=core-avx2 -fma -ftz -fomit-frame-pointer "
CLINK=" -lpthread -lm -ldl"

#CC=clang++
#CFLAG=" -m64 -I${MKLROOT}/include"
#CLINK=" -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"
#OFLAG=" -O3 -march=znver1 -mfma -fvectorize -mfma -mavx2 -m3dnow -fuse-ld=lld  "

INC="-DNDEBUG  -I../../include"
INC="  -I../../include"
SRC="../../src"
LIB=""

BUILD_DIR="../.."

$CC $CFLAG $INC $OFLAG -c $SRC/mkl_sparse_matrix.cpp -o $SRC/mkl_sparse_matrix.o
$CC $CFLAG $INC $OFLAG -c $SRC/sparse_matrix.cpp -o $SRC/sparse_matrix.o
$CC $CFLAG $INC $OFLAG -c $SRC/mkl_linear_algebra.cpp -o $SRC/linear_algebra.o
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_solver.cpp -o $SRC/chebyshev_solver.o
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_moments.cpp -o $SRC/chebyshev_moments.o
$CC $CFLAG $INC $OFLAG -c $SRC/chebyshev_vectors.cpp -o $SRC/chebyshev_vectors.o

$CC $CFLAG $INC -c inline_compute-kpm-nonEqOp.cpp  -o inline_compute-kpm-nonEqOp.o
$CC $CFLAG $OFLAG -o $BUILD_DIR/inline_compute-kpm-nonEqOp inline_compute-kpm-nonEqOp.o $SRC/chebyshev_solver.o $SRC/chebyshev_moments.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o  $CLINK

$CC $CFLAG $INC -c inline_compute-kpm-local_nonEqOp.cpp  -o inline_compute-kpm-local_nonEqOp.o
$CC $CFLAG $OFLAG -o $BUILD_DIR/inline_compute-kpm-local_nonEqOp inline_compute-kpm-local_nonEqOp.o $SRC/chebyshev_solver.o $SRC/chebyshev_moments.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o  $CLINK

#$CC $CFLAG $INC -c kuboBastinFromChebmom_parallel.cpp  -o inline_kuboBastinFromChebmom_parallel.o
#$CC $CFLAG $OFLAG -o $BUILD_DIR/inline_kuboBastinFromChebmom_parallel inline_kuboBastinFromChebmom_parallel.o  $SRC/chebyshev_solver.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o  $CLINK

#$CC $CFLAG $INC -c kuboGreenwoodFromChebmom_parallel.cpp  -o inline_kuboGreenwoodFromChebmom_parallel.o
#$CC $CFLAG $OFLAG -o $BUILD_DIR/inline_kuboGreenwoodFromChebmom_parallel inline_kuboGreenwoodFromChebmom_parallel.o  $SRC/chebyshev_solver.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o  $CLINK

rm -f ../*.o *.o
