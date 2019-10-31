CC=icpc
CFLAG=" -std=c++17 -mkl=parallel "
OFLAG=" -fopenmp -O3 -march=core-avx2 -fma -ftz -fomit-frame-pointer "
CLINK=" -liomp5 -lpthread -lm -ldl"

#CC=clang++
#CFLAG=" -m64 -I${MKLROOT}/include"
#CLINK=" -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"
#OFLAG=" -O3 -march=znver1 -mfma -fvectorize -mfma -mavx2 -m3dnow -fuse-ld=lld  "

INC="-I../../include"
SRC="../../src"
LIB="-lconfig++"

$CC $CFLAG $INC -c $SRC/mkl_sparse_matrix.cpp -o $SRC/mkl_sparse_matrix.o
$CC $CFLAG $INC -c $SRC/sparse_matrix.cpp -o $SRC/sparse_matrix.o
$CC $CGLAG $INC -c $SRC/linear_algebra.cpp -o $SRC/linear_algebra.o
$CC $CGLAG $INC -c $SRC/chebyshev_solver.cpp -o $SRC/chebyshev_solver.o
$CC $CFLAG $INC -c main.cpp -o main.o
$CC $CFLAG $OFLAG -o main main.o $SRC/chebyshev_solver.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o  $CLINK


