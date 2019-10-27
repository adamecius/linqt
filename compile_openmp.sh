CC=icpc
CFLAG=" -std=c++17"
OFLAG=" -fopenmp -O3 -march=core-avx2 -fma -ftz -fomit-frame-pointer "
INC="-Iinclude"
SRC="src"
LIB=" -lconfig++ -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -lgslcblas -lm -ldl"


$CC $CFLAG $OFLAG $INC -c $SRC/mkl_sparse_matrix.cpp -o $SRC/mkl_sparse_matrix.o
$CC $CFLAG $OFLAG $INC -c $SRC/sparse_matrix.cpp -o $SRC/sparse_matrix.o
$CC $CGLAG $OFLAG $INC -c $SRC/linear_algebra.cpp -o $SRC/linear_algebra.o
$CC $CGLAG $OFLAG $INC -c $SRC/chebyshev_solver.cpp -o $SRC/chebyshev_solver.o
$CC $CFLAG $OFLAG $INC -c $SRC/main.cpp -o $SRC/main.o
$CC $CFLAG $OFLAG -o main $SRC/main.o $SRC/chebyshev_solver.o $SRC/sparse_matrix.o $SRC/mkl_sparse_matrix.o $SRC/linear_algebra.o -mkl=sequential

