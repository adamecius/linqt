CC=icpc
CFLAG=" -std=c++17 -mkl=parallel "
OFLAG=" -fopenmp -O3 -march=core-avx2 -fma -ftz -fomit-frame-pointer "
CLINK=" -liomp5 -lpthread -lm -ldl"

#CC=clang++
#CFLAG=" -m64 -I${MKLROOT}/include"
#CLINK=" -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"
#OFLAG=" -O3 -march=znver1 -mfma -fvectorize -mfma -mavx2 -m3dnow -fuse-ld=lld  "

INC="-Iinclude"
SRC="src"

$CC $CFLAG $OFLAG  kuboBastinFromChebmom.cpp -o kuboBastinFromChebmom 
