To compile please use

mkdir build
cd build
cmake -DCMAKE_TOOLCHAIN_FILE=../toolchain-icc+mkl.cmake  ..
cmake --build .

or

cmake -DCMAKE_TOOLCHAIN_FILE=./toolchain-icc+mkl.cmake -S . -B build
cmake --build build -j 

#Requires Eigen
