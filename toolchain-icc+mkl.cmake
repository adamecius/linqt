set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_COMPILER_VENDOR "intel")

set(CMAKE_C_COMPILER icc)
set(CMAKE_C_FLAGS "-mkl=parallel -O3 --std=c++11 -Wall -Wextra -fopenmp " CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE " ")
set(CMAKE_C_FLAGS_DEBUG "-g -fopenmp")

set(CMAKE_CXX_COMPILER "icpc")
set(CMAKE_CXX_FLAGS "-mkl=parallel -O3 --std=c++11 -Wall -Wextra -fopenmp " CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")

set(INTEL_MKL "TRUE")
set(BLAS_LIB "-liomp5 -lpthread -lm -ldl")

