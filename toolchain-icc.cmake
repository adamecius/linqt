set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_COMPILER_VENDOR "intel")

set(CMAKE_C_COMPILER icc)
set(CMAKE_C_FLAGS "-O3 --std=c++11 -Wall -Wextra -fopenmp " CACHE STRING "" FORCE)
set(CMAKE_C_FLAGS_RELEASE " ")
set(CMAKE_C_FLAGS_DEBUG "-g -fopenmp")
set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?

set(CMAKE_CXX_COMPILER "icpc")
set(CMAKE_CXX_FLAGS " -O3 --std=c++11 -Wall -Wextra -fopenmp " CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fopenmp")
set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?

set(INTEL_MKL "TRUE")
