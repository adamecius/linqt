cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 11)
message( ${CMAKE_CURRENT_SOURCE_DIR}/include)
# Set the name of your project here
project(wannier2sparse)
# Set the version number of your project here (format is MAJOR.MINOR.PATCHLEVEL - e.g. 1.0.0)
set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

#CHECK REQUIRED DIRECTORIES
if(NOT DEFINED EIGEN_DIR)
set(EIGEN_DIR  "/usr/include/eigen3/")
endif()


include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/
    ${EIGEN_DIR}
)


set( PYTHON_MODULES "yes" )
if( DEFINED $PYTHON_MODULES )
 find_package(pybind11 REQUIRED)
 set(PYBIND11_CPP_STANDARD -std=c++11)
endif()

add_subdirectory(src)
