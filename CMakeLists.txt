#BASED ON https://github.com/Barthelemy/CppProjectTemplate
cmake_minimum_required(VERSION 2.0.0)


set(CMAKE_CXX_COMPILER "icpc")


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) # project specific cmake dir
set(CMAKE_CXX_FLAGS "--std=c++11 -Wall -Wextra ")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?

if( "${CMAKE_CXX_COMPILER}" STREQUAL "icpc")
    message("Using intel icc optimization options")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=core-avx2 -fma -ftz -fomit-frame-pointer -qopenmp")
endif()


# CMake useful variables
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib") 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Set the name of your project here
project("inline_kuboBastinFromChebmom")

# Set the version number of your project here (format is MAJOR.MINOR.PATCHLEVEL - e.g. 1.0.0)
set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# This sets the include directory for the reference project. This is the -I flag in gcc.
# All the includes should be in this variable To add new directories to the include path,
# just add them in this list (after a space), or add an include_directories statement in a
# subdirectory's CMakeLists.
message( ${PROJECT_SOURCE_DIR} )
include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

enable_testing()
add_subdirectory(include)
add_subdirectory(src)
#add_subdirectory(doc)
add_subdirectory(test)
#add_subdirectory(samples)
add_subdirectory(utilities)

