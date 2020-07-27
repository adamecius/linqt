#BASED ON https://github.com/Barthelemy/CppProjectTemplate
cmake_minimum_required(VERSION 2.4.0)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) # project specific cmake dir

#Configuration variables
set(CMAKE_INSTALL_LIBDIR "$ENV{HOME}/.local/usr/lib/")
set(CMAKE_INSTALL_INCLUDEDIR "$ENV{HOME}/.local/usr/include/")

#Relevant package for the projec
set( PYTHON_MODULES 1 ) #TRUE

if( ${PYTHON_MODULES} )
 find_package(PythonInterp 3.2 REQUIRED)
 string( REPLACE "/bin/python" "" PYTHON_USER_ROOT ${PYTHON_EXECUTABLE} )
 set( PYTHON_USER_LIB ${PYTHON_USER_ROOT}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
 set( PYTHON_USER_BIN ${PYTHON_USER_ROOT}/bin )
 set( PYTHON_USER_INC ${PYTHON_USER_ROOT}/include )
endif()



set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib") 
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Set the name of your project here
project("LinQT")

# Set the version number of your project here (format is MAJOR.MINOR.PATCHLEVEL - e.g. 1.0.0)
set(VERSION_MAJOR "2")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

# This sets the include directory for the reference project. This is the -I flag in gcc.
# All the includes should be in this variable To add new directories to the include path,
# just add them in this list (after a space), or add an include_directories statement in a
# subdirectory's CMakeLists.
include_directories(
    ${PROJECT_SOURCE_DIR}/include
)

#enable_testing()
add_subdirectory(src)
add_subdirectory(utilities)
#add_subdirectory(doc)
#add_subdirectory(test)
#add_subdirectory(samples)




# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)
# check if Doxygen is installed
find_package(Doxygen)

if (DOXYGEN_FOUND)
 set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/linqt.dox)
 set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/docs/linqt.dox)
 configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
 message("Doxygen build started")

 # note the option ALL which allows to build the docs together with the application
 add_custom_target( doc_doxygen ALL
     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
     COMMENT "Generating API documentation with Doxygen"
     VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
