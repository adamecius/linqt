cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 11)
# Set the name of your project here
set( SUBPROJECT_NAME wannier2sparse)
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


#BUILD FIRST THE WANNIER2SPARSE LIBRARY
set(LIBRARY_NAME ${SUBPROJECT_NAME})
set(LIB_HEADERS
  Version.h
  wannier_parser.hpp
  hopping_list.hpp
  tbmodel.hpp
  wannier2sparse.hpp
)
set(LIB_SOURCES src/wannier_parser.cpp src/hopping_list.cpp src/tbmodel.cpp src/wannier2sparse.cpp)
add_library(${LIBRARY_NAME} STATIC ${LIB_SOURCES})
set_target_properties(${LIBRARY_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(${LIBRARY_NAME} PROPERTIES VERSION ${VERSION})
set_target_properties(${LIBRARY_NAME} PROPERTIES PUBLIC_HEADER include/wannier2sparse.hpp)

#Create the install command for the library
include(GNUInstallDirs)
install(TARGETS ${LIBRARY_NAME}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


#BUILD THE EXECUTABLE USED TO CAST WANNIER2SPARSE FROM COMMAND LINE
set(APPLICATION_NAME main_${SUBPROJECT_NAME} )
set(APP_SOURCES src/main.cpp )
add_executable(${APPLICATION_NAME} ${APP_SOURCES})
target_link_libraries(${APPLICATION_NAME} ${LIBRARY_NAME} )
set_target_properties(${APPLICATION_NAME} PROPERTIES OUTPUT_NAME ${SUBPROJECT_NAME})


#BUIDL THE PYBINDING MODULES
if( ${PYTHON_MODULES} )

 find_package(pybind11 REQUIRED)
 set(PYBIND11_CPP_STANDARD -std=c++11)
 set( PYTHON_MODULE pybind11_${SUBPROJECT_NAME} )
 pybind11_add_module( ${PYTHON_MODULE} src/pybind11_wrapper.cpp)
 target_link_libraries(${PYTHON_MODULE} PRIVATE ${LIBRARY_NAME} )
 set_target_properties(${PYTHON_MODULE} PROPERTIES OUTPUT_NAME ${SUBPROJECT_NAME})

 install(TARGETS ${PYTHON_MODULE}
         ARCHIVE DESTINATION ${PYTHON_USER_LIB}/site-packages/
         LIBRARY DESTINATION ${PYTHON_USER_LIB}/site-packages/
         PUBLIC_HEADER DESTINATION ${PYTHON_USER_INC}/include)

 install(DIRECTORY lib/ DESTINATION ${PYTHON_USER_LIB}
        PATTERN "lib/*.py"
        PERMISSIONS OWNER_EXECUTE OWNER_READ
        GROUP_EXECUTE GROUP_READ)
endif()


# Require dot, treat the other components as optional
#find_package(Doxygen)
#doxygen_add_docs(wan)
#message( ${DOXYGEN_VERSION} )
