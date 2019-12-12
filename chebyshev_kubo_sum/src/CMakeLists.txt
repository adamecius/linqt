# @author  Jose H. Garcia


set(LIB_HEADERS
 "${CMAKE_CURRENT_BINARY_DIR}/../include/Version.h" # Generated header
 "${CMAKE_CURRENT_BINARY_DIR}/../include/chebyshev_coefficients.hpp" # Normal header
 "${CMAKE_CURRENT_BINARY_DIR}/../include/chebyshev_kubo.hpp" # Normal header
 "${CMAKE_CURRENT_BINARY_DIR}/../include/chebyshev_moments.hpp" # Normal header
)

set(LIB_SOURCES
 chebyshev_coefficients.cpp # Normal header
 chebyshev_kubo.cpp # Normal header
 chebyshev_moments.cpp # Normal header
)

set(APP_HEADERS
)

set(APP_SOURCES
  kuboBastinFromChebmom.cpp
)

# Define two variables in order not to repeat ourselves.
set(LIBRARY_NAME chebyshev_kubo_lib)

# Define the library
add_library(${LIBRARY_NAME} STATIC ${LIB_SOURCES})

# Set the build version. It will be used in the name of the lib, with corresponding
# symlinks created. SOVERSION could also be specified for api version. 
set_target_properties(${LIBRARY_NAME} PROPERTIES
  VERSION ${VERSION}  # ${VERSION} was defined in the main CMakeLists.
  FRAMEWORK FALSE
  PUBLIC_HEADER "${LIB_HEADERS}" # specify the public headers
)

# This is a make target, so you can do a "make projArunner"
set(APPLICATION_NAME inline_kuboBastinFromChebmom)
add_executable(${APPLICATION_NAME} ${APP_SOURCES})

# Links the target exe (projArunner) against the library (libprojA) and other required libraries
target_link_libraries(${APPLICATION_NAME} ${LIBRARY_NAME} )

# Says how and where to install software
# Targets:
#   * <prefix>/lib/<libraries>
#   * header location after install: <prefix>/include/<project>/*.h
#   * headers can be included by C++ code `#<project>/Bar.hpp>`
install(TARGETS ${LIBRARY_NAME} ${APPLICATION_NAME}
  EXPORT ${PROJECT_NAME}Targets            # for downstream dependencies
  ARCHIVE DESTINATION lib COMPONENT libs   # static lib
  RUNTIME DESTINATION bin COMPONENT libs   # binaries
  LIBRARY DESTINATION lib COMPONENT libs   # shared lib
  FRAMEWORK DESTINATION bin COMPONENT libs # for mac
  PUBLIC_HEADER DESTINATION include/${PROJECT_NAME} COMPONENT devel   # headers for mac (note the different component -> different package)
  INCLUDES DESTINATION include             # headers
)

