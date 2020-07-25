set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake) # project specific cmake dir
set(CMAKE_CXX_STANDARD 11) # tODO move up to a general cmake config for all sub projects ?

#Directories where make is going to look for new projects
#add_subdirectory(wannier2sparse)
add_subdirectory(python)
