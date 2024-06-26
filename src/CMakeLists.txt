# @author  Jose H. Garcia


#CREATES KPM_LIB
set(KPM_LIB_SOURCES
 "chebyshev_solver.cpp"
 "quantum_states.cpp"
 "chebyshev_vectors.cpp"
 "chebyshev_moments.cpp"
 "chebyshev_moments2D.cpp"
 "chebyshev_moments1D.cpp"
 "chebyshev_momentsTD.cpp"
 "sparse_matrix.cpp"
 "special_functions.cpp"
 "Kubo_solver_FFT.cpp"
 "Kubo_solver_FFT_postProcess.cpp"
 "Greenwood_FFTs.cpp"
 "Bastin_FFTs.cpp"
)

#IF INTEL BACKEND IS USED
if( ${INTEL_MKL} STREQUAL "TRUE" )
	set(KPM_LIB_SOURCES "${KPM_LIB_SOURCES}" "mkl_sparse_matrix.cpp" "mkl_linear_algebra.cpp" )
	message( ${KPM_LIB_SOURCES} )
endif()

#CREATION OF THE LIBRARY
set( KPM_LIB kpm_lib)
add_library(${KPM_LIB} STATIC ${KPM_LIB_SOURCES})
target_link_libraries(${KPM_LIB} ${BLAS_LIB})


#Create the  different executables
set(MAIN inline_compute-kpm-nonEqOp)
add_executable(${MAIN} inline_compute-kpm-nonEqOp.cpp )
target_link_libraries(${MAIN} ${KPM_LIB} )

set(MAIN inline_compute-kpm-CorrOp )
add_executable(${MAIN} inline_compute-kpm-CorrOp.cpp )
target_link_libraries(${MAIN} ${KPM_LIB} )

set(MAIN inline_compute-kpm-spectralOp )
add_executable(${MAIN} inline_compute-kpm-spectralOp.cpp )
target_link_libraries(${MAIN} ${KPM_LIB} )

set(MAIN inline_spectralFunctionFromChebmom)
add_executable(${MAIN} spectralFunctionFromChebmom.cpp )
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_kuboGreenwoodFromChebmom)
add_executable(${MAIN} kuboGreenwoodFromChebmom.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )

set(MAIN inline_kuboBastinIFromChebmom)
add_executable(${MAIN} kuboBastinIFromChebmom.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )

set(MAIN inline_kuboBastinIIFromChebmom)
add_executable(${MAIN} kuboBastinIIFromChebmom.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )

set(MAIN inline_kuboBastinFromChebmom)
add_executable(${MAIN} kuboBastinFromChebmom.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_kuboBastinSeaFromChebmom)
add_executable(${MAIN} kuboBastinSeaFromChebmom.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_kuboBastinSurfFromChebmom)
add_executable(${MAIN} kuboBastinSurfFromChebmom.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_timeCorrelationsFromChebmom )
add_executable(${MAIN} timeCorrelationsFromChebmom.cpp  )
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_compute-Kubo-kpm-FFT )
add_executable(${MAIN} inline_compute-Kubo-kpm-FFT.cpp  )
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_spectralFunctionFromChebmom_FFTgrid)
add_executable(${MAIN} spectralFunctionFromChebmom_FFTgrid.cpp )
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_kuboGreenwoodFromChebmom_FFTgrid)
add_executable(${MAIN} kuboGreenwoodFromChebmom_FFTgrid.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )


set(MAIN inline_kuboBastinFromChebmom_FFTgrid)
add_executable(${MAIN} kuboBastinFromChebmom_FFTgrid.cpp)
target_link_libraries(${MAIN} ${KPM_LIB} )

set(MAIN inline_compute-kpm-spectralOp-nonOrth )
add_executable(${MAIN} inline_compute-kpm-spectralOp-nonOrth.cpp )
target_link_libraries(${MAIN} ${KPM_LIB} )

