set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

if(TARGET lsquant)
    set_target_properties(lsquant PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/legacy_bin")
endif()

add_library(lsquant_core STATIC
    core/algebra/sparse_matrix.cpp
    core/algebra/composite_operators.cpp
    core/algebra/linalg.cpp
    core/algebra/state_generator.cpp
    core/algebra/spectral_bounds.cpp
    core/chebyshev/rescaled_operator.cpp
    core/chebyshev/chebyshev_recursion.cpp
    core/chebyshev/chebyshev_vectors.cpp
    core/chebyshev/chebyshev_moments.cpp
    core/chebyshev/moment_engine.cpp
    io/csr_text.cpp
    io/chebmom_file.cpp
)

target_include_directories(lsquant_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_features(lsquant_core PUBLIC cxx_std_17)
target_compile_options(lsquant_core PUBLIC -Wall -Wextra)
target_link_libraries(lsquant_core PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

add_executable(stage1_moments apps/stage1_moments.cpp)
target_compile_options(stage1_moments PRIVATE -Wall -Wextra)
target_link_libraries(stage1_moments PRIVATE lsquant_core)

add_executable(stage_1_check "${PROJECT_SOURCE_DIR}/tests/dev/stage_1_check.cpp")
target_compile_options(stage_1_check PRIVATE -Wall -Wextra)
target_link_libraries(stage_1_check PRIVATE lsquant_core)
