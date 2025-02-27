cmake_minimum_required(VERSION 3.10)
project(fem_solver C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

# Force GCC (MinGW)
if (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    message(FATAL_ERROR "MSVC detected! Run CMake with GCC: cmake -G Ninja -DCMAKE_C_COMPILER=gcc")
endif()

# Enable OpenMP (for GCC)
find_package(OpenMP REQUIRED)

# Set the source and include directories
set(FEM_SRC_DIR ${CMAKE_SOURCE_DIR}/src)
set(FEM_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)
set(MATRIX_SRC_DIR ${CMAKE_SOURCE_DIR}/../matrix/src)
set(MATRIX_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../matrix/include)

# Include directories (FEM and Matrix)
include_directories(${FEM_INCLUDE_DIR} ${MATRIX_INCLUDE_DIR})

# Define source files
# Add FEM source files
set(FEM_SOURCES
    ${FEM_SRC_DIR}/elliptic_solver.c
    ${FEM_SRC_DIR}/mesh.c
)

set(MATRIX_SOURCES
    ${MATRIX_SRC_DIR}/matrix_ops.c
    ${MATRIX_SRC_DIR}/qr.c
    ${MATRIX_SRC_DIR}/cg.c
    ${MATRIX_SRC_DIR}/lsqr.c
    ${MATRIX_SRC_DIR}/gmres.c
)

# Create FEM library
add_library(fem_solver ${FEM_SOURCES})

# Ensure OpenMP is used for GCC
if(OpenMP_C_FOUND)
    target_link_libraries(fem_solver PUBLIC OpenMP::OpenMP_C)
endif()

# Add the matrix module (explicit binary directory)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../matrix ${CMAKE_BINARY_DIR}/matrix)

# Include directories for external linking
target_include_directories(fem_solver PUBLIC ${MATRIX_INCLUDE_DIR})
target_link_directories(fem_solver PUBLIC ${CMAKE_BINARY_DIR}/../matrix)

# Link the matrix module (instead of individual files)
target_link_libraries(fem_solver PUBLIC matrix_module)

# Create test executables
add_executable(test_elliptic_solver tests/test_elliptic_solver.c)
add_executable(test_mesh tests/test_mesh.c)

# Link tests with FEM solver and matrix module
target_link_libraries(test_elliptic_solver PUBLIC fem_solver matrix_module m)
target_link_libraries(test_mesh PUBLIC fem_solver matrix_module m)