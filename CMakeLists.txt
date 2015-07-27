# CMakeLists.txt for CAFFE
# Author: Adam O'Brien

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(CAFFE)

# Version number

set(CAFFE_VERSION_MAJOR 1)
set(CAFFE_VERSION_MINOR 0)

# Compiler
set(CMAKE_C_COMPILER "mpicc")
set(CMAKE_CXX_COMPILER "mpic++")

# Compiler flags
set(CMAKE_CXX_FLAGS "-g -Wall -Wno-reorder -Wno-sign-compare -Wno-switch -std=c++11 -O3")
set(CMAKE_C_FLAGS "-Wall")

# Output paths
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)


# Find mpi
find_package(MPI REQUIRED)

if(MPI_FOUND)
    include_directories(${MPI_INCLUDE_DIRS})
    message("-- MPI_CXX includes: " ${MPI_CXX_INCLUDE_PATH})
    message("-- MPI_CXX libraries: " ${MPI_CXX_LIBRARIES})
endif()

# Find boost package
find_package(Boost 1.53.0 REQUIRED COMPONENTS program_options system date_time filesystem)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    message("-- Boost version: " ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION})
    message("-- Boost includes: " ${Boost_INCLUDE_DIRS})
    message("-- Boost libraries: " ${Boost_LIBRARY_DIRS})
endif()

# Find blas package
find_package(BLAS)

if(BLAS_FOUND)
	message("-- BLAS libraries: " ${BLAS_LIBRARIES})
endif()

# Find lapack package
find_package(LAPACK)

if(LAPACK_FOUND)
	message("-- LAPACK libraries: " ${LAPACK_LIBRARIES})
endif()

# Includes
set(INCLUDE_DIR src/DataStructures
		src/Math
		src/RunControl
                src/Domains/HexaFvmMesh
                src/Domains/Field
		src/Solvers
		src/Schemes
                src/Parallel)

include_directories(${INCLUDE_DIR})

# Libraries
add_subdirectory(src)

# Modules
add_subdirectory(Modules)

# Testing
add_subdirectory(test)
