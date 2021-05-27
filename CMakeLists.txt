cmake_minimum_required(VERSION 3.10.2)
project(4_point_correlation_calculation)

include (CTest)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
# Find the BLAS stuff
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules)

set(BLA_DEFINITIONS)

if ("$ENV{MKLROOT}" STREQUAL "")
    message(FATAL_ERROR "MKLROOT is not set. Please source the Intel MKL mklvars.sh file.")
endif()

# user defined options for MKL
option(MKL_USE_parallel "Use MKL parallel" True)
option(MKL_USE_sdl "Single Dynamic Library or static/dynamic" False)
set(MKL_USE_interface "lp64" CACHE STRING "for Intel(R)64 compatible arch: ilp64/lp64 or for ia32 arch: cdecl/stdcall")

if (BLA_VENDOR MATCHES "_seq")
    set(MKL_USE_parallel OFF)
else()
    set(MKL_USE_parallel ON)
endif()

find_package(MKL REQUIRED)

if (MKL_FOUND)
    set(BLA_INCLUDE_DIR ${MKL_INCLUDE_DIR})
    set(BLAS_LIBRARIES ${MKL_LIBRARIES})
    set(BLA_DEFINITIONS USE_MKL)
endif()

message(STATUS "BLAS Libraries: ${BLAS_LIBRARIES}")
# include the subdirectory containing our libs
include_directories(${MKL_INCLUDE_DIRS})

include_directories("/home/phyzch/CLionProjects/4_point_correlation_calculation/MPI_Code")

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER  mpicc)

add_executable(4_point_correlation_calculation main.cpp Detector.cpp Full_system.cpp System.cpp
        auxiliary_function.cpp util.h quotient_state.h
        ./MPI_Code/Read_MPI.cpp MPI_Code/save_load_detector_MPI.cpp
        ./MPI_Code/construct_detector_MPI.cpp ./MPI_Code/Communication_function_MPI.cpp
        ./MPI_Code/Detector_Evolve_MPI.cpp
        ./MPI_Code/System_Read_MPI.cpp ./MPI_Code/compute_detector_energy_window_MPI.cpp
        vanvlk_construct_detector_offdiag.cpp MPI_Code/OTOC_for_xp.cpp MPI_Code/Another_form_of_OTOC.cpp Eigenvalue_spectrum.cpp USE_MKL_MODULE.cpp Compute_Eigenstate_OTOC.cpp New_MKL_file.cpp)

target_link_libraries(4_point_correlation_calculation stdc++fs mpi)
# Indicate which libraries to include during the link process.
target_link_libraries(4_point_correlation_calculation "-Wl,--start-group" ${MKL_LIBRARIES} "-Wl,--end-group -lpthread -lm -ldl")