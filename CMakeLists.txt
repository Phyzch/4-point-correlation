cmake_minimum_required(VERSION 3.10.2)
project(4_point_correlation_calculation)

include (CTest)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

include_directories("/home/phyzch/CLionProjects/4_point_correlation_calculation/MPI_Code")

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER  mpicc)

add_executable(4_point_correlation_calculation main.cpp Detector.cpp Full_system.cpp System.cpp
        auxiliary_function.cpp util.h quotient_state.h
        ./MPI_Code/Read_MPI.cpp MPI_Code/save_load_detector_MPI.cpp
        ./MPI_Code/construct_detector_MPI.cpp ./MPI_Code/Communication_function_MPI.cpp
        ./MPI_Code/Detector_Evolve_MPI.cpp
        ./MPI_Code/System_Read_MPI.cpp ./MPI_Code/compute_detector_energy_window_MPI.cpp
        vanvlk_construct_detector_offdiag.cpp MPI_Code/OTOC_for_xp.cpp MPI_Code/Another_form_of_OTOC.cpp Eigenvalue_spectrum.cpp MPI_Code/Verify_TOC_factorization.cpp MPI_Code/Boltzmann_factor_Chebyshev.cpp MPI_Code/ladder_operator_operation.cpp)

target_link_libraries(4_point_correlation_calculation stdc++fs mpi)