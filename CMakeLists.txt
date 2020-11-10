cmake_minimum_required(VERSION 3.10.2)
project(4_point_correlation_calculation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

include_directories("/home/phyzch/CLionProjects/4_point_correlation_calculation/MPI_Code")

SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER  mpicc)

add_executable(4_point_correlation_calculation main.cpp Detector.cpp Full_system.cpp System.cpp  resource_estimate.cpp
        Compute_Detector_reduced_density_matrix.cpp Compute_sys_reduced_density_matrix.cpp auxiliary_function.cpp util.h quotient_state.h
        ./MPI_Code/Read_MPI.cpp MPI_Code/save_load_detector_MPI.cpp
        ./MPI_Code/construct_detector_MPI.cpp ./MPI_Code/Communication_function_MPI.cpp
        ./MPI_Code/Detector_Evolve_MPI.cpp ./MPI_Code/construct_full_system_MPI.cpp ./MPI_Code/construct_quotient_state_MPI.cpp
        ./MPI_Code/Compute_Off_diagonal_full_system_MPI.cpp  ./MPI_Code/Compute_Output_MPI.cpp ./MPI_Code/full_system_evolve_MPI.cpp
        ./MPI_Code/System_Read_MPI.cpp ./MPI_Code/compute_detector_energy_window_MPI.cpp ./MPI_Code/save_load_fullsystem_MPI.cpp
        vanvlk_construct_detector_offdiag.cpp Lanczos_and_MFT.cpp)

target_link_libraries(4_point_correlation_calculation stdc++fs mpi)