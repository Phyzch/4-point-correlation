//
// Created by phyzch on 6/18/20.
//
#pragma once
#ifndef QUANTUM_MEASUREMENT_UTIL_H
#define QUANTUM_MEASUREMENT_UTIL_H

#include<iostream>
#include<time.h>
#include<stdio.h>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <experimental/filesystem>
#include <random>
#include<iomanip>
#include <complex>
#include <assert.h>
#include <vector>
#include<ctime>
#include<algorithm>
#include<stdlib.h>
#include<mpi/mpi.h>
#include<sys/resource.h>
//using namespace concurrency;
#define pi2 3.141592653589793*2
using namespace std;
extern bool intra_detector_coupling;
extern double intra_detector_coupling_noise;
extern bool inter_detector_coupling;
extern double inter_detector_coupling_noise;
extern bool Continue_Simulation;
extern bool energy_window;
extern double energy_window_size;
extern double initial_energy;
extern double system_energy;
extern double noise_strength;
extern int Rmax;                // defined in compute_matrix_energy_window

extern bool detector_only;
extern bool Detector_Continue_Simulation;
extern bool Random_bright_state;
extern double detector_lower_bright_state_energy_window_shrink;

extern bool Turn_on_Vanvleck;
extern int ndegre;
extern int ndegrx2;

extern int distance_cutoff_for_4_piont_corre; // use to shrink simulation for 4-point correlation function.
extern double Energy_Range_4_point_corre_function_average;
extern int Distance_Range_4_point_corre_function_average;
extern bool turn_on_random_self_anharmonicity;
extern bool Sphere_cutoff_in_state_space;
extern bool read_Hamltonian_from_file;
extern bool save_state;
extern bool Evolve_dynamics;  // bool variable to decide if we run Detector_Evolve
extern bool compute_eigenvalue_spectrum; // bool variable to decide if use Lanczos algorithm to compute spectrum of syste
extern bool no_coupling;  // if this is Ture, we do not have off-diagonal coupling
extern bool compute_overlap_with_eigenstate; // If this is True, we will use MFD to compute overlap of initial state with eigenvalue
extern bool compute_state_space_and_coupling_suing_symmetry_bool;

// define function here
float ran2(long& idum);
void estimate_memory_cost(ofstream & resource_output);  // output resource cost to file at give time step.
void convert_dv(const vector<vector<int>> & vec_2d, vector <int>  & vec_1d , vector <int> & displacement , vector <int> & element_size );
// used for cnostruct buffer for communication between process for matrix multiplication.
int construct_send_buffer_index(int * remoteVecCount, int * remoteVecPtr, int * remoteVecIndex, int * tosendVecCount_element, int * tosendVecPtr_element, int * & tosendVecIndex_ptr);

int compar(const void * a, const void * b);

#endif //QUANTUM_MEASUREMENT_UTIL_H

