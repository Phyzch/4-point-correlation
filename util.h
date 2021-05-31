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
#include "mkl.h"

//using namespace concurrency;
#define pi2 3.141592653589793*2

// time is in unit of ps. energy is in unit of cm^{-1} as in spectroscopy.
#define  cf2  0.0299792458* pi2

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
extern bool compute_state_space_and_coupling_using_symmetry_bool;

extern bool compute_eigenvector_use_MKL_module ;
extern bool compute_Eigenstate_OTOC_bool ;
extern double Emin, Emax; // range to solve eigenvalue and eigenvector using Intel MKL.  We will read it from input file.
extern double Emin2, Emax2;
extern double Emin_for_core, Emax_for_core;

extern int Symmetry_Table[4][4];
extern bool  use_multiple_core_to_solve_eigenstate_spectrum;

// define function here
float ran2(long& idum);
void estimate_memory_cost(ofstream & resource_output);  // output resource cost to file at give time step.
void convert_dv(const vector<vector<int>> & vec_2d, vector <int>  & vec_1d , vector <int> & displacement , vector <int> & element_size );
// used for cnostruct buffer for communication between process for matrix multiplication.
int construct_send_buffer_index(int * remoteVecCount, int * remoteVecPtr, int * remoteVecIndex, int * tosendVecCount_element, int * tosendVecPtr_element, int * & tosendVecIndex_ptr);

void allocate_diagonalization_energy_range_for_diff_proc( vector<double> & sorted_dmat_diagonal_part ,const  vector<double> & Emin_energy_range_list,const vector<double> & Emax_energy_range_list,
                                                          vector<double>  & Emin_list , vector<double> & Emax_list);

int MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector_submodule(  int * dirow_list,  int * dicol_list,  double * dmat_list , int dmatsize  ,int dmatnum,
                                                             vector<double> & dmat_diagonal_part ,
                                                             vector<double> & E , vector<vector<double>> & Matrix_X,
                                                             double Emin_of_choice, double Emax_of_choice );

int MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector_divide_by_part(int * dirow_list,  int * dicol_list,  double * dmat_list , int dmatsize  ,int dmatnum,
                                                                          vector<double> & dmat_diagonal_part ,
                                                                          vector<double> & Eigenvalue_list , vector<vector<double>> & Eigenvector_list,
                                                                          vector<double> & Emin_list, vector<double> & Emax_list );

int MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector_given_energy_and_num(  int * dirow_list,  int * dicol_list,  double * dmat_list , vector<double> &  dmat_diagonal_list ,  int dmatsize  ,int dmatnum,
                                                                                  ofstream & Eigenvector_output,
                                                                                  double * &E , double ** &Matrix_X,
                                                                                  double energy_of_choice, int eigenstate_number);

void construct_energy_window_for_eigenstate(int nmode, double * mfreq, double eigen_energy, double energy_window_for_eigenstate,
                                            vector<double> & Energy_range_min_list,
                                            vector<double> & Energy_range_max_list);

void convert_COO_to_CSR(const int * dirow_list, const int * dicol_list, const double * dmat_list , int dmatsize, int dmatnum,
                        vector<int> & dirow_CSR_form_fortran, vector<int> & dicol_CSR_form_fortran, vector<double> & sorted_dmat );


int compar(const void * a, const void * b);


struct phi_operator_phi_tuple{
    // record phi_operator_phi information with value larger than criteria
    int eigenstate_index;  // state m linked to state l
    double phi_operator_phi_value;  // <phi_m | a_{i} | phi_l>
    phi_operator_phi_tuple(int state_index , double phi_operator_phi_value1){
        eigenstate_index = state_index;
        phi_operator_phi_value = phi_operator_phi_value1;
    }
};

struct phi_operator_phi_tuple_complex{
    // record phi_operator_phi information with value larger than criteria
    int eigenstate_index;  // state m linked to state l
    complex<double> phi_operator_phi_value;  // <phi_m | a_{i} | phi_l>
    phi_operator_phi_tuple_complex(int state_index , complex<double> phi_operator_phi_value1){
        eigenstate_index = state_index;
        phi_operator_phi_value = phi_operator_phi_value1;
    }
};

#endif //QUANTUM_MEASUREMENT_UTIL_H

