#include<iostream>
#include"util.h"//
// Created by phyzch on 6/22/20.
//

#ifndef QUANTUM_MEASUREMENT_QUOTIENT_STATE_H
#define QUANTUM_MEASUREMENT_QUOTIENT_STATE_H
using namespace std;

struct sys_quotient_state {
    vector<int> vmode1; // mode of detector 1
    vector<int> vmode2; // mode of detector 2
    // we can also store the detector index in it.
    int dindex1;  // index in detector 1
    int dindex2;  // index in detector 2
    vector <int> xindex;  // index in system + detector wave function. (x,y)
    vector <int> sysxindex;  // index in system reduced density matrix.
    sys_quotient_state(const vector<int> & vmode1_0, const vector<int> & vmode2_0, int dxindex1, int dxindex2){
        vmode1 = vmode1_0;
        vmode2 = vmode2_0;
        dindex1=dxindex1;
        dindex2=dxindex2;
    }
};

struct quotient_state {   // detector quotient_space_state.
    int sys_state;  // 0,1,2,3
    vector<int> vmode; // mode of detector

    vector <int> xindex;  // index in system + detector wave function. (x,y)
    vector <int> dxindex; // index in detector reduced density matrix basis (dx, dy).
    vector<vector<int>> q_index_list; // list of tuple (i,j,k,l,m): i: reduced density matrix index for rho_{i,j}, j: reduced density matrix index for rho_{i,j}
    // k index in xindex for i (use binary search at first time if we compute density matrix, then store result in k), l index in xindex for j, m: index in dmat
    vector<double> dmat_value_list;
    quotient_state(vector<int>  & vmode1, int sys_state1){
        sys_state= sys_state1;
        vmode=vmode1;
    }
};

int binary_insert_dxindex(vector <int> & list, int key);
int compare_quotient_state(int sys_state1, const vector <int>  & vmode1, int sys_state2, const vector <int> & vmode2);
int binary_search_dxindex(const vector <int> & list, int key, bool & exist);
void insert_quotient_state(vector <quotient_state> & list, int sys_state, vector<int> & vmode1,  int xindex, int dxindex);
void save_detector_quotient_state_data_for2( const vector <quotient_state> & d1list, const vector <quotient_state> & d2list, string path);
void load_detector_quotient_state_data_for2( vector<quotient_state> & d1list, vector <quotient_state> & d2list, string path);

void insert_sys_quotient_state(vector <sys_quotient_state> & list, const vector<int> vmode1_0, const vector<int> vmode2_0, int mod_dim, int xindex, int sysxindex, int dxindex1, int dxindex2);
vector<sys_quotient_state> merge_sort_list(vector<sys_quotient_state> & list);
void save_sys_quotient_state_data(const vector<sys_quotient_state> & list, string path);
void load_sys_quotient_state_data( vector <sys_quotient_state> &list, string path);


int find_position_for_insert_binary(const vector<vector<int>> & vmode, const vector<int> & ndetector, bool & exist);  // write in compute_matrix_energy_window.cpp
int find_location_binarysearch_sys_quotient_state(const vector<sys_quotient_state> & list, vector<int> & vmode1, const vector  <int> & vmode2, bool & exist);  // write in Compute_sys_reduced_density.cpp

#endif //QUANTUM_MEASUREMENT_QUOTIENT_STATE_H
