//
// Created by phyzch on 4/14/20.
//Used for putting non essential function in same file.
//
#include"system.h"
#include"util.h"
using namespace std;

void detector:: shift_detector_Hamiltonian(ofstream & log){
    int i;
    int irow_index, icol_index;
    // -----  compute state's energy and shift it before doing simulation -------------
    vector<complex<double>>  H_phi;
    H_phi.resize(dmatsize[0]);
    for(i=0;i<dmatsize[0];i++){
        H_phi[i] = 0;
    }

    double de;
    double de_all;


    // use ground electronic state specified in input.txt to compute energy and shift Hamiltonian according to it.
    for(i=0;i<dmatnum[0];i++){
        irow_index = local_dirow[0][i];
        icol_index = local_dicol[0][i]; // compute to point to colindex in
        H_phi[irow_index] = H_phi[irow_index] + dmat[0][i] * complex(xd[ initial_state_in_sampling_state_index_list[0] ][icol_index],
                                                                       yd[ initial_state_in_sampling_state_index_list[0] ][icol_index]);
    }
    de=0;
    for(i=0;i<dmatsize[0];i++){
        de= de+ real(H_phi[i] * complex(xd[initial_state_in_sampling_state_index_list[0] ][i],
                                        -yd[initial_state_in_sampling_state_index_list[0] ][i]));
    }
    MPI_Allreduce(&de,&de_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    // shift Hamiltonian by detector energy de_all:
    if(my_id == 0){
        cout <<"Shift Hamiltonian by energy of system  "<< de_all << endl;
        log << "Shift Hamiltonian by energy of system  " << de_all <<endl;
    }
    for(i=0;i<dmatsize[0];i++){
        dmat[0][i] = dmat[0][i] - de_all;
    }
}

// check some dimension parameter
void full_system::dimension_check() {
    double errflag = 0;
    if (!energy_window) {
        if (d.dmatnum[0] > d.dmatdim || d.dmatnum[1] > d.dmatdim) errflag = errflag + 1;
        if (d.dmatsize[0] > d.detdim || d.dmatsize[1] > d.detdim) errflag = errflag + 2;
    }
    if (errflag != 0) {
        log << " Dimension Problem: Error flag=" << errflag << endl;
        cout<<"Dimension error, Error flag="<<errflag;
        exit(-1);
    }

    if (s.tlnum == 1) {
        if (! Detector_Continue_Simulation) {
            output << "Global Matrix: 2*" << d.dmatsize[0] << " = " << matsize << endl;
        }
    }
    else if (s.tlnum == 2) {
        if (! Detector_Continue_Simulation ) {
            if (! energy_window) {
                output << "Global Matrix : 4*" << d.dmatsize[0] << " * " << d.dmatsize[1] << " = " << matsize << endl;
            }
            else {
                output << "Global Matrix: " << total_matsize << endl;
            }
        }
    }
    if (!Detector_Continue_Simulation) {
        output << "off-diagonal matrix number  " << total_offnum << endl;
        output << "Whole matrix element number  " << total_matnum << endl;
        log << "off-diagonal matrix number  " << total_offnum << endl;
        log << "Whole matrix element number  " << total_matnum << endl;
    }
};


