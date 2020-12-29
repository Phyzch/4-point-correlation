//
// Created by phyzch on 6/29/20.
//
#include "../system.h"
#include "../util.h"

using namespace std;
void full_system::construct_fullmatrix_with_energy_window_MPI() {
    int i;
    // compute diagonal part of full system Hamiltonian.
    compute_sstate_dstate_diagpart_dirow_dicol_MPI();
    // allocate memory for pointer sdnum, sdindex, sdmode
    sdnum = new int[s.tldim];
    total_sd_num = new int [s.tldim];
    sdnum_each_process = new int * [s.tldim];
    sdnum_displacement_each_process = new int * [s.tldim];
    for(i=0;i<s.tldim;i++){
        sdnum_each_process [i] = new int [num_proc];
        sdnum_displacement_each_process [i] = new int [num_proc];
    }

    sdindex = new  vector<int> [s.tldim];
    sdmode = new vector<int> [s.tldim];
    construct_quotient_state_all_MPI();
    compute_offdiagonal_part_MPI();
}

void full_system:: compute_sstate_dstate_diagpart_dirow_dicol_MPI(){
   // compute diagonal part and sstate, dstate for full_system.
    int i,j,k,l;
    int matnum=0;
    double energy;
    // each process responsible to couple detector 1's state with detector 2's state to form state in x,y.
    int vsize_d1 = d.total_dmat_size[0]/num_proc;
    int begin_index_d1 = my_id * vsize_d1;
    int end_index_d1;
    if(my_id!=num_proc-1){
        end_index_d1= (my_id +1) * vsize_d1;
    }
    else{
        end_index_d1 = d.total_dmat_size[0];
    }
    // sstate ,dstate is index for matrix element in full matrix to record the corresponding index in system and detector.
    dstate = new vector <int> [s.tldim];

    if(s.tlnum==1){
        for(j=begin_index_d1;j<end_index_d1;j++){
            for(i=0;i<s.tlmatsize;i++){
                energy= s.tlmat[i] + dmat0[j];
                if(abs(energy - initial_energy) < energy_window_size){
                    sstate.push_back(i);
                    dstate[0].push_back(j);
                    mat.push_back(energy);
                    irow.push_back(matnum);
                    icol.push_back(matnum);
                    matnum = matnum +1 ;
                }
            }
        }
    }
    else if(s.tlnum==2){
        for(j= begin_index_d1 ; j< end_index_d1;j++){
                for (k=0;k<d.total_dmat_size[1];k++) {
//                    for(i=0;i<s.tlmatsize;i++){
                        i=0;
                        energy = s.tlmat[i] + dmat0[j] + dmat1[k];
//                        if (abs(energy - initial_energy) < energy_window_size) {
                            sstate.push_back(i);
                            dstate[0].push_back(j); // dstate record detector global index
                            dstate[1].push_back(k);
                            mat.push_back(energy);
                            irow.push_back(
                                    matnum); //matnum is local, we have to re-compute it after all process compute its own irow, icol.
                            icol.push_back(matnum);
                            matnum = matnum + 1;
//                        }
//                    }
            }
        }
    }

    // you have to re-assign the irow, icol in each process:
    matsize= mat.size();
    matsize_each_process= new int [num_proc];
    matsize_offset_each_process = new int [num_proc];
    MPI_Allgather(&matsize,1,MPI_INT,&matsize_each_process[0],1,MPI_INT,MPI_COMM_WORLD);

    matsize_offset_each_process[0] = 0;
    for (i = 1; i < num_proc; i++) {
        matsize_offset_each_process[i] = matsize_offset_each_process[i - 1] + matsize_each_process[i-1];
    }
    total_matsize=0;
    for(i=0;i<num_proc;i++){
        total_matsize= total_matsize + matsize_each_process[i];
    }
    int offset;
    offset= matsize_offset_each_process[my_id];
    // each process re-assign irow, icol according to offset. Now irow, icol is index in global.
    for(i=0;i<matsize;i++){
        irow[i]= irow[i] + offset;
        icol[i]= icol[i] + offset;
    }
}


// initialize the wavefunction of our system.
// This code uses function  etot(). MPI_version of etot() is not finished yet.
void full_system::Initial_state_MPI() {
    // Construct initial wavefunction for full matrix
    // we should use information stored in vmode0, vmode1.
    double x1, y1, x2, y2, x3, y3;
    int m, i;
    norm = 0;
    x.reserve(matsize);
    y.reserve(matsize);

    // collect xd yd from all other process.
    d.gather_xd_yd();

    double value;
    for (i = 0; i < matsize; i++) {
        x1 = s.xtl[sstate[i]];
        y1 = s.ytl[sstate[i]];
        x2 = d.xd_all[0][dstate[0][i]];
        y2 = d.yd_all[0][dstate[0][i]];
        x3 = d.xd_all[1][dstate[1][i]];
        y3 = d.yd_all[1][dstate[1][i]];
        value = x1*(x2*x3 - y2*y3) - y1*(x2*y3 + x3*y2);
        x.push_back(value);
        value = y1*(x2*x3 - y2*y3) + x1*(x2*y3 + x3*y2);
        y.push_back(value);
        norm = norm + pow(x[i], 2) + pow(y[i], 2);
        if(isnan(norm)){
            cout<< "Norm is not a number"<<endl;
        }
    }

    MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    total_norm = 1 / sqrt(total_norm);
    for (i = 0; i < matsize; i++) {
        x[i] = x[i] * total_norm;
        y[i] = y[i] * total_norm;
    }
    // convert all matrix elements for ps-wavenumber units
    for (i = 0; i < matnum; i++) {
        mat[i] = cf * mat[i];
    }

}
void full_system::shift_mat(){
    int m;
    int i;
    // MPI_version of computing energy is not finished.
    double *hx, *hy;
    hx = new double[matsize];  // we have to rewrite this part. hx should be allocated space outside the function.
    hy = new double[matsize];
    etot_MPI(hx,hy);  // calculate total energy of system.  MPI_version is not finished yet.

    // Shift matrix for maximum SUR algorithm efficiency unless matflag ==2
    if (d.matflag != 2) {
        for (i = 0; i < matsize; i++) {
            mat[i] = mat[i] - total_energy;
        }
    }

//    for(i=0;i<matsize;i++) {
//        cout << mat[i] << "  ";
//        if(my_id == 0){
//            log << mat[i]<<" ";
//        }
//    }
//    cout <<endl;

    if(my_id==0) {
        if(! Detector_Continue_Simulation) {
            output << "Total energy before SUR algorithm shifting  " << total_energy / cf << endl;
            output << "Intra-detector coupling strength V set as: " << d.V_intra << endl;
            output << "Intra-detector coupling scaling parameter a set as  " << d.a_intra << endl;
            output << "We add noise in mode frequency, noise strength:  " << noise_strength << endl;
            output << "Initial Detector state is : " << endl;
            for (m = 0; m < s.tlnum; m++) {
                for (i = 0; i < d.nmodes[m]; i++) {
                    output << d.initial_detector_state[m][i] << " ";
                }
                output << endl;
            }

            if (energy_window) {
                output << "We apply energy window here:  " << energy_window_size << endl;
            }
        }
    }
    delete [] hx;
    delete [] hy;
}

