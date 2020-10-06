//
// Created by phyzch on 9/10/20.
// Save Hamiltonian and Load Hamiltonian, wave function.
//
#include "../system.h"
#include"../util.h"
using namespace std;

void full_system::save_wave_function_MPI(){
    double * x_all ;
    double * y_all ;
    int i,j;
    if(my_id == 0){
        save.open(path+"save_wavefunction.txt");
        if( ! save.is_open()){
            output << "Can not open the save_wavefunction.txt file"<<endl;
            MPI_Abort(MPI_COMM_WORLD,-10);
        }
    }
    if (my_id == 0) {
        x_all = new double[total_matsize];
        y_all = new double[total_matsize];
    }
    gather_x_y(x_all, y_all);
    if(my_id ==0) {
        for (i = 0; i < s.tlnum; i++) {
            for (j = 0; j < d.nmodes[i]; j++) {
                save << d.mfreq[i][j] << " ";
            }
        }
        save << endl;
        save << t << endl;
        save << "Real Part of Wavefunction:" << endl;
        for (i = 0; i < total_matsize; i++) {
            save << x_all[i] << " ";
        }
        save << endl;
        save << "Imaginary Part of Wavefunction:" << endl;
        for (i = 0; i < total_matsize; i++) {
            save << y_all[i] << " ";
        }
        save << endl;
        delete[] x_all;
        delete[] y_all;
    }

}

void full_system::load_wave_function_MPI() {
    string ss;
    int i,j ;
    double * x_all, * y_all;
    if(my_id == 0){
        load.open(path+"save_wavefunction.txt");
        if(! load.is_open()){
            output << "Can not open the save_wavefunction.txt file"<<endl;
            MPI_Abort(MPI_COMM_WORLD,-10);
        }
        x_all = new double[total_matsize];
        y_all = new double[total_matsize];
        for (i = 0; i < s.tlnum; i++) {
            for (j = 0; j < d.nmodes[i]; j++) {
                load >> d.mfreq[i][j];
            }
        }
        std::getline(load, ss);
        load >> t;
        std::getline(load, ss);
        std::getline(load, ss);
        for (i = 0; i < total_matsize; i++) {
            load >> x_all[i];
        }
        std::getline(load, ss);
        std::getline(load, ss);
        for (i = 0; i < total_matsize; i++) {
            load >> y_all[i];
        }
        load.close();
    }

    // Broadcast time, frequency to all other process
    MPI_Bcast(&t,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(i=0;i<s.tlnum;i++){
        MPI_Bcast(&d.mfreq[i][0],d.nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    scatter_x_y(x_all,y_all);

    for (i = 0; i < matnum; i++) {
        mat[i] = cf * mat[i];
    }
    if(my_id == 0){
        delete [] x_all;
        delete [] y_all;
    }
}

void full_system::save_Hamiltonian_MPI(){
    int m,i,j;
    double * mat_all;
    int * irow_all, * icol_all;
    int * sstate_all;
    int ** dstate_all;
    int ** sdindex_all;
    int ** sdmode_all;
    if(my_id==0) {
        save.open(path + "save_Hamiltonian.txt");
        if(! save.is_open()){
            log<<"Can not open file that store Hamiltonian for photon + detector"<<endl;
            MPI_Abort(MPI_COMM_WORLD,-9);
        }
    }
    dstate_all = new int * [2];
    sdindex_all = new int * [2];
    sdmode_all = new int *[2];
    if(my_id == 0){
        mat_all = new double [total_matnum];
        irow_all = new int [total_matnum];
        icol_all = new int [total_matnum];
        sstate_all = new int [total_matsize];
        for(i=0;i<s.tlnum;i++){
            dstate_all[i] = new int [total_matsize];
            sdindex_all[i] = new int [total_sd_num[i]];
            sdmode_all[i] = new int [total_sd_num[i]];
        }
    }
    gather_mat_irow_icol_sstate_dstate_sdmode_sdindex(mat_all,irow_all,icol_all,sstate_all,dstate_all,sdmode_all,sdindex_all);
    if(my_id == 0){
        save<< total_matsize << " "<<total_offnum <<" " << total_matnum <<"  ";
        for(m=0;m<s.tldim;m++){
            save << total_sd_num[m] <<"  ";
        }
        // record matrix size and number for each process
        for(i=0;i<num_proc;i++){
            save << matsize_each_process[i] << "  "<< matsize_offset_each_process[i] << "  "<< mat_offnum_each_process[i]<<"  "<<matnum_each_process[i] <<"  " << matnum_offset_each_process[i]<<"  ";
        }
        for(m=0;m<s.tldim;m++){
            for(i=0;i<num_proc;i++){
                save << sdnum_each_process[m][i] <<"  ";
                save << sdnum_displacement_each_process[m][i] <<" ";
            }
        }
        save << endl;

        for(i=0;i<total_matnum;i++){
            save << mat_all[i] <<"  "<< irow_all[i] <<"  "<< icol_all[i] << "  ";
        }
        save <<endl;

        for(i=0;i<total_matsize;i++){
            save << sstate_all[i] << " ";
            for(m=0;m<s.tlnum;m++){
                save << dstate_all[m][i] <<" ";
            }
        }
        save<<endl;

        for(m=0;m<s.tlnum;m++){
            for(i=0;i<total_sd_num[m];i++){
                save << sdindex_all[m][i] <<" "<< sdmode_all[m][i] <<" ";
            }
        }

        save.close();
    }

    if(my_id == 0){
        delete [] mat_all;
        delete [] irow_all;
        delete [] icol_all;
        for(m=0;m<s.tlnum;m++){
            delete [] dstate_all[m];
            delete [] sdindex_all[m];
            delete [] sdmode_all[m];
        }
        delete [] sstate_all;
    }
    delete [] dstate_all;
    delete [] sdindex_all;
    delete [] sdmode_all;
}

void full_system::load_Hamiltonian_MPI() {
    // in load function you have to allocate space for matsize_each_process, matnum_each_process etc.
    int i,m;
    string ss;

    double * mat_all;
    int * irow_all, * icol_all;
    int * sstate_all;
    int ** dstate_all;
    int ** sdindex_all;
    int ** sdmode_all;

    // allocate memory for pointer sdnum, sdindex, sdmode
    sdnum = new int[s.tldim];
    total_sd_num = new int [s.tldim];
    sdnum_each_process = new int * [s.tldim];
    sdnum_displacement_each_process = new int * [s.tldim];
    for(i=0;i<s.tldim;i++){
        sdnum_each_process [i] = new int [num_proc];
        sdnum_displacement_each_process [i] = new int [num_proc];
    }
    matsize_each_process= new int [num_proc];
    matsize_offset_each_process = new int [num_proc];
    matnum_each_process = new int [num_proc];
    matnum_offset_each_process = new int [num_proc];
    mat_offnum_each_process = new int [num_proc];

    dstate = new vector<int> [s.tlnum];
    sdmode = new vector<int> [s.tlnum];
    sdindex = new vector<int> [s.tlnum];
    if(my_id == 0){
        load.open(path + "save_Hamiltonian.txt");
        if(! load.is_open()){
            log<<"Can not open file that store Hamiltonian for photon + detector"<<endl;
            MPI_Abort(MPI_COMM_WORLD,-9);
        }
    }

    if(my_id == 0){
        load >> total_matsize >> total_offnum>> total_matnum;
        for(m=0;m<s.tldim;m++){
            load >> total_sd_num[m];
        }
        for(i=0;i<num_proc;i++){
           load >> matsize_each_process[i] >> matsize_offset_each_process[i] >> mat_offnum_each_process[i] >> matnum_each_process[i] >> matnum_offset_each_process[i];
        }
        for(m=0;m<s.tlnum;m++){
            for(i=0;i<num_proc;i++){
                load >> sdnum_each_process[m][i] >> sdnum_displacement_each_process[m][i] ;
            }
        }

        std::getline(load, ss);
    }
    // now we begin to allocate space for  mat_all, irow_all, icol_all, sstate_all, dstate_all, sdindex_all, sdmode_all
    dstate_all = new int * [2];
    sdindex_all = new int * [2];
    sdmode_all = new int *[2];
    if(my_id == 0){
        mat_all = new double [total_matnum];
        irow_all = new int [total_matnum];
        icol_all = new int [total_matnum];
        sstate_all = new int [total_matsize];
        for(i=0;i<s.tlnum;i++){
            dstate_all[i] = new int [total_matsize];
            sdindex_all[i] = new int [total_sd_num[i]];
            sdmode_all[i] = new int [total_sd_num[i]];
        }
    }

    if(my_id == 0){
        for(i=0;i<total_matnum;i++){
            load >> mat_all[i] >> irow_all[i] >> icol_all [i];
        }
        std::getline(load, ss);

        for(i=0;i<total_matsize;i++){
            load >> sstate_all[i];
            for(m=0;m<s.tlnum;m++){
                load >> dstate_all[m][i];
            }
        }
        std::getline(load, ss);

        for(m=0;m<s.tlnum;m++){
            for(i=0;i<total_sd_num[m];i++){
                load >> sdindex_all[m][i] >> sdmode_all[m][i];
            }
        }

        load.close();
    }

    MPI_Bcast(&total_matsize,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&total_matnum,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&total_offnum,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&matsize_each_process[0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&matsize_offset_each_process[0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&mat_offnum_each_process[0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&matnum_each_process[0],num_proc,MPI_INT, 0 , MPI_COMM_WORLD);
    MPI_Bcast(&matnum_offset_each_process[0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
    for(m=0;m<s.tlnum;m++){
        MPI_Bcast(&total_sd_num[m],1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&sdnum_each_process[m][0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&sdnum_displacement_each_process[m][0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
    }
    matsize = matsize_each_process[my_id];
    matnum = matnum_each_process[my_id];
    offnum = mat_offnum_each_process[my_id];
    for(m=0;m<s.tlnum;m++){
        sdnum[m] = sdnum_each_process[m][my_id];
    }
    // scatter matrix and irow, icol , sstate etc to all other process from process 0
    scatter_mat_irow_icol_sstate_dstate_sdmode_sdindex(mat_all,irow_all,icol_all,sstate_all,dstate_all,sdmode_all,sdindex_all);

    if(my_id == 0){
        delete [] mat_all;
        delete [] irow_all;
        delete [] icol_all;
        for(m=0;m<s.tlnum;m++){
            delete [] dstate_all[m];
            delete [] sdindex_all[m];
            delete [] sdmode_all[m];
        }
        delete [] sstate_all;
    }
    delete [] dstate_all;
    delete [] sdindex_all;
    delete [] sdmode_all;
}