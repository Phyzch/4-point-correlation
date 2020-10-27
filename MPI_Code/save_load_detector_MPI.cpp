//
// Created by phyzch on 6/23/20.
// This file contain function that save and load data using our previously program suitable for one process.
// we will wrap up these functions, load it from master process and send/ broadcast them to other process for computation.
//
#include"../system.h"
#include"../util.h"
using namespace std;

void detector:: save_detector_Hamiltonian_MPI(string path, ofstream & log){
    // we should record process number, dmat, dirow, dicol in each process (due to off-diagonal matrix element)
    int m,i,j;

    // now we have all dv, dirow, dicol in master process, we do I/O there.
    // first output  total_dmatsize, total_doffnum, total_dmatnum. Then for each process: dmatsize, doffnum, dmatnum
    // Then dmat_all, dirow_all, dicol_all data there.
    // finally save dv_all
    if(my_id==0){
        ofstream save;
        save.open(path+"save_detector_Hamiltonian.txt");
        if(save.is_open()){
            save << num_proc <<endl;
            for(m=0;m<stlnum;m++){
                save<< total_dmat_size[m] <<" "<< total_dmat_off_num[m] <<" "<<total_dmat_num[m]<<" ";
                for(i=0;i<num_proc;i++){
                    save<<dmatsize_each_process[m][i] << " "<< doffnum_each_process[m][i]<<" "<< dmatnum_each_process[m][i]<<" ";
                }
                save<<endl;
            }
            for(m=0;m<stlnum;m++){
                for(i=0;i<total_dmat_num[m];i++){
                    save<< total_dmat[m][i] << " "<< total_dirow[m][i]<<" "<<total_dicol[m][i]<<" ";
                }
                save<<endl;
            }
            for(m=0;m<stlnum;m++){
                for(i=0;i<total_dmat_size[m];i++){
                    for(j=0;j<nmodes[m];j++){
                        save<< dv_all[m][i][j]<<" ";
                    }
                }
                save<<endl;
            }
            save.close();
        }
        else{
            log<<" Error. We can not open save file for detector Hamiltonian."<<endl;
            MPI_Abort(MPI_COMM_WORLD,-11);
        }
    }

}


void detector::load_detector_Hamiltonian_MPI(string path, ofstream & log) {
    // load detector Hamiltonian from save_detector_Hamiltonian.txt
    int m, i, j;
    int value;
    int required_num_proc;
    total_dmat = new double * [stlnum];
    total_dirow = new int * [stlnum];
    total_dicol = new int * [stlnum];
    dv_all = new vector<vector<int>>[2];
    vector<int> temporary_dv;
    double ** temp_dmat = new double * [stlnum];
    int ** temp_dirow = new int * [stlnum];
    int ** temp_dicol = new int * [stlnum];
    if (my_id == 0) {
        ifstream load;
        string ss;
        load.open(path + "save_detector_Hamiltonian.txt");
        if (load.is_open()) {
            // check number of process is same as before.
            load >> required_num_proc;
            std::getline(load,ss);
            if(required_num_proc != num_proc){
                cout<<"Wrong number of process input. Previous simulation is run with num proc = "<< required_num_proc
                <<" Now it's   "<<num_proc<<endl;
                log <<"Wrong number of process input. Previous simulation is run with num proc = "<< required_num_proc
                    <<" Now it's   "<<num_proc<<endl;
                MPI_Abort(MPI_COMM_WORLD,-8);  // -8 for wrong number of process as input.
            }

            // load dmatsize, dmatnum, doff_num for each process.
            for (m = 0; m < stlnum; m++) {
                load >> total_dmat_size[m] >> total_dmat_off_num[m] >> total_dmat_num[m];
                for (i = 0; i < num_proc; i++) {
                    load >> dmatsize_each_process[m][i] >> doffnum_each_process[m][i] >> dmatnum_each_process[m][i];
                }
                std::getline(load, ss);
            }
            //use MPI_Gatherv function to gather dmat, dirow, dicol to process 0
            for (m = 0; m < stlnum; m++) {
                total_dmat[m] = new double [ total_dmat_num[m] ];
                total_dirow[m] = new int [ total_dmat_num[m] ];
                total_dicol[m] = new int [ total_dmat_num[m] ];
            }
            for (m = 0; m < stlnum; m++) {
                for (i = 0; i < total_dmat_num[m]; i++) {
                    load >> total_dmat[m][i] >> total_dirow[m][i] >> total_dicol[m][i];
                }
                std::getline(load, ss);
            }
            for (m = 0; m < stlnum; m++) {
                dv_all[m].reserve(total_dmat_size[m]);
                for (i = 0; i < total_dmat_size[m]; i++) {
                    dv_all[m][i].reserve(nmodes[m]);
                }
            }
            for (m = 0; m < stlnum; m++) {
                for (i = 0; i < total_dmat_size[m]; i++) {
                    temporary_dv.clear();
                    for (j = 0; j < nmodes[m]; j++) {
                        load >> value;
                        temporary_dv.push_back(value);
                    }
                    dv_all[m].push_back(temporary_dv);
                }
                std::getline(load, ss);
            }
            load.close();
        }
        else{
            log<< "Can not open load file for detector Hamiltonian."<<endl;
            MPI_Abort(MPI_COMM_WORLD,-11);
        }
    }
    // tell each process total matrix number and size
    for (m = 0; m < stlnum; m++) {
        MPI_Bcast(&total_dmat_size[m], 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&total_dmat_off_num[m], 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&total_dmat_num[m], 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    // tell each process the number of dmat_size, dmat_off_num, dmat_num in each of them
    for (m = 0; m < stlnum; m++) {
        MPI_Scatter(&dmatsize_each_process[m][0], 1, MPI_INT, &dmatsize[m], 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(&doffnum_each_process[m][0], 1, MPI_INT, &doffnum[m], 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(&dmatnum_each_process[m][0], 1, MPI_INT, &dmatnum[m], 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    // record dmatsize_each_process, doffnum_each_process and dmatnum_each_process;
    for(m=0;m<stlnum;m++){
        MPI_Bcast(&dmatsize_each_process[m][0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&doffnum_each_process[m][0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&dmatnum_each_process[m][0],num_proc,MPI_INT,0,MPI_COMM_WORLD);
    }
    // compute dmat_offset_each_process
    for(m=0;m<stlnum;m++){
        dmat_offset_each_process[m][0]=0;
        for(i=1;i<num_proc;i++){
            dmat_offset_each_process[m][i]= dmat_offset_each_process[m][i-1] + dmatnum_each_process[m][i-1];
        }
    }
    // scatter dmat, dirow, dicol to each process.
    for(m=0;m<stlnum;m++){
        dmat[m].reserve(dmatnum[m]);
        dirow[m].reserve(dmatnum[m]);
        dicol[m].reserve(dmatnum[m]);
    }
    for(m=0;m<stlnum;m++){
        temp_dmat[m] = new double [dmatnum[m]];
        temp_dirow[m] = new int [dmatnum[m]];
        temp_dicol[m] = new int [dmatnum[m]];
    }
    for(m=0;m<stlnum;m++){
        MPI_Scatterv(&total_dmat[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_DOUBLE,&temp_dmat[m][0],dmatnum[m],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Scatterv((void *)&total_dirow[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_INT,&temp_dirow[m][0],dmatnum[m],MPI_INT,0,MPI_COMM_WORLD);
        MPI_Scatterv( &total_dicol[m][0], dmatnum_each_process[m],dmat_offset_each_process[m],MPI_INT,&temp_dicol[m][0],dmatnum[m],MPI_INT,0,MPI_COMM_WORLD);
        for(i=0;i<dmatnum[m];i++){
            dmat[m].push_back(temp_dmat[m][i]);
            dirow[m].push_back(temp_dirow[m][i]);
            dicol[m].push_back(temp_dicol[m][i]);
        }
    }
    if(my_id!=0){ // allocate space for total_dmat
        for(m=0;m<stlnum;m++){
            total_dmat[m] = new double [total_dmat_num[m]];
            total_dirow[m] = new int [total_dmat_num[m]];
            total_dicol[m] = new int [total_dmat_num[m]];
        }
    }
    for(m=0;m<stlnum;m++){
        MPI_Bcast(&total_dmat[m][0],total_dmat_num[m],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&total_dirow[m][0],total_dmat_num[m],MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&total_dicol[m][0],total_dmat_num[m],MPI_INT,0,MPI_COMM_WORLD);
    }
    Scatter_dv(total_dmat_size);  // Scatter dv_all to dv.
    Broadcast_dv_all(); // broadcast dv_all from process 0 to other process.
    for(m=0;m<stlnum;m++){
        delete [] temp_dmat[m];
        delete [] temp_dirow[m];
        delete [] temp_dicol[m];
    }
    delete [] temp_dmat;
    delete [] temp_dirow;
    delete [] temp_dicol;

}

void detector:: save_detector_state_MPI(string path,double * final_time,ofstream & log,int initial_state_choice){
    int m,i;
    int nearby_state_list_size = nearby_state_index.size();
    int state_for_average_size = states_for_4_point_correlation_average.size();
    //---------------------------------------------------------------------------------
    gather_xd_yd(); // gather xd yd to xd_all, yd_all.
    if(my_id==0){
        ofstream save;
        save.open(path + "save_detector_state.txt");
        if(save.is_open()){
            save << nearby_state_list_size <<endl;
            for(i=0;i<nearby_state_list_size;i++){
                save << nearby_state_index[i] <<" ";
            }
            save << endl;

            save<< state_for_average_size << endl;
            for(i=0;i<state_for_average_size;i++){
                save << states_for_4_point_correlation_average[i] << " ";
            }
            save<<endl;

            save<<final_time[0]<<endl;
            save<<"Real part of Detector state :"<<endl;
            for(m=0;m<nearby_state_list_size;m++) {
                for (i = 0; i < total_dmat_size[0]; i++) {
                    save << xd_all[m][i] << " ";
                }
                save << endl;
            }
            save<<"Imaginary part of Detector state: "<<endl;
            for(m=0;m< nearby_state_list_size ;m++){
                for(i=0;i<total_dmat_size[0];i++){
                    save<<yd_all[m][i]<<" ";
                }
                save<<endl;
            }
            save.close();
        }
        else{
            log<<"Error.  We can not open file to save detector wave function."<<endl;
            MPI_Abort(MPI_COMM_WORLD,-12);
        }

    }
}

void detector:: load_detector_state_MPI(string path,double * start_time,ofstream & log, int initial_state_choice){
    int m,i,j;
    ifstream load;
    string ss;
    int state_index;
    int nearby_state_list_size;
    int state_for_average_size;
    int initial_state_index_in_total_dmatrix;
    initial_state_index_in_total_dmatrix = initial_state_index[0]
                                           + total_dmat_size[0]/num_proc * initial_state_pc_id[0] ;

    if(my_id==0) {
        //--------------Allocate space for detector state number for each process. -------------------
        load.open(path + "save_detector_state.txt");
        //-------------------Load data from file ----------------------------
        if (load.is_open()) {
            load >> nearby_state_list_size;
            std::getline(load,ss);
            for(m=0;m<nearby_state_list_size;m++){
                load >> state_index;
                nearby_state_index.push_back(state_index);
            }
            std::getline(load,ss);
            load >> state_for_average_size;
            std::getline(load,ss);
            for(m=0;m<state_for_average_size;m++){
                load >> state_index;
                states_for_4_point_correlation_average.push_back(state_index);
            }
            std::getline(load,ss);
        }
            //------------------------------------------------------------------
        else {
            log << "Error.  We can not open file to save detector wave function." << endl;
            MPI_Abort(MPI_COMM_WORLD, -12);
        }
    }
    // Broadcast nearby_state_index list
    MPI_Bcast(&nearby_state_list_size,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&state_for_average_size,1,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id !=0){
        nearby_state_index.resize(nearby_state_list_size);
        states_for_4_point_correlation_average.resize(state_for_average_size);
    }
    MPI_Bcast(&nearby_state_index[0],nearby_state_list_size,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&states_for_4_point_correlation_average[0],state_for_average_size,MPI_INT,0,MPI_COMM_WORLD);

    // compute initial_state_index_in_state_index_list
    for (i = 0; i < nearby_state_list_size; i++) {
        if (initial_state_index_in_total_dmatrix == nearby_state_index[i]) {
            initial_state_index_in_state_index_list = i;
            break;
        }
    }

    // compute states_for_average_in_nearby_state_index_list
    for(j=0;j<state_for_average_size;j++){
        for(i=0;i<nearby_state_list_size;i++){
            if(nearby_state_index[i] == states_for_4_point_correlation_average[j]){
                states_for_average_in_nearby_state_index_list.push_back(i);
                break;
            }
        }
    }

    xd_all = new double * [nearby_state_list_size];
    yd_all = new double * [nearby_state_list_size];
    for(i=0;i<nearby_state_list_size;i++){
        xd_all[i] = new double [total_dmat_size[0]];
        yd_all[i] = new double [total_dmat_size[0]];
    }
    // for 4 - point correlation function we will make xd , yd as N*N system.
    // we use first detector's state space
    xd = new vector <double> [nearby_state_list_size];
    yd = new vector<double> [nearby_state_list_size];
    for (i = 0; i < nearby_state_list_size ; i++) {
        xd[i].reserve(dmatsize[0]);
        yd[i].reserve(dmatsize[0]);
    }

    if(my_id == 0){
        load >> start_time[0];
        std::getline(load, ss);
        std::getline(load, ss);
        for(m=0;m<nearby_state_list_size;m++) {
            for (i = 0; i < total_dmat_size[0]; i++) {
                load >> xd_all[m][i];
            }
            std::getline(load, ss);
        }
        std::getline(load, ss);
        for(m=0;m< nearby_state_list_size ;m++) {
            for (i = 0; i < total_dmat_size[0]; i++) {
                load >> yd_all[m][i];
            }
        }
        load.close();
    }


    MPI_Bcast(&start_time[0],1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    Scatter_xd_yd();
}

