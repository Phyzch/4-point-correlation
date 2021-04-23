//
// Created by phyzch on 6/27/20.
//
#include "../system.h"
#include "../util.h"
using namespace std;

int  detector::construct_receive_buffer_index(int * remoteVecCount_element, int * remoteVecPtr_element, int * remoteVecIndex_element, int detector_index){
    // input: remoteVecCount: total number of element need to receive from each process.
    //        remoteVecPtr: displacement in remoteVecIndex for element in each process.
    //        remoteVecIndex: index for remote vector we need to receive. (they may allocate in different remote process.)
    // return: length of remoteVecIndex.

    int i,j;
    // range for element in process is [local_begin, local_end)
    int total_remoteVecCount=0;
    int vsize = total_dmat_size[detector_index] / num_proc;
    int local_begin= total_dmat_size[detector_index] / num_proc * my_id;
    int local_end;
    int remote_pc_id;
    if(my_id!=num_proc-1) {
        local_end = total_dmat_size[detector_index] / num_proc * (my_id + 1);
    }
    else{
        local_end = total_dmat_size[detector_index];
    }
    // ---------------------------------------------------------------
    vector <int> col_index_copy = dicol[detector_index];
    sort(col_index_copy.begin(),col_index_copy.end()); // sort vector.
    int col_array_size = col_index_copy.size();
    int prev_col=-1;
    j=0;
    for(i=0;i<col_array_size;i++){
        if( (col_index_copy[i]>prev_col)   and ( (col_index_copy[i]<local_begin) or (col_index_copy[i] >= local_end) )  ){
            // this matrix element is not in process.
            if (col_index_copy[i] >= vsize * (num_proc-1) ){
                remote_pc_id = num_proc-1;
            }
            else{
                remote_pc_id = col_index_copy[i] / vsize;
            }
            remoteVecCount_element[remote_pc_id] ++;
            remoteVecIndex_element [j] = col_index_copy[i];  // vector index need to receive. (global index , ordered)
            j++;
        }
        prev_col= col_index_copy[i];
    }
    remoteVecPtr_element[0]=0;   // displacement for remote vector from each process in remoteVecIndex.
    for(i=1;i<num_proc;i++){
        remoteVecPtr_element[i] = remoteVecPtr_element[i-1] + remoteVecCount_element[i-1];
    }
    for(i=0;i<num_proc;i++){
        total_remoteVecCount = total_remoteVecCount + remoteVecCount_element[i];
    }
    return total_remoteVecCount;
}

int construct_send_buffer_index(int * remoteVecCount_element, int * remoteVecPtr_element, int * remoteVecIndex_element,
                                int * tosendVecCount_element, int * tosendVecPtr_element, int* & tosendVecIndex_ptr){
    //  tosend_Vec_count record number of element to send to each process.
    // tosend_Vec_Index record the global index of vector the process have to send
    // tosend_Vec_Ptr record the offset of vector to send to each other process.
    // return to_send_buffer_len: lenfth of tosendVecIndex
    int i;
    int to_send_buffer_len;

    MPI_Alltoall(&remoteVecCount_element[0],1,MPI_INT,&(tosendVecCount_element[0]),1,MPI_INT,MPI_COMM_WORLD);

    // compute displacement for each process's data.
    tosendVecPtr_element[0]=0;
    for(i=1;i<num_proc;i++){
        tosendVecPtr_element[i]= tosendVecPtr_element[i-1] + tosendVecCount_element[i-1];
    }
    // compute total length of buffer to send
    to_send_buffer_len=0;
    for(i=0;i<num_proc;i++){
        to_send_buffer_len= to_send_buffer_len + tosendVecCount_element[i];
    }
    // Index (in global) of element to send. use MPI_Alltoallv to receive the index to send.
    tosendVecIndex_ptr = new int [to_send_buffer_len];
    MPI_Alltoallv(&remoteVecIndex_element[0],remoteVecCount_element,remoteVecPtr_element,MPI_INT,
                  & tosendVecIndex_ptr[0],tosendVecCount_element,tosendVecPtr_element,MPI_INT,MPI_COMM_WORLD);

    return to_send_buffer_len;
}

int compar(const void * a, const void * b){
    return *(int *) a - * (int *)b;
}

// this function is called every time we do evolution
void detector::prepare_evolution(){
    // compute buffer to receive and send for each process.
    // resize xd,yd to provide extra space for recv_buffer.
    // allocate space for send_xd , send_yd buffer.
    // Index for remoteVecIndex, tosendVecIndex are computed here.
    int m,i;
    int vsize;
    int sampling_state_list_size = sampling_state_index.size();
    // Index for vector to send and receive.
    // remoteVecCount: total number to receive. remoteVecPtr: displacement in remoteVecIndex for each process. remoteVecIndex: index in other process to receive.
    // tosendVecCount: total number to send to other process. tosendVecPtr: displacement in tosendVecIndex in each process.  tosendVecIndex: Index of element in itself to send. (it's global ,need to be converted to local index)
    remoteVecCount= new int * [1];
    remoteVecPtr= new int * [1];
    remoteVecIndex= new int * [1];
    to_recv_buffer_len = new int  [1];

    //------------------Allocate space for vector to receive ---------------------
    for (m=0;m<1;m++){
        remoteVecCount[m] = new int [num_proc];
        remoteVecPtr[m] = new int [num_proc];
        remoteVecIndex[m] = new int [dmatnum[m]];
        for(i=0;i<num_proc;i++){
            remoteVecCount[m][i] = 0;
        }
    }
    tosendVecCount= new int *[1];
    tosendVecPtr =  new int * [1];
    tosendVecIndex = new int * [1];
    to_send_buffer_len= new int [1];
    for(m=0;m<1;m++){
        tosendVecCount[m] = new int [num_proc];
        tosendVecPtr[m] = new int [num_proc];
    }

    int * search_Ind; // local variable, used for compute local_dicol;
    int col_index_to_search;
    // local column index used when we do H *x and H*y
    local_dirow= new vector<int> [1];
    local_dicol = new vector<int> [1]; // column index for computation in local matrix.
    // buffer to send and receive buffer to/from other process.
    recv_xd= new double * [sampling_state_list_size];
    recv_yd= new double * [sampling_state_list_size];
    send_xd= new double * [sampling_state_list_size];
    send_yd = new double *[sampling_state_list_size];

    vsize= total_dmat_size[0]/num_proc;
    to_recv_buffer_len[0] = construct_receive_buffer_index(remoteVecCount[0],remoteVecPtr[0],
            remoteVecIndex[0],0);  // construct buffer to receive.
    to_send_buffer_len[0]= construct_send_buffer_index(remoteVecCount[0],remoteVecPtr[0],remoteVecIndex[0],
                                                    tosendVecCount[0],tosendVecPtr[0], tosendVecIndex[0]);
    for(m=0;m< sampling_state_list_size;m++){
        xd[m].resize(dmatsize[0] + to_recv_buffer_len[0]);
        yd[m].resize(dmatsize[0] + to_recv_buffer_len[0]);
        recv_xd[m] = new double [to_recv_buffer_len[0]];
        recv_yd[m]= new double [to_recv_buffer_len[0]];
        send_xd[m]= new double [to_send_buffer_len[0]];
        send_yd[m] = new double [to_send_buffer_len[0]];
    }
    // construct local_dirow, local_dicol
    local_dirow[0].reserve(dmatnum[0]);
    local_dicol[0].reserve(dmatnum[0]);
    for(i=0;i<dmatnum[0];i++){
        local_dirow[0].push_back(dirow[0][i] - my_id * vsize);  // set local index for row index
        col_index_to_search= dicol[0][i];
        search_Ind=(int *) bsearch(&col_index_to_search,remoteVecIndex[0],to_recv_buffer_len[0],sizeof(int),compar);
        if(search_Ind!=NULL){
            // this column index is not in local matrix, and we should get it from other process (remoteVec)
            local_dicol[0].push_back(dmatsize[0] + (search_Ind-remoteVecIndex[0]) );
        }
        else{ // this column index is in local matrix.
            local_dicol[0].push_back (dicol[0][i] - my_id * vsize );
        }
    }


}

void detector::update_dx(int nearby_state_list_size){
    int i;
    int m;
    int vsize;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;

    for(m=0;m<nearby_state_list_size;m++){
        for (i = 0; i < to_send_buffer_len[0]; i++) {
            send_xd[m][i] = xd[m][tosendVecIndex[0][i] - my_id * vsize];
        }
        MPI_Alltoallv(&send_xd[m][0],tosendVecCount[0],tosendVecPtr[0],MPI_DOUBLE,
                      &recv_xd[m][0],remoteVecCount[0],remoteVecPtr[0],MPI_DOUBLE,MPI_COMM_WORLD);
        for(i=0;i<to_recv_buffer_len[0];i++){
            xd[m][i+ dmatsize[0]]= recv_xd[m][i];
        }
    }
}
void detector::update_dy(int nearby_state_list_size){
    int i;
    int vsize;
    int m;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;
    for(m=0;m<nearby_state_list_size;m++){
        for (i = 0; i < to_send_buffer_len[0]; i++) {
            send_yd[m][i] = yd[m][tosendVecIndex[0][i] - my_id * vsize];
        }
        MPI_Alltoallv(&send_yd[m][0],tosendVecCount[0],tosendVecPtr[0],MPI_DOUBLE,
                      &recv_yd[m][0],remoteVecCount[0],remoteVecPtr[0],MPI_DOUBLE,MPI_COMM_WORLD);
        for(i=0;i<to_recv_buffer_len[0];i++){
            yd[m][i+ dmatsize[0]]= recv_yd[m][i];
        }
    }

}


// call it every time we do SUR algorithm. do it for detector specified by detector_index.
void detector::SUR_onestep_MPI(double cf){
    int m, i;
    int irow,icol;
    int nearby_state_list_size = sampling_state_index.size();
    // do simulation starting from different state
    update_dy(nearby_state_list_size);
    for(i=0;i<dmatnum[0];i++){
        // make sure when we compute off-diagonal matrix, we record both symmetric and asymmetric part
        irow = local_dirow[0][i];
        icol = local_dicol[0][i]; // compute to point to colindex in
        for(m=0;m<nearby_state_list_size;m++) {
            xd[m][irow] = xd[m][irow] + dmat[0][i] * yd[m][icol] * cf;
        }
    }

    update_dx(nearby_state_list_size);
    for(i=0;i<dmatnum[0];i++){
        irow= local_dirow[0][i];
        icol= local_dicol[0][i];
        for(m=0;m<nearby_state_list_size;m++){
            yd[m][irow] = yd[m][irow] - dmat[0][i] * xd[m][icol] * cf;
        }
    }

}

void full_system::pre_coupling_evolution_MPI(int initial_state_choice){
    // we do not write output function now, try to make it as simple as possible.
    int irow_index, icol_index;
    int start_index;
    int m,i,j,k;


    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    int steps;
    double detector_tprint = tprint; // set detector_tprint == tprint specified in input file.
    int output_step= int(detector_tprint/delt); //Output every output_step.


    double * start_time = new double [s.tlnum];
    double * final_time = new double [s.tlnum];

    int sampling_state_index_size = d.sampling_state_index.size();

    // ----------------- variable for IPR --------------
    double IPR;

    double * inverse_IPR_in_one_process = new double [sampling_state_index_size]  ;
    double * inverse_IPR_all = new double [sampling_state_index_size];
    double * IPR_all = new double [sampling_state_index_size];

    double * normalization_in_one_process = new double [sampling_state_index_size];
    double * normalization_in_all_process = new double [sampling_state_index_size];

    for(i=0;i<s.tlnum;i++){
        start_time[i]=0;
    }
    // ------------ initialize state 's wave function ----------------
     d.initialize_detector_state_MPI(log); // initialize detector lower bright state


    int * initial_state_index_in_total_dmatrix = new int [2];
    for(i=0;i<2;i++){
        initial_state_index_in_total_dmatrix[i] = d.initial_state_index[i]
                                               + d.total_dmat_size[0]/num_proc * d.initial_state_pc_id[i] ;
    }

    // -----------------------------------------------------------------------------------------------
    // prepare sendbuffer and recv_buffer and corresponding index.
    d.prepare_evolution();

    // // use ground electronic state specified in input.txt to compute energy and shift Hamiltonian according to it.
    d.shift_detector_Hamiltonian(log);


    ofstream IPR_for_all_state_output;
    if(my_id == 0){
        IPR_for_all_state_output.open(path+"IPR_all_state.txt");
        IPR_for_all_state_output << sampling_state_index_size << endl;
        // output energy
        for(i=0;i<sampling_state_index_size ; i++ ){
            IPR_for_all_state_output << dmat0[ d.sampling_state_index[i] ] << " ";
        }
        IPR_for_all_state_output << endl;

        // output mode quanta
        for(i=0;i<sampling_state_index_size;i++){
            for(j=0;j<d.nmodes[0];j++){
                IPR_for_all_state_output << d.dv_all[0][d.sampling_state_index[i]][j] <<"  ";
            }
            IPR_for_all_state_output << endl;
        }
    }

    final_time[0] = 0;
    if(d.proptime[0]>0){
        if(my_id==0){
            log <<" Begin Propogation in Molecule Phase space "<<endl;
            cout << " Begin Propogation in Molecule Phase space "<<endl;
        }
        t=0;
        steps= d.proptime[0]/delt+1;
        if(Detector_Continue_Simulation){
            // update simulation time and t for detector to begin simulation.
            t=start_time[0];
            steps= (d.proptime[0]-start_time[0])/delt+2;
        }

        if(Detector_Continue_Simulation){
            start_index=1;
        }
        else{
            start_index = 0;
        }

        // Do simulation in loop
        for(k=start_index;k<steps;k++){
            // ------------- Output result ------------------------

            // --------- output IPR for all state --------------------------------------
            if(k% output_step == 0){
                for(i=0;i< sampling_state_index_size;i++){
                    inverse_IPR_in_one_process[i] = 0;
                    normalization_in_one_process[i] = 0;
                    for(j=0;j<d.dmatsize[0];j++){
                        inverse_IPR_in_one_process[i] = inverse_IPR_in_one_process[i] + pow(pow(d.xd[i][j],2) + pow(d.yd[i][j],2),2);
                        normalization_in_one_process[i] = normalization_in_one_process[i] +
                                pow(d.xd[i][j],2) + pow(d.yd[i][j],2);
                    }
                }

                MPI_Reduce(&inverse_IPR_in_one_process[0],&inverse_IPR_all[0],sampling_state_index_size,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
                MPI_Reduce(&normalization_in_one_process[0], &normalization_in_all_process[0], sampling_state_index_size, MPI_DOUBLE, MPI_SUM, 0 , MPI_COMM_WORLD);
                if(my_id == 0){
                    for(i=0;i<sampling_state_index_size;i++){
                        inverse_IPR_all[i] = inverse_IPR_all[i] / normalization_in_all_process[i];
                    }

                    for(i=0;i<sampling_state_index_size;i++){
                        IPR_all[i] = 1/ inverse_IPR_all[i];
                    }

                    // output result
                    IPR_for_all_state_output << t << endl;
                    for(i=0;i<sampling_state_index_size;i++){
                        IPR_for_all_state_output << IPR_all[i] <<" ";
                    }
                    IPR_for_all_state_output << endl;
                }

            }

            t= t+ delt;
            d.SUR_onestep_MPI(cf);
        }
        final_time[0] = t;
    }


    if(my_id==0){
        cout<<"Detector_pre_coupling simulation finished"<<endl;
        IPR_for_all_state_output.close();
    }
    // -------------- free remote_Vec_Count, remote_Vec_Index -------------------------
    for(i=0;i<1;i++){
        delete [] d.remoteVecCount[i];
        delete [] d.remoteVecPtr[i];
        delete []  d.remoteVecIndex[i];
        delete [] d.tosendVecCount[i];
        delete [] d.tosendVecPtr[i];
        delete [] d.tosendVecIndex[i];
    }
    for(i=0;i<sampling_state_index_size;i++){
        delete [] d.send_xd[i];
        delete [] d.send_yd[i];
        delete [] d.recv_xd[i];
        delete [] d.recv_yd[i];
    }
//
    delete [] d.to_recv_buffer_len;
    delete [] d.remoteVecCount;
    delete [] d.remoteVecPtr;
    delete [] d.remoteVecIndex;
    delete [] d.to_send_buffer_len;
    delete [] d.tosendVecCount;
    delete [] d.tosendVecPtr;
    delete [] d.tosendVecIndex;
    delete [] d.send_xd;
    delete [] d.send_yd;
    delete [] d.recv_xd;
    delete [] d.recv_yd;
    delete [] d.local_dirow;
    delete [] d.local_dicol;

    delete [] start_time;
    delete [] final_time;

    delete [] IPR_all;
    delete [] inverse_IPR_all;
    delete [] inverse_IPR_in_one_process;

    delete [] normalization_in_one_process;
    delete [] normalization_in_all_process;
};
