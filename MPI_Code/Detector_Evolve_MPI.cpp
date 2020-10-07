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
    recv_xd= new double * [total_dmat_size[0]];
    recv_yd= new double * [total_dmat_size[0]];
    send_xd= new double * [total_dmat_size[0]];
    send_yd = new double *[total_dmat_size[0]];

    vsize= total_dmat_size[0]/num_proc;
    to_recv_buffer_len[0] = construct_receive_buffer_index(remoteVecCount[0],remoteVecPtr[0],
            remoteVecIndex[0],0);  // construct buffer to receive.
    to_send_buffer_len[0]= construct_send_buffer_index(remoteVecCount[0],remoteVecPtr[0],remoteVecIndex[0],
                                                    tosendVecCount[0],tosendVecPtr[0], tosendVecIndex[0]);
    for(m=0;m< total_dmat_size[0];m++){
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

void detector::update_dx(){
    int i;
    int m;
    int vsize;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;
    for(m=0;m<total_dmat_size[0];m++){
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
void detector::update_dy(){
    int i;
    int vsize;
    int m;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;
    for(m=0;m<total_dmat_size[0];m++){
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
    update_dy();
        // go simulation starting from different state
    for(i=0;i<dmatnum[0];i++){
        // make sure when we compute off-diagonal matrix, we record both symmetric and asymmetric part
        irow = local_dirow[0][i];
        icol = local_dicol[0][i]; // compute to point to colindex in
        for(m=0;m<total_dmat_size[0];m++) {
            xd[m][irow] = xd[m][irow] + dmat[0][i] * yd[m][icol] * cf;
        }
    }

    update_dx();
    for(i=0;i<dmatnum[0];i++){
        irow= local_dirow[0][i];
        icol= local_dicol[0][i];
        for(m=0;m<total_dmat_size[0];m++){
            yd[m][irow] = yd[m][irow] - dmat[0][i] * xd[m][icol] * cf;
        }
    }

}

void detector:: compute_n_off_diag_element(int index_b, complex<double> * n_off_diag_element,
                                           int initial_state_index_in_total_dmatrix){
    // compute <a|n_{i}(t)|b> for i = 1, 2, .., N (number of mode) mode
    // here a is our state of interest: We want to compute <a | [n_{i}(t) , n_{i}]|^{2} a>.
    // this index of a is given as initial_state_index_in_total_dmatrix. (state a is chosen as initial state in detector 0)
    // result store in n_off_diag_element
    int i,k;
    complex<double> n_k_CC;
    complex<double> n_k_CC_i; // n_k_CC : n_{i}(k) * C^{b}_{k}(t) * (C^{a}_{k}(t))^{*}. Here C_{b}^{k}(t) = <k|b(t)>
    double nk_i; // n_i_k is n_i(k) which is number of quanta in mode i for state k
    complex<double> C_kb; //   C_{b}^{k}(t) = <k|b(t)>
    complex<double> C_ka_conjugate; //  (C^{a}_{k}(t))^{*}
    for(i=0;i<nmodes[0];i++) {
        n_off_diag_element[i] = 0;
    }

    for(k=0;k<dmatsize[0];k++){
        C_ka_conjugate = complex<double> (xd[initial_state_index_in_total_dmatrix][k],
                                -yd[initial_state_index_in_total_dmatrix][k]);
        C_kb = complex<double> (xd[index_b][k],
                                yd[index_b][k]);
        n_k_CC = C_kb * C_ka_conjugate;
        for(i=0;i<nmodes[0];i++){
            nk_i = dv[0][k][i];
            n_k_CC_i= n_k_CC * nk_i;
            n_off_diag_element[i] = n_off_diag_element[i] + n_k_CC_i;
        }
    }

}

void full_system::pre_coupling_evolution_MPI(int initial_state_choice){
    // we do not write output function now, try to make it as simple as possible.
    int irow_index, icol_index;
    int m,i,j,k;
    int b;
    double var;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    int steps;
    double detector_tprint = 0.01;
    int output_step= int(detector_tprint/delt); //Output every output_step.
    double * start_time = new double [s.tlnum];
    double * final_time = new double [s.tlnum];
    for(i=0;i<s.tlnum;i++){
        start_time[i]=0;
    }
    ofstream four_point_correlation_output;
    ofstream Detector_precoup_output;

    int initial_state_index_in_total_dmatrix;
    initial_state_index_in_total_dmatrix = d.initial_state_index[0]
            + d.total_dmat_size[0]/num_proc * d.initial_state_pc_id[0] ;
    // -----------Open 4_point_correlation_output ofstream -----------------------------
    if(my_id==0){
        if(Detector_Continue_Simulation){
            four_point_correlation_output.open(path + "4 point correlation.txt", ofstream::app);
            Detector_precoup_output.open(path + "detector_precoupling.txt", ofstream::app);
        }
        else {
            four_point_correlation_output.open(path + "4 point correlation.txt");
            Detector_precoup_output.open(path + "detector_precoupling.txt");
        }
    }
    // -------------- Allocate space for <a| |[n_{i}(t),n_{i}(0)]|^{2} |a> -------------
    double * four_point_correlation_function_at_initial_state = new double [d.nmodes[0]];

    //---------- Allocate space for <a| n_{i}(t) |b>  size: nmode * total_dmat_size[0] -------------------------------------------
    // each row is one n_{i}.  each column is one site |b>
    complex<double> ** n_offdiag_total;
    double ** n_offdiag_total_real;
    double ** n_offdiag_total_imag;
    n_offdiag_total = new complex<double> * [d.nmodes[0]];
    n_offdiag_total_real = new double * [d.nmodes[0]];
    n_offdiag_total_imag = new double * [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag_total[i] = new complex <double> [d.total_dmat_size[0]];
        n_offdiag_total_real[i] = new double [d.total_dmat_size[0]];
        n_offdiag_total_imag[i] = new double [d.total_dmat_size[0]];
    }
    complex<double> ** n_offdiag;
    double ** n_offdiag_real;
    double ** n_offdiag_imag;
    n_offdiag = new complex <double> * [d.nmodes[0]];
    n_offdiag_real = new double * [d.nmodes[0]];
    n_offdiag_imag = new double * [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag [i] = new complex <double> [d.total_dmat_size[0]];
        n_offdiag_real[i] = new double [d.total_dmat_size[0]];
        n_offdiag_imag[i] = new double [d.total_dmat_size[0]];
    }
    complex<double> * n_offdiag_element;
    n_offdiag_element = new complex <double> [d.nmodes[0]];
    // ------- For output for our state <a(t)|a> trajectory (check ergodicity)
    double special_state_x;
    double special_state_y;
    double survival_prob;
    // -------------Load detector state from save data if we want to continue simulation of detector.------------------
    if(Detector_Continue_Simulation){
        d.load_detector_state_MPI(path,start_time,log,initial_state_choice);
    }
    else{
        d.initialize_detector_state_MPI(log, 0); // initialize detector lower bright state
    }


    // -----------------------------------------------------------------------------------------------
    // prepare sendbuffer and recv_buffer and corresponding index.
    d.prepare_evolution();
    final_time[0] = 0;
    if(d.proptime[0]>0){
        if(my_id==0){
            log <<" Begin Propogation in Molecule Phase space "<<endl;
        }
        t=0;
        steps= d.proptime[0]/delt;
        if(Detector_Continue_Simulation){
            // update simulation time and t for detector to begin simulation.
            t=start_time[0];
            steps= (d.proptime[0]-start_time[0])/delt;
        }
        if(!Detector_Continue_Simulation) {
            if (my_id == 0) {
                four_point_correlation_output << "4- point correlation function for molecule " << endl;
                four_point_correlation_output << "total time: " << (d.proptime[0] - start_time[0]) << " "
                                              << delt * output_step << endl;
            }
        }
        // Do simulation in loop
        for(k=0;k<=steps;k++){
            //-------------------- output result ----------------------------
            if(k % output_step ==0) {
                // ---------------------------------output 4-point correlation function -------------------------------------------------
                for(b=0;b<d.total_dmat_size[0];b++){
                    // compute <b|n_{i}(t)|a>
                    // go through all state b
                    d.compute_n_off_diag_element(b,n_offdiag_element,initial_state_index_in_total_dmatrix);
                    for(i=0;i<d.nmodes[0];i++){
                        n_offdiag[i][b] = n_offdiag_element[i];
                        n_offdiag_real[i][b] = real(n_offdiag[i][b]);
                        n_offdiag_imag[i][b] = imag(n_offdiag[i][b]);
                    }
                }
                for(i=0;i<d.nmodes[0];i++){
                    MPI_Allreduce(&n_offdiag_real[i][0], & n_offdiag_total_real[i][0], d.total_dmat_size[0], MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                    MPI_Allreduce(&n_offdiag_imag[i][0], & n_offdiag_total_imag[i][0], d.total_dmat_size[0], MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
                }
                for(b=0;b<d.total_dmat_size[0];b++){
                    for(i=0;i<d.nmodes[0];i++){
                        n_offdiag_total[i][b] = complex<double> (n_offdiag_total_real[i][b], n_offdiag_total_imag[i][b]);
                    }
                }
                for(i=0;i<d.nmodes[0];i++){
                    four_point_correlation_function_at_initial_state[i] = 0;
                    for(b=0;b<d.total_dmat_size[0];b++){
                        // var = |<a|n_{i}(t)|b>|^{2} * (n_{i}(b) - n_{i}(a))^{2}
                       var = std::norm(n_offdiag_total[i][b]) * pow(d.dv_all[0][b][i] - d.dv_all[0][initial_state_index_in_total_dmatrix][i],2);
                       four_point_correlation_function_at_initial_state[i] =
                               four_point_correlation_function_at_initial_state[i] + var;
                    }
                }
                if(my_id == 0){
                    four_point_correlation_output << "Time:   " << t << endl;
                    for(i=0;i<d.nmodes[0];i++){
                        four_point_correlation_output << four_point_correlation_function_at_initial_state[i] << " ";
                    }
                    four_point_correlation_output<<endl;
                }
                // ------------- output state <a(t)|a>  ----------------------------------------
                if (my_id == d.initial_state_pc_id[0]){
                    special_state_x = d.xd[initial_state_index_in_total_dmatrix][d.initial_state_index[0]];
                    special_state_y = d.yd[initial_state_index_in_total_dmatrix][d.initial_state_index[0]];
                }
                MPI_Bcast(&special_state_x,1,MPI_DOUBLE,d.initial_state_pc_id[0],MPI_COMM_WORLD);
                MPI_Bcast(&special_state_y,1,MPI_DOUBLE,d.initial_state_pc_id[1],MPI_COMM_WORLD);
                survival_prob = pow(special_state_x,2) + pow(special_state_y,2);
                if(my_id == 0){
                    Detector_precoup_output << "Time:   " << t << endl;
                    Detector_precoup_output << survival_prob <<endl;
                }
            }
            t= t+ delt;
            d.SUR_onestep_MPI(cf);
        }
        final_time[0] = t;
    }

    d.save_detector_state_MPI(path,final_time,log,initial_state_choice);


    if(my_id==0){
        cout<<"Detector_pre_coupling simulation finished"<<endl;
        four_point_correlation_output.close();
        Detector_precoup_output.close();
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
    for(i=0;i<d.total_dmat_size[0];i++){
        delete [] d.send_xd[i];
        delete [] d.send_yd[i];
        delete [] d.recv_xd[i];
        delete [] d.recv_yd[i];
    }
    for(i=0;i<d.nmodes[0];i++){
        delete [] n_offdiag_total[i];
        delete [] n_offdiag_total_real[i];
        delete [] n_offdiag_total_imag[i];
        delete [] n_offdiag[i];
        delete [] n_offdiag_real[i];
        delete [] n_offdiag_imag[i];
    }
    delete [] n_offdiag_total;
    delete [] n_offdiag_total_real;
    delete [] n_offdiag_total_imag;
    delete [] n_offdiag;
    delete [] n_offdiag_real;
    delete [] n_offdiag_imag;
    delete [] n_offdiag_element;
    delete [] four_point_correlation_function_at_initial_state;

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

};
