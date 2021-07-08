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
    recv_xd= new double * [state_number_for_evolution];
    recv_yd= new double * [state_number_for_evolution];
    send_xd= new double * [state_number_for_evolution];
    send_yd = new double *[state_number_for_evolution];

    vsize= total_dmat_size[0]/num_proc;
    to_recv_buffer_len[0] = construct_receive_buffer_index(remoteVecCount[0],remoteVecPtr[0],
            remoteVecIndex[0],0);  // construct buffer to receive.
    to_send_buffer_len[0]= construct_send_buffer_index(remoteVecCount[0],remoteVecPtr[0],remoteVecIndex[0],
                                                    tosendVecCount[0],tosendVecPtr[0], tosendVecIndex[0]);
    for(m=0;m< state_number_for_evolution ;m++){
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

void detector::update_dx(int state_number_for_evolution){
    int i;
    int m;
    int vsize;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;

    for(m=0;m< state_number_for_evolution ;m++){
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
void detector::update_dy(int state_number_for_evolution){
    int i;
    int vsize;
    int m;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;
    for(m=0;m< state_number_for_evolution ;m++){
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
    // do simulation starting from different state
    update_dy(state_number_for_evolution);
    for(i=0;i<dmatnum[0];i++){
        // make sure when we compute off-diagonal matrix, we record both symmetric and asymmetric part
        irow = local_dirow[0][i];
        icol = local_dicol[0][i]; // compute to point to colindex in
        for(m=0;m< state_number_for_evolution ;m++) {
            xd[m][irow] = xd[m][irow] + dmat[0][i] * yd[m][icol] * cf;
        }
    }

    update_dx(state_number_for_evolution);
    for(i=0;i<dmatnum[0];i++){
        irow= local_dirow[0][i];
        icol= local_dicol[0][i];
        for(m=0;m< state_number_for_evolution ;m++){
            yd[m][irow] = yd[m][irow] - dmat[0][i] * xd[m][icol] * cf;
        }
    }

}

void detector:: compute_n_off_diag_element(int index_b, int index_a, complex<double> * n_off_diag_element){
    // compute <a|n_{i}(t)|b> for i = 1, 2, .., N (number of mode) mode.
    // Sum over all state k: <a(t)|k> n_{k}<k|b(t)>. We have to first sum it over in one process, then sum over all process.
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
        C_ka_conjugate = complex<double> (xd[index_a][k],
                                -yd[index_a][k]);
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

void full_system::compute_4_point_corre_for_single_state(int nearby_state_index_size, complex<double> * n_offdiag_element,
                                                         complex<double> ** n_offdiag,double ** n_offdiag_real, double ** n_offdiag_imag,
                                                         complex<double> **n_offdiag_total, double ** n_offdiag_total_real, double ** n_offdiag_total_imag,
                                                         int initial_state_index_in_total_dmatrix ,
                                                         double * four_point_correlation_function_at_initial_state){
    int i, b;
    double var;
    for(b=0;b<nearby_state_index_size;b++){
        // compute <b|n_{i}(t)|a>
        // go through all state b
        d.compute_n_off_diag_element(b,d.initial_state_index_in_nearby_state_index_list,n_offdiag_element);
        for(i=0;i<d.nmodes[0];i++){
            n_offdiag[i][b] = n_offdiag_element[i];
            n_offdiag_real[i][b] = real(n_offdiag[i][b]);
            n_offdiag_imag[i][b] = imag(n_offdiag[i][b]);
        }
    }
    for(i=0;i<d.nmodes[0];i++){
        MPI_Allreduce(&n_offdiag_real[i][0], & n_offdiag_total_real[i][0],nearby_state_index_size, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&n_offdiag_imag[i][0], & n_offdiag_total_imag[i][0], nearby_state_index_size, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    for(b=0;b<nearby_state_index_size;b++){
        for(i=0;i<d.nmodes[0];i++){
            n_offdiag_total[i][b] = complex<double> (n_offdiag_total_real[i][b], n_offdiag_total_imag[i][b]);
        }
    }
    for(i=0;i<d.nmodes[0];i++){
        four_point_correlation_function_at_initial_state[i] = 0;
        for(b=0;b<nearby_state_index_size;b++){
            // var = |<a|n_{i}(t)|b>|^{2} * (n_{i}(b) - n_{i}(a))^{2}
            // note there are two set of index: (index in nearby_state_list) <-> (index in dmat_all,dv_all)
            // (b, initial_state_index_in_nearby_state_index_list) <-> (nearby_state_index[b] , initial_state_index_in_total_dmatrix)
            // nearby_state_index record the index for all these independent state we evolve with time.

            var = std::norm(n_offdiag_total[i][b]) *
                  pow(d.dv_all[0][ d.nearby_state_index[b] ][i]
                      - d.dv_all[0][initial_state_index_in_total_dmatrix][i],2);
            four_point_correlation_function_at_initial_state[i] =
                    four_point_correlation_function_at_initial_state[i] + var;
        }
    }
}

void full_system::compute_4_point_corre_for_multiple_states(int state_for_average_size,int nearby_state_index_size,
                                                            complex<double> * n_offdiag_element,
                                                            complex<double> *** n_offdiag_for_states_ensemble, double *** n_offdiag_for_states_ensemble_real,
                                                            double *** n_offdiag_for_states_ensemble_imag,
                                                            complex<double> *** n_offdiag_total_for_states_ensemble,
                                                            double *** n_offdiag_total_for_states_ensemble_real, double *** n_offdiag_total_for_states_ensemble_imag,
                                                            double * four_point_correlation_function_average_over_states, double * four_point_correlation_function_variance_over_states,
                                                            double ** four_point_correlation_function_for_each_states){
    int a, b, i;
    double var;
    for(a=0;a<state_for_average_size;a++){
        for(b=0;b<nearby_state_index_size;b++){
            // compute <b|n_{i}(t)|a>
            // go through all state b
            d.compute_n_off_diag_element(b,d.states_for_average_in_nearby_state_index_list[a],n_offdiag_element);
            for(i=0;i<d.nmodes[0];i++){
                n_offdiag_for_states_ensemble[i][a][b] = n_offdiag_element[i];
                n_offdiag_for_states_ensemble_real[i][a][b] = real(n_offdiag_element[i]);
                n_offdiag_for_states_ensemble_imag[i][a][b] = imag(n_offdiag_element[i]);
            }
        }
        for(i=0;i<d.nmodes[0];i++){
            MPI_Allreduce(&n_offdiag_for_states_ensemble_real[i][a][0], & n_offdiag_total_for_states_ensemble_real[i][a][0],
                          nearby_state_index_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(& n_offdiag_for_states_ensemble_imag[i][a][0], & n_offdiag_total_for_states_ensemble_imag[i][a][0],
                          nearby_state_index_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

        for(i=0;i<d.nmodes[0];i++){
            for (b=0;b<nearby_state_index_size;b++){
                n_offdiag_total_for_states_ensemble[i][a][b] =
                        complex<double> (n_offdiag_total_for_states_ensemble_real[i][a][b],
                                         n_offdiag_total_for_states_ensemble_imag[i][a][b]);
            }
        }

        for(i=0;i<d.nmodes[0];i++){
            four_point_correlation_function_for_each_states[i][a] = 0;
            for (b=0;b<nearby_state_index_size;b++){
                var = std::norm(n_offdiag_total_for_states_ensemble[i][a][b]) *
                      pow(d.dv_all[0][  d.nearby_state_index[b]  ][i] -
                          d.dv_all[0][ d.states_for_4_point_correlation_average[a] ][i],2);
                four_point_correlation_function_for_each_states[i][a] =
                        four_point_correlation_function_for_each_states[i][a] + var;
            }
        }
    }

    // average over all states:
    for(i=0;i<d.nmodes[0];i++){
        four_point_correlation_function_average_over_states[i] = 0;
        for(a=0;a<state_for_average_size;a++){
            four_point_correlation_function_average_over_states[i] =
                    four_point_correlation_function_average_over_states[i] +
                    four_point_correlation_function_for_each_states[i][a];
        }
        four_point_correlation_function_average_over_states[i] =
                four_point_correlation_function_average_over_states[i]/ state_for_average_size;
    }

    // compute the variance over all states:
    for(i=0;i<d.nmodes[0];i++){
        four_point_correlation_function_variance_over_states [i] = 0;
        for(a=0;a<state_for_average_size; a++ ){
            four_point_correlation_function_variance_over_states[i] =
                    four_point_correlation_function_variance_over_states[i] +
                            pow(four_point_correlation_function_for_each_states[i][a] -four_point_correlation_function_average_over_states[i],2);
        }
        four_point_correlation_function_variance_over_states[i] = four_point_correlation_function_variance_over_states[i]/(max(state_for_average_size-1,1));
    }

}

void full_system::compute_Stability_Matrix( double ** Stability_Matrix,
                                           int nearby_state_index_size,
                                           complex<double> * n_offdiag_element, double * n_offdiag_element_real, double * n_offdiag_element_imag,
                                           complex<double> * n_offdiag_element_all_pc, double * n_offdiag_element_real_all_pc, double * n_offdiag_element_imag_all_pc,
                                           int initial_state_index_in_total_dmatrix ){
    int b;
    int i,j;
    double n_offdiag_element_norm_sum;  // \sum_{k} |<a|n_{k}(t)|b>|^{2}

    // clear previous result in Stability Matrix
    for(i=0;i<d.nmodes[0];i++){
        for(j=0;j<d.nmodes[0];j++){
            Stability_Matrix[i][j] = 0;
        }
    }

    for(b=0;b<nearby_state_index_size;b++){
        // b is state index
        d.compute_n_off_diag_element(b,d.initial_state_index_in_nearby_state_index_list,n_offdiag_element); // n_offdiag_element has <a|n_{k}(t)|b> there
        for(i=0;i<d.nmodes[0];i++){
            n_offdiag_element_real[i] = real(n_offdiag_element[i]);
            n_offdiag_element_imag[i] = imag(n_offdiag_element[i]);
        }
        MPI_Allreduce(&n_offdiag_element_real[0],&n_offdiag_element_real_all_pc[0],d.nmodes[0],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&n_offdiag_element_imag[0],&n_offdiag_element_imag_all_pc[0],d.nmodes[0],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        for(i=0;i<d.nmodes[0];i++){
            n_offdiag_element_all_pc[i] = complex<double>( n_offdiag_element_real_all_pc[i], n_offdiag_element_imag_all_pc[i] );
        }
        n_offdiag_element_norm_sum = 0;
        for(i=0;i<d.nmodes[0];i++){
            n_offdiag_element_norm_sum = n_offdiag_element_norm_sum +  std::norm( n_offdiag_element_all_pc[i] );
        }

        for(i=0;i<d.nmodes[0];i++){
            for(j=0;j<d.nmodes[0];j++){
                Stability_Matrix[i][j] = Stability_Matrix[i][j] + n_offdiag_element_norm_sum *
                        (d.dv_all[0][d.nearby_state_index[b]][i] - d.dv_all[0][initial_state_index_in_total_dmatrix][i])*
                      (d.dv_all[0][d.nearby_state_index[b]][j] - d.dv_all[0][initial_state_index_in_total_dmatrix][j]);
            }
        }

    }
}

void full_system::pre_coupling_evolution_MPI(int initial_state_choice){
    // we do not write output function now, try to make it as simple as possible.
    int irow_index, icol_index;
    int start_index;
    int m,i,j,k;
    int a;
    int b;
    double var;

    double start_clock;
    double end_clock;

    double Magnitude;
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    int steps;
    double detector_tprint = tprint; // set detector_tprint == tprint specified in input file.
    int output_step= int(detector_tprint/delt); //Output every output_step.

    int scale_of_output;
    double output_all_state_time_unit ;
    int output_step_for_all_state;
    double minimum_output_all_state_time_unit = pow(10,-2) ;

    double * start_time = new double [s.tlnum];
    double * final_time = new double [s.tlnum];
    for(i=0;i<s.tlnum;i++){
        start_time[i]=0;
    }
    int state_number_for_compputing_OTOC_for_xp = 0;
    ofstream four_point_correlation_output;
    ofstream four_point_correlation_average_output;
    ofstream four_point_correlation_variance_output;
    ofstream four_point_correlation_every_state_output;
    ofstream Detector_precoup_output;
    ofstream Detector_precoup_mode_quanta;
    ofstream Stability_matrix_output;
    ofstream Detector_precoup_all_state_output;

    ofstream n_offdiag_total_output; // output <j| n(t) |l> in our system

    ofstream IPR_output;
    ofstream IPR_for_all_state_output;

    ofstream overlap_with_initial_state_output;

    ofstream another_form_of_OTOC_output;

    ofstream Lyapunov_spectrum_for_xp_output;

    ofstream Every_states_contribution_to_OTOC_xp;

    ofstream TOC_output ;

    ofstream regularized_Thermal_OTOC_output;
    // -----------Open 4_point_correlation_output ofstream -----------------------------
    if(my_id==0){
        if(Detector_Continue_Simulation){
            four_point_correlation_output.open(path + "4 point correlation.txt", ofstream::app);
            four_point_correlation_average_output.open(path+ "4 point correlation average over states.txt",ofstream::app);
            four_point_correlation_variance_output.open(path+"4 point correlation function variance.txt",ofstream::app);
            four_point_correlation_every_state_output.open(path+"4 point correlation function all states.txt",ofstream::app);
            Detector_precoup_output.open(path + "detector_precoupling.txt", ofstream::app);
            Detector_precoup_all_state_output.open(path + "detector_precoupling_all_state.txt",ofstream::app);
            Detector_precoup_mode_quanta.open(path + "detector_precoup_mode_quanta.txt", ofstream::app);
            Stability_matrix_output.open(path+"Stability_Matrix.txt",ofstream::app);

            n_offdiag_total_output.open(path+"n_offdiag.txt",ofstream::app);
            IPR_output.open(path+"IPR.txt",ofstream::app);
            IPR_for_all_state_output.open(path+"IPR_all_state.txt",ofstream::app);
            overlap_with_initial_state_output.open(path+"other_state_overlap_with_initial_state.txt",ofstream::app);

            another_form_of_OTOC_output.open(path+ "another_OTOC.txt", ofstream::app);
            Lyapunov_spectrum_for_xp_output.open(path+ "Lyapunov_spectrum_for_xp.txt",ofstream::app);
            Every_states_contribution_to_OTOC_xp.open(path + "OTOC_xp_for_every_state.txt", ofstream::app);

            TOC_output.open(path + "TOC_factorization.txt", ofstream::app);

            regularized_Thermal_OTOC_output.open(path + "regularized_Thermal_OTOC.txt" , ofstream::app);
        }
        else {
            four_point_correlation_output.open(path + "4 point correlation.txt");
            four_point_correlation_average_output.open(path+ "4 point correlation average over states.txt");
            four_point_correlation_variance_output.open(path+"4 point correlation function variance.txt");
            four_point_correlation_every_state_output.open(path+"4 point correlation function all states.txt");
            Detector_precoup_output.open(path + "detector_precoupling.txt");
            Detector_precoup_all_state_output.open(path + "detector_precoupling_all_state.txt");
            Detector_precoup_mode_quanta.open(path + "detector_precoup_mode_quanta.txt");
            Stability_matrix_output.open(path+"Stability_Matrix.txt");

            n_offdiag_total_output.open(path+"n_offdiag.txt");
            IPR_output.open(path+"IPR.txt");
            IPR_for_all_state_output.open(path+"IPR_all_state.txt");
            overlap_with_initial_state_output.open(path+"other_state_overlap_with_initial_state.txt");

            another_form_of_OTOC_output.open(path+ "another_OTOC.txt");
            Lyapunov_spectrum_for_xp_output.open(path+"Lyapunov_spectrum_for_xp.txt");
            Every_states_contribution_to_OTOC_xp.open(path+"OTOC_xp_for_every_state.txt");
            TOC_output.open(path + "TOC_factorization.txt");

            regularized_Thermal_OTOC_output.open(path + "regularized_Thermal_OTOC.txt" );
        }
    }

    // ------------- prepare variable computing for Lyapunovian for xp --------------
    d.prepare_computing_Lyapunovian_for_xp();

    // -----------------------------------------------------------------------------------------------
    // prepare sendbuffer and recv_buffer and corresponding index.
    d.prepare_evolution();

    // ---------- prepare computing Boltzmann weighted state ------------
    d.prepare_compute_Boltzmann_factor_use_Chebyshev_polynomial(Boltzmann_beta/4, log);

    // ---- prepare computing state after ladder operator a_{j} operation ------
    d.prepare_ladder_operation();

    // -------------Load detector state from save data if we want to continue simulation of detector.------------------
    if(Detector_Continue_Simulation){
        d.load_detector_state_MPI(path,start_time,log,initial_state_choice);
    }
    else{
        d.initialize_detector_state_MPI(log); // initialize detector lower bright state
    }

    // generate Boltzmann weighted wave function for basis set.
    // also compute result : decorate Boltzmann weighted wave function with ladder operator and compute <{m} | a_{j} e^{-\beta H} | {n} >
    d.Boltzmann_factor_decorated_basis_set_and_with_ladder_operator( );

    // compute normalization factor: <Haar | e^{-\beta H } | Haar>
    d.compute_normalization_factor_for_Boltzmann_weighted_factor() ;

    // allocate space for Haar random variable calculation
    d.allocate_space_for_Haar_state_calculation();

    int initial_state_index_in_total_dmatrix;
    initial_state_index_in_total_dmatrix = d.initial_state_index[0]
            + d.total_dmat_size[0]/num_proc * d.initial_state_pc_id[0] ;





    // ---------- allocate space for mode quanta -----------------
    double * total_mode_quanta;
    total_mode_quanta= new double [d.nmodes[0]];


    double * mode_quanta;
    mode_quanta= new double [d.nmodes[0]];


    // -------------- Allocate space for <a| |[n_{i}(t),n_{i}(0)]|^{2} |a> -------------
    int state_for_average_size = d.states_for_4_point_correlation_average.size();
    int nearby_state_index_size = d.nearby_state_index.size();
    double * four_point_correlation_function_at_initial_state = new double [d.nmodes[0]];
    double * four_point_correlation_function_average_over_states = new double [d.nmodes[0]];
    double * four_point_correlation_function_variance_over_states = new double [d.nmodes[0]];
    double ** four_point_correlation_function_for_each_states = new double * [d.nmodes[0]]; // this is average over bunch of states
    for(i=0;i<d.nmodes[0];i++){
        four_point_correlation_function_for_each_states[i] = new double [state_for_average_size];
    }

    // ----------------- variable for IPR --------------
    double IPR;

    double * inverse_IPR_in_one_process = new double [nearby_state_index_size]  ;
    double * inverse_IPR_all = new double [nearby_state_index_size];
    double * IPR_all = new double [nearby_state_index_size];

    // ------------ Other states in nearby state index list's overlap with initial state ---------------------
    double * other_state_overlap_with_initial_state = new double [nearby_state_index_size];

    //---------- Allocate space for <a| n_{i}(t) |b>  size: nmode * total_dmat_size[0] -------------------------------------------
    // each row is one n_{i}.  each column is one site |b>
    complex<double> ** n_offdiag_total;
    double ** n_offdiag_total_real;
    double ** n_offdiag_total_imag;
    n_offdiag_total = new complex<double> * [d.nmodes[0]];
    n_offdiag_total_real = new double * [d.nmodes[0]];
    n_offdiag_total_imag = new double * [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag_total[i] = new complex <double> [nearby_state_index_size];
        n_offdiag_total_real[i] = new double [nearby_state_index_size];
        n_offdiag_total_imag[i] = new double [nearby_state_index_size];
    }
    complex<double> ** n_offdiag;
    double ** n_offdiag_real;
    double ** n_offdiag_imag;
    n_offdiag = new complex <double> * [d.nmodes[0]];
    n_offdiag_real = new double * [d.nmodes[0]];
    n_offdiag_imag = new double * [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag [i] = new complex <double> [nearby_state_index_size];
        n_offdiag_real[i] = new double [nearby_state_index_size];
        n_offdiag_imag[i] = new double [nearby_state_index_size];
    }


    //  Another form of OTOC
    double * another_OTOC = new double [d.nmodes[0]];
    complex<double> ** n_offdiag_total_one_mode_quanta_below;  // compute <m-| n_{i}(t) | l->  states one mode below in corresponding mode index
    double ** n_offdiag_total_one_mode_quanta_below_real;
    double ** n_offdiag_total_one_mode_quanta_below_imag;
    n_offdiag_total_one_mode_quanta_below = new complex<double> * [d.nmodes[0]];
    n_offdiag_total_one_mode_quanta_below_real = new double * [d.nmodes[0]];
    n_offdiag_total_one_mode_quanta_below_imag = new double * [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag_total_one_mode_quanta_below[i] = new complex<double> [nearby_state_index_size];
        n_offdiag_total_one_mode_quanta_below_real[i] = new double [nearby_state_index_size];
        n_offdiag_total_one_mode_quanta_below_imag[i] = new double [nearby_state_index_size];
    }

    complex<double> ** n_offdiag_one_mode_quanta_below;  // compute <m-| n_{i}(t) | l->  states one mode below in corresponding mode index
    double ** n_offdiag_one_mode_quanta_below_real;
    double ** n_offdiag_one_mode_quanta_below_imag;
    n_offdiag_one_mode_quanta_below = new complex<double> * [d.nmodes[0]];
    n_offdiag_one_mode_quanta_below_real = new double * [d.nmodes[0]];
    n_offdiag_one_mode_quanta_below_imag = new double * [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag_one_mode_quanta_below[i] = new complex<double> [nearby_state_index_size];
        n_offdiag_one_mode_quanta_below_real[i] = new double [nearby_state_index_size];
        n_offdiag_one_mode_quanta_below_imag[i] = new double [nearby_state_index_size];
    }


    complex<double> * n_offdiag_element;
    n_offdiag_element = new complex <double> [d.nmodes[0]];
    double * n_offdiag_element_real = new double [d.nmodes[0]];
    double * n_offdiag_element_imag = new double [d.nmodes[0]];
    complex<double> * n_offdiag_element_all_pc;
    n_offdiag_element_all_pc = new complex <double> [d.nmodes[0]];
    double * n_offdiag_element_real_all_pc = new double [d.nmodes[0]];
    double * n_offdiag_element_imag_all_pc = new double [d.nmodes[0]];

    // ---------- Allocate space for <a(states)|n_{i}(t)|b> for computing average over states  ---------------
    // -------------------  size: nmodes * state_for_average_num * total_dmat_size[0] ---------------------------------------------
    complex<double> *** n_offdiag_total_for_states_ensemble = new complex<double> ** [d.nmodes[0]];
    double *** n_offdiag_total_for_states_ensemble_real = new double ** [d.nmodes[0]];
    double *** n_offdiag_total_for_states_ensemble_imag = new double ** [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag_total_for_states_ensemble[i] = new complex <double> * [state_for_average_size];
        n_offdiag_total_for_states_ensemble_real [i] = new double * [state_for_average_size];
        n_offdiag_total_for_states_ensemble_imag [i] = new double * [state_for_average_size];
        for(j=0;j<state_for_average_size;j++){
            n_offdiag_total_for_states_ensemble[i][j] = new complex <double> [nearby_state_index_size];
            n_offdiag_total_for_states_ensemble_real[i][j] = new double [nearby_state_index_size];
            n_offdiag_total_for_states_ensemble_imag[i][j] = new double [nearby_state_index_size];
        }
    }

    complex<double> *** n_offdiag_for_states_ensemble = new complex <double> ** [d.nmodes[0]];
    double *** n_offdiag_for_states_ensemble_real = new double ** [d.nmodes[0]];
    double *** n_offdiag_for_states_ensemble_imag = new double ** [d.nmodes[0]];
    for(i=0;i<d.nmodes[0];i++){
        n_offdiag_for_states_ensemble[i] = new complex <double> * [state_for_average_size];
        n_offdiag_for_states_ensemble_real [i] = new double * [state_for_average_size];
        n_offdiag_for_states_ensemble_imag [i] = new double * [state_for_average_size];
        for(j=0;j<state_for_average_size;j++){
            n_offdiag_for_states_ensemble[i][j] = new complex<double> [nearby_state_index_size];
            n_offdiag_for_states_ensemble_real[i][j] = new double [nearby_state_index_size];
            n_offdiag_for_states_ensemble_imag[i][j] = new double [nearby_state_index_size];
        }
    }

    // ----------- Lyapunov spectrum for x,p operator: [a_{i}(t),a_{j}] or [a_{i}^{+}(t), a_{j}] or [a_{i}(t),a_{j}^{+}] or [a_{i}^{+}(t), a_{j}^{+}(t) ]-----------------

    complex<double> ** Lyapunov_spectrum_for_xp;
    Lyapunov_spectrum_for_xp = new complex<double> * [2 * d.nmodes[0]];
    for(i=0;i< 2*d.nmodes[0];i++){
        Lyapunov_spectrum_for_xp[i] = new complex<double> [2 * d.nmodes[0]];
    }

    complex<double> ** Lyapunov_spectrum_for_xp_from_single_state ;
    Lyapunov_spectrum_for_xp_from_single_state = new complex<double> * [2* d.nmodes[0]];
    for(i=0;i<2*d.nmodes[0];i++){
        Lyapunov_spectrum_for_xp_from_single_state[i] = new complex<double> [2 * d.nmodes[0]];
    }

    complex<double> ** Matrix_M; // L is Lyapunov spectrum,<l| L_{ij} |l> = \sum_{m}  [ \sum_{k} (M_{ki}^{ml})^{*} ( M_{kj}^{ml} )  ]
    // Here M_{kj}^{ml} = <m| M_{kj} | l>: M_{kj} = [c_{k}(t), c_{j}] where c_{k} = a_{k} if (0<=k<=N-1), c_{k} = (a_{k-N})^{+} if N<=k<= 2N-1
    Matrix_M = new complex<double> * [2 * d.nmodes[0]];
    for(i=0;i< 2*d.nmodes[0];i++){
        Matrix_M[i] = new complex<double> [ 2 * d.nmodes[0]];
    }

    vector<vector<double>> * xd_for_xp = new vector<vector<double>> [1 + 2*d.nmodes[0]]; // assum initial state l. Then we have [ l , l_{k}^{-} (k=1, cdots ,N) , l_{k}^{+}] (k=1, cdots, N)
    vector<vector<double>> * yd_for_xp = new vector<vector<double>> [1 + 2*d.nmodes[0]];
    // size [ 1 + 2*d.nmodes[0] , 2* d.nmodes[0], dmatsize[0] ]
    vector <double> v1 (d.dmatsize[0],0);
    vector<vector<double>> v2;
    for(j=0;j <2*d.nmodes[0];j++){ // this is for <j^{-} | l> or  <j^{+}| l>
        v2.push_back(v1);
    }
    for(i=0;i<1+2*d.nmodes[0];i++){ // first index i is for trajectory
        xd_for_xp[i] = v2;
        yd_for_xp[i] = v2;
    }

    // ----------------- Variable for Time ordered correlation function TOC -------------------
    double * Two_point_func1 = new double [2 * d.nmodes[0]];
    double * Two_point_func2 = new double [2* d.nmodes[0] ];
    double ** TOC_per_state = new double * [ 2* d.nmodes[0] ];
    double ** TOC = new double * [2* d.nmodes[0]];
    for(i=0;i<2*d.nmodes[0];i++){
        TOC_per_state[i] = new double [2 * d.nmodes[0] ];
        TOC[i] = new double [2 * d.nmodes[0] ];
    }


    //------ Allocate space for sparse version of xd_for_xp and yd_for_xp
    vector<vector<double>> * xd_for_xp_sparsify = new vector<vector<double>> [1 + 2 * d.nmodes[0]];
    vector<vector<double>> * yd_for_xp_sparsify = new vector<vector<double>> [1 + 2* d.nmodes[0] ];
    vector<vector<int>> * index_for_xp_sparsify = new vector<vector<int>> [1+ 2* d.nmodes[0] ];
// size [ 1 + 2*d.nmodes[0] , 2* d.nmodes[0], dmatsize[0] ]


    // ---------- Allocate space for stability matrix L:  ---------------------
    // L = sum_{b} (sum_{k}| <a |n_{k}(t)|b> |^{2}* (n_{i}^{b} - n_{i}^{a})* (n_{j}^{b} - n_{j}^{b}) ) here k,i,j is dof, a,b is state
    double ** Stability_Matrix;
    Stability_Matrix = new double * [d.nmodes[0]] ;
    for(i=0;i<d.nmodes[0];i++){
        Stability_Matrix[i] = new double [d.nmodes[0]];
    }


    // ------- For output for our state <a(t)|a> trajectory (check ergodicity)
    double special_state_x;
    double special_state_y;
    double survival_prob;


    // -----  compute state's energy and shift it before doing simulation -------------
    vector<complex<double>>  H_phi;
    H_phi.resize(d.dmatsize[0]);
    for(i=0;i<d.dmatsize[0];i++){
        H_phi[i] = 0;
    }
    double de;
    double de_all;
    for(i=0;i<d.dmatnum[0];i++){
        irow_index = d.local_dirow[0][i];
        icol_index = d.local_dicol[0][i]; // compute to point to colindex in
        H_phi[irow_index] = H_phi[irow_index] + d.dmat[0][i] * complex(d.xd[d.initial_state_index_in_nearby_state_index_list][icol_index],
                                                                       d.yd[d.initial_state_index_in_nearby_state_index_list][icol_index]);
    }
    de=0;
    for(i=0;i<d.dmatsize[0];i++){
        de= de+ real(H_phi[i] * complex(d.xd[d.initial_state_index_in_nearby_state_index_list][i],
                                        -d.yd[d.initial_state_index_in_nearby_state_index_list][i]));
    }
    MPI_Allreduce(&de,&de_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    // shift Hamiltonian by detector energy de_all:
    if(my_id == 0){
        cout <<"Shift Hamiltonian by energy of system  "<< de_all << endl;
        log << "Shift Hamiltonian by energy of system  " << de_all <<endl;
    }
    for(i=0;i<d.dmatsize[0];i++){
        d.dmat[0][i] = d.dmat[0][i] - de_all;
    }


    final_time[0] = 0;
    if(d.proptime[0]>0){
        if(my_id==0){
            log <<" Begin Propogation in Molecule Phase space "<<endl;
        }
        t=0;
        steps= d.proptime[0]/delt+1;
        if(Detector_Continue_Simulation){
            // update simulation time and t for detector to begin simulation.
            t=start_time[0];
            steps= (d.proptime[0]-start_time[0])/delt+2;
        }
        // output information to output file before we go into loop
        if(!Detector_Continue_Simulation) {
            if (my_id == 0) {
                four_point_correlation_output << "4- point correlation function for molecule " << endl;
                four_point_correlation_output << "total time: " << d.proptime[0]  << " "
                                              << delt * output_step << endl;
//                four_point_correlation_average_output << "4 - point correlation function for molecule average over states" <<endl;
//                four_point_correlation_average_output << "total time: " << d.proptime[0]  << " "
//                                              << delt * output_step << endl;
//                four_point_correlation_variance_output << "4 - point correlation function for molecule variance over states" <<endl;
//                four_point_correlation_variance_output << "total time: " << d.proptime[0]  << " "
//                                                      << delt * output_step << endl;
//
//                four_point_correlation_every_state_output << d.initial_state_index_in_states_for_4_point_correlation_list<<endl;
//                four_point_correlation_every_state_output <<   state_for_average_size <<endl;

                // output state mode number information for all states we compute for 4 point correlation function
                for (i=0;i<state_for_average_size;i++){
                    for(j=0;j<d.nmodes[0];j++){
                        four_point_correlation_every_state_output << d.dv_all[0][d.states_for_4_point_correlation_average[i]][j] <<" ";
                    }
                    four_point_correlation_every_state_output << endl;
                }
                four_point_correlation_every_state_output << "total time: " << d.proptime[0]  << " "
                                                       << delt * output_step << endl;

                // output state information for which we compute overlap and n_offdiag for b
                n_offdiag_total_output << d.initial_state_index_in_nearby_state_index_list << endl;
                n_offdiag_total_output << nearby_state_index_size << endl;
                for(i=0;i<nearby_state_index_size;i++){
                    for(j=0;j< d.nmodes[0]; j++){
                        n_offdiag_total_output << d.dv_all[0][d.nearby_state_index[i]][j] <<" ";
                    }
                    n_offdiag_total_output << endl;
                }

                Stability_matrix_output << d.nmodes[0] <<endl;

                for(i=0;i<d.nmodes[0];i++){
                    Detector_precoup_mode_quanta << d.mfreq[0][i] <<" ";
                }
                Detector_precoup_mode_quanta<<endl;

//                Detector_precoup_all_state_output << d.proptime[0] << endl;
//                Detector_precoup_all_state_output << initial_state_index_in_total_dmatrix << endl;
//                Detector_precoup_all_state_output << d.total_dmat_size[0] << endl;
//                for(i=0;i<d.total_dmat_size[0];i++){
//                    for(j=0;j<d.nmodes[0];j++){
//                        Detector_precoup_all_state_output << d.dv_all[0][i][j] <<" ";
//                    }
//                        Detector_precoup_all_state_output << endl;
//                }

//                IPR_for_all_state_output << nearby_state_index_size << endl;
//                for(i=0;i<nearby_state_index_size;i++){
//                    for(j=0;j<d.nmodes[0];j++){
//                        IPR_for_all_state_output << d.dv_all[0][d.nearby_state_index[i]][j] <<"  ";
//                    }
//                    IPR_for_all_state_output << endl;
//                }

//                overlap_with_initial_state_output<< nearby_state_index_size << endl;
//                for(i=0;i<nearby_state_index_size;i++){
//                    for(j=0;j<d.nmodes[0];j++){
//                        overlap_with_initial_state_output << d.dv_all[0][d.nearby_state_index[i]][j] <<"  ";
//                    }
//                    overlap_with_initial_state_output << endl;
//                }

                Lyapunov_spectrum_for_xp_output << d.nmodes[0] << endl;

                TOC_output << d.nmodes[0] << endl;

                regularized_Thermal_OTOC_output << d.nmodes[0] << endl;

//                Every_states_contribution_to_OTOC_xp << d.nmodes[0] << endl;
//                for(m=0;m<nearby_state_index_size;m++){
//                    if(d.bool_neighbor_state_all_in_nearby_state_index[m]){
//                        state_number_for_compputing_OTOC_for_xp = state_number_for_compputing_OTOC_for_xp + 1;
//                    }
//                }
//                Every_states_contribution_to_OTOC_xp << state_number_for_compputing_OTOC_for_xp << endl;
//                for(m=0;m<nearby_state_index_size;m++){
//                    if(d.bool_neighbor_state_all_in_nearby_state_index[m]){
//                        for(i=0;i<d.nmodes[0];i++){
//                            Every_states_contribution_to_OTOC_xp << d.dv_all[0][d.nearby_state_index[m]][i] <<" ";
//                        }
//                        Every_states_contribution_to_OTOC_xp << endl;
//                    }
//                }

            }
        }

        if(Detector_Continue_Simulation){
            start_index=1;
        }
        else{
            start_index = 0;
        }
        // Do simulation in loop
        for(k=start_index;k<steps;k++){
            //-------------------- output result ----------------------------
            if(t!=0){
                scale_of_output =  floor(std::log(t/delt)/std::log(10));  // 0 for 10^{-4} : 1 for 10^{-1}
            }
            else{
                scale_of_output = 0;
            }
            output_all_state_time_unit = min(delt * pow(10,scale_of_output), minimum_output_all_state_time_unit);
            output_step_for_all_state = int(output_all_state_time_unit / delt);

//            // ----- output state real and imaginary part for all state in simulation ----------------------
//            if(k%output_step_for_all_state==0){
//                MPI_Gatherv(&d.xd[d.initial_state_index_in_nearby_state_index_list][0],d.dmatsize[0],MPI_DOUBLE,
//                            &d.xd_all[d.initial_state_index_in_nearby_state_index_list][0],d.dmatsize_each_process[0],
//                            d.dmatsize_offset_each_process[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
//                MPI_Gatherv(&d.yd[d.initial_state_index_in_nearby_state_index_list][0], d.dmatsize[0], MPI_DOUBLE,
//                            &d.yd_all[d.initial_state_index_in_nearby_state_index_list][0], d.dmatsize_each_process[0],
//                            d.dmatsize_offset_each_process[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
//                if(my_id == 0){
//                    Detector_precoup_all_state_output << t << endl;
//                    for(i=0;i<d.total_dmat_size[0];i++){
//                        Detector_precoup_all_state_output << d.xd_all[d.initial_state_index_in_nearby_state_index_list][i] <<" ";
//                    }
//                    Detector_precoup_all_state_output << endl;
//                    for(i=0;i<d.total_dmat_size[0];i++){
//                        Detector_precoup_all_state_output << d.yd_all[d.initial_state_index_in_nearby_state_index_list][i] <<" ";
//                    }
//                    Detector_precoup_all_state_output << endl;
//                }
//            }
//            // --------- output IPR for all state --------------------------------------
//            if(k%output_step_for_all_state == 0){
//                for(i=0;i<nearby_state_index_size;i++){
//                    inverse_IPR_in_one_process[i] = 0;
//                    for(j=0;j<d.dmatsize[0];j++){
//                        inverse_IPR_in_one_process[i] = inverse_IPR_in_one_process[i] + pow(pow(d.xd[i][j],2) + pow(d.yd[i][j],2),2);
//                    }
//                }
//
//                MPI_Reduce(&inverse_IPR_in_one_process[0],&inverse_IPR_all[0],nearby_state_index_size,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
//
//                if(my_id == 0){
//                    for(i=0;i<nearby_state_index_size;i++){
//                        IPR_all[i] = 1/ inverse_IPR_all[i];
//                    }
//
//                    // output result
//                    IPR_for_all_state_output << t << endl;
//                    for(i=0;i<nearby_state_index_size;i++){
//                        IPR_for_all_state_output << IPR_all[i] <<" ";
//                    }
//                    IPR_for_all_state_output << endl;
//                }
//
//            }


            if(k % output_step ==0) {
                // ---------------------------------output 4-point correlation function -------------------------------------------------
               compute_4_point_corre_for_single_state(nearby_state_index_size,n_offdiag_element,n_offdiag,n_offdiag_real,n_offdiag_imag,
                                                      n_offdiag_total,n_offdiag_total_real,n_offdiag_total_imag,initial_state_index_in_total_dmatrix,
                                                      four_point_correlation_function_at_initial_state);
                if(my_id == 0){
                    four_point_correlation_output << "Time:   " << t << endl;
                    for(i=0;i<d.nmodes[0];i++){
                        four_point_correlation_output << four_point_correlation_function_at_initial_state[i] << " ";
                    }
                    four_point_correlation_output<<endl;
                }
                // ------------- output overlap with initial state and n_offdiag_total for all state |b> we choose:
//                if(my_id == 0){
//                    n_offdiag_total_output <<"Time:  "<<t << endl;
//                    for(i=0; i< nearby_state_index_size; i++){
//                        for(j=0;j<d.nmodes[0];j++){
//                            n_offdiag_total_output << std::norm(n_offdiag_total[j][i]) <<"  ";
//                        }
//                        n_offdiag_total_output << endl;
//                    }
//                }

                // -------- output another form of OTOC ------------------------------
//                compute_another_form_of_OTOC(nearby_state_index_size,n_offdiag_element,n_offdiag,n_offdiag_real,n_offdiag_imag,
//                                             n_offdiag_total,n_offdiag_total_real,n_offdiag_total_imag,
//                                             n_offdiag_one_mode_quanta_below,n_offdiag_one_mode_quanta_below_real,n_offdiag_one_mode_quanta_below_imag,
//                                             n_offdiag_total_one_mode_quanta_below,n_offdiag_total_one_mode_quanta_below_real,n_offdiag_total_one_mode_quanta_below_imag,
//                                             initial_state_index_in_total_dmatrix, another_OTOC);
//                if(my_id == 0){
//                    another_form_of_OTOC_output << "Time:  " << t <<endl;
//                    for(i=0;i<d.nmodes[0];i++){
//                        another_form_of_OTOC_output << another_OTOC[i] << "  ";
//                    }
//                    another_form_of_OTOC_output << endl;
//                }

//                 ----------- output Lyapunovian spectrum for xp ---------------------------------------
//                d.compute_Lyapunov_spectrum_for_xp(Lyapunov_spectrum_for_xp,Lyapunov_spectrum_for_xp_from_single_state,
//                                                   Matrix_M,
//                                                   xd_for_xp,yd_for_xp,
//                                                   xd_for_xp_sparsify, yd_for_xp_sparsify, index_for_xp_sparsify,
//                                                   t,Every_states_contribution_to_OTOC_xp);
//                if(my_id == 0){
//                    Lyapunov_spectrum_for_xp_output << t <<endl;
//                    for(i=0;i< 2 * d.nmodes[0];i++){
//                        for(j=0;j<2*d.nmodes[0];j++){
//                            Lyapunov_spectrum_for_xp_output << real(Lyapunov_spectrum_for_xp[i][j]) << " ";
//                        }
//                    }
//                    Lyapunov_spectrum_for_xp_output << endl;
//
//                    for(i=0;i<2*d.nmodes[0];i++){
//                        for(j=0;j<2*d.nmodes[0];j++){
//                            Lyapunov_spectrum_for_xp_output << imag(Lyapunov_spectrum_for_xp[i][j]) <<" ";
//                        }
//                    }
//                    Lyapunov_spectrum_for_xp_output << endl;
//
//                }

                // --------------- output information for Time ordered correlation function --------
//                d.output_TOC_factorization(Two_point_func1, Two_point_func2, TOC_per_state, TOC, xd_for_xp, yd_for_xp,
//                                           xd_for_xp_sparsify, yd_for_xp_sparsify, index_for_xp_sparsify, t, TOC_output);

                // -------------- compute regularized thermal OTOC --------------------
                start_clock = clock();
                d.compute_regularized_thermal_OTOC_Lyapunov_spectrum( );
                end_clock = clock();
                if(my_id == 0){
                    cout << " Finish computing one step of regularized thermal OTOC , t = " << (end_clock - start_clock)/CLOCKS_PER_SEC << endl;
                    regularized_Thermal_OTOC_output << t <<endl;
                    for(i=0;i< 2 * d.nmodes[0];i++){
                        for(j=0;j<2*d.nmodes[0];j++){
                            regularized_Thermal_OTOC_output << real(d.regularized_thermal_Lyapunov_spectrum[i][j]) << " ";
                        }
                    }
                    regularized_Thermal_OTOC_output << endl;

                    for(i=0;i<2*d.nmodes[0];i++){
                        for(j=0;j<2*d.nmodes[0];j++){
                            regularized_Thermal_OTOC_output << imag(d.regularized_thermal_Lyapunov_spectrum[i][j]) <<" ";
                        }
                    }
                        regularized_Thermal_OTOC_output << endl;

                }


                // ---------- output 4-point correlation function average over states and variance -------------------
//                compute_4_point_corre_for_multiple_states(state_for_average_size,nearby_state_index_size,n_offdiag_element,
//                                                          n_offdiag_for_states_ensemble,n_offdiag_for_states_ensemble_real,
//                                                          n_offdiag_for_states_ensemble_imag,n_offdiag_total_for_states_ensemble,
//                                                          n_offdiag_total_for_states_ensemble_real,n_offdiag_total_for_states_ensemble_imag,
//                                                          four_point_correlation_function_average_over_states, four_point_correlation_function_variance_over_states,
//                                                          four_point_correlation_function_for_each_states);
//
//                if(my_id == 0){
//                    four_point_correlation_average_output << "Time:   " << t << endl;
//                    for(i=0;i<d.nmodes[0];i++){
//                        four_point_correlation_average_output<< four_point_correlation_function_average_over_states[i] << " ";
//                    }
//                    four_point_correlation_average_output<<endl;
//
//                    // for variance.
//                    four_point_correlation_variance_output <<"Time:    " << t << endl;
//                    for(i=0;i<d.nmodes[0];i++){
//                        four_point_correlation_variance_output << four_point_correlation_function_variance_over_states[i]<< " ";
//                    }
//                    four_point_correlation_variance_output<<endl;
//
//                    // output all states
//                    four_point_correlation_every_state_output << "Time:   " <<t << endl;
//                    for(a=0;a<state_for_average_size;a++){
//                        for(i=0;i<d.nmodes[0];i++){
//                            four_point_correlation_every_state_output << four_point_correlation_function_for_each_states[i][a] <<" ";
//                        }
//                        four_point_correlation_every_state_output << endl;
//                    }
//                }

                // --------- output Stability_Matrix  -----------------------------------
//                compute_Stability_Matrix(Stability_Matrix,nearby_state_index_size,n_offdiag_element,n_offdiag_element_real,n_offdiag_element_imag,
//                                         n_offdiag_element_all_pc,n_offdiag_element_real_all_pc,n_offdiag_element_imag_all_pc,initial_state_index_in_total_dmatrix);
//
//                if(my_id == 0){
//                    Stability_matrix_output << t << endl;
//                    for(i=0;i<d.nmodes[0];i++){
//                        for(j=0;j<d.nmodes[0];j++){
//                            Stability_matrix_output << Stability_Matrix[i][j] <<" ";
//                        }
//                    }
//                    Stability_matrix_output << endl;
//                }
                // ------------- output state <a(t)|a>  ----------------------------------------
                if (my_id == d.initial_state_pc_id[0]){
                    special_state_x = d.xd[d.initial_state_index_in_nearby_state_index_list][d.initial_state_index[0]];
                    special_state_y = d.yd[d.initial_state_index_in_nearby_state_index_list][d.initial_state_index[0]];
                }
                MPI_Bcast(&special_state_x,1,MPI_DOUBLE,d.initial_state_pc_id[0],MPI_COMM_WORLD);
                MPI_Bcast(&special_state_y,1,MPI_DOUBLE,d.initial_state_pc_id[0],MPI_COMM_WORLD);
                survival_prob = pow(special_state_x,2) + pow(special_state_y,2);
                if(my_id == 0){
                    Detector_precoup_output << "Time:   " << t << endl;
                    Detector_precoup_output << survival_prob << endl;
                }

                // --------- output other state's overlap with initial state |a> ----------------------
//                if(my_id == d.initial_state_pc_id[0]){
//                    for(i=0;i<nearby_state_index_size;i++){
//                        other_state_overlap_with_initial_state[i] = pow(d.xd[i][d.initial_state_index[0]] ,2)
//                                                                    + pow(d.yd[i][d.initial_state_index[0]],2) ;
//                    }
//                }
//                if(d.initial_state_pc_id[0] != 0){
//                    // pass data to process with id 0 if data is not initially there
//                    if(my_id == d.initial_state_pc_id[0]){
//                        MPI_Send(&other_state_overlap_with_initial_state[0],nearby_state_index_size,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
//                    }
//                    if(my_id == 0){
//                        MPI_Recv(&other_state_overlap_with_initial_state[0],nearby_state_index_size,MPI_DOUBLE,d.initial_state_pc_id[0],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
//                    }
//                }
//
//                if(my_id == 0){
//                    overlap_with_initial_state_output << t << endl;
//                    for(i=0;i<nearby_state_index_size;i++){
//                        overlap_with_initial_state_output << other_state_overlap_with_initial_state[i] <<" ";
//                    }
//                    overlap_with_initial_state_output << endl;
//                }

                // ----  output IPR (inverse participation ratio) ----------------------------
                MPI_Gatherv(&d.xd[d.initial_state_index_in_nearby_state_index_list][0],d.dmatsize[0],MPI_DOUBLE,
                &d.xd_all[d.initial_state_index_in_nearby_state_index_list][0],d.dmatsize_each_process[0],d.dmatsize_offset_each_process[0],MPI_DOUBLE,0,MPI_COMM_WORLD);

                MPI_Gatherv(&d.yd[d.initial_state_index_in_nearby_state_index_list][0], d.dmatsize[0], MPI_DOUBLE,
                                             &d.yd_all[d.initial_state_index_in_nearby_state_index_list][0], d.dmatsize_each_process[0],
                                             d.dmatsize_offset_each_process[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
                if(my_id == 0){
                    IPR = 0;
                    IPR_output << t << endl;
                    for(i=0;i<d.total_dmat_size[0];i++){
                        Magnitude = pow(d.xd_all[d.initial_state_index_in_nearby_state_index_list][i] ,2) +
                                pow(d.yd_all[d.initial_state_index_in_nearby_state_index_list][i],2);
                        IPR = IPR + pow(Magnitude,2);
                    }
                    IPR = 1/ IPR;
                    IPR_output << IPR << endl;
                }

                // ----------- output mode quanta -----------------------------------
                for (j = 0; j < d.nmodes[0]; j++) {
                    mode_quanta[j] = 0;
                }
                for (i = 0; i < d.dmatsize[0]; i++) {
                    for (j = 0; j < d.nmodes[0]; j++) {
                        mode_quanta[j] =
                                mode_quanta[j] + (pow(d.xd[d.initial_state_index_in_nearby_state_index_list][i], 2) +
                                pow(d.yd[d.initial_state_index_in_nearby_state_index_list][i], 2)) * d.dv[0][i][j];
                    }
                }
                MPI_Reduce(&mode_quanta[0], &total_mode_quanta[0], d.nmodes[0], MPI_DOUBLE, MPI_SUM, 0,
                           MPI_COMM_WORLD);
                if (my_id == 0) {
                    Detector_precoup_mode_quanta << "Time:   " << t << endl;
                    for (j = 0; j < d.nmodes[0]; j++) {
                        Detector_precoup_mode_quanta << total_mode_quanta[j] << " ";
                    }
                    Detector_precoup_mode_quanta << endl;
                }

            }
            t= t+ delt;
            start_clock = clock();
            d.SUR_onestep_MPI(cf);
            end_clock = clock();
            if(my_id == 0 and t == delt ){
                cout << "Evolve one step t : " << (end_clock - start_clock)/CLOCKS_PER_SEC << endl;
            }
        }
        final_time[0] = t;
    }
    if(save_state){
        d.save_detector_state_MPI(path,final_time,log,initial_state_choice);
    }
    d.replace_4_point_corr_second_line(detector_tprint);

    if(my_id==0){
        cout<<"Detector_pre_coupling simulation finished"<<endl;
        four_point_correlation_output.close();
        four_point_correlation_average_output.close();
        four_point_correlation_variance_output.close();
        four_point_correlation_every_state_output.close();
        Detector_precoup_output.close();
        Detector_precoup_all_state_output.close();
        Detector_precoup_mode_quanta.close();

        n_offdiag_total_output.close();
        IPR_output.close();
        IPR_for_all_state_output.close();
        overlap_with_initial_state_output.close();
        another_form_of_OTOC_output.close();
        Lyapunov_spectrum_for_xp_output.close();
        Every_states_contribution_to_OTOC_xp.close();
        TOC_output.close();

        regularized_Thermal_OTOC_output.close();
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
    for(i=0;i<nearby_state_index_size;i++){
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

    for(i=0;i<d.nmodes[0];i++){
        delete [] n_offdiag_total[i];
        delete [] n_offdiag_total_real[i];
        delete [] n_offdiag_total_imag[i];
        delete [] n_offdiag[i];
        delete [] n_offdiag_real[i];
        delete [] n_offdiag_imag[i];
        delete [] four_point_correlation_function_for_each_states[i];

        for(j=0;j<state_for_average_size;j++){
            delete [] n_offdiag_total_for_states_ensemble[i][j];
            delete [] n_offdiag_total_for_states_ensemble_real[i][j];
            delete [] n_offdiag_total_for_states_ensemble_imag[i][j];
            delete [] n_offdiag_for_states_ensemble[i][j];
            delete [] n_offdiag_for_states_ensemble_real[i][j];
            delete [] n_offdiag_for_states_ensemble_imag[i][j];
        }
        delete [] n_offdiag_total_for_states_ensemble[i];
        delete [] n_offdiag_total_for_states_ensemble_real[i];
        delete [] n_offdiag_total_for_states_ensemble_imag[i];
        delete [] n_offdiag_for_states_ensemble[i];
        delete [] n_offdiag_for_states_ensemble_real[i];
        delete [] n_offdiag_for_states_ensemble_imag[i];
    }
    delete [] n_offdiag_total;
    delete [] n_offdiag_total_real;
    delete [] n_offdiag_total_imag;
    delete [] n_offdiag;
    delete [] n_offdiag_real;
    delete [] n_offdiag_imag;

    delete [] n_offdiag_element;
    delete [] n_offdiag_element_real;
    delete [] n_offdiag_element_imag;
    delete [] n_offdiag_element_all_pc;
    delete [] n_offdiag_element_real_all_pc;
    delete [] n_offdiag_element_imag_all_pc;

    delete [] n_offdiag_total_for_states_ensemble;
    delete [] n_offdiag_total_for_states_ensemble_real;
    delete [] n_offdiag_total_for_states_ensemble_imag;
    delete [] n_offdiag_for_states_ensemble;
    delete [] n_offdiag_for_states_ensemble_real;
    delete [] n_offdiag_for_states_ensemble_imag;

    // n_offdiag one mode quanta below.
    delete [] another_OTOC;
    for(i=0;i<d.nmodes[0];i++){
        delete [] n_offdiag_total_one_mode_quanta_below[i];
        delete [] n_offdiag_total_one_mode_quanta_below_real[i];
        delete [] n_offdiag_total_one_mode_quanta_below_imag[i];
        delete [] n_offdiag_one_mode_quanta_below[i];
        delete [] n_offdiag_one_mode_quanta_below_real[i];
        delete [] n_offdiag_one_mode_quanta_below_imag[i];
    }
    delete [] n_offdiag_one_mode_quanta_below;
    delete [] n_offdiag_one_mode_quanta_below_real;
    delete [] n_offdiag_one_mode_quanta_below_imag;
    delete [] n_offdiag_total_one_mode_quanta_below;
    delete [] n_offdiag_total_one_mode_quanta_below_real;
    delete [] n_offdiag_total_one_mode_quanta_below_imag;

    delete [] four_point_correlation_function_at_initial_state;
    delete [] four_point_correlation_function_average_over_states;
    delete [] four_point_correlation_function_variance_over_states;
    delete [] four_point_correlation_function_for_each_states;

    for(i=0;i<d.nmodes[0];i++){
        delete [] Stability_Matrix[i];
    }
    delete [] Stability_Matrix;



    delete [] inverse_IPR_in_one_process;
    delete [] inverse_IPR_all;
    delete [] IPR_all;

    delete [] other_state_overlap_with_initial_state;

    // delete variable for computing Lyapunovian for xp.

    d.delete_variable_for_computing_Lyapunovian_xp();
    delete [] xd_for_xp;
    delete [] yd_for_xp;
    for(i=0;i<2*d.nmodes[0];i++){
        delete [] Lyapunov_spectrum_for_xp[i];
        delete [] Lyapunov_spectrum_for_xp_from_single_state[i];
        delete [] Matrix_M[i];
    }
    delete [] Lyapunov_spectrum_for_xp;
    delete [] Lyapunov_spectrum_for_xp_from_single_state;
    delete [] Matrix_M;

    delete [] Two_point_func1;
    delete [] Two_point_func2;
    for(i=0;i<2*d.nmodes[0];i++){
        delete [] TOC_per_state[i];
        delete [] TOC[i];
    }
    delete [] TOC_per_state;
    delete [] TOC;


    delete [] mode_quanta;
    delete [] total_mode_quanta;
    delete [] start_time;
    delete [] final_time;



};
