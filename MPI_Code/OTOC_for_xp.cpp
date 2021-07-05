//
// Created by phyzch on 12/21/20.
#include "../system.h"
#include "../util.h"
using namespace std;

complex<double> detector:: compute_c_overlap(int state_m, int relative_position_to_initial_state, int mode_k,
                                             vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                             vector<vector<int>> * index_for_xp_sparsify){
    // compute <m| c_{k}(t) | l>   l is state of choice (could be l or l_{i}^{+}, not necessarily initial state)
    // xd_for_xp contain real part of state (for trajectory of initial state) after lowering or raising at mode k.
    // size : (2*d.nmodes[0] , dmatsize[0]). First N data is for lowering |j_{k}^{-}>, last N data is for raising |j_{k}^{+}>
    // yd_for_xp contain imaginary part of state. (for trajectory of initial state)
    // if our simulation do not contain that neighbor state , we set xd_for_xp and yd_for_xp for that state as 0
    // relative_position_to_initial_state is first index for xd_for_xp, yd_for_xp. If our state of choice l is initial state, relative_position == 0
    // else if it's state near initial state, it can be in range [0,nmodes[0]] if it has one quanta lower or in range[nmodes[0], 2 * nmodes[0]] if have one quanta higher
    int i, j ;
    complex<double> c_overlap = 0;  // c_overlap = <m| c_{k}(t) | l>
    complex<double> c_overlap_sum = 0;

    double real_c_overlap;
    double imag_c_overlap;

    double real_c_overlap_sum;
    double imag_c_overlap_sum;

    int xd_for_xp_sparsify_list_len ;
    int basis_set_index;
    double real_part_value;
    double imag_part_value;

    if(mode_k<nmodes[0]){
        // c_{k} = a_{k} : lowering operator. <m(t) | a_{k} | l(t)> = \sum_{j} <m(t)| | j> <j| a_{k} |l(t)> = \sum_{j} <m(t) | j> <j_{k}^{+} || l(t)> * sqrt(nj_{k} + 1)

        // Here because we have to compute <j_{k}^{+} | l(t)>, we choose [relative_position_to_initial_state] (stands for l(t))
        // [mode_k + nmodes[0]] stands for j_{k}^{+}.
        // we only compute <j_{k}^{+}| l(t)> which is larger than cutoff. So we do not have to spend time computing quantity near 0
        xd_for_xp_sparsify_list_len = xd_for_xp_sparsify[relative_position_to_initial_state][mode_k + nmodes[0]].size();
        for(j=0;j<xd_for_xp_sparsify_list_len;j++){
            // basis set index is |j> in c_{k} expression
            basis_set_index = index_for_xp_sparsify[relative_position_to_initial_state][mode_k + nmodes[0]][j];
            real_part_value = xd_for_xp_sparsify[relative_position_to_initial_state][mode_k + nmodes[0]][j];
            imag_part_value = yd_for_xp_sparsify[relative_position_to_initial_state][mode_k + nmodes[0]][j];
            c_overlap = c_overlap + sqrt(dv[0][basis_set_index][mode_k] + 1) *
                                    complex<double> (real_part_value,imag_part_value) *
                                    complex<double>(xd[state_m][basis_set_index], -yd[state_m][basis_set_index]);
        }


    }
    else{
        // c_{k} = a_{k}^{+} :raising operator: <m(t) | a_{k}^{+} | l(t)> = \sum_{j} <m(t)|j> <j| a_{k}^{+} |l(t)>
        // I have mode_k-nmodes[0] because we want to compute <j_{k}^{-}|l(t)> and mode_k > nmodes[0] corresponding to j_{k}^{+}

        xd_for_xp_sparsify_list_len = xd_for_xp_sparsify[relative_position_to_initial_state][mode_k - nmodes[0]].size();
        for(j=0;j<xd_for_xp_sparsify_list_len; j++){
            basis_set_index = index_for_xp_sparsify[relative_position_to_initial_state][mode_k - nmodes[0]][j];
            real_part_value = xd_for_xp_sparsify[relative_position_to_initial_state][mode_k-nmodes[0]][j];
            imag_part_value = yd_for_xp_sparsify[relative_position_to_initial_state][mode_k - nmodes[0] ][j];
            c_overlap = c_overlap + sqrt(dv[0][basis_set_index][mode_k - nmodes[0]]) *
                                    complex<double> (real_part_value , imag_part_value)
                                    *complex<double> (xd[state_m][basis_set_index] , -yd[state_m][basis_set_index]);
        }

    }

    real_c_overlap = real(c_overlap);
    imag_c_overlap = imag(c_overlap);

    MPI_Allreduce(&real_c_overlap,&real_c_overlap_sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&imag_c_overlap, &imag_c_overlap_sum,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    c_overlap_sum = complex<double>(real_c_overlap_sum, imag_c_overlap_sum);
    return c_overlap_sum;
}

void detector:: compute_M_matrix(int state_m, int state_l, complex<double> ** M_matrix,
                                 vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                 vector<vector<int>> * index_for_xp_sparsify){
    // M_{ki}(t) = [c_{k}(t) , c_{i} ]
    int i,j,k;
    int m_one_mode_quanta_below;
    int m_one_mode_quanta_above;

    complex<double> C1;   // C1 is <m| c_{k}(t) |l_{i}^{-}> if i<N, <m| c_{k}(t) |l_{i}^{+}> if l >= N
    complex<double> C2;   // C2 is <m_{i}^{+} | c_{k}(t) |l> if i<N, <m_{i}^{-}|c_{k}(t) | l> if l >= N

    for(k=0;k<2*nmodes[0];k++){
        for(i=0;i<2*nmodes[0];i++){
            if( i < nmodes[0] ){
                // c_{i} == a_{i}, lowering operator

                m_one_mode_quanta_above = neighbor_state_in_nearby_state_index_list[state_m][i+nmodes[0]];

                if(neighbor_state_in_nearby_state_index_list[initial_state_index_in_nearby_state_index_list][i] != -1){
                    C1 = compute_c_overlap(state_m,i + 1,k,
                                           xd_for_xp_sparsify,
                                           yd_for_xp_sparsify,
                                           index_for_xp_sparsify);  // i+1 stands for |l_{i}^{-}> state
                }
                else{
                    // state do not exist. set C1 to 0.
                    C1 = 0;
                }

                if(m_one_mode_quanta_above !=-1){
                    C2 = compute_c_overlap(m_one_mode_quanta_above,0,k,
                                           xd_for_xp_sparsify,
                                           yd_for_xp_sparsify,
                                           index_for_xp_sparsify);
                }
                else{
                    // state do not exist, set C2 = 0.
                    C2 = 0;
                }
                M_matrix[k][i] = sqrt( dv_all[0][nearby_state_index[state_l]][i] ) * C1 - sqrt(dv_all[0][nearby_state_index[state_m]][i] + 1) * C2;
            }
            else{
                // c_{i} == a_{i}^{+} , raising operator.  i>=nmodes[0]
                m_one_mode_quanta_below = neighbor_state_in_nearby_state_index_list[state_m][i-nmodes[0]];
                if(neighbor_state_in_nearby_state_index_list[initial_state_index_in_nearby_state_index_list][i] != -1){
                    C1 = compute_c_overlap(state_m,i + 1,k,
                                           xd_for_xp_sparsify,
                                           yd_for_xp_sparsify,
                                           index_for_xp_sparsify);  // i+1 stands for <l_{i}^{+}> state
                }
                else{
                    // state do not exist, set C1 = 0.
                    C1 = 0;
                }

                if(m_one_mode_quanta_below != -1){
                    C2 = compute_c_overlap(m_one_mode_quanta_below,0,k,
                                           xd_for_xp_sparsify,
                                           yd_for_xp_sparsify,
                                           index_for_xp_sparsify);
                }
                else{
                    // state do not exist, set C2 = 0.
                    C2 = 0;
                }
                M_matrix[k][i] = sqrt(dv_all[0][nearby_state_index[state_l]][i - nmodes[0] ]+1) * C1 - sqrt(dv_all[0][nearby_state_index[state_m]][i - nmodes[0] ]) * C2;
            }
        }
    }
}

void detector:: compute_Lyapunov_spectrum_for_xp(complex<double>  ** Lyapunov_spectrum_for_xp,
                                                 complex<double>  ** Lyapunov_spectrum_for_xp_from_single_state,
                                                 complex<double> ** M_matrix,
                                                 vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp,
                                                 vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                                 vector<vector<int>> * index_for_xp_sparsify,
                                                 double t , ofstream & Every_states_contribution_to_OTOC_xp){
    /*
    because after applying raising and lowering operator, our state may not in same process, so we have to create array to store
     states after we apply raising and lowering operator.
    xd_for_xp:  size (2*nmodes[0] + 1, nmodes[0], dmatsize[0] ).  2*nmodes[0] + 1 stands for trajectory of state |l, l^{-}_{1}, cdots, l^{+}_{1} )
    nmodes[0] stands for states one quanta below or above current state j: <j|l(t)>
     dmatsize[0]: stands for array size for states.
     */
    int i,j,k,m;
    int nearby_state_index_size = nearby_state_index.size();
    int   sparsify_M_matrix_element_len;
    int mode_index1, mode_index2;
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0];j++){
            Lyapunov_spectrum_for_xp[i][j] = 0;
        }
    }

    double cutoff_criteria = pow(10,-3);

    update_xd_yd_for_xp(xd_for_xp,yd_for_xp,xd_for_xp_sparsify,yd_for_xp_sparsify,index_for_xp_sparsify, cutoff_criteria);

    for(m=0;m<nearby_state_index_size;m++){

        // go through state m

        compute_M_matrix(m,initial_state_index_in_nearby_state_index_list,M_matrix,
                         xd_for_xp_sparsify, yd_for_xp_sparsify, index_for_xp_sparsify);

        if(my_id == 0){

            // sparsify M_matrix and speed up simulation.
            vector<vector<complex<double>>> M_matrix_sparsify ;
            vector<vector<int>> M_matrix_sparsify_index;
            for(i=0;i<2*nmodes[0];i++){
                vector<complex<double>> M_matrix_sparsify_element;
                vector<int> M_matrix_sparsify_index_element;
                for(j=0;j<2*nmodes[0];j++){

                    if(abs(M_matrix[i][j]) > cutoff_criteria){
                        M_matrix_sparsify_element.push_back(M_matrix[i][j]);
                        M_matrix_sparsify_index_element.push_back(j);
                    }

                }
                M_matrix_sparsify.push_back(M_matrix_sparsify_element);
                M_matrix_sparsify_index.push_back(M_matrix_sparsify_index_element);

            }

            for(k=0;k<2*nmodes[0];k++){
                sparsify_M_matrix_element_len = M_matrix_sparsify_index[k].size();
                for(i=0;i<sparsify_M_matrix_element_len; i++){
                    mode_index1 = M_matrix_sparsify_index[k][i];
                    for(j=0;j<sparsify_M_matrix_element_len;j++){
                        mode_index2 = M_matrix_sparsify_index[k][j];

                        Lyapunov_spectrum_for_xp[mode_index1][mode_index2] = Lyapunov_spectrum_for_xp[mode_index1][mode_index2] +
                                                                             std::conj(M_matrix[k][mode_index1]) * M_matrix[k][mode_index2];

                    }
                }

            }

        }



    }

    if(my_id == 0){
        cout << "finish computing one step " << endl;
    }
}

void detector:: update_xd_yd_for_xp(vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp,
                                    vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                    vector<vector<int>> * index_for_xp_sparsify, double cutoff_criteria
) {
    // xd_for_xp : size (2*nmodes[0] + 1, nmodes[0], dmatsize[0] )
    // yd_for_xp : size (2*nmodes[0] + 1, nmodes[0] , dmatsize[0])
    // first index is for states l of choice: |l, l_{k}^{-}, l_{k}^{+}>
    // second index is for raising and lowering operator we apply to state j of choice |j_{k}^{-} , j_{k}^{+} >
    // third index is for state j.
    // for example real part for  <j_{mode_k}^{+} | l_{mode_i}^{-}> = xd_for_xp[mode_i + 1] [mode_k + nmodes[0]] [j]

    // xd_for_xp_sparsify, yd_for_xp_sparsify , index_for_xp_sparsify : only record matrix element that is big to speed up simulation.
    // cutoff_criteria : criteria to be included in xd_for_xp_sparisfy

    int i, j, k;
    int state_l_index;
    int local_state_index; // state index in process
    int begin_index = dmatsize_offset_each_process[0][my_id];

    // clear xd_for_xp_sparsify etc.
    for (i = 0; i < 2 * nmodes[0] + 1; i++) {
        xd_for_xp_sparsify[i].clear();
        yd_for_xp_sparsify[i].clear();
        index_for_xp_sparsify[i].clear();
    }


    for (i = 0; i < 2 * nmodes[0] + 1; i++) {
        // i is index for trajectory
        if (i == 0) {
            state_l_index = initial_state_index_in_nearby_state_index_list;
        } else {
            // this is state index with one quanta difference in quantum number
            state_l_index = neighbor_state_in_nearby_state_index_list[initial_state_index_in_nearby_state_index_list][
                    i - 1];
        }

        if (state_l_index != -1) {
            for (j = 0; j < 2 * nmodes[0]; j++) {
                // j is index for move in quantum number space for state j.
                for (k = 0; k < to_send_buffer_len_list_for_xp[j]; k++) {
                    local_state_index = tosendVecIndex_for_xp[j][k] - begin_index;
                    send_xd_for_xp[i][j][k] = xd[state_l_index][local_state_index];
                    send_yd_for_xp[i][j][k] = yd[state_l_index][local_state_index];
                }
            }
            // send data.
            for (j = 0; j < 2 * nmodes[0]; j++) {
                MPI_Alltoallv(&send_xd_for_xp[i][j][0], tosendVecCount_for_xp[j], tosendVecPtr_for_xp[j], MPI_DOUBLE,
                              &receive_xd_for_xp[i][j][0], remoteVecCount_for_xp[j], remoteVecPtr_for_xp[j], MPI_DOUBLE,
                              MPI_COMM_WORLD);

                MPI_Alltoallv(&send_yd_for_xp[i][j][0], tosendVecCount_for_xp[j], tosendVecPtr_for_xp[j], MPI_DOUBLE,
                              &receive_yd_for_xp[i][j][0], remoteVecCount_for_xp[j], remoteVecPtr_for_xp[j], MPI_DOUBLE,
                              MPI_COMM_WORLD);
            }

            // now use Index_in_remoteVEcIndex_for_xp to construct one-to-one correspondence between receive_xd_for_xp and xd_for_xp.
            for (j = 0; j < 2 * nmodes[0]; j++) {
                for (k = 0; k < dmatsize[0]; k++) {

                    if (Index_in_remoteVecIndex_for_xp[j][k] == -1) {
                        // this state do not exist, we assign xd,yd = 0 in this case
                        xd_for_xp[i][j][k] = 0;
                        yd_for_xp[i][j][k] = 0;
                    } else {
                        xd_for_xp[i][j][k] = receive_xd_for_xp[i][j][Index_in_remoteVecIndex_for_xp[j][k]];
                        yd_for_xp[i][j][k] = receive_yd_for_xp[i][j][Index_in_remoteVecIndex_for_xp[j][k]];
                    }

                }
            }

            double sparsify_ratio;
            int state_num = 0;
            // sparsify xd_for_xp (real part of initial state's nearby state's trajectory) and yd_for_xp
            for (j = 0; j < 2 * nmodes[0]; j++) {
                vector<double> xd_for_xp_sparsify_one_mode;
                vector<double> yd_for_xp_sparsify_one_mode;
                vector<int> index_for_xp_sparsify_one_mode;
                for (k = 0; k < dmatsize[0]; k++) {

                    if (abs(xd_for_xp[i][j][k]) > cutoff_criteria or abs(yd_for_xp[i][j][k]) > cutoff_criteria) {
                        xd_for_xp_sparsify_one_mode.push_back(xd_for_xp[i][j][k]);
                        yd_for_xp_sparsify_one_mode.push_back(yd_for_xp[i][j][k]);
                        index_for_xp_sparsify_one_mode.push_back(k);

                        state_num = state_num + 1;
                    }

                }

                xd_for_xp_sparsify[i].push_back(xd_for_xp_sparsify_one_mode);
                yd_for_xp_sparsify[i].push_back(yd_for_xp_sparsify_one_mode);
                index_for_xp_sparsify[i].push_back(index_for_xp_sparsify_one_mode);

            }
            sparsify_ratio = double(state_num) / (2 * nmodes[0] * dmatsize[0]);

        } else {
            // state do not exist. We will not use xd_for_xp and yd_for_xp. Thus we do not have to distribute data there.
            ;
        }


    }
}

vector<int> detector::construct_receive_buffer_index_for_xp(int ** remoteVecCount_for_xp, int ** remoteVecPtr_for_xp, int ** remoteVecIndex_for_xp,
                                                            int ** Index_in_remoteVecIndex_for_xp){
    /*
     * remoteVecCount_for_xp : size [2 * nmodes[0] , num_proc]
     * remoteVecPtr_for_xp : size [2 * nmodes[0] , num_proc]
     * remoteVecIndex_for_xp: size [2* nmodes[0], dmatsize[0] ]
     * Index_in_remote_VecIndex_for_xp : record given state's index in remoteVecIndex. used for Transfer received state to xd_for_xp, yd_for_xp.
    */
    int i, j;
    int element_number;
    int remote_pc_id;
    int vsize = total_dmat_size[0]/num_proc;
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<num_proc;j++){
            remoteVecCount_for_xp[i][j] = 0;
        }
    }

    int * search_Ind;
    int to_receive_buffer_len;
    vector<int> to_receive_buffer_len_list; // record receive_buffer_len for different raising and lowering operator.

    vector<int> index_for_neighbor_state;
    for(i=0;i<2*nmodes[0];i++){
        // i is index for raising and lowering state
        index_for_neighbor_state.clear();
        for(j=0;j<dmatsize[0];j++){
            // find all state in process 's neighbor state (after operation of given operator).
            // Notice in neighbor_state_index_for_all_state_list we record global index, thus we need to use dirow[0][j]
            index_for_neighbor_state.push_back(neighbor_state_index_for_all_state_list[dirow[0][j]][i]);
        }
        vector<int> index_for_neighbor_state_copy = index_for_neighbor_state;
        sort(index_for_neighbor_state_copy.begin(), index_for_neighbor_state_copy.end()); // sort index for neighbor_state_copy
        // now states in different process is in order, we can compute remoteVecCount and remoteVecIndex now
        element_number = 0;
        for(j=0;j<dmatsize[0];j++){
            if(index_for_neighbor_state_copy[j] != -1){
                // this nearby state exists

                if(index_for_neighbor_state_copy[j] >= vsize * (num_proc-1)){
                    remote_pc_id = num_proc -1 ;
                }
                else{
                    remote_pc_id = index_for_neighbor_state_copy[j] / vsize;
                }
                remoteVecCount_for_xp[i][remote_pc_id] ++ ;
                remoteVecIndex_for_xp[i][element_number] = index_for_neighbor_state_copy[j];
                element_number ++ ;
            }

        }

        remoteVecPtr_for_xp[i][0] = 0;
        for(j=1;j<num_proc;j++){
            remoteVecPtr_for_xp[i][j] = remoteVecPtr_for_xp[i][j-1] + remoteVecCount_for_xp[i][j-1];
        }

        to_receive_buffer_len = 0;
        for(j=0;j<num_proc;j++){
            to_receive_buffer_len = to_receive_buffer_len + remoteVecCount_for_xp[i][j];
        }
        to_receive_buffer_len_list.push_back(to_receive_buffer_len);

        // compute Index_in_remoteVecIndex_for_xp
        for(j=0;j<dmatsize[0];j++){

            if(index_for_neighbor_state[j] != -1){
                search_Ind = (int *) bsearch( &index_for_neighbor_state[j],remoteVecIndex_for_xp[i], to_receive_buffer_len,sizeof(int),compar );
                Index_in_remoteVecIndex_for_xp[i][j] = search_Ind - remoteVecIndex_for_xp[i];
            }
            else{
                Index_in_remoteVecIndex_for_xp[i][j] = -1;
            }

        }

    }

    return to_receive_buffer_len_list;
}

vector<int> detector:: construct_send_buffer_index_for_xp(int ** remoteVecCount_for_xp, int ** remoteVecPtr_for_xp, int ** remoteVecIndex_for_xp,
                                               int ** tosendVecCount_for_xp, int ** tosendVecPtr_for_xp, int ** tosendVecIndex_for_xp){
    /*
 * remoteVecCount_for_xp : size [2 * nmodes[0] , num_proc]
 * remoteVecPtr_for_xp : size [2 * nmodes[0] , num_proc]
 * remoteVecIndex_for_xp: size [2* nmodes[0], dmatsize[0] ]
 * tosendVecCount_for_xp: size [2* nmodes[0], num_proc]
 * tosendVecPtr_for_xp : size [2*nmodes[0] , num_proc]
 * tosendVecIndex_for_xp : size [2*nmodes[0], dmatsize[0] ]
     */
     int i,j;
     int to_send_buffer_len;
     vector<int> to_send_buffer_len_list;
     for(i=0;i<2 * nmodes[0]; i++){
         MPI_Alltoall(&remoteVecCount_for_xp[i][0],1,MPI_INT,&tosendVecCount_for_xp[i][0],1,MPI_INT,MPI_COMM_WORLD);

         // compute length of data to send and record it in to_send_buffer_len_list
         to_send_buffer_len = 0;
         for(j=0;j<num_proc;j++){
             to_send_buffer_len = to_send_buffer_len + tosendVecCount_for_xp[i][j];
         }
         to_send_buffer_len_list.push_back(to_send_buffer_len);

         // compute displacement
         tosendVecPtr_for_xp[i][0] = 0;
         for(j=1;j<num_proc;j++){
             tosendVecPtr_for_xp[i][j] = tosendVecPtr_for_xp[i][j-1] + tosendVecCount_for_xp[i][j-1];
         }

         // compute index of element to send to other process.
         MPI_Alltoallv(&remoteVecIndex_for_xp[i][0],remoteVecCount_for_xp[i],remoteVecPtr_for_xp[i],MPI_INT,
                       &tosendVecIndex_for_xp[i][0],tosendVecCount_for_xp[i],tosendVecPtr_for_xp[i],MPI_INT,MPI_COMM_WORLD);
     }

     return to_send_buffer_len_list;
}

void detector::prepare_computing_Lyapunovian_for_xp() {
    int i,j;
    remoteVecCount_for_xp = new int * [2 * nmodes[0]];
    remoteVecPtr_for_xp = new int * [2* nmodes[0]];
    remoteVecIndex_for_xp = new int * [2*nmodes[0]];
    Index_in_remoteVecIndex_for_xp = new int * [2* nmodes[0]];

    tosendVecCount_for_xp = new int * [2*nmodes[0]];
    tosendVecPtr_for_xp = new int * [2*nmodes[0]];
    tosendVecIndex_for_xp = new int * [2*nmodes[0]];

    for(i=0;i<2*nmodes[0];i++){
        remoteVecCount_for_xp[i] = new int [num_proc];
        remoteVecPtr_for_xp[i] = new int [num_proc];
        remoteVecIndex_for_xp[i] = new int [dmatsize[0]];
        Index_in_remoteVecIndex_for_xp[i] = new int [dmatsize[0]];

        tosendVecCount_for_xp[i] =  new int [num_proc];
        tosendVecPtr_for_xp[i] =  new int [num_proc];
        tosendVecIndex_for_xp[i] = new int [dmatsize[0]];
    }


    to_receive_buffer_len_list_for_xp = construct_receive_buffer_index_for_xp(remoteVecCount_for_xp,
                                                                              remoteVecPtr_for_xp,remoteVecIndex_for_xp,
                                                                              Index_in_remoteVecIndex_for_xp);

    to_send_buffer_len_list_for_xp = construct_send_buffer_index_for_xp(remoteVecCount_for_xp,remoteVecPtr_for_xp,remoteVecIndex_for_xp,
                                                                 tosendVecCount_for_xp,tosendVecPtr_for_xp,tosendVecIndex_for_xp);

    send_xd_for_xp = new double **[2*nmodes[0] + 1];
    send_yd_for_xp = new double ** [2*nmodes[0] + 1];

    receive_xd_for_xp = new double ** [2*nmodes[0] + 1];
    receive_yd_for_xp = new double ** [2*nmodes[0] + 1];

    for(i=0;i<2*nmodes[0]+1;i++){
        send_xd_for_xp[i] = new double * [2*nmodes[0]];
        send_yd_for_xp[i] = new double * [2*nmodes[0]];

        receive_xd_for_xp[i] = new double * [2*nmodes[0]];
        receive_yd_for_xp[i] = new double * [2*nmodes[0]];

        for(j=0;j<2*nmodes[0];j++){
            send_xd_for_xp[i][j] = new double [dmatsize[0]];
            send_yd_for_xp[i][j] = new double [dmatsize[0]];

            receive_xd_for_xp[i][j] = new double[dmatsize[0]];
            receive_yd_for_xp[i][j] = new double [dmatsize[0]];
        }
    }

}

void detector::delete_variable_for_computing_Lyapunovian_xp(){
    int i,j;
    for(i=0;i<2*nmodes[0];i++){
        delete [] remoteVecCount_for_xp[i];
        delete [] remoteVecPtr_for_xp[i];
        delete [] remoteVecIndex_for_xp[i];

        delete [] tosendVecCount_for_xp[i];
        delete [] tosendVecPtr_for_xp[i];
        delete [] tosendVecIndex_for_xp[i];

        delete [] Index_in_remoteVecIndex_for_xp[i];
    }
    delete [] remoteVecCount_for_xp;
    delete [] remoteVecPtr_for_xp;
    delete [] remoteVecIndex_for_xp;

    delete [] tosendVecCount_for_xp;
    delete [] tosendVecPtr_for_xp;
    delete [] tosendVecIndex_for_xp;

    delete [] Index_in_remoteVecIndex_for_xp;

    for(i=0 ; i< 2*nmodes[0] + 1 ; i++){
        for(j=0 ; j < 2*nmodes[0] ; j++){
            delete [] send_xd_for_xp[i][j];
            delete [] send_yd_for_xp[i][j];
            delete [] receive_xd_for_xp[i][j];
            delete [] receive_yd_for_xp[i][j];
        }
        delete [] send_xd_for_xp[i];
        delete [] send_yd_for_xp[i];
        delete [] receive_xd_for_xp[i];
        delete [] receive_yd_for_xp[i];
    }

    delete [] send_xd_for_xp;
    delete [] send_yd_for_xp;
    delete [] receive_xd_for_xp;
    delete [] receive_yd_for_xp;

}