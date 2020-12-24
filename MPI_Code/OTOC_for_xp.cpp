//
// Created by phyzch on 12/21/20.
#include "../system.h"
#include "../util.h"
using namespace std;

complex<double> detector:: compute_c_overlap(int state_m, int relative_position_to_initial_state, int mode_k,
                                             vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp){
    // compute <m| c_{k}(t) | l>   l is state of choice (could be l or l_{i}^{+}, not necessarily initial state)
    // xd_for_xp contain real part of state (for trajectory of initial state) after lowering or raising at mode k.
    // size : (2*d.nmodes[0] , dmatsize[0]). First N data is for lowering |j_{k}^{-}>, last N data is for raising |j_{k}^{+}>
    // yd_for_xp contain imaginary part of state. (for trajectory of initial state)
    // if our simulation do not contain that neighbor state , we set xd_for_xp and yd_for_xp for that state as 0
    // relative_position_to_initial_state is first index for xd_for_xp, yd_for_xp. If our state of choice l is initial state, relative_position == 0
    // else if it's state near initial state, it can be in range [0,nmodes[0]] if it has one quanta lower or in range[0,range[0]] if have one quanta higher
    int i, j ;
    complex<double> c_overlap = 0;  // c_overlap = <m| c_{k}(t) | l>
    if(mode_k<nmodes[0]){
        // c_{k} = a_{k} : lowering operator. <m(t) | a_{k} | l(t)> = \sum_{j} <m(t) | j> <j| a_{k} |l(t)> = \sum_{j} <m(t) | j> <j_{k}^{+} | l(t)> * sqrt(nj_{k} + 1)
        for(j=0;j<dmatsize[0];j++){
            c_overlap = c_overlap + sqrt(dv[0][j][mode_k] + 1) *
                    complex<double>(xd_for_xp[relative_position_to_initial_state][mode_k + nmodes[0]][j] , yd_for_xp[relative_position_to_initial_state][mode_k + nmodes[0]][j])
                            * complex<double> (xd[state_m][j] , - yd[state_m][j]);
            // I have mode_k + nmodes[0] because we want to compute <j_{k}^{+} | l(t)> and mode_k <nmodes[0] corresponding to j_{k}^{-}
        }
    }
    else{
        // c_{k} = a_{k}^{+} :raising operator: <m(t) | a_{k}^{+} | l(t)> = \sum_{j} <m(t)|j> <j| a_{k}^{+} |l(t)>
        for(j=0;j<dmatsize[0];j++){
            c_overlap = c_overlap + sqrt(dv[0][j][mode_k]) *
                    complex<double> (xd_for_xp[relative_position_to_initial_state][mode_k-nmodes[0]][j] , yd_for_xp[relative_position_to_initial_state][mode_k - nmodes[0]][j])
                            * complex<double> (xd[state_m][j], -yd[state_m][j]);
            // I have mode_k-nmodes[0] because we want to compute <j_{k}^{-}|l(t)> and mode_k > nmodes[0] corresponding to j_{k}^{+}
        }
    }
    return c_overlap;
}

void detector:: compute_M_matrix(int state_m, int state_l, complex<double> ** M_matrix,
                                 vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp){
    // M_{ki}(t) = [c_{k}(t) , c_{i} ]
    int i,j,k;
    int m_one_mode_quanta_below;
    int m_one_mode_quanta_above;
    int l_one_mode_quanta_below;
    int l_one_mode_quanta_above;
    complex<double> C1;   // C1 is <m| c_{k}(t) |l_{i}^{-}> if i<N, <m| c_{k}(t) |l_{i}^{+}> if l >= N
    complex<double> C2;   // C2 is <m_{i}^{+} | c_{k}(t) |l> if i<N, <m_{i}^{-}|c_{k}(t) | l> if l >= N

    for(k=0;k<2*nmodes[0];k++){
        for(i=0;i<2*nmodes[0];i++){
            if( i < nmodes[0] ){
                // c_{i} == a_{i}, lowering operator
                l_one_mode_quanta_below = neighbor_state_in_nearby_state_index_list[state_l][i];
                m_one_mode_quanta_above = neighbor_state_in_nearby_state_index_list[state_m][i+nmodes[0]];
                C1 = compute_c_overlap(state_m,i + 1,k,xd_for_xp,yd_for_xp);  // i+1 stands for |l_{i}^{-}> state
                if(m_one_mode_quanta_above !=-1){
                    C2 = compute_c_overlap(m_one_mode_quanta_above,0,k,xd_for_xp,yd_for_xp);
                }
                else{
                    C2 = 0;
                }
                M_matrix[k][i] = sqrt( dv_all[0][nearby_state_index[state_l]][i] ) * C1 - sqrt(dv_all[0][nearby_state_index[state_m]][i] + 1) * C2;
            }
            else{
                // c_{i} == a_{i}^{+} , raising operator.  i>=nmodes[0]
                l_one_mode_quanta_above = neighbor_state_in_nearby_state_index_list[state_l][i];
                m_one_mode_quanta_below = neighbor_state_in_nearby_state_index_list[state_m][i-nmodes[0]];
                C1 = compute_c_overlap(state_m,i + 1,k,xd_for_xp,yd_for_xp);  // i+1 stands for <l_{i}^{+}> state
                if(m_one_mode_quanta_below != -1){
                    C2 = compute_c_overlap(m_one_mode_quanta_below,0,k,xd_for_xp,yd_for_xp);
                }
                else{
                    C2 = 0;
                }
                M_matrix[k][i] = sqrt(dv_all[0][nearby_state_index[state_l]][i]+1) * C1 - sqrt(dv_all[0][nearby_state_index[state_m]][i]) * C2;
            }
        }
    }
}

void detector:: compute_Lyapunov_spectrum_for_xp(complex<double>  ** Lyapunov_spectrum_for_xp, complex<double> ** M_matrix,
                                                 vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp){
    /*
    because after applying raising and lowering operator, our state may not in same process, so we have to create array to store
     states after we apply raising and lowering operator.
    xd_for_xp:  size (2*nmodes[0] + 1, nmodes[0], dmatsize[0] ).  2*nmodes[0] + 1 stands for trajectory of state |l, l^{-}_{1}, cdots, l^{+}_{1} )
    nmodes[0] stands for states one quanta below or above current state j: <j|l(t)>
     dmatsize[0]: stands for array size for states.
     */
    int i,j,k,m;
    int nearby_state_index_size = nearby_state_index.size();
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0];j++){
            Lyapunov_spectrum_for_xp[i][j] = 0;
        }
    }

    update_xd_yd_for_xp(xd_for_xp,yd_for_xp);

    for(m=0;m<nearby_state_index_size;m++){

        // go through state m
        if(bool_neighbor_state_all_in_nearby_state_index[m]){
            compute_M_matrix(m,initial_state_index_in_nearby_state_index_list,M_matrix,xd_for_xp,yd_for_xp);
            for(i=0;i<2*nmodes[0];i++){
                for(j=0;j<2*nmodes[0];j++){
                    for(k=0;k<2*nmodes[0];k++){
                        Lyapunov_spectrum_for_xp[i][j] = Lyapunov_spectrum_for_xp[i][j] +
                                                         std::conj(M_matrix[k][i]) * M_matrix[k][j];  // \sum_{k} (M_{ki})^{*} * M_{kj}
                    }
                }
            }
        }

    }
}

void detector:: update_xd_yd_for_xp(vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp){
    // xd_for_xp : size (2*nmodes[0] + 1, nmodes[0], dmatsize[0] )
    // yd_for_xp : size (2*nmodes[0] + 1, nmodes[0] , dmatsize[0])
    // first index is for states l of choice: |l, l_{k}^{-}, l_{k}^{+}>
    // second index is for raising and lowering operator we apply to state j of choice |j_{k}^{-} , j_{k}^{+} >
    // third index is for state j.
    // for example real part for  <j_{mode_k}^{+} | l_{mode_i}^{-}> = xd_for_xp[mode_i + 1] [mode_k + nmodes[0]] [j]

    int i,j,k;
    int state_l_index ;
    int local_state_index; // state index in process
    int begin_index = dmatsize_offset_each_process[0][my_id];
    for(i=0;i<2*nmodes[0]+1;i++){
        // i is index for trajectory
        if(i==0){
            state_l_index = initial_state_index_in_nearby_state_index_list;
        }
        else{
            // this is state index with one quanta difference in quantum number
            state_l_index = neighbor_state_in_nearby_state_index_list[initial_state_index_in_nearby_state_index_list][i-1];
        }

        for(j=0;j<2*nmodes[0];j++){
            // j is index for move in quantum number space for state j.
           for(k=0;k<to_send_buffer_len_list_for_xp[j];k++){
               local_state_index = tosendVecIndex_for_xp[j][k] - begin_index;
               send_xd_for_xp[i][j][k] = xd[state_l_index][local_state_index];
               send_yd_for_xp[i][j][k] = yd[state_l_index][local_state_index];
           }
        }

        for(j=0;j<2*nmodes[0];j++){
            MPI_Alltoallv(&send_xd_for_xp[i][j][0],tosendVecCount_for_xp[j],tosendVecPtr_for_xp[j],MPI_DOUBLE,
                          &receive_xd_for_xp[i][j][0],remoteVecCount_for_xp[j],remoteVecPtr_for_xp[j],MPI_DOUBLE,MPI_COMM_WORLD);
        }

        // now use Index_in_remoteVEcIndex_for_xp to construct one-to-one correspondence between receive_xd_for_xp and xd_for_xp.
        for(j=0;j<2*nmodes[0];j++){
            for(k=0;k<dmatsize[0];k++){

                if(Index_in_remoteVecIndex_for_xp[j][k] == -1){
                    // this state do not exist, we assign xd,yd = 0 in this case
                    xd_for_xp[i][j][k] = 0;
                    yd_for_xp[i][j][k] = 0;
                }
                else{
                    xd_for_xp[i][j][k] = receive_xd_for_xp[i][j][Index_in_remoteVecIndex_for_xp[j][k]];
                    yd_for_xp[i][j][k] = receive_yd_for_xp[i][j][Index_in_remoteVecIndex_for_xp[j][k]];
                }

            }
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

        // cmpute Index_in_remoteVecIndex_for_xp
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


    to_receive_buffer_len_list_for_xp = construct_receive_buffer_index_for_xp(remoteVecCount_for_xp,remoteVecPtr_for_xp,remoteVecIndex_for_xp,Index_in_remoteVecIndex_for_xp);

    to_send_buffer_len_list_for_xp = construct_send_buffer_index_for_xp(remoteVecIndex_for_xp,remoteVecPtr_for_xp,remoteVecIndex_for_xp,
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
        delete remoteVecCount_for_xp[i];
        delete remoteVecPtr_for_xp[i];
        delete remoteVecIndex_for_xp[i];

        delete tosendVecCount_for_xp[i];
        delete tosendVecPtr_for_xp[i];
        delete tosendVecIndex_for_xp[i];

        delete Index_in_remoteVecIndex_for_xp[i];
    }
    delete [] remoteVecCount_for_xp;
    delete [] remoteVecPtr_for_xp;
    delete [] remoteVecIndex_for_xp;

    delete [] tosendVecCount_for_xp;
    delete [] tosendVecPtr_for_xp;
    delete [] tosendVecIndex_for_xp;

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

}