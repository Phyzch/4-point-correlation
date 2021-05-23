//
// Created by phyzch on 3/31/21.
//
#include "util.h"
#include "system.h"
using namespace std;

void detector::compute_selected_eigenstate_index(){
    int i;
    for(i=0;i<eigenstate_num;i++){
        if(Eigenvalue_list[i] > Emin2 and Eigenvalue_list[i] < Emax2){
            selected_eigenstate_index.push_back(i);
        }
    }
    selected_eigenstate_num = selected_eigenstate_index.size();
}


void detector::construct_neighbor_state_index_list_for_all_state(){
    int i, j;
    // for neighbor_state_index_for_all_state_list
    vector<int> state_mode;
    vector<int> neighbor_state_mode;
    bool exist;
    int position;
    for(i=0;i<total_dmat_size[0];i++){
        state_mode = dv_all[0][i];
        vector<int> neighbor_state_index;
        for(j=0;j<nmodes[0];j++){
            neighbor_state_mode = state_mode;
            // one quanta lower in mode j
            neighbor_state_mode[j] = state_mode[j] - 1;
            position = find_position_for_insert_binary(dv_all[0],neighbor_state_mode,exist);
            if(exist){
                neighbor_state_index.push_back(position);
            }
            else{
                neighbor_state_index.push_back(-1);
            }

        }

        for(j=0;j<nmodes[0];j++){
            neighbor_state_mode = state_mode;
            // one quanta higher in mode j
            neighbor_state_mode[j] = state_mode[j] + 1;
            position = find_position_for_insert_binary(dv_all[0], neighbor_state_mode, exist);
            if(exist){
                neighbor_state_index.push_back(position);
            }
            else{
                neighbor_state_index.push_back(-1);
            }
        }

        neighbor_state_index_for_all_state_list.push_back(neighbor_state_index);
    }

// ----------- For debug ----------------------
    // Find maximum eigenstate index in all eigenstate.
    int * Maximum_eigenstate_index = new int [eigenstate_num];
    double maximum_eigenstate_component = 0;
    for(i=0;i<eigenstate_num;i++){
        maximum_eigenstate_component = 0;
        Maximum_eigenstate_index[i] = -1;
        for(j=0;j< total_dmat_size[0]; j++){
            if( abs(Eigenstate_list[i][j]) > maximum_eigenstate_component ){
                maximum_eigenstate_component = abs(Eigenstate_list[i][j]);
                Maximum_eigenstate_index[i] = j ;
            }
        }
    }

    i = 0;
    // ------------------ For debug ---------------------
}

void detector::Broadcast_eigenstate_and_eigenvalue(){
    int i;
    // broadcast total number of eigenstate found to other process
    MPI_Bcast(&eigenstate_num,1,MPI_INT,0,MPI_COMM_WORLD);

    // other process allocate space for eigenvalue and eigenvector
    if(my_id!=0){
        Eigenvalue_list = new double[eigenstate_num];
        Eigenstate_list = new double * [eigenstate_num];
        for(i=0;i<eigenstate_num ; i++){
            Eigenstate_list[i] = new double [total_dmat_size[0]];
        }
    }
    // Broadcast eigenvalue and eigenvector to other process
    MPI_Bcast(&Eigenvalue_list[0], eigenstate_num,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(i=0;i<eigenstate_num;i++){
        MPI_Bcast(&Eigenstate_list[i][0], total_dmat_size[0], MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    }
}

void detector::compute_eigenstate_energy_std(){
    int i , j , k;
    int index;
    int index1 ;
    Eigenstate_energy_std_list = new double [eigenstate_num];
    for(i = 0; i<eigenstate_num ; i++){
        Eigenstate_energy_std_list[i] = 0;
        // \sum (E_{n} - eigenvalue)^2 * |<n|\phi>|^2

        // As diagonal part is partitioned in different location in total_dmat. we have to find them in different place
        index1 = 0;
        for(j=0;j<num_proc;j++){
            for(k=0;k<dmatsize_each_process[0][j];k++){
              index = dmat_offset_each_process[0][j] + k;

              Eigenstate_energy_std_list[i] = Eigenstate_energy_std_list[i] +
                                                pow(Eigenvalue_list[i] - total_dmat[0][index] , 2) * pow(Eigenstate_list[i][index1] , 2);

              index1 = index1 + 1;

            }
        }

        Eigenstate_energy_std_list[i] = sqrt(Eigenstate_energy_std_list[i]);
    }
}

void detector:: compute_phi_ladder_operator_phi( ){
    int i, j, k ;
    int l, m;
    int local_eigenstate_num;
    if(my_id != num_proc - 1){
        local_eigenstate_num = eigenstate_num / num_proc ;
    }
    else{
        local_eigenstate_num =   eigenstate_num - (num_proc - 1) * int(eigenstate_num / num_proc) ;
    }

    int local_eigenstate_begin_index = int(eigenstate_num / num_proc)  * my_id  ;

    // size : [2 * nmodes[0] ,  local_eigenstate_num ,  eigenstate_num ]. First nmodes[0] index is for lowering operator and last nmodes[0] is for raising operator
    double *** local_phi_ladder_operator_phi = new double ** [2 * nmodes[0]];
    for(i=0;i< 2 * nmodes[0] ; i++){
        local_phi_ladder_operator_phi[i] = new double * [local_eigenstate_num];
        for(j=0;j< local_eigenstate_num ; j++){
            local_phi_ladder_operator_phi[i][j] = new double [ eigenstate_num];
        }
    }

    int eigenstate_m_index , eigenstate_l_index;
    double ladder_operator_energy_change ; // if ladder operator is a_{i} , energy change is -frequency[i]. if it's a_{i}^{+} , energy change is + frequency[i]
    // begin computing <\phi_m | a_i | \phi_l> where \phi_m and \phi_l is eigenstate
    double energy_difference;  // energy difference is |E_l + ladder_operator_energy_change - E_m|
    int basis_set_state_index;
    int nearby_basis_set_state_index;
    double local_ladder_operator_value;

    double Energy_sift_criteria = 0;

    for(i=0;i<2 * nmodes[0]; i++){
        for(m=0;m<local_eigenstate_num;m++){
            for(l=0;l<eigenstate_num;l++){
                if(i< nmodes[0]){
                    // lowering operator
                    ladder_operator_energy_change = - mfreq[0][i] ;
                }
                else{
                    // raising operator
                    ladder_operator_energy_change = + mfreq[0][i - nmodes[0]];
                }

                eigenstate_m_index = local_eigenstate_begin_index + m;
                eigenstate_l_index = l;

                energy_difference = abs( Eigenvalue_list[eigenstate_l_index] +
                        ladder_operator_energy_change - Eigenvalue_list[eigenstate_m_index] );


                Energy_sift_criteria =  5 * (Eigenstate_energy_std_list[eigenstate_l_index] + Eigenstate_energy_std_list[eigenstate_m_index]);
                if(Energy_sift_criteria < 10000){
                    Energy_sift_criteria = 10000;
                }
                if(energy_difference > Energy_sift_criteria ){
                    local_phi_ladder_operator_phi[i][m][l] = 0;
                }
                else {
                    // we need to compute <\phi_m | a_{i} | phi_l> or <\phi_m | a_{i}^{+} | \phi_l>
                    local_ladder_operator_value = 0;
                    for (j = 0; j < total_dmat_size[0]; j++) {
                        basis_set_state_index = j;
                        // for i < nmodes[0] , that's state_mode[i] -1. which corresponds to lowering operator.
                        // for i > nmodes[0], that's state_mode[i] + 1, corresponds to raising operator
                        nearby_basis_set_state_index = neighbor_state_index_for_all_state_list[j][i];
                        if (nearby_basis_set_state_index == -1) {
                            continue;
                        } else {
                            if (i < nmodes[0]) {
                                // lowering operator
                                local_ladder_operator_value = local_ladder_operator_value +
                                                              Eigenstate_list[eigenstate_l_index][basis_set_state_index] *
                                                              Eigenstate_list[eigenstate_m_index][nearby_basis_set_state_index] *
                                                              sqrt(dv_all[0][basis_set_state_index][i]);
                            } else {
                                // raising operator
                                local_ladder_operator_value = local_ladder_operator_value +
                                                              Eigenstate_list[eigenstate_l_index][basis_set_state_index] *
                                                              Eigenstate_list[eigenstate_m_index][nearby_basis_set_state_index] *
                                                              sqrt(dv_all[0][basis_set_state_index][i - nmodes[0]] + 1);
                            }

                        }

                    }
                    local_phi_ladder_operator_phi[i][m][l] = local_ladder_operator_value;

                }


            }
        }
    }

    // now we send local ladder operator to all process and free this space.
    double  ** local_phi_ladder_operator_send = new double * [2 * nmodes[0]];
    int index = 0;
    for(i=0;i< 2 * nmodes[0] ; i++ ){
        local_phi_ladder_operator_send[i] = new double[local_eigenstate_num * eigenstate_num] ;
        index = 0;
        for(j=0;j< local_eigenstate_num; j++ ){
            for(k=0;k<eigenstate_num ; k++){
                local_phi_ladder_operator_send[i][index] = local_phi_ladder_operator_phi[i][j][k];
                index ++ ;
            }
        }

    }

    double ** phi_ladder_operator_phi_1d = new double * [2 * nmodes[0]];
    for(i=0;i< 2* nmodes[0] ; i++ ){
        phi_ladder_operator_phi_1d [i] = new double [eigenstate_num * eigenstate_num ] ;
    }

    int * recv_count = new int [num_proc];
    int * displs = new int [num_proc];
    for(i=0;i<num_proc-1;i++){
        recv_count [i] = int(eigenstate_num / num_proc) * eigenstate_num ;
    }
    recv_count[num_proc - 1] = int(eigenstate_num - int(eigenstate_num / num_proc) * (num_proc - 1)) * eigenstate_num;
    displs[0] = 0;
    for(i=1;i<num_proc;i++){
        displs[i] = displs[i - 1] + recv_count[ i - 1 ];
    }

    for(i=0;i<2*nmodes[0];i++){
        MPI_Gatherv(& local_phi_ladder_operator_send[i][0], local_eigenstate_num * eigenstate_num , MPI_DOUBLE,
                    &phi_ladder_operator_phi_1d[i][0], &recv_count[0], &displs[0], MPI_DOUBLE,0,MPI_COMM_WORLD );

        MPI_Bcast(&phi_ladder_operator_phi_1d[i][0], eigenstate_num * eigenstate_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }


    // free space:
    for(i=0;i<2 * nmodes[0] ; i++){
        for(j=0;j<local_eigenstate_num;j++){
            delete [] local_phi_ladder_operator_phi[i][j];
        }
        delete [] local_phi_ladder_operator_phi[i];
    }
    delete [] local_phi_ladder_operator_phi;

    for(i=0;i<2*nmodes[0];i++){
        delete [] local_phi_ladder_operator_send[i];
    }
    delete [] local_phi_ladder_operator_send;
    delete [] recv_count;
    delete [] displs;

    // construct 3d array to store phi_ladder_operator_phi
    phi_ladder_operator_phi = new double ** [2 * nmodes[0]];
    for(i=0;i< 2* nmodes[0]; i++){
        phi_ladder_operator_phi[i] = new double * [eigenstate_num];
        for(j=0; j< eigenstate_num; j++ ){
            phi_ladder_operator_phi[i][j] = new double [ eigenstate_num ];
        }
    }

    for(i=0;i<2 * nmodes[0]; i++){
        index = 0;
        for(j=0;j<eigenstate_num; j++){
            for(k=0;k<eigenstate_num;k++){
                phi_ladder_operator_phi[i][j][k] = phi_ladder_operator_phi_1d[i][index];
                index ++ ;
            }
        }
    }

    Eigenstate_OTOC_sift_criteria = double(1) / (100) ;

    // construct phi_operator_phi_tuple_list : size : [2 * nmodes[0] , eigenstate_num ] (list)
    // phi_ladder_operator_phi_tuple_list record <l | a_{i} | m> for all state m. in this list, it will have state l which have value of operator larger than criteria.
    int Total_Element_num = 0;
    for(i=0; i< 2 * nmodes[0] ; i++ ){
        vector<vector<phi_operator_phi_tuple>> phi_operator_phi_tuple_for_dof ;
        for(m=0;m<eigenstate_num ; m++){
            vector< phi_operator_phi_tuple > phi_operator_phi_tuple_for_state;
            for(l=0;l<eigenstate_num;l++){

                if( abs(phi_ladder_operator_phi[i][l][m])  > Eigenstate_OTOC_sift_criteria){
                    Total_Element_num ++ ;
                    phi_operator_phi_tuple Tuple1 (l,phi_ladder_operator_phi[i][l][m]);
                    phi_operator_phi_tuple_for_state.push_back(Tuple1);
                }

            }

            phi_operator_phi_tuple_for_dof .push_back (phi_operator_phi_tuple_for_state );
        }
        phi_ladder_operator_phi_tuple_list.push_back(phi_operator_phi_tuple_for_dof);
    }

    double Non_zero_Ratio;
    Non_zero_Ratio = double(Total_Element_num) / ( 2 * nmodes[0] * eigenstate_num * eigenstate_num )  ;

    if(my_id == 0){
        cout << "for phi_operator_phi, nonzero ratio is :  " << Non_zero_Ratio << endl;
    }

    // free space
    for(i=0;i<2*nmodes[0];i++){
        delete [] phi_ladder_operator_phi_1d[i];
    }
    delete [] phi_ladder_operator_phi_1d;


}

int Binary_search_phi_operator_phi_tuple_complex(const vector<phi_operator_phi_tuple_complex> & List, int state_index){
    int Size = List.size();
    if(Size == 0){
        return -1;
    }
    int begin_index = 0;
    int end_index = Size - 1;

    int index = int(begin_index + end_index) / 2;
    while(begin_index < end_index){
        index = int(begin_index + end_index) / 2;
        if(state_index < List[index].eigenstate_index){
            end_index = index - 1;
            index = (begin_index + end_index) / 2;
        }
        else if (state_index > List[index].eigenstate_index){
            begin_index = index + 1;
            index = (begin_index + end_index ) / 2;
        }
        else{
            return index;
        }
    }

    if(state_index == List[index].eigenstate_index){
        return index;
    }
    else{
        return -1;
    }

}

void detector:: compute_Eigenstate_OTOC_submodule(ofstream & Eigenstate_OTOC_output, double time, double *** Eigenstate_OTOC ,
                                                  double *** local_Eigenstate_OTOC,
                                                  complex<double> **** l_M_m_overlap_value , int **** l_M_m_index_l,
                                                  vector<complex<double>> *** l_M_m_local_overlap_value , vector<int> *** l_M_m_local_index_l  ,
                                                  int * recv_count, int * displs ){
    int l, m,  p; // index for state
    int i, j, k; // index for dof
    int local_eigenstate_num = selected_eigenstate_num / num_proc;
    if(my_id == num_proc -1){
        local_eigenstate_num = selected_eigenstate_num -  int(selected_eigenstate_num / num_proc) * (num_proc - 1);
    }

    int local_eigenstate_begin_index = my_id * int (selected_eigenstate_num / num_proc);

    int state_l_index ;
    int state_m_index;
    complex<double> A;
    complex<double> B;



    int m_phi_ladder_phi_list_length = 0;
    int p_phi_ladder_phi_list_length = 0;
    int index_m;
    int index_p;
    double Value1;
    double Value2;
    complex<double> Value;

    vector<phi_operator_phi_tuple> * ptr1;
    vector<phi_operator_phi_tuple> * ptr2;

    complex <double> * temporary_list = new complex<double> [eigenstate_num];
    double l_M_m_abs_cutoff = Eigenstate_OTOC_sift_criteria ;

    // compute l_M_m_local
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0];j++){
            for(m=0;m< local_eigenstate_num ; m++){
                for(l=0;l<eigenstate_num;l++){
                    temporary_list[l] = 0;
                }


                state_m_index = selected_eigenstate_index[ local_eigenstate_begin_index + m ];

                // first compute part : <l | a_{i}(t) a_{j} | m>
                ptr1 = & phi_ladder_operator_phi_tuple_list[j][state_m_index];  // List for a_{j}|m>
                m_phi_ladder_phi_list_length = (*ptr1) .size();
                for(index_m = 0; index_m < m_phi_ladder_phi_list_length; index_m ++){
                    Value1 = (*ptr1)[index_m].phi_operator_phi_value;  // <p| a_{j} | m>
                    p = (*ptr1)[index_m].eigenstate_index;  // index for eigenstate p;

                    ptr2 = & phi_ladder_operator_phi_tuple_list[i][p];  // a_{i}|p>
                    p_phi_ladder_phi_list_length = (*ptr2).size();
                    for (index_p=0; index_p < p_phi_ladder_phi_list_length ; index_p ++){
                        Value2 = (*ptr2)[index_p].phi_operator_phi_value;  // <l| a_{i} | p>
                        l = (*ptr2)[index_p].eigenstate_index;

                        Value =   complex<double> (Value1 * Value2) *
                                  std::exp(complex<double>(1j * cf2 *  time * ( Eigenvalue_list[l] - Eigenvalue_list[p] )  )  );



                        temporary_list[l] = temporary_list[l] + Value;
                    }

                }

                // second compute part : <l|a_{j} a_{i}(t) | m>
                ptr1 = & phi_ladder_operator_phi_tuple_list[i][state_m_index];
                m_phi_ladder_phi_list_length = (*ptr1).size();
                for(index_m = 0; index_m < m_phi_ladder_phi_list_length; index_m ++ ){
                    Value1 = (*ptr1)[index_m].phi_operator_phi_value; // <p|a_{i}|m>
                    p = (*ptr1)[index_m].eigenstate_index;

                    ptr2 = & phi_ladder_operator_phi_tuple_list[j][p]; // a_{j}|p>
                    p_phi_ladder_phi_list_length = (*ptr2).size();
                    for(index_p = 0; index_p < p_phi_ladder_phi_list_length ; index_p ++){
                        Value2 = (*ptr2)[index_p].phi_operator_phi_value; // <l | a_{j} | p >
                        l = (*ptr2)[index_p].eigenstate_index;

                        Value =  complex<double> (Value1 * Value2) *
                                 std::exp(complex<double> ( 1j * cf2 * time * (Eigenvalue_list[p] - Eigenvalue_list[state_m_index])  ) );


                        temporary_list[l] = temporary_list[l] - Value;
                    }

                }

                // Sift result in temporary_list and only store result that is larger than cutoff.
                for(l=0;l<eigenstate_num; l++){
                    if ( abs(temporary_list[l]) > l_M_m_abs_cutoff ){
                        l_M_m_local_index_l[i][j][m].push_back(l);
                        l_M_m_local_overlap_value[i][j][m].push_back(temporary_list[l]);
                    }
                }

            }

        }
    }


    int * l_M_m_list_size = new int [selected_eigenstate_num];
    int * l_M_m_list_size_local = new int [local_eigenstate_num];

    double * l_M_m_bcast_real ;
    double * l_M_m_bcast_imag ;

    int send_pc_id ;
    int local_eigenstate_index ;
    // transfer l_M_m_local to l_M_m
    for(i=0;i<2 * nmodes[0] ; i++){
        for(j=0;j<2 * nmodes[0]; j++){

            // send number of state in l_M_m_local_index_l and l_M_m_overlap_value
            for (m=0; m< local_eigenstate_num; m++ ){
                l_M_m_list_size_local[m] = l_M_m_local_index_l[i][j][m].size();
            }
            MPI_Allgatherv(&l_M_m_list_size_local[0] , local_eigenstate_num, MPI_INT,
                           & l_M_m_list_size[0], &recv_count[0], &displs[0], MPI_INT, MPI_COMM_WORLD );


            for( m = 0 ; m < selected_eigenstate_num ; m++ ){
                l_M_m_overlap_value[i][j][m] = new complex<double> [ l_M_m_list_size[m] ];
                l_M_m_index_l[i][j][m] = new int [ l_M_m_list_size[m] ];
            }

            for(m=0; m< selected_eigenstate_num ; m++ ){
                send_pc_id = m / int(selected_eigenstate_num / num_proc) ;
                if(send_pc_id  >= num_proc ){
                    send_pc_id = num_proc - 1 ;
                }

                local_eigenstate_index = m - send_pc_id * int(selected_eigenstate_num / num_proc) ;
                l_M_m_bcast_real = new double [ l_M_m_list_size[m] ];
                l_M_m_bcast_imag = new double [ l_M_m_list_size[m] ];

                if(send_pc_id == my_id ){
                    for( l=0; l<l_M_m_list_size[m] ; l ++ ){
                        l_M_m_bcast_real [l] = real (l_M_m_local_overlap_value[i][j][local_eigenstate_index][l] ) ;
                        l_M_m_bcast_imag [l] = imag (l_M_m_local_overlap_value[i][j][local_eigenstate_index][l] ) ;
                        l_M_m_index_l[i][j][m][l] = l_M_m_local_index_l [i][j][ local_eigenstate_index ][l];
                    }
                }

                // bcast l_M_m_overlap[i][j][m]
                MPI_Bcast(&l_M_m_bcast_real[0], l_M_m_list_size[m], MPI_DOUBLE, send_pc_id , MPI_COMM_WORLD);
                MPI_Bcast(&l_M_m_bcast_imag[0] , l_M_m_list_size[m] , MPI_DOUBLE, send_pc_id , MPI_COMM_WORLD);
                for(l=0; l< l_M_m_list_size[m]; l++ ){
                    l_M_m_overlap_value[i][j][m][l] = complex<double> ( l_M_m_bcast_real[l] , l_M_m_bcast_imag[l] );
                }
                // bacst l_M_m_index_l[i][j][m]
                MPI_Bcast(&l_M_m_index_l[i][j][m][0], l_M_m_list_size[m] , MPI_INT, send_pc_id ,MPI_COMM_WORLD);

                delete [] l_M_m_bcast_real ;
                delete [] l_M_m_bcast_imag ;

            }


        }
    }


    // sparsify result:
    // l_M_m_nonzero (size : [ 2 * nmodes[0], 2*nmodes[0] , eigenstate_num, list]
    vector< vector< vector< vector<phi_operator_phi_tuple_complex> > > > l_M_m_nonzero;
    for(i=0;i<2*nmodes[0];i++){
        vector< vector< vector<phi_operator_phi_tuple_complex> > >  l_M_m_nonzero_1;
        for(j=0;j<2*nmodes[0];j++){
            vector< vector<phi_operator_phi_tuple_complex> > l_M_m_nonzero_2;

            // compute size for each l_M_m_overlap_value[i][j][m]
            for (m=0; m< local_eigenstate_num; m++ ){
                l_M_m_list_size_local[m] = l_M_m_local_index_l[i][j][m].size();
            }
            MPI_Allgatherv(&l_M_m_list_size_local[0] , local_eigenstate_num, MPI_INT,
                           & l_M_m_list_size[0], &recv_count[0], &displs[0], MPI_INT, MPI_COMM_WORLD );


            for(m=0;m<selected_eigenstate_num;m++) {
                vector<phi_operator_phi_tuple_complex> l_M_m_nonzero_3;
                for(l=0;l< l_M_m_list_size[m] ;l++){
                    // record eigenstate index and eigenstate overlap for <l|[a_{i}, a_{j}] |m>
                    phi_operator_phi_tuple_complex Tuple1( l_M_m_index_l[i][j][m][l], l_M_m_overlap_value[i][j][m][l] );
                    l_M_m_nonzero_3.push_back(Tuple1);
                }
                l_M_m_nonzero_2.push_back(l_M_m_nonzero_3);
            }

            l_M_m_nonzero_1.push_back(l_M_m_nonzero_2);
        }
        l_M_m_nonzero.push_back(l_M_m_nonzero_1);
    }

    // compute Eigenstate_OTOC in each process and stored in local_Eigenstate_OTOC
    int a_i_a_j_m_index  = 0;
    int l_M_m_nonzero_m_list_length = 0;
    vector<phi_operator_phi_tuple_complex> * ptr3;
    vector<phi_operator_phi_tuple_complex> * ptr4;
    complex<double> Value3;
    complex<double> Value4;
    complex<double> local_l_M_m_sum ;
    int position ;
    for(m=0;m<local_eigenstate_num;m++){
        for(i=0;i<2 * nmodes[0] ; i++){
            for(j=0; j<2 * nmodes[0] ; j++){
                local_l_M_m_sum = 0;

                state_m_index = local_eigenstate_begin_index + m  ;

                for(k=0;k<2 * nmodes[0] ; k++){
                    ptr3 = & l_M_m_nonzero[k][i][state_m_index];   // l_M_m[k][i][:][state_m]
                    l_M_m_nonzero_m_list_length = (*ptr3).size();
                    ptr4 = & l_M_m_nonzero[k][j][state_m_index]; // l_M_m[k][j][:][state_m]

                    // loop for a_i_a_j_m_index is equivalent to loop for state l
                    for(a_i_a_j_m_index = 0; a_i_a_j_m_index < l_M_m_nonzero_m_list_length ; a_i_a_j_m_index ++ ){
                        l = (*ptr3)[a_i_a_j_m_index].eigenstate_index;
                        Value3 = (*ptr3)[a_i_a_j_m_index].phi_operator_phi_value;  // l_M_m[k][i][l][state_m]

                        position = Binary_search_phi_operator_phi_tuple_complex((*ptr4) , l );
                        if(position != -1){
                            // state l have non-negligible contribution
                            Value4 = (*ptr4)[position].phi_operator_phi_value; // l_M_m[k][j][l][state_m]
                            local_l_M_m_sum = local_l_M_m_sum +  conj(Value3) * Value4;

                        }

                    }

                }

                local_Eigenstate_OTOC[i][j][m] = real(local_l_M_m_sum);
            }
        }
    }

    // Gather local_Eigenstate_OTOC to Eigenstate_OTOC
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0]; j++){
            MPI_Gatherv(&local_Eigenstate_OTOC[i][j][0], local_eigenstate_num, MPI_DOUBLE,
                        &Eigenstate_OTOC[i][j][0], &recv_count[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }

    // output to Eigenstate_OTOC_output
    if(my_id == 0){
        // output time
        Eigenstate_OTOC_output << time << endl;

        // output format:  each state, output real and imaginary part. So in total 2 *  eigenstate_num lines
        for(m=0;m<selected_eigenstate_num;m++){
            for(i=0;i< 2*nmodes[0] ; i++){
                for(j=0;j<2*nmodes[0];j++){
                    Eigenstate_OTOC_output << Eigenstate_OTOC[i][j][m] <<" ";
                }
            }
            Eigenstate_OTOC_output << endl;
        }

    }

    // clear the vector in l_M_m_local_***
    for(i = 0; i < 2 * nmodes[0] ; i++){
        for(j = 0 ; j < 2 * nmodes[0] ; j++) {
            for( m = 0 ; m < local_eigenstate_num ; m++ ){
                l_M_m_local_overlap_value[i][j][m].clear();
                l_M_m_local_index_l[i][j][m].clear();
            }
        }
    }
    // free space for l_M_m_index_l and l_M_m_overlap_value
    for(i = 0; i < 2 * nmodes[0] ; i++){
        for(j = 0 ; j < 2 * nmodes[0] ; j++) {
            for( m = 0 ; m < selected_eigenstate_num ; m++ ){
                delete [] l_M_m_overlap_value[i][j][m];
                delete [] l_M_m_index_l [i][j][m];
            }
        }
    }


    delete [] temporary_list ;
    delete [] l_M_m_list_size;
    delete [] l_M_m_list_size_local ;
}

void detector::compute_Eigenstate_OTOC(){
    int i,j,k;
    vector<double> Time_series ;
    double t = 0;

    while(t <= proptime[0]){
        Time_series.push_back(t);
        t = t + tprint;
    }
    int Time_series_len = Time_series.size();

    int local_eigenstate_num = selected_eigenstate_num / num_proc;
    if(my_id == num_proc -1){
        local_eigenstate_num = selected_eigenstate_num -  int(selected_eigenstate_num / num_proc) * (num_proc - 1);
    }


    // size: [eigenstate_num , 2* nmodes[0], 2* nmodes[0]] store final result
    double *** Eigenstate_OTOC;

    Eigenstate_OTOC = new double ** [2 * nmodes[0]];
    for(i=0;i<2 * nmodes[0];i++){
        Eigenstate_OTOC[i] = new double * [ 2 * nmodes[0]];
        for(j=0;j<2*nmodes[0];j++){
            Eigenstate_OTOC[i][j] = new double [selected_eigenstate_num];
        }
    }



    double *** local_Eigenstate_OTOC;
    local_Eigenstate_OTOC = new double ** [2*nmodes[0]];
    for(i=0;i<2*nmodes[0];i++){
        local_Eigenstate_OTOC[i] = new  double * [ 2 * nmodes[0]];
        for(j=0;j<2*nmodes[0]; j++){
            local_Eigenstate_OTOC[i][j] = new  double  [local_eigenstate_num];
        }
    }


    // size [ 2*nmodes[0] , 2*nmodes[0], eigenstate_num, eigenstate_num]   : <phi_l | [a_{i}(t) , a_{j}] | phi_m>
    // third dimension will be partitioned to let different process compute each part. [i, j ,m ,l ]
    // last dimension will only be allocated in submodule when we compute sparse <l| [a_{i}(t) , a_{j}]  | m> at time t.
    complex<double>  **** l_M_m_overlap_value = new complex<double> *** [2*nmodes[0]];
    int **** l_M_m_index_l = new int *** [ 2 * nmodes[0] ];
    for(i=0;i< 2* nmodes[0]; i++){
        l_M_m_overlap_value[i] = new complex<double> ** [2 * nmodes[0]];
        l_M_m_index_l[i] = new int ** [2 * nmodes[0] ];
        for(j=0;j<2*nmodes[0];j++){
            l_M_m_overlap_value[i][j] = new complex<double> * [selected_eigenstate_num];
            l_M_m_index_l[i][j] = new int * [ selected_eigenstate_num ];
        }
    }

    vector<complex<double> > *** l_M_m_local_overlap_value = new vector<complex<double>> ** [2*nmodes[0]];
    vector<int> *** l_M_m_local_index_l = new vector<int> ** [2 * nmodes[0]];
    for(i=0;i< 2* nmodes[0]; i++){
        l_M_m_local_overlap_value[i] = new vector<complex<double>> * [2 * nmodes[0]];
        l_M_m_local_index_l[i] = new vector<int> * [ 2 * nmodes[0] ];
        for(j=0;j<2*nmodes[0];j++){
            l_M_m_local_overlap_value[i][j] = new vector<complex<double>>  [local_eigenstate_num];
            l_M_m_local_index_l[i][j] = new vector<int>  [ local_eigenstate_num ] ;
        }
    }

    int * recv_count = new int [num_proc];
    int * displs = new int [num_proc];
    for(i=0;i<num_proc-1;i++){
        recv_count [i] = int(selected_eigenstate_num / num_proc) ;
    }
    recv_count[num_proc - 1] = selected_eigenstate_num -  int(selected_eigenstate_num / num_proc) * (num_proc - 1);
    displs[0] = 0;
    for(i=1;i<num_proc;i++){
        displs[i] = displs[i - 1] + recv_count[ i - 1 ];
    }

    ofstream Eigenstate_OTOC_output;
    if(my_id == 0){
        Eigenstate_OTOC_output.open(path + "Eigenstate_OTOC.txt");
        //output information about dof
        Eigenstate_OTOC_output << nmodes[0] << endl;

        // output information about eigenstate (eigenvalue)
        Eigenstate_OTOC_output << selected_eigenstate_num << endl;

        for(i=0;i<selected_eigenstate_num;i++){
            Eigenstate_OTOC_output << Eigenvalue_list[ selected_eigenstate_index[i] ] <<"  ";
        }
        Eigenstate_OTOC_output<<endl;

        // output totoal number of time step to output
        Eigenstate_OTOC_output << Time_series_len << endl;

    }

    // compute Eigenstate OTOC and output
    for(i=0;i<Time_series_len;i++){
        t = Time_series[i];
        compute_Eigenstate_OTOC_submodule(Eigenstate_OTOC_output,t,Eigenstate_OTOC, local_Eigenstate_OTOC, l_M_m_overlap_value, l_M_m_index_l, l_M_m_local_overlap_value  ,  l_M_m_local_index_l , recv_count,displs);
        if(my_id == 0){
            cout << " finish computing Eigenstate OTOC for time:  t =  " << t << endl;
        }
    }

    // free space:
    for(i=0;i< 2 * nmodes[0] ;i++){
        for(j=0;j<2*nmodes[0];j++){
            delete [] Eigenstate_OTOC[i][j];
        }
        delete [] Eigenstate_OTOC[i];
    }
    delete [ ] Eigenstate_OTOC;


    for(i=0;i< 2 * nmodes[0] ;i++){
        for(j=0;j<2*nmodes[0];j++){
            delete [] local_Eigenstate_OTOC[i][j];
        }
        delete [] local_Eigenstate_OTOC[i];
    }
    delete [ ] local_Eigenstate_OTOC;


    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0];j++){
            delete [] l_M_m_overlap_value[i][j];
            delete [] l_M_m_index_l[i][j];

            delete [] l_M_m_local_overlap_value[i][j];
            delete [] l_M_m_local_index_l[i][j];
        }
        delete [] l_M_m_overlap_value[i];
        delete [] l_M_m_index_l[i];

        delete [] l_M_m_local_overlap_value[i];
        delete [] l_M_m_local_index_l[i];
    }
    delete [] l_M_m_overlap_value;
    delete [] l_M_m_index_l ;

    delete [] l_M_m_local_overlap_value;
    delete [] l_M_m_local_index_l;

    delete [] recv_count;
    delete [] displs;

    // close files
    if(my_id == 0){
        Eigenstate_OTOC_output.close();
    }
}