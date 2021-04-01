//
// Created by phyzch on 3/31/21.
//
#include "util.h"
#include "system.h"
using namespace std;

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
    int i , j;
    Eigenstate_energy_std_list = new double [eigenstate_num];
    for(i = 0; i<eigenstate_num ; i++){
        Eigenstate_energy_std_list[i] = 0;
        // \sum (E_{n} - eigenvalue)^2 * |<n|\phi>|^2
        for(j=0;j<total_dmat_size[0];j++){
            Eigenstate_energy_std_list[i] = Eigenstate_energy_std_list[i] +
                    pow(Eigenvalue_list[i] - total_dmat[0][j] , 2) * pow(Eigenstate_list[i][j] , 2);
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

                if(energy_difference > 2 * (Eigenstate_energy_std_list[eigenstate_m_index] + Eigenstate_energy_std_list[eigenstate_l_index] )){
                    local_phi_ladder_operator_phi[i][m][l] = 0;
                }
                else{
                    // we need to compute <\phi_m | a_{i} | phi_l> or <\phi_m | a_{i}^{+} | \phi_l>
                    local_ladder_operator_value = 0;
                    for(j=0;j<total_dmat_size[0]; j++){
                        basis_set_state_index = j;
                        // for i < nmodes[0] , that's state_mode[i] -1. which corresponds to lowering operator.
                        // for i > nmodes[0], that's state_mode[i] + 1, corresponds to raising operator
                        nearby_basis_set_state_index = neighbor_state_index_for_all_state_list[j][i];
                        if(nearby_basis_set_state_index == -1){
                            continue;
                        }
                        else{
                            if(i<nmodes[0]){
                                // lowering operator
                                local_ladder_operator_value = local_ladder_operator_value + Eigenstate_list[eigenstate_l_index][basis_set_state_index] *
                                                    Eigenstate_list[eigenstate_m_index][ nearby_basis_set_state_index ] * sqrt(dv_all[0][basis_set_state_index][i]);
                            }
                            else{
                                // raising operator
                                local_ladder_operator_value = local_ladder_operator_value + Eigenstate_list[eigenstate_l_index][basis_set_state_index] *
                                            Eigenstate_list[eigenstate_m_index][ nearby_basis_set_state_index ] * sqrt(dv_all[0][basis_set_state_index][i - nmodes[0] ] + 1);
                            }

                        }

                    }
                    local_phi_ladder_operator_phi[i][m][l]  = local_ladder_operator_value;

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

    // free space
    for(i=0;i<2*nmodes[0];i++){
        delete [] phi_ladder_operator_phi_1d[i];
    }
    delete [] phi_ladder_operator_phi_1d;


}

void detector:: compute_Eigenstate_OTOC_submodule(ofstream & Eigenstate_OTOC_output, double time, complex<double> *** Eigenstate_OTOC ,
                                                  complex<double> **** l_M_m , complex<double> **** l_M_m_local , int * recv_count, int * displs ){
    int l, m,  p; // index for state
    int i, j, k; // index for dof
    int local_eigenstate_num = eigenstate_num / num_proc;
    if(my_id == num_proc -1){
        local_eigenstate_num = eigenstate_num -  int(eigenstate_num / num_proc) * (num_proc - 1);
    }

    int local_eigenstate_begin_index = my_id * int (eigenstate_num / num_proc);

    int state_l_index ;
    int state_m_index;
    complex<double> A;
    complex<double> B;

    double * send_local_l_M_m_real = new double [local_eigenstate_num];
    double * send_local_l_M_m_imag = new double [local_eigenstate_num];

    double * recv_l_M_m_real = new double [eigenstate_num];
    double * recv_l_M_m_imag = new double [eigenstate_num];

    // compute l_M_m_local
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0];j++){
            for(l=0;l<eigenstate_num;l++){
                for(m=0;m<local_eigenstate_num;m++){
                    state_m_index = local_eigenstate_begin_index + m;
                    state_l_index = l;
                    // compute <l [a_{i}(t) , a_{j}] | m>
                    // phi_ladder_operator_phi = <l | a_{i} | p>
                    l_M_m_local[i][j][l][m] = 0;
                    for(p=0;p<eigenstate_num;p++){
                        if(phi_ladder_operator_phi[i][state_l_index][p] != 0 and phi_ladder_operator_phi[j][p][state_m_index] != 0 ){
                            A = std::exp(complex<double>(1j * cf2 *( Eigenvalue_list[state_l_index] - Eigenvalue_list[p] ) * time )  )
                                   *  complex<double>( phi_ladder_operator_phi[i][state_l_index][p] * phi_ladder_operator_phi[j][p][state_m_index]  );
                        }
                        else{
                            A = 0 ;
                        }
                        if( phi_ladder_operator_phi[j][state_l_index][p] != 0 and phi_ladder_operator_phi[i][p][state_m_index] != 0 ){
                            B = std::exp( complex<double>(  1j * cf2 * time * (Eigenvalue_list[p] - Eigenvalue_list[state_m_index])  ) )
                                    * complex<double> ( phi_ladder_operator_phi[j][state_l_index][p] * phi_ladder_operator_phi[i][p][state_m_index]  );
                        }
                        else{
                            B = 0;
                        }

                        l_M_m_local[i][j][l][m] = l_M_m_local[i][j][l][m] +  (A - B );
                    }

                }
            }
        }
    }

    if(my_id == 0){
        cout <<"Finsh computing  l M m local " << endl;
    }

    // transfer l_M_m_local to l_M_m
    for(i=0;i<2 * nmodes[0] ; i++){
        for(j=0;j<2 * nmodes[0]; j++){
            for(l=0;l<eigenstate_num;l++){

                for(m=0;m<local_eigenstate_num;m++){
                    send_local_l_M_m_real[m] = real(l_M_m_local[i][j][l][m]);
                    send_local_l_M_m_imag [m] = imag(l_M_m_local[i][j][l][m]);
                }
                MPI_Gatherv(&send_local_l_M_m_real[0], local_eigenstate_num, MPI_DOUBLE,
                            & recv_l_M_m_real[0], &recv_count[0],&displs[0],MPI_DOUBLE, 0 ,MPI_COMM_WORLD);

                MPI_Bcast(&recv_l_M_m_real[0], eigenstate_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                MPI_Gatherv(&send_local_l_M_m_imag[0], local_eigenstate_num, MPI_DOUBLE,
                            &recv_l_M_m_imag[0], &recv_count[0], &displs[0], MPI_DOUBLE, 0, MPI_COMM_WORLD );

                MPI_Bcast(&recv_l_M_m_imag[0], eigenstate_num , MPI_DOUBLE, 0 , MPI_COMM_WORLD);

                for(m=0;m<eigenstate_num;m++){
                    l_M_m[i][j][l][m] = complex<double> ( recv_l_M_m_real[m] , recv_l_M_m_imag[m] );
                }

            }
        }
    }

    if(my_id == 0){
        cout <<"Finsh transfer l M m local to l_M_m" << endl;
    }

    complex<double> local_l_M_m_sum ;
    double local_l_M_m_sum_real;
    double local_l_M_m_sum_imag;

    double l_M_m_sum_real;
    double l_M_m_sum_imag;

    for (m=0;m<eigenstate_num;m++){
        for(i=0;i<2 * nmodes[0] ; i++){
            for(j=0; j<2 * nmodes[0] ; j++){
                local_l_M_m_sum = 0;
                // (M^{\dagger})_{ik}^{ml} *M_{kj}^{lm} here lm is index for state. i,j is index for dof
                // Here we first sum over state k in each process, then MPI_sum over state l.
                for(l=0;l<local_eigenstate_num;l++){
                        state_l_index = l + local_eigenstate_begin_index;
                        // sum over dof k
                        for(k=0 ; k< 2* nmodes[0] ; k++){
                            local_l_M_m_sum = local_l_M_m_sum + conj(l_M_m[k][i][state_l_index][m])
                                    * l_M_m[k][j][state_l_index][m];
                        }
                }

                local_l_M_m_sum_real = real(local_l_M_m_sum);
                local_l_M_m_sum_imag = imag(local_l_M_m_sum);

                MPI_Reduce(&local_l_M_m_sum_real, &l_M_m_sum_real, num_proc, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
                MPI_Reduce(&local_l_M_m_sum_imag , &l_M_m_sum_imag, num_proc, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

                if(my_id == 0){
                    Eigenstate_OTOC[i][j][m] = complex<double> ( l_M_m_sum_real , l_M_m_sum_imag);
                }

            }
        }
    }

    // output to Eigenstate_OTOC_output
    if(my_id == 0){
        // output time
        Eigenstate_OTOC_output << time << endl;

        // output format:  each state, output real and imaginary part. So in total 2 *  eigenstate_num lines
        for(m=0;m<eigenstate_num;m++){
            for(i=0;i< 2*nmodes[0] ; i++){
                for(j=0;j<2*nmodes[0];j++){
                    Eigenstate_OTOC_output << real(Eigenstate_OTOC[i][j][m]) <<" ";
                }
            }
            Eigenstate_OTOC_output << endl;
            for(i=0;i< 2*nmodes[0] ; i++){
                for(j=0;j<2*nmodes[0];j++){
                    Eigenstate_OTOC_output << imag(Eigenstate_OTOC[i][j][m]) <<" ";
                }
            }
            Eigenstate_OTOC_output << endl;
        }

    }


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

    int local_eigenstate_num = eigenstate_num / num_proc;
    if(my_id == num_proc -1){
        local_eigenstate_num = eigenstate_num -  int(eigenstate_num / num_proc) * (num_proc - 1);
    }


    // size: [eigenstate_num , 2* nmodes[0], 2* nmodes[0]] store final result
    complex<double> *** Eigenstate_OTOC;
    if(my_id == 0){

        Eigenstate_OTOC = new complex<double> ** [2 * nmodes[0]];
        for(i=0;i<2 * nmodes[0];i++){
            Eigenstate_OTOC[i] = new complex<double> * [ 2*nmodes[0]];
            for(j=0;j<2*nmodes[0];j++){
                Eigenstate_OTOC[i][j] = new complex<double> [eigenstate_num];
            }
        }

    }


    // size [ 2*nmodes[0] , 2*nmodes[0], eigenstate_num, eigenstate_num]   : <phi_l | [a_{i}(t) , a_{j}] | phi_m>
    // last dimension will be partitioned to let different process compute each part.
    complex<double> **** l_M_m = new complex<double> *** [2*nmodes[0]];
    for(i=0;i< 2* nmodes[0]; i++){
        l_M_m[i] = new complex<double> ** [2 * nmodes[0]];
        for(j=0;j<2*nmodes[0];j++){
            l_M_m[i][j] = new complex<double> * [eigenstate_num];
            for(k=0;k<eigenstate_num;k++){
                l_M_m[i][j][k] = new complex<double> [eigenstate_num];
            }
        }
    }

    complex<double> **** l_M_m_local = new complex<double> *** [2*nmodes[0]];
    for(i=0;i< 2* nmodes[0]; i++){
        l_M_m_local[i] = new complex<double> ** [2 * nmodes[0]];
        for(j=0;j<2*nmodes[0];j++){
            l_M_m_local[i][j] = new complex<double> * [eigenstate_num];
            for(k=0;k<eigenstate_num;k++){
                l_M_m_local[i][j][k] = new complex<double> [local_eigenstate_num];
            }
        }
    }

    int * recv_count = new int [num_proc];
    int * displs = new int [num_proc];
    for(i=0;i<num_proc-1;i++){
        recv_count [i] = int(eigenstate_num / num_proc) ;
    }
    recv_count[num_proc - 1] = int(eigenstate_num - int(eigenstate_num / num_proc) * (num_proc - 1));
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
        Eigenstate_OTOC_output << eigenstate_num << endl;

        for(i=0;i<eigenstate_num;i++){
            Eigenstate_OTOC_output << Eigenvalue_list[i] <<"  ";
        }
        Eigenstate_OTOC_output<<endl;

        // output totoal number of time step to output
        Eigenstate_OTOC_output << Time_series_len << endl;

    }

    if(my_id == 0){
        cout <<"Begin computing Eigenstate OTOC " << endl;
    }

    // compute Eigenstate OTOC and output
    for(i=0;i<Time_series_len;i++){
        t = Time_series[i];
        compute_Eigenstate_OTOC_submodule(Eigenstate_OTOC_output,t,Eigenstate_OTOC, l_M_m, l_M_m_local,recv_count,displs);
    }

    if(my_id == 0){
        cout <<"Finish computing Eigenstate OTOC " << endl;
    }


    // free space:
    if(my_id == 0){
        for(i=0;i< 2 * nmodes[0] ;i++){
            for(j=0;j<2*nmodes[0];j++){
                delete [] Eigenstate_OTOC[i][j];
            }
            delete [] Eigenstate_OTOC[i];
        }
        delete [ ] Eigenstate_OTOC;


    }

    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0];j++){
            for(k=0;k<eigenstate_num;k++){
                delete [] l_M_m[i][j][k];
                delete [] l_M_m_local[i][j][k];
            }
            delete [] l_M_m[i][j];
            delete [] l_M_m_local[i][j];
        }
        delete [] l_M_m[i];
        delete [] l_M_m_local[i];
    }
    delete [] l_M_m;
    delete [] l_M_m_local;

    delete [] recv_count;
    delete [] displs;

    // close files
    if(my_id == 0){
        Eigenstate_OTOC_output.close();
    }
}