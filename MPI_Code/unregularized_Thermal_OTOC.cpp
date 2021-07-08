//
// Created by phyzch on 7/8/21.
//
#include"../util.h"
#include"../system.h"

void detector::allocate_space_for_unregularized_Lyapunov_spectrum_calculation(  ){
    int i, j, k, l;
    int nearby_state_basis_size = nearby_state_index.size();

    // regularized_thermal_Lyapunov_spectrum
    for(i=0; i<2*nmodes[0]; i++ ){
        vector<complex<double>> v1 (2 * nmodes[0] , 0);
        unregularized_thermal_Lyapunov_spectrum.push_back(v1);
    }

    // regularized_thermal_Lyapunov_spectrum_each_Haar_state  :   [N_Haar] [ 2*nmodes[0] ] [2 * nmodes[0] ]
    for(i=0;i<N_Haar; i++ ){
        vector<vector<complex<double> >> v1 ;
        for(j=0;j < 2* nmodes[0] ; j++ ){
            vector<complex<double>> v2 (2 * nmodes[0] , 0);
            v1.push_back(v2);
        }
        unregularized_thermal_Lyapunov_spectrum_each_Haar_state.push_back(v1);
    }

    // unregularized_thermal_OTOC_overlap_Haar_state_basis_set  :  [N_Haar] [N_basis_set] [2 * nmodes[0] ] [2 * nmodes[0] ]
    for(i = 0; i<N_Haar; i++ ){
        vector<vector< vector< complex<double> > > > v1 ;
        for( j = 0 ; j < nearby_state_basis_size; j++ ){
            vector<vector< complex<double> > > v2 ;
            for(k = 0; k< 2 * nmodes[0]; k++ ){
                vector<complex<double>> v3 ( 2 * nmodes[0] , 0 );
                v2.push_back(v3);
            }
            v1.push_back(v2);
        }
        unregularized_thermal_OTOC_overlap_Haar_state_basis_set.push_back(v1);
    }

    // Haar_state_with_ladder_operator_x_sparsify , y_sparsify , basis_set_sparsify
    unregularized_Haar_state_with_ladder_operator_x_sparsify = new vector<vector<double>> * [N_Haar];
    unregularized_Haar_state_with_ladder_operator_y_sparsify = new vector<vector<double>> * [N_Haar] ;
    unregularized_Haar_state_with_ladder_operator_basis_set_sparsify =  new vector<vector<int>> * [N_Haar];
    for(i=0;i<N_Haar ; i++ ){
        unregularized_Haar_state_with_ladder_operator_x_sparsify[i] = new vector<vector<double>> [2 * nmodes[0] + 1 ];
        unregularized_Haar_state_with_ladder_operator_y_sparsify[i] = new vector<vector<double>> [2 * nmodes[0] + 1 ];
        unregularized_Haar_state_with_ladder_operator_basis_set_sparsify[i] = new vector<vector<int>> [2 * nmodes[0] + 1 ];
    }
}

void detector::compute_unregularized_Haar_random_state_with_ladder_operator (double sparsify_criteria ){
    // compute  a_{i} e^{-iHt} a_{j} e^{-\beta H/2} |\phi_{Haar}>   and   a_{i}e^{-iHt} e^{-\beta H/2} |\phi_{Haar}>
    // result store in array with size [N_{Haar}] [2 * nmodes[0] + 1 ] [2 * nmodes[0] ] [vector<double (int) >]
    // Here 2*nmodes[0]+1 (2nd dimension) stands for 1 or a_{j} in above expression.  2 * nmodes[0]  (3rd dimension) stands for a_{i}
    // Haar_state related variable should have been allocated space before calling this function.
    int i, j, k , l;
    int state_index;
    double normalization;
    double magnitude;
    double total_normalization;
    for(i=0;i<N_Haar;i++){
        for(j=0;j<2*nmodes[0] + 1 ;j++){
            unregularized_Haar_state_with_ladder_operator_x_sparsify[i][j].clear();
            unregularized_Haar_state_with_ladder_operator_y_sparsify[i][j].clear();
            unregularized_Haar_state_with_ladder_operator_basis_set_sparsify[i][j].clear();
        }
    }


    for(i=0;i<N_Haar ; i++){
        for(j=0;j<2*nmodes[0]+1; j++ ){
            state_index = unregularized_Haar_state_index_list[i] + j ;
            vector<vector<double>> xd_for_ladder_operation;
            vector<vector<double>> yd_for_ladder_operation;
            ladder_operator_operation(xd[state_index] , yd[state_index] , xd_for_ladder_operation, yd_for_ladder_operation);
            for(l=0;l<2*nmodes[0];l++){
                // after operation of a_{l} ladder operator
                vector<double> vx;
                vector<double> vy;
                vector<int> v_basis_set;

                // sparsify xd(yd)_for_ladder_operation[l]
                normalization = 0;
                for(k=0;k< dmatsize[0]; k++ ){
                    normalization = normalization + norm( xd_for_ladder_operation[l][k]) + norm(yd_for_ladder_operation[l][k]);
                }
                MPI_Allreduce(&normalization, &total_normalization, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
                normalization = total_normalization;
                normalization = sqrt(normalization);

                for(k=0;k<dmatsize[0];k++){
                    magnitude = sqrt ( norm(xd_for_ladder_operation[l][k]) + norm(yd_for_ladder_operation[l][k])  );
                    if(magnitude > sparsify_criteria * normalization / total_dmat_size[0] ){
                        vx.push_back(xd_for_ladder_operation[l][k]);
                        vy.push_back(yd_for_ladder_operation[l][k]);
                        v_basis_set.push_back(k);
                    }
                }


                unregularized_Haar_state_with_ladder_operator_x_sparsify[i][j].push_back(vx);
                unregularized_Haar_state_with_ladder_operator_y_sparsify[i][j].push_back(vy);
                unregularized_Haar_state_with_ladder_operator_basis_set_sparsify[i][j].push_back(v_basis_set);

            }

        }

    }

}

void detector:: compute_unregularized_thermal_OTOC_component(){
    // compute <m| [a_{i}(t) a_{j}] e^{-\beta H / 2 } | Haar>
    // result store in unregularized_thermal_OTOC_overlap_Haar_state_basis_set . dimension [N_Haar][N_basis_set ] [2 * N_dof] [2* N_dof]

    int i, j, k, l, m;

    int nearby_state_basis_size = nearby_state_index.size();

    int begin_index = total_dmat_size[0] / num_proc * my_id ;

    int list_size;

    int basis_state_index;

    int m_state_index_after_ladder_operator;
    int mode_index;
    double Coeff;

    double real_part_value_Haar_state;
    double imag_part_value_Haar_state;

    double real_part_value_basis_state;
    double imag_part_value_basis_state;

    complex<double> C1;
    complex<double> C2;

    double * overlap_real_list = new double [nearby_state_basis_size];
    double * overlap_imag_list = new double [nearby_state_basis_size];
    double * overlap_real_list_sum = new double [nearby_state_basis_size];
    double * overlap_imag_list_sum = new double [nearby_state_basis_size];

    compute_unregularized_Haar_random_state_with_ladder_operator();

    for(i=0;i<N_Haar; i++ ){
        for(j=0;j<2*nmodes[0]; j++ ){
            for(k=0;k<2*nmodes[0];k++){

                for(m=0; m<nearby_state_basis_size; m++ ){
                    // compute <m| [a_{j}(t) , a_{k}] y^2 | Haar>

                    // compute <m| a_{j}(t) a_{k} y^2 | Haar> . result store in C1
                    C1 = 0;
                    list_size = unregularized_Haar_state_with_ladder_operator_basis_set_sparsify[i][k+1][j].size();
                    for(l=0;l<list_size;l++){
                        basis_state_index = unregularized_Haar_state_with_ladder_operator_basis_set_sparsify[i][k+1][j][l];
                        real_part_value_Haar_state = unregularized_Haar_state_with_ladder_operator_x_sparsify[i][k+1][j][l];
                        imag_part_value_Haar_state = unregularized_Haar_state_with_ladder_operator_y_sparsify[i][k+1][j][l];

                        real_part_value_basis_state = xd[m][basis_state_index] ;
                        imag_part_value_basis_state = yd[m][basis_state_index] ;

                        C1 = C1 + complex<double> (real_part_value_basis_state, - imag_part_value_basis_state) *
                                complex<double> (real_part_value_Haar_state, imag_part_value_Haar_state) ;

                    }

                    // compute <m| a_{k} a_{j}(t) y^2 | Haar> . result store in C2
                    C2 = 0;
                    list_size = unregularized_Haar_state_with_ladder_operator_basis_set_sparsify[i][0][j].size();
                    if(k<nmodes[0]){
                        mode_index = k + nmodes[0];  // <m| a_{k} == conj (a_{k}^{+} |m>)
                        Coeff = sqrt(dv_all[0][m][k] + 1) ;  // \sqrt{n_{k} + 1 }
                    }
                    else{
                        mode_index = k - nmodes[0]; // <m|a_{k}^{+} == conj ( a_{k} |m> )
                        Coeff = sqrt(dv_all[0][m][k - nmodes[0]]) ;  // sqrt{n_{k}}
                    }
                    m_state_index_after_ladder_operator = neighbor_state_index_for_all_state_list[m][mode_index];

                    if(m_state_index_after_ladder_operator != -1){
                        for(l=0; l<list_size; l++ ){
                            basis_state_index = unregularized_Haar_state_with_ladder_operator_basis_set_sparsify[i][0][j][l];
                            real_part_value_Haar_state = unregularized_Haar_state_with_ladder_operator_x_sparsify[i][0][j][l];
                            imag_part_value_Haar_state = unregularized_Haar_state_with_ladder_operator_y_sparsify[i][0][j][l];

                            real_part_value_basis_state = xd[m_state_index_after_ladder_operator][basis_state_index];
                            imag_part_value_basis_state = yd [m_state_index_after_ladder_operator] [basis_state_index];

                            C2 = C2 + Coeff * complex<double> (real_part_value_basis_state, - imag_part_value_basis_state) *
                                              complex<double> (real_part_value_Haar_state, imag_part_value_Haar_state) ;

                        }
                    }

                    // sum result in different process
                    overlap_real_list[m] = real(C1 - C2 );
                    overlap_imag_list[m] = imag(C1- C2 );

                }

                MPI_Allreduce(&overlap_real_list[0] , &overlap_real_list_sum[0] , nearby_state_basis_size, MPI_DOUBLE, MPI_SUM,MPI_COMM_WORLD);
                MPI_Allreduce(&overlap_imag_list[0] , &overlap_imag_list_sum[0] , nearby_state_basis_size , MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);
                for(m=0;m<nearby_state_basis_size; m++ ){
                    unregularized_thermal_OTOC_overlap_Haar_state_basis_set[i][m][j][k] = complex<double> (overlap_real_list_sum[m] , overlap_imag_list_sum[m]) ;
                }

            }
        }
    }


    delete [] overlap_real_list;
    delete [] overlap_real_list_sum;
    delete [] overlap_imag_list;
    delete [] overlap_imag_list_sum;

}

void  detector:: compute_unregularized_thermal_OTOC_Lyapunov_spectrum(  ){
    int i,j, k, l, m;
    complex<double> Lyapunov_spectrum_component;
    int nearby_state_basis_size = nearby_state_index.size();
    compute_unregularized_thermal_OTOC_component();

    int Haar_index;


    for( Haar_index = 0; Haar_index < N_Haar; Haar_index ++ ){
        for(i=0;i<2 * nmodes[0]; i++){
            for(j = 0; j< 2*nmodes[0] ; j++ ){

                // compute L_{ij} for Haar state |\phi>
                Lyapunov_spectrum_component = 0 ;
                for(m=0;m< nearby_state_basis_size ; m++ ){
                    for(k=0; k < 2*nmodes[0]; k++ ){
                        Lyapunov_spectrum_component = Lyapunov_spectrum_component +
                                                      std::conj(unregularized_thermal_OTOC_overlap_Haar_state_basis_set[Haar_index][m][k][i]) *
                                                      unregularized_thermal_OTOC_overlap_Haar_state_basis_set[Haar_index][m][k][j] ;
                    }
                }

                unregularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j] = Lyapunov_spectrum_component ;
                unregularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j] =
                        unregularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j] / Haar_state_normalization_list[Haar_index];

            }
        }
    }

    // average result
    for(i=0;i<2*nmodes[0]; i++ ){
        for(j=0;j<2*nmodes[0]; j++){
            unregularized_thermal_Lyapunov_spectrum[i][j] = 0;
            for(Haar_index = 0; Haar_index < N_Haar; Haar_index ++ ){
                unregularized_thermal_Lyapunov_spectrum[i][j] = unregularized_thermal_Lyapunov_spectrum[i][j] +
                                                              unregularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j];
            }
            unregularized_thermal_Lyapunov_spectrum[i][j] = unregularized_thermal_Lyapunov_spectrum[i][j] / double(N_Haar);

        }
    }



}
