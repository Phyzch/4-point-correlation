//
// Created by phyzch on 7/7/21.
//
#include "../util.h"
#include "../system.h"

void detector::allocate_space_for_Haar_state_calculation(  ){
    int i, j, k, l;
    int nearby_state_basis_size = nearby_state_index.size();

    // regularized_thermal_Lyapunov_spectrum
    for(i=0; i<2*nmodes[0]; i++ ){
        vector<complex<double>> v1 (2 * nmodes[0] , 0);
        regularized_thermal_Lyapunov_spectrum.push_back(v1);
    }

     // regularized_thermal_Lyapunov_spectrum_each_Haar_state  :   [N_Haar] [ 2*nmodes[0] ] [2 * nmodes[0] ]
     for(i=0;i<N_Haar; i++ ){
         vector<vector<complex<double> >> v1 ;
         for(j=0;j < 2* nmodes[0] ; j++ ){
             vector<complex<double>> v2 (2 * nmodes[0] , 0);
             v1.push_back(v2);
         }
         regularized_thermal_Lyapunov_spectrum_each_Haar_state.push_back(v1);
     }

   // regularized_thermal_OTOC_overlap_Haar_state_basis_set  :  [N_Haar] [N_basis_set] [2 * nmodes[0] ] [2 * nmodes[0] ]
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
       regularized_thermal_OTOC_overlap_Haar_state_basis_set.push_back(v1);
   }

   // Haar_state_overlap_time_dependent_basis_set
    for(i = 0; i<N_Haar; i++ ){
        vector<vector< vector< complex<double> > > > v1 ;
        for( j = 0 ; j < nearby_state_basis_size; j++ ){
            vector<vector< complex<double> > > v2 ;
            for(k = 0; k< 2 * nmodes[0] + 1 ; k++ ){
                vector<complex<double>> v3 ( 2 * nmodes[0] , 0 );
                v2.push_back(v3);
            }
            v1.push_back(v2);
        }
        Haar_state_overlap_time_dependent_basis_set .push_back(v1);
    }

    // Haar_state_with_ladder_operator_x_sparsify , y_sparsify , basis_set_sparsify
    Haar_state_with_ladder_operator_x_sparsify = new vector<vector<double>> * [N_Haar];
     Haar_state_with_ladder_operator_y_sparsify = new vector<vector<double>> * [N_Haar] ;
    Haar_state_with_ladder_operator_basis_set_sparsify =  new vector<vector<int>> * [N_Haar];
    for(i=0;i<N_Haar ; i++ ){
        Haar_state_with_ladder_operator_x_sparsify[i] = new vector<vector<double>> [2 * nmodes[0] + 1 ];
        Haar_state_with_ladder_operator_y_sparsify[i] = new vector<vector<double>> [2 * nmodes[0] + 1 ];
        Haar_state_with_ladder_operator_basis_set_sparsify[i] = new vector<vector<int>> [2 * nmodes[0] + 1 ];
    }

}

void detector::compute_Haar_random_state_with_ladder_operator (double sparsify_criteria ){
    // compute  a_{i} e^{-iHt} a_{j} e^{-\beta H/4} |\phi_{Haar}>   and   a_{i}e^{-iHt} e^{-\beta H/4} |\phi_{Haar}>
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
            Haar_state_with_ladder_operator_x_sparsify[i][j].clear();
            Haar_state_with_ladder_operator_y_sparsify[i][j].clear();
            Haar_state_with_ladder_operator_basis_set_sparsify[i][j].clear();
        }
    }


    for(i=0;i<N_Haar ; i++){
        for(j=0;j<2*nmodes[0]+1; j++ ){
            state_index = Haar_state_index_list[i] + j ;
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
                    if(magnitude > sparsify_criteria * normalization){
                        vx.push_back(xd_for_ladder_operation[l][k]);
                        vy.push_back(yd_for_ladder_operation[l][k]);
                        v_basis_set.push_back(k);
                    }
                }


                Haar_state_with_ladder_operator_x_sparsify[i][j].push_back(vx);
                Haar_state_with_ladder_operator_y_sparsify[i][j].push_back(vy);
                Haar_state_with_ladder_operator_basis_set_sparsify[i][j].push_back(v_basis_set);

            }

        }

    }

}

void detector:: compute_Haar_random_state_with_ladder_operator_overlap_with_time_dependent_basis_set(   ){
    // compute  <m(t) | a_{i} e^{-iHt} a_{j} e^{-\beta H/4} |Haar>
    // results in Haar_state_overlap_time_dependent_basis_set . size [N_Haar] [N_basis_set ] [2 * N_dof + 1 ] [ 2 * N_dof ]
    // Here dimension:  N_Haar is for |Haar> , N_basis_set is for |m(t)> , 2 * N_dof + 1 is for a_{j}  2 * N_dof is for a_{i}
    // computational cost: N_basis * N_Haar * (N_dof) ^2 * O(sparsified basis set )
    int i,j , k, l, m;
    complex<double> overlap;
    complex<double> overlap_sum;
    double overlap_real;
    double overlap_imag;
    double overlap_real_sum;
    double overlap_imag_sum;

    int basis_index;
    double real_part_value_Haar_state;
    double imag_part_value_Haar_state;

    double real_part_value_basis_state;
    double imag_part_value_basis_state;


    int sparsify_basis_list_size;
    // compute a_{i} e^{-iHt/\hbar} a_{j} e^{-\beta H/4} |\phi_{Haar}>
    compute_Haar_random_state_with_ladder_operator();

    int nearby_state_basis_size = nearby_state_index.size();
    for(i=0;i<N_Haar; i++ ){
        for(m=0;m<nearby_state_basis_size; m++ ){
            for(j=0;j<2 * nmodes[0]+ 1 ; j++  ){
                for(k=0; k<2 *nmodes[0]; k++ ){

                    overlap = 0;
                    overlap_sum = 0 ;
                    sparsify_basis_list_size = Haar_state_with_ladder_operator_basis_set_sparsify[i][j][k].size();
                    for(l=0;l<sparsify_basis_list_size;l++){
                        basis_index = Haar_state_with_ladder_operator_basis_set_sparsify[i][j][k][l];
                        real_part_value_Haar_state = Haar_state_with_ladder_operator_x_sparsify[i][j][k][l];
                        imag_part_value_Haar_state = Haar_state_with_ladder_operator_y_sparsify[i][j][k][l];

                        real_part_value_basis_state = xd[m][basis_index];
                        imag_part_value_basis_state = yd[m][basis_index];

                        overlap = overlap + complex<double> (real_part_value_basis_state, - imag_part_value_basis_state) *
                                complex<double> (real_part_value_Haar_state , imag_part_value_Haar_state) ;

                    }

                    // sum result in different process
                    overlap_real = real(overlap);
                    overlap_imag = imag(overlap);
                    MPI_Allreduce(&overlap_real, &overlap_real_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(&overlap_imag, &overlap_imag_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    overlap_sum = complex<double> (overlap_real_sum, overlap_imag_sum);

                    Haar_state_overlap_time_dependent_basis_set[i][m][j][k] = overlap_sum ;

                }
            }
        }
    }


}

void detector:: compute_regularized_thermal_OTOC_component(){
    // compute <m| e^{-\beta H / 4} [a_{i}(t) , a_{j} ] e^{-\beta H / 4 } |Haar> .
    // result store in regularized_thermal_OTOC_overlap_Haar_state_basis_set . dimension [N_Haar][N_basis_set ] [2 * N_dof] [2* N_dof]
    // computational cost for this part :  N_basis * N_Haar * N_dof^2 * O(sparsified_state_number)
    int i, j, k, l , m;

    int nearby_state_basis_size = nearby_state_index.size();

    int begin_index = total_dmat_size[0] / num_proc * my_id ;  // actually here total_dmat_size[0] == nearby_state_basis_size

    int Boltzmann_weighted_list_size;
    int local_basis_set_index;
    int basis_set_index;

    int mode_index;

    complex<double> C1;
    complex<double> C2;
    complex<double> C1_sum;
    complex<double> C2_sum;
    double C1_real;
    double C1_imag;
    double C2_real;
    double C2_imag;
    double C1_real_sum;
    double C1_imag_sum;
    double C2_real_sum;
    double C2_imag_sum;



    complex<double> Boltzmann_weighted_basis_basis_overlap;
    complex<double> overlap_with_time_dependent_basis_set;

    // compute <m~|e^{iHt} a_{i} e^{-iHt} a_{j} e^{-\beta H/4} |Haar>
    compute_Haar_random_state_with_ladder_operator_overlap_with_time_dependent_basis_set();

    for(i=0;i<N_Haar;i++){
        for(m=0;m< nearby_state_basis_size;m++){
            for(j=0;j<2*nmodes[0] ; j++) {
                for(k=0;k<2*nmodes[0] ; k++ ){
                    // compute <m| y [a_{j}(t) , a_{k}] y | Haar>

                    // compute <m| y a_{j}(t) a_{k} y | Haar> . result store in C1
                    C1 = 0;
                    Boltzmann_weighted_list_size = Boltzmann_weighted_basis_index_sparsify[m].size();
                    for(l=0;l<Boltzmann_weighted_list_size; l++ ){
                        local_basis_set_index = Boltzmann_weighted_basis_index_sparsify[m][l] ;
                        basis_set_index = local_basis_set_index + begin_index;
                        // <m|e^{-\beta H/4} | m~>:  index l stands for |m~> . index m stands for |m>
                        Boltzmann_weighted_basis_basis_overlap = complex<double> ( Boltzmann_factor_weighted_x_sparsify[m][l]  ,
                                                                                   - Boltzmann_factor_weighted_y_sparsify[m][l] );
                        // <m~| a_{j}(t) a_{k}|Haar>
                        overlap_with_time_dependent_basis_set =  Haar_state_overlap_time_dependent_basis_set[i][basis_set_index][k+1][j]  ;
                        C1 = C1 + Boltzmann_weighted_basis_basis_overlap * overlap_with_time_dependent_basis_set ;
                    }
                    C1_real = real(C1);
                    C1_imag = imag(C1);
                    MPI_Allreduce(&C1_real, &C1_real_sum , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(&C1_imag, &C1_imag_sum , 1,  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    C1_sum = complex<double> (C1_real_sum, C1_imag_sum);


                    // compute <m | y a_{k} a_{j}(t) y | Haar> , result store in C2
                    C2 = 0;
                    // This is for <m|y a_{k} | m~>  = conjugate (<m~| a_{k}^{+} y | m>)
                    if(k<nmodes[0]){
                        mode_index = k + nmodes[0];
                    }
                    else{
                        mode_index = k - nmodes[0];
                    }

                    Boltzmann_weighted_list_size = ladder_operator_Boltzmann_weighted_basis_index_sparsify[m][mode_index].size();
                    for(l=0;l<Boltzmann_weighted_list_size;l++){
                        // this basis set index is local index in process.
                        local_basis_set_index = ladder_operator_Boltzmann_weighted_basis_index_sparsify[m][mode_index][l];
                        basis_set_index = local_basis_set_index + begin_index;
                        // <m| e^{-\beta H/4} a_{k} | m~> : index l stands for |m~> , index m stands for |m>
                        Boltzmann_weighted_basis_basis_overlap = complex<double> ( ladder_operator_Boltzmann_weighted_x_sparsify[m][mode_index][l] ,
                                                                                   - ladder_operator_Boltzmann_weighted_y_sparsify[m][mode_index][l]  ) ;
                        // <m~|a_{j}(t) y | Haar>
                        overlap_with_time_dependent_basis_set = Haar_state_overlap_time_dependent_basis_set[i][basis_set_index][0][j];
                        C2 = C2 + Boltzmann_weighted_basis_basis_overlap * overlap_with_time_dependent_basis_set ;
                    }

                    C2_real = real(C2);
                    C2_imag = imag(C2);
                    MPI_Allreduce(&C2_real, &C2_real_sum , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    MPI_Allreduce(&C2_imag, &C2_imag_sum , 1,  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    C2_sum = complex<double> (C2_real_sum, C2_imag_sum);

                     // <m| y [a_{j}(t) , a_{k}] y | Haar>
                    regularized_thermal_OTOC_overlap_Haar_state_basis_set [i][m][j][k] = C1_sum - C2_sum ;

                }
            }
        }
    }



}

void detector:: compute_regularized_thermal_OTOC_Lyapunov_spectrum( ){
    // result sotre in regularized_thermal_Lyapunov_spectrum_each_Haar_state : size : [N_Haar][ 2* nmodes[0]] [2 * nmodes[0] ]
    // computational cost for this part : O(N_dof ^3 ) * N_basis_set * N_Haar
    int i, j , k, l, m;
    complex<double> Lyapunov_spectrum_component;
    int nearby_state_basis_size = nearby_state_index.size();
    compute_regularized_thermal_OTOC_component() ;

    int Haar_index;
    for(Haar_index = 0; Haar_index <N_Haar; Haar_index ++ ){
        for(i=0;i< 2 * nmodes[0]; i++){
            for(j=0;j<2 * nmodes[0] ; j++ ){
                // compute L_{ij} for Haar state |\phi>
                Lyapunov_spectrum_component = 0 ;
                for(m=0;m< nearby_state_basis_size ; m++ ){
                    for(k=0; k < 2*nmodes[0]; k++ ){
                        Lyapunov_spectrum_component = Lyapunov_spectrum_component +
                                std::conj(regularized_thermal_OTOC_overlap_Haar_state_basis_set[Haar_index][m][k][i]) *
                                regularized_thermal_OTOC_overlap_Haar_state_basis_set[Haar_index][m][k][j] ;
                    }
                }

                regularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j]  = Lyapunov_spectrum_component ;
                regularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j] =
                        regularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j] / Haar_state_normalization_list[Haar_index] ;

            }
        }
    }

    // average result
    for(i=0;i<2*nmodes[0]; i++ ){
        for(j=0;j<2*nmodes[0]; j++){
            regularized_thermal_Lyapunov_spectrum[i][j] = 0;
            for(Haar_index = 0; Haar_index < N_Haar; Haar_index ++ ){
                regularized_thermal_Lyapunov_spectrum[i][j] = regularized_thermal_Lyapunov_spectrum[i][j] +
                                                                regularized_thermal_Lyapunov_spectrum_each_Haar_state[Haar_index][i][j];
            }
            regularized_thermal_Lyapunov_spectrum[i][j] = regularized_thermal_Lyapunov_spectrum[i][j] / double(N_Haar);

        }
    }

}
