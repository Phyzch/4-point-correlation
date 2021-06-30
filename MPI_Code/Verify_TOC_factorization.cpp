//
// Created by phyzch on 6/27/21.
//
#include "../util.h"
#include "../system.h"

void detector:: compute_Time_ordered_correlation_func_per_state(int state_m, int state_l, double ** M_matrix,
                                 vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                 vector<vector<int>> * index_for_xp_sparsify){
    // M_{ij}(t) = <m| z_{i}(t) z_{j} |l> . Here z is x if ( i < N) and z is p if ( i > N )
    // x == a + a^{+} (we ignore square root 2 here)   p == a - a^{+} we ignore square root 2 here.
    int i,j,k;
    complex<double> C1;  // C1 == a_{i}(t) a_{j}
    complex<double> C2;  // C2 == a_{i}^{+}(t) a_{j}
    complex<double> C3;  // C3 == a_{i}(t) a_{j}^{+}
    complex<double> C4;  // C4 == a_{i}^{+}(t) a_{j}^{+}

    for(i=0; i<nmodes[0];i++ ){
        for(j=0;j<nmodes[0];j++){
            if(neighbor_state_in_nearby_state_index_list[initial_state_index_in_nearby_state_index_list][j] != -1){
                C1 = compute_c_overlap(state_m, j+1, i, xd_for_xp_sparsify, yd_for_xp_sparsify, index_for_xp_sparsify);
                C2 = compute_c_overlap(state_m,j+1, i+nmodes[0], xd_for_xp_sparsify, yd_for_xp_sparsify, index_for_xp_sparsify);
            }
            else{
                C1 = 0;
                C2 = 0;
            }
            if(neighbor_state_in_nearby_state_index_list[initial_state_index_in_nearby_state_index_list][j + nmodes[0]] != -1){
                C3 = compute_c_overlap(state_m, j+1+nmodes[0],i,xd_for_xp_sparsify,yd_for_xp_sparsify,index_for_xp_sparsify);
                C4 = compute_c_overlap(state_m,j+1+nmodes[0],i+nmodes[0], xd_for_xp_sparsify, yd_for_xp_sparsify, index_for_xp_sparsify);
            }
            else{
                C3 = 0;
                C4 = 0;
            }


            // ladder operator \sqrt{n} or \sqrt{n+1}
            C1 = C1 *  sqrt( dv_all[0][ nearby_state_index[state_l] ][j]) ;
            C2 = C2 * sqrt(dv_all[0][ nearby_state_index[state_l] ][j] );
            C3 = C3 * sqrt(dv_all[0][ nearby_state_index[state_l] ][j] + 1 );
            C4 = C4 * sqrt(dv_all[0][ nearby_state_index[state_l] ][j] + 1 );

            // [i][j] == |<m|x_{i}(t) x_{j}|l>|^2 , [i+nmodes[0]][j] == |<m|p_{i}(t) x_{j}|l>|^2 , [i][j+nmodes[0]] == |<m|x_{i}(t) p_{j}|l>|^2
            M_matrix [i] [j] = norm(C1 + C2 + C3 + C4) ;
            M_matrix [i+ nmodes[0] ] [j] = norm( C1 - C2 + C3 - C4 ) ;
            M_matrix [i] [j+ nmodes[0] ] = norm( C1 + C2 - C3 - C4 ) ;
            M_matrix [i+ nmodes[0] ] [j+ nmodes[0] ] = norm( C1 - C2 - C3 + C4 ) ;
        }
    }

}

double detector:: compute_two_point_func_zt_zt ( int mode_k, vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                            vector<vector<int>> * index_for_xp_sparsify , vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp ){
    // compute <m| z_{k}(t) z_{k}(t) |m> . If k>=dof, then z_{k} = p_{k}. Else z_{k} == x_{k}
    int i,j;
    int xd_for_xp_sparsify_list_len ;
    int basis_set_index;
    double real_part_value;
    double imag_part_value;

    double real_part_value2;
    double imag_part_value2;

    complex<double> C1;
    complex<double> C2;
    double C;

    int mode_k_index;
    if(mode_k < nmodes[0]){
        mode_k_index = mode_k;
    }
    else{
        mode_k_index = mode_k - nmodes[0];
    }

    double sign;
    if(mode_k<nmodes[0]){
        sign = 1;
    }
    else{
        sign = -1 ;
    }

    double norm_sum = 0;
    double norm_sum_all_proc = 0 ;

    // z_{k} = x_{k} == (a_{k} + a_{k}^{+})
    xd_for_xp_sparsify_list_len = xd_for_xp_sparsify[0][mode_k_index].size();  // for <n|a_{k}^{+}|phi(t)> = \sqrt{n_{k}} <n_{k}^{-}|\phi(t)>
    for(j=0;j<xd_for_xp_sparsify_list_len;j++){

        // \sqrt{n_{k}} <n_{k}^{-}|\phi(t)>
        basis_set_index = index_for_xp_sparsify[0][mode_k_index][j];
        real_part_value = xd_for_xp_sparsify[0][mode_k_index][j];
        imag_part_value = yd_for_xp_sparsify[0][mode_k_index][j];
        C1 = sign * complex<double> (real_part_value, imag_part_value) * sqrt(dv[0][basis_set_index][mode_k_index]);

        //  \sqrt{n_{k} + 1 } <n_{k}^{+}|\phi(t)>
        real_part_value2 = xd_for_xp[0][mode_k_index + nmodes[0] ][basis_set_index];
        imag_part_value2 = yd_for_xp[0][mode_k_index + nmodes[0] ][basis_set_index];
        C2 = complex<double> (real_part_value2 , imag_part_value2) * sqrt( dv[0][basis_set_index][mode_k_index] + 1 );

        C = norm(C1 + C2);
        norm_sum = norm_sum + C;
    }

    MPI_Allreduce(&norm_sum, &norm_sum_all_proc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return norm_sum_all_proc;
}

void detector::output_TOC_factorization(double * Two_point_func1 , double * Two_point_func2, double ** TOC_per_state,
                                        double ** TOC,
                                        vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp,
                                        vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                        vector<vector<int>> * index_for_xp_sparsify,
                                        double t,  ofstream & TOC_output ){
    int i,j,k,m;
    int nearby_state_index_size = nearby_state_index.size();
    double cutoff_criteria = pow(10,-2);

    for(i=0; i<2 *nmodes[0]; i++ ){
        for( j=0 ; j< 2* nmodes[0]; j++ ){
            TOC[i][j] = 0;
        }
    }

    update_xd_yd_for_xp(xd_for_xp,yd_for_xp,xd_for_xp_sparsify,yd_for_xp_sparsify,index_for_xp_sparsify, cutoff_criteria);

    // compute Time ordered correlation function (TOC) :
    for(m=0;m<nearby_state_index_size;m++){

        // go through state m to compute TOC per state
        compute_Time_ordered_correlation_func_per_state(m,initial_state_index_in_nearby_state_index_list, TOC_per_state,
                                                        xd_for_xp_sparsify, yd_for_xp_sparsify,index_for_xp_sparsify);

        for(i=0; i<2 *nmodes[0]; i++ ){
            for( j=0 ; j< 2* nmodes[0]; j++ ){
                // x_{i} = 1/\sqrt{2} (a_{i} + a_{i}^{+}). Thus for 4 point func we have to divide by 4
                TOC[i][j] = TOC[i][j] + TOC_per_state[i][j] / 4 ;
            }
        }

    }

    // Two point function1 is <\phi| z_{i}(t) z_{i}(t) |\phi>
    for(i=0; i<2 *nmodes[0]; i++){
        Two_point_func1[i] = compute_two_point_func_zt_zt(i, xd_for_xp_sparsify, yd_for_xp_sparsify, index_for_xp_sparsify, xd_for_xp, yd_for_xp);
        // x_{i} = 1/\sqrt{2} (a_{i} + a_{i}^{+}). Thus for 2 point func we have to divide by 2
        Two_point_func1[i] =  Two_point_func1[i] / 2;
    }

    // Two point function2 is <\phi|z_{i} z_{i}|\phi> . z_{i} = x_{i} or p_{i}
    for(i=0;i<nmodes[0]; i++ ){
        Two_point_func2[i] = 2 * dv_all[0][ nearby_state_index[initial_state_index_in_nearby_state_index_list] ][i] + 1 ;
        // x_{i} = 1/\sqrt{2} (a_{i} + a_{i}^{+}). Thus for 2 point func we have to divide by 2
        Two_point_func2[i] = Two_point_func2[i] / 2 ;
        Two_point_func2[i + nmodes[0]] = Two_point_func2[i] ;
    }

    TOC_output << t << endl;
    // First output <\phi|z_{i}(t) z_{i}(t) |\phi>
    for(i=0;i<2*nmodes[0];i++){
        TOC_output << Two_point_func1[i] << " ";
    }
    TOC_output << endl;
    TOC_output << endl;

    // Second: output <\phi|z_{i} z_{i} | \phi>
    for(i=0;i<2*nmodes[0];i++){
        TOC_output << Two_point_func2[i] <<" ";
    }
    TOC_output << endl;
    TOC_output << endl;

    // Third : output Time ordered correlation function : <\phi| z_{i}(t) z_{j} z_{j} z_{i}(t) |\phi>
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0];j++){
            TOC_output << TOC[i][j] <<" ";
        }
    }
    TOC_output << endl;
    TOC_output << endl;

    // Fourth: output epsilon == TOC - Two_point_func1 * Two_point_func2 . This is error for factorization.
    for(i=0;i<2*nmodes[0];i++){
        for(j=0;j<2*nmodes[0]; j++){
            TOC_output << TOC[i][j] - Two_point_func1[i] * Two_point_func2[j] << " ";
        }
    }
    TOC_output << endl;
    TOC_output << endl;

}
