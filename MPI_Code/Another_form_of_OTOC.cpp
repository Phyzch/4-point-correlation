//
// Created by phyzch on 12/21/20.
//
#include "../system.h"
#include "../util.h"

void detector:: compute_n_off_diag_element_one_mode_quanta_below(int index_b, int index_a, complex<double> * n_off_diag_element){
    // compute <a_{i} -  | n_{i}(t) | b_{i} - > : here | b_{i} - > is the state one mode quanta below coordinate i from state b. same for |a_{i}->
    int i,j,k;
    int a_one_mode_quanta_below;
    int b_one_mode_quanta_below;

    complex<double> C_kb; //   C_{b}^{k}(t) = <k|b(t)>
    complex<double> C_ka_conjugate; //  (C^{a}_{k}(t))^{*}
    double n_k_i;

    for(i=0;i<nmodes[0];i++){
        n_off_diag_element[i] = 0;
    }

    for(i=0;i<nmodes[0];i++){
        a_one_mode_quanta_below = neighbor_state_in_nearby_state_index_list[index_a][i+nmodes[0]]; // one mode up for first dof number element, then one mode below
        b_one_mode_quanta_below = neighbor_state_in_nearby_state_index_list[index_b][i+nmodes[0]];

        for(k=0;k<dmatsize[0];k++){
            C_ka_conjugate = complex<double> (xd[a_one_mode_quanta_below][k],
                                              -yd[a_one_mode_quanta_below][k]);
            C_kb = complex<double> (xd[b_one_mode_quanta_below][k],yd[b_one_mode_quanta_below][k]);
            n_k_i = dv[0][ k ][i];

            n_off_diag_element[i] = n_off_diag_element[i] + C_ka_conjugate * C_kb * n_k_i;
        }

    }

}

void full_system :: compute_another_form_of_OTOC(int nearby_state_index_size, complex<double> * n_offdiag_element,
                                                 complex<double> ** n_offdiag,double ** n_offdiag_real, double ** n_offdiag_imag,
                                                 complex<double> **n_offdiag_total, double ** n_offdiag_total_real, double ** n_offdiag_total_imag,
                                                 complex<double> ** n_offdiag_one_mode_quanta_below, double ** n_offdiag_one_mode_quanta_below_real, double ** n_offdiag_one_mode_quanta_below_imag,
                                                 complex<double> ** n_offdiag_total_one_mode_quanta_below, double ** n_offdiag_total_one_mode_quanta_below_real, double ** n_offdiag_total_one_mode_quanta_below_imag,
                                                 int initial_state_index_in_total_dmatrix ,
                                                 double * another_OTOC){
    int i,j;
    int b;
    vector<int> nearby_state_1;
    vector<int> nearby_state_2;

    complex<double> variable1;


    for(b=0;b<nearby_state_index_size;b++) {
        if (d.neighbor_state_all_in_nearby_state_index_list[b]) {
            d.compute_n_off_diag_element(b, d.initial_state_index_in_nearby_state_index_list, n_offdiag_element);
            for (i = 0; i < d.nmodes[0]; i++) {
                n_offdiag[i][b] = n_offdiag_element[i];
                n_offdiag_real[i][b] = real(n_offdiag[i][b]);
                n_offdiag_imag[i][b] = imag(n_offdiag[i][b]);
            }

            d.compute_n_off_diag_element_one_mode_quanta_below(b, d.initial_state_index_in_nearby_state_index_list,
                                                               n_offdiag_element);
            for (i = 0; i < d.nmodes[0]; i++) {
                n_offdiag_one_mode_quanta_below[i][b] = n_offdiag_element[i];
                n_offdiag_one_mode_quanta_below_real[i][b] = real(n_offdiag_element[i]);
                n_offdiag_one_mode_quanta_below_imag[i][b] = imag(n_offdiag_element[i]);
            }

        }

        else {
            for (i = 0; i < d.nmodes[0]; i++) {
                n_offdiag[i][b] = 0;
                n_offdiag_real[i][b] = 0;
                n_offdiag_imag[i][b] = 0;
                n_offdiag_one_mode_quanta_below[i][b] = 0;
                n_offdiag_one_mode_quanta_below_real[i][b] = 0;
                n_offdiag_one_mode_quanta_below_imag[i][b] = 0;
            }
        }

    }
    for(i=0;i<d.nmodes[0];i++){
        MPI_Allreduce(&n_offdiag_real[i][0], & n_offdiag_total_real[i][0],nearby_state_index_size, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&n_offdiag_imag[i][0], & n_offdiag_total_imag[i][0], nearby_state_index_size, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        MPI_Allreduce(&n_offdiag_one_mode_quanta_below_real[i][0], & n_offdiag_total_one_mode_quanta_below_real[i][0] , nearby_state_index_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(& n_offdiag_one_mode_quanta_below_imag[i][0], & n_offdiag_total_one_mode_quanta_below_imag[i][0], nearby_state_index_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }

    for(b=0;b<nearby_state_index_size;b++){
        for(i=0;i<d.nmodes[0];i++){
            n_offdiag_total[i][b] = complex<double> (n_offdiag_total_real[i][b], n_offdiag_total_imag[i][b]);
            n_offdiag_total_one_mode_quanta_below [i][b] = complex<double> (n_offdiag_total_one_mode_quanta_below_real[i][b],
                                                                            n_offdiag_total_one_mode_quanta_below_imag[i][b]);
        }
    }

    // Now compute | \sqrt{n_{i}^{l} } <m| n_{i}(t) |l> - \sqrt{n_{i}^{m}} <m^{-} | n_{i}(t) | l-> |^{2}
    for(i=0;i<d.nmodes[0];i++){
        another_OTOC[i] = 0;
        for(b=0;b<nearby_state_index_size;b++){
            variable1 = sqrt(d.dv_all[0][initial_state_index_in_total_dmatrix][i]) * n_offdiag_total[i][b] -
                        sqrt(d.dv_all[0][d.nearby_state_index[b]][i]) * n_offdiag_total_one_mode_quanta_below[i][b];
            another_OTOC[i] = another_OTOC[i] + d.dv_all[0][initial_state_index_in_total_dmatrix][i] *
                                                std::norm(variable1);
        }

    }

}
