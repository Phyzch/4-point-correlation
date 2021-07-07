//
// Created by phyzch on 7/5/21.
//

// Use Chebyshev polynomial to expand e^{-\beta H}.   exp(Rx) [x in [-1 ,1] ] \approx  \sum_{k=0}^{n} a_{k}(R) T_{k}(x)
// a_{k}(R) = C_{k} I_{k}(R)  .  C_{k} = 1 if k=-=0 , C_{k} = 2 otherwise.   I_{k}(R) is modified Bessel function of first kind.
// When k>=1.5 R.  I_{k}(R) will be very small and decrease exponentially.

#include "../util.h"
#include "../system.h"
int N_Harr = 5 ; // Haar random state number.
double Boltzmann_beta; // 1/T.

void detector::prepare_compute_Boltzmann_factor_use_Chebyshev_polynomial(double one_fourth_beta , ofstream & log ){
    // except to call this function, we also have to call prepare_evolution function to prepare parallel computation : remote_Vec_count etc.
    // this function should be called after prepare_evolution() function
    // beta == 1/T
    int i, j;
    double max_Hamiltonian, min_Hamiltonian;
    max_Hamiltonian = total_dmat[0][0];
    min_Hamiltonian = total_dmat[0][0];
    for(i=0;i<total_dmat_size[0];i++){
        if (max_Hamiltonian < total_dmat[0][i]){
            max_Hamiltonian = total_dmat[0][i];
        }
        if (min_Hamiltonian > total_dmat[0][i]){
            min_Hamiltonian = total_dmat[0][i];
        }
    }

    Chebyshev_e0 = ( max_Hamiltonian + min_Hamiltonian )/ 2;
    Chebyshev_R  = (max_Hamiltonian - min_Hamiltonian ) * 0.55 ;

    shifted_total_dmat = new double [total_dmat_num[0]];
    shifted_dmat = new double [dmatnum[0]];

    for(i=0;i<total_dmat_num[0];i++){
        if(total_dirow[0][i] == total_dicol[0][i]){
            shifted_total_dmat[i] = total_dmat[0][i] - Chebyshev_e0;
        }
        else{
            shifted_total_dmat[i] = total_dmat[0][i];
        }
        shifted_total_dmat[i] =  shifted_total_dmat[i] / Chebyshev_R ;
    }

    for(i=0;i<dmatnum[0];i++){
        if(dirow[0][i] == dicol[0][i]){
            shifted_dmat[i] = dmat[0][i] - Chebyshev_e0;
        }
        else{
            shifted_dmat[i] = dmat[0][i];
        }
        shifted_dmat[i] = shifted_dmat[i] / Chebyshev_R;
    }

    Chebyshev_prefactor = std::exp( - one_fourth_beta * Chebyshev_e0);
    Chebyshev_R_beta = Chebyshev_R * one_fourth_beta;

    N_chebyshev = ceil(1.5 * Chebyshev_R);
    if(my_id == 0){
        cout <<"Using Chebyshev method to compute Boltzmann factor. beta * R / 4 = " << Chebyshev_R_beta << endl;
        log << "Using Chebyshev method to compute Boltzmann factor. beta * R / 4 = " << Chebyshev_R_beta << endl;
        cout << "Using Chebychev method.  order of Chebychev polynomial cutoff = " << N_chebyshev << endl;
        log << "Using Chebychev method.  order of Chebychev polynomial cutoff = " << N_chebyshev << endl;
    }

    Bessel_function_array = new double [N_chebyshev + 1] ;
    for(i=0; i<= N_chebyshev; i++){
        Bessel_function_array[i] = std::cyl_bessel_i(i,Chebyshev_R_beta) ;  // modified Bessel function of first kind I_{i}(x)
    }

    // Used for recursive solving T_{n}(x) : T_{n+1}(x) = 2* x * T_{n}(x) - T_{n-1}(x)
    Chebyshev_polyn = new vector<double> [6];
    for(i=0;i<6;i++){
        vector <double> v (dmatsize[0] + to_recv_buffer_len[0] , 0 );
        Chebyshev_polyn[i] = v;
    }

    send_polyn = new double * [6];
    recv_polyn = new double * [6];
    for(i=0; i<6 ;i++ ){
        send_polyn[i] = new double [to_send_buffer_len[0]];
        recv_polyn[i] = new double [to_recv_buffer_len[0]];
    }

}

void detector::update_polyn23() {
    // update Chebychev_polyn[2] and Chebychev_polyn[3]
    int i;
    int vsize;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;
    for(i=0;i<to_send_buffer_len[0];i++){
        send_polyn[2][i] = Chebyshev_polyn[2][tosendVecIndex[0][i] - my_id * vsize];
        send_polyn[3][i] = Chebyshev_polyn[3][tosendVecIndex[0][i] - my_id * vsize];
    }
    MPI_Alltoallv(&send_polyn[2][0],tosendVecCount[0],tosendVecPtr[0],MPI_DOUBLE,
                  &recv_polyn[2][0],remoteVecCount[0],remoteVecPtr[0],MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Alltoallv(&send_polyn[3][0],tosendVecCount[0],tosendVecPtr[0],MPI_DOUBLE,
                  &recv_polyn[3][0],remoteVecCount[0],remoteVecPtr[0],MPI_DOUBLE,MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len[0];i++){
        Chebyshev_polyn[2][ i + dmatsize[0]]  = recv_polyn[2][i];
        Chebyshev_polyn[3][ i + dmatsize[0]]  = recv_polyn[3][i];
    }
}


void detector:: Chebyshev_method_Boltzmann_factor(const  vector<double> & wave_func_x ,const vector<double> & wave_func_y,
                                                  vector<double> & Boltzmann_factor_weighted_wave_func_x, vector<double> & Boltzmann_factor_weighted_wave_func_y ){
    // compute Bolzmann weighted wave function .
    // input :: one_fourth_beta,  wave_func_x,  wave_func_y
    // output :: Boltzmann_factor_weighted_wave_func_x ,  Boltzmann_factor_weighted_wave_func_y
    // caution: before call this function, make sure update wave function component from other process.
    int i, j, k , m;
    int irow , icol;
    double bess;
    double prefactor;

    Boltzmann_factor_weighted_wave_func_x = wave_func_x;
    Boltzmann_factor_weighted_wave_func_y = wave_func_y;


    vector<double> & creal =   Boltzmann_factor_weighted_wave_func_x;
    vector<double> & cimag = Boltzmann_factor_weighted_wave_func_y  ;

    // Chebyshev_polyn is Chebyshev polynomial of Hamiltonian act upon wave function
    // 0,2,4 for real part. 1,3,5 for imag part
    for(i=0;i<dmatsize[0] + to_recv_buffer_len[0]; i++) {
        Chebyshev_polyn[0][i] = creal[i];  // now it's just xd
        Chebyshev_polyn[1][i] = cimag[i];  // now it's just yd
        Chebyshev_polyn[2][i] = 0;
        Chebyshev_polyn[3][i] = 0;
        Chebyshev_polyn[4][i] = 0;
        Chebyshev_polyn[5][i] = 0;
    }
    // zeroth order
    bess = Bessel_function_array[0];
    prefactor = Chebyshev_prefactor * bess;  // Chebyshev prefactor : e^{-\beta E_{0}} . bess = I_{0}(R)
    for(i=0;i<dmatsize[0];i++){
        creal[i] = prefactor * Chebyshev_polyn[0][i];
        cimag[i] = prefactor * Chebyshev_polyn[1][i];
    }

    // first order
    bess = Bessel_function_array[1];
    prefactor = 2 * bess * Chebyshev_prefactor;
    for(i=0;i<dmatnum[0];i++){  // \omega =  - shifted_dmat
        irow = local_dirow[0][i];
        icol = local_dicol[0][i];
        Chebyshev_polyn[2][irow ] = Chebyshev_polyn[2][irow]  + (-shifted_dmat[i])  * Chebyshev_polyn[0][icol];
        Chebyshev_polyn[3][irow ] = Chebyshev_polyn[3][irow] +  (-shifted_dmat[i]) * Chebyshev_polyn[1][icol];
    }

    // used for communication between different process
    update_polyn23();
    for(i=0;i<dmatsize[0];i++){
        creal[i] = creal[i] + prefactor * Chebyshev_polyn[2][i];
        cimag[i] = cimag[i] + prefactor * Chebyshev_polyn[3][i];
    }

    // Remaining terms :
    for(k=2;k<=N_chebyshev;k++){
        bess = Bessel_function_array[k];
        prefactor = Chebyshev_prefactor * bess * 2 ;

        // Use Chebyshev polynomial relationship  T_{k+2}(x) = 2 * x *  T_{k+1}(x) - T_{k}(x)
        for(i=0;i<dmatnum[0]; i++ ){
            irow = local_dirow[0][i];
            icol = local_dicol[0][i];
            Chebyshev_polyn[4][irow] = Chebyshev_polyn[4][irow] + 2 * (-shifted_dmat[i]) * Chebyshev_polyn[2][icol];
            Chebyshev_polyn[5][irow] = Chebyshev_polyn[5][irow] + 2 * (-shifted_dmat[i]) * Chebyshev_polyn[3][icol];
        }
        for(i=0;i<dmatsize[0];i++){
            Chebyshev_polyn[4][i] = Chebyshev_polyn[4][i] + Chebyshev_polyn[0][i];
            Chebyshev_polyn[5][i] = Chebyshev_polyn[5][i] + Chebyshev_polyn[1][i];
        }

        // update dx , dy
        for(i=0;i<dmatsize[0];i++){
            creal[i] = creal[i]  + prefactor * Chebyshev_polyn[4][i];
            cimag[i] = cimag[i] + prefactor * Chebyshev_polyn[5][i];
        }

        for(i=0;i<dmatsize[0];i++){
            Chebyshev_polyn[0][i] = Chebyshev_polyn[2][i];
            Chebyshev_polyn[1][i] = Chebyshev_polyn[3][i];
            Chebyshev_polyn[2][i] = Chebyshev_polyn[4][i];
            Chebyshev_polyn[3][i] = Chebyshev_polyn[5][i];
            Chebyshev_polyn[4][i] = 0;
            Chebyshev_polyn[5][i] = 0;
        }
        update_polyn23();

    }



}

void detector::Boltzmann_factor_decorated_basis_set_and_with_ladder_operator(){
    int i, j, k, m;
    int nearby_state_basis_size = nearby_state_index.size();
    update_dx(nearby_state_basis_size);
    update_dy(nearby_state_basis_size);

    vector<vector<double>> xd_for_ladder_operator;
    vector<vector<double>> yd_for_ladder_operator;

    double sparsify_criteria =  pow(10,-2);
    double normalization;
    double normalization_tot;
    double magnitude;


    for(m=0; m<nearby_state_basis_size; m++){
        vector<double> Boltzmann_weighted_x;
        vector<double> Boltzmann_weighted_y;

        vector<double> x_sparsify;
        vector<double> y_sparsify;
        vector<int> basis_set_index_sparsify;

        Chebyshev_method_Boltzmann_factor(xd[m] , yd[m], Boltzmann_weighted_x, Boltzmann_weighted_y);

        // normalization
        normalization = 0;
        for(k=0;k<dmatsize[0];k++){
            normalization = normalization + std::norm(Boltzmann_weighted_x[k]) + std::norm( Boltzmann_weighted_y[k] );
        }
        MPI_Allreduce(&normalization , &normalization_tot , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        normalization = normalization_tot;
        normalization = sqrt(normalization);

        //sparsify
        for(k=0;k<dmatsize[0];k++){
            magnitude =  sqrt( norm(Boltzmann_weighted_x[k]) + norm(Boltzmann_weighted_y[k]) );
            if( magnitude > sparsify_criteria * normalization ){
                x_sparsify.push_back(Boltzmann_weighted_x[k]);
                y_sparsify.push_back(Boltzmann_weighted_y[k]);
                basis_set_index_sparsify.push_back(k);
            }
        }
        Boltzmann_factor_weighted_x_sparsify.push_back(x_sparsify);
        Boltzmann_factor_weighted_y_sparsify.push_back(y_sparsify);
        Boltzmann_weighted_basis_index_sparsify.push_back(basis_set_index_sparsify);

        // compute a_{j} e^{-\beta H} | \phi>  ladder operator Boltzmann weighted basis set
        vector<vector<double>> wave_func_x_sparsify_ladder;
        vector<vector<double>> wave_func_y_sparsify_ladder;
        vector<vector<int>> basis_set_index_sparsify_ladder;

        xd_for_ladder_operator.clear();
        yd_for_ladder_operator.clear();
        ladder_operator_operation(Boltzmann_weighted_x, Boltzmann_weighted_y ,
                                  xd_for_ladder_operator, yd_for_ladder_operator);

        for(j=0;j<2*nmodes[0];j++){

            vector<double> x_sparsify_ladder;
            vector<double> y_sparsify_ladder;
            vector<int> index_sparsify_ladder;

            // compute magnitude of a_{j} e^{-\beta H} | {n} >
            normalization = 0;
            for(k=0;k<dmatsize[0];k++){
                normalization = normalization + std::norm(xd_for_ladder_operator[j][k]) + std::norm(yd_for_ladder_operator[j][k]);
            }
            MPI_Allreduce(&normalization , &normalization_tot , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            normalization = normalization_tot;
            normalization = sqrt(normalization);

            for(k=0;k<dmatsize[0];k++){
                magnitude = sqrt( norm(xd_for_ladder_operator[j][k]) + norm(yd_for_ladder_operator[j][k]) );
                if(   magnitude > sparsify_criteria * normalization  ){
                    x_sparsify_ladder.push_back(xd_for_ladder_operator[j][k]);
                    y_sparsify_ladder.push_back(yd_for_ladder_operator[j][k]);
                    index_sparsify_ladder.push_back(k);
                }
            }

            wave_func_x_sparsify_ladder.push_back(x_sparsify_ladder);
            wave_func_y_sparsify_ladder.push_back(y_sparsify_ladder);
            basis_set_index_sparsify_ladder.push_back(index_sparsify_ladder);

        }

        ladder_operator_Boltzmann_weighted_x_sparsify.push_back(wave_func_x_sparsify_ladder);
        lader_operator_Boltzmann_weighted_y_sparsify.push_back(wave_func_y_sparsify_ladder);
        ladder_operator_Boltzmann_weighted_basis_index_sparsify.push_back(basis_set_index_sparsify_ladder);

    }

}