//
// Created by phyzch on 10/10/21.
//
#include "../util.h"
# include"../system.h"

void detector::read_rotation_parameter(){
    int i, j;
    int vib_index1, vib_index2;
    int rotation_axis_index;
    double Coriolis_strength;
    string Coriolis_coefficient_folder = "/home/phyzch/CLionProjects/4_point_correlation_calculation/";
    string Coriolis_coefficient_file_name =  Coriolis_coefficient_folder +  "Coriolis_coeff.txt" ;
    ifstream coriolis_input;
    if(my_id == 0){
        coriolis_input.open(Coriolis_coefficient_file_name);
        for(i=0;i<3;i++){
            coriolis_input >> rotational_constant[i] ;
        }
        coriolis_input >> Coriolis_coupling_term_num ;

        Coriolis_coupling_rot_axis = new int [Coriolis_coupling_term_num];
        Coriolis_coupling_vib_index = new int * [Coriolis_coupling_term_num];
        for(i=0 ; i < Coriolis_coupling_term_num ; i++ ){
            Coriolis_coupling_vib_index[i] = new int [2];
        }
        Coriolis_coupling_strength = new double [Coriolis_coupling_term_num];

        for(i=0;i<Coriolis_coupling_term_num; i++){
            coriolis_input >> rotation_axis_index;
            Coriolis_coupling_rot_axis[i] = rotation_axis_index;

            coriolis_input >>vib_index1 >> vib_index2;
            vib_index1 = vib_index1 - 1 ;  // have to -1 to make index start from 0.
            vib_index2 = vib_index2 - 1;
            Coriolis_coupling_vib_index[i][0] = vib_index1;
            Coriolis_coupling_vib_index[i][1] = vib_index2;

            coriolis_input >> Coriolis_strength;
            Coriolis_coupling_strength[i] = Coriolis_strength;

        }

        coriolis_input.close();
    }

    MPI_Bcast(&rotational_constant[0] , 3, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    MPI_Bcast(&Coriolis_coupling_term_num , 1, MPI_INT, 0 , MPI_COMM_WORLD);

    if(my_id != 0){
        Coriolis_coupling_rot_axis = new int [Coriolis_coupling_term_num];
        Coriolis_coupling_vib_index = new int * [Coriolis_coupling_term_num];
        for(i=0 ; i < Coriolis_coupling_term_num ; i++ ){
            Coriolis_coupling_vib_index[i] = new int [2];
        }
        Coriolis_coupling_strength = new double [Coriolis_coupling_term_num];
    }

    MPI_Bcast(&Coriolis_coupling_rot_axis[0] , Coriolis_coupling_term_num , MPI_INT, 0 , MPI_COMM_WORLD);
    for ( i = 0 ; i < Coriolis_coupling_term_num ; i++ ){
        MPI_Bcast(&Coriolis_coupling_vib_index[i][0] , 2, MPI_INT, 0 , MPI_COMM_WORLD);
    }
    MPI_Bcast(&Coriolis_coupling_strength[0] , Coriolis_coupling_term_num , MPI_DOUBLE, 0, MPI_COMM_WORLD);



}

complex<double>  detector:: compute_rotational_offdiag_part_MPI(   vector<vector <int>> * vmode_ptr , int i , int j , int m ){
    // i is state index, j is state index.
    // m = 0 or 1. detector index
    int k;
    int begin_index;
    int ntot;
    complex<double> value;
    double lij;
    bool exist;
    int position;
    int rot_term_index;
    int coriolis_vib_term1;
    int coriolis_vib_term2;
    int rot_term_index2;
    int coriolis_vib_term3;
    int coriolis_vib_term4;

    int i_1, i_2, i_3, i_4;
    double coeff_1 , coeff_2, coeff_3, coeff_4;
    vector<int> quanta_move = {-1 , 1 };

    vector<int> quantum_num_new;

    complex<double> rot_sign_1;
    complex<double> rot_sign_2;
    int v_l , v_k;
    int rotation_quantum_m;

    int higher_rot_quanta;
    int rot_diff;


    if(i == j ) return 0 ;

    // check vib term diff for vibrational term
    value = 0;
    // difference in rotational quantum number M.
    rot_diff = abs ((*vmode_ptr)[i][nmodes[m]] - (*vmode_ptr)[j][nmodes[m]]);

    ntot = 0;
    for(k=0;k<nmodes[m];k++){
        deln[k]= abs( (*vmode_ptr)[i][k] - (*vmode_ptr)[j][k] ); // same as deln[k] = abs(dv[m][i][k] - dv[m][j][k]);
        ntot = ntot + deln[k] ;
    }

    if(ntot == 0){
        if( (*vmode_ptr)[i][nmodes[m]] == (*vmode_ptr)[j][nmodes[m]] ){
            cout  << "Error! The off-diagonal element must differ in q.n." << endl;
            MPI_Abort(MPI_COMM_WORLD,-8);
        }
        // J_x ^ 2 , J_y^2  coupling term.
        if(  rot_diff == 2  ){
            higher_rot_quanta = max((*vmode_ptr)[i][nmodes[m]] , (*vmode_ptr)[j][nmodes[m]]);
            value = value + (rotational_constant[0] + rotational_constant[1]) * 1/4 * pow (
                    (angular_momentum_J * (angular_momentum_J + 1 ) - higher_rot_quanta * (higher_rot_quanta - 1) )  * ( angular_momentum_J * (angular_momentum_J + 1 )  -  (higher_rot_quanta-2) * (higher_rot_quanta - 1) )
                    , 1/2) ;
        }

    }
    else{

        // Coriolis coupling
        if(ntot == 2 ){
            for(rot_term_index = 0; rot_term_index < Coriolis_coupling_term_num; rot_term_index ++ ){
                coriolis_vib_term1 = Coriolis_coupling_vib_index[rot_term_index][0];
                coriolis_vib_term2 = Coriolis_coupling_vib_index[rot_term_index][1];
                // check if this is coupling term.
                if(deln[coriolis_vib_term1] == 1 and deln[coriolis_vib_term2] == 1){
                    if( (*vmode_ptr)[i][coriolis_vib_term2] - (*vmode_ptr)[j][coriolis_vib_term2] == 1 ){
                        rot_sign_1 = complex<double>(0,1); // i
                    }
                    else{
                        rot_sign_1 = complex<double>(0,-1); // -i
                    }
                    v_k = max( (*vmode_ptr)[i][coriolis_vib_term1] , (*vmode_ptr)[j][coriolis_vib_term1] );
                    v_l = max( (*vmode_ptr)[i][coriolis_vib_term2] , (*vmode_ptr)[j][coriolis_vib_term2] );
                    if( Coriolis_coupling_rot_axis[rot_term_index] == 0 and rot_diff == 1){
                        // Jx
                        rot_sign_2 = 1;
                        rotation_quantum_m = min( (*vmode_ptr)[i][nmodes[m]] ,  (*vmode_ptr)[j][nmodes[m]] ) ;
                        value = value + (-2 * rotational_constant[0]) * Coriolis_coupling_strength[rot_term_index] * sqrt( double(v_k * v_l ) ) / 2 * rot_sign_1
                                        *  rot_sign_2 * double(1)/double(2) * sqrt( double(angular_momentum_J * (angular_momentum_J + 1 ) - rotation_quantum_m * (rotation_quantum_m + 1 ) ) );
                    }
                    if (Coriolis_coupling_rot_axis[rot_term_index] == 1 and rot_diff == 1){
                        // Jy
                        rot_sign_2 = complex<double>(1, -1) ; // -1j
                        rotation_quantum_m = min( (*vmode_ptr)[i][nmodes[m]] ,  (*vmode_ptr)[j][nmodes[m]] ) ;
                        value = value +  (-2 * rotational_constant[1] ) * Coriolis_coupling_strength[rot_term_index] * sqrt( double(v_k * v_l) ) /2 * rot_sign_1 *
                                         rot_sign_2 * double(1)/double(2) * sqrt( double(angular_momentum_J * (angular_momentum_J + 1 ) - rotation_quantum_m * (rotation_quantum_m + 1 ) ) );
                    }
                    if(Coriolis_coupling_rot_axis[rot_term_index] == 2 and rot_diff == 0){
                        //Jz term
                        rot_sign_2 = 1;
                        rotation_quantum_m = (*vmode_ptr)[i][nmodes[m]] ;
                        value = value +  (-2*rotational_constant[2]) * Coriolis_coupling_strength[rot_term_index] * sqrt( double(v_k * v_l) ) /2 * rot_sign_1 *
                                         rot_sign_2 * double(rotation_quantum_m);
                    }

                }

            }

        }

        // Rotational energy P *B * P : here P is vibrational rotation angular momentum, B is rotational constant
        if(rot_diff == 0 and (ntot == 4 or ntot == 2 )   ){
            for(rot_term_index = 0; rot_term_index < Coriolis_coupling_term_num; rot_term_index ++ ){
                for(rot_term_index2 = 0; rot_term_index2 < Coriolis_coupling_term_num; rot_term_index2 ++ ){
                    // here rotational constant B is diagonal.
                    if(Coriolis_coupling_rot_axis[rot_term_index] == Coriolis_coupling_rot_axis[rot_term_index2]){

                        coriolis_vib_term1 = Coriolis_coupling_vib_index[rot_term_index][0];
                        coriolis_vib_term2 = Coriolis_coupling_vib_index[rot_term_index][1];
                        coriolis_vib_term3 = Coriolis_coupling_vib_index[rot_term_index2][0];
                        coriolis_vib_term4 = Coriolis_coupling_vib_index[rot_term_index2][1];

                        for(i_1 = 0; i_1 < 2; i_1 ++){ // index 0 means a, index 1 means a^{+}
                            for(i_2 = 0; i_2 <2 ; i_2 ++ ){
                                if(i_2 == 0)  rot_sign_1 = 1; else rot_sign_1 = -1;

                                for(i_3 = 0; i_3 <2 ; i_3 ++ ){
                                    for(i_4 = 0; i_4 < 2; i_4 ++ ){
                                        if(i_4 == 0)  rot_sign_2 = 1; else rot_sign_2 = -1; // for a^{+} , we have - sign.

                                        quantum_num_new = (*vmode_ptr)[j];

                                        coeff_1 = sqrt ( quantum_num_new[coriolis_vib_term1] + i_1 ); // raising and lowering quanta term
                                        quantum_num_new[coriolis_vib_term1] = quantum_num_new[coriolis_vib_term1]  + quanta_move [i_1];
                                        coeff_2 = sqrt(quantum_num_new[coriolis_vib_term2] + i_2) ;
                                        quantum_num_new[coriolis_vib_term2] = quantum_num_new[coriolis_vib_term2] + quanta_move [i_2];
                                        coeff_3 = sqrt(quantum_num_new[coriolis_vib_term3] + i_3);
                                        quantum_num_new[coriolis_vib_term3] = quantum_num_new[coriolis_vib_term3] + quanta_move [i_3];
                                        coeff_4 = sqrt(quantum_num_new[coriolis_vib_term4] + i_4) ;
                                        quantum_num_new[coriolis_vib_term4] = quantum_num_new[coriolis_vib_term4] + quanta_move [i_4];

                                        if(quantum_num_new == (*vmode_ptr)[i]){
                                            value = value + (-1/double(4)) * rot_sign_1 * rot_sign_2 * coeff_1 * coeff_2 * coeff_3 * coeff_4 *
                                                            rotational_constant[ Coriolis_coupling_rot_axis[rot_term_index] ] * Coriolis_coupling_strength[rot_term_index] * Coriolis_coupling_strength[rot_term_index2];
                                        }

                                    }
                                }
                            }
                        }


                    }

                }
            }
        }

    }

    return value;

}