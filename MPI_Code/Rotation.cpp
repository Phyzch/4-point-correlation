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

void detector:: compute_rotational_offdiag_part_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1){
    int i,j,m,k;
    int begin_index;
    int ntot;
    double value, lij;
    vector<vector <int>> * vmode_ptr;
    vector<double> * dmat_ptr;
    bool exist;
    int position;


    int higher_rot_quanta;
    int rot_diff;
    for(m=0; m< stlnum; m++ ){
        if(m==0){
            vmode_ptr = &(vmode0);
            dmat_ptr= &(dmat0);
        }
        else {
            vmode_ptr= &(vmode1);
            dmat_ptr= &(dmat1);
        }
        begin_index= total_dmat_size[m]/num_proc * my_id;
        // compute off diagonal matrix :
        for(i=begin_index ; i< begin_index + dmatsize[m] ; i++ ){
            for(j=0;j<total_dmat_size[m]; j++ ){
                if(i == j ) continue;
                // check vib term diff for vibrational term
                value = 0;
                rot_diff = abs ((*vmode_ptr)[i][nmodes[m]] - (*vmode_ptr)[j][nmodes[m]]);

                ntot = 0;
                for(k=0;k<nmodes[m];k++){
                    deln[k]= abs( (*vmode_ptr)[i][k] - (*vmode_ptr)[j][k] ); // same as deln[k] = abs(dv[m][i][k] - dv[m][j][k]);
                    ntot = ntot + deln[k] ;
                }

                if(ntot == 0){
                    if( (*vmode_ptr)[i][nmodes[m]] == (*vmode_ptr)[j][nmodes[m]] ){
                        log << "Error! The off-diagonal element must differ in q.n." << endl;
                        MPI_Abort(MPI_COMM_WORLD,-8);
                    }
                    // J_x ^ 2 , J_y^2  coupling term.
                    if(  rot_diff == 2  ){
                        higher_rot_quanta = max((*vmode_ptr)[i][nmodes[m]] , (*vmode_ptr)[j][nmodes[m]]);
                        value = (rotational_constant[0] + rotational_constant[1]) * 1/4 * pow (
                                (angular_momentum_J * (angular_momentum_J + 1 ) - higher_rot_quanta * (higher_rot_quanta - 1) )  * ( angular_momentum_J * (angular_momentum_J + 1 )  -  (higher_rot_quanta-2) * (higher_rot_quanta - 1) )
                                , 1/2) ;
                    }

                }
                else{
                    // Coriolis coupling
                    if(ntot == 2 and rot_diff == 1 ){

                    }


                }


            }
        }

    }
}