//
// Created by phyzch on 7/5/21.
//
#include "../util.h"
#include "../system.h"

void detector::prepare_ladder_operation(){
    int j, k;
    send_xd_ladder_operator = new double * [2 * nmodes[0] ];
    send_yd_ladder_operator = new double * [2 * nmodes[0] ];
    recv_xd_ladder_operator = new double * [2 * nmodes[0] ];
    recv_yd_ladder_operator = new double * [2 * nmodes[0] ];
    for(j=0;j<2*nmodes[0]; j++){
        send_xd_ladder_operator[j] = new double [dmatsize[0]];
        send_yd_ladder_operator[j] = new double [dmatsize[0]];
        recv_xd_ladder_operator[j] = new double [dmatsize[0]];
        recv_yd_ladder_operator[j] = new double [dmatsize[0]];
    }
}

void detector::ladder_operator_operation( const vector<double> & wave_func_x , const vector<double> & wave_func_y ,
                                          vector<vector<double>> & xd_for_ladder_operator ,
                                          vector<vector<double>> & yd_for_ladder_operator
                                          ){
    // state after operate ladder opearator. We already take \sqrt{n_{j}} coefficient into account here.
    // Before calling this function, we have to make sure a_{i} | \phi> is incorporated in basis set.
    int i, j ,k;
    int local_state_index;
    int begin_index = dmatsize_offset_each_process[0][my_id];
    for(j=0; j< 2 * nmodes[0]; j++ ){
        // j is index for move in quantum number space for mode j.
        for(k=0;k<to_send_buffer_len_list_with_ladder_operator[j]; k++ ){
            local_state_index = tosendVecIndex_for_xp[j][k] - begin_index;  //tosendVecIndex_for_xp store index for basis set after ladder opeartor operate upon state
            send_xd_ladder_operator[j][k] = wave_func_x[local_state_index];
            send_yd_ladder_operator[j][k] = wave_func_y[local_state_index];
        }
    }

    for(j=0;j<2*nmodes[0]; j++ ){
        MPI_Alltoallv( &send_xd_ladder_operator[j][0] , tosendVecCount_for_xp[j] , tosendVecPtr_for_xp[j] , MPI_DOUBLE,
                       &recv_xd_ladder_operator[j][0] , remoteVecCount_for_xp[j] , remoteVecPtr_for_xp[j] , MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(&send_yd_ladder_operator[j][0] , tosendVecCount_for_xp[j] , tosendVecPtr_for_xp[j] , MPI_DOUBLE,
                      &recv_yd_ladder_operator[j][0] , remoteVecCount_for_xp[j] , remoteVecPtr_for_xp[j] , MPI_DOUBLE, MPI_COMM_WORLD );
    }

    // now put variable received into right place for wave function.
    for(j=0;j<2*nmodes[0];j++){
        vector<double> v (dmatsize[0],0);
        xd_for_ladder_operator.push_back(v);
        yd_for_ladder_operator.push_back(v);
    }

    // a_{j} | \phi>  lowering operator.   <{n} | a_{j} | \phi> = \sqrt{n_{j} + 1 } <{n}_{j} + 1| \phi>
    for(j=0;j<nmodes[0];j++){

        for(k=0;k<dmatsize[0];k++){

            if(Index_in_remoteVecIndex_for_xp[j + nmodes[0]][k] == -1){
                // this state do not exist, we assign xd,yd = 0 in this case
                xd_for_ladder_operator[j][k] = 0;
                yd_for_ladder_operator [j][k] = 0;
            }
            else{
                // add \sqrt{n} to coefficient.
                xd_for_ladder_operator[j][k] = recv_xd_ladder_operator [j + nmodes[0] ][Index_in_remoteVecIndex_for_xp[j + nmodes[0]  ][k]]  * sqrt(dv[0][k][j]  + 1 );
                yd_for_ladder_operator [j][k] = recv_yd_ladder_operator [j + nmodes[0] ][Index_in_remoteVecIndex_for_xp[j + nmodes[0] ][k]]  * sqrt(dv[0][k][j]  + 1 ) ;
            }

        }
    }
    // a_{j}^{+} raising operator.  <{n}| a_{j}^{+} | \phi> =\sqrt{n}_{j}  <{n}_{j}-1 | \phi>
    for(j=nmodes[0]; j< 2*nmodes[0]; j++){

        for(k=0;k<dmatsize[0];k++){
            if(Index_in_remoteVecIndex_for_xp[j - nmodes[0]][k] == -1){
                // this state do not exist, we assign xd,yd = 0 in this case
                xd_for_ladder_operator[j][k] = 0;
                yd_for_ladder_operator [j][k] = 0;
            }
            else{
                xd_for_ladder_operator[j][k] = recv_xd_ladder_operator[ j - nmodes[0] ] [Index_in_remoteVecIndex_for_xp [j - nmodes[0]][k] ] * sqrt(dv[0][k][j - nmodes[0]]);
                yd_for_ladder_operator[j][k] = recv_yd_ladder_operator[ j - nmodes[0] ] [Index_in_remoteVecIndex_for_xp [j - nmodes[0]][k] ] * sqrt(dv[0][k][j - nmodes[0]]);
            }

        }

    }



}

