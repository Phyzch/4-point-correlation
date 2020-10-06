//
// Created by phyzch on 7/19/20.
//
#include"../util.h"
#include"../system.h"
using namespace std;

// function is almost the same for the one Detector_Evolve_MPI
int full_system::construct_recvbuffer_index(){
    int i,j,k;
    int total_recv_count=0;
    int local_begin= matsize_offset_each_process[my_id];
    int local_end;
    if(my_id!=num_proc-1){
        local_end = matsize_offset_each_process[my_id+1];
    }
    else{
        local_end= total_matsize;
    }

    // --------- Allocate space for vector to receive -----------------------
    remoteVecCount = new int [num_proc];
    remoteVecPtr= new int [num_proc];
    remoteVecIndex = new int [offnum];
    for(i=0;i<num_proc;i++){
        remoteVecCount[i] = 0;
    }
    // ----------------------------------------------------------------


    vector<int> icol_copy = icol;
    sort(icol_copy.begin(),icol_copy.end());
    j=0;
    int prev_col=-1;
    int col_index;
    int col_pc_id;
    for(i=0;i<matnum;i++){
        col_index= icol_copy[i];
        if(col_index>prev_col and (col_index<local_begin or col_index >= local_end) ){
            // find process id this icol belong to
            for(k=1;k<num_proc;k++){
                if(col_index < matsize_offset_each_process[k]){
                    break;
                }
            }
            col_pc_id= k-1;
            remoteVecCount[col_pc_id]++;
            remoteVecIndex[j]= col_index;
            j++;
        }
        prev_col= col_index;
    }

    remoteVecPtr[0]=0;
    for(i=1;i<num_proc;i++){
        remoteVecPtr[i]= remoteVecPtr[i-1] + remoteVecCount[i-1];
    }
    total_recv_count=0;
    for(i=0;i<num_proc;i++){
        total_recv_count = total_recv_count + remoteVecCount[i];
    }
    return total_recv_count;
}


void full_system:: prepare_evolution(){

    int m,i;
    int col_index_for_search;
    int * search_Ind;
    // You have to allocate space for tosendVecCount, tosendVecPtr.
    tosendVecCount = new int [num_proc];
    tosendVecPtr = new int [num_proc];
    to_recv_buffer_len = construct_recvbuffer_index();
    to_send_buffer_len= construct_send_buffer_index(remoteVecCount,remoteVecPtr,remoteVecIndex,
                                                    tosendVecCount,tosendVecPtr,tosendVecIndex);
    x.resize(matsize + to_recv_buffer_len);
    y.resize(matsize+ to_recv_buffer_len);
    recv_x = new double [to_recv_buffer_len];
    recv_y= new double [to_recv_buffer_len];
    send_x = new double[to_send_buffer_len];
    send_y= new double [to_send_buffer_len];
    local_irow.reserve(matnum);
    local_icol.reserve(matnum);
    for(i=0;i<matnum;i++){
        local_irow.push_back(irow[i] - matsize_offset_each_process[my_id]);
        col_index_for_search = icol[i];
        search_Ind = (int *)bsearch(&col_index_for_search,remoteVecIndex,to_recv_buffer_len,sizeof(int),compar);
        if (search_Ind!=NULL){
            // this column element is in remote vector:
            local_icol.push_back(matsize + (search_Ind - remoteVecIndex));
        }
        else{
            // this column element is in local vector:
            local_icol.push_back(icol[i] - matsize_offset_each_process[my_id]);
        }
    }
    sleep(1);
}

void full_system:: full_system_SUR_one_step()
{
    int i;
    int irow_index, icol_index;
    update_x_y();
    // SUR_algorithm
    for(i=0;i<matnum;i++){
        irow_index = local_irow[i];
        icol_index= local_icol[i];
        x[irow_index] = x[irow_index] + mat[i] * y[icol_index];
    }
    for(i=0;i<matnum;i++){
        irow_index=local_irow[i];
        icol_index= local_icol[i];
        y[irow_index] = y[irow_index] - mat[i] * x[icol_index];
    }
}

void full_system:: update_x_y(){
    int i;
    int local_index;
    // update tosend_xd for sending  to other process
    for(i=0;i<to_send_buffer_len;i++){
        local_index = tosendVecIndex[i] - matsize_offset_each_process[my_id];
        send_x[i] = x[local_index];
        send_y[i] = y[local_index];
    }

    MPI_Alltoallv(&send_x[0],tosendVecCount,tosendVecPtr,MPI_DOUBLE,
                  &recv_x[0],remoteVecCount,remoteVecPtr,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Alltoallv(&send_y[0],tosendVecCount,tosendVecPtr,MPI_DOUBLE,
                  &recv_y[0],remoteVecCount,remoteVecPtr,MPI_DOUBLE,MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len;i++){
        x[matsize+i] = recv_x[i];
        y[matsize+i] = recv_y[i];
    }
}

void full_system:: update_x(){
    int i;
    int local_index;
    for(i=0;i<to_send_buffer_len;i++){
        local_index = tosendVecIndex[i] - matsize_offset_each_process[my_id];
        send_x[i] = x[local_index];
    }

    MPI_Alltoallv(&send_x[0],tosendVecCount,tosendVecPtr,MPI_DOUBLE,
                  &recv_x[0],remoteVecCount,remoteVecPtr,MPI_DOUBLE,MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len;i++){
        x[matsize+i] = recv_x[i];
    }
}

void full_system::update_y(){
    int i;
    int local_index;
    // update tosend_xd for sending  to other process
    for(i=0;i<to_send_buffer_len;i++){
        local_index = tosendVecIndex[i] - matsize_offset_each_process[my_id];
        send_y[i] = y[local_index];
    }

    MPI_Alltoallv(&send_y[0],tosendVecCount,tosendVecPtr,MPI_DOUBLE,
                  &recv_y[0],remoteVecCount,remoteVecPtr,MPI_DOUBLE,MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len;i++){
        y[matsize+i] = recv_y[i];
    }
}

void full_system::Normalize_wave_function(){
    int i;
    norm = 0;
    total_norm = 0;
    for(i=0;i<matsize;i++) {
        norm = norm + x[i] * x[i] + y[i] * y[i];
    }
    MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<matsize;i++){
        x[i] = x[i] / sqrt(total_norm);
        y[i] = y[i] / sqrt(total_norm);
    }
}