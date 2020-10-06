//
// Created by phyzch on 6/30/20.
//
# include"../system.h"
#include "../util.h"
using namespace std;
double inter_detector_coupling_scaling=10;

void full_system::compute_offdiagonal_part_MPI(){
    vector<double> d_off_mat;
    vector<int> d_off_irow;
    vector<int> d_off_icol;
    compute_dmat_off_diagonal_matrix_in_full_matrix_MPI(d_off_mat,d_off_irow,d_off_icol);

    vector < double > * sys_detector_mat = new vector <double> [2];
    vector  <int> * sys_detector_irow = new vector <int> [2];
    vector<int> * sys_detector_icol = new vector <int> [2];
    compute_sys_detector_coupling_MPI(sys_detector_mat,sys_detector_irow,sys_detector_icol);

    vector <double>  d_d_mat;
    vector<int>  d_d_irow;
    vector<int>  d_d_icol;
    compute_detector_detector_coupling_MPI(d_d_mat,d_d_irow,d_d_icol);  // unfortunately we have to share slist among all process to find out detector detector coupling.

    // we have to rearrange off-diagonal_matrix in full_system to make sure irow is in  corresponding process.
    //Also we have to recompute offnum, matnum
    combine_offdiagonal_term(sys_detector_mat,sys_detector_irow,sys_detector_icol,
                             d_off_mat,d_off_irow,d_off_icol,d_d_mat,d_d_irow,d_d_icol);

    offnum= matnum-matsize;
    // compute total_matnum, total_offnum, matnum_each_process, offnum_each_process.
    mat_offnum_each_process= new int [num_proc];
    matnum_each_process= new int [num_proc];
    matnum_offset_each_process = new int [num_proc];
    MPI_Allgather(&matnum,1,MPI_INT,&matnum_each_process[0],1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&offnum,1,MPI_INT,&mat_offnum_each_process[0],1,MPI_INT,MPI_COMM_WORLD);
    int i;
    matnum_offset_each_process[0]=0;
    for(i=1;i<num_proc;i++){
        matnum_offset_each_process[i]= matnum_offset_each_process[i-1] + matnum_each_process[i-1];
    }

    MPI_Allreduce(&matnum,&total_matnum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&offnum,&total_offnum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    delete [] sys_detector_mat;
    delete [] sys_detector_irow;
    delete [] sys_detector_icol;
}

void full_system:: combine_offdiagonal_term(vector <double> * sys_detector_mat, vector<int> * sys_detector_irow, vector<int> * sys_detector_icol,
        vector<double> & d_off_mat, vector<int> & d_off_irow, vector<int> & d_off_icol,
        vector<double> & d_d_mat, vector<int> & d_d_irow, vector<int> & d_d_icol){
    // combine 3 part of off-diagonal term together and construct sdindex, sdmode, sdnum
    int m,i;
    matnum=matsize;
    int size;
    for(m=0;m<s.tlnum;m++){
        matnum = matnum + sys_detector_mat[m].size();
    }
    matnum = matnum + d_off_mat.size() + d_d_mat.size();
    mat.reserve(matnum);
    irow.reserve(matnum);
    icol.reserve(matnum);
    int matindex=matsize;
    for(m=0;m<s.tlnum;m++){
        size = sys_detector_mat[m].size();
        for(i=0;i<size;i++) {
            mat.push_back(sys_detector_mat[m][i]);
            irow.push_back(sys_detector_irow[m][i]);
            icol.push_back(sys_detector_icol[m][i]);

            sdindex[m].push_back(matindex);  // this sdindex is local index in process.
            sdmode[m].push_back(0);
            sdnum[m]++;
            matindex++;
        }
    }
    for(m=0;m<s.tlnum;m++){
        MPI_Allgather(&sdnum[m],1,MPI_INT,&sdnum_each_process[m][0],1,MPI_INT,MPI_COMM_WORLD);
        total_sd_num[m] = 0;
        for(i=0;i<num_proc;i++){
            total_sd_num[m] = total_sd_num[m] + sdnum_each_process[m][i] ;
        }
        sdnum_displacement_each_process[m][0] = 0;
        for(i=1;i<num_proc;i++){
            sdnum_displacement_each_process[m][i] = sdnum_displacement_each_process[m][i-1] + sdnum_each_process[m][i-1];
        }
    }

    size= d_off_mat.size();
    for(i=0;i<size;i++){
        mat.push_back(d_off_mat[i]);
        irow.push_back(d_off_irow[i]);
        icol.push_back(d_off_icol[i]);
        matindex++;
    }
    size = d_d_mat.size();
    for(i=0;i<size;i++){
        mat.push_back(d_d_mat[i]);
        irow.push_back(d_d_irow[i]);
        icol.push_back(d_d_icol[i]);
        matindex++;
    }
}

void full_system::compute_dmat_off_diagonal_matrix_in_full_matrix_MPI(vector < double > & d_off_mat,vector  <int> & d_off_irow, vector<int> & d_off_icol){
    for(int i=0;i<s.tlnum;i++){
        compute_dmat_off_diagonal_matrix_in_full_matrix_one_dmat_MPI(i,d_off_mat,d_off_irow,d_off_icol);
    }
    rearrange_off_diagonal_term(d_off_mat,d_off_irow,d_off_icol);
}
void full_system::compute_dmat_off_diagonal_matrix_in_full_matrix_one_dmat_MPI(int index,vector < double > & d_off_mat,vector  <int> & d_off_irow, vector<int> & d_off_icol){
    // we just use q_index easily compute the result.
    vector<quotient_state> * dlist_ptr;
    quotient_state * d_ptr;
    if(index ==0){
        dlist_ptr = &(d1list);
    }
    else{
        dlist_ptr = &(d2list);
    }
    int i,j,k,l,p;
    int m,n;
    int size= (*dlist_ptr).size();
    int q_index_list_size;
    for(m=0;m<size;m++){
        d_ptr =&((*dlist_ptr)[m]);
        q_index_list_size = (*d_ptr).q_index_list.size();
        for(n=0;n<q_index_list_size;n++){
            i= (*d_ptr).q_index_list[n][0];   // index in dmat  (diag)
            j = (*d_ptr).q_index_list[n][1];  // index in dmat  (diag)
            k=  (*d_ptr).q_index_list[n][2];  // index in mat (diagonal)
            l= (*d_ptr).q_index_list[n][3];   // index in mat (diagonal)
            p= (*d_ptr).q_index_list[n][4];   // index in dmat (off_diag)  p in dmat link i,j. or dirow[p] = i, dicol[p]=j
            if(i!=j){
                // this is off-diagonal element
                d_off_mat.push_back((*d_ptr).dmat_value_list[n]);
                d_off_irow.push_back(k);
                d_off_icol.push_back(l);

                // we have to count a_{ji} to enable parallel version of SUR algorithm
                d_off_mat.push_back( (*d_ptr).dmat_value_list[n] );
                d_off_irow.push_back(l);
                d_off_icol.push_back(k);
            }
        }
    }
}

void full_system::compute_sys_detector_coupling_MPI(vector < double > * sys_detector_mat, vector  <int> * sys_detector_irow, vector<int> * sys_detector_icol){
    sdnum[0]=0;
    sdnum[1]=0;
    int i,j,l;
    int k,m;
    int slist_size= slist.size();
    int sysx_index_size;
    int sysx_index_size2;
    vector <int> * coup_vmode1_ptr;  // pointer for vmode1
    vector <int> * coup_vmode2_ptr;  // pointer for vmode2
    vector <int> coup_vmode1_near; // adjacent vmode 1
    vector <int> coup_vmode2_near; // adjacent vmode 2
    int position1;
    int position2;
    bool exist1;
    bool exist2;
    //------------- for communication betweeen processes (serialization) --------------------------------
    vector<vector<int>> * vmode_tosend_list = new vector<vector<int>> [2];
    //--------------------------------------------------------------------
    for(i=0;i<slist_size;i++){
        sysx_index_size= slist[i].sysxindex.size();
        for(j=0;j<sysx_index_size;j++){
            if(slist[i].sysxindex[j]==0){// we found |00000>
                // check coupling between detector 0 and photon:  |00>  <-> |10>
                if(slist[i].vmode1[0] > 0 and d.modtype[0][0] == 1){
                    coup_vmode1_ptr= &(slist[i].vmode1);
                    coup_vmode1_near.clear();
                    coup_vmode1_near.push_back((*coup_vmode1_ptr)[0] -1); // vmode1[0] -1
                    for(k=1;k<d.nmodes[0];k++){
                        coup_vmode1_near.push_back((*coup_vmode1_ptr)[k]);
                    }
                    // search in detector subspace . May cause confusion, vmode0 actually is vmode 1 here.
                    find_position_for_insert_binary(vmode0,coup_vmode1_near,exist1);
                    if(exist1){
                        vector <int> vmode_to_send;
                        vmode_to_send.insert(vmode_to_send.end(),coup_vmode1_near.begin(),coup_vmode1_near.end());
                        vmode_to_send.insert(vmode_to_send.end(),slist[i].vmode2.begin(), slist[i].vmode2.end());
                        vmode_to_send.push_back(slist[i].xindex[j]);
                        vmode_tosend_list[0].push_back(vmode_to_send);
                    }
                }
                if(slist[i].vmode2[0] > 0 and d.modtype[1][0] == 1){
                    // check coupling between detector 0 and photon:  |00>  <-> |01>
                    coup_vmode2_ptr= &(slist[i].vmode2);
                    coup_vmode2_near.clear();
                    coup_vmode2_near.push_back((*coup_vmode2_ptr)[0] -1);
                    for(k=1;k<d.nmodes[1];k++){
                        coup_vmode2_near.push_back((*coup_vmode2_ptr)[k]);
                    }
                    // search in detector subspace. (May  cause confusion: vmode1 actually is vmode2 in second detector here)
                    find_position_for_insert_binary(vmode1,coup_vmode2_near,exist1);
                    if(exist1){
                        vector<int> vmode_to_send;
                        vmode_to_send.insert(vmode_to_send.end(),slist[i].vmode1.begin(),slist[i].vmode1.end());
                        vmode_to_send.insert(vmode_to_send.end(),coup_vmode2_near.begin(),coup_vmode2_near.end());
                        vmode_to_send.push_back(slist[i].xindex[j]);
                        vmode_tosend_list[1].push_back(vmode_to_send);
                    }
                }
            }
        }
    }
    //------------- for communication betweeen processes (serialization) --------------------------------
    int ** vmode_list_1d = new int * [2];  // 1d version of vmode
    int list_size;
    int list_size_1d;
    int index;
    for(i=0;i<2;i++){
        list_size = vmode_tosend_list[i].size();
        list_size_1d = list_size * (d.nmodes[0] + d.nmodes[1] + 1); // record 2 vmode + xindex.
        vmode_list_1d[i] = new int [list_size_1d];
        index=0;
        for(j=0;j<list_size;j++){
            for(k=0;k< d.nmodes[0] + d.nmodes[1] +1;k++){
                vmode_list_1d[i][index] = vmode_tosend_list[i][j][k];
                index ++;
            }
        }
    }
    // ----- tell each process how many (vmode1,vmode2,xindex) tuple they have to receive from others.
    int * vmode_send_amount = new int [2];
    for(i=0;i<2;i++){
        vmode_send_amount[i]= vmode_tosend_list[i].size();
    }
    int * vmode_send_data_amount = new int [2];  // amount of int data process will receive.
    for(i=0;i<2;i++){
        vmode_send_data_amount[i] = vmode_send_amount[i] * (d.nmodes[0] + d.nmodes[1] +1);
    }
    // --------place to receive vmode_send_amount ----------------
    int ** vmode_send_amount_each_process = new int * [2];
    for(i=0;i<2;i++){
        vmode_send_amount_each_process[i] = new int [num_proc];
    }
    for(i=0;i<2;i++){
        MPI_Allgather(&vmode_send_amount[i],1,MPI_INT,
                      &vmode_send_amount_each_process[i][0],1,MPI_INT,MPI_COMM_WORLD);
    }
    //---------------- compute amount of data to receive and its displacement for Allgatherv API to send data to all process --------------------------
    int ** data_amount_each_process = new int * [2];
    int ** data_amount_displacement_each_process= new int *[2];
    for(i=0;i<2;i++){
        data_amount_each_process[i] = new int [num_proc];
        for(j=0;j<num_proc;j++){
            data_amount_each_process[i][j] = vmode_send_amount_each_process[i][j] * ( d.nmodes[0] + d.nmodes[1] +1);
        }
    }
    for(i=0;i<2;i++){
        data_amount_displacement_each_process[i] = new int [num_proc];
        data_amount_displacement_each_process[i][0]=0;
        for(j=1;j<num_proc;j++){
            data_amount_displacement_each_process[i][j] = data_amount_displacement_each_process[i][j-1] + data_amount_each_process[i][j-1];
        }
    }
    // ------------- total amount of (vmode1,vmode2,xindex) data they will receive-----------------------------
    int * total_vmode_to_receive = new int [2];
    for(i=0;i<2;i++){
        total_vmode_to_receive[i]=0;
        for(j=0;j<num_proc;j++){
            total_vmode_to_receive[i]= total_vmode_to_receive[i] + vmode_send_amount_each_process[i][j];
        }
    }
    // ------------- receive data from vmode_list_1d to store in vmode_list_recv_1d
    int ** vmode_list_recev_1d = new int * [2];
    for(i=0;i<2;i++){
        vmode_list_recev_1d[i]= new int [total_vmode_to_receive[i] * (d.nmodes[0] + d.nmodes[1] +1)];
    }

    //------------------- send data to vmode_list_recv_1d --------------------------
    for(i=0;i<2;i++) {
        MPI_Allgatherv(&vmode_list_1d[i][0], vmode_send_data_amount[i],MPI_INT,
                       &vmode_list_recev_1d[i][0],data_amount_each_process[i],data_amount_displacement_each_process[i],MPI_INT,MPI_COMM_WORLD);
    }
    // ----------- convert data in vmode_list_recv_1d to sysindex list and vmode1, vmode2 list -------------------------
    vector<vector<int>> * vmode1_recv = new vector<vector<int>> [2];
    vector<vector<int>> * vmode2_recv = new vector<vector<int>> [2];
    vector<int> * sys_xindex_recv = new vector<int> [2];
    for(i=0;i<2;i++){
        index=0;
        for(j=0;j<total_vmode_to_receive[i];j++){
            vector<int> vmode1_temp;
            vector<int> vmode2_temp;
            for(k=0;k<d.nmodes[0];k++){
                vmode1_temp.push_back(vmode_list_recev_1d[i][index]);
                index++;
            }
            for(k=0;k<d.nmodes[1];k++){
                vmode2_temp.push_back(vmode_list_recev_1d[i][index]);
                index++;
            }
            sys_xindex_recv[i].push_back(vmode_list_recev_1d[i][index]);
            index++;
            vmode1_recv[i].push_back(vmode1_temp);
            vmode2_recv[i].push_back(vmode2_temp);
        }
    }

    // ----------------- search these data in our sys_quotient_state_list --------------
    for(i=0;i<2;i++) {
        for (j = 0; j < total_vmode_to_receive[i]; j++) {
            l = find_location_binarysearch_sys_quotient_state(slist, vmode1_recv[i][j], vmode2_recv[i][j], exist2);
            if (exist2) {
                sysx_index_size2 = slist[l].sysxindex.size();
                for (k = 0; k < sysx_index_size2; k++) {
                    if (i==0){
                        if(slist[l].sysxindex[k]==1){  // coupling between |00> <-> |10>
                            sys_detector_irow[0].push_back(sys_xindex_recv[i][j]);  // push xindex into irow
                            sys_detector_icol[0].push_back(slist[l].xindex[k]);  // push xindex into jrow
                            sys_detector_mat[0].push_back(d.premodcoup[0][0]);

                            // record the transpose element here for Parallel version of SUR algorithm
                            sys_detector_irow[0].push_back(slist[l].xindex[k]);
                            sys_detector_icol[0].push_back(sys_xindex_recv[i][j]);
                            sys_detector_mat[0].push_back(d.premodcoup[0][0]);
                        }
                    }
                    else if (i==1){
                        if(slist[l].sysxindex[k]==2){  // coupling between |00> <-> |01>
                            sys_detector_irow[1].push_back(sys_xindex_recv[i][j]);
                            sys_detector_icol[1].push_back(slist[l].xindex[k]);
                            sys_detector_mat[1].push_back(d.premodcoup[1][0]);

                            // record the transpose element here for Parallel version of SUR algorithm.
                            sys_detector_irow[1].push_back(slist[l].xindex[k]);
                            sys_detector_icol[1].push_back(sys_xindex_recv[i][j]);
                            sys_detector_mat[1].push_back(d.premodcoup[1][0]);
                        }
                    }
                }
            }
        }
    }
    for(m=0;m<s.tlnum;m++){
        rearrange_off_diagonal_term(sys_detector_mat[m],sys_detector_irow[m],sys_detector_icol[m]);
    }
    // -----------------free the space ----------------------------
    for(i=0;i<2;i++){
        delete [] vmode_list_1d[i];
        delete [] vmode_send_amount_each_process[i];
        delete [] data_amount_each_process[i];
        delete [] data_amount_displacement_each_process[i];
        delete [] vmode_list_recev_1d[i];
    }
    delete [] vmode_tosend_list;
    delete [] vmode_list_1d;
    delete [] vmode_send_amount;
    delete [] vmode_send_data_amount;
    delete [] vmode_send_amount_each_process;
    delete [] data_amount_each_process;
    delete [] data_amount_displacement_each_process;
    delete [] total_vmode_to_receive;
    delete [] vmode_list_recev_1d;
    delete [] vmode1_recv;
    delete [] vmode2_recv;
    delete [] sys_xindex_recv;
}



 void full_system::compute_detector_detector_coupling_MPI(vector <double> & d_d_mat, vector<int> & d_d_irow, vector<int> & d_d_icol){
    int detector_detector_coupling_term=0;
    int total_detector_detector_coupling_term=0;
    int i,j,k,l,m;
    int n;
    int a,b;
    vector<int> direction = {1, -1};
    int sys_list_size = slist.size();
    int sys_xindex_size1, sys_xindex_size2;
    int sys_xindex_1, sys_xindex_2;
    vector<int> vmode1_near;
    vector<int> vmode2_near;
    bool exist1, exist2,exist3;
    double value, lij;
    vector <vector<int>> vmode_to_send_list; // [vmode1, vmode2, sysindex, xinidex]
    vector<double> d_d_coupling_strength_list;

    int first_pass=0;

    for(i=0;i<sys_list_size;i++){
        for(k=0;k< d.nmodes[0] ;k++){
            for(l=0;l< d.nmodes[1] ;l++){
                if(d.modtype[0][k]!=0 or d.modtype[1][l]!=0) continue;
                value= d.modcoup[0][k] * d.modcoup[1][l] / inter_detector_coupling_scaling;
                if (inter_detector_coupling) {
                    value = value * (1+ inter_detector_coupling_noise * 2* (double(rand()) / RAND_MAX -0.5 )) ;
                }
                if(d.mfreq[0][k]!=d.mfreq[1][l]){
                    lij = abs(value) / abs(d.mfreq[0][k] - d.mfreq[1][l]);
                }
                else{
                    if(value!=0){
                        lij=10000;
                    }
                    else {
                        lij = 0;
                    }
                }
                if(lij>d.cutoff2){
                    for(m=0;m<2;m++){
                        vmode1_near.clear();
                        vmode1_near= slist[i].vmode1;
                        vmode1_near[k] = vmode1_near[k] + direction[m];
                        if (vmode1_near[k] < 0 or vmode1_near[k]> d.nmax[0][k]) continue;
                        find_position_for_insert_binary(vmode0,vmode1_near,exist1);  // confusion may come here, vmode0 is actually first detector, vmode1_near may be also in first detector.
                        if(exist1){
                            vmode2_near.clear();
                            vmode2_near= slist[i].vmode2;
                            vmode2_near[l] = vmode2_near[l] - direction[m];
                            if (vmode2_near[l] < 0 or vmode2_near[l] > d.nmax[1][l]) continue;
                            find_position_for_insert_binary(vmode1, vmode2_near, exist2);  // confusion may come here: vmode1 is actually second detector, vmode2_near is also mode in 2rd detector.
                            if(exist2){
                                first_pass++;
                                // we have to record sysxindex, xindex here.
                                sys_xindex_size1= slist[i].sysxindex.size();
                                for(a=0;a<sys_xindex_size1;a++){
                                    vector<int> vmode_to_send;  // form : [vmode1, vmode2, sysxindex, xindex]
                                    vmode_to_send.insert(vmode_to_send.end(),vmode1_near.begin(),vmode1_near.end());
                                    vmode_to_send.insert(vmode_to_send.end(),vmode2_near.begin(),vmode2_near.end());
                                    vmode_to_send.push_back(slist[i].sysxindex[a]);
                                    vmode_to_send.push_back(slist[i].xindex[a]);
                                    vmode_to_send_list.push_back(vmode_to_send);
                                    d_d_coupling_strength_list.push_back(value); // record value of coupling strength.
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    int * vmode_list_1d; // 1d array of vmode_to_send_list
    double * d_d_coupling_strength_to_send;
    int vmode_list_size; // size of vmode_to_send_list
    int vmode_list_size_1d; // size of vmode_list_1d
    vmode_list_size= vmode_to_send_list.size();
    vmode_list_size_1d= vmode_list_size * (d.nmodes[0] + d.nmodes[1] + 2);
    vmode_list_1d= new int [vmode_list_size_1d];
    int index=0;
    for(i=0;i<vmode_list_size;i++){
        for (j=0;j< d.nmodes[0] + d.nmodes[1] + 2;j++){
            vmode_list_1d[index] = vmode_to_send_list[i][j];
            index++;
        }
    }
    d_d_coupling_strength_to_send= new double [vmode_list_size];
    for(i=0;i<vmode_list_size;i++){
        d_d_coupling_strength_to_send[i]= d_d_coupling_strength_list[i];
    }
    // ------------------- record vmode list size gonna to receive from all process -------------------------
    int * vmode_list_size_each_process = new int [num_proc];
    int * vmode_list_displacement_each_process= new int [num_proc];
    int * vmode_list_size_1d_each_process= new int [num_proc];
    int * vmode_list_displacement_1d_each_process= new int [num_proc];
    MPI_Allgather(&vmode_list_size,1,MPI_INT,&vmode_list_size_each_process[0],1,MPI_INT,MPI_COMM_WORLD);

    vmode_list_displacement_each_process[0]=0;
    for(i=1;i<num_proc;i++){
        vmode_list_displacement_each_process[i] = vmode_list_displacement_each_process[i-1] + vmode_list_size_each_process[i-1];
    }

    for(i=0;i<num_proc;i++){
        vmode_list_size_1d_each_process[i]= vmode_list_size_each_process[i] * (d.nmodes[0] + d.nmodes[1] +2);
    }

    vmode_list_displacement_1d_each_process[0]=0;
    for(i=1;i<num_proc;i++){
        vmode_list_displacement_1d_each_process[i]= vmode_list_displacement_1d_each_process[i-1] + vmode_list_size_1d_each_process[i-1];
    }

     int total_vmode_list_size=0;
     int total_vmode_list_size_1d;
     for(i=0;i<num_proc;i++){
         total_vmode_list_size= total_vmode_list_size + vmode_list_size_each_process[i];
     }
     total_vmode_list_size_1d= total_vmode_list_size * ( d.nmodes[0] + d.nmodes[1] + 2);

    // -------------------- use MPI_Allgatherv to gather vmode_to_send_list and d_d_coupling_strength_list -------------
    int * vmode_list_to_recv =new int [total_vmode_list_size_1d];
    double * d_d_coupling_strength_list_recv= new double [total_vmode_list_size];

    MPI_Allgatherv(&vmode_list_1d[0],vmode_list_size_1d,MPI_INT,
                   &vmode_list_to_recv[0],vmode_list_size_1d_each_process,vmode_list_displacement_1d_each_process,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgatherv(&d_d_coupling_strength_to_send[0],vmode_list_size,MPI_DOUBLE,
                   &d_d_coupling_strength_list_recv[0],vmode_list_size_each_process,vmode_list_displacement_each_process,MPI_DOUBLE,MPI_COMM_WORLD);

    // ------- construct detector detecotr coupling -------------------------
    int vmode_element_index=0;
    int d_d_coupling_element_index=0;
    int sysxindex , xindex;
    for(i=0;i<total_vmode_list_size;i++){
        vector<int> vmode1_recv;
        vector<int> vmode2_recv;

        for(j=0;j<d.nmodes[0];j++){
            vmode1_recv.push_back(vmode_list_to_recv[vmode_element_index]);
            vmode_element_index++;
        }

        for(j=0;j<d.nmodes[1];j++){
            vmode2_recv.push_back(vmode_list_to_recv[vmode_element_index]);
            vmode_element_index++;
        }

        sysxindex=vmode_list_to_recv[vmode_element_index];
        vmode_element_index++;

        xindex= vmode_list_to_recv[vmode_element_index];
        vmode_element_index++;

        value=d_d_coupling_strength_list_recv[d_d_coupling_element_index];
        d_d_coupling_element_index++;

        n=find_location_binarysearch_sys_quotient_state(slist,vmode1_recv,vmode2_recv,exist3);
        if(exist3) {
            // we have found space that may have detector detector coupling with our slist[i].
            sys_xindex_size2=slist[n].sysxindex.size();
            for(b=0;b<sys_xindex_size2;b++){
                if(slist[n].sysxindex[b] == sysxindex){
                    d_d_irow.push_back(xindex);
                    d_d_icol.push_back(slist[n].xindex[b]);
                    d_d_mat.push_back(value);
                    detector_detector_coupling_term++;
                }
            }
        }
    }

    // rearrange off_diagpma_matrix
    rearrange_off_diagonal_term(d_d_mat,d_d_irow,d_d_icol);

    MPI_Reduce(&detector_detector_coupling_term,&total_detector_detector_coupling_term,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(my_id==0){
        cout<<"Total number of detector coupling in simulation:  "<<total_detector_detector_coupling_term<<endl;
        output<< "Total number of detector coupling in simulation:  "<<total_detector_detector_coupling_term<<endl;
    }
    //--------------- free the space ------------------------------
    delete [] vmode_list_1d;
    delete [] d_d_coupling_strength_to_send;
    delete [] vmode_list_size_each_process;
    delete [] vmode_list_displacement_each_process;
    delete [] vmode_list_size_1d_each_process;
    delete [] vmode_list_displacement_1d_each_process;
    delete [] vmode_list_to_recv;
    delete [] d_d_coupling_strength_list_recv;
}


void full_system:: rearrange_off_diagonal_term(vector < double > & sys_detector_mat,vector  <int> & sys_detector_irow,
        vector<int> & sys_detector_icol){
    /*
     rearrange term for off-diagonal part of mat, irow, icol. to make sure for each element, its irow is in
     corresponding process. Although called sys_detector_mat, it can be used for all off-diagonal part for full_matrix.
    */

    int i,j;
    int pc_index;
    int irow_index;
    int size;
    int index;
    int size_per_pc;
    size = sys_detector_irow.size();
    int * sys_detector_coup_to_send_pc_number = new int  [num_proc];  // number of element to send each pc
    int * sys_detector_coup_to_send_pc_displ= new int  [num_proc];  // displs for each pc.

    vector <int> * sys_detector_irow_to_send_each_pc = new vector <int> [num_proc];  // record irow to send for each pc
    vector <int> * sys_detector_icol_to_send_each_pc = new vector <int> [num_proc];  // record icol to send for each pc
    vector <double> * sys_detector_mat_to_send_each_pc = new vector <double> [num_proc]; // record mat to send for each pc

    int * sys_detector_irow_to_send = new int  [size];  //irow ordered by pc_index
    int * sys_detector_icol_to_send= new int  [size]; // icol ordered by pc_index
    double * sys_detector_mat_to_send= new double  [size];

    for(i=0;i<num_proc;i++){
        sys_detector_coup_to_send_pc_number[i] = 0;
    }
    for(i=0;i<size;i++){
        irow_index=sys_detector_irow[i];
        for(j=1;j<num_proc;j++){
            if(irow_index < matsize_offset_each_process[j]){
                break;
            }
        }
        pc_index= j-1;
        sys_detector_coup_to_send_pc_number[pc_index]++;
        sys_detector_irow_to_send_each_pc[pc_index].push_back(irow_index);
        sys_detector_icol_to_send_each_pc[pc_index].push_back(sys_detector_icol[i]);
        sys_detector_mat_to_send_each_pc[pc_index].push_back(sys_detector_mat[i]);
    }

    sys_detector_coup_to_send_pc_displ[0]=0;
    for(i=1;i<num_proc;i++){
        sys_detector_coup_to_send_pc_displ[i] = sys_detector_coup_to_send_pc_displ[i-1] +
                                                sys_detector_coup_to_send_pc_number[i-1];
    }

    index=0;
    for(i=0;i<num_proc;i++){
        size_per_pc = sys_detector_irow_to_send_each_pc[i].size();
        for(j=0;j<size_per_pc;j++){
            sys_detector_irow_to_send[index] = sys_detector_irow_to_send_each_pc[i][j] ;
            sys_detector_icol_to_send[index] = sys_detector_icol_to_send_each_pc[i][j];
            sys_detector_mat_to_send[index] = sys_detector_mat_to_send_each_pc[i][j];
            index++;
        }
    }

    int * sys_detector_coup_to_recv_num ;
    int * sys_detector_coup_to_recv_displ;
    int * sys_detector_coup_irow_to_recv;
    int * sys_detector_coup_icol_to_recv;
    double * sys_detector_coup_mat_to_recv;
    int  total_recv_num ;


    sys_detector_coup_to_recv_num = new int [num_proc];
    for(i=0;i<num_proc;i++){
        sys_detector_coup_to_recv_num[i]=0;
    }
    sys_detector_coup_to_recv_displ= new int [num_proc];
    MPI_Alltoall(sys_detector_coup_to_send_pc_number, 1,MPI_INT,
                 sys_detector_coup_to_recv_num, 1 , MPI_INT , MPI_COMM_WORLD);
    total_recv_num=0;
    for(i=0;i<num_proc;i++){
        total_recv_num = total_recv_num + sys_detector_coup_to_recv_num[i];
    }
    sys_detector_coup_to_recv_displ[0]=0;
    for(i=1;i<num_proc;i++){
        sys_detector_coup_to_recv_displ[i]= sys_detector_coup_to_recv_displ[i-1] + sys_detector_coup_to_recv_num[i-1];
    }
    sys_detector_coup_irow_to_recv = new int [total_recv_num];
    sys_detector_coup_icol_to_recv = new int [total_recv_num];
    sys_detector_coup_mat_to_recv = new double [total_recv_num];

    //  send irow, icol, mat
    MPI_Alltoallv(&sys_detector_irow_to_send[0],sys_detector_coup_to_send_pc_number,
                  sys_detector_coup_to_send_pc_displ,MPI_INT,&sys_detector_coup_irow_to_recv[0],
                  sys_detector_coup_to_recv_num, sys_detector_coup_to_recv_displ,MPI_INT,MPI_COMM_WORLD);

    MPI_Alltoallv(&sys_detector_icol_to_send[0],sys_detector_coup_to_send_pc_number,
                  sys_detector_coup_to_send_pc_displ,MPI_INT,&sys_detector_coup_icol_to_recv[0],
                  sys_detector_coup_to_recv_num,sys_detector_coup_to_recv_displ,MPI_INT,MPI_COMM_WORLD);

    MPI_Alltoallv(&sys_detector_mat_to_send[0],sys_detector_coup_to_send_pc_number,
                  sys_detector_coup_to_send_pc_displ, MPI_DOUBLE, &sys_detector_coup_mat_to_recv[0],
                  sys_detector_coup_to_recv_num,sys_detector_coup_to_recv_displ,MPI_DOUBLE,MPI_COMM_WORLD);

    // assign new value to mat, irow, icol
    sys_detector_mat.clear();
    sys_detector_mat.reserve(total_recv_num);
    sys_detector_irow.clear();
    sys_detector_irow.reserve(total_recv_num);
    sys_detector_icol.clear();
    sys_detector_icol.reserve(total_recv_num);

    for(i=0;i<total_recv_num;i++){
        sys_detector_mat.push_back(sys_detector_coup_mat_to_recv[i]);
        sys_detector_irow.push_back(sys_detector_coup_irow_to_recv[i]);
        sys_detector_icol.push_back(sys_detector_coup_icol_to_recv[i]);
    }


    // free the space
    delete [] sys_detector_coup_to_send_pc_number;
    delete [] sys_detector_coup_to_send_pc_displ;

    delete [] sys_detector_irow_to_send_each_pc;
    delete [] sys_detector_icol_to_send_each_pc;
    delete [] sys_detector_mat_to_send_each_pc;

    delete [] sys_detector_irow_to_send;
    delete [] sys_detector_icol_to_send;
    delete [] sys_detector_mat_to_send;


    delete [] sys_detector_coup_to_recv_num;
    delete [] sys_detector_coup_to_recv_displ;

    delete [] sys_detector_coup_irow_to_recv;
    delete [] sys_detector_coup_icol_to_recv;
    delete [] sys_detector_coup_mat_to_recv;

}