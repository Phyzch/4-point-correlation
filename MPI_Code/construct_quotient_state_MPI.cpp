//
// Created by phyzch on 6/29/20.
//
#include"../system.h"
#include"../util.h"
using namespace std;
vector<quotient_state> merge_d1list(vector<quotient_state> & qlist1, vector<quotient_state> & qlist2);
void full_system::construct_quotient_state_all_MPI(){
    int i;
    vector<sys_quotient_state> unsorted_slist;
    for(i=0;i<matsize;i++) {
        insert_quotient_state(d1list, sstate[i], d.dv_all[1][dstate[1][i]],  irow[i], dstate[0][i]);
        // d2list is already sorted because each process are gauranteed to be sorted.
        insert_quotient_state(d2list, sstate[i], d.dv_all[0][dstate[0][i]],  irow[i], dstate[1][i]);
    }
    //------------------------------------------------------------------
    // For d1list, it is marked by (detector2 state, system state). because each process only have a fraction of detector 1 state but all detector 2 state, the d1list constructed is incomplete and only contain part of detector 1 in dxindex.
    rearrange_d1list();

    for(i=0;i<matsize;i++){
        //---------------------------------------------------------------------------
        // slist: quotient system state list for photon system.
        // we will use slist to compute energy of photon and construct off diagonal matrix.
        //insert_sys_quotient_state(slist,vmode0[dstate[0][i]],vmode1[dstate[1][i]],d.moddim,irow[i],sstate[i],dstate[0][i], dstate[1][i]);
        sys_quotient_state sq(d.dv_all[0][dstate[0][i]],d.dv_all[1][dstate[1][i]],dstate[0][i],dstate[1][i]);
        sq.xindex.push_back(irow[i]);
        sq.sysxindex.push_back(sstate[i]);
        unsorted_slist.push_back(sq);
    }
    // now sort slist using merge_sort
    slist= merge_sort_list(unsorted_slist);
    unsorted_slist.clear();
    unsorted_slist.shrink_to_fit();

    // construct MPI version of q_index is easy to do. just let every process search the index in their local dlist.
    // q_index is detector Hamiltonian's element relation to location in full matrix.
    construct_q_index_MPI();
}

vector<vector<quotient_state>>  full_system::Gather_quotient_state_vec(){
    // return d1list collected from each process: a list of d1list: d1list each process.
    int i,j,k;
    // for d1list, we still need to re-sort d1list to make it sorted among all process to enable binary search.
    // Gather d1list to process 0 for merge: Process 0 will use O(matsize * log(matsize)) to do sorting.
    //------------  resize d1list in process 0 --------------------------------
    int d1list_size = d1list.size();
    int total_d1list_size; // total number of d1list in each process.
    int * total_d1list_vmode;
    int * total_d1list_sys_state;
    MPI_Reduce(&d1list_size,&total_d1list_size,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    vector< vector <quotient_state>> d1list_each_process; // contain d1list gathered from each process.
    if(my_id==0){
        d1list_each_process.reserve(num_proc);
        total_d1list_vmode = new int [total_d1list_size * d.nmodes[1]];
        total_d1list_sys_state = new int [total_d1list_size];
    }
    //------------ prepare size and displacement to send vmode and sys_state. Need in MPI_Gather. --------------------------------
    int * d1list_size_each_process;
    int * d1list_displacement ;
    int * d1list_vmode_size_each_process;
    int * d1list_vmode_displacement;

    d1list_size_each_process= new int [num_proc];
    d1list_displacement= new int [num_proc];
    d1list_vmode_size_each_process= new int [num_proc]; // number of element to receive from each process for vmode.
    d1list_vmode_displacement = new int [num_proc];

    MPI_Gather(&d1list_size,1,MPI_INT,&d1list_size_each_process[0],1,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id==0) {
        d1list_displacement[0] = 0;
        for (i = 1; i < num_proc; i++) {
            d1list_displacement[i] = d1list_displacement[i - 1] + d1list_size_each_process[i - 1];
        }
        for (i = 0; i < num_proc; i++) {
            d1list_vmode_size_each_process[i] = d1list_size_each_process[i] * d.nmodes[1];
        }
        d1list_vmode_displacement[0] = 0;
        for (i = 1; i < num_proc; i++) {
            d1list_vmode_displacement[i] = d1list_vmode_displacement[i - 1] + d1list_vmode_size_each_process[i - 1];
        }
    }

    //----------------------------Convert vmode and sys_state from d1list to vector.------------------------------------------------------
    int * sys_state_list = new int [d1list_size];     // 1d array sys_state_list contain d1list's system state
    for(i=0;i<d1list_size;i++){
        sys_state_list[i]= d1list[i].sys_state;
    }
    int * vmode_list= new int [d1list_size * d.nmodes[1]]; // 1d array vmode_list contain d1list's vmode
    int index=0;
    for(i=0;i<d1list_size;i++){
        for(j=0;j<d.nmodes[1];j++){
            vmode_list[index] = d1list[i].vmode[j];
            index++;
        }
    }
    //-------------------------------------Gather sys_state and vmode to process 0 ---------------------------------
    MPI_Gatherv(&sys_state_list[0],d1list_size,MPI_INT,
            &total_d1list_sys_state[0],d1list_size_each_process,d1list_displacement,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gatherv(&vmode_list[0],d1list_size*d.nmodes[1],MPI_INT,
            &total_d1list_vmode[0],d1list_vmode_size_each_process,d1list_vmode_displacement,MPI_INT,0,MPI_COMM_WORLD);

    //--------------------------------- send xlist, dxlist --------------------------
    int * xindex_size_per_state_single_process = new int [d1list_size];  // record xindex_size in each quotient state.
    int total_xindex_size_single_process=0;  // total_xindex_size in process.
    for(i=0;i<d1list_size;i++){
        xindex_size_per_state_single_process[i] = d1list[i].xindex.size();
        total_xindex_size_single_process = total_xindex_size_single_process + xindex_size_per_state_single_process[i];
    }
    int * xindex_all_single_process = new int [total_xindex_size_single_process];  // all xindex in this process stored in xindex_all_single_process
    int * dxindex_all_single_process = new int [total_xindex_size_single_process];// all dxindex in this process stored in dxindex_All_single_process
    index=0;
    for(i=0;i<d1list_size;i++){
        for(j=0;j<xindex_size_per_state_single_process[i];j++){
            xindex_all_single_process[index] = d1list[i].xindex[j];
            dxindex_all_single_process[index] = d1list[i].dxindex[j];
            index++;
        }
    }

    // ----------------------process 0 record each quotient state 's size of xinsdex, dxindex list----------
    int * xindex_size_each_quotient_state; // xindex size for each quotient state in each process.
    if(my_id ==0) {
        xindex_size_each_quotient_state = new int[total_d1list_size];  // total quotient state size.
    }
    MPI_Gatherv(&xindex_size_per_state_single_process[0],d1list_size,MPI_INT,
            &xindex_size_each_quotient_state[0],d1list_size_each_process,d1list_displacement,MPI_INT,0,MPI_COMM_WORLD);
//--------------------------  record xindex size each process  ---------------------------------------------------------------------------------------
    // to record xindex size for each process and use to gather all xindex list to process 0.
    int * total_xindex_size_each_process; // total_xindex_size in each process (gather in process id 0)
    int * total_xindex_displacement_each_process; // displacement for total_xindex_size
    total_xindex_size_each_process = new int [num_proc];
    total_xindex_displacement_each_process = new int [num_proc];

    MPI_Gather(&total_xindex_size_single_process,1,MPI_INT,&total_xindex_size_each_process[0],1,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id==0){
        total_xindex_displacement_each_process[0]=0;
        for(i=1;i<num_proc;i++){
            total_xindex_displacement_each_process[i]= total_xindex_displacement_each_process[i-1] + total_xindex_size_each_process[i-1];
        }
    }
//------------------------------------------------------------------------------------------------
    int total_xindex_list_size; // Sum of total_xindex_size_each_process
    int * xindex_all_process;  // store all xindex from different process.
    int * dxindex_all_process; // store all dxindex from different process
    if(my_id==0){
        total_xindex_list_size=0;
        for(i=0;i<num_proc;i++) {
            total_xindex_list_size = total_xindex_list_size +  total_xindex_size_each_process[i];
        }
        xindex_all_process = new int [total_xindex_list_size];
        dxindex_all_process = new int [total_xindex_list_size];
    }
    MPI_Gatherv(&xindex_all_single_process[0],total_xindex_size_single_process,MPI_INT,
            &xindex_all_process[0],total_xindex_size_each_process,total_xindex_displacement_each_process,MPI_INT,
            0,MPI_COMM_WORLD);
    MPI_Gatherv(&dxindex_all_single_process[0],total_xindex_size_single_process,MPI_INT,
            &dxindex_all_process[0],total_xindex_size_each_process,total_xindex_displacement_each_process,MPI_INT,
            0,MPI_COMM_WORLD);
//-------------------------------------------------------------------------------------------------------------------
    vector<vector<int>> xindex_list_all_vec;  // all xindex ready to reconstruct quotient state
    vector<vector<int>> dxindex_list_all_vec; // all dxindex ready to reconstruct quotient state
    if (my_id==0) {
        index = 0;
        for (i = 0; i < total_d1list_size; i++) {
            vector<int> xindex_vec;
            vector<int> dxindex_vec;
            for (j = 0; j < xindex_size_each_quotient_state[i]; j++) {
                xindex_vec.push_back(xindex_all_process[index]);
                dxindex_vec.push_back(dxindex_all_process[index]);
                index++;
            }
            xindex_list_all_vec.push_back(xindex_vec);
            dxindex_list_all_vec.push_back(dxindex_vec);
        }
    }

    // -------------convert total_d1list_vmode and total_d1list_sys_state and xindex, dxindex to d1list there ----------
    if(my_id==0){
        d1list_each_process.push_back(d1list);
        int sys_state_index= d1list_size ;  // we will not reconstruct d1list originally in process 0
        int vmode_index= d1list_size * d.nmodes[1]; // we will not reconstruct d1list originally in process 0
        int sys_state;
        for(k=1;k<num_proc;k++) {
            vector<quotient_state> temp_d1list;  // reconstruct d1list from other process.
            for (i = d1list_displacement[k];i < d1list_displacement[k] + d1list_size_each_process[k]; i++) {
                // we will not reconstruct d1list originally in process 0
                vector<int> vmode_temp;  // vmode characterize d1list
                for (j = 0; j < d.nmodes[1]; j++) {
                    vmode_temp.push_back(total_d1list_vmode[vmode_index]);
                    vmode_index++;
                }
                sys_state = total_d1list_sys_state[sys_state_index];  // system state characterize d1list
                sys_state_index++;
                quotient_state qs(vmode_temp, sys_state);  // quotient state qs constructed in d1list
                qs.xindex = xindex_list_all_vec[i];
                qs.dxindex = dxindex_list_all_vec[i];
                temp_d1list.push_back(qs);
            }
            d1list_each_process.push_back(temp_d1list);
        }
    }

    //---------------------free space -----------------------------------
    if (my_id ==0){
        delete [] total_d1list_vmode;
        delete [] total_d1list_sys_state;
        delete [] xindex_all_process;
        delete [] dxindex_all_process;
        delete [] xindex_size_each_quotient_state;
    }
    delete [] d1list_size_each_process;
    delete [] d1list_displacement;
    delete [] d1list_vmode_size_each_process;
    delete [] d1list_vmode_displacement;

    delete [] total_xindex_size_each_process;
    delete [] total_xindex_displacement_each_process;

    delete [] sys_state_list;
    delete [] vmode_list;
    delete [] xindex_size_per_state_single_process;
    delete [] xindex_all_single_process;
    delete [] dxindex_all_single_process;
    //-----------------------------------------------------------

    return (d1list_each_process);
}

vector<quotient_state> full_system::sort_d1list(vector<vector<quotient_state>> & d1list_each_process){
    // use d1list_each_process collected from  full_system::Gather_quotient_state_vec().
    int i;
    // now you sort d1list using d1list_each_process (each member contain part of d1list):
    // we use function merge_d1list to merge quotient list in process 0;
    int merge_list_size;  // size of d1list after merge.
    vector < quotient_state> sorted_d1list;
    if(my_id==0) {
        merge_list_size = num_proc;
        vector<vector<quotient_state>> * new_list_ptr;
        vector<vector<quotient_state>> * old_list_ptr;
        vector<vector<quotient_state>> *list_ptr_3;
        old_list_ptr = & (d1list_each_process);
        vector<vector<quotient_state>> new_list;
        new_list_ptr= & (new_list);
        while (merge_list_size != 1) {
            vector <quotient_state> merge_result;
            for (i = 0; i+1 < merge_list_size; i = i + 2) {
                merge_result = merge_d1list( (*old_list_ptr)[i], (*old_list_ptr)[i + 1]);
                (*new_list_ptr).push_back(merge_result);
            }
            if (merge_list_size % 2 == 1) {
                (* new_list_ptr).push_back( (*old_list_ptr) [merge_list_size - 1]);
            }

            (*old_list_ptr).clear();
            (*old_list_ptr).shrink_to_fit();  // free the space.
            // exchange old_list and new list
            list_ptr_3 = & (*old_list_ptr);
            old_list_ptr = &(* new_list_ptr);
            new_list_ptr = & (*list_ptr_3);

            merge_list_size= (merge_list_size+1)/2;  // odd number of merge_list will left the final one.
        }
        sorted_d1list =(*old_list_ptr)[0];
    }
    // Now we use merge sort to sort d1list. next we have to scatter them to all other process.
    return sorted_d1list;
}

quotient_state merge_quotient_state(quotient_state s1, quotient_state s2){
    int i;
    if( s1.sys_state != s2.sys_state or s1.vmode != s2.vmode){
        cout<<"Error, the two quotient_state do not equal to each other!"<<endl;
        MPI_Abort(MPI_COMM_WORLD,-14);
    }
    quotient_state s3(s1.vmode,s1.sys_state);
    int xindex_size1= s1.xindex.size();
    int xindex_size2= s2.xindex.size();
    int dx_position;
    for(i=0;i<xindex_size1;i++){
        dx_position=binary_insert_dxindex(s3.dxindex,s1.dxindex[i]);
        s3.xindex.insert(s3.xindex.begin()+ dx_position, s1.xindex[i]);
    }
    for(i=0;i<xindex_size2;i++){
        dx_position = binary_insert_dxindex(s3.dxindex,s2.dxindex[i]);
        s3.xindex.insert(s3.xindex.begin() + dx_position, s2.xindex[i]);
    }
    return s3;
}

vector<quotient_state> merge_d1list(vector<quotient_state> & qlist1, vector<quotient_state> & qlist2){
    // Algorithm used here: merge sort.  merge d1list in different process and return a new d1list.
    // if we found two quotient state in d1list is the same, we combine them together using merge_quotient_state() function.
    // input: qlist1, qlist2:  two quotient state list (d1list) in different process to merge
    // output: qlist3: quotient state list after merge.
    int mark;
    int qlist1_size = qlist1.size();
    int qlist2_size= qlist2.size();
    int q1_index=0;
    int q2_index=0;
    if(qlist1_size==0){
        return qlist2;
    }
    if(qlist2_size==0){
        return qlist1;
    }
    vector<quotient_state> qlist3;
    quotient_state s3 = qlist1[0];  // quotient state used to push into qlist3.
    while(q1_index<qlist1_size and q2_index < qlist2_size){
        mark = compare_quotient_state(qlist1[q1_index].sys_state,qlist1[q1_index].vmode,qlist2[q2_index].sys_state,qlist2[q2_index].vmode);
        if(mark>0){
            // qlist2[q2_index] is smaller.
            qlist3.push_back(qlist2[q2_index]);
            q2_index++;
        }
        else if (mark<0){
            //qlist1 [q1_index] is smaller
            qlist3.push_back(qlist1[q1_index]);
            q1_index++;
        }
        else{
            // two quotient_state is the same. merge them.
            s3= merge_quotient_state(qlist1[q1_index],qlist2[q2_index]);
            qlist3.push_back( s3 );
            q1_index++;
            q2_index++;
        }
    }
    if(q1_index== qlist1_size){
        // we push back the qlist2 left over
        while(q2_index<qlist2_size){
            qlist3.push_back(qlist2[q2_index]);
            q2_index++;
        }
    }
    if(q2_index== qlist2_size){
        // we push back the qlist1 left over
        while(q1_index<qlist1_size){
            qlist3.push_back(qlist1[q1_index]);
            q1_index++;
        }
    }
    return qlist3;
}

void full_system:: Scatter_sorted_d1list(vector<quotient_state> & sorted_d1list){
    // Scatter sorted_d1list to each process.
    int i,j;
    int * vsize_each_process;  // size of d1list scatter to each process.
    int * vsize_displacement_each_process; // displacement of d1list to scatter in each process.
    int * total_sys_state_list;

    int *vsize_vmode_each_process; // size of vmode to scatter to each process.
    int *vsize_vmode_displacement_each_process; // displacement for vmode to scatter to each process.
    int * total_vmode_list; // 1d array that record all vmode

    int *vsize_xindex_size_each_state;  // xindex size each quotient state
    int * vsize_xindex_size_each_process;  // xindex size in each process.
    int * vsize_xindex_displacement_each_process;

    int * total_xindex_list;  // list record all xindex
    int * total_dxindex_list;  // list record all dxindex

    vsize_each_process = new int[num_proc];  // quotient state each process
    vsize_displacement_each_process = new int[num_proc];  // displacement for vsize_each_process

    vsize_vmode_each_process = new int[num_proc];  // quotient state's mode for each process.
    vsize_vmode_displacement_each_process = new int[num_proc];  // quotient state's mode's displacement for each process.

    vsize_xindex_size_each_process = new int [num_proc];  // size of 1d array for xindex.
    vsize_xindex_displacement_each_process = new int [num_proc];  // displacement for 1d array for xindex

    if(my_id==0) {
        int index;
        int d1list_total_size = sorted_d1list.size();
        int vsize = d1list_total_size / num_proc;
        int vsize2 = d1list_total_size - (num_proc - 1) * vsize;

        // for sys_state
        for (i = 0; i < num_proc - 1; i++) {
            vsize_each_process[i] = vsize;
        }
        vsize_each_process[num_proc - 1] = vsize2;

        vsize_displacement_each_process[0] = 0;
        for (i = 1; i < num_proc; i++) {
            vsize_displacement_each_process[i] = vsize_displacement_each_process[i - 1] + vsize_each_process[i - 1];
        }

        total_sys_state_list = new int [d1list_total_size];
        for(i=0;i<d1list_total_size;i++){
            total_sys_state_list[i] = sorted_d1list[i].sys_state;
        }

        // for vmode
        for (i = 0; i < num_proc; i++) {
            vsize_vmode_each_process[i] = vsize_each_process[i] * d.nmodes[1];
        }

        vsize_vmode_displacement_each_process[0] = 0;
        for (i = 1; i < num_proc; i++) {
            vsize_vmode_displacement_each_process[i] =
                    vsize_vmode_displacement_each_process[i - 1] + vsize_vmode_each_process[i - 1];
        }
        total_vmode_list = new int [d1list_total_size*d.nmodes[1]];
        index=0;
        for(i=0;i<d1list_total_size;i++){
            for(j=0;j<d.nmodes[1];j++){
                total_vmode_list[index] = sorted_d1list[i].vmode[j];
                index++;
            }
        }

        // for xindex, dxindex
        vsize_xindex_size_each_state = new int[d1list_total_size];    // xindex size for each quotient state
        int total_vsize_xindex_size = 0;       // total xindex size for all quotient state
        for (i = 0; i < d1list_total_size; i++) {
            vsize_xindex_size_each_state[i] = sorted_d1list[i].xindex.size();
            total_vsize_xindex_size = total_vsize_xindex_size + vsize_xindex_size_each_state[i];
        }


        for(i=0;i<num_proc;i++){
            vsize_xindex_size_each_process[i]=0;
            for(j=vsize_displacement_each_process[i];j< vsize_displacement_each_process[i] + vsize_each_process[i];j++){
                vsize_xindex_size_each_process[i] = vsize_xindex_size_each_process[i] + vsize_xindex_size_each_state[j];
            }
        }

        vsize_xindex_displacement_each_process[0]=0;
        for(i=1;i<num_proc;i++){
            vsize_xindex_displacement_each_process[i] =
                    vsize_xindex_displacement_each_process[i-1] + vsize_xindex_size_each_process[i-1];
        }

        total_xindex_list = new int [total_vsize_xindex_size];  // array prepare to send xindex
        total_dxindex_list = new int [total_vsize_xindex_size]; // array prepare to send dxindex
        index =0;
        int xindex_size;
        for(i=0;i<d1list_total_size;i++){
            xindex_size= sorted_d1list[i].xindex.size();
            for(j=0;j<xindex_size;j++){
                total_xindex_list[index] = sorted_d1list[i].xindex[j];
                total_dxindex_list[index] = sorted_d1list[i].dxindex[j];
                index++;
            }
        }
    }
    // scatter sys_state
    int quotient_state_number;  // number of quotient contain in d1list.
    MPI_Scatter(&vsize_each_process[0],1,MPI_INT,&quotient_state_number,1,MPI_INT,0,MPI_COMM_WORLD);
    d1list.clear();
    d1list.reserve(quotient_state_number);
    int * quotient_state_sys_1d = new int [quotient_state_number];
    MPI_Scatterv(&total_sys_state_list[0],vsize_each_process,vsize_displacement_each_process,MPI_INT,
            &quotient_state_sys_1d[0],quotient_state_number,MPI_INT,0,MPI_COMM_WORLD);

    // scatter vmode
    int quotient_state_vmode_number;
    quotient_state_vmode_number = quotient_state_number  * d.nmodes[1];
    int * quotient_state_vmode_1d = new int [quotient_state_vmode_number];
    MPI_Scatterv(&total_vmode_list[0],vsize_vmode_each_process,vsize_vmode_displacement_each_process,MPI_INT,
            &quotient_state_vmode_1d[0],quotient_state_vmode_number,MPI_INT,0,MPI_COMM_WORLD);


    //scatter xindex, dxindex
    int total_xindex_size_single_process;  // total number of xindex
    MPI_Scatter(&vsize_xindex_size_each_process[0],1,MPI_INT,
            &total_xindex_size_single_process,1,MPI_INT,0,MPI_COMM_WORLD);
    int * xindex_size_each_quotient_state = new int [quotient_state_number];  // xindex_size for each quotient state in list
    MPI_Scatterv(&vsize_xindex_size_each_state[0],vsize_each_process,vsize_displacement_each_process,MPI_INT,
            &xindex_size_each_quotient_state[0],quotient_state_number,MPI_INT,0,MPI_COMM_WORLD);

    int * xindex_each_quotient_state_1d = new int [total_xindex_size_single_process];
    MPI_Scatterv(&total_xindex_list[0],vsize_xindex_size_each_process,vsize_xindex_displacement_each_process,MPI_INT,
            &xindex_each_quotient_state_1d[0],total_xindex_size_single_process,MPI_INT,0,MPI_COMM_WORLD);
    int * dxindex_each_quotient_state_1d = new int [total_xindex_size_single_process];
    MPI_Scatterv(&total_dxindex_list[0],vsize_xindex_size_each_process,vsize_xindex_displacement_each_process,MPI_INT,
            &dxindex_each_quotient_state_1d[0],total_xindex_size_single_process,MPI_INT,0,MPI_COMM_WORLD);

    // ----------------------------convert it to quotient state list ------------------------------------
    int sys_state;
    int vmode_index=0;
    int xindex_index=0;
    for(i=0;i<quotient_state_number;i++){
        sys_state = quotient_state_sys_1d[i];
        vector<int> vmode;
        for(j=0;j<d.nmodes[1];j++){
            vmode.push_back(quotient_state_vmode_1d[vmode_index]);
            vmode_index++;
        }
        quotient_state sq (vmode,sys_state);
        // insert xindex, dxindex
        for(j=0;j<xindex_size_each_quotient_state[i];j++){
            sq.xindex.push_back(xindex_each_quotient_state_1d[xindex_index]);
            sq.dxindex.push_back(dxindex_each_quotient_state_1d[xindex_index]);
            xindex_index++;
        }
        d1list.push_back(sq);
    }

    //--------------------free the space --------------------------------
    if(my_id==0){
        delete [] total_sys_state_list;

        delete [] total_vmode_list;

        delete [] vsize_xindex_size_each_state;

        delete [] total_xindex_list;
        delete [] total_dxindex_list;
    }
    delete [] vsize_each_process;
    delete [] vsize_displacement_each_process;

    delete [] vsize_vmode_each_process;
    delete [] vsize_vmode_displacement_each_process;

    delete [] vsize_xindex_size_each_process;
    delete []  vsize_xindex_displacement_each_process;

    delete []  quotient_state_sys_1d;
    delete []  quotient_state_vmode_1d;
    delete [] xindex_size_each_quotient_state;
    delete [] xindex_each_quotient_state_1d;
    delete [] dxindex_each_quotient_state_1d;
}

void full_system:: rearrange_d1list(){
    vector <vector<quotient_state>> d1list_each_process = Gather_quotient_state_vec();
    vector<quotient_state> sorted_d1list= sort_d1list(d1list_each_process);
    Scatter_sorted_d1list(sorted_d1list);
}

void full_system::construct_q_index_MPI(){
    // Before we compute Reduced density matrix of detector to compute detector Energy,
    // we have to compute q_index to speed up calculation
    // This q_index is also used to construct off diagonal matrix for full system.
    // q_index compose detector's Hamiltonian and its corresponding element in full matrix.
    int m,i,j,p,n;
    int k1,l1;
    int k,l;
    double value;
    bool exist1, exist2;
    vector <quotient_state> * dlist_ptr;
    int list_size;
    for(m=0;m<s.tlnum;m++){
        if(m==0) dlist_ptr = &(d1list);  // d1list: quotient state list for detector 1
        else dlist_ptr = &(d2list);  // d2list: quotient state list for detector 2
        list_size= (*dlist_ptr).size();
        for(p=0;p<d.total_dmat_num[m];p++){
            i= d.total_dirow[m][p];
            j= d.total_dicol[m][p];
            if(i>j) continue; // we only record half the result.
            value= d.total_dmat[m][p];
            for(n=0;n<list_size;n++){
                // go through quotient space for specific detector.
                k1=binary_search_dxindex((*dlist_ptr)[n].dxindex,i,exist1);
                if(i==j){
                    exist2=exist1;
                    l1=k1;
                }
                else{
                    l1= binary_search_dxindex((*dlist_ptr)[n].dxindex,j,exist2);
                }
                if(exist1 and exist2){
                    k= (*dlist_ptr)[n].xindex[k1];
                    l=(*dlist_ptr)[n].xindex[l1];
                    vector <int> qindex = {i,j,k,l,p};
                    (*dlist_ptr)[n].q_index_list.push_back(qindex);
                    (*dlist_ptr)[n].dmat_value_list.push_back(value);
                }
            }
        }
    }
}

