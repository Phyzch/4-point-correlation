//
// Created by phyzch on 7/13/20.
//
#include "../system.h"
#include"../util.h"
using namespace std;

void full_system::evaluate_system_output_MPI(double *hx, double * hy, double &se, double &s0, double &s1, double &s2,
        double &trsr2, double * de,  double ** mode_quanta, complex<double> ** sr,
        complex<double> ** dr, complex<double> ** total_dr){
    // evaluate the output of whole system
    etot_MPI(hx,hy);
    compute_sys_energy_MPI(se,s0,s1,s2);
    detenergy_MPI(de,dr,total_dr);
    if(my_id==0) {
        output << setprecision(3);
        output << setw(4) << t << "  " << setw(4) << s1 << " " << setw(4) << s2 << "  " << setw(4) << trsr2
               << "  " << setw(7) << std::setprecision(7) << se << "  " << setw(4) << std::setprecision(7) << de[0]
               << "  " << setw(4) << std::setprecision(7) << de[1] << "  " << setw(4) << total_energy / cf
               << "  " << setw(4) << total_norm << "  " << setw(4) << s0 << endl;
    }
    average_vibrational_mode_quanta_MPI(total_dr,mode_quanta);
}


void full_system:: etot_MPI(double * hx, double * hy){
    int i;
    int irow_index, icol_index;
    double local_test=0;
    double total_test=0;
    double local_e=0;
    double local_e_add = 0;
    norm=0;
    total_norm=0;
    // update tosend_xd, tosend_yd to send data to other process
    update_x_y();

    for(i=0;i<matsize;i++){
        hx[i]=0;
        hy[i]=0;
    }

    for(i=0;i<matnum;i++){
        irow_index = local_irow[i];
        icol_index = local_icol[i];
        hx[irow_index] = hx[irow_index] + mat[i] * x[icol_index];
        hy[irow_index] = hy[irow_index] + mat[i] * y[icol_index];
    }

    for(i=0;i<matsize;i++){
        local_e = local_e + x[i] * hx[i] + y[i] * hy[i];
        norm = norm + x[i] * x[i] + y[i] * y[i];
        local_test= local_test+ hx[i]*y[i]- hy[i]*x[i];
    }
    total_energy=0;
    MPI_Allreduce(&local_e,&total_energy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&local_test,&total_test,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(abs(total_test)>0.01){
        if(my_id==0){
            log << "Warning: The Hermiticity  is not conserved in Etot " << endl;
            cout << "Warning: The Hermiticity  is not conserved in Etot " << endl;
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
    }
    if(abs(total_norm -1)>0.01){
        cout<<"Wrong"<<endl;
        if(my_id == 0){
            log <<" Warning the norm of whole system is not conserved "<<endl;
            cout <<"Warning the norm of whole system is not conserved"<<endl;
            cout<<total_norm<<endl;
        }
    }
    total_energy= total_energy/ total_norm;
}


void full_system:: compute_sys_energy_MPI(double & se, double &s0, double &s1, double &s2){
    int i,j;
    int size= slist.size();
    int sys_xindex; // index in system state: 0,1,2,3
    int xindex;  // index in x
    int local_xindex;
    int xindex_list_size;
    double sr_norm = 0;
    sys_quotient_state * sq;

    double * sr = new double [s.tlmatnum];
    double  * total_sr = new double [s.tlmatnum];
    for(i=0;i<s.tlmatnum;i++){
        sr[i]=0;
        total_sr[i]=0;
    }
    for(i=0;i<size;i++) {
        sq = &(slist[i]);
        xindex_list_size = (*sq).xindex.size();
        for (j = 0; j < xindex_list_size; j++) {
            xindex = (*sq).xindex[j];
            // all xindex here are local in process. This is because the way we construct slist is pushing
            // diagonal matrix term local in process to construct slist.
            local_xindex = xindex - matsize_offset_each_process[my_id];
            sys_xindex=(*sq).sysxindex[j];
            sr[sys_xindex] = sr[sys_xindex] + pow(x[local_xindex],2) + pow(y[local_xindex],2);
        }
    }
    MPI_Allreduce(&sr[0],&total_sr[0],s.tlmatnum,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for(i=0;i<s.tlmatnum;i++){
        sr_norm = sr_norm + total_sr[i];
    }

    for(i=0;i<s.tlmatnum;i++){
        total_sr[i] = total_sr[i] / sr_norm ;
    }
//    if(my_id == 0){
//        cout <<"sr_norm: "<<endl;
//        cout << sr_norm <<endl;
//    }
    if( abs( sr_norm -1)> 0.01 ){
        if(my_id==0){
            log << "Warning: Norm of system reduced density matrix is not conserved."<<endl;
            log<<sr_norm<<endl;
            cout << "Warning: Norm of system reduced density matrix is not conserved." << endl;
            cout << sr_norm <<endl;
        }
    }
    se=0;
    for(i=0;i<s.tlmatnum;i++){
        se = se + total_sr[i] * s.tlmat[i];
    }

    delete [] sr;
    delete [] total_sr;
}

void full_system::detenergy_MPI(double * de, complex <double> ** dr, complex <double> ** total_dr){
    /*  input: de: record energy of two detector.
 *         dr: density matrix of detector, size: [s.tlnum] [dmatnum] same size as dmat
 */
    int m,n;
    int i;
    int k,l,p;
    int index;
    int dr_index_list_size = dr_index_list.size();

    update_x_y_for_detenergy();

    for(m=0;m<s.tlnum;m++) {
        // initialize dr
        for (n = 0; n < d.total_dmat_num[m]; n++) {
            dr[m][n] = 0;
        }
    }
    index=0;
    for(i=0;i<dr_index_list_size;i++){
        k = local_vector_index_for_detenergy[index];
        index++;
        l = local_vector_index_for_detenergy[index];
        index++;
        p= dr_index_list[i][2];  // index in detector Hamiltonian
        m= dr_index_list[i][3];  // index for detector.
        dr[m][p] = dr[m][p] + complex<double> (x_for_detenergy[k] , -y_for_detenergy[k]) *
                              complex<double> ( x_for_detenergy[l] , y_for_detenergy[l] );
    }

    for(m=0;m<s.tlnum;m++){
        MPI_Allreduce(&dr[m][0], & total_dr[m][0], d.total_dmat_num[m],MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD);
    }
    
    //Normalize result.
    double total_de_norm = 0;
    for(m=0;m<s.tlnum;m++){
        total_de_norm = 0;

        for(i=0;i<num_proc;i++){
            for(n=d.dmat_offset_each_process[m][i];n<d.dmat_offset_each_process[m][i] + d.dmatsize_each_process[m][i]; n++){
                // diagonal term
                total_de_norm = total_de_norm + real(total_dr[m][n]);
            }
        }

        for(i=0;i<d.total_dmat_num[m];i++){
            total_dr[m][i] = total_dr[m][i] / total_de_norm ;
        }
    }
    // dr is not organlized, we want to multiply 2 to the off-diag terms. How to identify off-diag term?
    for(m=0;m<s.tlnum;m++){
        de[m]=0;
        for(i=0;i<num_proc;i++){
            for(n=d.dmat_offset_each_process[m][i];n<d.dmat_offset_each_process[m][i] + d.dmatsize_each_process[m][i]; n++){
                // diagonal part
                de[m] = de[m] + real(d.total_dmat[m][n] * total_dr[m][n]);
            }
            for(n=d.dmat_offset_each_process[m][i] + d.dmatsize_each_process[m][i]; n<d.dmat_offset_each_process[m][i] + d.dmatnum_each_process[m][i];n++){
                // off diagonal part
                de[m] = de[m] + real(d.total_dmat[m][n] * total_dr[m][n])*2;
            }
        }
    }
}

void full_system::prepare_detenergy_computation_MPI(){
    int m,n,r;
    int i,j;
    int size, qlist_size;
    vector<quotient_state> * dlist;
    vector <vector <int>> * qlist;
    for(m=0;m<s.tlnum;m++){
        if(m==0) dlist= &(d1list);
        else if (m==1) dlist= &(d2list);
        size= (*dlist).size();
        for(n=0;n<size;n++){
            qlist = &((*dlist)[n].q_index_list);
            qlist_size= (*qlist).size();
            for(r=0;r<qlist_size;r++){
                vector <int> dr_index;
                dr_index.reserve(4);
                dr_index.push_back((*qlist)[r][2]);  // index stored in xindex,
                // corresponding photon+detector state used for compute density matrix state |i>
                dr_index.push_back( (*qlist)[r][3] );  // index stored in xindex,
                // corresponding photon + detector state used for compute density matrix state <j|
                dr_index.push_back( (*qlist)[r][4] );  // p is corresponding index of |i><j| in dr (or dmat)
                dr_index.push_back( m );  // m is detector index
                dr_index_list.push_back(dr_index);
            }
        }
    }

    // Collect all xindex together for communication in dr_xindex_list
    int dr_index_list_size = dr_index_list.size();
    vector <int> dr_xindex_list;  // collect all xindex stored in dlist1, dlist2 's qlist.
    dr_xindex_list.reserve(2*dr_index_list_size);
    // will sort dr_xindex_list and use it to construct vector to communicate between process.
    for(i=0;i<dr_index_list_size;i++){
        dr_xindex_list.push_back(dr_index_list[i][0]);
        dr_xindex_list.push_back(dr_index_list[i][1]);    // record all xindex we need to recv vector from other process
    }
    int dr_xindex_list_size = dr_xindex_list.size();
    sort(dr_xindex_list.begin(),dr_xindex_list.end());


    remote_vec_count_for_detenergy =  new int [num_proc];
    remote_vec_ptr_for_detenergy = new int [num_proc];
    vector <int> * remote_vec_index_for_detenergy_each_process = new vector <int> [num_proc];

    for(i=0;i<num_proc;i++){
        remote_vec_count_for_detenergy[i] = 0;
    }

    int prev_xindex =-1;
    int xindex;
    int xindex_pc_id;  // process id this xindex belong to.

    for(i=0;i<dr_xindex_list_size;i++){
        xindex = dr_xindex_list[i];
        if(xindex > prev_xindex) {
            for (j = 1; j < num_proc; j++) {
                if (xindex < matsize_offset_each_process[j]) {
                    break;
                }
            }
            xindex_pc_id = j-1;
            remote_vec_count_for_detenergy[xindex_pc_id] ++;
            remote_vec_index_for_detenergy_each_process[xindex_pc_id].push_back(xindex);
        }
        prev_xindex = xindex;
    }

    total_remote_vec_num_for_detenergy = 0;
    for(i=0;i<num_proc;i++){
        total_remote_vec_num_for_detenergy = total_remote_vec_num_for_detenergy + remote_vec_count_for_detenergy[i];
    }

    remote_vec_ptr_for_detenergy[0]=0;
    for(i=1;i<num_proc;i++){
        remote_vec_ptr_for_detenergy[i] = remote_vec_ptr_for_detenergy [i-1] + remote_vec_count_for_detenergy[i-1];
    }

    remote_vec_index_for_detenergy = new int [total_remote_vec_num_for_detenergy];
    int index = 0;
    for(i=0;i<num_proc;i++){
        for(j=0;j<remote_vec_count_for_detenergy[i];j++){
            remote_vec_index_for_detenergy[index] = remote_vec_index_for_detenergy_each_process[i][j];
            index++;
        }
    }

    int * search_Ind;

    local_vector_index_for_detenergy.reserve(dr_xindex_list_size);
    index=0;
    for(i=0;i<dr_index_list_size;i++){
        for(m=0;m<2;m++) {
            xindex = dr_index_list[i][m];
            search_Ind = (int *) bsearch(&xindex, remote_vec_index_for_detenergy, total_remote_vec_num_for_detenergy,
                                         sizeof(int), compar);
            local_vector_index_for_detenergy[index] = search_Ind - remote_vec_index_for_detenergy;
            index++;
        }
    }

    x_for_detenergy = new double [total_remote_vec_num_for_detenergy];
    y_for_detenergy = new double [total_remote_vec_num_for_detenergy];

    to_send_vec_count_for_detenergy = new int [num_proc];
    to_send_vec_ptr_for_detenergy = new int [num_proc];
    total_to_send_vec_num_for_detenergy =
            construct_send_buffer_index(remote_vec_count_for_detenergy,remote_vec_ptr_for_detenergy,
            remote_vec_index_for_detenergy,to_send_vec_count_for_detenergy,to_send_vec_ptr_for_detenergy,
            to_send_vec_index_for_detenergy);

    send_x_for_detenergy = new double [total_to_send_vec_num_for_detenergy];
    send_y_for_detenergy = new double [total_to_send_vec_num_for_detenergy];


    // MPI_ Gather dr.

}

void full_system:: update_x_y_for_detenergy(){
    int i;
    int local_index;
    for(i=0;i<total_to_send_vec_num_for_detenergy;i++){
        local_index = to_send_vec_index_for_detenergy[i] - matsize_offset_each_process[my_id];
        send_x_for_detenergy [i] = x[local_index];
        send_y_for_detenergy [i] = y[local_index];
    }

    MPI_Alltoallv(&send_x_for_detenergy[0],to_send_vec_count_for_detenergy,to_send_vec_ptr_for_detenergy,MPI_DOUBLE,
                  &x_for_detenergy[0],remote_vec_count_for_detenergy,remote_vec_ptr_for_detenergy,MPI_DOUBLE, MPI_COMM_WORLD);

    MPI_Alltoallv(&send_y_for_detenergy[0],to_send_vec_count_for_detenergy,to_send_vec_ptr_for_detenergy,MPI_DOUBLE,
                  &y_for_detenergy[0],remote_vec_count_for_detenergy,remote_vec_ptr_for_detenergy,MPI_DOUBLE,MPI_COMM_WORLD);

}

void full_system:: average_vibrational_mode_quanta_MPI(complex <double> ** total_dr, double ** mode_quanta){
    int m,i,j;
    vector<vector<int>> * vmode_ptr;
    if(my_id==0) {
        for (m = 0; m < s.tlnum; m++) {
            if (m == 0) {
                vmode_ptr = &d.dv_all[0];
            } else {
                vmode_ptr = &d.dv_all[1];
            }
            for (j = 0; j < d.nmodes[m]; j++) {
                mode_quanta[m][j] = 0;
                for (i = 0; i < d.total_dmat_size[m]; i++) {
                    mode_quanta[m][j] = mode_quanta[m][j] + (*vmode_ptr)[i][j] * real(total_dr[m][i]);
                }
            }
        }
        for(m=0;m<s.tlnum;m++){
            Detector_mode_quanta<<"t=  "<< t << "  Detector "<< m << "  vibrational mode quanta  "<<endl;
            for(j=0;j<d.nmodes[m];j++){
                Detector_mode_quanta << mode_quanta[m][j]<<"  ";
            }
            Detector_mode_quanta<<endl;
        }
    }
}