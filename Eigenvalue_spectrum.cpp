//
// Created by phyzch on 1/6/21.
//
#include "util.h"
#include "system.h"

double pythag(double a, double b){
    // compute a^2 + b^2 with high precision
    double absa = abs(a);
    double absb = abs(b);
    double pythag ;
    if( absa > absb){
        pythag = absa * sqrt(1+pow(absb/absa,2));
    }
    else{
        if(absb == 0){
            pythag = 0;
        }
        else{
            pythag = absb * sqrt(1+pow(absa/absb, 2) );
        }
    }
    return pythag;
}

void detector::prepare_computation_for_Lanczos(){
    // This function should be called before we call diagonalize()

    // compute buffer to receive and send for each process.
    // Index for remoteVecIndex, tosendVecIndex are computed here.
    int m,i;
    int vsize;
    int nearby_state_list_size = nearby_state_index.size();
// Index for vector to send and receive.
// remoteVecCount: total number to receive. remoteVecPtr: displacement in remoteVecIndex for each process. remoteVecIndex: index in other process to receive.
// tosendVecCount: total number to send to other process. tosendVecPtr: displacement in tosendVecIndex in each process.  tosendVecIndex: Index of element in itself to send. (it's global ,need to be converted to local index)
    remoteVecCount= new int * [1];
    remoteVecPtr= new int * [1];
    remoteVecIndex= new int * [1];
    to_recv_buffer_len = new int  [1];

//------------------Allocate space for vector to receive ---------------------
    for (m=0;m<1;m++){
        remoteVecCount[m] = new int [num_proc];
        remoteVecPtr[m] = new int [num_proc];
        remoteVecIndex[m] = new int [dmatnum[m]];
        for(i=0;i<num_proc;i++){
            remoteVecCount[m][i] = 0;
        }
    }
    tosendVecCount= new int *[1];
    tosendVecPtr =  new int * [1];
    tosendVecIndex = new int * [1];
    to_send_buffer_len= new int [1];
    for(m=0;m<1;m++){
        tosendVecCount[m] = new int [num_proc];
        tosendVecPtr[m] = new int [num_proc];
    }

    int * search_Ind; // local variable, used for compute local_dicol;
    int col_index_to_search;
// local column index used when we do H *x and H*y
    local_dirow= new vector<int> [1];
    local_dicol = new vector<int> [1]; // column index for computation in local matrix.

    vsize= total_dmat_size[0]/num_proc;
    to_recv_buffer_len[0] = construct_receive_buffer_index(remoteVecCount[0],remoteVecPtr[0],
                                                           remoteVecIndex[0],0);  // construct buffer to receive.
    to_send_buffer_len[0]= construct_send_buffer_index(remoteVecCount[0],remoteVecPtr[0],remoteVecIndex[0],
                                                       tosendVecCount[0],tosendVecPtr[0], tosendVecIndex[0]);

    // construct local_dirow, local_dicol
    local_dirow= new vector<int> [1];
    local_dicol = new vector<int> [1]; // column index for computation in local matrix.

    local_dirow[0].reserve(dmatnum[0]);
    local_dicol[0].reserve(dmatnum[0]);
    for(i=0;i<dmatnum[0];i++){
        local_dirow[0].push_back(dirow[0][i] - my_id * vsize);  // set local index for row index
        col_index_to_search= dicol[0][i];
        search_Ind=(int *) bsearch(&col_index_to_search,remoteVecIndex[0],to_recv_buffer_len[0],sizeof(int),compar);
        if(search_Ind!=NULL){
            // this column index is not in local matrix, and we should get it from other process (remoteVec)
            local_dicol[0].push_back(dmatsize[0] + (search_Ind-remoteVecIndex[0]) );
        }
        else{ // this column index is in local matrix.
            local_dicol[0].push_back (dicol[0][i] - my_id * vsize );
        }
    }
}

void detector:: allocate_space_for_vector_in_Lanczos(double * &local_v, double * &local_vold, double * &local_w, double * &send_v, double * &recv_v){
    local_v = new double [dmatsize[0] + to_recv_buffer_len[0]];
    local_vold = new double [dmatsize[0]];
    local_w = new double [dmatsize[0]];
    // prepare for w = H * v
    send_v = new double [to_send_buffer_len[0]];
    recv_v = new double [to_recv_buffer_len[0]];
}

void detector:: update_v_in_Lanczos(double * local_v, double * send_v, double * recv_v){
    int i;
    int vsize;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;
    for(i=0;i<to_send_buffer_len[0];i++){
        send_v[i] = local_v[tosendVecIndex[0][i] - my_id * vsize];
    }
    MPI_Alltoallv(&send_v[0],tosendVecCount[0],tosendVecPtr[0], MPI_DOUBLE,
                  &recv_v[0], remoteVecCount[0],remoteVecPtr[0],MPI_DOUBLE,MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len[0];i++){
        local_v[i + dmatsize[0]] = recv_v[i];
    }
}

void detector::compute_Lanczos( double * alpha, double * beta,
                               double * local_v, double * local_v_old, double * local_w, double * send_v, double * recv_v,
                               int nlev, int maxit){
    /*
     * maxit : dimension of matrix (alpha, beta) we finally compute in Lanczos algorithm
     * v is vector we use to compute Lanczos algorithm
     * vold is old vector
     * alpha[maxit], beta[maxit] : alpha is diagonal part in tridiagonal matrix, beta is subdiagonal part in tri-diagonal matrix
     * results store in alpha, beta
     */
    int i,j;
    int irow, icol;
    double local_alpha, local_beta;
    double ibet;  // inverse of beta

    beta[0] = 0;

    for(i=0;i<dmatsize[0];i++){
        local_v_old[i]  = 0;
    }

    for(i=0;i<maxit;i++){
        // maxit is number of eigenvalue we want to compute with Lanczos algorithm
        // results store in alpha[0~maxit-1] , beta[0~maxit-1]. beta[maxit] is not used. maxit <= lancdim- 1 to make sure we do not have overflow
        for(j=0;j<dmatsize[0];j++){
            local_w [j] = 0;
        }
        // matrix multiplication below: w1 = A * v1.  parallel version:
        update_v_in_Lanczos(local_v,send_v,recv_v);
        for(j=0;j<dmatnum[0];j++){
            irow = local_dirow[0][j];
            icol = local_dicol[0][j];
            local_w[irow] = local_w[irow] + real(dmat[0][j]) * local_v[icol];
        }

        // alpha = w^{T} * v
        local_alpha = 0;
        for(j=0;j<dmatsize[0];j++){
            local_alpha = local_alpha + local_w[j] * local_v[j];
        }
        MPI_Allreduce(&local_alpha,&alpha[i],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        // wi = wi' -alphai * vi - betai * v_{i-1}: here v_{i} = v, v_{i-1} = vold
        for(j=0;j<dmatsize[0];j++){
            local_w[j] = local_w[j] - alpha[i] * local_v[j] - beta[i] * local_v_old[j];
        }


        // beta[i+1] is Euclidean norm of w[j]
        local_beta = 0;
        for(j=0;j<dmatsize[0];j++){
            local_beta = local_beta + local_w[j] * local_w[j];
        }
        MPI_Allreduce(&local_beta,&beta[i+1],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        beta[i+1] = sqrt(beta[i+1]);

        // inverse of beta
        ibet = 1/beta[i+1];

        // vold = v.  v= w / beta[i+1]
        for(j=0;j<dmatsize[0];j++){
            local_v_old[j] = local_v[j];
            local_v[j] = local_w[j] * ibet;
        }

    }

}

double SIGN(double a,double b){
    if(b>=0) return abs(a);
    else return -abs(a);
}

void tqli( double * d, double *e , int n, ofstream & eigenvalue_log_output) {
    /*
    Implicit QL algorithm implemented here
     d is diagonal part of tri-diagonal matrix : 1~n
    e is subdiagonal part of tri-diagonal matrix : 1~n
    n: size of matrix

     for tqli: https://www.cimat.mx/~posada/OptDoglegGraph/DocLogisticDogleg/projects/adjustedrecipes/tqli.cpp.html
    */
    int i, l, m, iter;
    double g, r, s, c, p, f, b, dd;

    if (n > 1) {
        for (i = 2; i <= n; i++) {
            e[i - 1] = e[i];  // shift off-diagonal element
        }
        e[n] = 0;

        for (l = 1; l <= n; l++) {
            iter = 0;
            do {
                m = l ;
                for (m = l; m <= n - 1; m++) {
                    dd = fabs(d[m]) + fabs(d[m + 1]);
                    if ((double) (fabs(e[m]) + dd) == dd) break; // e[m] : offdiagonal part is relatively small to diagonal part
                }

                if (m != l) {
                    if (iter > 50) {
                        // too many iteration
                        eigenvalue_log_output << d[l] << "  > 50 iteration" << endl;
                        cout << d[l] << "  > 50 iteration" << endl;
                        goto label2;
                    }
                    iter = iter + 1;

                    // compute eigenvalue in submatrix as ks to shift matrix,
                    g = (d[l + 1] - d[l]) / (2.0 * e[l]); // Form shift.
                    r = pythag(g, 1.0);
                    g = d[m] - d[l] + e[l] / (g + SIGN(r, g)); // This is dm - ks.

                    s = c = 1.0;
                    p = 0.0;

                    for (i = m - 1; i >= l; i--) { // A plane rotation as in the original QL,
                        // followed by Givens rotations to restore tridiagonal form.
                        f = s * e[i];
                        b = c * e[i];
                        r = pythag(f, g);
                        e[i + 1] = r;
                        if (r == 0.0) {   // Recover from underflow.
                            d[i + 1] = d[i + 1] - p;
                            e[m] = 0.0;
                            break;
                        }
                        s = f / r;
                        c = g / r;
                        g = d[i + 1] - p;
                        r = (d[i] - g) * s + 2.0 * c * b;
                        p = s * r;
                        d[i + 1] = g + p;
                        g = c * r - b;
                    }
                    if (r == 0.0 && i >= l) continue;  // this is for if(r==0.0 condition)
                    d[l] = d[l] - p;
                    e[l] = g;
                    e[m] = 0.0;
                }
            } while (m != l);
            label2:;
        }
    }
}

int  partition(double * eigenvalue, int lo, int hi){
    double pivot = eigenvalue[hi];
    int i , j;
    double temp;
    i = lo;
    for(j=lo;j<hi;j++){
        if(eigenvalue[j] < pivot ){
            // swap eigenvalue[i] , eigenvalue[j]
            temp = eigenvalue[i];
            eigenvalue[i] = eigenvalue[j];
            eigenvalue[j] = temp;
            i ++ ;
        }
    }

    temp = eigenvalue[hi];
    eigenvalue[hi] = eigenvalue[i];
    eigenvalue[i] = temp;
    return i;
}

void quicksort(double * eigenvalue, int lo, int hi){
    int p ;
    if(lo < hi){
        p = partition(eigenvalue,lo,hi);
        quicksort(eigenvalue,lo,p-1);
        quicksort(eigenvalue,p+1,hi);
    }
}

void sortev(int begin_index, int end_index, double * eigenvalue){
    // use quick sort algorithm
    quicksort(eigenvalue,begin_index,end_index);
}

int partition_lambda (double * eigenvalue, int * repfind, int lo, int hi){
    double pivot = eigenvalue[hi];
    int i , j;
    double temp;
    int temp1;
    i = lo;
    for(j=lo;j<hi;j++){
        if(eigenvalue[j] < pivot ){
            // swap eigenvalue[i] , eigenvalue[j]
            temp = eigenvalue[i];
            eigenvalue[i] = eigenvalue[j];
            eigenvalue[j] = temp;

            temp1 = repfind[i];
            repfind[i] = repfind[j];
            repfind[j] = temp1;
            i ++ ;
        }
    }

    temp = eigenvalue[hi];
    eigenvalue[hi] = eigenvalue[i];
    eigenvalue[i] = temp;

    temp1 = repfind[hi];
    repfind[hi] = repfind[i];
    repfind[i] = temp1;

    return i;
}

void quicksort_lambda(double * eigenvalue, int * repfind, int lo, int hi){
    int p;
    if(lo < hi){
        p = partition_lambda(eigenvalue,repfind,lo,hi);
        quicksort_lambda(eigenvalue,repfind,lo,p-1);
        quicksort_lambda(eigenvalue,repfind,p+1,hi);
    }
}

void sortlambda(int begin_index, int end_index, double * lambda, int * repfind){
    // still use quick sort
    quicksort_lambda(lambda,repfind,begin_index,end_index);
}

void detector :: diagonalize(double * eigenvalue_list, int & numlam,  ofstream & eigenvalue_log_file){
    /*
    eigenvalue_list : size: total_dmat_size, store eigenvalue computed using Lanczos algorithm
     maxit: dimension of vector in Lanczos algorithm
     numlam: number of eigenvalue obtained
     lb: lower bound of eigenvalue want to compute
     ub: upper bound of eigenvalue want to compute
    */
    int i, j , k;
    int stopflag = 0;
    int iter = 0; // iteratioin
    numlam = 0;
    int numold;
    int iold;
    int verified = 0; // number of verified eigenvalue
    int nlev = total_dmat_size[0] ; // size of matrix we want to compute eigenvalue
    int lanc_matrix_size_ratio = 4;  // lanczos algorithm matrix size over input matrix size
    int lancdim = nlev * lanc_matrix_size_ratio; // maximum number of eigenvalue computed in Lanczos algorithm
    int maxit = lancdim - 1;  // maximum iteration in Lanczos algorithm (also size of Lanczos matrix size - 1 )
    int actit ;

    double minimum ;
    double diff ;
    int maximum = 4;  // maximum is maximum iteration we do to compute eigenvalue

    double vsum;
    double * lambda = new double [lancdim + 1];  // store potential eigenvalue. Results store in [1, lancdim] ( 0 is left as blank)
    int * repfind = new int [lancdim + 1] ; //
    double * v = new double [nlev]; // store vector v used in Lanczos algorithm
    double * alpha = new double [lancdim]; // diagonal part.
    double * beta = new double [lancdim]; // subddiagonal part. code for computing Lanczos
    double * d = new double[lancdim + 1];  // d is diagonal part of matrix after Lanczos algorithm and after we filter out part that loss orthogonality
    double * e = new double[lancdim + 1]; // e is subdiagonal part of matrix after we filter the matrix produced by Lanczos algorithm

    // for Lanczos MPI computation
    double * local_v;
    double * local_vold;
    double * local_w;
    double * send_v;
    double * recv_v;

    // prepare receive and send buffer for Lanczos algorithm.
    prepare_computation_for_Lanczos();

    allocate_space_for_vector_in_Lanczos(local_v, local_vold,local_w,send_v,recv_v);

    double error_spacing_between_eigenvalue = pow(10,-7); // we recognize eigenvalue's relative difference within this spacing as the same

    if(my_id == 0){
        eigenvalue_log_file << "Size of Hamiltonian to be diagonalized  " << nlev << endl;
        cout << "Size of Hamiltonian to be diagonalized  " << nlev << endl;
    }

    for(i=0;i<lancdim;i++){
        lambda[i] = 0;
    }

    // random number generator
    std::default_random_engine generator(time(NULL));
    std::uniform_real_distribution<double> distribution(0,1);  // random number in range [0,1]

    while(stopflag == 0){
         iter = iter + 1;
         if(my_id == 0){
             eigenvalue_log_file << "iteration:  "<<iter << endl;
             cout <<"iteration:  "<<iter<< endl;
         }
         //
         // compute initial Lanczos vector
         //
        if(my_id == 0){
            vsum = 0;
            for(i=0;i<nlev;i++){
                v[i] = distribution(generator) - 0.5;
                vsum = vsum + pow(v[i],2);
            }
            vsum = sqrt(vsum);
            for(i=0;i<nlev;i++){
                v[i] = v[i] / vsum;
            }
        }
        MPI_Scatterv(&v[0],dmatsize_each_process[0],dmatsize_offset_each_process[0],MPI_DOUBLE,
                     &local_v[0],dmatsize[0],MPI_DOUBLE,0,MPI_COMM_WORLD);
        //
        //  Call Lanczos algorithm then diagonalize tri-diagonal using tqli: implicit QL algorithm
        //
        // result store in alpha, beta: alpha: diagonal element, beta: subdiagonal element
        compute_Lanczos(alpha,beta,
                        local_v, local_vold, local_w,send_v,recv_v,
                        nlev,maxit);

        if(my_id == 0){
            cout << "Finish Lanczos run for this iteration: " << iter << endl;
            eigenvalue_log_file << "Finish Lanczos run for this iteration: " << iter << endl;
            i =0 ;
            // we make d[0], e[0] left as blank not to use. Record data in d[1], e[1]. This way we do not confuse when transfer fortran code to C++ code
            d[0] = 0;
            e[0] = 0;

            while(i<maxit and alpha[i] < pow(10,5)){
                d[i+1] = alpha[i];  // d , alpha is diagonal part tridiagonal matrix
                e[i+1] = beta[i];   // e , beta is off diagonal part of tridiagonal matrix
                i = i + 1;
            }
            // .. we filter possible Lanczos matrix according to size of alpha[i]
            actit = i;
            if(actit != maxit){
                cout << "Warning: Lanczos iteration truncated at  " << actit << " due to orthogonality" << endl;
                eigenvalue_log_file << "Warning: Lanczos iteration truncated at  " << actit << " due to orthogonality" << endl;
            }

            // call tqli : Implicit QL algorithm to get eigenvalue
            tqli(d,e,actit,eigenvalue_log_file);

            eigenvalue_log_file << "Completed first tqli this iteration:  " << iter << endl;
            cout <<  "Completed first tqli this iteration:  " << iter << endl;

            // sort eigenvalue
            sortev(1,actit,d);

            //
            // Find eigenvalue and multiplicity and save them
            //

            i = 1;
            numold = numlam;
            while(i<=actit){

                iold = i;  // i is index for newly found eigenvalue
                if(i < actit){
                    while( abs(d[i+1] - d[i])/abs(d[i]) < pow(10,-9)){
                        i = i + 1;
                    }
                }

                k = 1 ; // k is index for already found eigenvalue in lambda
                while(k <= numlam and abs(lambda[k] - d[i])/ abs(d[i]) > error_spacing_between_eigenvalue  ){
                    k = k + 1;
                }

                if(k<=numlam){
                    // found repeat eigenvalule
                    repfind[k] = 2;
                }
                else{
                    // this is newly found eigenvalue
                    numlam = numlam + 1;
                    lambda[numlam] = d[i];
                    if(i == iold){
                        repfind[numlam] = 1;
                    }
                    else{
                        repfind[numlam] = 2;
                    }
                }
                i = i + 1;
            }


            //
            // Run tqli with first alpha and beta removed to find eigenvalue
            //

            i = 1;
            while(i<actit){
                d[i] = alpha[i]; // d[i] = alpha[i]  : previously d[1] = alpha[0], so we omit first element
                e[i] = beta[i];
                i = i + 1;
            }
            d[actit] = 0;
            e[actit] = 0;

            tqli(d,e,actit,eigenvalue_log_file);

            eigenvalue_log_file << "Completed second tqli this iteration:  " << iter << endl;
            cout << "Complete second tqli this iteration: " << iter << endl;

            // sort env
            sortev(1,actit-1,d);
            for(i=1;i<=actit-1;i++){
                for(j=numold + 1; j<=numlam; j ++ ){
                    if(repfind[j] < 2 and abs(d[i] - lambda[j]) / abs(d[i]) < error_spacing_between_eigenvalue ){
                        repfind[j] = 0 ;
                    }
                }
            }

            // eliminate eigenvalue with repfind[k] == 0 from lambda list
            i = 1;
            j= 0 ;
            while( i <= numlam){
                if(repfind[i+j] == 0){
                    j = j + 1;
                    numlam = numlam - 1;
                }
                else{
                    repfind[i] = repfind[i+j];
                    lambda[i] = lambda[i+j];
                    i = i + 1;
                }
            }

            cout << "Removed suspicious eigenvalue from this iteration" << endl;
            eigenvalue_log_file << "Removed suspicious eigenvalue from this iteration" << endl;
            //
            // compute number of potential and verified eigenvalue and continue simulation if necessary
            //
            verified = 0;
            for(i=1;i<=numlam;i++){
                if(repfind[i] == 2) verified = verified + 1;
            }

            eigenvalue_log_file << numlam << " potential eigenvalue found so far " << endl;
            eigenvalue_log_file << verified << " out of " << nlev <<" eigenvalue verified so far " << endl;
            eigenvalue_log_file << endl;
            cout << numlam << " potential eigenvalue found so far" << endl;
            cout << verified << " out of " << nlev <<" eigenvalue verified so far " << endl;
            cout << endl;

            sortlambda(1,numlam,lambda,repfind);

            if(verified >= nlev or iter > maximum or numlam == nlev){
                stopflag = 1;
            }
        }

        // broadcase stopflag to let other process stop if stopflag == 1
        MPI_Bcast(&stopflag,1,MPI_INT,0,MPI_COMM_WORLD);

    }

    if(my_id == 0){
        if(verified >= nlev) {
            // eliminate eigenvalue with repfind[k] < 2, only verified eigenvalue left
            i = 1;
            j = 0;
            while (i <= numlam) {
                if (repfind[i + j] < 2) {
                    j = j + 1;
                    numlam = numlam - 1;
                } else {
                    repfind[i] = repfind[i + j];
                    lambda[i] = lambda[i + j];
                    i = i + 1;
                }
            }
            // now numlam == vereified
            if (verified > nlev) {
                // eliminate excess found eigenvalue
                eigenvalue_log_file << "Warning :  excess eigenvalue found " << endl;
                cout << "Warning:  excess eigenvalue found " << endl;


                minimum = pow(10, 6);
                while (verified > nlev) {

                    for (i = 2; i <= numlam; i++) {
                        diff = abs(lambda[i] - lambda[i - 1]) / abs(lambda[i]);
                        if (diff < minimum) {
                            minimum = diff;
                            k = i;
                        }
                    }
                    // eliminate element k
                    numlam = numlam - 1;
                    verified = verified - 1;

                    for (i = k; i <= numlam; i++) {
                        lambda[i] = lambda[i + 1];
                    }
                    lambda[numlam + 1] = 0;
                }

                eigenvalue_log_file << "Error criteria has to be raised to " << minimum
                                        << "to eliminate excess eigenvalue found " << endl;
                cout  << "Error criteria has to be raised to " << minimum
                          << "to eliminate excess eigenvalue found " << endl;

            }
        }
        else if (verified < nlev){

            eigenvalue_log_file <<"Warning:  only " << verified <<" eigenvalue verified" << endl;
            cout << "Warning:  only " << verified <<" eigenvalue verified" << endl;
            if(numlam < nlev){
                eigenvalue_log_file <<" Warning:  only " << numlam <<"eigenvalue found." << endl;
                cout << " Warning:  only " << numlam <<"eigenvalue found." << endl;
            }


            else if (numlam > nlev){
                minimum = pow(10,6);

                while(numlam > nlev){

                    for(i=2;i<=numlam;i++){
                        diff = abs(lambda[i] - lambda[i-1])/abs(lambda[i]);
                        if(diff < minimum){
                            minimum = diff;
                            k = i;
                        }
                    }
                    numlam = numlam -1 ;

                    if(repfind[k] == 2){
                        // switch element of k and k-1
                        diff = lambda[k-1];
                        i = repfind[k-1];
                        lambda[k-1] = lambda[k];
                        repfind[k-1] = repfind[k];
                        lambda[k] = diff;
                        repfind[k] = i;
                    }

                    if (repfind[k] == 2){
                        verified = verified - 1;
                    }

                    for(i=k;i<=numlam;i++){
                        lambda[i] = lambda[i + 1];
                        repfind[i] = repfind[ i + 1 ];
                    }

                    lambda[numlam + 1] = 0;

                }

                eigenvalue_log_file << "Error criteria has to be raised to " <<
                                    minimum << "to eliminate necessary number of unverified eigenvalue found " << endl;
                eigenvalue_log_file <<"Now verified eigeenvalue is:  " << verified << endl;

                cout << "Error criteria has to be raised to " <<
                     minimum << "to eliminate necessary number of unverified eigenvalue found " << endl;
                cout << "Now verified eigeenvalue is:  " << verified << endl;

            }
        }
    }

    if(my_id == 0){
        for(i=1;i<=numlam;i++){
            eigenvalue_list[i-1] = lambda[i];   // store found eigenvalue in eigenvalue_list.
        }
        for(i=numlam+1;i<=nlev;i++){
            eigenvalue_list[i-1] = 0;
        }
    }

    delete [] lambda;
    delete [] repfind;
    delete [] v;
    delete [] alpha;
    delete [] beta;
    delete [] d;
    delete [] e;

    // delete variable used in MPI
    delete [] local_v;
    delete [] local_vold;
    delete [] local_w;
    delete [] send_v;
    delete [] recv_v;

    delete [] to_recv_buffer_len;
    delete [] remoteVecCount;
    delete [] remoteVecPtr;
    delete [] remoteVecIndex;
    delete [] to_send_buffer_len;
    delete [] tosendVecCount;
    delete [] tosendVecPtr;
    delete [] tosendVecIndex;
    delete [] local_dirow;
    delete [] local_dicol;
}