//
// Created by phyzch on 5/27/21.
//
#include "util.h"
#include "system.h"

void compute_energy_range( double & Emin_return, double & Emax_return, double Energy, vector<double> & sorted_energy_list , int eigenstate_num){
    // return energy range : [Emin_return, Emax_return] , which estimate energy range within which we have number of state given by state_num around Energy.
    double min_energy_diff = 10000 ;
    double energy_diff ;
    int index_for_state;
    int begin_index;
    int end_index;
    int i ;
    int sorted_energy_list_len = sorted_energy_list.size();
    for(i=0; i< sorted_energy_list_len ; i++ ){
        energy_diff = abs (sorted_energy_list[i] - Energy) ;
        if(energy_diff < min_energy_diff){
            index_for_state = i;
            min_energy_diff = energy_diff ;
        }
    }
    begin_index = index_for_state - int (eigenstate_num / 2);
    end_index = index_for_state + int(eigenstate_num / 2);

    Emin_return = sorted_energy_list[begin_index];
    Emax_return = sorted_energy_list[end_index];

}


int MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector_given_energy_and_num(  int * dirow_list,  int * dicol_list,  double * dmat_list , vector<double> &  sorted_dmat_diagonal_list ,  int dmatsize  ,int dmatnum,
                                                                                  ofstream & Eigenvector_output,
                                                                                  double * &E , double ** &Matrix_X,
                                                                                  double energy_of_choice, int eigenstate_number){
    // E :[M] eigenvalue found .  Matrix_X: eigenvector found [M, dmatsize]. M is returned by function.

    // for CSR format, row index should be row_num + 1, which last element should be row_num + 1. (be careful about that)
    int i, j , k;
    vector<int>  dirow_CSR_form_fortran;
    vector<int>  dicol_CSR_form_fortran;
    vector<double>  sorted_dmat;

    // convert COO format which is the one we use in this program to CSR format:
    convert_COO_to_CSR(dirow_list,dicol_list,dmat_list,dmatsize,dmatnum,dirow_CSR_form_fortran,dicol_CSR_form_fortran,sorted_dmat);


    // to use dfeast_scsrev, we have to stall value in MKL_INT format
    // row format should be special for CSR form , pay attention to that.
    int * rows = new int   [dmatsize + 1];
    int * cols = new int   [dmatnum];
    double * val = new double  [dmatnum];
    for (i=0;i<dmatsize;i++){
        rows[i] = dirow_CSR_form_fortran[i];
    }
    rows[dmatsize] = dmatnum + 1;
    for (i=0;i<dmatnum;i++){
        cols[i] = dicol_CSR_form_fortran[i];
        val[i] = sorted_dmat[i];
    }


/* Declaration of FEAST variables */
    char          UPLO = 'F'; /* Type of matrix: (F means full matrix, L/U - lower/upper triangular part of matrix) */

    MKL_INT      fpm[128];      /* Array to pass parameters to Intel MKL Extended Eigensolvers */

    double       epsout;        /* Relative error on the trace */
    MKL_INT      loop;          /* Number of refinement loop */
    MKL_INT      info;          /* Errors */

    epsout = 0.0;
    loop = 0;
    info = 0;

    // Use Intel MKL extended eigensolver

    /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
    feastinit(
            fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
    );


    double dmat_list_min = 10000;
    double dmat_list_max = 0;

    // Below we estimate number of all eigenstates within energy range.
    // However, for large matrix, we only want sample small number of eigenvalues in matrix according to Boltzmann dist.
    MKL_INT      M0 = eigenstate_number * 2;  /* Initial guess for subspace dimension to be used */
    MKL_INT      M = eigenstate_number;             /* Total number of eigenvalues found in the interval */
    double Emin_of_choice;
    double Emax_of_choice;
    compute_energy_range(Emin_of_choice, Emax_of_choice, energy_of_choice, sorted_dmat_diagonal_list,eigenstate_number);

    if(M0 > dmatsize){
        M0 = dmatsize - 100 ;
    }





    E = new double [M0];         /* Eigenvalues */
    double * X = new double [M0 * dmatsize]; /* Eigenvectors */

    double * res = new double [M0];       /* Residual */


    /* Initialize matrix X */
    for (i=0; i<M0*dmatsize; i++)
    {
        X[i] = 0.0;
    }

    for(i=0;i<M0;i++){
        E[i] = 0;
    }


    printf("Sparse matrix size for solving eigenvector %i\n", (int)dmatsize);
    /* Search interval [Emin,Emax] */
    // Emin, Emax is generated from energy we input and number of state we want
    printf("number of states try to solve:  %d \n" , eigenstate_number);
    printf("Energy we want to estimate eigenstate: %.5e \n" ,energy_of_choice );
    printf("M value estimated %d \n",M0);
    printf("From energy and eigenstate_number: Search interval [ %.5e, %.5e  ]  \n", Emin_of_choice, Emax_of_choice);




    //fpm[0] =  1; /* Extended Eigensolver routines print runtime status to the screen. */

    /* Step 2. Solve the standard Ax = ex eigenvalue problem. */
    printf("Testing dfeast_scsrev routine:\n");
    dfeast_scsrev(
            &UPLO,   /* IN: UPLO = 'F', stores the full matrix */
            &dmatsize,      /* IN: Size of the problem */
            val,     /* IN: CSR matrix A, values of non-zero elements */
            rows,    /* IN: CSR matrix A, index of the first non-zero element in row */
            cols,    /* IN: CSR matrix A, columns indices for each non-zero element */
            fpm,     /* IN/OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
            &epsout, /* OUT: Relative error of on the trace */
            &loop,   /* OUT: Contains the number of refinement loop executed */
            &Emin_of_choice,   /* IN: Lower bound of search interval */
            &Emax_of_choice,   /* IN: Upper bound of search interval */
            &M0,     /* IN: The initial guess for subspace dimension to be used. */
            E,       /* OUT: The first M entries of Eigenvalues */
            X,       /* IN/OUT: The first M entries of Eigenvectors */
            &M,      /* OUT: The total number of eigenvalues found in the interval */
            res,     /* OUT: The first M components contain the relative residual vector */
            &info    /* OUT: Error code */
    );
    // X eigenvector have the form: M(row) * dmatsize(column). Here M is number of eigenvector found.
    printf("FEAST OUTPUT INFO %d \n",info);



    printf("Number of eigenvalues found %d \n", M);

    if(M<0 or M > dmatsize){
        printf("Bad M value returned from dfeast_scsrev:  Exit.");
        Eigenvector_output.close();
        MPI_Abort(MPI_COMM_WORLD,-22);
    }

    if ( info != 0 )
    {
        printf("Routine dfeast_scsrev returns code of ERROR: %i", (int)info);
        Eigenvector_output.close();
        MPI_Abort(MPI_COMM_WORLD, -21);
    }

    // check normalization and orthogonality of eigenvector by computing Y =  X^{transpose} * X - I
    /* Y=(X')*X-I */
    Matrix_X = new double * [M];
    for(i=0;i<M;i++){
        Matrix_X[i] = new double [dmatsize];
    }
    int index = 0;
    for(i=0;i<M;i++){
        for(j=0;j<dmatsize;j++){
            Matrix_X[i][j] = X[index];
            index ++ ;
        }
    }
    printf("Finish constructing Matrix_X \n");

    double ** Y = new double * [M];
    for(i=0;i<M;i++){
        Y[i] = new double [M];
    }
    for(i=0;i<M;i++){
        for(j=0;j<M;j++){
            Y[i][j] = 0;
            for (k = 0; k<dmatsize;k++){
                Y[i][j] = Y[i][j] + Matrix_X[i][k] * Matrix_X[j][k];
            }
        }
    }

    for(i=0;i<M;i++){
        Y[i][i] = Y[i][i] - 1;
    }
    printf("Finsih computing Y \n");
    double small_value = 0;
    for(i=0;i<M;i++){
        for(j=0;j<M;j++){
            if (abs(Y[i][j]) > small_value){
                small_value = abs(Y[i][j]);
            }
        }
    }
    printf("Maximum value in X*X - I is %f \n" , small_value);

    for(i=0;i<M;i++){
        // output eigenvalue
        Eigenvector_output << E[i] << endl;
        for(j=0;j<dmatsize;j++){
            //output eigenvector
            Eigenvector_output << Matrix_X[i][j] << " ";
        }
        Eigenvector_output << endl;
    }

    // delete dynamics array X. (if not used)
    delete [] X ;
    delete [] rows;
    delete [] cols;
    delete [] val;


    for(i=0;i<M;i++){
        delete [] Y[i];
    }
    delete [] Y;
    delete [] res;

    return M;
}

void allocate_diagonalization_energy_range_for_diff_proc( vector<double> & sorted_dmat_diagonal_part ) {
    // we divide energy range into different part and let different process work in parallel to solve eigenvalue and eigenvector.
    int i, j , k ;
    int index1;
    int total_state_num_in_range = 0 ;
    int dmat_diagonal_num = sorted_dmat_diagonal_part.size();


    bool initial_index_bool = false;
    bool end_index_bool = false;
    int initial_index; // initial index in sorted diagonal energy part within energy range
    int end_index;  // end index in sorted diagonal energy part within energy range
    for(index1=0; index1 < dmat_diagonal_num ; index1++){
        if(sorted_dmat_diagonal_part[index1] > Emin and sorted_dmat_diagonal_part[index1] < Emax){
            total_state_num_in_range = total_state_num_in_range + 1 ;
        }
        if(sorted_dmat_diagonal_part[index1] > Emin and !initial_index_bool ){
            initial_index_bool = true;
            initial_index = index1;
        }
        if(sorted_dmat_diagonal_part[index1] >= Emax and !end_index_bool){
            end_index_bool = true;
            end_index = index1 - 1 ;
        }
    }

    int num_of_state_for_solving_eigenvec_per_proc = int (total_state_num_in_range / num_proc ) ;
    vector<int> block ;
    double block_Emin ;
    double block_Emax;

    for(i=0;i< num_proc   ; i++ ){
        block.push_back( i * num_of_state_for_solving_eigenvec_per_proc );
    }
    block.push_back(total_state_num_in_range);

    vector<double> block_energy_range ;
    block_energy_range.push_back(Emin);

    vector<double> sorted_dmat_diagonal_part_slice ;
    for(i=initial_index ; i <= end_index; i++ ){
        sorted_dmat_diagonal_part_slice.push_back(sorted_dmat_diagonal_part[i]);
    }

    for(i=1 ; i< num_proc  ; i++ ){
        block_energy_range.push_back (  sorted_dmat_diagonal_part_slice[ block[i] ] ) ;
    }
    block_energy_range.push_back(Emax);

    Emin_for_core = block_energy_range[ my_id ];
    Emax_for_core = block_energy_range[ my_id + 1 ];
    if(my_id == 0){
        cout << "total number of eigenstate estimated:  " << total_state_num_in_range << endl;
        cout << "each process estimate to solve number of eigenstate:  " << num_of_state_for_solving_eigenvec_per_proc << endl;
        cout <<" Each process solve block index below : " << endl;
        for(i=0;i<num_proc + 1 ;i++){
            cout << block[i] << " ";
        }
        cout << endl;

        cout <<"energy range for block " << endl;
        for(i=0;i<num_proc + 1 ; i++){
            cout << block_energy_range[i] << "  ";
        }
        cout << endl;
    }
}

int MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector_divide_by_part(int * dirow_list,  int * dicol_list,  double * dmat_list , int dmatsize  ,int dmatnum,
                                                                          vector<double> & dmat_diagonal_part ,
                                                                          vector<double> & Eigenvalue_list , vector<vector<double>> & Eigenvector_list,
                                                                          double Emin_of_choice, double Emax_of_choice){
    int i, j , k ;
    int index1;
    int total_state_num_in_range = 0 ;
    int dmat_diagonal_num = dmat_diagonal_part.size();

    int num_of_state_for_solving_eigenvec = 20;

    int eigenstate_num_solved  =  0;
    int total_eigenstate_num_solved = 0;


    vector<double> sorted_dmat_diagonal_part = dmat_diagonal_part ;
    sort(sorted_dmat_diagonal_part.begin() , sorted_dmat_diagonal_part.end());

    bool initial_index_bool = false;
    bool end_index_bool = false;
    int initial_index; // initial index in sorted diagonal energy part within energy range
    int end_index;  // end index in sorted diagonal energy part within energy range
    for(index1=0; index1 < dmat_diagonal_num ; index1++){
        if(sorted_dmat_diagonal_part[index1] > Emin_of_choice and sorted_dmat_diagonal_part[index1] < Emax_of_choice){
            total_state_num_in_range = total_state_num_in_range + 1 ;
        }
        if(sorted_dmat_diagonal_part[index1] > Emin_of_choice and !initial_index_bool ){
            initial_index_bool = true;
            initial_index = index1;
        }
        if(sorted_dmat_diagonal_part[index1] >= Emax_of_choice and !end_index_bool){
            end_index_bool = true;
            end_index = index1;
        }
    }

    int block_number = int(total_state_num_in_range / num_of_state_for_solving_eigenvec) + 1;
    vector<int> block ;
    double block_Emin ;
    double block_Emax;

    for(i=0;i<block_number -1 ; i++ ){
        block.push_back( i * num_of_state_for_solving_eigenvec );
    }
    block.push_back(total_state_num_in_range);

    vector<double> block_energy_range ;
    block_energy_range.push_back(Emin_of_choice);

    vector<double> sorted_dmat_diagonal_part_slice ;
    for(i=initial_index ; i < end_index; i++ ){
        sorted_dmat_diagonal_part_slice.push_back(sorted_dmat_diagonal_part[i]);
    }

    for(i=1 ; i< block_number - 1 ; i++ ){
        block_energy_range.push_back (  sorted_dmat_diagonal_part_slice[ block[i] ] ) ;
    }
    block_energy_range.push_back(Emax_of_choice);

//    cout << "Block_number :  " << block_number << endl;
//    cout << "Block" << endl;
//    for(i=0;i<block_number;i++){
//        cout << block[i] << " ";
//    }
//    cout << endl;
//    cout << "Block energy: " << endl;
//    for(i=0;i<block_number;i++){
//        cout << block_energy_range[i] << " ";
//    }
//    cout << endl;

    for (i=0; i<block_number-1 ; i++ ){
        block_Emin = block_energy_range[i];
        block_Emax = block_energy_range[i + 1];
        vector<double> Eigenvalue_temp ;
        vector<vector<double>> Eigenvector_temp ;
        eigenstate_num_solved = MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector_submodule(dirow_list ,dicol_list ,dmat_list , dmatsize,  dmatnum ,
                                                                                       dmat_diagonal_part , Eigenvalue_temp , Eigenvector_temp , block_Emin, block_Emax );
        total_eigenstate_num_solved = total_eigenstate_num_solved + eigenstate_num_solved ;

        for(j = 0;j < eigenstate_num_solved; j++ ){
            Eigenvalue_list.push_back(Eigenvalue_temp[j]);
            vector<double> v = Eigenvector_temp[j] ;
            Eigenvector_list.push_back(v);
        }
    }

    printf("total number of eigenstate solved in proc %d :   %d \n" , my_id ,total_eigenstate_num_solved );


    return total_eigenstate_num_solved ;
}

void detector::Broadcast_eigenstate_and_eigenvalue(vector<double> & Eigenvalue_temp, vector<vector<double>> & Eigenvector_temp){
    int i,j , k;
    // broadcast total number of eigenstate found to other process
    int eigenstate_num_in_all_proc;
    int * eigenstate_num_in_each_proc = new int [num_proc];
    int * eigenstate_num_disp = new int [num_proc];
    MPI_Allgather(&eigenstate_num, 1, MPI_INT, &eigenstate_num_in_each_proc[0], 1, MPI_INT, MPI_COMM_WORLD );

    eigenstate_num_disp[0] = 0;
    for(i=1;i<num_proc;i++){
        eigenstate_num_disp[i] = eigenstate_num_disp[i - 1] + eigenstate_num_in_each_proc[i-1];
    }

    MPI_Allreduce(&eigenstate_num,&eigenstate_num_in_all_proc,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    eigenstate_num = eigenstate_num_in_all_proc;


    // other process allocate space for eigenvalue and eigenvector
    Eigenvalue_list = new double[eigenstate_num];
    Eigenstate_list = new double * [eigenstate_num];
    for(i=0;i<eigenstate_num ; i++){
        Eigenstate_list[i] = new double [total_dmat_size[0]];
    }

    double * Eigenstate_recv_list = new double [total_dmat_size[0]];

    // AllGather Eigenvalue_list
    MPI_Allgatherv(&Eigenvalue_temp[0], eigenstate_num_in_each_proc[my_id], MPI_DOUBLE, &Eigenvalue_list[0], eigenstate_num_in_each_proc, eigenstate_num_disp, MPI_DOUBLE, MPI_COMM_WORLD);

    // AllGather Eigenstate_list
    int eigenstate_index = 0;
    for(i=0;i<num_proc;i++){
        for(j=0;j<eigenstate_num_in_each_proc[i];j++){
            if(my_id == i){
                for(k=0;k<total_dmat_size[0]; k++){
                    Eigenstate_recv_list[k] = Eigenvector_temp[j][k];
                }
            }
            MPI_Bcast(&Eigenstate_recv_list[0],total_dmat_size[0], MPI_DOUBLE, i, MPI_COMM_WORLD);
            for(k=0;k<total_dmat_size[0];k++){
                Eigenstate_list[eigenstate_index][k] = Eigenstate_recv_list[k] ;
            }

            eigenstate_index = eigenstate_index + 1;
        }
    }

    delete [] eigenstate_num_in_each_proc;
    delete [] eigenstate_num_disp;
    delete [] Eigenstate_recv_list;

}

void construct_energy_window_for_eigenstate(int nmode, double * mfreq, double eigen_energy, double energy_window_for_eigenstate,
                                            vector<double> & Energy_range_min_list,
                                            vector<double> & Energy_range_max_list){
    // To compute OTOC for particular eigenstate, we do not have to resolve all eigenstate and eigenvalue.
    vector<double> Energy_range_list;
    int energy_range_list_size;
    int index;
    int i , j;
    double connected_eigen_state_energy;
    for (i=0; i<nmode; i++){
        connected_eigen_state_energy = eigen_energy - mfreq[i] ;
        if(connected_eigen_state_energy > 0){
            Energy_range_list.push_back(connected_eigen_state_energy);
        }
    }
    for(i=0; i<nmode; i++ ){
        connected_eigen_state_energy = eigen_energy + mfreq[i] ;
        Energy_range_list.push_back(connected_eigen_state_energy) ;
    }

    sort(Energy_range_list.begin(), Energy_range_list.end());

    // construct an energy window around these energy to solve for possible eigenstate .
    energy_range_list_size = Energy_range_list.size();
    for(i=0;i<energy_range_list_size; i++ ){
        Energy_range_min_list.push_back(Energy_range_list[i] - energy_window_for_eigenstate);
        Energy_range_max_list.push_back(Energy_range_list[i] + energy_window_for_eigenstate);
    }

    i = 0;
    while(1){
        if(i == energy_range_list_size - 1 ){
            break;
        }
        if(Energy_range_max_list[i] > Energy_range_min_list[i+1]){
            Energy_range_max_list[i] = Energy_range_max_list[i+1];
            Energy_range_min_list.erase(Energy_range_min_list.begin() + i + 1);
            Energy_range_max_list.erase(Energy_range_max_list.begin() + i + 1);
            energy_range_list_size = energy_range_list_size - 1 ;
        }
        else{
            i = i + 1;
        }
    }


}