//
// Created by phyzch on 5/28/21.
//

//
// Created by phyzch on 3/29/21.
//
#include"system.h"
#include"util.h"
using namespace std;

// Interval for solving eigenvector and eigenvalue.
/* Lower/upper bound of search interval [Emin,Emax] */
double Emin = 0; // default value
double Emax = 100; // default value
double Emin2 = 0;
double Emax2 = 100;

struct Matrix_COO_element{
    int row;
    int col;
    double value;
    Matrix_COO_element(int row1, int col1, double value1){
        value = value1;
        row = row1;
        col = col1;
    }
};

vector<Matrix_COO_element> merge_sort(const vector<Matrix_COO_element> & v1, const vector<Matrix_COO_element> & v2 ){
    int size1= v1.size();
    int size2= v2.size();
    if(size1 == 0) return v2;
    if(size2 == 0) return v1;
    int v1_index = 0;
    int v2_index = 0;

    int i;
    int row_index_difference;
    int column_index_difference ;
    vector <Matrix_COO_element> v3;
    while(v1_index<size1 and v2_index<size2){
        row_index_difference = v1[v1_index].row - v2[v2_index].row;
        if(row_index_difference < 0){
            v3.push_back(v1[v1_index]);
            v1_index++;
        }
        else if (row_index_difference > 0){
            v3.push_back(v2[v2_index]);
            v2_index++;
        }
        else{
            // compare column index
            column_index_difference = v1[v1_index].col - v2[v2_index].col;
            if(column_index_difference < 0){
                v3.push_back(v1[v1_index]);
                v1_index ++ ;
            }
            else{
                v3.push_back(v2[v2_index]);
                v2_index++;
            }
        }
    }
    if(v1_index<size1){
        for(i=v1_index; i <size1; i++){
            v3.push_back(v1[i]);
        }
    }
    if(v2_index<size2){
        for(i=v2_index;i<size2;i++){
            v3.push_back(v2[i]);
        }
    }
    return v3;
}

void convert_COO_to_CSR(const int * dirow_list, const int * dicol_list, const double * dmat_list , int dmatsize, int dmatnum,
                        vector<int> & dirow_CSR_form_fortran, vector<int> & dicol_CSR_form_fortran, vector<double> & sorted_dmat ){
    // convert COO (row , column, element) form of sparse matrix to CSR form of matrix.
    //CSR form matrix is used in MKL dfeast_scsrev function to solve eigenvalue and eigenvector in a given range:
    // See https://software.intel.com/content/www/us/en/develop/articles/introduction-to-the-intel-mkl-extended-eigensolver.html
    // use sort algorithm to sort Matrix according to their
    int i, j;

    // construct list of Matrix_COO_element for merge sort
    vector<Matrix_COO_element> list_to_sort ;
    for (i=0;i<dmatnum;i++){
        Matrix_COO_element element(dirow_list[i], dicol_list[i], dmat_list[i]);
        list_to_sort.push_back(element);
    }

    // start merge sorting
    vector <vector<Matrix_COO_element>> List_for_list_old;
    vector<vector<Matrix_COO_element>> * old_ptr = & List_for_list_old;
    vector<vector<Matrix_COO_element>> List_for_list_new;
    vector<vector<Matrix_COO_element>> * new_ptr = & List_for_list_new;
    vector<vector<Matrix_COO_element>> * list_ptr_3;
    vector<Matrix_COO_element> v3;

    int list_size = dmatnum;
    for(i=0;i<dmatnum;i++){
        vector<Matrix_COO_element> small_list ;
        small_list.push_back(list_to_sort[i]);
        List_for_list_old.push_back(small_list);
    }

    // sort with order as : first compare row_index to have min_row, then compare column index to have min_column
    while(list_size > 1){
        for(i=0; i+1 < list_size; i=i+2){
            v3 = merge_sort( (*old_ptr)[i], (*old_ptr)[i+1] );
            (*new_ptr).push_back(v3);
        }
        if(list_size % 2 ==1){
            (*new_ptr).push_back((*old_ptr)[list_size - 1]);
        }
        list_size = (list_size + 1) /2;
        //exchange two list
        (*old_ptr).clear();
        (*old_ptr).shrink_to_fit();

        list_ptr_3 = old_ptr;
        old_ptr = new_ptr;
        new_ptr = list_ptr_3;
    }

    list_to_sort.clear();
    list_to_sort = (*old_ptr)[0];  // sorting result store in *old_ptr

    // see dexample_sparse_c.c for input of matrix row and column. To use fortran package, we have to take initial index for row and column to 1.

    int old_row = -1;
    int row_index;
    for (i=0;i<dmatnum;i++){
        dicol_CSR_form_fortran.push_back(list_to_sort[i].col + 1); // for fortran input for  dfeast_scsrev, index start with 1
        sorted_dmat.push_back(list_to_sort[i].value);
        row_index = list_to_sort[i].row;
        if(row_index != old_row){
            dirow_CSR_form_fortran.push_back( i + 1);
            old_row = row_index;
        }
    }

    // dirow we get should have size dmatsize.
    if(dirow_CSR_form_fortran.size()!= dmatsize){
        cout<<"Error in COO to CSR, CSR row does not have size dmatsize. " << endl;

        MPI_Abort(MPI_COMM_WORLD,-20);
    }

}




int MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector(  int * dirow_list,  int * dicol_list,  double * dmat_list , int dmatsize  ,int dmatnum,
                                                             vector<double> & dmat_diagonal_part ,
                                                             ofstream & Eigenvector_output,
                                                             vector<double> & E , vector<vector<double>> & Matrix_X,
                                                             double Emin_of_choice, double Emax_of_choice ){
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

    MKL_INT      L = 0;
    double dmat_list_min = 10000;
    double dmat_list_max = 0;

    int index1;
    int dmat_diagonal_num = dmat_diagonal_part.size();
    for( index1 = 0 ; index1 < dmat_diagonal_num ; index1++ ){
        if(dmat_diagonal_part[index1] < dmat_list_min){
            dmat_list_min = dmat_diagonal_part [index1];
        }
        if(dmat_diagonal_part[index1] > dmat_list_max){
            dmat_list_max = dmat_diagonal_part[index1];
        }
    }


    cout << "Min energy for system:  " << dmat_list_min << endl;
    cout <<" Max energy for system:   "<< dmat_list_max << endl;
    cout <<"dmatsize :   " << dmatsize << "   dmatnum:   " << dmatnum << endl;

    // Below we estimate number of all eigenstates within energy range.
    // However, for large matrix, we only want sample small number of eigenvalues in matrix according to Boltzmann dist.
    for(index1=0; index1 < dmat_diagonal_num ; index1++){
        if(dmat_diagonal_part[index1] > Emin_of_choice and dmat_diagonal_part[index1] < Emax_of_choice){
            L = L + 1 ;
        }
    }

    L = int( L * 1.5 );
    if(L > dmatsize){
        L = dmatsize - 100 ;
    }

    MKL_INT      M0 = 0;            /* Initial guess for subspace dimension to be used */
    MKL_INT      M = 0;             /* Total number of eigenvalues found in the interval */

    double * E_array = new double [L];         /* Eigenvalues */
    double * X = new double [L * dmatsize]; /* Eigenvectors */
//    double       X[dmatsize * dmatsize];        /* Eigenvectors */

    double * res = new double [dmatsize];       /* Residual */
    MKL_INT      info;          /* Errors */

    /* Initialize matrix X */
    for (i=0; i< L *dmatsize; i++)
    {
        X[i] = 0.0;
    }

    for(i=0;i<L ;i++){
        E_array[i] = 0;
    }

    M0   = L;
    M    = L;

    printf("Sparse matrix size for solving eigenvector %i\n", (int)dmatsize);
    /* Search interval [Emin,Emax] */
    // Emin, Emax is read from input file
    printf("Search interval [ %.5e, %.5e  ]  \n", Emin_of_choice, Emax_of_choice);
    printf("M value estimated %d \n",M0);


    loop = 0;
    info = 0;
    epsout = 0.0;

    // Use Intel MKL extended eigensolver

    /* Step 1. Call  FEASTINIT to define the default values for the input FEAST parameters */
    feastinit(
            fpm /* OUT: Array is used to pass parameters to Intel MKL Extended Eigensolvers */
    );

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
            &Emax_of_choice ,   /* IN: Upper bound of search interval */
            &M0,     /* IN: The initial guess for subspace dimension to be used. */
            E_array,       /* OUT: The first M entries of Eigenvalues */
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

    // construct Matrix X
    for(i=0;i<M;i++){
        vector<double> v;
        for(j=0;j<dmatsize; j++ ){
            v.push_back(0);
        }
        Matrix_X.push_back(v);
    }
    int index = 0;
    for(i=0;i<M;i++){
        for(j=0;j<dmatsize;j++){
            Matrix_X[i][j] = X[index];
            index ++ ;
        }
    }
    printf("Finish constructing Matrix_X \n");

    // construct eigenvalue E
    for(i=0; i<M; i++){
        E.push_back(E_array[i]);
    }


    // check normalization and orthogonality of eigenvector by computing Y =  X^{transpose} * X - I
    /* Y=(X')*X-I */
//    double ** Y = new double * [M];
//    for(i=0;i<M;i++){
//        Y[i] = new double [M];
//    }

//    for(i=0;i<M;i++){
//        for(j=0;j<M;j++){
//            Y[i][j] = 0;
//            for (k = 0; k<dmatsize;k++){
//                Y[i][j] = Y[i][j] + Matrix_X[i][k] * Matrix_X[j][k];
//            }
//        }
//    }
//
//    for(i=0;i<M;i++){
//        Y[i][i] = Y[i][i] - 1;
//    }
//    printf("Finsih computing Y \n");
//    double small_value = 0;
//    for(i=0;i<M;i++){
//        for(j=0;j<M;j++){
//            if (abs(Y[i][j]) > small_value){
//                small_value = abs(Y[i][j]);
//            }
//        }
//    }
//    printf("Maximum value in X*X - I is %f \n" , small_value);

    for(i=0;i<M;i++){
        // output eigenvalue
        Eigenvector_output << E[i] << endl;
//        for(j=0;j<dmatsize;j++){
//            //output eigenvector
//            Eigenvector_output << Matrix_X[i][j] << " ";
//        }
//        Eigenvector_output << endl;
    }

    // delete dynamics array X. (if not used)
    delete [] X ;
    delete [] rows;
    delete [] cols;
    delete [] val;
    delete [] res;
    delete [] E_array ;

//    for(i=0;i<M;i++){
//        delete [] Y[i];
//    }
//    delete [] Y;

    return M;
}