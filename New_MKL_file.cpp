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