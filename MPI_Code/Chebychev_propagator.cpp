//
// Created by phyzch on 1/26/21.
//
#include "../util.h"
#include "../system.h"

// instead of using SUR algorithm, we use Chebychev algorithm
void detector::prepare_evolution_Chebychev_method(ofstream & log){
    // See : https://doi.org/10.1063/1.448136
    // except to call this function, we also have to call prepare_evolution function to prepare parallel computation : remote_Vec_count etc.
    // this function should be called after prepare_evolution() function
    int i, j;
    double max_Hamiltonian, min_Hamiltonian;
    max_Hamiltonian = total_dmat[0][0];
    min_Hamiltonian = total_dmat[0][0];
    for(i=0;i<total_dmat_size[0];i++){
        if (max_Hamiltonian < total_dmat[0][i]){
            max_Hamiltonian = total_dmat[0][i];
        }
        if (min_Hamiltonian > total_dmat[0][i]){
            min_Hamiltonian = total_dmat[0][i];
        }
    }

    Chebychev_e0 = ( max_Hamiltonian + min_Hamiltonian )/ 2;
    Chebychev_R = (max_Hamiltonian - min_Hamiltonian ) * 0.55 ;

    shifted_total_dmat = new double * [stlnum];
    shifted_dmat = new double * [stlnum];
    for(i=0;i<stlnum;i++){
        shifted_total_dmat[i] = new double [total_dmat_num[i]];
        shifted_dmat[i] = new double [dmatnum[i]];
    }

    for(i=0;i<total_dmat_num[0];i++){
        if(total_dirow[0][i] == total_dicol[0][i]){
            shifted_total_dmat[0][i] = total_dmat[0][i] - Chebychev_e0;
        }
        else{
            shifted_total_dmat[0][i] = total_dmat[0][i];
        }
        shifted_total_dmat[0][i] = shifted_total_dmat[0][i] / Chebychev_R;
    }

    for(i=0;i<dmatnum[0];i++){
        if(dirow[0][i] == dicol[0][i]){
            shifted_dmat[0][i] = dmat[0][i] - Chebychev_e0;
        }
        else{
            shifted_dmat[0][i] = dmat[0][i];
        }
        shifted_dmat[0][i] = shifted_dmat[0][i] / Chebychev_R;
    }

    Chebychev_Rt = Chebychev_R * cf;
    Chebychev_expr = std::cos(cf * Chebychev_e0);
    Chebychev_expi = - std::sin(cf * Chebychev_e0); // Here it is imaginary part of exp(-ie_{0} dt )

    Nchev = int(Chebychev_Rt) + 5; // J_{k}(R * dt) will decrease exponentially when k >= (R * dt)
    if(my_id == 0){
        cout << "Using Chebychev method.  R * dt=   " << Chebychev_Rt << endl;
        log << "Using Chebychev method.  R * dt=   " << Chebychev_Rt << endl;
        cout << "Using Chebychev method. order of Chebychev polynomial:    " << Nchev << endl;
        log << "Using Chebychev method. order of Chebychev polynomial:     " << Nchev << endl;
    }
    Bessel_function_array = new double [Nchev + 1];
    for(i=0;i<=Nchev;i++){
        Bessel_function_array[i] = std::cyl_bessel_j(i,Chebychev_Rt);  // J_{k}(R*t)
    }

    Chebychev_polyn = new vector<double> [6];
    for(i=0;i<6;i++){
        vector<double> v (dmatsize[0] + to_recv_buffer_len[0],0);
        Chebychev_polyn[i]=v;
    }

    send_polyn = new double * [6];
    recv_polyn = new double * [6];
    for(i=0;i<6;i++){
        send_polyn[i] = new double [to_send_buffer_len[0]];
        recv_polyn[i] = new double [to_recv_buffer_len[0]];
    }
}


void detector::update_polyn23() {
    // update Chebychev_polyn[2] and Chebychev_polyn[3]
    int i;
    int vsize;
    // collect data for send_buffer.
    vsize = total_dmat_size[0]/num_proc;
    for(i=0;i<to_send_buffer_len[0];i++){
        send_polyn[2][i] = Chebychev_polyn[2][tosendVecIndex[0][i] - my_id * vsize];
        send_polyn[3][i] = Chebychev_polyn[3][tosendVecIndex[0][i] - my_id * vsize];
    }
    MPI_Alltoallv(&send_polyn[2][0],tosendVecCount[0],tosendVecPtr[0],MPI_DOUBLE,
                  &recv_polyn[2][0],remoteVecCount[0],remoteVecPtr[0],MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Alltoallv(&send_polyn[3][0],tosendVecCount[0],tosendVecPtr[0],MPI_DOUBLE,
                  &recv_polyn[3][0],remoteVecCount[0],remoteVecPtr[0],MPI_DOUBLE,MPI_COMM_WORLD);
    for(i=0;i<to_recv_buffer_len[0];i++){
        Chebychev_polyn[2][ i + dmatsize[0]]  = recv_polyn[2][i];
        Chebychev_polyn[3][ i + dmatsize[0]]  = recv_polyn[3][i];
    }
}


void detector:: Chebychev_method_one_step(){
    int i,j,k,m;
    int irow,icol;
    double bess;
    double air, aii;
    int nearby_state_list_size = nearby_state_index.size();
    update_dx(nearby_state_list_size);
    update_dy(nearby_state_list_size);


    for(m=0;m<nearby_state_list_size;m++){
        vector<double> & creal = xd[m];
        vector<double> & cimag = yd[m];

        // Initialize Chebychev polynomial * |phi>. Here Chebychev_polyn is Chebychev polynomial of Hamiltonian act upon wave function
        // THis part also have to be changed in parallel code
        for(i=0;i<dmatsize[0] + to_recv_buffer_len[0];i++){
            Chebychev_polyn[0][i] = creal[i];
            Chebychev_polyn[1][i] = cimag[i];
            Chebychev_polyn[2][i] = 0;
            Chebychev_polyn[3][i] = 0;
            Chebychev_polyn[4][i] = 0;
            Chebychev_polyn[5][i] = 0;
        }

        //zeroth order: C0 = 1 , J0(Rt) * exp(-i* e0 * dt) * wavefunction
        bess = Bessel_function_array[0];
        air = Chebychev_expr * bess;
        aii = Chebychev_expi * bess;

        for(i=0;i<dmatsize[0];i++){
            creal[i] = air * Chebychev_polyn[0][i] - aii * Chebychev_polyn[1][i];
            cimag[i] = air * Chebychev_polyn[1][i] + aii * Chebychev_polyn[0][i];
        }

        // first order
        bess = Bessel_function_array[1];
        air = 2 * bess * Chebychev_expr;
        aii = 2 * bess * Chebychev_expi;

        for(i=0;i<dmatnum[0];i++){
            // make sure when we compute off-diagonal matrix, we record both symmetric and asymmetric part
            irow = local_dirow[0][i];
            icol = local_dicol[0][i]; // compute to point to colindex in
            Chebychev_polyn[3][irow] = Chebychev_polyn[3][irow] - shifted_dmat[0][i] * Chebychev_polyn[0][icol];
            Chebychev_polyn[2][irow] = Chebychev_polyn[2][irow] + shifted_dmat[0][i] *Chebychev_polyn[1][icol];
        }
        update_polyn23();
        for(i=0;i<dmatsize[0];i++){
            creal[i] = creal[i] +  air * Chebychev_polyn[2][i] - aii * Chebychev_polyn[3][i];
            cimag[i] = cimag[i] + air * Chebychev_polyn[3][i] + aii * Chebychev_polyn[2][i];
        }

        // Remaining terms:
        for(k=2;k<=Nchev;k++){
            bess = Bessel_function_array[k];
            air = 2 * bess * Chebychev_expr;
            aii = 2 * bess * Chebychev_expi;

            // use Chebychev polynomial recursion relationship. T_{k+2}(x) = T_{k+1}(x) *2x + T_{k}(x)
            for(i=0;i<dmatnum[0];i++){
                irow = local_dirow[0][i];
                icol = local_dicol[0][i];
                Chebychev_polyn[5][irow] = Chebychev_polyn[5][irow] - 2 * shifted_dmat[0][i] * Chebychev_polyn[2][icol];
                Chebychev_polyn[4][irow] = Chebychev_polyn[4][irow] + 2 * shifted_dmat[0][i] * Chebychev_polyn[3][icol];
            }
            for(i=0;i<dmatsize[0];i++){
                Chebychev_polyn[4][i] = Chebychev_polyn[4][i] + Chebychev_polyn[0][i];
                Chebychev_polyn[5][i] = Chebychev_polyn[5][i] + Chebychev_polyn[1][i];
            }
            // update dx, dy
            for(i=0;i<dmatsize[0];i++){
                creal[i] = creal[i] + air * Chebychev_polyn[4][i] - aii * Chebychev_polyn[5][i];
                cimag[i] = cimag[i] + aii * Chebychev_polyn[4][i] + air * Chebychev_polyn[5][i];
            }
            for(i=0;i<dmatsize[0];i++){
                Chebychev_polyn[0][i] = Chebychev_polyn[2][i];
                Chebychev_polyn[1][i] = Chebychev_polyn[3][i];
                Chebychev_polyn[2][i] = Chebychev_polyn[4][i];
                Chebychev_polyn[3][i] = Chebychev_polyn[5][i];
                Chebychev_polyn[4][i] = 0;
                Chebychev_polyn[5][i] = 0;
            }
            update_polyn23();
        }

    }
}

void detector:: delete_variable_for_Chebychev_method(){
    int i ,j;
    for(i=0;i<stlnum;i++){
        delete [] shifted_total_dmat[i];
        delete [] shifted_dmat[i];
    }
    delete [] shifted_dmat;
    delete [] shifted_total_dmat;
    delete [] Bessel_function_array;
    delete [] Chebychev_polyn;
    for(i=0;i<6;i++){
        delete [] send_polyn[i];
        delete [] recv_polyn[i];
    }
    delete [] send_polyn;
    delete [] recv_polyn;

}