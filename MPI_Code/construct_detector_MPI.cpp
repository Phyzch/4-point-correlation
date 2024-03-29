//
// Created by phyzch on 6/23/20.
// This file contain function that use our previously coded programm to construct matrix.
//
#include"../system.h"
#include"../util.h"
using namespace std;
int distance_cutoff_for_4_piont_corre ;
double Energy_Range_4_point_corre_function_average = 0; // average over bunch of states within energy window for computation of 4 point correlation function
int Distance_Range_4_point_corre_function_average = 0;   // average over bunch of states within distance cutoff for computation of 4 point correlation function

void Broadcast_dmat_vmode(int stlnum, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);

// initialize parameters for detector and allocate space there.
void detector::allocate_space(int tlnum) {
    int i, j;
    dmatdim = detdim*detdim / fillfrac;

    nmodes = new int[tlnum];  //  number of modes in each detector

    proptime = new double[tlnum];  //  pre-coupling propagation time

    nmax = new int *[tlnum];  // maximum number of quanta in eqch mode.

    // if mode is bright or dark
    modtype = new int *[tlnum];

    mfreq = new double *[tlnum];  // mfreq: frequency of each mode here.

    aij = new double * [tlnum];

    modcoup = new double *[tlnum];  // coupling of matrix mod here.

    premodcoup = new double *[tlnum]; // coupling before system and detector contact with each other.

    dmat = new vector <double> [tlnum];

    for (i=0;i<tlnum;i++){
        vector<int> v1 ;
        dirow.push_back(v1);
        dicol.push_back(v1);
        vector<vector<int>> v2 ;
        dv.push_back(v2);
    }


    // matrix element number for detector matrix
    dmatnum = new int[tlnum];
    // off diagonal matrix element number for detector matrix
    doffnum = new int[tlnum];

    dmatsize = new int[tlnum];   // size of detector matrix

    // tell detector total matrix size..
    total_dmat_size.reserve(2);
    total_dmat_num.reserve(2);
    total_dmat_off_num.reserve(2);
    dmatsize_each_process = new int * [2];
    dmatsize_offset_each_process = new int * [2];
    doffnum_each_process= new int * [2];
    dmatnum_each_process= new int * [2];;  // record detector matrix element number in each process.
    dmat_offset_each_process= new int * [2];; // record local first detector matrix's index in global matrix.
    for(i=0;i<tlnum;i++){
        dmatsize_each_process[i]= new int [num_proc];
        dmatsize_offset_each_process[i] = new int [num_proc];
        doffnum_each_process[i]= new int [num_proc];
        dmatnum_each_process[i] = new int [num_proc];
        dmat_offset_each_process[i] = new int [num_proc];
    }

};

void detector:: allocate_space_single_detector(int detector_index){
    int i= detector_index;
    nmax[i] = new int[nmodes[i]];
    modtype[i] = new int[nmodes[i]];
    mfreq[i] = new double[nmodes[i]];
    aij[i] = new double[nmodes[i]];
    modcoup[i] = new double[nmodes[i]];
    premodcoup[i] = new double[nmodes[i]];
}

// read parameters for detector matrix construction.
void detector::read_MPI(ifstream & input, ofstream & output, ofstream & log, int tlnum, int tldim, string path) {
    int i, j;
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);

    a_intra= 0.2;
    stlnum=tlnum;
    stldim=tldim;
    allocate_space(tlnum);
    if (my_id==0){
        // matflag: mark used to see if we shift energy window.
        // maxdis: maximum allowed distance.
        // cutoff: cutoff strength for intra-detector coupling strength.  cutoff2: cutoff strength for inter-detector coupling.
        // kelvin: temperature used to initialize detector state.
        input >> matflag >> maxdis >> cutoff >> cutoff2 >> kelvin;
        if (! Detector_Continue_Simulation) {
            output << "Detector  " << matflag << " " << maxdis << " " << cutoff << " " << cutoff2 << " " << kelvin << endl;
        }
        for (i = 0; i < tlnum; i++) {
            input >> nmodes[i] >> proptime[i];
            allocate_space_single_detector(i);
            for (j = 0; j < nmodes[i]; j++) {
                input >> mfreq[i][j] >> nmax[i][j] >> modtype[i][j] >> premodcoup[i][j] >> modcoup[i][j];
            }
        }
        // -----------------------------------------------------------------------------------------------------
        // add noise to frequency.
        std::default_random_engine generator(time(NULL));
        std::uniform_real_distribution<double> distribution(0,1);  // random number in range [0,1]
        if (!Continue_Simulation) {
            for(i=0;i<tlnum;i++) {
                for(j=0;j<nmodes[0];j++) {
                    if (j == 0) {
                        mfreq[i][j] = mfreq[i][j] * (1 + noise_strength * (distribution(generator)-0.5) * 2 );
                        cout << mfreq[i][j] << " ";
                    }
                    else {
                        mfreq[i][j] = mfreq[i][j] * (1 + noise_strength * (distribution(generator)-0.5) * 2 );
                        cout << mfreq[i][j] << " ";
                    }
                    // set precision of mfreq to 0.01. Convenient to restart simulation.
                    mfreq[i][j] = floor(mfreq[i][j] *100) /100;
                }

                if (noise_strength not_eq 0){
                    cout << "Nonzero noise strength" << endl;
                }

            }
        }
        else{
            ifstream load;
            load.open(path+"save_wavefunction.txt");
            if (load.is_open()) {
                for (i = 0; i < tlnum; i++) {
                    for (j = 0; j < nmodes[0]; j++) {
                        load >> mfreq[i][j];
                    }
                }
                load.close();
            }
            else{
                log<<"Can not load Detector frequency from save_wave function.txt"<<endl;
                MPI_Abort(MPI_COMM_WORLD,-6); // error code -6 : fail to load detector frequency
            }
        }
        //-----------------------------------------------------------------------------------------------------
        for(i=0;i<tlnum;i++){
            if (!Detector_Continue_Simulation) {
                output << nmodes[i] << " " << proptime[i] << endl;
            }
            for(j=0;j<nmodes[i];j++){
                // aij[i][j] = f^{0.5} / 270
                aij[i][j] = pow(double(mfreq[i][j]),0.5) / 270 ;
                // aij corresponding to scaling factor for f= f_{bright}/2 cm^{-1}.
                if (! Detector_Continue_Simulation) {
                    output << mfreq[i][j] << " " << nmax[i][j] << " " << modtype[i][j] << " " << premodcoup[i][j]
                    << " " << modcoup[i][j] << endl;
                }
            }
        }
    }

    MPI_Bcast(&nmodes[0],tlnum,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&proptime[0],tlnum,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(my_id!=0){
        for(i=0;i<stlnum;i++){
            allocate_space_single_detector(i);
        }
    }

    // temporary value for calculation of coupling.
    deln = new int[max(nmodes[0],nmodes[1])];
    nbar = new double[max(nmodes[0],nmodes[1])];

    MPI_Bcast(&matflag,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&maxdis,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&cutoff,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&cutoff2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&kelvin,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // Boltzmann_beta used for compute regularized thermal OTOC is 1 / T. T given by kelvin
    double kelvin_in_wavenumber_unit = kelvin * 0.6945 ;
    Boltzmann_beta = 1 / kelvin_in_wavenumber_unit ;

    for(i=0;i<stlnum;i++){
        MPI_Bcast(&mfreq[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&nmax[i][0],nmodes[i],MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&modtype[i][0],nmodes[i],MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&premodcoup[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&modcoup[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
        MPI_Bcast(&aij[i][0],nmodes[i],MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    ndegre = nmodes[0];
    ndegrx2 = ndegre * 2;
};

void detector:: construct_dmatrix_MPI(ifstream & input, ofstream & output, ofstream & log,  vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1) {
    int m, i , j;
    construct_dv_dirow_dicol_dmatrix_MPI(log, dmat0, dmat1, vmode0, vmode1);
    compute_important_state_index();
    // -------------------------- Two different way of constructing off-diagonal term for detector  -----------------------------
    // 1:  traditional way of constructing detector off-diagonal part of matrix
    if(not read_Hamltonian_from_file){
        if(not no_coupling){
            compute_detector_offdiag_part_MPI(log,dmat0,dmat1,vmode0,vmode1);
        }
    }
    else{
//  2:  applying Van Vleck transformation:
        if (Turn_on_Vanvleck) {
            if (my_id == 0) {
                log << "Use Van Vleck transformation" << endl;
            }
        }
        if(!no_coupling){
            construct_state_coupling_vanvlk(dmat[0], dmat0, vmode0, dirow[0], dicol[0],output);
        }
    }


    //--------------------------------------------------------------------------------------------------
    update_initial_and_bright_detector_energy();

    output_state_density(dmat0,dmat1);

    broadcast_dmatnum_doffnum();
    broadcast_total_dmat();

    if(my_id==0){
        for(m=0;m<stlnum;m++){
            if (! Detector_Continue_Simulation) {
                output << "Matrix for detector " << m << " : " << total_dmat_size[m] << "  " << total_dmat_off_num[m] << endl;
                if (intra_detector_coupling) {
                    output << "Turn on function of rescale intra-detector coupling strength:  We scale coupling strength by:  " << intra_detector_coupling_noise << endl;
                }
                if (inter_detector_coupling ) {
                    output << "Turn on function of rescale inter-detector coupling strength:  We scale coupling strength between detector by:  " << inter_detector_coupling_noise << endl;
                }
            }
        }

    }
    // compute density of state
    compute_local_density_of_state(output,dmat0);

    //prepare variable for 4 piont correlation function

    prepare_variable_for_4_point_correlation_function(dmat0,dmat1,log);
}


void detector:: construct_dv_dirow_dicol_dmatrix_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1){
    int i,m;
    int vsize,vsize2;
    int displacement;
    if(my_id==0){
        total_dmat_size[0]= dmat0.size();
        if(stlnum==2){
            total_dmat_size[1] = dmat1.size();
        }
        else{
            total_dmat_size[1]=0;
        }
    }
    MPI_Bcast(&total_dmat_size[0],2,MPI_INT,0,MPI_COMM_WORLD);
    vector <int> * dirow_all, *dicol_all;
    vector< double> * dmat_all;
    int ** vector_size, ** displacement_list;

    for(i=0;i<2;i++){
        vector<vector<int>> v1 ;
        dv_all.push_back(v1);
    }
    dmat_all = new vector<double>[2];  // do not confuse dmat_all with total_dmat. dmat_all only contain diagonal term
    dirow_all = new vector<int>[2];
    dicol_all = new vector<int>[2];
    vector_size = new int *[2];
    displacement_list = new int *[2];

    Broadcast_dmat_vmode(stlnum, dmat0,dmat1,vmode0,vmode1); // broadcast dmat0, dmat1, vmode0, vmode1 to all process to compute off-diagonal matrix.
    dv_all[0] = vmode0;
    dv_all[1] = vmode1;

    if(my_id==0) {
        dmat_all[0] = dmat0;
        dmat_all[1] = dmat1;
        // prepare dirow, dicol for broadcasting.
        for (m = 0; m < stlnum; m++) {
            int size = dmat_all[m].size();
            for (i = 0; i < size; i++) {
                dirow_all[m].push_back(i);
                dicol_all[m].push_back(i);
            }
        }
        // prepare vector size and vector displacement:
        for (m = 0; m < stlnum; m++) {
            vsize = total_dmat_size[m] / num_proc;
            vsize2 = total_dmat_size[m] - (num_proc - 1) * vsize;
            vector_size[m] = new int[num_proc];  // size of vector to scatter to each process.
            displacement_list[m] = new int[num_proc]; // displacement of vector to scatter to each process.
            displacement = 0;
            for (i = 0; i < num_proc - 1; i++) {
                vector_size[m][i] = vsize;
                displacement_list[m][i] = displacement;
                displacement = displacement + vsize;
            }
            vector_size[m][num_proc - 1] = vsize2;
            displacement_list[m][num_proc - 1] = displacement;
        }
    }
    Scatter_dirow_dicol_dmatrix(dmat_all, dirow_all, dicol_all, vector_size, displacement_list, log);
    Scatter_dv(total_dmat_size);  // scatter dv_all
    // construct detector matrix size for each process.
    if(my_id != num_proc-1){
        for(m=0;m<stlnum;m++) {
            vsize= total_dmat_size[m]/num_proc;
            dmatsize[m] = vsize;
        }
    }
    else{
        for(m=0;m<stlnum;m++){
            vsize= total_dmat_size[m] - total_dmat_size[m]/num_proc *(num_proc-1);
            dmatsize[m]=vsize;
        }
    }
    for(m=0;m<stlnum;m++){
        MPI_Allgather(&dmatsize[m],1, MPI_INT,&dmatsize_each_process[m][0],1,MPI_INT,MPI_COMM_WORLD);
    }
    for(m=0;m<stlnum;m++){
        dmatsize_offset_each_process[m][0] = 0;
        for(i=1;i<num_proc;i++){
            dmatsize_offset_each_process[m][i] = dmatsize_offset_each_process[m][i-1] + dmatsize_each_process[m][i-1];
        }
    }

    if(my_id == 0) {
        for (m = 0; m < stlnum; m++) {
            delete[] vector_size[m];
            delete[] displacement_list[m];
        }
    }
    delete[] vector_size;
    delete[] displacement_list;
    delete[] dirow_all;
    delete[] dicol_all;
}

void Broadcast_dmat_vmode(int stlnum, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1){
    int i,j,m;
    int dmat0_size, dmat1_size;
    // we broadcast dmat0, dmat1, .. to all other process. This is need for us to compute off diagonal matrix
    // first allocate space for dmat0 , dmat1.
    if(my_id==0){
        dmat0_size= dmat0.size();
        dmat1_size= dmat1.size();
    }
    MPI_Bcast(&dmat0_size,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&dmat1_size,1,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id!=0) {
        dmat0.resize(dmat0_size);
        dmat1.resize(dmat1_size);
    }
    // Broadcast dmat0, dmat1
    MPI_Bcast((void *) dmat0.data(),dmat0_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast((void *) dmat1.data(),dmat1_size,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // ----------------Broadcast vmode0, vmode1 -------------------------------
    vector<vector<int>> * v_ptr;
    // ------------- For sending vmode0,vmode1 --------------------------
    vector<int> vmode_1d;
    vector<int> displacement;
    vector <int> element_size;
    int vmode_1d_size;
    int element_number;
//-------------------------------------------
    for (m = 0; m < stlnum; m++) {
        vmode_1d.clear();
        displacement.clear();
        element_size.clear();
        if (m == 0) v_ptr = &(vmode0);
        else v_ptr = &(vmode1);
        if(my_id==0) {
            convert_dv(*v_ptr, vmode_1d, displacement, element_size);
            vmode_1d_size = vmode_1d.size();
            MPI_Bcast(&vmode_1d_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
            element_number = element_size.size();
            MPI_Bcast(&element_number, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast((void *) element_size.data(), element_number, MPI_INT, 0,
                      MPI_COMM_WORLD); // send size of each 2d element
            MPI_Bcast((void *) vmode_1d.data(), vmode_1d_size, MPI_INT, 0,
                      MPI_COMM_WORLD); // send converted 1d data to other process.
        }
        else{
            MPI_Bcast(&vmode_1d_size,1,MPI_INT,0,MPI_COMM_WORLD);
            vmode_1d.reserve(vmode_1d_size);
            MPI_Bcast(&element_number,1,MPI_INT,0,MPI_COMM_WORLD);
            element_size.reserve(element_number);
            MPI_Bcast( (void *) element_size.data(), element_number,MPI_INT,0,MPI_COMM_WORLD);
            MPI_Bcast((void *) vmode_1d.data(),vmode_1d_size,MPI_INT,0,MPI_COMM_WORLD);

            int index=0;
            // convert vmode_1d to vmode0 / vmode1
            for(i=0;i<element_number;i++){
                vector<int> dv_vmode;
                dv_vmode.reserve(element_size[i]);
                for(j=0;j<element_size[i];j++){
                    dv_vmode.push_back(vmode_1d[index]);
                    index++;
                }
                (*v_ptr).push_back(dv_vmode);
            }
        }
    }
}




void detector::compute_detector_offdiag_part_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1){
    int i,j,m,k;
    int begin_index;
    int ntot;
    double value, lij;
    vector<vector <int>> * vmode_ptr;
    vector<double> * dmat_ptr;
    bool exist;
    int position;
    double random_number;
    // different process do different amount of work.
    for(m=0;m<stlnum;m++){
        if(m==0){
            vmode_ptr = &(vmode0);
            dmat_ptr= &(dmat0);
        }
        else {
            vmode_ptr= &(vmode1);
            dmat_ptr= &(dmat1);
        }
        begin_index= total_dmat_size[m]/num_proc * my_id;
        // compute off diagonal matrix element, using dmat0, dmat1.
        for(i=begin_index;i<begin_index + dmatsize[m];i++){  // for my_id==0 , O(dmatsize * dmatsize/ proc_num)
            for(j=0;j<total_dmat_size[m];j++){ // j is different from serial program. Now we record both symmetric Hamiltonian element
                if (i==j) continue;
                ntot=0;
                for(k=0;k<nmodes[m];k++){
                    deln[k]= abs( (*vmode_ptr)[i][k] - (*vmode_ptr)[j][k] ); // same as deln[k] = abs(dv[m][i][k] - dv[m][j][k]);
                    nbar[k]= sqrt(sqrt(double(max(1, (*vmode_ptr)[i][k])) * double(max(1, (*vmode_ptr)[j][k]  ))));
                    ntot= ntot+ deln[k];
                }
                // this is because we don't have 1st order ladder operator in Harmonic oscillator's expression
                if (ntot == 2) {
                    for (k = 0; k < nmodes[m]; k++) {
                        if (deln[k] == 2) deln[k] = 4;
                        if (deln[k] == 1) deln[k] = 2;
                    }
                } else if (ntot == 1) {
                    for (k = 0; k < nmodes[m]; k++) {
                        if (deln[k] == 1) deln[k] = 3;
                    }
                } else if (ntot == 0) {
                    log << "Error! The off-diagonal element must differ in q.n." << endl;
                    MPI_Abort(MPI_COMM_WORLD,-8);
                }
                // check ntot with maxdis in quantum number space:
                if (ntot < maxdis) {

                    if (ntot % 2 == 0) {
                        value = V_intra;  // V_intra == 3050 as in model
                    } else {
                        value = -V_intra;
                    }

                    // noise term
                    if (intra_detector_coupling) {
                        do (random_number = 2*((double(rand())/RAND_MAX)-0.5)  ); while (random_number==0) ;
                        value = value * (1+intra_detector_coupling_noise * random_number);
                    }

                    for (k = 0; k < nmodes[m]; k++) {
                        // aij  = f^{1/2} / 270
                        value = value * pow(aij[m][k]* nbar[k], deln[k]);
                    }
                    if ( (*dmat_ptr)[i] != (*dmat_ptr)[j] ) {
                        lij = abs(value / ((*dmat_ptr)[i] - (*dmat_ptr)[j]));
                        if (lij > cutoff) {
                            dmat[m].push_back(value);
                            dirow[m].push_back(i);
                            dicol[m].push_back(j);
                        }
                    } else {
                        dmat[m].push_back(value);
                        dirow[m].push_back(i);
                        dicol[m].push_back(j);
                    }
                }
            }
        }
    }
}

void detector:: broadcast_dmatnum_doffnum(){
    int m,i;
    // compute dmatnum and doffnum and dmatnum_each_process, total_dmatnum, total_doffnum etc.
    for(m=0;m<stlnum;m++){
        dmatnum[m]= dmat[m].size();
        doffnum[m]=dmatnum[m] - dmatsize[m];
    }
    // compute toatl_dmatnum, total_dmatoffnum
    for(m=0;m<stlnum;m++){
        MPI_Allreduce(&dmatnum[m],&total_dmat_num[m],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(&doffnum[m],&total_dmat_off_num[m],1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    }
    // broadcast detector matrix number and detector matrix off-diagonal number in each process.
    for(m=0;m<stlnum;m++){
        MPI_Allgather(&dmatnum[m],1,MPI_INT,&dmatnum_each_process[m][0],1,MPI_INT,MPI_COMM_WORLD);
        MPI_Allgather(&doffnum[m], 1, MPI_INT, &doffnum_each_process[m][0], 1, MPI_INT,MPI_COMM_WORLD);
    }
    // compute offset of detector matrix for each process.
    for(m=0;m<stlnum;m++){
        dmat_offset_each_process[m][0]=0;
        for(i=1;i<num_proc;i++){
            dmat_offset_each_process[m][i] = dmat_offset_each_process[m][i-1] + dmatnum_each_process[m][i-1];
        }
    }
}


void detector:: broadcast_total_dmat(){
    /*
     construct total_dmat, total_dirow, total_dicol
     use dmatnum_each_process,  dmat_offset_each_process\
    */
    total_dmat= new double * [stlnum];
    total_dirow= new int * [stlnum];
    total_dicol= new int * [stlnum];
    total_dmat_diagonal = new double * [stlnum];
    int m;
    for(m=0;m<stlnum;m++){
        total_dmat_diagonal[m] = new double [total_dmat_size[m]];
        total_dmat[m] = new double [total_dmat_num[m]];
        total_dirow[m] = new int [total_dmat_num[m]];
        total_dicol[m] = new int [total_dmat_num[m]];
        MPI_Allgatherv(&dmat[m][0], dmatsize[m] , MPI_DOUBLE,
                       &total_dmat_diagonal[m][0], dmatsize_each_process[m], dmatsize_offset_each_process[m], MPI_DOUBLE, MPI_COMM_WORLD);

        MPI_Allgatherv(&dmat[m][0],dmatnum[m],MPI_DOUBLE,
                &total_dmat[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_DOUBLE,MPI_COMM_WORLD);
        MPI_Allgatherv(&dirow[m][0],dmatnum[m],MPI_INT,
                &total_dirow[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_INT,MPI_COMM_WORLD);
        MPI_Allgatherv(&dicol[m][0],dmatnum[m],MPI_INT,
                &total_dicol[m][0],dmatnum_each_process[m],dmat_offset_each_process[m],MPI_INT,MPI_COMM_WORLD);
    }
}

//-------------------------------- Construct detector wave function ------------------------------
void detector::initialize_detector_state_MPI(ofstream & log) {
    // as we have to compute Boltzmann weighted wave function here. make sure prepare_evolution function is called before this function to prepare for update_dx, update_dy here.
    // also make sure prepare Boltzmann computation and prepare ladder operator computation is called before this function.
    // this is the code for initializing detector state's wavefunction.
    int m, i,j , k;
    double norm;
    double total_norm;
    int local_index = 0;
    int dmat0_offset = my_id * int(total_dmat_size[0]/num_proc);
    // initialize basis set state |{n}> to compute |{n}(t)>
    int nearby_state_list_size = nearby_state_index.size();
    for(m=0;m<nearby_state_list_size;m++){
        norm = 0;
        for(i=0;i<dmatsize[0];i++){
            xd[m].push_back(0);
            yd[m].push_back(0);
        }
        if(nearby_state_index[m] >= dmat0_offset and nearby_state_index[m] < dmat0_offset + dmatsize[0]){
            // state m is in this range
            local_index = nearby_state_index[m] - dmat0_offset;
            xd[m][local_index] = 1;
            norm = 1;
        }
        MPI_Allreduce(&norm,&total_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        if(total_norm==0){
            if(my_id==0) {
                cout << "Norm for detector state "<< m <<" is 0" << endl;
                log << "Norm for detector state "<<m<<" is 0" << endl;
                MPI_Abort(MPI_COMM_WORLD,-10);
            }
        }
    }

    // initialize Haar random state . and then decorated it with Boltzmann factor.
    std::default_random_engine  generator;
    std::normal_distribution<double> distribution(0,1);
    int state_index ;

    for(i=0;i<N_Haar; i++ ){
        state_index = regularized_Haar_state_index_list[i];

        // initialize Haar random state  xd [state_index] , yd[state_index]
        // draw distribution from Gaussian random distribution. Then normalize it.
        norm = 0;
        for(k=0;k<dmatsize[0];k++){
            xd[state_index][k] = distribution(generator);
            yd[state_index][k] = distribution(generator);

            norm = norm + std::norm(xd[state_index][k]) + std::norm(yd[state_index][k]) ;
        }

        MPI_Allreduce(&norm, &total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        total_norm = sqrt(total_norm);
        for(k=0;k<dmatsize[0];k++){
            xd[state_index][k] = xd[state_index][k] / total_norm;
            yd[state_index][k] = yd[state_index][k] / total_norm;
        }

    }


    // decorate Haar random state with Boltzmann factor e^{-\beta H/4}
    update_dx(state_number_for_evolution);
    update_dy(state_number_for_evolution);

    for(i = 0; i< N_Haar; i++ ){
        state_index = regularized_Haar_state_index_list[i] ;

        vector<double> Boltzmann_factor_weighted_wave_func_x ;
        vector<double> Boltzmann_factor_weighted_wave_func_y;
        Chebyshev_method_Boltzmann_factor( xd[state_index] , yd[state_index],
                                          Boltzmann_factor_weighted_wave_func_x, Boltzmann_factor_weighted_wave_func_y);
        xd[ state_index ] = Boltzmann_factor_weighted_wave_func_x ;
        yd[ state_index ] = Boltzmann_factor_weighted_wave_func_y ;

    }


    // compute wave function after ladder operator a_{j} is operated on Boltzmann weighted wave function:
    for( i = 0; i<N_Haar; i++ ){
        state_index = regularized_Haar_state_index_list[i] ;  // index for Haar random variable (after multiply e^{-\beta H/4})
        vector<vector<double>>  xd_for_ladder_operator;
        vector<vector<double>>  yd_for_ladder_operator;
        ladder_operator_operation(xd[state_index] , yd[state_index], xd_for_ladder_operator, yd_for_ladder_operator );
        for(j=0;j<2*nmodes[0];j++){
            state_index = state_index + 1;  // state_index for a_{j} e^{-beta H/4} | \phi>
            for(k=0;k<dmatsize[0];k++){
                xd[state_index][k] = xd_for_ladder_operator[j][k];
                yd[state_index][k] = yd_for_ladder_operator[j][k];
            }
        }
    }


    update_dx(state_number_for_evolution);
    update_dy(state_number_for_evolution);

    // for comparision, we use same Haar random state to compute unregularized Lyapunov spectrum.
    // we initialize it as e^{-\beta H /2 } |Haar>
    int unregularized_Haar_state_index;
    int regularized_Haar_state_index;
    for(i=0; i<N_Haar; i++ ){
        unregularized_Haar_state_index = unregularized_Haar_state_index_list[i];
        regularized_Haar_state_index = regularized_Haar_state_index_list[i];

        xd[unregularized_Haar_state_index] = xd[regularized_Haar_state_index];
        yd[unregularized_Haar_state_index] = yd[regularized_Haar_state_index];

        vector<double> Boltzmann_factor_weighted_wave_func_x ;
        vector<double> Boltzmann_factor_weighted_wave_func_y;

        Chebyshev_method_Boltzmann_factor(xd[unregularized_Haar_state_index] , yd[unregularized_Haar_state_index] ,
                                          Boltzmann_factor_weighted_wave_func_x, Boltzmann_factor_weighted_wave_func_y) ;
        // e^{-\beta H /2 } |Haar>
        xd[unregularized_Haar_state_index] = Boltzmann_factor_weighted_wave_func_x;
        yd[unregularized_Haar_state_index] = Boltzmann_factor_weighted_wave_func_y;

    }

    // apply ladder operator to e^{-\beta H /2 } |Haar>
    for(i=0;i<N_Haar; i++ ){
        unregularized_Haar_state_index = unregularized_Haar_state_index_list[i] ;  // index for Haar random variable (after multiply e^{-\beta H/4})
        vector<vector<double>>  xd_for_ladder_operator;
        vector<vector<double>>  yd_for_ladder_operator;
        ladder_operator_operation(xd[unregularized_Haar_state_index] , yd[unregularized_Haar_state_index], xd_for_ladder_operator, yd_for_ladder_operator );
        for(j=0;j<2*nmodes[0];j++){
            unregularized_Haar_state_index = unregularized_Haar_state_index + 1;  // state_index for a_{j} e^{-beta H/4} | \phi>
            for(k=0;k<dmatsize[0];k++){
                xd[unregularized_Haar_state_index][k] = xd_for_ladder_operator[j][k];
                yd[unregularized_Haar_state_index][k] = yd_for_ladder_operator[j][k];
            }
        }
    }

}

void detector::construct_bright_state_MPI(ifstream & input, ofstream & output){
    // MPI version of construct bright state
    /*
     *  Initial detector state: detector state populated at beginning of simulation
     *  Bright  detector state: detector state with bright mode set to 1 ,all other mode equal to initial detector state
     *  Initial_Detector_energy: Energy of detector at initial state.
     */
    int m,i,j;
    double norm;
    double random_number;
    double prob;
    bright_state=new int * [stlnum];
    initial_detector_state= new int * [stlnum];
    initial_Detector_energy= new double [stlnum];
    bright_state_energy= new double [stlnum];
    for(m=0;m<stlnum;m++){
        bright_state[m]= new int [nmodes[m]];
        initial_detector_state[m]= new int [nmodes[m]];
        initial_Detector_energy[m]=0;
        bright_state_energy[m]=0;
    }
    if(my_id==0){
        for(m=0;m<stlnum;m++) {
            // initialize our bright state from input.txt
            if (!Random_bright_state) {
                for (i = 0; i < nmodes[m]; i++) {
                    input >> initial_detector_state[m][i];
                }
            }
            else {
                // construct bright state according to Boltzmann distribution.
                for (i = 0; i < nmodes[m]; i++) {
                    norm = 0;
                    for (j = 0; j < nmax[m][i]; j++) {
                        norm = norm + exp(-mfreq[m][i] * j / (0.7 * kelvin)); // 300 K ~ 200 cm^{-1}
                    }
                    do (random_number = double(rand()) / RAND_MAX); while (random_number == 0);
                    prob = 0;
                    for (j = 0; j < nmax[m][i]; j++) {
                        prob = prob + exp(-mfreq[m][i] * j / (0.7 * kelvin)) / norm;
                        if (prob > random_number) {
                            initial_detector_state[m][i] = j;
                            break;
                        }
                    }
                }
            }
        }
    }
    for(m=0;m<stlnum;m++){ // Broad cast initial detector state.
        MPI_Bcast(&initial_detector_state[m][0],nmodes[m],MPI_INT,0,MPI_COMM_WORLD);
    }

    for(m=0;m<stlnum;m++){
        // initialize our initial detector state.  set dark mode's quanta equal to bright_state.
        bright_state[m][0]= initial_detector_state[m][0] + 1;
        for(i=1;i<nmodes[m];i++){   // other dark mode
            bright_state[m][i]=initial_detector_state[m][i];
        }
        for(i=0;i<nmodes[m];i++){
            initial_Detector_energy[m]= initial_Detector_energy[m] + initial_detector_state[m][i] * mfreq[m][i];
            bright_state_energy[m] = bright_state_energy[m] + bright_state[m][i] * mfreq[m][i];
        }
    }
    initial_energy= initial_energy + initial_Detector_energy[0] + initial_Detector_energy[1]; // update initial energy
    // as system + detector.

    if(my_id==0){  // output initial detector state to output.txt
        cout <<"Initial detector state:"<<endl;
        for(m=0;m<stlnum;m++){
            for(i=0;i<nmodes[m];i++){
                cout<< initial_detector_state[m][i] <<" ";
            }
            cout<<endl;
        }
    }
}

void detector:: update_initial_and_bright_detector_energy(){
    // we update energy of initial and bright state of detector. since in Van Vleck transformation, the energy level is shhifted.
    int m,i;
    for(m=0;m<stlnum;m++){
        initial_Detector_energy[m] = 0;
        bright_state_energy[m] = 0;
        if(my_id == initial_state_pc_id[m]) {
            initial_Detector_energy[m] = dmat[m][initial_state_index[m]];
        }
        if(my_id == bright_state_pc_id[m] and bright_state_index[m]<dmatsize[m]){
            bright_state_energy[m] = dmat[m][bright_state_index[m]];
        }
        MPI_Bcast(&initial_Detector_energy[m],1,MPI_DOUBLE,initial_state_pc_id[m],MPI_COMM_WORLD);
        MPI_Bcast(&bright_state_energy[m],1,MPI_DOUBLE,bright_state_index[m],MPI_COMM_WORLD);
    }
    initial_energy= system_energy + initial_Detector_energy[0] + initial_Detector_energy[1]; // update initial energy as system + detector.
}
void detector:: compute_important_state_index(){
    // compute bright state index and initial state index for two detector.
    int m,i,j;
    int index;
    bool check_mark1, check_mark2;
    vector<int> special_state_vec ;
    int ** special_state;
    int * special_state_index;
    int * special_state_pc_id;
    int position;
    bool exist;
    bool * exist_bool_for_pc = new bool [num_proc];
    bright_state_index = new int [2];
    initial_state_index = new int [2];
    bright_state_pc_id = new int [2];
    initial_state_pc_id = new int [2];
    // we initialize bright_state_index and initial_state_index.
    // We do this because sometimes we may not include bright state in our simulation, then this index need to be initialized.
    for(m=0;m<stlnum;m++){
        bright_state_index[m] = 0;
        initial_state_index[m]=0;
        bright_state_pc_id[m] = 0;
        initial_state_pc_id[m] = 0;
    }
    for(index =0 ; index <1; index ++) {    // when you want to consider bright state, set index<=1
        // loop for bright state and initial state
        if(index == 0) {
            special_state = initial_detector_state;
            special_state_index = initial_state_index;
            special_state_pc_id = initial_state_pc_id;
        }
        else{
            special_state = bright_state;
            special_state_index = bright_state_index;
            special_state_pc_id = bright_state_pc_id;
        }
        // loop for two detector.
        for (m = 0; m < stlnum; m++) {
            special_state_vec.clear();
            for (i = 0; i < nmodes[m]; i++) {
                special_state_vec.push_back(special_state[m][i]);
            }
            position=find_position_for_insert_binary(dv[m],special_state_vec,exist);
            MPI_Allgather(&exist,1,MPI_C_BOOL,&exist_bool_for_pc[0],1,MPI_C_BOOL,MPI_COMM_WORLD);
            special_state_pc_id [m] = -1;
            for(i=0;i<num_proc;i++){
                if(exist_bool_for_pc[i]){
                    special_state_pc_id[m] = i;
                }
            }
            if(special_state_pc_id[m] == -1){
                if(my_id==0){
                    cout<<" Can not find initial state or brigth state in all vmode. Must have bug here."<<endl;
                    MPI_Abort(MPI_COMM_WORLD,-7);
                }
            }
            MPI_Bcast(&position,1, MPI_INT,special_state_pc_id[m],MPI_COMM_WORLD);
            special_state_index [m] = position;
        }
    }
    delete [] exist_bool_for_pc;
}


void compute_detector_state_density(vector <int> & state_density, vector <double> & energy_distribution,
                                    vector<double> & energy_range){
    /*
     input: state_density: D_{l}(E) we are going to calculate
            energy_distribution: energy of different state in energy window
            energy_range: record range of energy.
    */
    int i;
    double min_energy= *min_element(energy_distribution.begin(),energy_distribution.end());
    double max_energy= *max_element(energy_distribution.begin(),energy_distribution.end());
    int block_number= 20;  // use to specify number of energy blocks we want to consider.
    double energy_step = (max_energy- min_energy)/block_number;
    int state_number = energy_distribution.size();
    int index;
    state_density.assign(block_number,0);
    for(i=0;i<state_number;i++){
        index = (energy_distribution[i] - min_energy) / energy_step;
        if(index == block_number ) index= index -1;  // in case we meet max_energy element
        state_density[index]++;
    }
    double energy=min_energy;
    for(i=0;i<block_number+1;i++){
        energy_range.push_back(energy);
        energy= energy + energy_step;
    }
}

void detector:: output_state_density(vector<double> & dmat0,  vector<double> & dmat1){
    int i;
    // compute_detector state_density
    vector<int> state_density0;
    vector <int> state_density1;
    vector <double> energy_range0;
    vector<double> energy_range1;
    int block_number;

    vector<double> dmat_energy_level0;
    vector<double> dmat_energy_level1;
    if(my_id == 0) {
        ofstream state_density_output(path + "detector_state_density");
        for (i = 0; i < total_dmat_size[0]; i++) {
            dmat_energy_level0.push_back(dmat0[i]);  // dmat0 is diagonal part in all matrix.
        }
        for (i = 0; i < total_dmat_size[1]; i++) {
            dmat_energy_level1.push_back(dmat1[i]);
        }
        compute_detector_state_density(state_density0, dmat_energy_level0, energy_range0);
        compute_detector_state_density(state_density1, dmat_energy_level1, energy_range1);
        block_number = state_density0.size();
        state_density_output << "Detector 1 Range of energy block: " << endl;
        for (i = 0; i <= block_number; i++) {
            state_density_output << energy_range0[i] << " ";
        }
        state_density_output << endl;
        state_density_output << "Detector 1 density of state " << endl;
        for (i = 0; i < block_number; i++) {
            state_density_output << state_density0[i] << " ";
        }
        state_density_output << endl;
        state_density_output << " Detector 2 Range of energy block" << endl;
        for (i = 0; i <= block_number; i++) {
            state_density_output << energy_range1[i] << " ";
        }
        state_density_output << endl;
        state_density_output << "Detector 2 density of state " << endl;
        for (i = 0; i < block_number; i++) {
            state_density_output << state_density1[i] << " ";
        }
        state_density_output << endl;
        // output initial state's energy
        for (i = 0; i < 2; i++) {
            state_density_output << initial_Detector_energy[i] << "  " << bright_state_energy[i] << " ";
        }
        state_density_output.close();
    }
}

void detector:: prepare_variable_for_4_point_correlation_function(vector<double> & dmat0, vector<double> & dmat1,ofstream & log){
    // for 4 - point correlation function we will make xd , yd as N*N system.
    // we will only include state near our initial state.
    int i,j,k,l;
    int Haar_state_index;
    int max_state_quanta_diff = 0;
    int state_distance = 0;
    double maximum_freq = 0;
    double energy_cutoff_for_nearby_state_index;

    double state_energy_difference;
    int initial_state_index_in_total_dmatrix;
    int nearby_state_index_size;
    int state_for_4_point_correlation_average_list_size ;
    initial_state_index_in_total_dmatrix = initial_state_index[0]
                                           + total_dmat_size[0]/num_proc * initial_state_pc_id[0] ;

    for(i=0;i< nmodes[0]; i++){
        if(mfreq[0][i] > maximum_freq){
            maximum_freq = mfreq[0][i];
        }
    }
    energy_cutoff_for_nearby_state_index = 4 / Boltzmann_beta + maximum_freq * 2 ;
    if(my_id == 0){
        cout << "energy cutoff for nearby state index " << energy_cutoff_for_nearby_state_index << endl;
    }

    // compute nearby_state_index and initial_state_index_in_nearby_state_index_list for computing 4-point correlation function
    for(i=0;i<total_dmat_size[0];i++){
        if( total_dmat_diagonal[0][i] < energy_cutoff_for_nearby_state_index ){
            nearby_state_index.push_back(i);
        }
    }

    // incorporate random Haar state to compute thermal average. For information for random Haar state, see: https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.112.120601
    state_number_for_evolution = nearby_state_index.size();
    Haar_state_index = state_number_for_evolution ;
    for(i=0;i<N_Haar; i++){
        regularized_Haar_state_index_list.push_back(Haar_state_index);
        Haar_state_index = Haar_state_index + 2 * nmodes[0] + 1 ; // we have to incorporate a_{j} e^{-\beta H/4} | Haar>
    }
    state_number_for_evolution = state_number_for_evolution + ( 2 * nmodes[0] + 1 ) * N_Haar;

    // for computing unregularized Lyapunov spectrum
    Haar_state_index = state_number_for_evolution;
    for(i=0; i< N_Haar; i++ ){
        unregularized_Haar_state_index_list.push_back(Haar_state_index);
        Haar_state_index = Haar_state_index + 2 * nmodes[0] + 1 ;
    }
    state_number_for_evolution = state_number_for_evolution + (2*nmodes[0] + 1 ) * N_Haar;

    nearby_state_index_size = nearby_state_index.size();
    for(i=0;i<nearby_state_index_size;i++){
        if(nearby_state_index[i] == initial_state_index_in_total_dmatrix){
            initial_state_index_in_nearby_state_index_list = i;
            break;
        }
    }


    log<< "Nearby state number:   "<< nearby_state_index_size <<endl;
    log<< "Total state number:    "<< total_dmat_size[0] <<endl;
    log<< "Total state for evolution " << state_number_for_evolution << endl;

    for (i = 0; i < state_number_for_evolution ; i++) {
        vector<double> v1 ;
        v1.reserve(dmatsize[0]);
        xd.push_back(v1);
        yd.push_back(v1);
    }
    xd_all = new double * [state_number_for_evolution];
    yd_all = new double * [state_number_for_evolution];
    for(i=0;i<state_number_for_evolution;i++){
        xd_all[i] = new double [total_dmat_size[0]];
        yd_all[i] = new double [total_dmat_size[0]];
    }


    // construct neighbor_state_index : These are states with one quantum number difference in state space. a_{j} |{n}>
    vector<int> neighbor_state_index;  // index in all dv_all[0]
    vector<int> neighbor_state_in_nearby_state_index; // index in nearby_state_index
    vector<int> state_mode;
    vector<int> neighbor_state_mode ;
    int position;
    bool exist;


    for(j=0;j<nearby_state_index_size;j++){
        state_mode = dv_all[0][ nearby_state_index[j] ];  // states we simulate trajectories.
        neighbor_state_index.clear();  // record states near to our states of choice for simulating trajectory.
        for(k=0;k<nmodes[0];k++){
            // one mode up
            neighbor_state_mode = state_mode;
            neighbor_state_mode[k] = state_mode[k] - 1;
            position = find_position_for_insert_binary(dv_all[0],neighbor_state_mode,exist);
            if(exist){
                neighbor_state_index.push_back(position);
            }
            else{
                neighbor_state_index.push_back(-1);
            }
        }

        for(k=0;k<nmodes[0];k++){
            neighbor_state_mode = state_mode;
            neighbor_state_mode[k] = state_mode[k] + 1;
            position = find_position_for_insert_binary(dv_all[0],neighbor_state_mode,exist);
            if(exist){
                neighbor_state_index.push_back(position);
            }
            else{
                neighbor_state_index.push_back(-1);
            }
        }
        neighbor_state_index_list.push_back(neighbor_state_index);

        neighbor_state_in_nearby_state_index.clear();
        // compute neighbor_state_in_nearby_state_index
        for(k=0;k<2*nmodes[0];k++){

            if(neighbor_state_index[k] != -1) {
                for (l = 0; l < nearby_state_index_size; l++) {
                    if (nearby_state_index[l] == neighbor_state_index[k]) {
                        neighbor_state_in_nearby_state_index.push_back(l);
                        break;
                    }
                }
                if (l == nearby_state_index_size) {
                    // not found
                    neighbor_state_in_nearby_state_index.push_back(-1);
                }
            }

            else {
                neighbor_state_in_nearby_state_index.push_back(-1);
            }

        }
        neighbor_state_in_nearby_state_index_list.push_back(neighbor_state_in_nearby_state_index);

        // compute bool_state_one_mode_quanta_below_all_in_nearby_state_index
        exist = true;
        for(k = 0 ;k< nmodes[0];k++){
            if (neighbor_state_in_nearby_state_index[k] == -1){
                exist = false;
            }
        }

        if(exist == true ){
            bool_state_one_mode_quanta_below_all_in_nearby_state_index.push_back(true);
        }
        else{
            bool_state_one_mode_quanta_below_all_in_nearby_state_index.push_back(false) ;
        }

        // compute bool_neighbor_state_all_in_nearby_state_index
        exist = true;
        for(k=0;k<nmodes[0];k++){
            if( neighbor_state_in_nearby_state_index[k] == -1 and state_mode[k]!=0 ){
                exist = false;
            }
        }
        for(k=nmodes[0];k<2*nmodes[0];k++){
            if (neighbor_state_in_nearby_state_index[k] == -1){
                exist = false;
            }
        }
        if(exist == false){
            bool_neighbor_state_all_in_nearby_state_index.push_back(false);
        }
        else{
            bool_neighbor_state_all_in_nearby_state_index.push_back(true);
        }

    }
    // construct_nearby_state_index_for_all_states
    for(j=0;j<total_dmat_size[0];j++){
        state_mode = dv_all[0][j];
        neighbor_state_index.clear();
        for(k=0;k<nmodes[0];k++){
            neighbor_state_mode = state_mode;
            neighbor_state_mode[k] = state_mode[k] - 1 ;
            position = find_position_for_insert_binary(dv_all[0],neighbor_state_mode,exist);
            if(exist){
                neighbor_state_index.push_back(position);
            }
            else{
                neighbor_state_index.push_back(-1);
            }
        }
        for(k=0;k<nmodes[0];k++){
            neighbor_state_mode = state_mode;
            neighbor_state_mode[k] = state_mode[k] + 1 ;
            position = find_position_for_insert_binary(dv_all[0],neighbor_state_mode,exist);
            if(exist){
                neighbor_state_index.push_back(position);
            }
            else{
                neighbor_state_index.push_back(-1);
            }
        }
        neighbor_state_index_for_all_state_list.push_back(neighbor_state_index);
    }

}

void detector:: compute_local_density_of_state(ofstream & output,vector<double> & dmat0){
    // using eq.(2) in https://doi.org/10.1063/1.476070 to compute density of states: sum Lij
    int i,l;
    double energy_difference;
    int begin_index = total_dmat_size[0]/num_proc * my_id ;
    double density_of_states = 0;
    int number_of_state= 0 ;
    if(my_id == initial_state_pc_id[0]){
        for(i=dmatsize[0];i<dmatnum[0];i++){
            if(dirow[0][i] == begin_index + initial_state_index[0]){
                number_of_state = number_of_state + 1;
                energy_difference = abs(dmat[0][initial_state_index[0]] - dmat0[dicol[0][i]]);
                density_of_states = density_of_states + 1 / (1+ pow(  energy_difference/dmat[0][i] , 2 ) );
            }
        }
    }
    MPI_Bcast(&density_of_states,1,MPI_DOUBLE,initial_state_pc_id[0],MPI_COMM_WORLD);
    MPI_Bcast(&number_of_state,1,MPI_INT,initial_state_pc_id[0],MPI_COMM_WORLD);
    if(my_id == 0){
        output<< "Density of states for detector 0 at initial state:   " << density_of_states << endl;
        output<< "number of states connected to it is   " << number_of_state << endl;
    }


    double  max_coupling_strength_in_all_process;
    double max_coupling_strength = 0;
    int num_terms=0;
    int num_terms_in_all_process=0;
    double  coupling_strength_sum ;
    double  coupling_strength_sum_all_process ;
    double  coupling_strength_average_all_process;
    double ratio_to_select_for_average = double(1)/5; // select value larger than 1/10 of maximum coupling strength

    for(i=dmatsize[0];i<dmatnum[0];i++){
        if(abs(dmat[0][i]) > max_coupling_strength){
            max_coupling_strength = abs(dmat[0][i]);
        }
    }
    MPI_Allreduce(&max_coupling_strength,&max_coupling_strength_in_all_process,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    coupling_strength_sum = 0;
    for(i=dmatsize[0];i<dmatnum[0];i++){
        if(abs(dmat[0][i]) > max_coupling_strength_in_all_process * ratio_to_select_for_average ){
            coupling_strength_sum = coupling_strength_sum + abs(dmat[0][i]);
            num_terms ++ ;
        }
    }
    MPI_Allreduce(&coupling_strength_sum,&coupling_strength_sum_all_process,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&num_terms,&num_terms_in_all_process,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    coupling_strength_average_all_process = coupling_strength_sum_all_process / (num_terms_in_all_process);
    if(my_id == 0){
        output <<"Largest coupling strength is "<<max_coupling_strength_in_all_process << endl;
        output <<"ratio for counting as states is "<< ratio_to_select_for_average << endl;
        output<< "number of terms for coupling strength average is: " << num_terms_in_all_process << endl;
        output << "Avereage coupling strength V is:   "<<coupling_strength_average_all_process << endl;
    }
}