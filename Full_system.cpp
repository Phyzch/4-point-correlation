#include"system.h"
#include"util.h"
using namespace std;
// advise: When you get lost in all these functions. Keep eye on fullsystem::fullsystem() and full_system::Quantum_evolution. Because
//they are functions do most of works and call other functions here.


//choice to turn on intra-detector coupling and scale it
bool intra_detector_coupling;  // mark if to scale intra-detector coupling strength.
double intra_detector_coupling_noise;
// We set this value to explore the coarse-graining of our system  i.e. Compare system with different dof

// choice to turn on inter_detector coupling and scale it.
bool inter_detector_coupling;
double inter_detector_coupling_noise;

// choice to continue our simulation: read wavefunction in  save_wavefunction.txt
bool Continue_Simulation;

// choice to only simulate detector
bool detector_only;

double noise_strength;

//choice to continue simulation of detector precoupling state: (Used to investigate decoherence of detector)
bool Detector_Continue_Simulation;

// energy window choice: We usually turn this on for large dof simulation. This will dramatically decrease computational cost for simulation
bool energy_window ;
double energy_window_size;  // size of energy window, we will only include whole system state whose energy difference with

double initial_energy;  // energy of system + detector.
double system_energy;  // energy of photon
bool Random_bright_state;

// initialization of parameters and do some pre-coupling set up
full_system::full_system(string path1, string cvpt_path1) {

	path = path1;
    d.path = path;
    d.cvpt_path = cvpt_path1;
    // read hyper parameter and time step from input.txt
    read_input_with_MPI();

    //  store energy of photon in system_energy
    system_energy=initial_energy;

	s.read_MPI(input, output, log);
	d.read_MPI(input, output, log, s.tlnum, s.tldim,path);
    d.construct_bright_state_MPI(input,output);

    // use symmetry to construct allowed mode combination.
    d.construct_Mode_combination_list();

    detector_only = true;
	if(energy_window) {
	    if(detector_only){ // construct matrix for detector only.
	        if(Detector_Continue_Simulation){  // load detector matrix.
                d.load_detector_Hamiltonian_MPI(path,log);
                vmode0 = d.dv_all[0];
                vmode1 = d.dv_all[1];
                d.compute_important_state_index();
                if(my_id == 0){
                    cout <<"Successfully load detector Hamiltonian" <<endl;
                }
	        }
	        else {

	            if(not compute_state_space_and_coupling_suing_symmetry_bool){
                    if(Sphere_cutoff_in_state_space){
                        compute_detector_matrix_size_MPI_sphere();
                    }
                    else{
                        compute_detector_matrix_size_MPI_cubed();
                    }
	            }
	            else{
	                construct_state_space_using_symmetry();
	            }

                d.construct_dmatrix_MPI(input,output,log,dmat0,dmat1,vmode0,vmode1);
	            if(save_state){
                    d.save_detector_Hamiltonian_MPI(path,log);
	            }
	        }
	    }
	    else {  // construct matrix for detector + photon.
	        cout << "Mode do not support. We only simulate detector in this code." << endl;
	        MPI_Abort(MPI_COMM_WORLD, -15);
        }
    }
	if(my_id ==0){
        cout<<"Finish constructing Matrix"<<endl;
        if( ! energy_window){
            cout<<" This mode is now not supported."<<endl;
            log<<" This mode is now not supported."<<endl;
            exit(-1);
        }
        // dimension_check(); // check if all matrix's dimension is right.
	}
}

// Doing Quantum Simulation with SUR algorithm, parallelized version.
void full_system::Quantum_evolution() {
    // -----------------------------------------------------------------------------------
	// Now we construct our wavefunction /phi for our detector and full_system. (For system it is already constructed in s.read())
    int i,j;
	clock_t start_time, end_time, duration;


    if(compute_eigenvalue_spectrum){
        start_time = clock();
        ofstream eigenvalue_log_file;
        ofstream eigenvalue_output_file;
        ofstream diagonal_state_mode_flie;
        if(my_id == 0){
            eigenvalue_log_file.open(path + "spectrum_log.txt");
            eigenvalue_output_file.open(path + "spectrum.txt");
            diagonal_state_mode_flie.open(path + "diagonal_state.txt");
        }
        double * eigenvalue_list = new double [d.total_dmat_size[0]];
        int numlam = 0;

        // diagonalize matrix using Lanczos algorithm. Eigenvalue store in eigenvalue list. nonzero number of eigenvalue is numlam
        if(not no_coupling){
            d.diagonalize(eigenvalue_list,numlam, eigenvalue_log_file);
        }


        if(my_id == 0){
            if(not no_coupling){
                eigenvalue_output_file << numlam << endl;
                for(i=0;i<numlam;i++){
                    eigenvalue_output_file << eigenvalue_list[i] <<" ";
                }
                eigenvalue_output_file << endl;
            }
            else{
                // no off-diagonal coupling. Thus we have states whose energy is exactly same. just output diagonal term
                eigenvalue_output_file << d.total_dmat_size[0] << endl;
                for(i=0;i<d.total_dmat_size[0];i++){
                    eigenvalue_output_file << dmat0[i] <<" ";
                }
                eigenvalue_output_file << endl;

                diagonal_state_mode_flie << d.total_dmat_size[0] << endl;
                for(i=0;i<d.total_dmat_size[0];i++){
                    for(j=0;j<d.nmodes[0];j++){
                        diagonal_state_mode_flie << vmode0[i][j] <<" ";
                    }
                    diagonal_state_mode_flie << endl;
                }
            }
        }
        if(my_id == 0){
            eigenvalue_log_file.close();
            eigenvalue_output_file.close();
            diagonal_state_mode_flie.close();
        }

        delete [] eigenvalue_list;
    }

    // compute eigenstate of system using MFD
    if(compute_overlap_with_eigenstate){
        compute_eigenstate_overlap_with_initial_state();
    }

    if(Evolve_dynamics){
        start_time = clock();

        pre_coupling_evolution_MPI(0); // pre-coupling evolution of detector state (lower bright state)

        end_time = clock();
        duration = end_time - start_time;
        if(my_id == 0) {
            log << "The total run time for parallel computing is " << (double(duration) /CLOCKS_PER_SEC)/60 << " minutes  for simulation time  " << d.proptime[0] << endl;
            cout << "The total run time for parallel computing is " << (double(duration)/CLOCKS_PER_SEC)/60 << " minutes  for simulation time  " << d.proptime[0] << endl;
        }
    }

    if(compute_eigenvector_use_MKL_module){

        vector<double> dmat0_copy = dmat0;
        sort(dmat0_copy.begin() , dmat0_copy.end());

        // compute Emin_for_core,  Emax_for_core
        allocate_diagonalization_energy_range_for_diff_proc(dmat0_copy);

        vector<double> Eigenvalue_temp ;
        vector<vector<double>> Eigenvector_temp ;

        if(my_id == 0 ){
            cout << "Min energy for system:  " << dmat0_copy[0] << endl;
            cout <<" Max energy for system:   "<< dmat0_copy[dmat0_copy.size() - 1 ] << endl;
            cout <<"dmatsize :   " << d.total_dmat_size[0] << "   dmatnum:   " << d.total_dmat_num[0] << endl;
        }

        d.eigenstate_num = MKL_Extended_Eigensolver_dfeast_scsrev_for_eigenvector_divide_by_part(d.total_dirow[0],d.total_dicol[0],d.total_dmat[0],d.total_dmat_size[0], d.total_dmat_num[0],
                                                                                                 dmat0 ,Eigenvalue_temp , Eigenvector_temp , Emin_for_core , Emax_for_core );

        // Broadcast all eigenvalue and eigenvector to all process.
        d.Broadcast_eigenstate_and_eigenvalue(Eigenvalue_temp, Eigenvector_temp);

        if(my_id == 0){
            ofstream Eigenvector_solver(path + "MKL_eigenvalue_Eigenvector.txt");


            cout <<"Finish solving eigenvector and eigenvalue " << endl;
            Eigenvector_solver.close();
        }

        if (compute_Eigenstate_OTOC_bool){

            // Use eigenstate and eigenvalue to compute OTOC.
            // first for every state construct neighbor state index.
            d.construct_neighbor_state_index_list_for_all_state();

            // compute eigenstate energy spectrum on basis set
            d.compute_eigenstate_energy_std();

            // only eigenstate whose OTOC may converge to 1 at initial time is computed. This implies eigenstate selected to compute OTOC should be in small energy window.
            d.compute_selected_eigenstate_index();

            // comput <phi_m | a_{i} | phi_l> and <phi_m | a^{+}_{i} | phi_l>. result in 3d array: phi_ladder_operator_phi
            d.compute_phi_ladder_operator_phi();

//
//
//            // compute OTOC for eigenstate solved and output to file.
//            d.compute_Eigenstate_OTOC();

            if(my_id == 0){
                cout <<" Finish computing Eigenstate OTOC ." << endl;
            }

            delete [] d.Eigenstate_energy_std_list;

            // free space for phi_ladder_operator_phi
            for(i=0;i<2*d.nmodes[0];i++){
                for(j=0;j< d.eigenstate_num;j++){
                    delete [] d.phi_ladder_operator_phi[i][j];
                }
                delete [] d.phi_ladder_operator_phi[i];
            }
            delete [] d.phi_ladder_operator_phi;

        }


        // all process free the space
        delete [] d.Eigenvalue_list;
        for (i=0;i<d.eigenstate_num;i++){
            delete [] d.Eigenstate_list[i];
        }
        delete [] d.Eigenstate_list;

    }

}

void full_system::compute_eigenstate_overlap_with_initial_state(){
    int i,j;
    int initial_state_index_in_total_dmatrix;
    int matrix_element_number;
    int numlam1 ;
    int numlam2 ;
    int min_num_lam;
    int index1 , index2;
    int last_index;
    double initial_state_energy = d.initial_Detector_energy[0];
    double hnn;
    double eigenvalue_diff1;
    double eigenvalue_diff2;
    double eigenvalue_diff3;
    double mean_eigenvalue;
    double * overlap = new double [d.total_dmat_size[0]];
    double * eigenvalue_list1 = new double [d.total_dmat_size[0]];
    double * eigenvalue_list2 = new double [d.total_dmat_size[0]];
    int * matrix_element_number_in_each_proces = new int [num_proc];
    double shift_mat_scale = 0.001;

    vector<double> eigenvalue1_after_sift;
    vector<double> eigenvalue2_after_sift;

    ofstream eigenvalue_log_file;
    ofstream eigenvalue_overlap;
    if(my_id == 0){
        eigenvalue_log_file.open(path + "MFD_log_file.txt");
        eigenvalue_overlap.open(path+"eigenvalue_overlap.txt");
    }

    initial_state_index_in_total_dmatrix = d.initial_state_index[0] + d.total_dmat_size[0] / num_proc * d.initial_state_pc_id[0];
    // record matrix element couple to bright state and its position
    vector<int> matrix_element_position;
    vector<double> matrix_element_value;
    matrix_element_number = 0;
    for(i=d.dmatsize[0];i<d.dmatnum[0];i++){
        if(d.dirow[0][i] == initial_state_index_in_total_dmatrix or d.dicol[0][i] == initial_state_index_in_total_dmatrix){
            matrix_element_number ++ ;
            matrix_element_position.push_back(i);
            matrix_element_value.push_back(d.dmat[0][i]); // coupling strength value
        }
    }

    MPI_Gather(&matrix_element_number,1,MPI_INT,&matrix_element_number_in_each_proces[0],1,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id == 0){
        for(i=0;i<num_proc;i++){
            eigenvalue_log_file << matrix_element_number_in_each_proces[i] <<" ";
        }
        eigenvalue_log_file << endl;
    }

    // increase coupling strength by 0.001
    for(i=0;i<matrix_element_number;i++){
        d.dmat[0][matrix_element_position[i]] = (1+shift_mat_scale) * matrix_element_value[i];
    }
    if(my_id == 0){
        cout << "First diagonalization: " << endl;
    }
    d.diagonalize(eigenvalue_list1,numlam1, eigenvalue_log_file);

    // decrease coupling strength by 0.001
    for(i=0;i<matrix_element_number;i++){
        d.dmat[0][matrix_element_position[i]] = (1-shift_mat_scale) * matrix_element_value[i];
    }
    if(my_id == 0){
        cout << "Second diagonalization: " << endl;
    }
    d.diagonalize(eigenvalue_list2,numlam2,eigenvalue_log_file);

    // revert dmat to original value
    for(i=0;i<matrix_element_number;i++){
        d.dmat[0][matrix_element_position[i]] = matrix_element_value[i];
    }

    if(my_id == 0){

        // we may have mismatch. sift these eigenvalue found to match them.
        index1 = 0;
        index2 = 0;
        while(index1 < numlam1 -1  and index2 < numlam2 -1){
            eigenvalue_diff1 = abs(eigenvalue_list1[index1] - eigenvalue_list2[index2]);
            eigenvalue_diff2 = abs(eigenvalue_list1[index1 + 1] - eigenvalue_list2[index2]);
            eigenvalue_diff3 = abs(eigenvalue_list1[index1] - eigenvalue_list2[index2 + 1]);
            if(eigenvalue_diff1 > eigenvalue_diff2){
                // mismatch. index1 should +1. which means list1 has one more eigenvalue found not found in list2
                index1 ++ ;
                continue;
            }
            else if (eigenvalue_diff1 > eigenvalue_diff3){
                // mismatch. index2 should +1. which means list2 has one more eigenvalue found not found in list1
                index2 ++;
                continue;
            }
            else{
                // match
                eigenvalue1_after_sift.push_back(eigenvalue_list1[index1]);
                eigenvalue2_after_sift.push_back(eigenvalue_list2[index2]);
                index1++;
                index2++;
            }

        }

        // deal with last element
        eigenvalue_diff1 = 1000;
        if(index1 == numlam1-1){
            for(j=index2;j<numlam2;j++){
                if( abs(eigenvalue_list1[index1] - eigenvalue_list2[j]) < eigenvalue_diff1 ){
                    eigenvalue_diff1 = abs(eigenvalue_list1[index1] - eigenvalue_list2[j]);
                    last_index = j;
                }
            }
            eigenvalue1_after_sift.push_back(eigenvalue_list1[numlam1 -1 ]);
            eigenvalue2_after_sift.push_back(eigenvalue_list2[last_index]);
        }
        else if (index2 == numlam2 -1){
            for(j = index1 ; j<numlam1; j++){
                if( abs(eigenvalue_list1[j] - eigenvalue_list2[index2]) < eigenvalue_diff1 ){
                    eigenvalue_diff1 = abs(eigenvalue_list1[j] - eigenvalue_list2[index2]);
                    last_index = j;
                }
            }
            eigenvalue1_after_sift.push_back(eigenvalue_list1[last_index]);
            eigenvalue2_after_sift.push_back(eigenvalue_list2[numlam2-1]);
        }

        min_num_lam = eigenvalue1_after_sift.size();
        for(i=0;i<min_num_lam;i++){
            hnn = (eigenvalue1_after_sift[i]- eigenvalue2_after_sift[i])/(2*shift_mat_scale);
            mean_eigenvalue = ( eigenvalue1_after_sift[i] + eigenvalue2_after_sift[i])/2 ;
            overlap[i] = hnn / (2 * (mean_eigenvalue - initial_state_energy) );
        }

    }

    if(my_id == 0){
        eigenvalue_overlap << initial_state_energy << endl;

        eigenvalue_overlap << min_num_lam << endl;
        for(i=0;i<numlam2;i++){
            eigenvalue_overlap << (eigenvalue1_after_sift[i] - eigenvalue2_after_sift[i])/(2*shift_mat_scale) <<" ";
        }
        eigenvalue_overlap << endl;

        eigenvalue_overlap << min_num_lam << endl;
        for(i=0;i<min_num_lam;i++){
            mean_eigenvalue = ( eigenvalue1_after_sift[i] + eigenvalue2_after_sift[i])/2 ;
            eigenvalue_overlap << mean_eigenvalue <<" ";
        }
        eigenvalue_overlap << endl;
        for(i=0;i<min_num_lam;i++){
            eigenvalue_overlap << overlap[i] << " ";
        }
        eigenvalue_overlap << endl;
    }

    delete [] eigenvalue_list1;
    delete [] eigenvalue_list2;
    delete [] overlap;
    delete [] matrix_element_number_in_each_proces;
    eigenvalue_log_file.close();
    eigenvalue_overlap.close();
}