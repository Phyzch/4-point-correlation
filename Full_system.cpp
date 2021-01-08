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
	            if(Sphere_cutoff_in_state_space){
                    compute_detector_matrix_size_MPI_sphere();
	            }
	            else{
                    compute_detector_matrix_size_MPI_cubed();
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
        dimension_check(); // check if all matrix's dimension is right.
	}
}

// Doing Quantum Simulation with SUR algorithm, parallelized version.
void full_system::Quantum_evolution() {
    // -----------------------------------------------------------------------------------
	// Now we construct our wavefunction /phi for our detector and full_system. (For system it is already constructed in s.read())
    int i;
	clock_t start_time, end_time, duration;

    if(Evolve_dynamics){
        start_time = clock();

        pre_coupling_evolution_MPI(0); // pre-coupling evolution of detector state (lower bright state)

        end_time = clock();
        duration = end_time - start_time;
        if(my_id == 0) {
            log << "The total run time for parallel computing is " << (double(duration) /CLOCKS_PER_SEC)/60 << " minutes  for simulation time  " << tmax << endl;
            cout << "The total run time for parallel computing is " << (double(duration)/CLOCKS_PER_SEC)/60 << " minutes  for simulation time  " << tmax << endl;
        }
    }

    if(compute_eigenvalue_spectrum){
        start_time = clock();
        ofstream eigenvalue_log_file;
        ofstream eigenvalue_output_file;
        if(my_id == 0){
            eigenvalue_log_file.open(path + "spectrum_log.txt");
            eigenvalue_output_file.open(path + "spectrum.txt");
        }
        double * eigenvalue_list = new double [d.total_dmat_size[0]];
        int numlam = 0;
        // diagonalize matrix using Lanczos algorithm. Eigenvalue store in eigenvalue list. nonzero number of eigenvalue is numlam
        if(my_id == 0){
            if(not no_coupling){
                d.diagonalize(eigenvalue_list,numlam, eigenvalue_log_file);
            }
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
            }
        }
        if(my_id == 0){
            eigenvalue_log_file.close();
            eigenvalue_output_file.close();
        }

        delete [] eigenvalue_list;
    }

}
