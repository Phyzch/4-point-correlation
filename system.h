#pragma once
#include<iostream>
#include<iomanip>
#include <fstream>
#include <complex>
#include <vector>
#include<ctime>
#include"quotient_state.h"
using namespace std;

// define function here
void estimate_memory_cost(ofstream & resource_output);  // output resource cost to file at give time step.

// Information for MPI program
extern int my_id;
extern int num_proc;

class system {
public:
	const static int tldim = 2;
	friend class full_system;
	friend class detector;
	double *xtl, *ytl, *tle, *tlmat, *tlrho; //
	int tlnum, *tlirow, *tlicol, tlmatsize, tlmatnum;
	void read_MPI(ifstream &input, ofstream &output, ofstream &log);
    void initialize_energy_level(ifstream & input, ofstream & output);
    void initialize_wavefunction(ifstream & input, ofstream & output);
    void initialize_state_energy();
	system();
	~system();
};

class detector {
private:
    const int fillfrac = 10;
    const static int detdim =80;  // maximum n*n dimension of our detector matrix
    int dmatdim; // detector matrix size
    int stlnum;
    int stldim;
    string path;
    string cvpt_path;

    int * bright_state_index;
    int * initial_state_index;
    int * bright_state_pc_id;
    int * initial_state_pc_id;
public:
	// for mode:  modtype: =0 dark or =1 bright state.  nmax: maximum state number for each mode. nmodes: total mode numbers.
	friend class full_system;
	friend class system;
	double cf;
	int *nmodes, **nmax, **modtype;
	int *dmatsize;
    int *dmatnum , *doffnum;  // detector matrix elemetn array
	vector<int> total_dmat_size; // size of whole matrix across various process.
	vector <int> total_dmat_num, total_dmat_off_num; // total matrix number and off-diagonal number across various process.
    vector<vector<vector <int>>>  dv_all;
	int ** dmatsize_each_process;
	int ** doffnum_each_process;
	int ** dmatnum_each_process;  // record detector matrix element number in each process.
    int ** dmatsize_offset_each_process;
	int ** dmat_offset_each_process; // record local first detector matrix's index in global matrix.
	double ** total_dmat;
	int ** total_dirow, ** total_dicol; // dirow, dicol, dmat in all process.

	// == dmat0. energy of state in all process.
	vector<double> dmat_energy;

    double average_coupling_strength;

	vector<vector<vector<int>>> dv;  //dv: the q.n. for states (m,i) at coordinate j.  [m][i][j]: [detector][state][mode]
	vector<vector<int>> dirow;
	vector<vector<int>> dicol;
	int *deln;  // deln= |n_{i} - n_{j}| at coordinate k
	double *nbar;
	double **mfreq, ** modcoup, **premodcoup; // frequency of mode
	double **aij;
	vector<vector<double>>  xd,  yd; // wavefunction of detector state
    double ** xd_all, ** yd_all;
	vector<double> * dmat; // matrix
	double *proptime; // we set proptime for two different mode to be the same


    int ** remoteVecCount,  ** remoteVecPtr,  **  remoteVecIndex;
    int ** tosendVecCount,  **tosendVecPtr,  ** tosendVecIndex;
    int * to_send_buffer_len, * to_recv_buffer_len;
    double ** recv_xd,  ** recv_yd,  ** send_xd,  ** send_yd;
    vector <int > * local_dirow;
    vector <int > * local_dicol;


    //matflag	: = 1 compute detector matrixusing Madsen scaling Hamiltonian
    //maxdis : largest distance in q.n.space for which matrix element is calculated
    //cutoff : perturbation cutoff criterion V / delta - E(typ. 0.05) for states within a single detector
    //cutoff2 : same for states between different detectors : a(i)a(j) / delta - E must be greater than cutoff2
    //kelvin : detector temperature in kelvin
	int  matflag, maxdis;
	double cutoff, cutoff2, kelvin;

	double V_intra, a_intra; // intra detector coupling strength.  a_intra = coupling strength for mfreq = 50.
    double detector_energy_window_size;
	int ** bright_state, ** initial_detector_state; // record bright mode for two detectors when we try to see decoherence in our model.
	double * initial_Detector_energy;
	double * bright_state_energy;  // energy of detector's bright state.

	vector<int> sampling_state_index; // states that we simulate dynamics, |b> for computing <a|n(t)|b>
    vector<int> initial_state_in_sampling_state_index_list;
    // For states in states_for_4_point_correlation_average, we compute its index in sampling_state_index


	detector();
	~detector();
	void allocate_space(int tlnum);
    void allocate_space_single_detector(int detector_index);
	void read_MPI(ifstream & input, ofstream & output, ofstream & log, int stlnum, int stldim, string path);


    // MPI version of function.
    void construct_dv_dirow_dicol_dmatrix_MPI(ofstream & log, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void construct_dmatrix_MPI(ifstream & input, ofstream & output, ofstream & log, vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void compute_detector_offdiag_part_MPI(ofstream & log,vector<double> & dmat0,  vector<double> & dmat1,  vector<vector<int>> & vmode0, vector<vector<int>> & vmode1);
    void broadcast_total_dmat();  // broadcast dmat, dirow , dicol to form total_dirow, total_dicol, total_dmat
    void broadcast_dmatnum_doffnum(); // broadcast dmatnum, doffnum, dmatnum_each_process, doffnum_each_process, total_dmatnum etc.
    void gather_xd_yd();
    void Scatter_xd_yd();
    void Scatter_dv(vector<int> & total_mat_num);
    void Scatter_dirow_dicol_dmatrix(vector <double> * dmat_all, vector<int> * dirow_data, vector<int> * dicol_data, int ** vector_size,int **vector_displs , ofstream & log);
    int construct_receive_buffer_index(int * remoteVecCount_element, int * remoteVecPtr_element, int * remoteVecIndex_element, int detector_index);
    void prepare_evolution();
    // MPI version of SUR for one detector for each timestep.
    void update_dx(int nearby_state_list_size);
    void update_dy(int nearby_state_list_size);
    void SUR_onestep_MPI(double cf);
    void construct_initial_state_MPI(ifstream & input, ofstream & output);
    void initialize_detector_state_MPI(ofstream & log);
    void save_detector_Hamiltonian_MPI(string path, ofstream & log);
    void load_detector_Hamiltonian_MPI(string path, ofstream & log);
    // used to broadcast dv_all , vmode0, vmode1 , dmat0, dmat1
    void Broadcast_dv_all();
    // for Van Vleck transformation
    void output_detector_Hamiltonian(vector<double> & state_energy, vector<vector<int>> & dv);
    void construct_state_coupling_vanvlk(vector<double> & state_energy_local, vector<double> & state_energy, vector<vector<int>> & dv,
                                         vector <int> & dirow, vector<int> & dicol, ofstream & output);
    // use hybrid method for van vleck Hamiltonian
    void update_initial_and_bright_detector_energy();
    void compute_important_state_index();

    void output_state_density(vector<double> & dmat0,  vector<double> & dmat1);

    // compute density of states:
    void compute_local_density_of_state(ofstream & output,vector<double> & dmat0);


    void prepare_computation_for_Lanczos();
    void allocate_space_for_vector_in_Lanczos(double * &local_v, double * &local_vold,
                                              double * &local_w, double * &send_v, double * &recv_v);
    void update_v_in_Lanczos(double * local_v, double * send_v, double * recv_v);

    void compute_Lanczos(double * alpha, double * beta,
                         double * local_v, double * local_vold, double * local_w, double * send_v, double * recv_v,
                         int nlev, int maxit);

    void diagonalize(double * eigenvalue_list, int & numlam,  ofstream & eigenvalue_log_file);

    void construct_sampling_state_index(vector<double> & dmat0);

    void shift_detector_Hamiltonian(ofstream & log);
};

class full_system {
	// detector+ system
private:
	int matdim;  // matdim is maximum size of full matrix
	int matnum, offnum, matsize; // matnum is total matrix element number, it should be smaller than matdim.
	                                // in our program, usually we set matdim=matnum and use these variable interchangably.
								 //offnum: off diagonal matrix element number. matsize: matrix size for full matrix (number of diagonal element)
    int total_matnum, total_offnum, total_matsize;
    int * matsize_each_process, *mat_offnum_each_process, *matnum_each_process;
    int * matsize_offset_each_process, * matnum_offset_each_process;
    vector <double> x;
    vector <double> y; // x[matsize], y[matsize]
	vector <double> mat; // full system matrix, size: matdim
	vector <int> irow, icol; // row index and column index of matrices, size:matdim
	vector <int>sstate;
	vector<int> * dstate;  // sstate[matsize], dstate[matsize].  index of system state and detector for that matrix element. (This could be replaced by function to save space)

	int * sdnum;
	int ** sdnum_each_process;
    int * total_sd_num;
	int ** sdnum_displacement_each_process;
	vector<int> * sdindex; // index in mat for system-detector coupling elements
	vector<int> * sdmode; // mode number k for coupling
	double total_energy;
	double norm; // used to check normalization
    double total_norm;
				 // below are some variables for density matrix

    vector <quotient_state> d1list;  // state in quotient Hilbert space for detector 1
    vector <quotient_state> d2list;  // state in quotient Hilbert space for detector 2

    vector <sys_quotient_state> slist;  // state in quotient Hilbert space for system.

    // vmode,dmat for each detector.
    vector<vector<int>> vmode0;
    vector<vector<int>> vmode1;
    vector<double> dmat0;
    vector<double> dmat1;


    // used for receiving and sending vedtor x , y from/to other process
    int * remoteVecCount,  *remoteVecPtr,  * remoteVecIndex,
     * tosendVecCount, * tosendVecPtr,  * tosendVecIndex;
    double * recv_x,  * recv_y, * send_x, * send_y;
    int  to_send_buffer_len,  to_recv_buffer_len;
    vector<int>  local_irow;
    vector <int>  local_icol;

    // used for recving and sending vector x,y for computing system_energy

public:
	class system s;
	class detector d;

	// output and input of file
	ofstream output; // output result we record.
	ofstream log;  // output status (problem with code)
	ofstream resource_output;
	ifstream input;
	ofstream save;
	ifstream load;
    ofstream Detector_output;
    ofstream Detector_mode_quanta;
	// timestep variable
	double delt, tstart, tmax, tprint;
	double t; // check the time

	double cf; // energy scale

	string path;

	full_system(string path1,string cvpt_path1);
	void dimension_check();
	void Quantum_evolution();;
	void replace_first_line(); // just ignore this code, this code do the very dumb work..


    // MPI version of code:
    void read_input_with_MPI();
    void compute_detector_matrix_size_MPI_sphere();
    void pre_coupling_evolution_MPI(int initial_state_choice);



    void compute_eigenstate_overlap_with_initial_state();
};



