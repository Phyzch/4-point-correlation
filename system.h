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

	double ** total_dmat_diagonal;
	double ** total_dmat;
	int ** total_dirow, ** total_dicol; // dirow, dicol, dmat in all process.

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

	vector<int> nearby_state_index; // states that we simulate dynamics, |b> for computing <a|n(t)|b>
	vector<int> states_for_4_point_correlation_average; // nearby states index for average over all these states 4 point correlation function.
    vector<int> states_for_average_in_nearby_state_index_list;
    // For states in states_for_4_point_correlation_average, we compute its index in nearby_state_index

    vector<vector<int>>  neighbor_state_index_list; // state_index with 1 quanta difference from states in nearby_state_index. [2*dof] , (1up,2up,3up, etc, 1down, 2down, 3down, etc). If not exist, choose -1.
    vector<vector<int>>  neighbor_state_in_nearby_state_index_list;
    vector<bool> bool_state_one_mode_quanta_below_all_in_nearby_state_index; // bool variable: indicating if all

    vector<bool> bool_neighbor_state_all_in_nearby_state_index;
    vector<vector<int>> neighbor_state_index_for_all_state_list;  // to compute OTOC for x and p , I need to record all states' nearby states.

	int initial_state_index_in_nearby_state_index_list;
    int initial_state_index_in_states_for_4_point_correlation_list;

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
    void update_dx(int state_number_for_evolution);
    void update_dy(int state_number_for_evolution);
    void SUR_onestep_MPI(double cf);
    void construct_bright_state_MPI(ifstream & input, ofstream & output);
    void initialize_detector_state_MPI(ofstream & log);
    void save_detector_Hamiltonian_MPI(string path, ofstream & log);
    void load_detector_Hamiltonian_MPI(string path, ofstream & log);
    void save_detector_state_MPI(string path,double * final_time,ofstream & log,int initial_state_choice);
    void load_detector_state_MPI(string path,double * start_time,ofstream & log,int initial_state_choice);

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
    void compute_n_off_diag_element(int index_b, int index_a, complex <double> * n_off_diag_element);

    void compute_n_off_diag_element_one_mode_quanta_below(int index_b, int index_a, complex<double> * n_off_diag_element);

    void replace_4_point_corr_second_line(double detector_tprint);

    // prepare variable for 4 point correlation function:
    void prepare_variable_for_4_point_correlation_function(vector<double> & dmat0, vector<double> & dmat1,ofstream & log);

    // compute density of states:
    void compute_local_density_of_state(ofstream & output,vector<double> & dmat0);

    // compute OTOC for x, p variable:
    int ** remoteVecCount_for_xp ;
    int ** remoteVecPtr_for_xp ;
    int ** remoteVecIndex_for_xp ;
    int ** Index_in_remoteVecIndex_for_xp ;
    vector<int> to_receive_buffer_len_list_with_ladder_operator;

    int ** tosendVecCount_for_xp;
    int ** tosendVecPtr_for_xp ;
    int ** tosendVecIndex_for_xp ;
    vector<int> to_send_buffer_len_list_with_ladder_operator;

    double *** send_xd_for_xp;
    double *** send_yd_for_xp;

    double *** receive_xd_for_xp;
    double *** receive_yd_for_xp;

    complex<double> compute_c_overlap(int state_m, int relative_position_to_initial_state, int mode_k,
                                      vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                      vector<vector<int>> * index_for_xp_sparsify);


    void compute_M_matrix(int state_m, int state_l, complex<double> ** M_matrix ,
                          vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                          vector<vector<int>> * index_for_xp_sparsify);
    void compute_Time_ordered_correlation_func_per_state(int state_m, int state_l, double ** M_matrix,
                                                  vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                                  vector<vector<int>> * index_for_xp_sparsify);

    double compute_two_point_func_zt_zt (int mode_k, vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                                  vector<vector<int>> * index_for_xp_sparsify , vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp );

    void output_TOC_factorization(double * Two_point_func1 , double * Two_point_func2, double ** TOC_per_state,
                                            double ** TOC,
                                            vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp,
                                            vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                            vector<vector<int>> * index_for_xp_sparsify,
                                            double t,  ofstream & TOC_output );

    void compute_Lyapunov_spectrum_for_xp(complex<double>  ** Lyapunov_spectrum_for_xp, complex<double>  ** Lyapunov_spectrum_for_xp_from_single_state,
                                          complex<double> ** M_matrix,
                                          vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp,
                                          vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                                          vector<vector<int>> * index_for_xp_sparsify,
                                          double t ,ofstream & Every_states_contribution_to_OTOC_xp);
    void prepare_computing_Lyapunovian_for_xp();

    void delete_variable_for_computing_Lyapunovian_xp();

    vector<int>  construct_receive_buffer_index_for_xp(int ** remoteVecCount_for_xp, int ** remoteVecPtr_for_xp, int ** remoteVecIndex_for_xp,
                                                       int ** Index_in_remoteVecIndex_for_xp);

    vector<int> construct_send_buffer_index_for_xp(int ** remoteVecCount_for_xp, int ** remoteVecPtr_for_xp, int ** remoteVecIndex_for_xp,
                                    int ** tosendVecCount_for_xp, int ** tosendVecPtr_for_xp, int ** tosendVecIndex_for_xp);

    void update_xd_yd_for_xp(vector<vector<double>> * xd_for_xp, vector<vector<double>> * yd_for_xp,
                             vector<vector<double>> * xd_for_xp_sparsify, vector<vector<double>> * yd_for_xp_sparsify,
                             vector<vector<int>> * index_for_xp_sparsify, double cutoff_criteria);

    void prepare_computation_for_Lanczos();
    void allocate_space_for_vector_in_Lanczos(double * &local_v, double * &local_vold,
                                              double * &local_w, double * &send_v, double * &recv_v);
    void update_v_in_Lanczos(double * local_v, double * send_v, double * recv_v);

    void compute_Lanczos(double * alpha, double * beta,
                         double * local_v, double * local_vold, double * local_w, double * send_v, double * recv_v,
                         int nlev, int maxit);

    void diagonalize(double * eigenvalue_list, int & numlam,  ofstream & eigenvalue_log_file);


    // For Boltzmann weight e^{-\beta H/4}
    double Chebyshev_e0;
    double Chebyshev_R;
    double Chebyshev_R_beta ;
    int N_chebyshev;
    double * Bessel_function_array ;

    double * shifted_total_dmat;
    double * shifted_dmat;

    vector<double> * Chebyshev_polyn;
    double ** recv_polyn;
    double ** send_polyn;

    double Chebyshev_prefactor ;
    vector<vector<double>> Boltzmann_factor_weighted_x_sparsify;  // store Boltzmann weighted basis set e^{-\beta H} | {n}>
    vector<vector<double>> Boltzmann_factor_weighted_y_sparsify;
    vector<vector<int>> Boltzmann_weighted_basis_index_sparsify;

    vector<vector<vector<double>>> ladder_operator_Boltzmann_weighted_x_sparsify; // use list to record a_{i} e^{-\beta H/4} |\phi>
    vector<vector<vector<double>>> ladder_operator_Boltzmann_weighted_y_sparsify;
    vector<vector<vector<int>>> ladder_operator_Boltzmann_weighted_basis_index_sparsify ;

    void update_polyn23();
    void prepare_compute_Boltzmann_factor_use_Chebyshev_polynomial(double one_fourth_beta , ofstream & log );
    void Chebyshev_method_Boltzmann_factor(const  vector<double> & wave_func_x ,const vector<double> & wave_func_y,
                                           vector<double> & Boltzmann_factor_weighted_wave_func_x, vector<double> & Boltzmann_factor_weighted_wave_func_y );
    void Boltzmann_factor_decorated_basis_set_and_with_ladder_operator(double sparsify_criteria =  pow(10,-2) ) ;


    // for ladder operator operation
    double ** send_xd_ladder_operator;
    double ** send_yd_ladder_operator;
    double ** recv_xd_ladder_operator;
    double ** recv_yd_ladder_operator;


    void prepare_ladder_operation();
    void ladder_operator_operation( const vector<double> & wave_func_x , const vector<double> & wave_func_y ,
                                    vector<vector<double>> & xd_for_ladder_operator ,
                                    vector<vector<double>> & yd_for_ladder_operator );

    // state_number_for evolution is number of state we have to evolve (including Haar random state.)
    int state_number_for_evolution;

    // Haar_random_vector should also be evolved in state space. Thus we need to incorporate it into nearby_state_index list.
    // We use vector<int> Haar_state_index to indicate their location. j ~ j+2N (N : dof) stores e^{-\beta H/4} |Haar> and  a_{j} e^{-\beta H/4} |Haar>
    vector<int> regularized_Haar_state_index_list;
    vector<double> Haar_state_normalization_list;

    // for regularized Lyapunov spectrum
    vector<vector<complex<double>>> regularized_thermal_Lyapunov_spectrum;
    vector<vector<vector< complex <double> >>>  regularized_thermal_Lyapunov_spectrum_each_Haar_state;
    // <m| e^{-\beta H / 4} [a_{i}(t) , a_{j} ] e^{-\beta H / 4 } |Haar> .
    vector<vector<vector<vector<  complex<double>   >>>>   regularized_thermal_OTOC_overlap_Haar_state_basis_set  ;
    // <m(t) | a_{i} e^{-iHt} a_{j} e^{-\beta H/4} |Haar>
    vector<vector<vector<vector<  complex<double>   >>>>  Haar_state_overlap_time_dependent_basis_set ;
    // a_{i} e^{-iHt} a_{j} e^{-\beta H/4} |\phi_{Haar}>   and   a_{i}e^{-iHt} e^{-\beta H/4} |\phi_{Haar}>
    vector<vector<double>> ** Haar_state_with_ladder_operator_x_sparsify;
    vector<vector<double>> ** Haar_state_with_ladder_operator_y_sparsify;
    vector<vector<int>> ** Haar_state_with_ladder_operator_basis_set_sparsify;



    // compute <Haar | e^{-\beta H } | Haar>
    void compute_normalization_factor_for_Boltzmann_weighted_factor();

    void allocate_space_for_regularized_thermal_Lyapunov_spectrum_calculation(  );

    void compute_Haar_random_state_with_ladder_operator(  double sparsify_criteria =  pow(10,-2)  );
    void compute_Haar_random_state_with_ladder_operator_overlap_with_time_dependent_basis_set(  );
    void  compute_regularized_thermal_OTOC_component(  ) ;

    void  compute_regularized_thermal_OTOC_Lyapunov_spectrum(  );


    // unregularized Lyapunov spectrum.
    vector<int> unregularized_Haar_state_index_list;

    // a_{i} e^{-iHt} a_{j} e^{-\beta H/2} |\phi_{Haar}>   and   a_{i}e^{-iHt} e^{-\beta H/2} |\phi_{Haar}>
    vector<vector<double>> ** unregularized_Haar_state_with_ladder_operator_x_sparsify;
    vector<vector<double>> ** unregularized_Haar_state_with_ladder_operator_y_sparsify;
    vector<vector<int>> ** unregularized_Haar_state_with_ladder_operator_basis_set_sparsify;

    //  <m| [a_{i}(t) a_{j}] e^{-\beta H / 2 } | Haar>
    vector<vector<vector<vector<  complex<double>   >>>>   unregularized_thermal_OTOC_overlap_Haar_state_basis_set  ;

    // for unregularized Lyapunov spectrum
    vector<vector<vector< complex <double> >>>  unregularized_thermal_Lyapunov_spectrum_each_Haar_state;
    vector<vector<complex<double>>> unregularized_thermal_Lyapunov_spectrum;

    void allocate_space_for_unregularized_Lyapunov_spectrum_calculation(  ) ;
    void compute_unregularized_Haar_random_state_with_ladder_operator (double sparsify_criteri = pow(10,-2) );
    void  compute_unregularized_thermal_OTOC_component();
    void  compute_unregularized_thermal_OTOC_Lyapunov_spectrum(  );

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
    void compute_detector_matrix_size_MPI_cubed();
    void compute_detector_matrix_size_MPI_sphere();
    void pre_coupling_evolution_MPI(int initial_state_choice);
    void construct_fullmatrix_with_energy_window_MPI();
    void compute_sstate_dstate_diagpart_dirow_dicol_MPI( );
    void construct_quotient_state_all_MPI();
    vector<vector<quotient_state>>  Gather_quotient_state_vec();
    vector<quotient_state> sort_d1list(vector<vector<quotient_state>> & d1list_each_proces);
    void Scatter_sorted_d1list(vector<quotient_state> & sorted_d1list);
    void construct_q_index_MPI();
    void rearrange_d1list();
    void compute_offdiagonal_part_MPI();
    void compute_dmat_off_diagonal_matrix_in_full_matrix_one_dmat_MPI(int index,vector < double > & mat,
            vector<int> & irow, vector<int> & icol);
    void compute_dmat_off_diagonal_matrix_in_full_matrix_MPI(vector < double > & mat,vector  <int> & irow, vector<int> & icol);
    void compute_sys_detector_coupling_MPI(vector < double > * sys_detector_mat, vector  <int> * sys_detector_irow,
            vector<int> * sys_detector_icol);
    void rearrange_off_diagonal_term(vector < double > & mat,vector  <int> & irow, vector<int> & icol);
    void  combine_offdiagonal_term(vector <double> * sys_detector_mat, vector<int> * sys_detector_irow, vector<int> * sys_detector_icol,
                                   vector<double> & d_off_mat, vector<int> & d_off_irow, vector<int> & d_off_icol,
                                   vector<double> & d_d_mat, vector<int> & d_d_irow, vector<int> & d_d_icol);
    void compute_detector_detector_coupling_MPI(vector <double> & d_d_mat, vector<int> & d_d_irow, vector<int> & d_d_icol);

    void Initial_state_MPI();

    // function to prepare and evolve photon+detector system:
    int construct_recvbuffer_index();
    void prepare_evolution();
    void  update_x_y();
    void update_x();
    void update_y();
    void full_system_SUR_one_step();
    // Output function MPI version

    void evaluate_system_output_MPI(double *hx, double * hy, double &se, double &s0, double &s1, double &s2,
                                    double &trsr2, double * de,  double ** mode_quanta, complex<double> ** sr,
                                    complex<double> ** dr, complex<double> ** total_dr);

    void etot_MPI(double * hx, double * hy);
    void compute_sys_energy_MPI(double & s0, double & s1, double &s2, double & se);

    // variable for compute detenergy_MPI:
    vector<vector<int>> dr_index_list;
    int * remote_vec_count_for_detenergy, * remote_vec_ptr_for_detenergy, * remote_vec_index_for_detenergy;
    int total_remote_vec_num_for_detenergy;
    int * to_send_vec_count_for_detenergy, * to_send_vec_ptr_for_detenergy, *to_send_vec_index_for_detenergy;
    int total_to_send_vec_num_for_detenergy;

    vector <int> local_vector_index_for_detenergy;

    double * x_for_detenergy, *y_for_detenergy;
    double * send_x_for_detenergy , * send_y_for_detenergy;

    void prepare_detenergy_computation_MPI();
    void update_x_y_for_detenergy();
    void detenergy_MPI(double * de, complex <double> ** dr, complex <double> ** total_dr);

    void average_vibrational_mode_quanta_MPI(complex <double> ** total_dr, double ** mode_quanta);
    void shift_mat();

    void gather_x_y(double * x_all, double * y_all); // gather x,y to save the state or load the state to file.
    void scatter_x_y(double * x_all, double * y_all); // scatter x_all, y_all to x,y.
    void gather_mat_irow_icol_sstate_dstate_sdmode_sdindex(double * mat_all, int * irow_all, int * icol_all ,
                                                                        int * sstate_all, int ** dstate_all,  int ** sdmode_all, int ** sdindex_all);
    void scatter_mat_irow_icol_sstate_dstate_sdmode_sdindex(double * mat_all, int * irow_all, int * icol_all ,
                                                       int * sstate_all, int ** dstate_all, int ** sdmode_all, int ** sdindex_all);
    // save and load data
    void save_Hamiltonian_MPI();
    void load_Hamiltonian_MPI();
    void save_wave_function_MPI();
    void load_wave_function_MPI();

    void Normalize_wave_function();

    // compute 4-point correlation function
    void compute_4_point_corre_for_single_state(int nearby_state_index_size, complex<double> * n_offdiag_element,
                                           complex<double> ** n_offdiag,double ** n_offdiag_real, double ** n_offdiag_imag,
                                           complex<double> **n_offdiag_total, double ** n_offdiag_total_real, double ** n_offdiag_total_imag,
                                           int initial_state_index_in_total_dmatrix ,
                                           double * four_point_correlation_function_at_initial_state);
    void compute_4_point_corre_for_multiple_states(int state_for_average_size,int nearby_state_index_size,
                                                   complex<double> * n_offdiag_element,
                                                   complex<double> *** n_offdiag_for_states_ensemble, double *** n_offdiag_for_states_ensemble_real,
                                                   double *** n_offdiag_for_states_ensemble_imag,
                                                   complex<double> *** n_offdiag_total_for_states_ensemble,
                                                   double *** n_offdiag_total_for_states_ensemble_real, double *** n_offdiag_total_for_states_ensemble_imag,
                                                   double * four_point_correlation_function_average_over_states, double * four_point_correlation_function_variance_over_states,
                                                   double ** four_point_correlation_function_for_each_states);
    void compute_Stability_Matrix( double ** Stability_Matrix,
                                                int nearby_state_index_size,
                                                complex<double> * n_offdiag_element, double * n_offdiag_element_real, double * n_offdiag_element_imag,
                                                complex<double> * n_offdiag_element_all_pc, double * n_offdiag_element_real_all_pc, double * n_offdiag_element_imag_all_pc,
                                                int initial_state_index_in_total_dmatrix );
    void compute_another_form_of_OTOC(int nearby_state_index_size, complex<double> * n_offdiag_element,
                                      complex<double> ** n_offdiag,double ** n_offdiag_real, double ** n_offdiag_imag,
                                      complex<double> **n_offdiag_total, double ** n_offdiag_total_real, double ** n_offdiag_total_imag,
                                      complex<double> ** n_offdiag_one_mode_quanta_below, double ** n_offdiag_one_mode_quanta_below_real, double ** n_offdiag_one_mode_quanta_below_imag,
                                      complex<double> ** n_offdiag_total_one_mode_quanta_below, double ** n_offdiag_total_one_mode_quanta_below_real, double ** n_offdiag_total_one_mode_quanta_below_imag,
                                      int initial_state_index_in_total_dmatrix ,
                                      double * another_OTOC);
    void compute_eigenstate_overlap_with_initial_state();
};



