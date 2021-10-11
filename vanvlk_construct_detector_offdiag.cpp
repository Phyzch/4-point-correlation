//
// Created by phyzch on 9/12/20.
//

#include"system.h"
#include"util.h"
using namespace std;
int ndegre;
int ndegrx2;

int bin_number = 1000;
void Broadcast_vanvlk_data(int & norder, vector<double> & frequency, vector<double> & energy_decrease_list,
                           vector<double> & Coeff, vector<vector<int>> & Normal_form);
struct sort_element{
    double energy;
    double Coeff;
    vector<int> Coordinate_index;
    sort_element(double energy1,double Coeff1, vector<int> & Coordinate_index1){
        energy= energy1;
        Coeff= Coeff1;
        Coordinate_index = Coordinate_index1;
    }
};


void Read_ham4x_out( string  path, int & norder, vector<double> & frequency, vector<double> & energy_decrease_list,
                     vector<double> & Coeff, vector<vector<int>> & Normal_form){
    // ham4x.out contain Hamiltonian which is the one before VaN Vleck perturbation.
    // By comparing the result using Hamiltonian in ham4x.out and the one in vanvlk.out,
    // we can tell if van vleck transformation is successful or not.
    ifstream ham4x_out (path + "ham4x.out");
    ifstream param_input (path+"param.input");
    assert(! param_input.fail() and !ham4x_out.fail());
    int i,j;
    int iorder;
    string ss;
    double f;
    double coefficient;
    int normal_form_number;
    double energy_decrease;
    int numterm_in_order;
    int number_of_coupling=0;

    getline(param_input,ss);
    param_input >> norder;
    param_input.close();

    frequency.clear();
    for(i=0;i<ndegre;i++){
        ham4x_out>> f;
        frequency.push_back(f);
        getline(ham4x_out,ss);
    }
    // read 0 order term
    ham4x_out >> numterm_in_order;
    getline(ham4x_out,ss);
    for(i=0;i<numterm_in_order;i++){
        getline(ham4x_out,ss);
    }

    for (iorder = 1; iorder<= norder; iorder++){
        ham4x_out >> numterm_in_order;
        number_of_coupling = number_of_coupling + numterm_in_order;
        getline(ham4x_out,ss);

        for(i=0;i<numterm_in_order;i++){
            ham4x_out >>  coefficient;
            Coeff.push_back(coefficient);
            vector <int> normal_form_element;
            for(j=0;j<ndegrx2;j++){
                ham4x_out >> normal_form_number;
                normal_form_element.push_back(normal_form_number);
            }
            Normal_form.push_back(normal_form_element);  // get normal form
            getline(ham4x_out,ss);
        }
    }
    cout<<" Total Number of coupling in ham4x.out is "<< number_of_coupling<<endl;

    for(i=0;i<number_of_coupling;i++){
        energy_decrease = 0;
        for(j=0;j<ndegre;j++){
            energy_decrease = energy_decrease + frequency[j] * (Normal_form[i][j+ndegre] - Normal_form[i][j]);
        }
        energy_decrease_list.push_back(energy_decrease);
    }
    ham4x_out.close();
}

void Read_Vanvleck_output( string  path, int & norder, vector<double> & frequency, vector<double> & energy_decrease_list,
                           vector<double> & Coeff, vector<vector<int>> & Normal_form){
    // read data from vanvlk.out
    // Coeff store coefficient for normal form operator
    // Normal_form store index for normal form
    // energy decrease list is energy decrease caused by operator.
    ifstream Vanvleck_out(path + "vanvlk.out");
    assert(!Vanvleck_out.fail());
    int i,j;
    string ss;
    double f;
    double coefficient;
    int normal_form_number;
    int iorder;
    int numterm_in_order;  // number of coupling term in each order
    int number_of_coupling = 0;  // total number of coupling
    double energy_decrease;
    Vanvleck_out >> norder;
    for(i=0;i<9;i++){
        getline(Vanvleck_out,ss);
    }

    frequency.clear();
    for(i=0;i<ndegre;i++){
        Vanvleck_out >> f ;
        frequency.push_back(f);
    }
    for(i=0;i<2;i++) getline(Vanvleck_out,ss);  // readlines in vanvleck.out

    for(iorder = 1; iorder<= norder; iorder++){
        Vanvleck_out >> numterm_in_order; // number of coupling term in each order
        number_of_coupling = number_of_coupling + numterm_in_order;
        getline(Vanvleck_out,ss);

        for(i=0;i<numterm_in_order;i++){
            Vanvleck_out>> coefficient;
            Coeff.push_back(coefficient); // get coefficient
            vector <int> normal_form_element;
            for(j=0;j<ndegrx2;j++){
                Vanvleck_out >> normal_form_number;
                normal_form_element.push_back(normal_form_number);
            }
            Normal_form.push_back(normal_form_element);  // get normal form
            getline(Vanvleck_out,ss);
        }
        getline(Vanvleck_out,ss);
    }

    // Compute energy decrease for each coupling term
    for(i=0;i<number_of_coupling;i++){
        energy_decrease = 0;
        for(j=0;j<ndegre;j++){
            energy_decrease = energy_decrease + frequency[j] * (Normal_form[i][j+ndegre] - Normal_form[i][j]);
        }
        energy_decrease_list.push_back(energy_decrease);
    }
    Vanvleck_out.close();
}


vector<sort_element> merge_sort(const vector<sort_element> & v1, const vector<sort_element> & v2){
    int size1= v1.size();
    int size2= v2.size();
    if(size1 == 0) return v2;
    if(size2 == 0) return v1;
    int v1_index = 0;
    int v2_index = 0;
    double energy_difference;
    int i;
    vector <sort_element> v3;
    while(v1_index<size1 and v2_index<size2){
        energy_difference = v1[v1_index].energy - v2[v2_index].energy;
        if(energy_difference < 0){
            v3.push_back(v1[v1_index]);
            v1_index++;
        }
        else if (energy_difference>0){
            v3.push_back(v2[v2_index]);
            v2_index++;
        }
        else{
            // energy change is the same
            v3.push_back(v1[v1_index]);
            v3.push_back(v2[v2_index]);
            v1_index++;
            v2_index++;
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

void order_coupling_term(vector<double> & energy_decrease_list, vector<double> & Coeff, vector<vector<int>> & Normal_Form){
    // order coupling term by its energy change
    // use merge sort algorithm
    int i,j;
    int coupling_num = Normal_Form.size();
    vector<sort_element> list_to_sort;
    // construct sort_element for merge_sort
    for(i=0;i<coupling_num;i++){
        sort_element element(energy_decrease_list[i],Coeff[i],Normal_Form[i]);
        list_to_sort.push_back(element);
    }

    // start merge sorting
    vector <vector<sort_element>> List_for_list_old;
    vector<vector<sort_element>> * old_ptr = & List_for_list_old;
    vector<vector<sort_element>> List_for_list_new;
    vector<vector<sort_element>> * new_ptr = & List_for_list_new;
    vector<vector<sort_element>> * list_ptr_3;
    vector<sort_element> v3;
    int list_size = coupling_num;
    for(i=0;i<coupling_num;i++){
        vector<sort_element> small_list;
        small_list.push_back(list_to_sort[i]);
        List_for_list_old.push_back(small_list);
    }

    while(list_size>1){
        for(i=0;i+1<list_size;i=i+2){
            v3 = merge_sort( (*old_ptr)[i], (*old_ptr)[i+1] );
            (*new_ptr).push_back(v3);
        }
        if(list_size % 2 == 1){
            (*new_ptr).push_back( (*old_ptr)[list_size -1] );
        }
        list_size = (list_size + 1 ) /2;
        //exchange two list
        (*old_ptr).clear();
        (*old_ptr).shrink_to_fit();

        list_ptr_3 = old_ptr;
        old_ptr = new_ptr;
        new_ptr = list_ptr_3;
    }
    list_to_sort.clear();
    energy_decrease_list.clear();
    Coeff.clear();
    Normal_Form.clear();
    list_to_sort = (*old_ptr)[0];  // sorting result store in *old_ptr
    for(i=0;i<coupling_num;i++){
        energy_decrease_list.push_back(list_to_sort[i].energy);
        Coeff.push_back(list_to_sort[i].Coeff);
        Normal_Form.push_back(list_to_sort[i].Coordinate_index);
    }
}

int * Bin_coupling_term(vector<double> & energy_decrease_list, vector<double> & Coeff, vector<vector<int>> & Normal_Form){
    // Binning the coupling term according to energy_change
    // return mcount: begin index for binning in ordered list
    order_coupling_term(energy_decrease_list,Coeff,Normal_Form);
    int i,j;
    int bin_index;
    int coupling_number = energy_decrease_list.size();
    double min_energy = energy_decrease_list[0] - 0.001;
    double max_energy = energy_decrease_list[coupling_number - 1] + 0.001;
    int * mcount = new int [bin_number+1];
    for(i=0;i<bin_number+1;i++){
        mcount[i] = 0;
    }
    for(i=0;i<coupling_number;i++){
        bin_index = 1 + (energy_decrease_list[i] - min_energy)/(max_energy - min_energy) * bin_number;
        mcount[bin_index] ++;
    }
    for(i=1;i<=bin_number;i++){
        mcount[i] = mcount[i] +mcount[i-1];
    }

    // Now element in bin has the index in ordered list from range mcount[i-1]+1 ~ mcount[i]
    return mcount;
}

double factorial(int i, int j, int n){
    // compute sqrt(i! * j!) /(j-n)!   here n is order of lowering operator.
    // Which is coefficient we get when apply normal form operator to connect two states <i| (a^{+})^{i-j+n} a^{n}|j>

    double prod =1;
    int k;
    int m;
    for(k=0;k<n;k++){
        prod = prod * (j-k);
    }
    m= i - (j-n);
    for(k=0;k<m;k++){
        prod = prod * (i-k);
    }
    prod = sqrt(prod);
    return prod;
}

void vmat(vector<double> & state_energy_change,vector<double> & state_energy_local, vector<double> & state_energy, vector<vector<int>> & dv,
          vector <int> & dirow, vector<int> & dicol,
          vector<double> & energy_decrease_list, vector<double> & Coeff,
          vector<vector<int>> & Normal_Form, int * mcount, int i, int j, int bin_index,
          int matrix_size, double cutoff, double random_number){
    // <i| operator |j>
    int k,l;
    double Prod;
    double lij;
    double Vmat= 0;
    double index_diff;
    // local index is index of i in state_energy_local
    int local_index;
    int begin_index, end_index;
    local_index = i - (matrix_size/num_proc)*my_id ;

    if(i==j){
        // diagonal term
        for(k=mcount[bin_index -1] ; k<mcount[bin_index]; k++){  // k is index for coupling operator
            for(l=0;l<ndegre;l++){
                if(Normal_Form[k][l] != Normal_Form[k][l+ndegre]){
                    goto label1;  // for contribution to diagonal term, order of raising operator should equal to lowering operator.
                }
            }

            Prod = 1;
            for(l=0;l<ndegre;l++){
                // here Normal_Form[k][l+ndegre] is order of lowering operator for operator in list with index k
                Prod = Prod * factorial(dv[i][l], dv[j][l], Normal_Form[k][l+ndegre]);
            }
            Vmat = Vmat + Prod * Coeff[k];
            label1:;
        }

 // we use  random number as anharmonicity, this way, we get Logan Wolynes model: energy difference between nearby state is random
 // Turn on random anharmonicity
        if (turn_on_random_self_anharmonicity){
           Vmat = random_number;
        }
    }
    else{
        // off-diagonal term
        begin_index = mcount[bin_index - 1];
        end_index = mcount[bin_index];

        // rotational dof index should equal to each other for vib term to have effect.
        if(dv[i][ndegre] != dv[j][ndegre] ) goto label3;

        for(k=begin_index;k<end_index;k++){
            if (k == mcount[bin_number]) break;
            for(l=0;l<ndegre;l++){
                index_diff = Normal_Form[k][l] - Normal_Form[k][l+ndegre];
                if( dv[i][l] - dv[j][l] != index_diff ) goto label2;
                if( dv[j][l] < Normal_Form[k][l+ndegre]) goto label2; // lowering operator first reach |0> here
            }

            Prod=1;
            for(l=0;l<ndegre;l++){
                Prod = Prod * factorial( dv[i][l], dv[j][l], Normal_Form[k][l+ndegre] );
            }
            Vmat = Vmat + Prod * Coeff[k];
            label2:;
        }

        label3:;
    }

    if(i!=j) {
        // off-diagonal term
        lij = abs(Vmat) / abs(state_energy[i] - state_energy[j]);
        if (lij > cutoff) {
            state_energy_local.push_back(Vmat);
            dirow.push_back(i);
            dicol.push_back(j);
        }
    }
    else{
        // diagonal term
        // state_energy represent dmat0/ dmat1, which is detector matrix diagonal form across all process.
        // It will be used in compute_sstate_dstate_diagpart_dirow_dicol_MPI() to construct full matrix so we have to update it.
        state_energy_change[i] = state_energy_change[i] + Vmat;
        state_energy_local[local_index] = state_energy_local[local_index] + Vmat;
    }
}


void  construct_state_coupling_subroutine(vector<double> & state_energy_local ,vector<double> & state_energy, vector<vector<int>> & dv,
                                          vector <int> & dirow, vector<int> & dicol,
                                          vector<double> & energy_decrease_list, vector<double> & Coeff,
                                          vector<vector<int>> & Normal_Form, int * mcount, double cutoff, ofstream & output){
    // state_energy is energy of state (diagonal term)
    // dv is coordinate in state .
    // dirow and dicol is row index and column index for matrix.  matrix is sparse here
    // energy_decrease_list (ordered) is decrease of energy caused by operator. This is also energy difference between original state and state it couple to.
    // Coeff is coefficient of operator.
    // Normal_Form is normal form for operator.
    // state_energy_local is dmat[0] \ dmat[1] , which is local detector matrix in each process.
    // state_energy is dmat0 or dmat1, which is detector matrix (diagonal part) shared by all process

    int matrix_size = state_energy.size();
    int i,j;
    double energy_difference;
    int bin_index;
    int coupling_operator_number = energy_decrease_list.size();
    double min_energy = energy_decrease_list[0] - 0.001;
    double max_energy = energy_decrease_list[coupling_operator_number-1] + 0.001 ;
    vector<double> original_state_energy = state_energy;  // original state energy is energy level before perturbed by Van Vleck  transformation (dmat0 or dmat1)
    vector<double> state_energy_change;
    vector<double> state_energy_change_in_all_process;
    state_energy_change.resize(matrix_size);
    state_energy_change_in_all_process.resize(matrix_size);
    int begin_index, end_index;  // begin index is beginning index for this process.
    begin_index = matrix_size / num_proc * my_id ;
    end_index = begin_index + state_energy_local.size();

   double random_number;
   double timee = time(nullptr);
   double self_anharmonicity_strength = 5000 ;
   std::default_random_engine generator(timee);
   std::normal_distribution<double> distribution(0,self_anharmonicity_strength);   // used to generate random anharmonicity

   if(turn_on_random_self_anharmonicity){
       output<<"Self anharmonicity on energy level:   "<<self_anharmonicity_strength << endl;
   }

    // first compute state energy shift: (state_energy_change)
    for(i=begin_index;i<end_index;i++){
        j=i;
        energy_difference = 0;
        bin_index = 1 + (energy_difference - min_energy) / (max_energy - min_energy) * bin_number;
        random_number = distribution(generator);
        if(bin_index >=1 and bin_index < bin_number+1) {
            vmat(state_energy_change,state_energy_local,state_energy, dv, dirow, dicol, energy_decrease_list, Coeff, Normal_Form, mcount, i, j,
                 bin_index,matrix_size,cutoff,random_number);  // compute coupling strength and update off-diagonal part and diagonal part correction
        }
    }
    MPI_Allreduce(&state_energy_change[0],&state_energy_change_in_all_process[0],matrix_size,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<matrix_size;i++){
        state_energy[i] = state_energy[i] + state_energy_change_in_all_process[i];
    }
    // now we compute coupling between different state.
    for(i=begin_index;i<end_index;i++){
        for(j=0;j<matrix_size;j++){
            if(j!=i) {
                energy_difference = original_state_energy[j] - original_state_energy[i];   // <i| operator | j>
                // to find the right bin for energy change of opeartor, we should use original state energy.
                bin_index = 1 + (energy_difference - min_energy) / (max_energy - min_energy) * bin_number;
                random_number = distribution(generator);
                if (bin_index >= 1 and bin_index < bin_number + 1) {
                    vmat(state_energy_change, state_energy_local, state_energy, dv, dirow, dicol, energy_decrease_list,
                         Coeff, Normal_Form, mcount, i, j,
                         bin_index, matrix_size,
                         cutoff,random_number);  // compute coupling strength and update off-diagonal part and diagonal part correction
                }
            }
        }
    }
}

void detector::construct_state_coupling_vanvlk(vector<double> & state_energy_local, vector<double> & state_energy, vector<vector<int>> & dv,
                                               vector <int> & dirow, vector<int> & dicol,ofstream & output){
    // state_energy_local is dmat[0] \ dmat[1] , which is local detector matrix in each process.
    // state_energy is dmat0 or dmat1, which is detector matrix (diagonal part) shared by all process
    int i ;
    vector <double> frequency;
    vector <double> Coeff;
    vector <double> energy_decrease_list;
    vector<vector<int>> Normal_Form;
    int norder;
    int * mcount;
    bool match = true;
    if(my_id == 0) {
        if (Turn_on_Vanvleck) {
            Read_Vanvleck_output(cvpt_path, norder, frequency, energy_decrease_list, Coeff, Normal_Form);
        }
        else {
            Read_ham4x_out(cvpt_path, norder, frequency, energy_decrease_list, Coeff, Normal_Form);
        }
                // check if frequency match
        for(i=0;i<ndegre;i++){
            if(mfreq[0][i] != frequency[i] ) match = false;
        }
        if(match == false){
            cout <<"Frequency doesn't match!" << endl;
            MPI_Abort(MPI_COMM_WORLD,-13);
        }

        ifstream pot_dat(cvpt_path + "pot.dat");
        assert(!pot_dat.fail());
        pot_dat >> a_intra >> V_intra;
        pot_dat.close();
    }
    MPI_Bcast(&a_intra,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&V_intra,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    Broadcast_vanvlk_data(norder,frequency,energy_decrease_list,Coeff,Normal_Form);
    // Bin the coupling term
    mcount = Bin_coupling_term(energy_decrease_list,Coeff,Normal_Form);
    // compute diagonal correction to matrix and off-diagonal term there
    construct_state_coupling_subroutine(state_energy_local, state_energy,dv,dirow,dicol,energy_decrease_list,Coeff,Normal_Form,mcount, cutoff,output);
}

void Broadcast_vanvlk_data(int & norder, vector<double> & frequency, vector<double> & energy_decrease_list,
                           vector<double> & Coeff, vector<vector<int>> & Normal_form){
    int operator_number=0;
    int * normal_form_1d;
    int index = 0;
    int i,j;
    vector<int> normal_form_element;
    if(my_id == 0){
        operator_number = Coeff.size();
    }
    MPI_Bcast(&operator_number, 1, MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&norder,1,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id != 0) {
        frequency.clear();
        frequency.resize(ndegre);
        energy_decrease_list.resize(operator_number);
        Coeff.resize(operator_number);
    }
    MPI_Bcast(&frequency[0],ndegre,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&energy_decrease_list[0],operator_number,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Coeff[0],operator_number,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // broadcast 2d vector Normal form:
    normal_form_1d = new int [ndegrx2*operator_number];
    if(my_id == 0){
        // copy Normal form to normal_form_1d
        index = 0;
        for(i=0;i<operator_number;i++){
            for(j=0;j<ndegrx2;j++){
                normal_form_1d[index] = Normal_form[i][j];
                index ++ ;
            }
        }
    }
    MPI_Bcast(&normal_form_1d[0],ndegrx2*operator_number,MPI_INT,0,MPI_COMM_WORLD);
    if(my_id!=0){
        index = 0;
        for(i=0;i<operator_number;i++){
            normal_form_element.clear();
            for(j=0;j<ndegrx2;j++){
                normal_form_element.push_back(normal_form_1d[index]);
                index ++;
            }
            Normal_form.push_back(normal_form_element);
        }
    }
    delete [] normal_form_1d;

}

double compute_edge_criteria(vector<int> dv_element, int nmode){
    // we use criteria from paper:  https://doi.org/10.1063/1.3105989 : es
    int i;
    double es = 0;
    double dv_average = 0;
    for(i=0;i<nmode;i++){
        dv_average = dv_average + double(dv_element[i])/nmode;
    }
    for(i=0;i<nmode;i++){
        es = es + pow(dv_element[i]/dv_average - 1, 2);
    }
    es = es / (nmode * (nmode -  1));
    es = sqrt(es);
    return es;
}

vector<bool> check_edge_state(const vector<vector<int>> & dv ){
    int i;
    int state_number = dv.size();
    int nmode = dv[0].size();
    double es; // criteria for deciding if it's edge state
    vector<bool> edge_state_mark; // for edge state, its mark will be true
    for(i=0;i<state_number;i++){
        es = compute_edge_criteria(dv[i],nmode);
        if(es > 0.8){  // criteria for edge state
            edge_state_mark.push_back(true); // edge state
        }
        else{
            edge_state_mark.push_back(false);
        }
    }
    return  edge_state_mark;
}


void detector:: output_detector_Hamiltonian(vector<double> & state_energy, vector<vector<int>> & dv){
    ofstream Hamiltonian(path+"Hamiltonian.txt");
    int m,i,j;
    int dmat_size = state_energy.size();
    Hamiltonian << "Hamiltonian for detector  "<<endl;
    for(i=0;i<dmat_size;i++){
        Hamiltonian <<i<<"  "<< state_energy[i] <<"   ";
        for(j=0;j<ndegre;j++){
            Hamiltonian << dv[i][j] <<"  ";
        }
        Hamiltonian <<endl;
    }
    for(i=0;i<10;i++){
        Hamiltonian << endl;
    }
}