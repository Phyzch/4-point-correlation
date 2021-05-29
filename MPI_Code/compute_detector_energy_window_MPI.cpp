//
// Created by phyzch on 7/22/20.
//
# include"../util.h"
# include "../system.h"

using namespace std;
int Rmax; // maximum distance allowed in detector state space.

double detector_lower_bright_state_energy_window_shrink;
int * compute_state_space(const vector <double> & dmat0, const vector <double> & dmat1, double * initial_state_energy, double * bright_state_energy){
    /* input:  dmat0, dmat1: energy for detector state
               initial_state_energy: energy for detector lower bright state, [0] for detector1,  [1] for detector2
        output: state_number: 0-3.  0 for detector0 lower bright state, 1 for detector0 higher bright state
        2 for detector 1 lower bright state, 3 for detector 1 higher bright state
    */
    int * state_number = new int [4];
    int i;
    for(i=0;i<4;i++){
        state_number[i]=0;
    }
    int dmat0_size=dmat0.size();
    int dmat1_size= dmat1.size();
    for(i=0;i<dmat0_size;i++){
        if(dmat0[i]<initial_state_energy[0] + energy_window_size and dmat0[i] > initial_state_energy[0] - energy_window_size){
            state_number[0]++;
        }
        if(dmat0[i] < bright_state_energy[0] + energy_window_size and  dmat0[i] > bright_state_energy[0] - energy_window_size){
            state_number[1] ++;
        }
    }
    for(i=0;i<dmat1_size;i++){
        if(dmat1[i]<initial_state_energy[1] + energy_window_size and dmat1[i] > initial_state_energy[1] - energy_window_size){
            state_number[2] ++;
        }
        if(dmat1[i]<bright_state_energy[1] + energy_window_size and dmat1[i] > bright_state_energy[1] - energy_window_size){
            state_number[3]++;
        }
    }
    return state_number;
}


int state_distance(const vector<int> & ndetector, int * state1, int moddim){
    // compute distance between two state : ndetector and state1. Dimension of mod is moddim
    int i;
    int distance=0;
    for(i=0;i<moddim;i++){
        distance= distance + abs(state1[i] - ndetector[i]);
    }
    return distance;
}

int state_max_quanta_diff(const vector<int> & ndetector, int * state1, int moddim){
    int i;
    int max_quanta_diff = 0;
    for(i=0;i<moddim;i++){
        if(max_quanta_diff < abs(state1[i] - ndetector[i])){
            max_quanta_diff = abs(state1[i] - ndetector[i]);
        }
    }
    return max_quanta_diff;
}

int vmode_compare(const vector <int> & lhs, const vector <int> & rhs){
    // compare function used for find_position_for_insert_binary.
    int i,j;
    int size= lhs.size();
    if(size != rhs.size()){
        cout<<"size of vector do not match."<<endl;
        exit(-5);
    }
    for(i=0;i<size;i++){
        if(lhs[i] > rhs[i]) return 1;
        else if(lhs[i] < rhs[i]) return -1;
    }
    return 0;
}

//called in compute_matrix_size() function.
void copy_array(vector<int> & temporary_vector, const vector<int> & ndetector,int modnumber) {
    // code used in compute_matrix_size() function
    // modnumber here is the size of ndetector (the two detectors have the same size)
    int i;
    for (i = 0; i < modnumber; i++) {
        temporary_vector[i] = ndetector[i];
    }
}


// find position to insert the corresponding detector state, called in compute_matrix_size() function.  Binary search.
int find_position_for_insert_binary(const vector<vector<int>> & vmode, const vector<int> & ndetector, bool & exist) {
    // code used in compute_matrix_size() function
    // first compare first mode, 0 is in front of 1. So finally our mode is |00000> , |00001>,|00010>, etc, |10000>
    if (vmode.size() ==0){
        exist = false;
        return 0;
    }
    int left_flag=0;
    int right_flag=vmode.size() -1;
    int position = (left_flag + right_flag) /2;
    exist = false;
    int mark;
    // binary search using compare function
    while(right_flag>left_flag) {
        mark = vmode_compare(ndetector, vmode[position]);
        if (mark == 1) {
            // ndetector should be after position
            left_flag = position+1;
            position = (left_flag + right_flag) / 2;
        }
        else if ( mark == -1) {
            right_flag = position-1;
            position = (left_flag + right_flag) / 2;
        }
        else if (mark== 0) {
            exist = true;
            return position;
        }
    }
    // now left == right. Decide position to insert and return it.
    mark = vmode_compare(ndetector, vmode[position]);
    if( mark == 0){
        exist=true;
        return position;  // we will not insert.
    }
    else if ( mark == -1){
        exist=false;
        return position;   // we will insert before this position.
    }
    else if ( mark == 1){
        exist=false;
        return position +1; // we will insert after this position.
    }

    //
}


void full_system:: compute_detector_matrix_size_MPI_sphere( ){
    // use this function we construct detector state in a cube.
    if(my_id==0){
        int i, j, k;
        int i1, i2;
        double  detector0_energy, detector1_energy;
        // ndetector0 and ndetector1 indicate current detector mode index we are calculating energy.
        vector<int> ndetector0(d.nmodes[0]);
        vector <int> ndetector1(d.nmodes[1]);
        // record size of total matrix
        int location;
        bool exist=false;

        int excited_bright_state_distance;
        int lower_bright_state_distance;

        double middle_state_energy = (d.initial_Detector_energy[0] + d.initial_Detector_energy[1])/2;
        double high_initial_state_energy = max(d.initial_Detector_energy[0] , d.initial_Detector_energy[1]);
        double low_initial_state_energy = min(d.initial_Detector_energy[0],d.initial_Detector_energy[1]);

        ndetector0[0] = -1; // this is for:  when we go into code: ndetector0[i]= ndetector0[i]+1, our first state is |000000>
        while (1) {
            label2:;  // label2 is for detector1 to jump out of while(1) loop (this is inner layer of while(1))
            detector0_energy = 0;
            for (i1 = 0; i1 < d.nmodes[0]; i1++) {  // loop through detector0
                // define the way we loop through detector0:
                ndetector0[i1] = ndetector0[i1] + 1;
                if (ndetector0[i1] <= d.nmax[0][i1]) break;
                if (ndetector0[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
                    ndetector0[d.nmodes[0] - 1] = 0;
                    goto label1;  // use goto to jump out of nested loop
                }
                ndetector0[i1] = 0;
            }
            // calculate detector 0 energy
            for (i = 0; i < d.nmodes[0]; i++) {
                detector0_energy = detector0_energy + ndetector0[i] * d.mfreq[0][i];
            }

            //--------------------------------------------------------------------------------------------
            // criteria below make sure detector 0 's energy is reasonable.
            if (detector0_energy > d.initial_Detector_energy[0] + d.detector_energy_window_size) {
                // detector 0's energy can not be larger than its initial energy + photon energy
                // jump to next detector state.
                k = 0;
                while (ndetector0[k] == 0) {
                    ndetector0[k] = d.nmax[0][k];
                    k++;
                    if (k >= d.nmodes[0]) {
                        break;
                    }
                }
                if (k < d.nmodes[0]) {
                    ndetector0[k] = d.nmax[0][k];
                }
                goto label2;
            }

            //------------------------------------------------------------------------------------------------
            // criteria below make sure detector 1 can not be too far away from bright state and lower bright state.

            lower_bright_state_distance = state_distance(ndetector0, d.initial_detector_state[0], d.nmodes[0]);
            // we do not use distance constraint for state whose energy is between two
            if ( lower_bright_state_distance > Rmax) {
                goto label2;
            }


            //--------------------------------------insert this state in detector's state.-----------------------------------------------------------
            location=find_position_for_insert_binary(vmode0, ndetector0, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                vmode0.insert(vmode0.begin() + location, ndetector0);
                dmat0.insert(dmat0.begin() + location, detector0_energy);
            }
        }
        label1:;

        ndetector1[0] = -1; // this is when we go into code: ndetector1[i] = ndetector1[i]+1. our first state is |000000>
        while (1) { // loop through detector 1
            label3:;
            detector1_energy = 0;
            for (i2 = 0; i2 < d.nmodes[1]; i2++) {
                // define the way we loop through detector1
                ndetector1[i2] = ndetector1[i2] + 1;
                if (ndetector1[i2] <= d.nmax[1][i2]) break;
                if (ndetector1[d.nmodes[1] - 1] > d.nmax[1][d.nmodes[1] - 1]) {
                    ndetector1[d.nmodes[1] - 1] = 0;
                    goto label4;
                }
                ndetector1[i2] = 0;
            }
            // calculate detector 1 energy
            for (i = 0; i < d.nmodes[1]; i++) {
                detector1_energy = detector1_energy + ndetector1[i] * d.mfreq[1][i];
            }
            // --------------------------------------------------------------
            //  criteria below make sure detector 1's energy is reasonable:
            if (detector1_energy > d.initial_Detector_energy[1] + d.detector_energy_window_size) {
                // initial energy is system energy.
                // detector 1 's energy can not be larger than its initial energy + photon energy
                j = 0;
                while (ndetector1[j] == 0) { // go to first mode whose n!=0;
                    ndetector1[j] = d.nmax[1][j];
                    j++;
                    if (j >= d.nmodes[1]) {
                        break;
                    }
                }
                if (j < d.nmodes[1]) {
                    ndetector1[j] = d.nmax[1][j];
                }
                goto label3;
            }

            //------------------------------------------------------------------------------------------------
            // criteria below make sure detector 1 can not be too far away from bright state and lower bright state.
            lower_bright_state_distance = state_distance(ndetector1, d.initial_detector_state[1], d.nmodes[1]);
            if ( lower_bright_state_distance > Rmax ) {
                goto label3;
            }
            location = find_position_for_insert_binary(vmode1, ndetector1, exist);
            if (!exist) {
                vmode1.insert(vmode1.begin() + location, ndetector1);
                dmat1.insert(dmat1.begin() + location, detector1_energy);
            }
        }
        label4:;
        // add function here to count state number in dmat1 and dmat0 to get state number.
        int * state_space_size = compute_state_space(dmat0,dmat1,d.initial_Detector_energy, d.bright_state_energy);
        log << "detector_lower_bright_state_energy_window_shrink   " << detector_lower_bright_state_energy_window_shrink
            << endl;
        log<< "lower bright state number for detector 1:  "<< state_space_size[0]<<endl;
        log<<"higher bright state number for detector 1:  "<<state_space_size[1]<<endl;
        log<<"lower bright state number for detector 2:  "<<state_space_size[2] <<endl;
        log<< "higher bright state number for detector 2:  "<< state_space_size[3]<<endl;

    }
}

void full_system:: compute_detector_matrix_size_MPI_cubed( ){
    if(my_id==0){
        int i, j, k;
        int i1, i2;
        double  detector0_energy, detector1_energy;
        // ndetector0 and ndetector1 indicate current detector mode index we are calculating energy.
        vector<int> ndetector0(d.nmodes[0]);
        vector <int> ndetector1(d.nmodes[1]);
        // record size of total matrix
        int location;
        bool exist=false;

        int excited_bright_state_distance;
        int lower_bright_state_distance;
        int quanta_diff_max;

        double middle_state_energy = (d.initial_Detector_energy[0] + d.initial_Detector_energy[1])/2;
        double high_initial_state_energy = max(d.initial_Detector_energy[0] , d.initial_Detector_energy[1]);
        double low_initial_state_energy = min(d.initial_Detector_energy[0],d.initial_Detector_energy[1]);

        ndetector0[0] = -1; // this is for:  when we go into code: ndetector0[i]= ndetector0[i]+1, our first state is |000000>
        while (1) {
            label2:;  // label2 is for detector1 to jump out of while(1) loop (this is inner layer of while(1))
            detector0_energy = 0;
            for (i1 = 0; i1 < d.nmodes[0]; i1++) {  // loop through detector0
                // define the way we loop through detector0:
                ndetector0[i1] = ndetector0[i1] + 1;
                if (ndetector0[i1] <= d.nmax[0][i1]) break;
                if (ndetector0[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
                    ndetector0[d.nmodes[0] - 1] = 0;
                    goto label1;  // use goto to jump out of nested loop
                }
                ndetector0[i1] = 0;
            }
            // calculate detector 0 energy
            for (i = 0; i < d.nmodes[0]; i++) {
                detector0_energy = detector0_energy + ndetector0[i] * d.mfreq[0][i];
            }

            //--------------------------------------------------------------------------------------------
            // criteria below make sure detector 0 's energy is reasonable.
            if (detector0_energy > d.bright_state_energy[0] + d.detector_energy_window_size) {
                // detector 0's energy can not be larger than its initial energy + photon energy
                // jump to next detector state.
                k = 0;
                while (ndetector0[k] == 0) {
                    ndetector0[k] = d.nmax[0][k];
                    k++;
                    if (k >= d.nmodes[0]) {
                        break;
                    }
                }
                if (k < d.nmodes[0]) {
                    ndetector0[k] = d.nmax[0][k];
                }
                goto label2;
            }

            //------------------------------------------------------------------------------------------------
            // we make state space be a cubic not sphere, so There are more states with difference in quantum number from our original state
            quanta_diff_max = state_max_quanta_diff(ndetector0,d.initial_detector_state[0],d.nmodes[0]);
            // We constraint state we construct in a cubic cell.
            if (    quanta_diff_max > Rmax) {
                goto label2;
            }
            //--------------------------------------insert this state in detector's state.-----------------------------------------------------------
            location=find_position_for_insert_binary(vmode0, ndetector0, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                vmode0.insert(vmode0.begin() + location, ndetector0);
                dmat0.insert(dmat0.begin() + location, detector0_energy);
            }
        }
        label1:;

        ndetector1[0] = -1; // this is when we go into code: ndetector1[i] = ndetector1[i]+1. our first state is |000000>
        while (1) { // loop through detector 1
            label3:;
            detector1_energy = 0;
            for (i2 = 0; i2 < d.nmodes[1]; i2++) {
                // define the way we loop through detector1
                ndetector1[i2] = ndetector1[i2] + 1;
                if (ndetector1[i2] <= d.nmax[1][i2]) break;
                if (ndetector1[d.nmodes[1] - 1] > d.nmax[1][d.nmodes[1] - 1]) {
                    ndetector1[d.nmodes[1] - 1] = 0;
                    goto label4;
                }
                ndetector1[i2] = 0;
            }
            // calculate detector 1 energy
            for (i = 0; i < d.nmodes[1]; i++) {
                detector1_energy = detector1_energy + ndetector1[i] * d.mfreq[1][i];
            }
            // --------------------------------------------------------------
            //  criteria below make sure detector 1's energy is reasonable:
            if (detector1_energy > d.bright_state_energy[1] + d.detector_energy_window_size) {
                // initial energy is system energy.
                // detector 1 's energy can not be larger than its initial energy + photon energy
                j = 0;
                while (ndetector1[j] == 0) { // go to first mode whose n!=0;
                    ndetector1[j] = d.nmax[1][j];
                    j++;
                    if (j >= d.nmodes[1]) {
                        break;
                    }
                }
                if (j < d.nmodes[1]) {
                    ndetector1[j] = d.nmax[1][j];
                }
                goto label3;
            }

            //------------------------------------------------------------------------------------------------
            quanta_diff_max = state_max_quanta_diff(ndetector1,d.initial_detector_state[1],d.nmodes[1]);
            if (excited_bright_state_distance > Rmax and lower_bright_state_distance > Rmax
                and (quanta_diff_max > Rmax)
                    ) {
                goto label3;
            }
            location = find_position_for_insert_binary(vmode1, ndetector1, exist);
            if (!exist) {
                vmode1.insert(vmode1.begin() + location, ndetector1);
                dmat1.insert(dmat1.begin() + location, detector1_energy);
            }
        }
        label4:;
        // add function here to count state number in dmat1 and dmat0 to get state number.
        int * state_space_size = compute_state_space(dmat0,dmat1,d.initial_Detector_energy, d.bright_state_energy);
        log << "detector_lower_bright_state_energy_window_shrink   " << detector_lower_bright_state_energy_window_shrink
            << endl;
        log<< "lower bright state number for detector 1:  "<< state_space_size[0]<<endl;
        log<<"higher bright state number for detector 1:  "<<state_space_size[1]<<endl;
        log<<"lower bright state number for detector 2:  "<<state_space_size[2] <<endl;
        log<< "higher bright state number for detector 2:  "<< state_space_size[3]<<endl;

    }
}


void full_system:: construct_state_space_using_symmetry_submodule(int distance_cutoff, vector<vector<int>> * old_layer_mode_ptr,
                                                                  vector<vector<int>> * new_layer_mode_ptr,
                                                                  vector<double> * old_layer_energy_ptr,
                                                                  vector<double> * new_layer_energy_ptr
                                                                  ){
    // as input: new_ptr point to empty vector. old_ptr point to initial mode index.
    // as output: These should all be empty. Result store in vmode0, dmat0
    int i, j, k , l , m;
    int old_layer_state_num;
    int order_3_coupling_num = d.Mode_combination_list3[0].size();
    int order_4_coupling_num = d.Mode_combination_list4[0].size();
    int order_coupling_num;
    int end_coupling_num; // because coupling criteria is ordered from large to small. When coupling criteria doesn't meet, we jump other coupling
    bool change_end_coupling_num ;

    vector<vector<int>> * layer3_mode_ptr;
    vector<double> * layer3_energy_ptr;

    vector<vector<int>> * Mode_combination_list;
    vector<vector<int>> * Mode_raising_lowering_list;

    int norder;
    int mode_index;
    bool allowed;
    double energy_change;
    double sign;
    double coupling_strength;
    double coupling_nquanta;
    double criteria;

    double start_state_energy;
    double new_state_energy;

    bool exist;
    int location;

    for(i=0;i<distance_cutoff;i++){
        // Rmax is upper bound for number of layers.
        old_layer_state_num = (*old_layer_mode_ptr).size();

        for(j=0;j<old_layer_state_num;j++){
            vector<int> Start_state = (*old_layer_mode_ptr)[j];
            start_state_energy = (*old_layer_energy_ptr)[j];

            for(m=0;m<2;m++){ // m is index for 3rd order or 4th order coupling. ==0 3rd order. == 1 : 4th order
                if(m==0){
                    Mode_combination_list = & d.Mode_combination_list3[0];
                    Mode_raising_lowering_list = & d.mode_raising_lowering_list3[0];
                    norder = 3;
                    order_coupling_num = order_3_coupling_num;
                }
                else{
                    Mode_combination_list = & d.Mode_combination_list4[0];
                    Mode_raising_lowering_list = & d.mode_raising_lowering_list4[0];
                    norder = 4;
                    order_coupling_num = order_4_coupling_num;
                }

                // check 3rd or 4th order coupling
                end_coupling_num = order_coupling_num;
                change_end_coupling_num = true;

                for(k=0; k<end_coupling_num ;k++){
                    vector<int> New_state = Start_state;
                    vector<int> coupling_mode_index = (*Mode_combination_list)[k];
                    vector<int> raising_lowering_operator = (*Mode_raising_lowering_list)[k];

                    allowed = true;
                    for(l=0;l < norder ; l++ ){
                        mode_index = coupling_mode_index[l];
                        if(raising_lowering_operator[l] == 0){
                            New_state[mode_index] = New_state[mode_index] - 1;
                        }
                        else{
                            New_state[mode_index] = New_state[mode_index] + 1;
                        }
                        if(New_state[mode_index] < 0 or New_state[mode_index] > d.nmax[0][mode_index]){
                            allowed = false;
                            break;
                        }

                    }
                    // newly constructed state is not allowed, continue
                    if(!allowed){
                        continue;
                    }

                    // check V/ \Delta E criteria
                    // compute \Delta E
                    energy_change = 0;
                    for(l=0;l<norder;l++){
                        mode_index = coupling_mode_index[l];
                        if(raising_lowering_operator[l] == 0){
                            sign = -1;
                        }
                        else{
                            sign = + 1;
                        }
                        energy_change = energy_change + sign *  d.mfreq[0][mode_index];
                    }

                    // compute V
                    coupling_strength = 3050;
                    for(l=0;l<norder;l++){
                        mode_index = coupling_mode_index[l];
                        coupling_nquanta = double(Start_state[mode_index]);
                        if(raising_lowering_operator[l] == 1){
                            coupling_nquanta = coupling_nquanta + 1;
                        }

                        coupling_strength = coupling_strength * ( sqrt(double(d.mfreq[0][mode_index]) * coupling_nquanta) / 270 );
                    }

                    criteria = abs( coupling_strength / energy_change);
                    if( criteria > d.cutoff  ){
                        // check if Newly constructed state already exists.
                        location = find_position_for_insert_binary(vmode0, New_state, exist);
                        if(!exist){
                            vmode0.insert(vmode0.begin() + location, New_state);
                            new_state_energy = start_state_energy + energy_change ;
                            dmat0.insert(dmat0.begin() + location, new_state_energy);

                            // construct new generation
                            (*new_layer_mode_ptr).push_back(New_state);
                            (*new_layer_energy_ptr).push_back(new_state_energy);

                        }
                    }
                    else{

                        if(change_end_coupling_num){
                            // we will end this early.
                            change_end_coupling_num = false;
                            end_coupling_num = min(k + int(order_coupling_num / 10) , order_coupling_num);
                        }

                    }

                }



            }

        }

        // Change new_layer as old_layer so we can construct new state on top of it
        (*old_layer_mode_ptr).clear();
        (*old_layer_energy_ptr).clear();

        layer3_mode_ptr =  old_layer_mode_ptr;
        old_layer_mode_ptr = new_layer_mode_ptr;
        new_layer_mode_ptr = layer3_mode_ptr;

        layer3_energy_ptr = old_layer_energy_ptr;
        old_layer_energy_ptr = new_layer_energy_ptr;
        new_layer_energy_ptr = layer3_energy_ptr;

    }

    (*old_layer_energy_ptr).clear();
    (*old_layer_mode_ptr).clear();
}

void full_system :: include_nearby_state(vector<int> & initial_state_mode, double initial_state_energy ){
    int i , j;
    int sign;
    int location;
    bool exist;
    // put nearby state into list
    for(i=0; i< d.nmodes[0];i++){
        for(j=0;j<2;j++){
            if(j==0) {
                sign = -1;
            }
            else{
                sign = 1;
            }
            vector<int> nearby_state_mode = initial_state_mode;
            double nearby_state_energy;
            nearby_state_mode[i] = initial_state_mode[i]  + sign ;
            nearby_state_energy = initial_state_energy + sign *  d.mfreq[0][i];
            if(nearby_state_mode[i] >=0 and nearby_state_mode[i] <= d.nmax[0][i]){
                location = find_position_for_insert_binary(vmode0, nearby_state_mode, exist);
                if(!exist){
                    // if this state not exist in list. start constructing nearby state.
                    vmode0.insert(vmode0.begin() + location, nearby_state_mode);
                    dmat0.insert(dmat0.begin() + location,  nearby_state_energy);
                }

            }
        }
    }

}

void full_system:: find_resonant_state_to_initial_state( double energy_window_cutoff , vector<vector<int> > & resonant_state_mode_list , vector<double> & resonant_state_energy ){
    int i, j , k;
    int i1;
    double  detector0_energy;
    vector<int> ndetector0(d.nmodes[0]);
    int location;
    bool exist=false;
    ndetector0[0] = -1; // this is for:  when we go into code: ndetector0[i]= ndetector0[i]+1, our first state is |000000>
    double Cutoff_criteria = 0.05;
    double energy_diff ;
    while (1) {
        label2:;  // label2 is for detector1 to jump out of while(1) loop (this is inner layer of while(1))
        detector0_energy = 0;
        for (i1 = 0; i1 < d.nmodes[0]; i1++) {  // loop through detector0
            // define the way we loop through detector0:
            ndetector0[i1] = ndetector0[i1] + 1;
            if (ndetector0[i1] <= d.nmax[0][i1]) break;
            if (ndetector0[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
                ndetector0[d.nmodes[0] - 1] = 0;
                goto label1;  // use goto to jump out of nested loop
            }
            ndetector0[i1] = 0;
        }
        // calculate detector 0 energy
        for (i = 0; i < d.nmodes[0]; i++) {
            detector0_energy = detector0_energy + ndetector0[i] * d.mfreq[0][i];
        }

        //--------------------------------------------------------------------------------------------
        // criteria below make sure detector 0 's energy is reasonable.
        if (detector0_energy > d.initial_Detector_energy[0] + energy_window_cutoff ) {
            // detector 0's energy can not be larger than its initial energy + photon energy
            // jump to next detector state.
            k = 0;
            while (ndetector0[k] == 0) {
                ndetector0[k] = d.nmax[0][k];
                k++;
                if (k >= d.nmodes[0]) {
                    break;
                }
            }
            if (k < d.nmodes[0]) {
                ndetector0[k] = d.nmax[0][k];
            }
            goto label2;
        }

        if(detector0_energy > d.initial_Detector_energy[0] - energy_window_cutoff and detector0_energy < d.initial_Detector_energy[0] + energy_window_cutoff){
            // criteria here will be coupling strength / energy_diff
            energy_diff = abs ( detector0_energy - d.initial_Detector_energy[0] );

            //--------------------------------------insert this state in detector's state.-----------------------------------------------------------
            location=find_position_for_insert_binary(vmode0, ndetector0, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                vmode0.insert(vmode0.begin() + location, ndetector0);
                dmat0.insert(dmat0.begin() + location, detector0_energy);
                resonant_state_mode_list.push_back(ndetector0);
                resonant_state_energy.push_back(detector0_energy);
            }
        }
    }
    label1:;

}

void full_system:: construct_state_space_using_symmetry(){

     if(my_id == 0){
         // Number of layers is given by Rmax.
         // For each layer, we have to record newly constructed states' mode and its energy.
         // each time when constructing new layer, we use information for previous layer. Tree-type structure

         // These two list will contain mode quanta in last generation and new generation and interchange at the end of each generation
         vector<vector<int>> Layer1_mode;
         vector<vector<int>> Layer2_mode;

         vector<vector<int>> * old_layer_mode_ptr;
         vector<vector<int>> * new_layer_mode_ptr;

         // These two list will contain old generation's states' energy and new generations' states energy and will interchange at the end of geenration
         vector<double> Layer1_energy;
         vector<double> Layer2_energy;

         vector<double> * old_layer_energy_ptr;
         vector<double> * new_layer_energy_ptr;


         int i, j, k , l , m;
         double sign;

         bool exist;
         int location;
         old_layer_mode_ptr = & Layer1_mode;
         new_layer_mode_ptr = & Layer2_mode;

         old_layer_energy_ptr = & Layer1_energy;
         new_layer_energy_ptr = & Layer2_energy;

         vector<int> initial_state_mode;
         for(i=0;i<d.nmodes[0];i++){
             initial_state_mode.push_back( d.initial_detector_state[0][i] );
         }

         // first push initial state and energy into list
         vmode0.push_back(initial_state_mode);
         (*old_layer_mode_ptr).push_back(initial_state_mode);

         dmat0.push_back(d.initial_Detector_energy[0]);
         (*old_layer_energy_ptr).push_back(d.initial_Detector_energy[0]);


         // use submodule. As return. old_layer_ptr and new_layer_ptr are all empty
         construct_state_space_using_symmetry_submodule(Rmax,old_layer_mode_ptr,new_layer_mode_ptr,old_layer_energy_ptr,new_layer_energy_ptr);

         // put state with one quanta difference into list
         for(i=0; i< d.nmodes[0];i++){
             for(j=0;j<2;j++){
                 if(j==0) {
                     sign = -1;
                 }
                 else{
                     sign = 1;
                 }
                 vector<int> nearby_state_mode = initial_state_mode;
                 double nearby_state_energy;
                 nearby_state_mode[i] = initial_state_mode[i]  + sign ;
                 nearby_state_energy = d.initial_Detector_energy[0] + sign *  d.mfreq[0][i];
                 if(nearby_state_mode[i] >=0 and nearby_state_mode[i] <= d.nmax[0][i]){
                     location = find_position_for_insert_binary(vmode0, nearby_state_mode, exist);
                     if(!exist){
                         // if not exist. start constructing states from nearby state.
                         vmode0.insert(vmode0.begin() + location, nearby_state_mode);
                         dmat0.insert(dmat0.begin() + location,  nearby_state_energy);
                     }

                     (*old_layer_mode_ptr).push_back(nearby_state_mode);
                     (*old_layer_energy_ptr).push_back(nearby_state_energy);

                 }
             }
         }

         // add one more layer of neraby state
         for(i=0;i<(*old_layer_energy_ptr).size(); i++ ){
             include_nearby_state((*old_layer_mode_ptr)[i], (*old_layer_energy_ptr)[i] );
         }

         int vmode_size = vmode0.size();
         if(my_id == 0){
             cout << "number of states around initial state " << vmode_size << endl;
         }

         int distance_cutoff_for_nearby_state = 1 ;
         construct_state_space_using_symmetry_submodule(distance_cutoff_for_nearby_state, old_layer_mode_ptr, new_layer_mode_ptr,old_layer_energy_ptr,new_layer_energy_ptr);

         int vmode_size2 = vmode0.size();
         int vmode_diff = vmode_size2 - vmode_size;
         if(my_id == 0){
             cout << "number of states in layers starting from states next to initial state   " << vmode_diff << endl;
         }

//        // include state resonant to initial state into list. We also have to include state with one ladder operator difference into list
//        double energy_window = 200;
//         vector<vector<int>> resonant_state_mode_list ;
//         vector<double> resonant_state_energy_list;
//         find_resonant_state_to_initial_state(energy_window, resonant_state_mode_list, resonant_state_energy_list );
//
//         int list_len = resonant_state_energy_list.size();
//         cout <<"resonant state : " << endl;
//         for(i=0;i<list_len;i++){
//             cout << "energy:   " << resonant_state_energy_list[i] <<"      mode :    ";
//             for(j=0;j<d.nmodes[0];j++){
//                 cout << resonant_state_mode_list[i][j] << "  ";
//             }
//             cout << endl;
//         }
//         // push state with one mode quanta into list
//         for(i=0;i< list_len ; i++ ){
//             include_nearby_state(resonant_state_mode_list[i], resonant_state_energy_list[i]);
//         }

     }


}