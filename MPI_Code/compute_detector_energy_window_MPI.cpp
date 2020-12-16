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


void full_system:: compute_detector_matrix_size_MPI( ){
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

            // criteria for energy window around bright_state and lower bright state for detector 0
            if ((detector0_energy > d.initial_Detector_energy[0] - d.detector_energy_window_size / detector_lower_bright_state_energy_window_shrink and
                 detector0_energy < d.initial_Detector_energy[0] + d.detector_energy_window_size / detector_lower_bright_state_energy_window_shrink)
//                or (detector0_energy > d.bright_state_energy[0] - d.detector_energy_window_size and
//                    detector0_energy < d.bright_state_energy[0] + d.detector_energy_window_size)
                    )
//            if(detector0_energy>=d.initial_Detector_energy[0] and detector0_energy<= d.bright_state_energy[0] + energy_window_size )
            { ;
            } else {
                goto label2;
            }
            //------------------------------------------------------------------------------------------------
            // criteria below make sure detector 1 can not be too far away from bright state and lower bright state.
            excited_bright_state_distance = state_distance(ndetector0, d.bright_state[0], d.nmodes[0]);
            lower_bright_state_distance = state_distance(ndetector0, d.initial_detector_state[0], d.nmodes[0]);
            if (excited_bright_state_distance > Rmax and lower_bright_state_distance > Rmax
                 and ( detector0_energy > high_initial_state_energy or detector0_energy<low_initial_state_energy )
            ) {
                goto label2;
            }
            //-------------------------------------------------------------------------------------------------
            location=find_position_for_insert_binary(vmode0, ndetector0, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                // this detector state has not been recorded yet. We should record them.
                // copy detector 0 mode information into temporary_vector and store this vector in vmode0
                vector<int> temporary_vector(d.nmodes[0]);
                copy_array(temporary_vector, ndetector0, d.nmodes[0]);
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                vmode0.insert(vmode0.begin() + location, temporary_vector);
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

            // criteria for energy window around bright_state and lower bright state for detector 1
            if ((detector1_energy > d.initial_Detector_energy[1] - d.detector_energy_window_size / detector_lower_bright_state_energy_window_shrink and
                 detector1_energy < d.initial_Detector_energy[1] + d.detector_energy_window_size / detector_lower_bright_state_energy_window_shrink)
//                or (detector1_energy > d.bright_state_energy[1] - d.detector_energy_window_size and
//                    detector1_energy < d.bright_state_energy[1] + d.detector_energy_window_size)
                    )
//                if(detector1_energy>=d.initial_Detector_energy[1]  and detector1_energy<= d.bright_state_energy[1]+energy_window_size)
            {  // criteria here means we only consider detector state whose energy is within small energy window
                ;
            } else {
                goto label3;
            }
            //------------------------------------------------------------------------------------------------
            // criteria below make sure detector 1 can not be too far away from bright state and lower bright state.
            excited_bright_state_distance = state_distance(ndetector1, d.bright_state[1], d.nmodes[1]);
            lower_bright_state_distance = state_distance(ndetector1, d.initial_detector_state[1], d.nmodes[1]);
            if (excited_bright_state_distance > Rmax and lower_bright_state_distance > Rmax
             and (detector1_energy< low_initial_state_energy or detector1_energy > high_initial_state_energy )
            ) {
                goto label3;
            }
            location = find_position_for_insert_binary(vmode1, ndetector1, exist);
            if (!exist) {
                vector<int> temporary_vector(d.nmodes[1]);
                copy_array(temporary_vector, ndetector1, d.nmodes[1]);
                vmode1.insert(vmode1.begin() + location, temporary_vector);
                dmat1.insert(dmat1.begin() + location, detector1_energy);
            }
        }
        label4:;
        // add function here to count state number in dmat1 and dmat0 to get state number.
        int * state_space_size = compute_state_space(dmat0,dmat1,d.initial_Detector_energy, d.bright_state_energy);
        log<< "lower bright state number for detector 1:  "<< state_space_size[0]<<endl;
        log<<"higher bright state number for detector 1:  "<<state_space_size[1]<<endl;
        log<<"lower bright state number for detector 2:  "<<state_space_size[2] <<endl;
        log<< "higher bright state number for detector 2:  "<< state_space_size[3]<<endl;

    }
}

void full_system:: compute_detector_matrix_size_MPI_new( ){
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
//            if (detector0_energy > d.bright_state_energy[0] + d.detector_energy_window_size) {
//                // detector 0's energy can not be larger than its initial energy + photon energy
//                // jump to next detector state.
//                k = 0;
//                while (ndetector0[k] == 0) {
//                    ndetector0[k] = d.nmax[0][k];
//                    k++;
//                    if (k >= d.nmodes[0]) {
//                        break;
//                    }
//                }
//                if (k < d.nmodes[0]) {
//                    ndetector0[k] = d.nmax[0][k];
//                }
//                goto label2;
//            }

            // criteria for energy window around bright_state and lower bright state for detector 0
//            if ((detector0_energy > low_initial_state_energy - d.detector_energy_window_size / detector_lower_bright_state_energy_window_shrink and
//                 detector0_energy < high_initial_state_energy + d.detector_energy_window_size / detector_lower_bright_state_energy_window_shrink)
//                    )
//                ;
//            else {
//                goto label2;
//            }
            //------------------------------------------------------------------------------------------------
            // criteria below make sure detector 1 can not be too far away from bright state and lower bright state.
//            excited_bright_state_distance = max(state_distance(ndetector0, d.bright_state[0], d.nmodes[0]),
//                                                state_distance(ndetector0, d.bright_state[1], d.nmodes[1]));
//            lower_bright_state_distance = max(state_distance(ndetector0, d.initial_detector_state[0], d.nmodes[0]),
//                                              state_distance(ndetector0, d.initial_detector_state[1], d.nmodes[1]));

            // we make state space be a cubic not sphere, so There are more states with difference in quantum number from our original state
               quanta_diff_max = state_max_quanta_diff(ndetector0,d.initial_detector_state[0],d.nmodes[0]);

            // We constraint state we construct in a cubic cell.
            if (    quanta_diff_max > Rmax
                    ) {
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