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
    for(i=0;i<dmat0_size;i++){
        if(dmat0[i]<initial_state_energy[0] + energy_window_size and dmat0[i] > initial_state_energy[0] - energy_window_size){
            state_number[0]++;
        }
        if(dmat0[i] < bright_state_energy[0] + energy_window_size and  dmat0[i] > bright_state_energy[0] - energy_window_size){
            state_number[1] ++;
        }
    }
    return state_number;
}


int state_distance(const vector<int> & ndetector, int * state1, int moddim){
    // compute distance between two state : ndetector and state1. Dimension of mod is moddim
    int i;
    int distance=0;
    int initial_vibrational_state_index = 1 ;
    for(i=  initial_vibrational_state_index ;i<moddim;i++){
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
        int i1;
        double  energy;

        double energy_for_vibrational_state;
        // ndetector0 and ndetector1 indicate current detector mode index we are calculating energy.
        vector<int> ndetector0(d.nmodes[0]);
        // record size of total matrix
        int location;
        bool exist=false;

        int lower_bright_state_distance;

        double middle_state_energy = (d.initial_Detector_energy[0] + d.initial_Detector_energy[1])/2;
        double high_initial_state_energy = max(d.initial_Detector_energy[0] , d.initial_Detector_energy[1]);
        double low_initial_state_energy = min(d.initial_Detector_energy[0],d.initial_Detector_energy[1]);

        // set ground electronic state energy = 0. excited electronic state energy = (\lambda_{down}^{2}^{2} - \lambda_{up}^{2}) / \hbar \omega
        // ground_electronic_state_initial_energy_shift : energy shift for ground state
        double ground_electronic_state_initial_energy_shift = 0 ;
        // excited_electronic_state_initial_energy_shift : energy shift for excited state.
        // this energy shift is due to different coupling strength in different electronic state.
        double excited_electronic_state_initial_energy_shift = pow(coupling_strength_to_mode0[0],2) / d.mfreq[0][1] - pow(coupling_strength_to_mode0[1] , 2) / d.mfreq[0][1] ;

        double crossing_region_energy = d.mfreq[0][1]/2 *
                pow( d.mfreq[0][0] /(sqrt(2) *(coupling_strength_to_mode0[0] + coupling_strength_to_mode0[1]) ) +
                sqrt(2) * coupling_strength_to_mode0[0] / d.mfreq[0][1] , 2);

        cout << "Crossing region energy " <<  crossing_region_energy << endl;
        output << "Crossing region energy " <<  crossing_region_energy << endl;
        log << "Crossing region energy " <<  crossing_region_energy << endl;

        ndetector0[0] = -1; // this is for:  when we go into code: ndetector0[i]= ndetector0[i]+1, our first state is |000000>
        while (1) {
            label2:;  // label2 is for detector1 to jump out of while(1) loop (this is inner layer of while(1))

            // ground electronic state
            ndetector0[0] = 0;
            for (i1 = 1; i1 < d.nmodes[0]; i1++) {  // loop through detector0
                // define the way we loop through detector0:
                ndetector0[i1] = ndetector0[i1] + 1;
                if (ndetector0[i1] <= d.nmax[0][i1]) break;
                if (ndetector0[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
                    ndetector0[d.nmodes[0] - 1] = 0;
                    goto label1;  // use goto to jump out of nested loop
                }
                ndetector0[i1] = 0;
            }

            energy = 0;
            energy_for_vibrational_state = 0;
            // calculate detector 0 energy
            for (i = 0; i < d.nmodes[0]; i++) {
                energy = energy + ndetector0[i] * d.mfreq[0][i];
            }
            for( i = 1; i< d.nmodes[0]; i++){
                energy_for_vibrational_state = energy_for_vibrational_state + ndetector0[i] * d.mfreq[0][i];
            }

            //--------------------------------------------------------------------------------------------
            // criteria below make sure detector 0 's energy is reasonable.
            // d.initial_Detector_energy[0] only include vibrational energy of electronic state.
            if ( energy_for_vibrational_state > d.initial_Detector_energy[0] + d.detector_energy_window_size) {
                // detector 0's energy can not be larger than its initial energy + photon energy
                // jump to next detector state.

                // start with first vibrational mode
                k = 1;
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
            // this distance only apply to vibrational state space.
            lower_bright_state_distance = state_distance(ndetector0, d.initial_detector_state[0], d.nmodes[0]);
            if ( lower_bright_state_distance > Rmax) {
                goto label2;
            }

            // shift energy by ground state energy of each electronic state

            energy = ground_electronic_state_initial_energy_shift + energy ;

            //--------------------------------------insert this state in detector's state.-----------------------------------------------------------
            location=find_position_for_insert_binary(vmode0, ndetector0, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                vmode0.insert(vmode0.begin() + location, ndetector0);
                dmat0.insert(dmat0.begin() + location, energy);
            }
        }
        label1:;

        for(i=0;i<d.nmodes[0];i++){
            ndetector0[i] = 0;
        }

        while(1){
            label4:;  // label2 is for detector1 to jump out of while(1) loop (this is inner layer of while(1))

            // excited electronic state
            ndetector0[0] = 1;
            for (i1 = 1; i1 < d.nmodes[0]; i1++) {  // loop through detector0
                // define the way we loop through detector0:
                ndetector0[i1] = ndetector0[i1] + 1;
                if (ndetector0[i1] <= d.nmax[0][i1]) break;
                if (ndetector0[d.nmodes[0] - 1] > d.nmax[0][d.nmodes[0] - 1]) {
                    ndetector0[d.nmodes[0] - 1] = 0;
                    goto label3;  // use goto to jump out of nested loop
                }
                ndetector0[i1] = 0;
            }

            energy = 0;
            energy_for_vibrational_state = 0;
            // calculate detector 0 energy
            for (i = 0; i < d.nmodes[0]; i++) {
                energy = energy + ndetector0[i] * d.mfreq[0][i];
            }
            for( i = 1; i< d.nmodes[0]; i++){
                energy_for_vibrational_state = energy_for_vibrational_state + ndetector0[i] * d.mfreq[0][i];
            }

            //--------------------------------------------------------------------------------------------
            // criteria below make sure detector 0 's energy is reasonable.
            // d.initial_Detector_energy[0] only include vibrational energy of electronic state.
            if ( energy_for_vibrational_state > d.initial_Detector_energy[1] + d.detector_energy_window_size) {
                // detector 0's energy can not be larger than its initial energy + photon energy
                // jump to next detector state.

                // start with first vibrational mode
                k = 1;
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
                goto label4;
            }

            //------------------------------------------------------------------------------------------------
            // criteria below make sure detector 1 can not be too far away from bright state and lower bright state.
            // this distance only apply to vibrational state space.
            lower_bright_state_distance = state_distance(ndetector0, d.initial_detector_state[1], d.nmodes[0]);
            if ( lower_bright_state_distance > Rmax) {
                goto label4;
            }

            // shift energy by ground state energy of each electronic state

            energy = excited_electronic_state_initial_energy_shift + energy ;

            //--------------------------------------insert this state in detector's state.-----------------------------------------------------------
            location=find_position_for_insert_binary(vmode0, ndetector0, exist);  // we check if this mode exist and the location we have to insert this state at the same time.
            if (!exist) {
                // when we push back we should consider arrange them in order. We compute location to insert in find_position_for_insert_binary() function:
                vmode0.insert(vmode0.begin() + location, ndetector0);
                dmat0.insert(dmat0.begin() + location, energy);
            }
        }
        label3:;
    }
}
