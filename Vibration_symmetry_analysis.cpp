//
// Created by phyzch on 4/12/21.
// This file include function about how to utilize symmetry of vibrational mode to constrain states constructed around initial states
//
# include "util.h"
# include "system.h"
void order_coupling_term(vector<vector<int>> & Mode_combination_list, vector<vector<int>> & mode_raising_lowering_list, vector<double> & criteria_list );
// Here A1 == 0, A2 == 1, B1 == 2, B2 == 3
int Symmetry_Table [4][4] = {
        {0, 1 ,2 ,3},
        {1, 0, 3, 2},
        {2, 3, 0, 1},
        {3, 2, 1, 0}
};

struct sort_coupling_element{
    double criteria;
    vector<int> mode_combination;
    vector<int> mode_raising_lowering;
    sort_coupling_element(double criteria1, const vector<int> & mode_combination_1, const vector<int> & mode_raising_lowering_1){
        criteria = criteria1;
        mode_combination = mode_combination_1;
        mode_raising_lowering = mode_raising_lowering_1;
    }
};


void construct_allowed_mode_combination( int norder, vector<vector<int>> & Mode_combination_list,
                                         vector<vector<int>> & mode_raising_lowering_list ,
                                         int nmodes, const int * Mode_symmetry ,
                                         double nquanta_estimate, const double * mfreq , double cutoff
                                         ){
    /*
     *  norder: order of anharmonicity , ==3 or == 4
     *  nmodes: number of modes for molecules
     *  Mode_symmetry: symmetry of modes. size: [nmodes]
     *  Mode_combination_list. share by norder ==3 and norder == 4. combination of raising and lowering operator allowed by symmetry and also pass V/\Delta E criteria
     * mode_raising_lowering_list : symbol for raising and lowering operator. 0 for lowering operator, 1 for raising operator

     *  (this is computationally easy , so each process can call this function and compute their own.)
     */
    int i, j, k, l, m ;
    int index1, index2, index_combine;
    if(norder != 3 and norder != 4 ){
        if(my_id == 0){
            cout <<"Norder for computing mode combination is wrong. " << endl;
        }
    }

    // Large list without considering symmetry.
    vector<vector<int>> candidate_Mode_combination_list ;
    if(norder ==3 ){
        for(i=0;i<nmodes; i++){
            for(j=i ; j< nmodes; j++ ){
                for(k = j; k< nmodes; k++){
                    vector<int> candidate_mode_combination { i, j, k};
                    candidate_Mode_combination_list.push_back(candidate_mode_combination);
                }
            }
        }
    }

    if(norder == 4){
        for(i = 0;i < nmodes;i++){
            for(j = i;j < nmodes;j++){
                for(k = j;k < nmodes; k++){
                    for( l = k ; l<nmodes; l++ ){
                        vector<int> candidate_mode_combination {i,j,k,l} ;
                        candidate_Mode_combination_list.push_back(candidate_mode_combination);
                    }
                }
            }
        }
    }
    unsigned long candidate_list_size = candidate_Mode_combination_list.size();
    if(my_id == 0){
        cout << "order of coupling:  " << norder <<"  candidate combination:   " << candidate_list_size << endl;
    }
    vector<vector<int>> symmetry_allowed_candidate_Mode_combination_list ;

    for(i=0;i< candidate_list_size; i++){
        vector<int> candidate_mode_combination = candidate_Mode_combination_list[i];
        index1 = Mode_symmetry [ candidate_mode_combination[0] ];
        index2 = Mode_symmetry[ candidate_mode_combination[1] ] ;
        for(j=0; j < norder -2 ; j++){
            index_combine = Symmetry_Table [index1][index2]; //  Read symmetry when combining two modes
            index1 = index_combine;
            index2 = Mode_symmetry[ candidate_mode_combination[j + 2] ];
        }
        index_combine = Symmetry_Table[index1][index2];
        // now index_combine == product of symmetry of modes.
        // For example, norder ==3 , 3 mode have symmetry:  A1, A2, B1, then index_combine give 3, which means B2.
        if(index_combine == 0){
            symmetry_allowed_candidate_Mode_combination_list.push_back(candidate_mode_combination);
        }

    }
    unsigned long symmetry_allowed_list_size = symmetry_allowed_candidate_Mode_combination_list.size();
    if(my_id == 0){
        cout <<" order of coupling : " << norder <<"  allowed combination by symmetry  " << symmetry_allowed_list_size << endl;
    }

    // Now sift allowed coupling using criteria V / \Delta E.
    // V: 3050 * \prod (omega^{0.5} * nquanta^{0.5} / 270)^{nk}
    double energy_change ;
    double coupling_strength;
    int sign ;
    bool continue_bool ;
    double criteria;

    vector<double> criteria_list ;
    for(i=0;i<symmetry_allowed_list_size;i++){
        vector<int> Mode_combination = symmetry_allowed_candidate_Mode_combination_list[i];

        vector<int> raising_lowering_tuple ;
        for(j=0; j < norder ; j++){
            raising_lowering_tuple.push_back(0); // start with all lowering operator.
        }
        raising_lowering_tuple[0] = -1 ;

        while(1){
            // go through all combination of raising and lowering operator. [000] -> [111] or [0000]->[1111]
            for(j=0;j < norder ; j++){
                raising_lowering_tuple[j]  = raising_lowering_tuple[j] + 1;
                if( raising_lowering_tuple[j] < 2){
                    break;
                }
                if (raising_lowering_tuple[norder -1 ] > 1 ){
                    raising_lowering_tuple[norder - 1] = 0;
                    goto label1;
                }
                raising_lowering_tuple[j] = 0;
            }

            // We have to eliminate repetition: for example: mode tuple (1 1 1). it can generate
            // a_{1} a_{1}^{*} a_{1} and a_{1} a_{1} a_{1}^{*}. But they are the same.
            // if mode_combination[i][j] == mode_combination[i][j+1] and raising_lowering_tuple[j] > raising_and_lowering_tuple[j+1].  continue
            continue_bool = false;
            for(j=0;j<norder - 1; j++){
                if(Mode_combination[j] == Mode_combination[j+1]
                and  raising_lowering_tuple[j] > raising_lowering_tuple[j+1] ){
                    // then we have case not follow normal order: have the form : a_{i}^{+} a_{i}
                    // we don't include this term in our Mode_combination_list
                    continue_bool = true;
                }
            }
            if( continue_bool ){
                continue;
            }

            // compute \Delta E
            energy_change = 0;
            for(j=0; j<norder; j++){
                if(raising_lowering_tuple[j] == 0){
                    sign = -1 ;
                }
                else{
                    sign = 1;
                }
                energy_change = energy_change + sign  * mfreq[Mode_combination[j]];
            }
            // compute V:
            coupling_strength = 3050;
            for(j=0;j<norder; j++){
                coupling_strength = coupling_strength * (   sqrt(mfreq[Mode_combination[j]]  * nquanta_estimate ) / 270
                        );
            }

            criteria = abs( double(coupling_strength) / energy_change ) ;
            if( criteria > cutoff ){
                Mode_combination_list.push_back(Mode_combination);
                mode_raising_lowering_list.push_back(raising_lowering_tuple);
                criteria_list.push_back(criteria);
            }

        }
        label1:;

    }

    int Mode_combination_list_size = Mode_combination_list.size();
    if(my_id == 0){
        cout << "order coupling :  " << norder <<"  allowed combination by symmetry and cutoff criteria:  " << Mode_combination_list_size << endl;
    }

    // In each order, we should order the coupling respective to their resonance. (criteria)
    order_coupling_term(Mode_combination_list,mode_raising_lowering_list, criteria_list);

}

vector<sort_coupling_element> merge_sort(const vector<sort_coupling_element> & v1, const vector<sort_coupling_element> & v2){
    int size1= v1.size();
    int size2= v2.size();
    if(size1 == 0) return v2;
    if(size2 == 0) return v1;
    int v1_index = 0;
    int v2_index = 0;
    double criteria_difference;
    int i;
    vector <sort_coupling_element> v3;
    while(v1_index<size1 and v2_index<size2){
        criteria_difference = v1[v1_index].criteria - v2[v2_index].criteria;
        if(criteria_difference > 0){
            v3.push_back(v1[v1_index]);
            v1_index++;
        }
        else if (criteria_difference < 0){
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

// Sort coupling according to their criteria
void order_coupling_term(vector<vector<int>> & Mode_combination_list, vector<vector<int>> & mode_raising_lowering_list, vector<double> & criteria_list ){
    // order coupling term by its energy change
    // use merge sort algorithm
    int i,j;
    int coupling_num = Mode_combination_list.size();
    vector<sort_coupling_element> list_to_sort;
    // construct sort_element for merge_sort
    for(i=0;i<coupling_num;i++){
        sort_coupling_element element(criteria_list[i],Mode_combination_list[i],mode_raising_lowering_list[i]);
        list_to_sort.push_back(element);
    }

    // start merge sorting
    vector <vector<sort_coupling_element>> List_for_list_old;
    vector<vector<sort_coupling_element>> * old_ptr = & List_for_list_old;
    vector<vector<sort_coupling_element>> List_for_list_new;
    vector<vector<sort_coupling_element>> * new_ptr = & List_for_list_new;
    vector<vector<sort_coupling_element>> * list_ptr_3;
    vector<sort_coupling_element> v3;
    int list_size = coupling_num;
    for(i=0;i<coupling_num;i++){
        vector<sort_coupling_element> small_list;
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
    Mode_combination_list.clear();
    mode_raising_lowering_list.clear();
    criteria_list.clear();
    list_to_sort = (*old_ptr)[0];  // sorting result store in *old_ptr
    for(i=0;i<coupling_num;i++){
        Mode_combination_list.push_back(list_to_sort[i].mode_combination);
        mode_raising_lowering_list.push_back(list_to_sort[i].mode_raising_lowering);
        criteria_list.push_back(list_to_sort[i].criteria);
    }
}



void detector::construct_Mode_combination_list() {
    int norder ;
    int i;

    double nquanta_mean_estimate = 3; // used to compute coupling strength for mode combination.
    if(my_id == 0){
        cout << " nquanta for estimating allowed operator :  " << nquanta_mean_estimate << endl;
    }
    norder = 3 ;
    for(i=0;i<stlnum; i ++){
        construct_allowed_mode_combination(norder, Mode_combination_list3[i], mode_raising_lowering_list3[i], nmodes[i], Mode_Symmetry[i], nquanta_mean_estimate, mfreq[i] , cutoff);
    }

    norder = 4;
    for(i=0;i<stlnum;i++){
        construct_allowed_mode_combination(norder, Mode_combination_list4[i],mode_raising_lowering_list4[i] ,nmodes[i], Mode_Symmetry[i] , nquanta_mean_estimate, mfreq[i], cutoff);
    }
}