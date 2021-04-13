//
// Created by phyzch on 4/12/21.
// This file include function about how to utilize symmetry of vibrational mode to constrain states constructed around initial states
//
# include "util.h"
# include "system.h"

// Here A1 == 0, A2 == 1, B1 == 2, B2 == 3
int Symmetry_Table [4][4] = {
        {0, 1 ,2 ,3},
        {1, 0, 3, 2},
        {2, 3, 0, 1},
        {3, 2, 1, 0}
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
            }

        }
        label1:;

    }

    int Mode_combination_list_size = Mode_combination_list.size();
    if(my_id == 0){
        cout << "order coupling :  " << norder <<"  allowed combination by symmetry and cutoff criteria:  " << Mode_combination_list_size << endl;
    }

}

void detector::construct_Mode_combination_list() {
    int norder ;
    int i;

    double nquanta_mean_estimate = nmax[0][0]; // used to compute coupling strength for mode combination.
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