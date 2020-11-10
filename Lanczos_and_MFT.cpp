//
// Created by phyzch on 11/8/20.
//
#include "system.h"
#include "util.h"
void diagonalize( vector<double> & v, int lanc_size, int & num_eigen,
                  vector<double>  & mat, vector<int> & irow, vector<int> & icol, int matsize, int matnum );
void MFT_for_density_of_state(vector<double>  & mat, vector<int> & irow, vector<int> & icol, int matsize, int matnum, int ibright){
    // use Matrix Fluctuation Theorem to compute density of state: https://doi.org/10.1021/jp960442q
    int i, j, k;
    int bnum = 0;
    vector<int> bright_coup_pos;
    vector<int> bright_coup_strength;
    double bright_state_energy;
    int lanc_size = 5 * matsize;
    double perturbation_strength = 0.001;
    int num_eigen1 , num_eigen2;
    double itotal;
    if(my_id == 0){
        vector <double> e1;
        vector<double> e2;
        vector<double> e3;
        double emean;
        double * hnn;
        double * ipred;
        ofstream MFT_output_file("Matrix Fluctuation Theorem output.txt");
        ofstream MFT_log_file("Matrix Fluctuation Theorem log.txt");
        // we apply MFT theory to bright state of left detector.
        for(i=matsize;i<matnum;i++){
            if(irow[i] == ibright or icol[i] == ibright){
                bright_coup_pos.push_back(i);
                bright_coup_strength.push_back(mat[i]);
                bnum = bnum + 1;
            }
        }
        bright_state_energy = mat[ibright];
        MFT_log_file << mat[ibright] << ":  energy of bright state; " << endl;
        MFT_log_file << ibright << ": position of bright state;  " <<endl;
        for(i=0;i<bnum;i++){
            if (irow[bright_coup_pos[i]] == ibright) j = icol[bright_coup_pos[i]];
            if (icol[bright_coup_pos[i]] == ibright) j = irow[bright_coup_pos[i]];
            MFT_log_file << j <<"  " << bright_coup_strength[i] << endl;
        }

        MFT_log_file <<"Calculating spectrum..." << endl;
        MFT_log_file<<"Lanc_size is:  "<<lanc_size;

        for(i=0;i<bnum;i++){
            mat[bright_coup_pos[i]] = (1 + perturbation_strength) * bright_coup_strength[i];
        }
        MFT_log_file << "First Diagonalization:   "<<endl;
        diagonalize(e2,lanc_size,num_eigen2,mat,irow,icol,matsize,matnum);

        for(i=0;i<bnum;i++){
            mat[bright_coup_pos[i]] = (1 - perturbation_strength) * bright_coup_strength[i];
        }

        MFT_log_file << "Second Diagonalization   "<<endl;
        diagonalize(e1,lanc_size,num_eigen1,mat,irow,icol,matsize,matnum);

        if(num_eigen1 != matsize or num_eigen2 != matsize){
            cout <<" Only total number of: "<< num_eigen1 <<"   and    " << num_eigen2<<"   out of  "<<matsize <<"  Found."<<endl;
            MFT_log_file <<" Only total number of: "<< num_eigen1 <<"   and    " << num_eigen2<<"   out of  "<<matsize <<"  Found."<<endl;
        }

        i=0;
        hnn = new double[min(num_eigen1,num_eigen2)];
        ipred = new double[min(num_eigen1, num_eigen2)];
        while(i<num_eigen1 and i<num_eigen2){
            emean = (e2[i] + e1[i])/2;
            hnn[i] = (e2[i] - e1[i])/(2*perturbation_strength); // ( En(1+\Delta) - En(1-\Delta) )/ (2*Delta)
            ipred[i] = 1/(2*(emean - mat[ibright])) * hnn[i];  // |<n|i>|^{2}
            itotal = itotal + ipred[i];  // sum of |<n|i>|^{2}
        }

        MFT_output_file<<"Total intensity (should be 1)  "<< itotal <<endl;
        for(i=0;i<min(num_eigen1,num_eigen2); i++){
            if(ipred[i] > 0.001){
                MFT_output_file << e1[i] <<"  "<<e2[i] << "  " << ipred[i] << endl;
            }
        }
        delete [] hnn;
        delete [] ipred;
    }

}

void diagonalize( vector<double> & v, int lanc_size, int & num_eigen ,
                  vector<double>  & mat, vector<int> & irow, vector<int> & icol, int matsize, int matnum ){
    // v is array contain eigenvalue.
    // lanc_size is size of Krylov space
    // num_eigen is total number of eigenvalue we obtain

}