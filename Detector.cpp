#include"system.h"
#include"util.h"
using namespace std;
// add comment to test git.

// if we are going to use energy_window, we don't have to pay attention to detdim and dmatdim, because we will calculate these number in our own functioin
// however, if we are not going to add energy window, we have to make sure these values are right before begin our simulation.

detector::detector(){
    ; // do nothing here, I write this code because I don't know what will compiler do if I don't initialize class.
}

detector::~detector(){
    int i;
    int sampling_state_index_size = sampling_state_index.size();

    delete [] dmatsize;
    delete [] dmatnum;
    delete [] doffnum;
    delete [] deln;
    delete [] nbar;
    delete [] mfreq;
    delete [] modcoup;
    delete [] premodcoup;

    delete [] proptime;
//    delete [] bright_state;
    delete [] initial_detector_state;
    delete [] initial_Detector_energy;
//    delete [] bright_state_energy;

    for(i=0;i<1;i++){
        delete [] dmatsize_each_process[i];
        delete [] dmatsize_offset_each_process[i];
        delete [] doffnum_each_process[i];
        delete [] dmatnum_each_process[i];
        delete [] dmat_offset_each_process[i];
        delete [] total_dmat[i];
        delete [] total_dirow[i];
        delete [] total_dicol[i];
    }
   for(i=0;i<sampling_state_index_size;i++){
            delete [] xd_all[i];
            delete [] yd_all[i];
    }
    delete [] dmat;
    delete [] dmatsize_each_process;
    delete [] dmatsize_offset_each_process;
    delete [] doffnum_each_process;
    delete [] dmatnum_each_process;
    delete [] dmat_offset_each_process;
    delete [] total_dmat;
    delete [] total_dirow;
    delete [] total_dicol;

    delete [] xd_all;
    delete [] yd_all;
}
