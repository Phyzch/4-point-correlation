//
// Created by phyzch on 7/22/20.
//
#include"../util.h"
#include"../system.h"

// read parameters for system matrix set up
void system::read_MPI(ifstream &input, ofstream &output, ofstream &log) {
    // read energy level tle , initial wavefunction xtl, ytl
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    if(my_id==0) {
        // only process 0 will do I/O here.
        // read number of detectors. typically 1 or 2.
        input >> tlnum;
        if (tlnum != 1 && tlnum != 2) {
            log << "TLM NUMBER NOT SUPPORTED" << endl;
            input.close();
            log.close();
            output.close();
            exit(-5 ); // tlnum is not right.
        }
        if (!Detector_Continue_Simulation) {
            output << "system  " << tlnum << " ";
        }
    }
    MPI_Bcast(&tlnum, 1, MPI_INT,0,MPI_COMM_WORLD);
    tlmatsize = pow(2, tlnum);  // system wave function array size.

    initialize_energy_level(input,output);

    initialize_wavefunction(input,output);
    // do not have to parallelize initializes_state_energy
    initialize_state_energy();

    tlmatnum = tlmatsize;
};


void system::initialize_energy_level(ifstream & input, ofstream & output){
    // initialize energy of state
    int i;
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if(my_id==0) {
        for (i = 0; i < tlnum; i++) {
            input >> tle[i];
            if (!Detector_Continue_Simulation) {
                output << tle[i] << " ";
            }
        }
        output << endl;
    }
    // broadcast tle to all other process. (Variable in class system is all very small, do not have to data decomposition).
    MPI_Bcast(tle,tlnum,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void system::initialize_wavefunction(ifstream & input, ofstream & output){
    //initialize wavefunction of state & normalize it
    int i;
    double norm = 0;
    int my_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    if(my_id==0) {
        for (i = 0; i < tlmatsize; i++) {
            input >> xtl[i] >> ytl[i];
            if (!Detector_Continue_Simulation) {
                output << xtl[i] << " " << ytl[i] << endl;
            }
            norm = norm + pow(xtl[i], 2) + pow(ytl[i], 2);
        }
        norm = 1 / sqrt(norm);
        for (i = 0; i < tlmatsize; i++) {
            xtl[i] = xtl[i] * norm;
            ytl[i] = ytl[i] * norm;
        }
    }
    MPI_Bcast(xtl,tlmatsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(ytl,tlmatsize,MPI_DOUBLE,0,MPI_COMM_WORLD);
}

void system::initialize_state_energy(){
    // initialize energy of system state.
    int i,j;
    int * marray = new int[tlnum];  // intermediate variable used to initialize matrix element
    for (i = 0; i < tlnum; i++) {
        marray[i] = 0;
    }
    for (i = 0; i < tlmatsize; i++) {
        tlmat[i] = 0;
        for (j = 0; j < tlnum; j++) {
            tlmat[i] = tlmat[i] + marray[j] * tle[j];
        }
        tlirow[i] = i;
        tlicol[i] = i;
        // marray will be [0,0],[1,0],[0,1],[1,1] for tlnum=2
        for (j = 0; j < tlnum; j++) {
            marray[j] = marray[j] + 1;
            if (marray[j] <= 1)  break;
            if (marray[tlnum - 1] > 1) break;
            marray[j] = 0;
        }
    }
    delete[] marray;
}