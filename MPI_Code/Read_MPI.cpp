//
// Created by phyzch on 6/23/20.
// Wrap up function to read data from input.txt and do output.
//
#include"../util.h"
#include"../system.h"

void full_system:: read_input_with_MPI(){
    if(my_id==0) {
        input.open(path + "input.txt"); // information recorded in input.txt
        if (!input.is_open()) {
            cout << "THE INFILE FAILS TO OPEN!" << endl;
            log<< "THE INFILE FAILS TO OPEN!" << endl;
            MPI_Abort(MPI_COMM_WORLD, -2);  // Abort the process and return error code -2. (input file can't open)
        }
        input >> intra_detector_coupling >> inter_detector_coupling >> Continue_Simulation >> energy_window
              >> detector_only >> Detector_Continue_Simulation >> Random_bright_state;
        input >> intra_detector_coupling_noise >> inter_detector_coupling_noise >> energy_window_size >> initial_energy
              >> noise_strength >> Rmax >> d.V_intra >> d.detector_energy_window_size
              >>detector_lower_bright_state_energy_window_shrink >> distance_cutoff_for_4_piont_corre
              >> Energy_Range_4_point_corre_function_average >> Distance_Range_4_point_corre_function_average ;
        // read time used for simulation.  delt: time step. tstart: time start to turn on coupling. tmax: maximum time for simulation.   tprint: time step to print result.
        input >> delt >> tstart >> tmax >> tprint;
        // check if input is valid
        if ( !Detector_Continue_Simulation) {
            log.open(path + "log.txt");  // log to record the error information
        }
        else{
            log.open(path + "log.txt", ios::app);
        }
        if (! log.is_open()){
            cout << "THE LOG FILE FAILS TO OPEN" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);      // Abort the process and return error code -1.  (log file can't open)
        }
        log <<"Number of Process running:  "<< num_proc <<endl;
        if ( ! Detector_Continue_Simulation) {  // start from beginning.
            output.open(path+"output.txt");
        }
        else {
            output.open(path+"output.txt", ios::app);
        }
        if( ! output.is_open()){
            cout<<"OUTPUT FILE FAILS TO OPEN."<<endl;
            log<< "OUTPUT FILE FAILS TO OPEN."<<endl;
            MPI_Abort(MPI_COMM_WORLD,-3); // Abort the process and return error code -3 (output file can't open).
        }

        if (delt <= 0 || tstart > tmax || delt > tstart) {
            cout << "Wrong time variable." << endl;
            log << "Wrong time variable." << "The simulation is cancelled." << endl;
            log.close();
            input.close();
            output.close();
            MPI_Abort(MPI_COMM_WORLD,-4);  // Abort with error code -4: Wrong time variable.
        }
        if (! Detector_Continue_Simulation) {
            output << delt << " " << tstart << " " << tmax << " " << tprint << endl;
        }
    }
    // Broadcast hyper parameter to all process.
    // function:  int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
    MPI_Bcast(&intra_detector_coupling,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
    MPI_Bcast(&inter_detector_coupling,1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Continue_Simulation,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
    MPI_Bcast(&energy_window,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
    MPI_Bcast(&detector_only,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
    MPI_Bcast(&Detector_Continue_Simulation,1,MPI_C_BOOL,0,MPI_COMM_WORLD);
    MPI_Bcast(&Random_bright_state,1,MPI_C_BOOL,0,MPI_COMM_WORLD);

    MPI_Bcast(&intra_detector_coupling_noise,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&inter_detector_coupling_noise,1, MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&energy_window_size,1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    MPI_Bcast(&initial_energy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&noise_strength,1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Rmax,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&d.V_intra,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&d.detector_energy_window_size,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&detector_lower_bright_state_energy_window_shrink,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&distance_cutoff_for_4_piont_corre,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&Energy_Range_4_point_corre_function_average,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&Distance_Range_4_point_corre_function_average,1,MPI_INT,0,MPI_COMM_WORLD);

    // Bcast delt tstart tmax tprint to other process.
    MPI_Bcast(&delt, 1, MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&tstart, 1, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
    MPI_Bcast(&tmax,1 ,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(&tprint,1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // used to rescale the matrix element amplitude.
    cf = 0.0299792458*delt * pi2;
}


