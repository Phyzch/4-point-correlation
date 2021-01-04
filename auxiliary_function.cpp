//
// Created by phyzch on 4/14/20.
//Used for putting non essential function in same file.
//
#include"system.h"
#include"util.h"
using namespace std;

// code used to reform the final output if we re-start our simulation
void detector::replace_4_point_corr_second_line(double detector_tprint){
    // This code is quite dumb... I find if I continue my simulation, I have to correct the endtime I give. so I have to rewrite whole output.txt for the first line.
    ifstream fin;
    ofstream temp;
    int x;
    if(my_id==0) {
        string old_path = path + "4 point correlation.txt";
        fin.open(old_path);
        string new_path = path + "4 point correlation_new.txt";
        temp.open(new_path);
        string line;
        getline(fin, line); // read first line
        temp << line <<endl;
        getline(fin,line); // read second line
        temp << "total time: "  << " " << proptime[0] << " " << detector_tprint << endl;
        while (std::getline(fin, line)) {
            if (line != "") {
                temp << line << endl;
            }
        }
        temp.close();
        fin.close();
        if (remove(old_path.c_str()) == 0) {
            rename(new_path.c_str(), old_path.c_str());
        } else {
            cout << "Remove file failed" << endl;
            cin >> x;
        }
    }

    if(my_id==0) {
        string old_path = path + "4 point correlation average over states.txt";
        fin.open(old_path);
        string new_path = path + "4 point correlation_new.txt";
        temp.open(new_path);
        string line;
        getline(fin, line); // read first line
        temp << line <<endl;
        getline(fin,line); // read second line
        temp << "total time: "  << " " << proptime[0] << " " << detector_tprint << endl;
        while (std::getline(fin, line)) {
            if (line != "") {
                temp << line << endl;
            }
        }
        temp.close();
        fin.close();
        if (remove(old_path.c_str()) == 0) {
            rename(new_path.c_str(), old_path.c_str());
        } else {
            cout << "Remove file failed" << endl;
            cin >> x;
        }
    }
}

// check some dimension parameter
void full_system::dimension_check() {
    double errflag = 0;
    if (!energy_window) {
        if (d.dmatnum[0] > d.dmatdim || d.dmatnum[1] > d.dmatdim) errflag = errflag + 1;
        if (d.dmatsize[0] > d.detdim || d.dmatsize[1] > d.detdim) errflag = errflag + 2;
    }
    if (errflag != 0) {
        log << " Dimension Problem: Error flag=" << errflag << endl;
        cout<<"Dimension error, Error flag="<<errflag;
        exit(-1);
    }

    if (s.tlnum == 1) {
        if (! Detector_Continue_Simulation) {
            output << "Global Matrix: 2*" << d.dmatsize[0] << " = " << matsize << endl;
        }
    }
    else if (s.tlnum == 2) {
        if (! Detector_Continue_Simulation ) {
            if (! energy_window) {
                output << "Global Matrix : 4*" << d.dmatsize[0] << " * " << d.dmatsize[1] << " = " << matsize << endl;
            }
            else {
                output << "Global Matrix: " << total_matsize << endl;
            }
        }
    }
    if (!Detector_Continue_Simulation) {
        output << "off-diagonal matrix number  " << total_offnum << endl;
        output << "Whole matrix element number  " << total_matnum << endl;
        log << "off-diagonal matrix number  " << total_offnum << endl;
        log << "Whole matrix element number  " << total_matnum << endl;
    }
};


