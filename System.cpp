#include"util.h"
#include"system.h"
using namespace std;

#define pi2 3.141592653589793*2 

//// allocate space for matrix and array
system::system() {
    // xtl: real part of system wave function
    // ytl: image part of system wave function
    // tle: energy level of system's eigen state, also diagonal part of our Hamiltonian matrix
    // tlrho: density matrix of system
    // tlirow: row index for matrix element in tlmat
    // tlicol: column index for matrix element in tlmat
	xtl = new double[int(pow(2, tldim))];
	ytl = new double[int(pow(2, tldim))];
	tle = new double[tldim];
	tlmat = new double[int(pow(2, tldim))];
	tlrho = new double[int(pow(2, tldim)), int(pow(2, tldim))];
	tlirow = new int[int(pow(2, tldim))];
	tlicol = new int[int(pow(2, tldim))];
}



system::~system(){
    delete [] xtl;
    delete [] ytl;
    delete [] tle;
    delete [] tlmat;
    delete [] tlrho;
    delete [] tlirow;
    delete [] tlicol;
}
