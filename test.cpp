// Testing PPPC4DMID_Reader.hh
// Omer Tzuk

#include <iostream>
#include "PPPC4DMID_Reader.hh"

using namespace std ;

int main () {
	try {DmDataReader r("data/AtProduction_neutrinos_e.dat") ;
	//r.print_data_for_wimp_mass(5);
	r.print_f_gamma_partial_for_wimp_mass(1000) ;
	}
	catch (ReaderException cre) { 
		cout << cre << endl ; 
	}
	return 0 ;
}
