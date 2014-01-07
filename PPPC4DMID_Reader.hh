#ifndef PPPC4DMID_READER
#define PPPC4DMID_READER

#include <fstream>
#include <string>
#include <iomanip>
#include "WimpData.hh"


class ReaderException {
public :
	ReaderException(const char* error_message):_what(error_message) {
		}
	
	const char* what() const {
		return _what.c_str() ;
	}

private :
	std::string _what ;	
} ;







std::ostream& operator<<(std::ostream& os, const ReaderException& cre) {
	return (os << "ERROR: File open failed : " << cre.what() ) ;
}

class DmDataReader {
public  :
	DmDataReader(const char* input_file):
	_file(input_file) { init() ;  }
	
	~DmDataReader () {}

	// Printers
	void print_buf(ostream& os = cout) {
		os << "buf is " << buf << endl ;
	}
	
	void print_data_for_wimp_mass(int _mass) {
		_wimp[_mass]->print_line_by_line() ;
	}
	
	void print_f_gamma_partial_for_wimp_mass(int _mass) {
		_wimp[_mass]->set_f_gamma_partial() ;
		_wimp[_mass]->print_f_gamma_partial() ;
	}
	
	
	
	

private :
	//WimpData* _wimp;	
	ifstream _file ;
	std::string word ;
	string buf ;
	vector<string> def;
	vector<WimpData> _mDM ;
	map<int, WimpData*> _wimp ;
	int _mass ;
	
	void init() {
		if (_file.fail()) { 
					ReaderException cre("error in opening the file") ;
					throw cre ; 
		} else {
			cout << "File was opened" << endl ;
			read_first_line() ;
			//_file >> _mass ;
			//read_mass() ;
			//read_mass() ;
			read_data();
		}	
	}
	
	void read_first_line () {
		
		int raw_length = 30 ;
		for (int i=0 ; i < raw_length ; i++) {
			_file >> buf ;
			def.push_back(buf) ;
			cout << def[i] << " ";
		}
		cout << endl << "So far for definitions." << endl ;
	}
	
	void read_data () {
		cout << "Starting reading data from file" ;
		int j = 0 ;
		_file >> _mass ;
		//for (int j = 0 ; j < 10 ; j++)
		while(!_file.fail()) {
			if (_file.eof()) break ;
			// Loop over the same masses and enter the values to the same WimpData object till the mass 
			// changes
			read_mass() ;
			
			j++;
			
		}
		cout << endl << "Finished reading data from file" << endl ;
		
	}
	
	void read_mass() {
		cout << "." ;
		vector<double> value(29);
		double _value ;
		
		int _next_mass ;
		
		WimpData* mDM = new WimpData(_mass) ;
		_next_mass = _mass ;
		//cout << "Reading data for mass " << _mass << endl;
		int j=0;
		do {
			for (int i=0 ; i < 29 ; i++) {
				_file >> _value ;
				j++;
				value[i]=_value;
				//cout << value[i] << " ";
			}
			//cout << endl ;
			mDM->add_line(value[0] ,value[23],value[24],value[25]);
			if (_file.eof()) break ;
			_file >> _next_mass ;
			if (_file.fail()) break ;
			//cout << "Next mass is " << _next_mass << endl ;
			j++;
		} while (_next_mass == _mass) ;
		_wimp.insert(pair<int, WimpData*>(_mass, mDM));
		//cout << "Finished reading data for mass " << _mass << endl ;
		_mass = _next_mass ;
		
	}
	
	
} ;



#endif
