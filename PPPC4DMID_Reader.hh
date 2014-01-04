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
	
	
	

private :
	//WimpData* _wimp;	
	ifstream _file ;
	std::string word ;
	string buf ;
	vector<string> def;
	vector<WimpData> _mDM ;
	map<double, WimpData*> _wimp ;
	
	void init() {
		if (_file.fail()) { 
					ReaderException cre("error in opening the file") ;
					throw cre ; 
		} else {
			cout << "File was opened" << endl ;
			read_first_line() ;
			read_data() ;
		}	
	}
	
	void read_first_line () {
		
		int raw_length = 30 ;
		for (int i=0 ; i < raw_length ; i++) {
			_file >> buf ;
			def.push_back(buf) ;
			cout << def[i] << " ";
		}
	}
	
	void read_data () {
		int j = 0 ;
		//for (int j = 0 ; j < 10 ; j++)
		while(!_file.eof()) {
			if (_file.eof()) break ;
			
			// Loop over the same masses and enter the values to the same WimpData object till the mass 
			// changes
			read_mass() ;
			
			j++;
			//cout << "For mass " << _mass << endl ;
			//cout << setw(10) << def[1] << setw(10) << def[24] << setw(10) << def[25] << setw(10) <<def[26] << endl ;
			//cout << setw(10) << value[1] << setw(10) <<value[24] << setw(10) <<value[25] << setw(10) <<value[26] << endl ;			
		}
		
	}
	
	void read_mass() {
		vector<double> value(29);
		double _value ;
		double _mass ;
		double _next_mass ;
		_file >> _mass ;
		WimpData* mDM = new WimpData(_mass) ;
		_next_mass = _mass ;
		do {
			if (_file.eof()) break ;
			for (int i=0 ; i < 29 ; i++) {
				_file >> _value ;
				value[i]=_value;
			}
			mDM->add_line(value[1] ,value[24],value[25],value[26]);
			_file >> _next_mass ;
		} while (_next_mass == _mass) ;
		_wimp.insert(pair<double, WimpData*>(_mass, mDM));
		
	}
	
	
} ;



#endif
