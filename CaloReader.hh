#ifndef CALOREADER_HH
#define CALOREADER_HH

#include <fstream>
#include <string>
#include <iomanip>
#include "Calorimeter.hh"

class CaloReaderException {
public :
	CaloReaderException(const char* error_message):_what(error_message) {
		
		}
	
	const char* what() const {
		return _what.c_str() ;
	}

private :
	std::string _what ;	
} ;







std::ostream& operator<<(std::ostream& os, const CaloReaderException& cre) {
	return (os << "ERROR: File open failed : " << cre.what() ) ;
}

class Caloreader {
public  :
	Caloreader(const char* input_file):
	_calo(0),
	_file(input_file) { init(input_file) ; }
	
	~Caloreader () {}

	// Accessors
	Calorimeter& calo() {
		// accessor that return a reference to the _calo pointer, that is, a reference to Calorimeter instance
		return *_calo ;
	}
	
	bool readEvent() {
		bool event = true ;
		CaloCell* _cell ;
		double _energy ;
		int _readoutID ;
		if (_file.fail()) { event = false ; }
		_calo->clear() ;
		_file >> word ;
		if (word != "BEGIN_EVENT") { event = false ;}
		_file >> word ;
		while ( word == "ENERGY" ) {
			_file >> _readoutID >> _energy ;
			_cell = _calo->findCellByID(_readoutID) ;
			_cell->setEnergy(_energy) ;
			_file >> word ;
		}
		if (word != "END_EVENT") { event = false ;}
		return event ;
	}
	
	
	

private :
	Calorimeter* _calo ;	
	ifstream _file ;
	std::string word ;
	int size_x, size_y ;
	
	// Setup functions
	void position_reader() {
		int readoutID, ix, iy ;
		CaloCell* cell;
		_file >> word ;
		while (word == "POSITION" ) {
			_file >> readoutID >> ix >> iy ;
			cell = _calo->grid().cell(ix, iy) ;
			cell->setReadoutID(readoutID) ;
			_file >> word ;			
		}	
	}
	
	void init(const char* input_file) {
				if (_file.fail()) { 
					CaloReaderException cre("error in opening the file") ;
					throw cre ; 
				} else { 
					_file >> word ;	
					if (word != "BEGIN_CALO_DEF") {
						CaloReaderException cre("No BEGIN_CALO_DEF in file") ;
						throw cre ;
					} else {
						_file >> word ;
						if (word != "SIZE" ) {
							CaloReaderException cre("No SIZE in file") ;
							throw cre ;
						} else {
							_file >> size_x >> size_y ;
							Calorimeter* calorimeter = new Calorimeter(size_x, size_y ) ;
							_calo = calorimeter ;
							position_reader();
						}
						
				if (word != "END_CALO_DEF" ) {
					CaloReaderException cre("No END_CALO_DEF in file") ;
					throw cre ;
				} 
			}								
		}
	}
	
} ;



#endif
