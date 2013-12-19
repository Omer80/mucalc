#ifndef WIMP_DATA
#define WIMP_DATA

#include <iostream>
#include <string>
#include <map>
#include <iterator>
#include <vector>

using namespace std ;

class WimpData {
public :	
	WimpData(double _mass) : mass(_mass) {}
	
	~WimpData () {}
	
	// Accessors
	double Nu_e_for_x(double _log_x) {
		return Nu_e[_log_x] ;
	}
	
	double Nu_Mu_for_x(double _log_x) {
		return Nu_Mu[_log_x] ;
	}
	
	double Nu_Tau_for_x(double _log_x) {
		return Nu_Tau[_log_x] ;
	}
	
	// Printers
	void print_line_by_line() {
		vector<double>::iterator iter ;
		iter = log_x.begin() ;
		//while (iter != log_x.end()) {
			cout << *iter << " " << Nu_e_for_x(*iter) << " " ;
			cout << Nu_Mu_for_x(*iter) << " " << Nu_Tau_for_x(*iter) << endl ;
		//}
	}
	
	// Modifiers
	void add_line(double _log_x,double _Nu_e,double _Nu_Mu,double _Nu_Tau) {
		log_x.push_back(_log_x) ;
		Nu_e.insert(pair<double, double> (_log_x, _Nu_e)) ;
		Nu_Mu.insert(pair<double, double> (_log_x, _Nu_Mu)) ;
		Nu_Tau.insert(pair<double, double> (_log_x, _Nu_Tau)) ;
	}
	
	
private :
	double mass;
	vector<double> log_x ;
	map<double, double>  Nu_e;
	map<double, double>  Nu_Mu;
	map<double, double>  Nu_Tau;
	
	
} ;


#endif
