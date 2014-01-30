import numpy as np
from astropy import units as u
from astropy.units import imperial
from astropy import cosmology
import scipy.integrate
from scipy.interpolate import interp1d
import const
import math

cosmology.core.set_current(cosmology.Planck13)

#print cosmology.get_current()
#print cosmology.H(10**6)

#Physical constants 
const.Sigma_T = 6.65246e-25; # Compton scattering cross-sections 
# const.Sigma_x = Sigma_T * 3/4 //1.+x//x^3 //2 x /1+x////1+2 x/-Log[1+2 x]/+1//2 x/ Log[1+2 x]-/1+3 x///1+2 x/^2/;
const.G = 6.6742e-8;#cgs
const.a = 7.5658e-15;#radiation constant cgs
const.c = 2.99792458e10;#cm/s
const.parsec = 3.0856776e18;#cm
const.mparsec = const.parsec*10**6;#cm
const.Alpha_fs= 1.0/137.036;# Fine structure constant 
const.kb=1.3806505e-16;#Boltzmann cgs
const.hb=1.05457168e-27;#Planck/2\[Pi] cgs
const.b = (2.404)/(math.pi**2)*((const.kb/(const.hb*const.c))**3); #radiation density constant cgs
# const.h=2 * \[Pi]*hb;
const.me= 9.1093826e-28;#mass of electron g
const.ev=1.60217646e-12;
const.GeV= const.ev * 10**9;
#Helium fraction Y=0.24
const.Nnu = 3.046; # Number of neutrinos 

#Cosmological  Parameters 

const.TCMB = 2.725; # CMB Temperature in K 
const.Omega_b = 4.8999e-2; # Baryon density 
const.Omega_cdm = 2.6709e-1; # Cold Dark matter density 
const.h0 = 0.6711; # Hubble parameter 
const.fmHe = 0.24; #Helium mass fraction 
const.zmin = 200.0; # Minimum redshift 

const.Omega_m = const.Omega_cdm + const.Omega_b; # Total Matter density 
const.H0 = const.h0*100*3.24077649*(10**(-20)) ;#s^-1
const.fnu = const.Nnu * (7.0/8) *  ((4/11)**(4/3));
const.zeq = (3 *((const.H0 * const.c)**2) * const.Omega_m)/(8 * math.pi * const.G * const.a * (const.TCMB**4) * (1+const.fnu))-1;
const.Omega_Lambda = 0.6825;
const.Omega_r = 1 - (const.Omega_Lambda + const.Omega_m);#Omega_m/(1+zeq);

const.Rho_cr = (3 * (const.H0**2)) / (8 * math.pi * const.G);


const.zmax = 5.0e6; # Maximum redshift 



# Dark matter annihilation example 

const.mdm = (10 * const.GeV)/(const.c**2); # Dark matter mass - !!!!! This will be eventually a function taking value from micOMEGAs 
const.f_Gamma = 1; # Fraction of energy that goes into particles with electromagnetic interaction and is deposited in the plasma 
const.sigma_v = (3e-27)/(const.Omega_cdm * (const.h0**2));

# Chemical potential calculation parameters from arXiv: 1203.2601v2 
const.C = 0.7768;
const.B = 1.803;
const.z_dC = 1.96e6;
const.z_br = 1.05e7;
const.z_eps = 3.67e5;
const.eps = 0.0151;
const.z_dC_prime = 7.11e6;
const.z_br_prime = 5.41e11;

class Wimp(object):
	def __init__(self, mass, definitions):
		self.mass = mass
		self.log_x = []
		self.dN_dx = {}
		self.definitions = definitions
		self.num_columns = len(definitions)
		self.int_x_dN_dx = {}
		for j in range(2, self.num_columns):
			self.dN_dx[definitions[j]] = []
			self.int_x_dN_dx[definitions[j]] = []
			#print definition, "added"
		#print "Wimp of mass", self.mass, "is added"
			
	def integration_of_columns(self):
		#print "Integration of columns"
		sum_integrations = 0.
		for j in range(2, self.num_columns):
			# Creating an array of the sample for integration \int E * dN/dx
			sampling_integrand = []
			for i in range(len(self.log_x)):
				sampling_integrand.append(self.mass*(10**(self.log_x[i]))*self.dN_dx[self.definitions[j]][i])
			# Integrating using scipy.integrate.simps
			integration_result = scipy.integrate.simps(sampling_integrand, self.log_x)
			#for char in self.definitions[j]:
				#if char in 'LRT':
					#sum_integrations = sum_integrations + integration_result
			#print "Integration results", integration_result
			# Storing the value of integration in the dictionary created int_x_dN_dx
			self.int_x_dN_dx[self.definitions[j]].append(integration_result)
			#print "The integration of column", self.definitions[j], "is", integration_result
		#difference = (2 * self.mass) - sum_integrations
		#print "For mDM ", self.mass, "the sum of the integration of columns is", sum_integrations
		#print "2* mDM - (sum of integrated columns) = ", difference
		
		
class CmbDistortion(object):
	""" (Wimp) -> CMB distortion type and magnitude
	This class implements all of the functions that are used to calculate the type and magnitude
	of distortion that a certain WIMP particle will produce at the frequency spectrum of the CMB
	in the early universe (redshift 10^6 > z > 10^4)
	"""
	def __init__(self, wimp):
		self.wimp = wimp
		# In order to decide whether the distortion is y-type or mu-type
		# The Compton y parameter should be calculated (eq. 2.5 in 1203.2601)
		self.calculate_compton_y_parameter()
		if self.compton_y_parameter > 1:
			self.calculate_mu_distortion()
		elif self.compton_y_parameter <= 1 and self.compton_y_parameter >= 0.01:
			self.calculate_y_distortion()
		else:
			pass
		
	def calculate_compton_y_parameter(self):
		self.compton_y_parameter = 0.5
		
	def calculate_mu_distortion(self):
		print "Calculating mu distortion.."
		
	def calculate_y_distortion(self):
		print "Calculating y distortion.."
			
					
"""
/*Cosmological parameters functions */
extern double ndm_z(double z)
{
	double ndm = (Rho_cr * Omega_cdm / mdm) * pow((1+z),3);
	return ndm;
}

extern double H(double z)
{
	double H_z = H0 * pow((pow((1+z),3)*Omega_m + pow((1+z),4)*Omega_r + (1-(Omega_m+Omega_r))),0.5);
	return H_z;
}

/*WIMP cross section*/ 
extern double Sigma_v(double sigma_v)
{
	//double sigma_v = (3 * pow(10,-27))/(Omega_cdm * pow(h0,2));
	return sigma_v;// /(Omega_cdm * pow(h0,2));
}

extern int photons_per_annihilation()
{
	return 2;
}
/*Energy injection calculation from eq. (5.5) in arXiv: 1203.2601v2 */
extern double dQdz(double z) /* Energy injection rate */
{
	return f_Gamma * (mdm * pow(c,2) * pow(ndm_z(z),2) * sigma_v)/(a * pow((TCMB * (1+z)),4));
}

extern double dNdz(double z)
{
	return photons_per_annihilation() * (pow(ndm_z(z),2) * sigma_v)/(b * pow((TCMB * (1+z)),3) );
}
extern double dummy_energy_injection(double z) /* Energy injection for a fixed WIMP candidate as in arXiv: 1203.2601v2 */
{
	double epsilon_dot = 1.4 * pow(10,-29) * pow((1+z),2) * f_Gamma *(10/*GeV*/ / 10 /* WIMP mass in GeV*/) * ((Omega_cdm * pow(h0,2))/0.105 );
	return epsilon_dot;
}


/*Blackbody optical depth calculation from eq. (4.5) in arXiv: 1203.2601v2 */
extern double tau(double z)
{
	double first_part_tau_1, second_part_tau_1, tau_1,tau_2,tau_precise;
	first_part_tau_1 = pow(((1 + z)/(1 + z_dC)),5);
	second_part_tau_1 = pow(((1 + z)/(1 + z_br)),(5/2));
	tau_1 = pow((pow(((1 + z)/(1 + z_dC)),5) + pow(((1 + z)/(1 + z_br)),(5/2))),0.5);
	tau_2 = eps * log((pow(((1 + z)/(1 + z_eps)),(5/4)))+(pow((1 + pow(((1 + z)/(1 + z_eps)),(5/2))),(1/2))));
	tau_precise = (pow(((1+z)/(1+z_dC_prime)),3) + pow(((1+z)/(1+z_br_prime)),0.5));
	return 1.007 * (tau_1+tau_2) + tau_precise;
}

/*Chemical potential calculation from eq. (3.6) in arXiv: 1203.2601v2 */
extern double mu(double z)
{
	double first_part = 7.43 * pow(10,-5) * ((1 + z)/(2 * pow(10,6)));
	double second_part = 1.07 * pow(10,-6) * pow(((1 + z)/(2 * pow(10,6))), (-3/2));
	double value = pow((first_part + second_part),(1/2));
	
	return value;
}

extern double mu_integrand(double z)
{
	return (1/((1+z)*H(z)))* (dQdz(z)-(4/3)*(dNdz(z))) * exp(-tau(z));
}

extern double mu_0(double z_i, double z_min, int subdivisions)
{
	double first_part = mu(z_i) * exp(-tau(z_i));
	double second_part = C*B*integrate_trapezoid((&mu_integrand),z_min, z_i, subdivisions);
	return first_part+second_part;
}

extern double mdm_calc(double mass_in_GeV)
{
	return ( mass_in_GeV * GeV)/pow(c,2);
}


"""		
		
