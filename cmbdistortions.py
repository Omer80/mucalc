"""
Omer Tzuk, February 2014
Version 0.1
This script 

"""

import numpy as np
from astropy import units as u
from astropy.units import imperial
from astropy import cosmology
from scipy import integrate
from scipy.interpolate import interp1d
import const
import math

cosmology.core.set_current(cosmology.Planck13)


#Physical constants 
const.Sigma_T = 6.65246e-25 # Compton scattering cross-sections 
# const.Sigma_x = Sigma_T * 3/4 //1.+x//x^3 //2 x /1+x////1+2 x/-Log[1+2 x]/+1//2 x/ Log[1+2 x]-/1+3 x///1+2 x/^2/
const.G = 6.6742e-8#cgs
const.a = 7.5658e-15#radiation constant cgs
const.c = 2.99792458e10#cm/s
const.parsec = 3.0856776e18#cm
const.mparsec = const.parsec*10**6#cm
const.Alpha_fs= 1.0/137.036# Fine structure constant 
const.kb=1.3806505e-16#Boltzmann cgs
const.hb=1.05457168e-27#Planck/2\[Pi] cgs
const.b = (2.404)/(math.pi**2)*((const.kb/(const.hb*const.c))**3) #radiation density constant cgs
# const.h=2 * \[Pi]*hb
const.me= 9.1093826e-28#mass of electron g
const.ev=1.60217646e-12
const.GeV= const.ev * 10**9
#Helium fraction Y=0.24
const.Nnu = 3.046 # Number of neutrinos 

#Cosmological  Parameters 

const.TCMB = 2.725 # CMB Temperature in K 
const.Omega_b = 4.8999e-2 # Baryon density 
const.Omega_cdm = 2.6709e-1 # Cold Dark matter density 
const.h0 = 0.6711 # Hubble parameter 
const.fmHe = 0.24 #Helium mass fraction 
const.zmin = 200.0 # Minimum redshift 

const.Omega_m = const.Omega_cdm + const.Omega_b # Total Matter density 
const.H0 = const.h0*100*3.24077649*(10**(-20)) #s^-1
const.fnu = const.Nnu * (7.0/8) *  ((4/11)**(4/3))
const.zeq = (3 *((const.H0 * const.c)**2) * const.Omega_m)/(8 * math.pi * const.G * const.a * (const.TCMB**4) * (1+const.fnu))-1
const.Omega_Lambda = 0.6825
const.Omega_r = 1 - (const.Omega_Lambda + const.Omega_m)#Omega_m/(1+zeq)

const.Rho_cr = (3 * (const.H0**2)) / (8 * math.pi * const.G)


const.zmax = 5.0e6 # Maximum redshift 



# Dark matter annihilation example 

const.mdm = (10 * const.GeV)/(const.c**2) # Dark matter mass - !!!!! This will be eventually a function taking value from micOMEGAs 
const.f_Gamma = 1 # Fraction of energy that goes into particles with electromagnetic interaction and is deposited in the plasma 
const.sigma_v = (3e-27)/(const.Omega_cdm * (const.h0**2))

# Chemical potential calculation parameters from arXiv: 1203.2601v2 
const.C = 0.7768
const.B = 1.803
const.z_dC = 1.96e6
const.z_br = 1.05e7
const.z_eps = 3.67e5
const.eps = 0.0151
const.z_dC_prime = 7.11e6
const.z_br_prime = 5.41e11

class Wimp(object):
	def __init__(self, mass):
		self.mass = mass
		self.f_gamma = 1.
		self.mu_distortion = 0
		self.y_distortion = 0
		self.compton_y_parameter = 1
		self.z_annihilation = 10e5
		self.sigma_v = const.sigma_v
		self.CMBdistortion = CmbDistortion() # Need to think what I should pass  to calculation
		
	def ndm_z(self, z):
		"""
		"""
		return (const.Rho_cr * const.Omega_cdm / const.mdm) * ((1+z)**3)
		
	def dQdz(self, z):
		"""
		Energy injection calculation from eq. (5.5) in arXiv: 1203.2601v2 
		"""
		return self.f_gamma * (self.mass * const.c**2 * (self.ndm_z(z)**2) * self.sigma_v)/(const.a * ((const.TCMB * (1+z))**4))
		
	def dNdz(self, z):
		"""
		"""
		photons_per_annihilation = 2
		return photons_per_annihilation * ((self.ndm_z(z)**2) * self.sigma_v)/(const.b * ((const.TCMB * (1+z))**3) )
		
	def dummy_energy_injection(self, z):
		"""
		Energy injection for a fixed WIMP candidate as in arXiv: 1203.2601v2
		"""
		epsilon_dot = (1.4e-29)*((1+z)**2)*self.f_gamma*(10 / self.mass ) * ((const.Omega_cdm * (const.h0**2))/0.105 )
		return epsilon_dot
			
			
		
class CmbDistortion(object):
	""" (Wimp) -> CMB distortion type and magnitude
	This class implements all of the functions that are used to calculate the type and magnitude
	of distortion that a certain WIMP particle will produce at the frequency spectrum of the CMB
	in the early universe (redshift 10^6 > z > 10^4)
	"""
	def __init__(self):
		# In order to decide whether the distortion is y-type or mu-type
		# The Compton y parameter should be calculated (eq. 2.5 in 1203.2601)
		self.calculate_compton_y_parameter(10**5)
		if self.compton_y_parameter > 1:
			self.calculate_mu_distortion()
		elif self.compton_y_parameter <= 1 and self.compton_y_parameter >= 0.01:
			self.calculate_y_distortion()
		else:
			pass
		
	def calculate_compton_y_parameter(self, z):
		"""
		Calculation of Compton y parameter based on Eq. 2.5 in 1203.2601
		"""
		H_z = cosmology.H(z).to(1/u.s)
		integrand = lambda z : (const.kb*const.Sigma_T)/(const.me*const.c) *((n_e * TCMB_z)/(H_z * (1+z)))
		self.compton_y_parameter = 0.5
		
	def calculate_mu_distortion(self):
		print "Calculating mu distortion.."
		
	def calculate_y_distortion(self):
		print "Calculating y distortion.."
		
	def tau(self, z):
		"""
		Blackbody optical depth calculation from eq. (4.5) in arXiv: 1203.2601v2
	    """
		z_dC = const.z_dC
		z_br = const.z_br
		z_eps = const.z_eps
		eps = const.eps
		z_dC_prime = const.z_dC_prime
		z_br_prime = const.z_br_prime
		# Calculation of tau
		first_part_tau_1 = (((1 + z)/(1 + z_dC))**5)
		second_part_tau_1 = (((1 + z)/(1 + z_br))**(5/2))
		tau_1 = ((((1 + z)/(1 + z_dC))**5) + (((1 + z)/(1 + z_br))**(5/2)))**0.5
		tau_2 = eps * np.log((((1 + z)/(1 + z_eps))**(5/4))+((1 + (((1 + z)/(1 + z_eps))**(5/2)))**(1/2)))
		tau_precise = ((((1+z)/(1+z_dC_prime))**3) + (((1+z)/(1+z_br_prime))**0.5))
		return 1.007 * (tau_1+tau_2) + tau_precise
	
	def mu(self, z):
		"""
		Chemical potential calculation from eq. (3.2) in arXiv: 1203.2601v2
		"""
		first_part = 7.43e-5 * ((1 + z)/(2e6))
		second_part = 1.07e-6 * (((1 + z)/(2e6))**(-3/2))
		return (first_part + second_part)**(1/2)
		
	def mu_0(self, z_i, z_min):
		"""
		Chemical potential calculation from eq. (3.6) in arXiv: 1203.2601v2
		"""
		first_part = self.mu(z_i) * np.exp(-self.tau(z_i))
		H_z = cosmology.H(z).to(1/u.s)
		integrand = lambda z: (1/((1+z)*H_z))* (self.wimp.dQdz(z)-(4/3)*(self.wimp.dNdz(z))) * exp(-self.tau(z))
		second_part = C*B*integrate.romberg(integrand,z_min,z_i)
		return first_part+second_part

			
if __name__ == "__main__":
	print cosmology.get_current()
	test_wimp = Wimp(10)
	print test_wimp.dummy_energy_injection(10**4)
					
"""


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


extern double mdm_calc(double mass_in_GeV)
{
	return ( mass_in_GeV * GeV)/pow(c,2);
}


"""		
		
