"""
Omer Tzuk, February 2014
[\mu] type distortion calculator Version 0.1
This script 

"""

import numpy as np
from astropy import units as u
from astropy import cosmology
from scipy import integrate
from scipy.interpolate import interp1d
import const
import math
import argparse
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

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

def wimp_mucalc(mDM, sigma_v, f_gamma):
	z_i = 2.0e6
	z_min = 5.0e4
	return mu_0(z_i, z_min,f_gamma, mDM, sigma_v)
	

H_z = lambda z : cosmology.H(z).to(1/u.s)

def ndm_z(z):
	"""
	"""
	return (const.Rho_cr * const.Omega_cdm / const.mdm) * ((1+z)**3)
	
def dQdz(f_gamma, mDM, sigma_v, z):
	"""
	Energy injection calculation from eq. (5.5) in arXiv: 1203.2601v2 
	"""
	mDM = (mDM * const.GeV)/(const.c**2)
	dQdz_numerator = (mDM * const.c**2 * (ndm_z(z)**2) * sigma_v)
	dQdz_denominator = (const.a * ((const.TCMB * (1+z))**4))
	return f_gamma * (dQdz_numerator/dQdz_denominator)
	
def dNdz(sigma_v,z):
	"""
	"""
	photons_per_annihilation = 2.
	return photons_per_annihilation * ((ndm_z(z)**2) * sigma_v)/(const.b * ((const.TCMB * (1+z))**3) )

def tau(z):
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

def mu(z): # I have problem in interpreting this equation from the paper, there must be mistake here..
	"""
	Chemical potential calculation from eq. (3.2) in arXiv: 1203.2601v2
	"""
	# What is mu_c and x_e ??
	mu_c = 0.
	x_e = 1.
	first_part = 7.43e-5 * ((1 + z)/(2e6))
	second_part = 1.07e-6 * (((1 + z)/(2e6))**(-3/2))
	x_c = (first_part + second_part)**(1/2)
	return mu_c * np.exp(-x_c/x_e)
	
def mu_0(z_i, z_min,f_gamma, mDM, sigma_v):
	"""
	Chemical potential calculation from eq. (3.6) in arXiv: 1203.2601v2
	I do not consider any initial mu, so the first part of the equation is not taken into account
	"""
	#first_part = mu(z_i) * np.math.exp(-tau(z_i))
	
	# With dN/dz
	#integrand = lambda z: ((1/((1+z)*H_z(z)))* (dQdz(z)-(4/3)*(dNdz(z))) * math.exp(-tau(z))) 
	# Without dN/dz
	integrand = lambda z: ((1/((1+z)*H_z(z)))* (dQdz(z)) * np.exp(-tau(z)))
	second_part = const.C*const.B*integrate.romberg(integrand,z_min,z_i)
	#return first_part+second_part
	return second_part
	

def main():
	sigma_v = 3e-27 / (const.Omega_cdm * (const.h0**2)) 
	dot_epsilon = lambda z: dQdz(1,10,sigma_v, z)
	y = lambda z: dot_epsilon(z) * (np.exp(-tau(z))/H_z(z).value)
	#print "sigma_v",sigma_v
	#print "dot_epsilon", dot_epsilon(10**5)
	#print "H", H_z(10**5)
	#print "tau",tau(10**5)
	#print "y axis of Figure 10 for z= 10^5 is :  ", y(10**5)
	#print "y axis of Figure 10 for z= 100 is :  ", y(10)
	#print "y axis of Figure 10 for z= 2.5*10^6 is :  ", y(2.5e6)
	plot(y,100, 2.5e6)

def plot(function, min_x, max_x):
	t = np.logspace(np.log10(min_x), np.log10(max_x),100000)
	s = function(t)
	plt.loglog(t, s, 'b-', lw=3)
	
	plt.xlabel(r'$z$')
	plt.ylabel(r'$(1+z)G d \varepsilon / dz$')
	plt.title(r'Energy injection from dark matter annihilation for 10 GeV WIMP with $f_{\gamma} = 1$')
	plt.grid(True, which="both")
	plt.show()

	



# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')

args = parser.parse_args()	
if __name__ == "__main__":
	main()
