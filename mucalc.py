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
import math
import argparse
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

cosmology.core.set_current(cosmology.Planck13)

from mucalc_constants import *
import ucmh

def wimp_mucalc(mDM, sigma_v, f_gamma):
	z_i = 2.0e6
	z_min = 5.0e4
	mDM = (mDM * const.GeV)/(const.c**2)
	return mu_0(z_i, z_min,f_gamma, mDM, sigma_v)
	

H_z = lambda z : cosmology.H(z).to(1/u.s).value

def ndm(mDM, z):
	"""
	"""
	return (const.Rho_cr * const.Omega_cdm / mDM) * ((1+z)**3)
	
def dQdz(f_gamma, mDM, sigma_v, z):
	"""
	Energy injection calculation from eq. (5.5) in arXiv: 1203.2601v2 
	"""
	
	dQdz_numerator = (mDM * const.c**2 * (ndm(mDM,z)**2) * sigma_v)
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
	dQdz_z = lambda z : dQdz(f_gamma, mDM, sigma_v, z)
	integrand = lambda z: ((1/((1+z)*H_z(z)))* (dQdz_z(z)) * np.exp(-tau(z)))
	second_part = const.C*const.B*integrate.romberg(integrand,z_min,z_i)
	#return first_part+second_part
	return second_part
	

def main():
	sigma_v = 3.0e-27 / (const.Omega_cdm * (const.h0**2)) 
	f_gamma = 1.
	#mDM = (10 * const.GeV)/(const.c**2)
	#dot_epsilon = lambda z: dQdz(f_gamma,mDM,sigma_v, z)
	#y = lambda z: dot_epsilon(z) * (np.exp(-tau(z))/H_z(z))
	#print "sigma_v",sigma_v
	#print "dot_epsilon", dot_epsilon(10**5)
	#print "H", H_z(10**5)
	#print "tau",tau(10**5)
	#print "y axis of Figure 10 for z= 10^5 is :  ", y(10**5)
	#print "y axis of Figure 10 for z= 100 is :  ", y(10)
	#print "y axis of Figure 10 for z= 2.5*10^6 is :  ", y(2.5e6)
	#dQdz_z = lambda z : dQdz(f_gamma, mDM, sigma_v, z)
	#integrand = lambda z: ((1/((1+z)*H_z(z)))* (dQdz_z(z)) * np.exp(-tau(z)))
	#print "integrand",integrand(10**5)
	#z_i = 2.0e6
	#z_min = 5.0e4
	#print "mu", integrate.romberg(integrand,z_min,z_i)
	#print "mu_0", mu_0(z_i, z_min, f_gamma, mDM, sigma_v)
	print "mu", wimp_mucalc(10, sigma_v, f_gamma)
	#print "mu distortion magnitude", wimp_mucalc(10,sigma_v, 1) 
	#plot(y,100, 2.5e6)

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
