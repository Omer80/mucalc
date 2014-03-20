"""
Omer Tzuk, March 2014
Version 0.1
Based on arXiv:0908.4082v5 this code is building the function related to Ultracompact Dark Matter Minihalos
Those equations will be used to estimate the amplification in WIMP annihilation at the early Universe
"""
import numpy as np
from astropy import units as u
from astropy import constants as const
from astropy import cosmology
from scipy import integrate
from mucalc_constants import *
import argparse


cosmology.core.set_current(cosmology.Planck13)

# delta_m as calculated from eq. (2) at arXiv:0908.4082v5
delta_m_EW  = 5.6e-19*u.solMass
delta_m_QCD = 1.1e-9*u.solMass
delta_m_ee  = 0.33*u.solMass

delta_m = delta_m_ee

z_eq = (2.32e4 *  cosmology.Planck13.Om(0) *((cosmology.H(0).value)/100.)**2.) -1.
f_chi = 0.834 # as in arXiv:0908.4082v5

def M_h(z):
	'''(float) -> float in units.solMass (from astropy)
	M_h is the total mass of the UCMH at redshift z
	equation (1) at arXiv:0908.4082v5
	'''
	M_h = delta_m * ((1. + z_eq)/(1. + z)) # delta_m is in solMass
	return M_h # M_h is in solMass !
	
def R_h(z):
	'''(float) -> float (in cm)
	R_h is the maximum extent of the UCMH at redshift z
	as in eq. (6) at arXiv:09084.082v5
	'''
	R_h = 0.019 * (1000./(z + 1.)) * ((M_h(z).value)**(1./3)) # gives R_h in parsec
	return R_h * u.pc # const.parsec

t = lambda z: cosmology.Planck13.age(z).to(u.s).value

rho_max = lambda mDM, sigma_v, t_i, z : mDM / (sigma_v *(t(z) - t_i))

def r_cut(mDM, sigma_v, t_i,z):
	''' 
	'''
	Rh_z = R_h(z).to(u.cm).value
	r_cut = ((3.*f_chi*(M_h(z).to(u.g).value)) / (16.*np.pi*(Rh_z**(3./4))*(rho_max(z))))**(4./9)
	return r_cut
	
		
def rho_chi(mDM, sigma_v, t_i,r,z):
	'''(float,float) -> float
	rho_chi as in eq. (5) at arXiv:0908.4082v5
	'''
	if r < r_cut(z):
		rho_chi = rho_max(z)
	elif r != 0.:
		Rh_z = R_h(z).to(u.cm).value
		rho_chi = (3.*f_chi*M_h(z).to(u.g)) / (16.*np.pi*(Rh_z**(3./4))*(r**(9./4)))
	else:
		rho_chi = rho_max(z)
	return rho_chi

sigma_ucmh_pert = lambda n: 9.5e-5 * (M_H_z_X/10e56)**((1. - n)/4.) # eq. (9) at arXiv:0908.4082v5

def Sigma_ucmh(z):
	'''(float) -> (float)
	eq. (8) from arXiv:0908.4082v5
	'''
	integrand = lambda delta : (1/(np.sqrt(2*np.pi)*sigma_ucmh_pert)) * np.exp(-(delta**2)/(2*(sigma_ucmh_pert**2)))
	Omega_CMD_z = f_chi * cosmology.Planck13.Om(z)
	Sigma_ucmh_z =  Omega_CMD_z * integrate.romberg(integrand,10-3,0.3)
	return	Sigma_ucmh_z

n_ucmh = lambda z: (const.Rho_cr * Sigma_ucmh(z)) / M_h(z)

def avg_rho_ucmh_squared(z):
	'''
	'''
	rho_z = lambda r: ''' '''
	integrand = lambda r: 4*np.pi*(r**2)*((rho_z(r))**2)
	return n_ucmh(z) *  integrate.romberg(integrand,0,R_h(z))


def main():
	a = 1. * u.solMass
	print a.to(u.g)
	print M_h(10**5)
	print R_h(10**5)
	
	
# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')

args = parser.parse_args()	
if __name__ == "__main__":
	main()
