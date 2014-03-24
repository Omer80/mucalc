"""
Omer Tzuk, March 2014
Version 0.1
Based on arXiv:0908.4082v5 this code is building the function related to Ultracompact Dark Matter Minihalos
Those equations will be used to estimate the amplification in WIMP annihilation at the early Universe
"""
import numpy as np
from astropy import units as u
from astropy import constants as c
from astropy import cosmology
from scipy import integrate
from scipy import special
from mucalc_constants import *
import decimal
import argparse
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

cosmology.core.set_current(cosmology.Planck13)

# delta_m as calculated from eq. (2) at arXiv:0908.4082v5
delta_m_EW  = 5.6e-19*u.solMass
delta_m_QCD = 1.1e-9*u.solMass
delta_m_ee  = 0.33*u.solMass

delta_m = delta_m_ee

z_eq = (2.32e4 *  cosmology.Planck13.Om(0) *((cosmology.H(0).value)/100.)**2.) -1.
f_chi = 0.834 # as in arXiv:0908.4082v5

T_EW = 200e9  * u.eV
T_QCD = 200e6 * u.eV
T_ee = 0.51e3   * u.eV

# Setting up the phase transition time to be at the positron-electron annihilation phase
T_phase_transition = T_ee
# Calculating time of phase transition
t_i = cosmology.Planck13.age(T_phase_transition.value/(cosmology.Planck13.Tcmb0.value * (8.617e-5))-1).to(u.s).value

# Perturbation spectrum index
n=1.

# Horizon mass parameters
M_H_Teq = 3.5e17 * u.solMass
T_eq = 0.75 * u.eV

g_s_eq = 3.91

def g_s(T):
	if T == T_EW.value:
		_g_s = 107.
	elif T == T_QCD.value:
		_g_s = 55.
	elif T == T_ee.value:
		_g_s = 10.8
	return _g_s

M_H = lambda T: M_H_Teq * ((((g_s_eq**(1./3.))*T_eq.value)/((g_s(T)**(1./3.))*T))**2.)

M_H_z_X = M_H(T_phase_transition.value)

def M_h(z):
	'''(float) -> float in units.solMass (from astropy)
	M_h is the total mass of the UCMH at redshift z
	equation (1) at arXiv:0908.4082v5
	'''
	M_h = delta_m * ((1. + z_eq)/(1. + z)) # delta_m is in solMass
	return M_h.to(u.g).value # M_h is in solMass !

solarmass = 1. * u.solMass
	
def R_h(z):
	'''(float) -> float (in cm)
	R_h is the maximum extent of the UCMH at redshift z
	as in eq. (6) at arXiv:09084.082v5
	'''
	R_h = 0.019 * (1000./(z + 1.)) * ((M_h(z)/(solarmass.to(u.g).value))**(1./3)) * u.pc # gives R_h in parsec
	return R_h.to(u.cm).value # const.parsec

t = lambda z: cosmology.Planck13.age(z).to(u.s).value


rho_max = lambda mDM,sigma_v, z: (1. / (sigma_v *(t(z) - t_i)))


rho_chi = lambda r,z:(3. * f_chi * M_h(z)) / (16. * np.pi * (R_h(z)**(3./4)) * (r**(9./4))) # eq. (5) at arXiv:0908.4082v5
sigma_ucmh_pert = lambda n: (9.5e-5 * (M_H_z_X/10e56)**((1. - n)/4.)).value # eq. (9) at arXiv:0908.4082v5

def Sigma_ucmh(z):
	'''(float) -> (float)
	eq. (8) from arXiv:0908.4082v5
	'''
	x = lambda delta: delta/sigma_ucmh_pert(n)
	print "sigma_ucmh_pert(n)",sigma_ucmh_pert(n)
	Omega_CMD_z = f_chi * cosmology.Planck13.Om(z)
	print "Sigma_ucmh first part",(Omega_CMD_z / np.sqrt(2.* np.pi))
	int_f = lambda delta: (1./2) * special.erf(x(delta)/np.sqrt(2))
	print "int_f(0.3)",int_f(0.3),"int_f(10e-3)",int_f(10e-3)
	Sigma_ucmh_z =  Omega_CMD_z * float(int_f(0.3) - int_f(10e-3))
	print "Sigma_ucmh_z",Sigma_ucmh_z
	return	Sigma_ucmh_z

n_ucmh = lambda z: (const.Rho_cr * Sigma_ucmh(z)) / M_h(z)

def avg_n_ucmh_squared(mDM, sigma_v, z):
	'''
	'''
	n_maz_z = rho_max(mDM,sigma_v, z)/mDM
	first_part = ((R_h(z)**3.)/3.)* (n_maz_z**2.)
	f = (3. * f_chi * M_h(z)) / (mDM * 16. * np.pi * (R_h(z)**(3./4)))
	second_part = (2./3) * (f**2) * (1./(R_h(z)**(3./2)))
	avg_n_ucmh_squared_z = 4 * n_ucmh(z) *np.pi*(first_part +  second_part)
	return avg_n_ucmh_squared_z


def main():
	sigma_v = 3.0e-27 / (const.Omega_cdm * (const.h0**2))
	mDM = (10 * const.GeV)/(const.c**2)
	z = 10e5
	print "Halo size", R_h(z)
	#print "n_ucmh",n_ucmh(z)
	print avg_n_ucmh_squared(mDM, sigma_v, z)
	#n_z = lambda r: rho_chi(mDM, sigma_v,r,z)/mDM
	#print avg_n_ucmh_squared(mDM, sigma_v, 10**5)

# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')

args = parser.parse_args()	
if __name__ == "__main__":
	main()
