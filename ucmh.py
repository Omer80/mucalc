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
from scipy.stats import norm
from mucalc_constants import *
import decimal
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
t_i = cosmology.Planck13.age(T_phase_transition.value/(cosmology.Planck13.Tcmb0.value * (8.617e-5))-1.).to(u.s).value

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
	return M_h.to(u.g).value 

solarmass = 1. * u.solMass
	
def R_h(z):
	'''(float) -> float (in cm)
	R_h is the maximum extent of the UCMH at redshift z
	as in eq. (6) at arXiv:09084.082v5
	'''
	R_h = 0.019 * (1000./(z + 1.)) * ((M_h(z)/(solarmass.to(u.g).value))**(1./3)) * u.pc # gives R_h in parsec
	return R_h.to(u.cm).value # const.parsec

t = lambda z: cosmology.Planck13.age(z).to(u.s).value


rho_max = lambda mDM,sigma_v, z: (mDM / (sigma_v *(t(z) - t_i)))


rho_chi = lambda r,z:(3. * f_chi * M_h(z)) / (16. * np.pi * (R_h(z)**(3./4)) * (r**(9./4))) # eq. (5) at arXiv:0908.4082v5
r_rho = lambda rho_chi,z:((3. * f_chi * M_h(z)) / (16. * np.pi * (R_h(z)**(3./4)) * (rho_chi)))**(4./9) # inverse eq. (5) at arXiv:0908.4082v5
sigma_ucmh_pert = lambda n: (9.5e-5 * (M_H_z_X.to(u.g).value/10e56)**((1. - n)/4.)) # eq. (9) at arXiv:0908.4082v5

def Sigma_ucmh(n, z):
	'''(float) -> (float)
	eq. (8) from arXiv:0908.4082v5
	'''
	#print "n", n
	#x = lambda delta: delta/sigma_ucmh_pert(n)
	#print "sigma_ucmh_pert(n)",sigma_ucmh_pert(n)
	Omega_CMD_z = f_chi * cosmology.Planck13.Om(z)
	#print "Sigma_ucmh first part",(Omega_CMD_z / np.sqrt(2.* np.pi))
	#int_f = lambda delta: (1./2) * special.erf(x(delta)/np.sqrt(2))
	#print "int_f(0.3)",int_f(0.3),"int_f(10e-3)",int_f(1e-3)
	#Sigma_ucmh_z =  Omega_CMD_z * float(int_f(0.3) - int_f(1.e-3))
	
	# Analytical integration using scipy.stats
	x0 = 1.e-3
	x1 = 0.3
	s = sigma_ucmh_pert(n)
	
	Sigma_ucmh_z = norm.sf(x0, scale = s) - norm.sf(x1, scale = s)
	
	# Numerical integration
	#integrand = lambda delta: (1./(np.sqrt(2 * np.pi)*sigma_ucmh_pert(n))) * np.exp(-(delta**2)/(2. * sigma_ucmh_pert(n)**2))
	#Sigma_ucmh_z = Omega_CMD_z * integrate.romberg(integrand,(10.**(-3)),0.3,tol=1.e-41, divmax=120)
	return	Sigma_ucmh_z

ndm_0 = lambda mDM, z: (const.Rho_cr * const.Omega_cdm / mDM) * ((1+z)**3)

n_ucmh = lambda n, z: (cosmology.critical_density(z).value * Sigma_ucmh(n, z)) / M_h(z)

def avg_n_ucmh_squared(mDM, sigma_v, z, n):
	'''
	'''
	#print "input", mDM, sigma_v, z, n
	r_cut = r_rho(rho_max(mDM, sigma_v,z),z)
	#print "r_cut", r_cut
	n_max_z = rho_max(mDM,sigma_v, z)/mDM
	#print "n_max_z",n_max_z
	ndm_0_z = ndm_0(mDM, z)
	R_h_z = R_h(z)
	#print "R_h_z",R_h_z
	first_part = ((r_cut**3.)/3.)* (n_max_z**2.)
	#print "n_ucmh(n,z)",n_ucmh(n,z)
	#print first_part
	f = (3. * f_chi * M_h(z)) / ( 16. * np.pi * (R_h_z**(3./4)))
	second_part = (2./3) * ((f/mDM)**2) * (1./(r_cut**(3./2)))
	#print second_part
	avg_n_ucmh_squared_z = 4. *np.pi * n_ucmh(n,z) *(first_part +  second_part)
	
	# Adding cross term
	first_part = ((r_cut**3.)/3.)* (n_max_z* ndm_0_z)
	second_part = (4./3) * (f/mDM) * ndm_0_z * (R_h_z**(3./4) - r_cut**(3./4))
	avg_n_ucmh_squared_z = avg_n_ucmh_squared_z + 8. * n_ucmh(n,z) *np.pi*(first_part +  second_part)
	return avg_n_ucmh_squared_z


def main():
	sigma_v = 3.0e-27 / (const.Omega_cdm * (const.h0**2))
	mDM = (10 * const.GeV)/(const.c**2)
	z = 2.e6
	n = 1.3
	#r_rho = lambda rho_chi,z:((3. * f_chi * M_h(z)) / (16. * np.pi * (R_h(z)**(3./4)) * (rho_chi)))**(4./9)
	print "rho_max", rho_max(mDM, sigma_v, z)
	print "f_chi", f_chi
	print "M_h", M_h(z)
	print "R_h", R_h(z)
	print "r_cut", r_rho(rho_max(mDM, sigma_v, z), z)
	print "parameters", sigma_v, mDM, z, n
	ndm_0_sq_n = lambda n: (ndm_0(mDM, z)**2)
	ndm_ucmh_sq_n = lambda n: avg_n_ucmh_squared(mDM, sigma_v, z, n)
	print "ndm_0_sq_n(1.3)",ndm_0_sq_n(1.3)
	print "ndm_ucmh_sq_n(1.3)",ndm_ucmh_sq_n(1.3)
	#print "Sigma_ucmh(z)", Sigma_ucmh(n, z)
	#print "avg_n_ucmh_squared(mDM, sigma_v, z, n)", avg_n_ucmh_squared(mDM, sigma_v, z, n)
	#print "Halo size", R_h(z)
	#print "n_ucmh",n_ucmh(z)
	#print avg_n_ucmh_squared(mDM, sigma_v, z)
	#n_z = lambda r: rho_chi(mDM, sigma_v,r,z)/mDM
	#print avg_n_ucmh_squared(mDM, sigma_v, 10**5)


if __name__ == "__main__":
	main()
