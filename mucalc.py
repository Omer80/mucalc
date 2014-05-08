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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib import cm
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

cosmology.core.set_current(cosmology.Planck13)

from mucalc_constants import *
import ucmh
import argparse
from model import Wimp_model

def wimp_mucalc(WIMP_model,spectral_index = 1., with_ucmh = False):
	z_i = 2.0e6
	z_min = 5.0e4
	return mu_0(z_i, z_min,WIMP_model.f_gamma, WIMP_model.mDM, 
	            WIMP_model.sigma_v,spectral_index, with_ucmh)
	
H_z = lambda z : cosmology.H(z).to(1/u.s).value

def ndm_0(mDM, z):
	"""(float, float) -> float
	Giving the average number density of Darm Matter, without considering UCMH
	"""
	return (const.Rho_cr * const.Omega_cdm / mDM) * ((1+z)**3)
	
def ndm_squared(mDM, sigma_v, z,n, with_ucmh):
	'''(float, float, float) -> float
	Choosing between Darm matter number density with and without UCMH
	'''
	ndm_squared_z = (ndm_0(mDM, z)**2)
	if with_ucmh == True:
		ndm_squared_z = ndm_squared_z + ucmh.avg_n_ucmh_squared(mDM, sigma_v, z, n)
	return ndm_squared_z
	
def dQdz(f_gamma, mDM, sigma_v, z,n, with_ucmh):
	"""
	Energy injection calculation from eq. (5.5) in arXiv: 1203.2601v2 
	"""
	
	dQdz_numerator = (mDM * const.c**2 * (ndm_squared(mDM, sigma_v,z,n, with_ucmh)) * sigma_v)
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

def mu(z):
	"""
	Chemical potential calculation from eq. (3.2) in arXiv: 1203.2601v2
	"""
	mu_c = 0.
	x_e = 1.
	first_part = 7.43e-5 * ((1 + z)/(2e6))
	second_part = 1.07e-6 * (((1 + z)/(2e6))**(-3/2))
	x_c = (first_part + second_part)**(1/2)
	return mu_c * np.exp(-x_c/x_e)
	
def mu_0(z_i, z_min,f_gamma, mDM, sigma_v,n, with_ucmh):
	"""
	Chemical potential calculation from eq. (3.6) in arXiv: 1203.2601v2
	I do not consider any initial mu, so the first part of the equation is not taken into account
	"""
	#first_part = mu(z_i) * np.math.exp(-tau(z_i))
	# With dN/dz
	#integrand = lambda z: ((1/((1+z)*H_z(z)))* (dQdz(z)-(4/3)*(dNdz(z))) * math.exp(-tau(z)))
	
	# Without dN/dz
	dQdz_z = lambda z : dQdz(f_gamma, mDM, sigma_v, z,n, with_ucmh)
	integrand = lambda z: ((1/((1+z)*H_z(z)))* (dQdz_z(z)) * np.exp(-tau(z)))
	second_part = const.C*const.B*integrate.romberg(integrand,z_min,z_i)
	#return first_part+second_part
	return second_part
	
def mu0_n_mDM(mDM, spectral_index):
	WIMP = Wimp_model(mDM,1.,3.0e-27)
	return wimp_mucalc(WIMP,spectral_index, with_ucmh = True)

def main(args):
	WIMP = Wimp_model(10.,1.,3.0e-27)
	z = 2.e6
	if args.print_density_squared_plot:
		plot_density_squared(WIMP)
	if args.print_energy_injection:
		WIMP.mDM = (float(args.print_energy_injection[0]) * const.GeV)/(const.c**2)
		if args.spectral_index:
			spectral_index = float(args.spectral_index[0])
			plot_energy_injection(WIMP, spectral_index, True)
		else:
			plot_energy_injection(WIMP)
	if args.print_mu_to_n:
		plot_mu_to_spectral_index(WIMP)
	if args.print_contours_mu_to_n_mDM:
		plot_mu_contours()
	else:
		print "T_phase_transition", ucmh.T_phase_transition
		print "mu0 with n=1",wimp_mucalc(WIMP)
		print "mu0 with n=1.3",wimp_mucalc(WIMP,1.3,True)
		ucmh.T_phase_transition = ucmh.T_QCD
		ucmh.M_H_z_X = ucmh.M_H(ucmh.T_phase_transition.value)
		print "T_phase_transition", ucmh.T_phase_transition
		print "mu0 with n=1.3",wimp_mucalc(WIMP,1.3,True)
		ucmh.T_phase_transition = ucmh.T_EW
		ucmh.M_H_z_X = ucmh.M_H(ucmh.T_phase_transition.value)
		print "T_phase_transition", ucmh.T_phase_transition	
		print "mu0 with n=1.3",wimp_mucalc(WIMP,1.3,True)	
		
	
# Plotting functions definitions
def plot_mu_to_spectral_index(model):
	with_ucmh = True
	n_range = np.linspace(1.,1.3,100)
	mu_n = lambda spectral_index: wimp_mucalc(model,spectral_index, with_ucmh)
	# Calculation for electron positron annihilation
	#ucmh.T_phase_transition = ucmh.T_ee
	mu_ee = []
	for n in n_range:
		mu_ee.append(mu_n(n))
	mu_ee = np.asarray(mu_ee)
	# Calculation for QCD phase transition
	ucmh.T_phase_transition = ucmh.T_QCD
	ucmh.M_H_z_X = ucmh.M_H(ucmh.T_phase_transition.value)
	mu_QCD = []
	for n in n_range:
		mu_QCD.append(mu_n(n))
	mu_QCD = np.asarray(mu_QCD)
	# Calculation for electroweak phase transition
	ucmh.T_phase_transition = ucmh.T_EW
	ucmh.M_H_z_X = ucmh.M_H(ucmh.T_phase_transition.value)
	mu_EW = []
	for n in n_range:
		mu_EW.append(mu_n(n))
	mu_EW = np.asarray(mu_EW)
	
	plt.yscale('log')
	plt.plot(n_range, mu_ee, 'b', n_range, mu_QCD, 'g',n_range, mu_EW, 'm')
	plt.xlabel(r'$n$')
	plt.ylabel(r'$ | \mu |$')
	plt.title(r'Change of $\mu$-type distortion as function of spectral index $n$')
	plt.legend({r"$e^+ e^-$ annihilation epoch ",r"QCD phase transition",r"Electroweak phase transition"},loc = 4)
	plt.axhline(y=9e-5, linewidth=2, ls='--',color='r',label="Limit from COBE observations")
	plt.axvline(x=1.25, linewidth=2, ls='--',color='r',label="Constraints from PBH")
	plt.grid(True, which="both")
	#plt.figure(1)
	#plt.yscale('log')
	#plt.plot(n_range, mu_ee, 'b', lw=1, label = r"$e^+ e^-$ annihilation epoch ")
	#plt.xlabel(r'$n$')
	#plt.ylabel(r'$ | \mu |$')
	#plt.title(r'Change of $\mu$-type distortion as function of spectral index $n$')
	#plt.legend(loc = 4)
	#plt.axhline(y=9e-5, linewidth=2, ls='--',color='r',label="Limit from COBE observations")
	#plt.axvline(x=1.25, linewidth=2, ls='--',color='r',label="Constraints from PBH")
	##plt.grid(True, which="both")
	#plt.figure(2)
	#plt.yscale('log')
	#plt.plot(n_range, mu_QCD, 'g', lw=1, label = r"QCD phase transition")
	#plt.xlabel(r'$n$')
	#plt.ylabel(r'$ | \mu |$')
	#plt.title(r'Change of $\mu$-type distortion as function of spectral index $n$')
	#plt.legend(loc = 4)
	#plt.axhline(y=9e-5, linewidth=2, ls='--',color='r',label="Limit from COBE observations")
	#plt.axvline(x=1.25, linewidth=2, ls='--',color='r',label="Constraints from PBH")
	##plt.grid(True, which="both")
	#plt.figure(3)
	#plt.yscale('log')
	#plt.plot(n_range, mu_EW, 'm', lw=1, label = r"Electroweak phase transition")
	#plt.xlabel(r'$n$')
	#plt.ylabel(r'$ | \mu |$')
	#plt.title(r'Change of $\mu$-type distortion as function of spectral index $n$')
	#plt.legend(loc = 4)
	#plt.axhline(y=9e-5, linewidth=2, ls='--',color='r',label="Limit from COBE observations")
	#plt.axvline(x=1.25, linewidth=2, ls='--',color='r',label="Constraints from PBH")
	##plt.grid(True, which="both")
	plt.show()
		
def plot_density_squared(model):
	z = 2.e6
	ndm_0_n = lambda n: (ndm_0(model.mDM, z)**2) * (n**0)
	ndm_ucmh_n = lambda n: ucmh.avg_n_ucmh_squared(model.mDM, model.sigma_v, z, n)
	n_min = 1.
	n_max = 1.30

	t = np.linspace(n_min, n_max,100)
	s1 = ndm_0_n(t)
	s2 = ndm_ucmh_n(t)
	s3 = ndm_0_n(t) + ndm_ucmh_n(t)
	plt.yscale('log')
	plt.plot(t, s1, 'b--', lw=1, label = r"Only normal DM $density^2$")
	plt.plot(t, s2, 'r--', lw=1, label = r"$Density^2$ only UCMH")
	plt.plot(t, s3, 'g-', lw=2, label = r"$Density^2$ normal plus UCMH")
	
	plt.xlabel(r'spectral index $n$')
	plt.ylabel(r'$< {density^{2}} > $')
	plt.title(r'Change of average of $density^2$ as function of spectral index')
	plt.legend(loc = 4)
	#plt.grid(True, which="both")
	plt.show()	

def plot_energy_injection(model, spectral_index = 1., with_ucmh = False):
	min_z = 100
	max_z = 2.5e6
	t = np.logspace(np.log10(min_z), np.log10(max_z),1000)
	function = lambda z: dQdz(model.f_gamma, model.mDM, model.sigma_v, z,spectral_index, with_ucmh)*(np.exp(-tau(z))/H_z(z))
	s = function(t)
	plt.loglog(t, s, 'b-', lw=3)
	plt.xlabel(r'$z$')
	plt.ylabel(r'$(1+z)G d \varepsilon / dz$')
	plt.title(r'Energy injection from dark matter annihilation')
	#plt.grid(True, which="both")
	plt.show()

def plot_mu_contours():
	print "Processing things you know..",
	import matplotlib
	import matplotlib.cm as cm
	import matplotlib.mlab as mlab
	matplotlib.rcParams['xtick.direction'] = 'out'
	matplotlib.rcParams['ytick.direction'] = 'out'
	n_range = np.linspace(1.,1.25,100)
	mDM_range = np.logspace(1,6,100)
	X = n_range
	Y = mDM_range
	mu = []
	for mDM in Y:
		for n in X:
			print ".",
			mu.append(np.log10(mu0_n_mDM(mDM,n)))
	X, Y = np.meshgrid(n_range, mDM_range)
	mu = np.asarray(mu).reshape(100,100)
	plt.figure()
	CS = plt.contourf(X,Y, mu)
	plt.clabel(CS, inline = 1, fontsize = 3)
	CB = plt.colorbar(CS, shrink=0.8, extend='both')
	plt.yscale('log')
	plt.xlabel(r'Spectral Index $n$')
	plt.ylabel(r'Wimp mass in $GeV$')
	plt.title(r'Contour of $\log_{10}(\mu)$')
	plt.show()
	
	
# Parser setup
parser = argparse.ArgumentParser(description='Passing some arguments')

parser.add_argument('-d','--print_density_squared_plot', 
					help="Print plot", action='store_true')

parser.add_argument('-m','--print_mu_to_n', 
					help="Print plot of \mu type distortion for range of spectral index", 
					action='store_true')

parser.add_argument('-c','--print_contours_mu_to_n_mDM', 
					help="Print contour maps of mu distortions as function of n and mDM", 
					action='store_true')
					
parser.add_argument('-e','--print_energy_injection',nargs='+', 
					help="Print energy injection plot for a WIMP with the given mass in GeV")
parser.add_argument('-n','--spectral_index',nargs='+', 
					help="Setting spectral index n")

					
args = parser.parse_args()	
if __name__ == "__main__":
	main(args)
