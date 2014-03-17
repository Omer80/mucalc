"""
Omer Tzuk, March 2014
Version 0.1
Based on arXiv:09084082v5 this code is building the function related to Ultracompact Dark Matter Minihalos
Those equations will be used to estimate the amplification in WIMP annihilation at the early Universe
"""
import numpy as np
from PPPC4DMID_Reader import *
from mssm_data_Reader import *
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.table import Table
from astropy import units as u
from astropy import constants as const
from astropy import cosmology
import argparse
from progressbar import Bar, ETA, Percentage, ProgressBar
import mucalc
import h5py
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

cosmology.core.set_current(cosmology.Planck13)

# delta_m as calculated from eq. (2) at arXiv:09084082v5
delta_m_EW  = 5.6e-19*u.solMass
delta_m_QCD = 1.1e-9*u.solMass
delta_m_ee  = 0.33*u.solMass
z_eq = (2.32e4 *  cosmology.Planck13.Om(0) *(cosmology.H(0).value/100)**2) -1
f_chi = 0.834 # as in arXiv:09084082v5

def M_h(z):
	'''(float) -> float in units.solMass (from astropy)
	M_h is the total mass of the UCMH at redshift z
	equation (1) at arXiv:09084082v5
	'''
	M_h = delta_m * ((1+z_eq)/(1+z))
	return M_h
	
def R_h(z):
	'''(float) -> float in units.pc
	R_h is the maximum extent of the UCMH at redshift z
	as in eq. (6) at arXiv:09084082v5
	'''
	R_h = 0.019 * (1000/(z + 1)) * ((M_h(z).value)**(1/3))
	return R_h

rho_max = lambda mDM, sigma_v, t_i, t : mDM / (sigma_v *(t * t_i))

def r_cut():
	''' 
	'''
	
		
def rho_chi(r,z):
	'''(float,float) -> float
	rho_chi as in eq. (5) at arXiv:09084082v5
	'''
	if r < r_cut:
		rho_chi = rho_max
	else:
		rho_chi = (3*f_chi*M_h(z).to(u.g)) / (16*np.PI*(R_h(z)**(3/4))*(r**(9/4)))
	return rho_chi


def main():
	delta_m = delta_m_ee
	a = 1 * u.solMass
	print a.to(u.g)
	
	
# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')

args = parser.parse_args()	
if __name__ == "__main__":
	main()
