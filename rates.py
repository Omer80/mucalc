"""
Omer Tzuk, February 2014
Calculation of K_C, K_dC, K_br
rates of three dominant physical processes between photons and matter particle
in the early Universe
Version 0.1
"""
import numpy as np
from astropy import units as u
from astropy import cosmology
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

Y_He_term = ((1-(const.fmHe/2))/0.88)
Omega_b_term = ((const.Omega_b*(const.h0**2))/0.0226)

def K_C(z):
	K_C_z = (2.045e-30)*((1+z)**4)*Omega_b_term*Y_He_term
	return K_C_z
def K_dC(z):
	K_dC_z = (7.561e-41)*((1+z)**5)*(const.g_dC)*Omega_b_term*Y_He_term
	return K_dC_z
def K_br(z):
	K_br_z = (2.074e-27)*((1+z)**(5./2.))*(const.g_br)*(Omega_b_term**2.)*Y_He_term
	return K_br_z

H = lambda z : cosmology.H(z).to(1/u.s).value

def main():
	#z = 2e6
	#print "K_C(z)", K_C(z)
	#print "K_dC(z)", K_dC(z)
	#print "K_br(z)", K_br(z)
	#print "H(z)", H(z)
	plot_rates()

def plot_rates():
	z = np.logspace(3,7,10)
	K_C_z = K_C(z)
	K_dC_z = K_dC(z)
	K_br_z = K_br(z)
	H_z = H(z)
	plt.xscale('log')
	plt.yscale('log')
	plt.plot(z,K_C_z , 'm--', lw=1, label = r"Compton scattering rate")
	plt.plot(z,K_dC_z , 'b--', lw=1, label = r"Double Compton scattering rate")
	plt.plot(z,K_br_z , 'g--', lw=1, label = r"Bremsstrahlung rate")
	plt.plot(z, H_z, 'r-', lw=2, label = r"Hubble expansion rate")
	
	plt.xlabel(r'redshift $z$')
	plt.ylabel(r'rate $(s^{-1})$')
	plt.title(r'Rates of physical processes in respect to Hubble rate')
	plt.legend(loc = 4)
	plt.grid(True, which="both")
	plt.show()	
	
# Parser setup
parser = argparse.ArgumentParser(description='Passing some arguments')


					
args = parser.parse_args()	
if __name__ == "__main__":
	main()
