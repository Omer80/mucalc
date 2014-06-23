"""
Omer Tzuk, February 2014
Version 0.1

"""
import numpy as np
from PPPC4DMID_Reader import *
from mssm_data_Reader import *
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.table import Table
import argparse
from progressbar import Bar, ETA, Percentage, ProgressBar
import mucalc
import h5py
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

from model import Wimp_model


def produce_matrix(mass_range, cross_section_range):
	final_SM_nu_e = ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
	final_SM_nu_mu = ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
	final_SM_nu_tau = ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
	
	nu_e_int = final_SM_nu_e.interp_integrated_column('b')
	nu_mu_int = final_SM_nu_mu.interp_integrated_column('b')
	nu_tau_int = final_SM_nu_tau.interp_integrated_column('b')
	
	X = mass_range
	
	Y = cross_section_range
	
	mu = []
	log_mu = []
	detection_levels = []
	
	for sigma_v in Y:
		for mDM in X:
			print ".",
			nu_e_int_for_mass = float(2 * nu_e_int(mDM))
			nu_mu_int_for_mass = float(2 * nu_mu_int(mDM))
			nu_tau_int_for_mass = float(2 * nu_tau_int(mDM))
						
			nu_energy = (nu_e_int_for_mass + nu_mu_int_for_mass + nu_tau_int_for_mass)
			
			f_gamma = 1. - (nu_energy / (2. * mDM))
			WIMP = Wimp_model(mDM, f_gamma ,sigma_v)		
			mu.append(mucalc.wimp_mucalc(WIMP))
			log_mu.append(np.log10(mucalc.wimp_mucalc(WIMP)))
			detection_levels.append(mucalc.WIMP_detection(WIMP))
			
	X, Y = np.meshgrid(mass_range,cross_section_range)
	
	mu = np.asarray(mu).reshape(3,3)	
	log_mu = np.asarray(mu).reshape(3,3)	
	detection_levels = np.asarray(mu).reshape(3,3)	

			
	
	import matplotlib as mpl
	import matplotlib.cm as cm
	import matplotlib.mlab as mlab
	mpl.rcParams['xtick.direction'] = 'out'
	mpl.rcParams['ytick.direction'] = 'out'
	
	plt.figure()
	
	#cmap = mpl.colors.ListedColormap(['w', 'c', 'b'])
	#levels = [np.log10(5e-8),np.log10(9e-5)]
	#CS = plt.contourf(X, Y,log_mu , levels,
                        #colors = ('w', 'c', 'b'),
                        #extend='both')
	#proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) for pc in CS.collections]
	#plt.legend(proxy, [r"$\mu < 5 \times 10^{-8}$", "Detectable by future experiments", "Detectable by COBE FIRAS"])
	
	CS = plt.contour(X, Y,log_mu)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel(r'')
	plt.ylabel(r'')
	plt.show()


		
def main(args):
	
	mass_range = np.logspace(0.7,3.3,3)
	cross_section_range = np.logspace(-27,-23,3)
	produce_matrix(mass_range, cross_section_range)


# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')

parser.add_argument('-n','--spectral_index',nargs='+', 
					help="Setting spectral index n")

parser.add_argument('--T_EW', 
					help="Calculate results for phase transitions at the electroweak scale", 
					action='store_true')

parser.add_argument('--T_QCD', 
					help="Calculate results for phase transitions at the QCD", 
					action='store_true')

args = parser.parse_args()	
if __name__ == "__main__":
	main(args)
