"""
Omer Tzuk, February 2014
Version 0.1
This script imports data from PPPC 4 DM ID files in order to reproduce 
the precentage of energy that is injected in form of electromagnetic interacting particles
and the mu-type distortion for every WIMP candidate, and create a table which is based on 
pp_log_mssm1_cdm_mup_CTA_sigmavXBR, with added two columns of f_gamma and \mu type distortion magnitude
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

def load_mssm_data():
	mssm = Read_pp_log_mssm_data("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR")
	return mssm.data

def produce_table(mssm_data, channels):	
	final_SM_nu_e = ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
	final_SM_nu_mu = ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
	final_SM_nu_tau = ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
	
	f_gamma = []
	#i = 0
	column_length = len(mssm_data)
	print "Calculating f_gamma for each row.",
	#progress_bar = ProgressBar(widgets = ['Progress: ', Percentage(), ' ',
	                                      #Bar(marker='X'), ' ', ETA(), ' ']).start()
	for row in 	mssm_data:
		mass = row["m_{\chi_1^0} (GeV)"]
		print ".",

		nu_fractions = {}
		for channel in channels:
			#print channel
			nu_e_int = final_SM_nu_e.interp_integrated_column(channel)
			nu_mu_int = final_SM_nu_mu.interp_integrated_column(channel)
			nu_tau_int = final_SM_nu_tau.interp_integrated_column(channel)
			
			nu_fraction = []
		
			nu_e_int_for_mass = float(2 * nu_e_int(mass))
			nu_mu_int_for_mass = float(2 * nu_mu_int(mass))
			nu_tau_int_for_mass = float(2 * nu_tau_int(mass))
			
			
			nu_energy = (nu_e_int_for_mass + nu_mu_int_for_mass + nu_tau_int_for_mass)
			
			nu_fractions[channel] = nu_energy / (2 * mass)
			
		f_gamma_for_mass = 	calc_f_gamma(row, nu_fractions)
		#print "mDM", mass, "f_gamma ",f_gamma_for_mass
		f_gamma.append(f_gamma_for_mass)
		#i = i+1
		#progress_bar.update(i/column_length)
			
	#progress_bar.finish()
	print " "
	mssm_data["f_gamma"] = np.array(f_gamma)				
			
	return mssm_data

def calc_f_gamma(row, nu_fractions):
	#print nu_fractions
	mass = row["m_{\chi_1^0} (GeV)"]
	sigma_v = row["<\sigma v> (cm^3 s^{-1})"]
	br_Z = row["<\sigma v> (Z0 Z0) (cm^3 s^{-1})"] / sigma_v
	br_W = row["<\sigma v> (W+ W-) (cm^3 s^{-1})"]/ sigma_v
	br_tau = row["<\sigma v> (tau+ tau-) (cm^3 s^{-1})"]/ sigma_v
	br_b = row["<\sigma v> (b bar) (cm^3 s^{-1})"]/ sigma_v
	br_gam = row["<\sigma v> (gam gam) (cm^3 s^{-1})"]/ sigma_v
	br_Z_gam = row["<\sigma v> (Z gam) (cm^3 s^{-1})"]/ sigma_v
	total =  br_Z+br_W+br_tau+br_b+br_gam
	nu_Z = nu_fractions['Z']
	nu_W = nu_fractions['W']
	nu_tau = nu_fractions['t']
	nu_b = nu_fractions['b']
	nu_gam = nu_fractions['\[Gamma]']
	#print "Each br",br_Z,br_W,br_tau,br_b,br_gam
	#print "total", total
	f_nu = 	(br_Z*nu_Z + br_W*nu_W + br_tau*nu_tau + br_b*nu_b + br_gam*nu_gam)/total
	#print "f_gamma", (1. - f_nu)
	return (1. - f_nu)
	
		
def main():
	mssm_data = load_mssm_data()
	table = produce_table(mssm_data,['Z','W','t','b','\[Gamma]'])
	#print table
	table.write("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR.hdf5", format='hdf5',path='data')	


# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')

args = parser.parse_args()	
if __name__ == "__main__":
	main()
	
