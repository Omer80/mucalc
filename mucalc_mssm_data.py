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

from model import Wimp_model

def load_mssm_data(mssm_filename):
	mssm = Read_pp_log_mssm_data("data/"+mssm_filename)
	return mssm.data

def produce_table(mssm_data, channels, spectral_index = 1., with_ucmh = False):	
	final_SM_nu_e = ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
	final_SM_nu_mu = ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
	final_SM_nu_tau = ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
	
	f_gamma = []
	mu = []
	i = 0
	column_length = len(mssm_data)
	print "number of mssm models", column_length
	print "Calculating f_gamma for each row.",
	progress_bar = ProgressBar(widgets = ['Progress: ', Percentage(), ' ',
	                                      Bar(marker='X'), ' ', ETA(), ' ']).start()
	for row in 	mssm_data:
		mass = row["m_{\chi_1^0} (GeV)"]
		sigma_v = row["<\sigma v> (cm^3 s^{-1})"]
		log10_sigma_v = np.log10(row["<\sigma v> (cm^3 s^{-1})"])
		row["<\sigma v> (cm^3 s^{-1})"] = log10_sigma_v
		#print ".",

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
		WIMP = Wimp_model(mass,f_gamma_for_mass,sigma_v)
		mu_distortion_for_mass = mucalc.wimp_mucalc(WIMP,spectral_index,with_ucmh)
		#print "mDM", mass, "f_gamma ",f_gamma_for_mass
		f_gamma.append(f_gamma_for_mass)
		mu.append(np.log10(mu_distortion_for_mass))
		i = i+1
		progress_bar.update(i/column_length)
			
	progress_bar.finish()
	print " "
	mssm_data["f_{\gamma}"] = np.array(f_gamma)
	mssm_data["\mu distortion"]	 = 	np.array(mu)
	mssm_data["f_{\gamma}"].format = '%.19e'
	mssm_data["\mu distortion"].format = '%.19e'
	mssm_data["<\sigma v> (cm^3 s^{-1})"].format = '%.19e'
	mssm_data.format = '%.19e'			
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
	
		
def main(args):
	mssm_filename = "pp_log_mssm1_cdm_mup_CTA_sigmavXBR"
	
	n=1
	if args.T_EW:
		# Calculation for electroweak phase transition
		print 'Setting phase transition to be in the EW scale'
		mucalc.ucmh.T_phase_transition = mucalc.ucmh.T_EW
		mucalc.ucmh.M_H_z_X = mucalc.ucmh.M_H(mucalc.ucmh.T_phase_transition.value)
	if args.T_QCD:
		# Calculation for QCD phase transition
		print 'Setting phase transition to be in the QCD scale'
		mucalc.ucmh.T_phase_transition = mucalc.ucmh.T_QCD
		mucalc.ucmh.M_H_z_X = mucalc.ucmh.M_H(mucalc.ucmh.T_phase_transition.value)
	
	if args.spectral_index:
		with_ucmh = True
		n = args.spectral_index[0]
		print 'Calculating for UCMHs with small-scale spectral index n=', n
		mssm_data = load_mssm_data(mssm_filename)
		table = produce_table(mssm_data,['Z','W','t','b','\[Gamma]'], float(n) , with_ucmh)
	else:
		mssm_data = load_mssm_data(mssm_filename)
		table = produce_table(mssm_data,['Z','W','t','b','\[Gamma]'])		
	#print table["multip","chisq","m_{\chi_1^0} (GeV)","f_gamma","\mu distortion"]
	table.write("results/mu_f_gamma_with_"+mssm_filename+".txt", 
	            format='ascii.fixed_width_no_header', delimiter=' ')
	         #   format = "latex")
	


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
	
