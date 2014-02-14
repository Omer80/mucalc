"""
Omer Tzuk, February 2014
Version 0.1
This script imports data from PPPC 4 DM ID files in order to reproduce 
the precentage of energy that is injected in form of electromagnetic interacting particles
"""
import numpy as np
import PPPC4DMID_Reader as pppc
import mssm_data_Reader as mssm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from astropy.table import Table
import argparse
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def load_mssm_data():
	data = mssm.Read_pp_log_mssm_data("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR")
	data_dict = {}
	data_dict["mass"] = data.data["m_{\chi_1^0} (GeV)"]
	data_dict["sigma_v"] = data.data["<\sigma v> (cm^3 s^{-1})"]
	data_dict["sigma_v_Z"] = data.data["<\sigma v> (Z0 Z0) (cm^3 s^{-1})"]
	data_dict["sigma_v_W"] = data.data["<\sigma v> (W+ W-) (cm^3 s^{-1})"]
	data_dict["sigma_v_tau"] = data.data["<\sigma v> (tau+ tau-) (cm^3 s^{-1})"]
	data_dict["sigma_v_b"] = data.data["<\sigma v> (b bar) (cm^3 s^{-1})"]
	data_dict["sigma_v_gam"] = data.data["<\sigma v> (gam gam) (cm^3 s^{-1})"]
	data_dict["sigma_v_Z_gam"] = data.data["<\sigma v> (Z gam) (cm^3 s^{-1})"]
	return data_dict

def load_pppc_data(masses, channels):	
	final_SM_nu_e = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
	final_SM_nu_mu = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
	final_SM_nu_tau = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
	#final_SM_positrons = pppc.ReadPPPC4DMID_data("data/AtProduction_positrons.dat")
	#final_SM_gammas = pppc.ReadPPPC4DMID_data("data/AtProduction_gammas.dat")
	#final_SM_antiprotons = pppc.ReadPPPC4DMID_data("data/AtProduction_antiprotons.dat")
	#final_SM_antideuterons = pppc.ReadPPPC4DMID_data("data/AtProduction_antideuterons.dat")
	
	#print masses
	#print channels
	table = Table([masses],names=('mDM',), meta={'name': 'f_gamma per mDM'})
	#print table
	
	for channel in channels:
		print channel
		nu_e_int = final_SM_nu_e.interp_integrated_column(channel)
		nu_mu_int = final_SM_nu_mu.interp_integrated_column(channel)
		nu_tau_int = final_SM_nu_tau.interp_integrated_column(channel)
		#positrons_int = final_SM_positrons.interp_integrated_column(channel)
		#gammas_int = final_SM_gammas.interp_integrated_column(channel)
		#antiprotons_int = final_SM_antiprotons.interp_integrated_column(channel)
		#antideuterons_int = final_SM_antideuterons.interp_integrated_column(channel)
		
		nu_fraction = []
		
		for mass in masses:
			#print mass				
			nu_e_int_for_mass = float(2 * nu_e_int(mass))
			#total_energy_for_channel = nu_e_int_for_mass
			nu_mu_int_for_mass = float(2 * nu_mu_int(mass))
			#total_energy_for_channel = total_energy_for_channel + nu_mu_int_for_mass
			nu_tau_int_for_mass = float(2 * nu_tau_int(mass))
			#total_energy_for_channel = total_energy_for_channel + nu_tau_int_for_mass
			#positrons_int_for_mass = float(2 * positrons_int(mass))
			#total_energy_for_channel = total_energy_for_channel + positrons_int_for_mass
			#gammas_int_for_mass = float(gammas_int(mass))
			#total_energy_for_channel = total_energy_for_channel + gammas_int_for_mass
			#antiprotons_int_for_mass = float(2 * antiprotons_int(mass))
			#total_energy_for_channel = total_energy_for_channel + antiprotons_int_for_mass
			#antideuterons_int_for_mass = float(2 * antideuterons_int(mass))
			#total_energy_for_channel = total_energy_for_channel + antideuterons_int_for_mass
			
			#print "Total energy for channel", channel, "is", total_energy_for_channel
			
			nu_energy = (nu_e_int_for_mass + nu_mu_int_for_mass + nu_tau_int_for_mass)
			#d_plus_p_energy = (antiprotons_int_for_mass + antideuterons_int_for_mass)
			#gamma_energy = gammas_int_for_mass
			#e_energy = positrons_int_for_mass
			
			nu_fraction.append(nu_energy/(2 * mass))
			
			#print [nu_energy , d_plus_p_energy, e_energy , gamma_energy]
			
			
		table[channel] = np.array(nu_fraction)
			
			
			
			
				
			
	return table
	
	
def main(channels = ['eL'], mass = 200, save_figures = False, detailed_plots = False):
	br_data = load_mssm_data()
	table = load_pppc_data(br_data["mass"],['Z','W','t','b','\[Gamma]'])
	table.write("table.tex", format='latex')
	
	
	
	
	


# Parser setup
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-c','--channels', nargs='+',
                   help='Choose specific channel for plotting')
                   
parser.add_argument('-a','--all_channels', help="Producing pie charts for all channels",
					action='store_true')
parser.add_argument('-d','--detailed_plots', help="Producing pie charts",
					action='store_true')
parser.add_argument('-s','--save_figures', help="Saves the pie charts to files",
					action='store_true')
parser.add_argument('-m','--mass', help="Specifing the mass for the calculations",
					 action='count', default = 200)
parser.add_argument('-p','--print_possible_channels', help="Prints all possible channels for calculation, exit",
					action='store_true')

args = parser.parse_args()	
if __name__ == "__main__":
	if args.print_possible_channels:
		with open("data/AtProduction_neutrinos_e.dat",'r') as definitions:
			definitions = definitions.readline()
			definitions = definitions.split()
			print definitions
			
	else:	
		if args.all_channels:
			with open("data/AtProduction_neutrinos_e.dat",'r') as definitions:
				definitions = definitions.readline()
				definitions = definitions.split()
			channels = definitions[2:]
		else:
			if args.channels:
				channels = args.channels
			else:
				channels = ['eL']
		
		if args.mass != None:
			mass = int(args.mass)
		else:
			mass = 200	
		main(channels, mass, args.save_figures, args.detailed_plots)
	
