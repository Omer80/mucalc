"""
Omer Tzuk, February 2014
Version 0.6
This script imports data from PPPC 4 DM ID files and reproduce
Figure 4 from 1012.4515v4
"""

import PPPC4DMID_Reader as pppc
import mssm_data_Reader as mssm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def main(channels = ['eL'], mass = 200, save_figures = False, detailed_plots = False):
	
	print "Producing plots for channels:",channels
	print "For mass - ",mass
	
	
	channels = calculate_percentages_per_channel(mass, channels)
	#print channels
	
	for channel in channels:
		plot_channel(channel,channels[channel], mass, save_figures,detailed_plots)


def calculate_percentages_per_channel(mass, channels):	
	final_SM_nu_e = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
	final_SM_nu_mu = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
	final_SM_nu_tau = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
	final_SM_positrons = pppc.ReadPPPC4DMID_data("data/AtProduction_positrons.dat")
	final_SM_gammas = pppc.ReadPPPC4DMID_data("data/AtProduction_gammas.dat")
	final_SM_antiprotons = pppc.ReadPPPC4DMID_data("data/AtProduction_antiprotons.dat")
	final_SM_antideuterons = pppc.ReadPPPC4DMID_data("data/AtProduction_antideuterons.dat")
	
	#print channels
	channels_dict = {}
	
	for channel in channels:
		nu_e_int = final_SM_nu_e.interp_integrated_column(channel)
		nu_mu_int = final_SM_nu_mu.interp_integrated_column(channel)
		nu_tau_int = final_SM_nu_tau.interp_integrated_column(channel)
		positrons_int = final_SM_positrons.interp_integrated_column(channel)
		gammas_int = final_SM_gammas.interp_integrated_column(channel)
		antiprotons_int = final_SM_antiprotons.interp_integrated_column(channel)
		antideuterons_int = final_SM_antideuterons.interp_integrated_column(channel)
		
		nu_e_int_for_mass = nu_e_int(mass)
		total_energy_for_channel = nu_e_int_for_mass
		nu_mu_int_for_mass = nu_mu_int(mass)
		total_energy_for_channel = total_energy_for_channel + nu_mu_int_for_mass
		nu_tau_int_for_mass = nu_tau_int(mass)
		total_energy_for_channel = total_energy_for_channel + nu_tau_int_for_mass
		positrons_int_for_mass = 2 * positrons_int(mass)
		total_energy_for_channel = total_energy_for_channel + positrons_int_for_mass
		gammas_int_for_mass = gammas_int(mass)
		total_energy_for_channel = total_energy_for_channel + gammas_int_for_mass
		antiprotons_int_for_mass = antiprotons_int(mass)
		total_energy_for_channel = total_energy_for_channel + antiprotons_int_for_mass
		antideuterons_int_for_mass = antideuterons_int(mass)
		total_energy_for_channel = total_energy_for_channel + antideuterons_int_for_mass
		
		print "Total energy for channel", channel, "is", total_energy_for_channel
		
		nu_energy = (nu_e_int_for_mass + nu_mu_int_for_mass + nu_tau_int_for_mass)/ total_energy_for_channel
		d_plus_p_energy = (antiprotons_int_for_mass + antideuterons_int_for_mass)/ total_energy_for_channel
		gamma_energy = gammas_int_for_mass/ total_energy_for_channel
		e_energy = positrons_int_for_mass / total_energy_for_channel
		
		print [nu_energy , d_plus_p_energy, e_energy , gamma_energy]
		channels_dict[channel] = [nu_energy , d_plus_p_energy, e_energy , gamma_energy]
		
	return channels_dict
		
			
			

def plot_channel(channel,percentages, mass, save_figures,detailed_plots):
	# build a rectangle in axes coords
	left, width = .25, .5
	bottom, height = .25, .5
	right = left + width
	top = bottom + height
	
	channel = ''.join(i for i in channel if i in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
	nu_energy , d_plus_p_energy, e_energy, gamma_energy  = percentages  
	
	# The slices will be ordered and plotted counter-clockwise.
	if detailed_plots == False:
		labels = [r'$E_{\nu} / E_{tot}$ '  ,
	              r'$E_{d+p} / E_{tot}$ ' ,
	              r'$E_{e} / E_{tot}$ ' ,
	              r'$E_{\gamma} / E_{tot}$']
	else:
		labels = [r'$E_{\nu} / E_{tot}$ = '+str(nu_energy) ,
	              r'$E_{d+p} / E_{tot}$ = '+str(d_plus_p_energy),
	              r'$E_{e} / E_{tot}$ = '+str(e_energy),
	              r'$E_{\gamma} / E_{tot}$ = '+str(gamma_energy)]		
		
	sizes = [nu_energy , d_plus_p_energy, e_energy, gamma_energy]
	labels_a = [r'${\nu}$',r'${d+p}$',r'${e}$',r'${\gamma}$']
	colors = ['gold','red','green', 'lightskyblue']
	explode = (0.1, 0,0,0)
	patches, texts = plt.pie(sizes, explode = explode, colors=colors) 
	plt.legend(patches, labels, loc = "best" )
	E_gamma_e = e_energy + gamma_energy
	#plt.pie(sizes, explode=explode, labels=labels, colors=colors,
	        #autopct='%1.1f%%', shadow=True)
	# Set aspect ratio to be equal so that pie is drawn as a circle.
	plt.axis('equal')
	plt.title(r'DM DM $\rightarrow$ $%s$ + $%s$'%(channel,channel),position=(0.5,1),bbox=dict(facecolor='0.8',), fontsize=30)
	plt.text(-0.4,-0.76, r'$E_{\gamma + e} / E_{tot}$ = %.3f'%(E_gamma_e)
	         , bbox=dict(facecolor='white', alpha=0.5), fontsize=25) 
	plt.text(-0.4,-0.98, r'$E_{p + d} / E_{\gamma + e}$ = %.3f'%(d_plus_p_energy/E_gamma_e)
	         , bbox=dict(facecolor='white', alpha=0.5), fontsize=25)  
	          
	plt.tight_layout()
	
	if save_figures:
		plt.savefig("./figures/energy_distribution_for_channel_"+channel+".png")
	else:
		plt.show()
	plt.close()



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

args = parser.parse_args()	
if __name__ == "__main__":
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
