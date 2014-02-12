"""
Omer Tzuk, February 2014
Version 0.1
This script imports data from PPPC 4 DM ID files in order to reproduce 
the precentage of energy that is injected in form of electromagnetic interacting particles
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
	pass


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
	
