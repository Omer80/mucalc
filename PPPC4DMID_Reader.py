import numpy as np
import scipy
from math import *

class Wimp(object):
	def __init__(self, mass, definitions):
		self.mass = mass
		self.dN_dx = {}
		self.num_columns = len(definitions)
		self.int_x_dN_dx = {}
		for definition in definitions:
			self.dN_dx[definition] = []
			self.int_x_dN_dx[definition] = []
			#print definition, "added"
		print "Wimp of mass", self.mass, "is added"
			
	def integration_of_columns():
		print "Integration of columns"
		for j in range(2, self.num_columns):
			print self.dN_dx[definitions[j]]
		
		
		
wimp_data = {}

def read_data():
	'''
	Read the data tables of PPPC 4 DM ID in order to calculate f_gamma
	for the calculation of mu distortion
	'''
	
	#Reading definitions from one of the files
	definitions = open("data/AtProduction_neutrinos_e.dat",'r')
	definitions = definitions.readline()
	definitions = definitions.split()
	
	print definitions

	nu_e = np.loadtxt("data/AtProduction_neutrinos_e.dat", skiprows=1 , unpack=True)
	#nu_mu = np.loadtxt("data/AtProduction_neutrinos_mu.dat", skiprows=1 , unpack=True)
	#nu_tau = np.loadtxt("data/AtProduction_neutrinos_tau.dat", skiprows=1 , unpack=True)
	
	nu_e_mass = nu_e[0]
	#nu_e_logx = nu_e[1]
	
	
	i=0
	print "Starting reading file"
	for mass in nu_e_mass:
		if (i <(len(nu_e_mass)-1)) and (mass == nu_e_mass[i+1]) :
			#Saving the columns into arrays
			if mass not in wimp_data:
				wimp_data[mass] = Wimp(mass, definitions)
			for j in range(1, len(definitions)):
				wimp_data[mass].dN_dx[definitions[j]].append(float(nu_e[j][i]))
				#print nu_e[j][i]
			#print mass, "Same mass", nu_e_mass[i+1]
			i = i+1
		else:
			for j in range(1, len(definitions)):
				wimp_data[mass].dN_dx[definitions[j]].append(float(nu_e[j][i]))
			#Integrating the columns and storing into Wimp class
			wimp_data[mass].integration_of_columns()					
			if i==(len(nu_e_mass)-1):
				print "End of file"
			else:
				#print mass, "Skipping to next mass"
				i = i+1
		
	
	#last_column = nu_e[len(nu_e)-1]
	#print last_column
	#print definitions
	
	
	
	
read_data()
#print wimp_data[6].dN_dx['Log[10,x]']
#print len(wimp_data[6].dN_dx['Log[10,x]'])
#print len(wimp_data[6].dN_dx['e'])
	
	
	

