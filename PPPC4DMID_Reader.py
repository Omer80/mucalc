import numpy as np
import scipy.integrate
from scipy.interpolate import interp1d


from math import *

class Wimp(object):
	def __init__(self, mass, definitions):
		self.mass = mass
		self.log_x = []
		self.dN_dx = {}
		self.definitions = definitions
		self.num_columns = len(definitions)
		self.int_x_dN_dx = {}
		for j in range(2, self.num_columns):
			self.dN_dx[definitions[j]] = []
			self.int_x_dN_dx[definitions[j]] = []
			#print definition, "added"
		#print "Wimp of mass", self.mass, "is added"
			
	def integration_of_columns(self):
		#print "Integration of columns"
		sum_integrations = 0.
		for j in range(2, self.num_columns):
			# Creating an array of the sample for integration \int E * dN/dx
			sampling_integrand = []
			for i in range(len(self.log_x)):
				sampling_integrand.append(self.mass*(10**(self.log_x[i]))*self.dN_dx[self.definitions[j]][i])
			# Integrating using scipy.integrate.simps
			integration_result = scipy.integrate.simps(sampling_integrand, self.log_x)
			#for char in self.definitions[j]:
				#if char in 'LRT':
					#sum_integrations = sum_integrations + integration_result
			#print "Integration results", integration_result
			# Storing the value of integration in the dictionary created int_x_dN_dx
			self.int_x_dN_dx[self.definitions[j]].append(integration_result)
			#print "The integration of column", self.definitions[j], "is", integration_result
		#difference = (2 * self.mass) - sum_integrations
		#print "For mDM ", self.mass, "the sum of the integration of columns is", sum_integrations
		#print "2* mDM - (sum of integrated columns) = ", difference
			
					
		

class ReadPPPC4DMID_data(object):
	def __init__(self, filename):
		'''
		Read the data tables of PPPC 4 DM ID in order to calculate f_gamma
		for the calculation of mu distortion
		'''
		
		self.wimp_data = {}
		
		#Reading definitions from one of the files
		with open(filename,'r') as definitions:
			definitions = definitions.readline()
			definitions = definitions.split()
			self.definitions = definitions
		
		#print self.definitions
		
		data_table = np.loadtxt(filename, skiprows=1 , unpack=True)
		
		self.mass_array = data_table[0]
		
		
		i=0
		print "Starting reading file", filename
		for mass in self.mass_array:
			if (i <(len(self.mass_array)-1)) and (mass == self.mass_array[i+1]) :
				#Saving the columns into arrays
				if mass not in self.wimp_data:
					self.wimp_data[mass] = Wimp(mass, self.definitions)
				self.wimp_data[mass].log_x.append(float(data_table[1][i]))
				for j in range(2, len(self.definitions)):
					self.wimp_data[mass].dN_dx[self.definitions[j]].append(float(data_table[j][i]))
					#print nu_e[j][i]
				#print mass, "Same mass", self.mass_array[i+1]
				i = i+1
			else:
				self.wimp_data[mass].log_x.append(float(data_table[1][i]))
				for j in range(2, len(self.definitions)):
					self.wimp_data[mass].dN_dx[self.definitions[j]].append(float(data_table[j][i]))
				#Integrating the columns and storing into Wimp class
				self.wimp_data[mass].integration_of_columns()					
				if i==(len(self.mass_array)-1):
					print "End of file", filename
				else:
					#print mass, "Skipping to next mass"
					i = i+1
		
		# Arranging the integrated values of the columns into arrays
		self.filtered_mass_array = []
		self.integrated_columns_values = {}
		for j in range(2, len(self.definitions)):
			self.integrated_columns_values[self.definitions[j]] = []
		for mass in self.mass_array:
			if mass not in self.filtered_mass_array:
				#print "Mass is not in filtered then.."
				self.filtered_mass_array.append(mass)
				for j in range(2, len(self.definitions)):
					#print "adding", self.definitions[j]
					value = self.wimp_data[mass].int_x_dN_dx[self.definitions[j]]
					#print mass, value
					self.integrated_columns_values[self.definitions[j]].append(value)
			else:
				#print "Mass is already in filtered"
				pass
		
		#print "len(filtered_mass)", len(self.filtered_mass_array) , "len(integrated_column)", len(self.integrated_columns_values['e'])
		
	def interp_integrated_column(self, definition):
		'''  (string) -> interpolated_function(mass)
		This functions output the interpolated value of the integrated columns
		as function of the mass of the WIMP (mDM)
		'''
		#print self.filtered_mass_array
		#print "len(filtered_mass)", len(self.filtered_mass_array) , "len(integrated_column)", len(self.integrated_columns_values[definition])
		#print self.integrated_columns_values[definition]
		x = np.asarray(self.filtered_mass_array)
		y = np.asarray(self.integrated_columns_values[definition]).squeeze()
		for i in range(len(x)):
			print x[i], y[i]
		interpolated_values = interp1d(x,y)
		return interpolated_values
	
		
	

	
	
	

