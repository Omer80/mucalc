import numpy as np
import scipy.integrate
from scipy.interpolate import interp1d
import PPPC4DMID_Reader as pppc
from astropy.table import Table


class Read_pp_log_mssm_data(object):
	'''
	This class is built to handle data from the files 
	data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR.info
	data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR.txt
	'''
	def __init__(self, filename):
		self.load_info(filename)
		self.load_data(filename)
		
	def load_info(self, filename):			
		info_file = filename+".info"
		#"data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR.info"
		
		with open(info_file,'r') as lines:
					lines = lines.readlines()
					#definitions = definitions.split()
					#index_labs = lines.index("lab1")
					#print index_lab
					self.columns_defs = {}
					for line in lines:
						if line.find("lab") == 0:
							zero_index = line.find("=")
							self.columns_defs[int(line[3:zero_index])+2] = line[zero_index+2:-1]
					#labs = [s for s in lines if "lab" in s]
					#print self.columns_defs
					
	
	def load_data(self, filename):
		txt_file = filename+".txt"
		#data_file = "data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR.txt"
		table = np.loadtxt(txt_file, skiprows=0 , unpack=True)
		
		#print len(table)
		#print len(table[0])
		
		#self.data = {}
		
		self.data = Table([table[0],table[1]],names=("multip","chisq"), meta={'name': 'mssm data table'})
		#self.data["multip"] = table[0]
		#self.data["chisq"] = table[1]
		
		for key in self.columns_defs:
			#print key, self.columns_defs[key]
			self.data[self.columns_defs[key]] = table[key-1]
		
if __name__ == "__main__":
	data = Read_pp_log_mssm_data("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR")
	print data.data["<\sigma v> (Z0 Z0) (cm^3 s^{-1})"][0:10]
