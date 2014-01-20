import numpy as np
import scipy.integrate
from scipy.interpolate import interp1d



#info_file = "data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR.info"

#with open(info_file,'r') as definitions:
			#definitions = definitions.readlines()
			##definitions = definitions.split()
			#print definitions
			
data_file = "data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR.txt"

data_table = np.loadtxt(data_file, skiprows=0 , unpack=True)

print data_table[0][10:]
