import PPPC4DMID_Reader as pppc
import mssm_data_Reader as mssm

def loaddata():
	mssm = mssm.Read_pp_log_mssm_data("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR")
	return mssm

mssm = loaddata()
print mssm.data["<\sigma v> (Z0 Z0) (cm^3 s^{-1})"][0:10]



		
