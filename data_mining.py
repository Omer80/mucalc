import PPPC4DMID_Reader as pppc
import mssm_data_Reader as mssm

mssm = mssm.Read_pp_log_mssm_data("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR")

print mssm.data["m_1 (GeV)"][0:10]



		
