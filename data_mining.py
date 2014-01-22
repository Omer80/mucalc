import PPPC4DMID_Reader as pppc
import mssm_data_Reader as mssm

mssm = mssm.Read_pp_log_mssm_data("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR")

#print mssm.data["m_1 (GeV)"][0:10]

#final_SM_nu_e = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
#final_SM_nu_mu = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
#final_SM_nu_tau = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
#final_SM_positrons = pppc.ReadPPPC4DMID_data("data/AtProduction_positrons.dat")
#final_SM_gammas = pppc.ReadPPPC4DMID_data("data/AtProduction_gammas.dat")
#final_SM_antiprotons = pppc.ReadPPPC4DMID_data("data/AtProduction_antiprotons.dat")
#final_SM_antideuterons = pppc.ReadPPPC4DMID_data("data/AtProduction_antideuterons.dat")

#e_int = nu_e.interp_integrated_column('e')
		
