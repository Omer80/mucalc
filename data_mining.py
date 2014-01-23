import PPPC4DMID_Reader as pppc
import mssm_data_Reader as mssm

#mssm = mssm.Read_pp_log_mssm_data("data/pp_log_mssm1_cdm_mup_CTA_sigmavXBR")

#print mssm.data["m_1 (GeV)"][0:10]

final_SM_nu_e = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
final_SM_nu_mu = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
final_SM_nu_tau = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
final_SM_positrons = pppc.ReadPPPC4DMID_data("data/AtProduction_positrons.dat")
final_SM_gammas = pppc.ReadPPPC4DMID_data("data/AtProduction_gammas.dat")
final_SM_antiprotons = pppc.ReadPPPC4DMID_data("data/AtProduction_antiprotons.dat")
final_SM_antideuterons = pppc.ReadPPPC4DMID_data("data/AtProduction_antideuterons.dat")

nu_e_int = final_SM_nu_e.interp_integrated_column('b')
nu_mu_int = final_SM_nu_mu.interp_integrated_column('b')
nu_tau_int = final_SM_nu_tau.interp_integrated_column('b')
positrons_int = final_SM_positrons.interp_integrated_column('b')
gammas_int = final_SM_gammas.interp_integrated_column('b')
antiprotons_int = final_SM_antiprotons.interp_integrated_column('b')
antideuterons_int = final_SM_antideuterons.interp_integrated_column('b')


total_energy_for_b_channel = nu_e_int(50)
total_energy_for_b_channel = total_energy_for_b_channel + nu_mu_int(50)
total_energy_for_b_channel = total_energy_for_b_channel + nu_tau_int(50)
total_energy_for_b_channel = total_energy_for_b_channel + positrons_int(50)
total_energy_for_b_channel = total_energy_for_b_channel + gammas_int(50)
total_energy_for_b_channel = total_energy_for_b_channel + antiprotons_int(50)
total_energy_for_b_channel = total_energy_for_b_channel + antideuterons_int(50)

print "Total energy for mDM 50 is", total_energy_for_b_channel


		
