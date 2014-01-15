import PPPC4DMID_Reader as pppc


nu_e = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
#nu_mu = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
#nu_tau = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")

e_int = nu_e.interp_integrated_column('e')
		
