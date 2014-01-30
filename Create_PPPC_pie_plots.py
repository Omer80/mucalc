"""
Omer Tzuk, January 2014
Version 0.1
This script imports data from PPPC 4 DM ID files and reproduce
Figure 4 from 1012.4515v4
"""

import PPPC4DMID_Reader as pppc
import mssm_data_Reader as mssm
import matplotlib.pyplot as plt

final_SM_nu_e = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_e.dat")
final_SM_nu_mu = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_mu.dat")
final_SM_nu_tau = pppc.ReadPPPC4DMID_data("data/AtProduction_neutrinos_tau.dat")
final_SM_positrons = pppc.ReadPPPC4DMID_data("data/AtProduction_positrons.dat")
final_SM_gammas = pppc.ReadPPPC4DMID_data("data/AtProduction_gammas.dat")
final_SM_antiprotons = pppc.ReadPPPC4DMID_data("data/AtProduction_antiprotons.dat")
final_SM_antideuterons = pppc.ReadPPPC4DMID_data("data/AtProduction_antideuterons.dat")

nu_e_int = final_SM_nu_e.interp_integrated_column('eL')
nu_mu_int = final_SM_nu_mu.interp_integrated_column('eL')
nu_tau_int = final_SM_nu_tau.interp_integrated_column('eL')
positrons_int = final_SM_positrons.interp_integrated_column('eL')
gammas_int = final_SM_gammas.interp_integrated_column('eL')
antiprotons_int = final_SM_antiprotons.interp_integrated_column('eL')
antideuterons_int = final_SM_antideuterons.interp_integrated_column('eL')


total_energy_for_b_channel = nu_e_int(200)
print nu_e_int(200)
total_energy_for_b_channel = total_energy_for_b_channel + nu_mu_int(200)
print nu_mu_int(200)
total_energy_for_b_channel = total_energy_for_b_channel + nu_tau_int(200)
print nu_tau_int(200)
total_energy_for_b_channel = total_energy_for_b_channel + positrons_int(200)
print positrons_int(200)
total_energy_for_b_channel = total_energy_for_b_channel + gammas_int(200)
print gammas_int(200)
total_energy_for_b_channel = total_energy_for_b_channel + antiprotons_int(200)
print antiprotons_int(200)
total_energy_for_b_channel = total_energy_for_b_channel + antideuterons_int(200)
print antideuterons_int(200)

gamma_plus_e = gammas_int(200) + positrons_int(200)



print "Total energy for mDM 200 is", total_energy_for_b_channel
print "The ratio is",gamma_plus_e / total_energy_for_b_channel

# The slices will be ordered and plotted counter-clockwise.
labels = 'gamma + e', 'total energy for eL channel'
sizes = [(total_energy_for_b_channel - gamma_plus_e)/total_energy_for_b_channel, gamma_plus_e/total_energy_for_b_channel]
colors = ['gold', 'lightskyblue']
explode = (0.1, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

plt.pie(sizes, explode=explode, labels=labels, colors=colors,
        autopct='%1.1f%%', shadow=True)
# Set aspect ratio to be equal so that pie is drawn as a circle.
plt.axis('equal')
plt.show()
		
