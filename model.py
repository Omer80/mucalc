"""
Omer Tzuk, April 2014
Script for class Model
with the following parameters:
mDM - Mass of Wimp in GeV
f_gamma - ratio of energy that is injected as EM particles
sigma_v
with_ucmh - if True, Ultracompact dark matter Minihalos are considered
n - spectral index

"""
from mucalc_constants import *

class Wimp_model(object):
	def __init__(self, mDM, f_gamma, sigma_v):
		self.mDM = ( mDM * const.GeV)/(const.c**2)
		self.f_gamma = f_gamma
		self.sigma_v = sigma_v / (const.Omega_cdm * (const.h0**2)) 
	def __str__(self):
		print_WIMP = "Wimp of mass {0} with f_gamma {1} and cross-section {2}".format(str(self.mDM),str(self.f_gamma), str(self.sigma_v))
		return print_WIMP
