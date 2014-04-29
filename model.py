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
	def __init__(self, mDM, f_gamma, sigma_v, spectral_index = 1. , with_ucmh = False):
		self.mDM = ( mDM * const.GeV)/(const.c**2)
		self.f_gamma = f_gamma
		self.sigma_v = sigma_v / (const.Omega_cdm * (const.h0**2)) 
		self.spectral_index = spectral_index
		self.with_ucmh = with_ucmh
