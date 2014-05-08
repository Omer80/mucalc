import const
import numpy as np

#Physical constants 
const.Sigma_T = 6.65246e-25 # Compton scattering cross-sections 
# const.Sigma_x = Sigma_T * 3/4 //1.+x//x^3 //2 x /1+x////1+2 x/-Log[1+2 x]/+1//2 x/ Log[1+2 x]-/1+3 x///1+2 x/^2/
const.G = 6.6742e-8#cgs
const.a = 7.5658e-15#radiation constant cgs
const.c = 2.99792458e10#cm/s
const.parsec = 3.0856776e18#cm
const.mparsec = const.parsec*10**6#cm
const.Alpha_fs= 1.0/137.036# Fine structure constant 
const.kb=1.3806505e-16#Boltzmann cgs
const.hb=1.05457168e-27#Planck/2\[Pi] cgs
const.b = (2.404)/(np.pi**2)*((const.kb/(const.hb*const.c))**3) #radiation density constant cgs
# const.h=2 * \[Pi]*hb
const.me= 9.1093826e-28#mass of electron g
const.ev=1.60217646e-12
const.GeV= const.ev * 10**9
#Helium fraction Y=0.24
const.Nnu = 3.046 # Number of neutrinos 

#Cosmological  Parameters 

const.TCMB = 2.725 # CMB Temperature in K 
const.Omega_b = 4.8999e-2 # Baryon density 
const.Omega_cdm = 2.6709e-1 # Cold Dark matter density 
const.h0 = 0.6711 # Hubble parameter 
const.fmHe = 0.24 #Helium mass fraction 
const.zmin = 200.0 # Minimum redshift 

const.Omega_m = const.Omega_cdm + const.Omega_b # Total Matter density 
const.H0 = const.h0*100*3.24077649*(10**(-20)) #s^-1
const.fnu = const.Nnu * (7.0/8) *  ((4/11)**(4/3))
const.zeq = (3 *((const.H0 * const.c)**2) * const.Omega_m)/(8 * np.pi * const.G * const.a * (const.TCMB**4) * (1+const.fnu))-1
const.Omega_Lambda = 0.6825
const.Omega_r = 1 - (const.Omega_Lambda + const.Omega_m)#Omega_m/(1+zeq)

const.Rho_cr = (3 * (const.H0**2)) / (8 * np.pi * const.G)


const.zmax = 5.0e6 # Maximum redshift 



# Dark matter annihilation example 

const.mdm = (10 * const.GeV)/(const.c**2) # Dark matter mass - !!!!! This will be eventually a function taking value from micOMEGAs 
const.f_Gamma = 1 # Fraction of energy that goes into particles with electromagnetic interaction and is deposited in the plasma 
const.sigma_v = (3e-27)/(const.Omega_cdm * (const.h0**2))

# Chemical potential calculation parameters from arXiv: 1203.2601v2 
const.C = 0.7768
const.B = 1.803
const.z_dC = 1.96e6
const.z_br = 1.05e7
const.z_eps = 3.67e5
const.eps = 0.0151
const.z_dC_prime = 7.11e6
const.z_br_prime = 5.41e11

# Rate coefficient calculations
const.g_dC = 1.005
const.g_br = 2.99
