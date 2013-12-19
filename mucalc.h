#ifndef  __MIUCALC__
#define  __MIUCALC__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

// /*Physical constants */
// double PI;
// double Sigma_T; /* Compton scattering crossections */
// // double Sigma_x = Sigma_T * 3/4 //1.+x//x^3 //2 x /1+x////1+2 x/-Log[1+2 x]/+1//2 x/ Log[1+2 x]-/1+3 x///1+2 x/^2/;
// double G;/*cgs*/
// double a;/*radiation constant cgs*/
// double c;/*cm/s*/
// double parsec;/*cm*/
// double mparsec;/*cm*/
// double Alpha_fs;/* Fine structure constant */
// double kb;/*Boltzmann cgs*/
// double hb;/*Planck/2\[Pi] cgs*/
// // double h=2 * \[Pi]*hb;
// double me;/*mass of electron g*/
// double ev;
// double GeV;
// /*Helium fraction Y=0.24*/
// double Nnu; /* Number of neutrinos */
// 
// /*Cosmological  Parameters */
// 
// double TCMB; /* CMB Temperature in K */
// double Omega_b; /* Baryon density */
// double Omega_cdm; /* Cold Dark matter density */
// double h0; /* Hubble parameter */
// double fmHe; /*Helium mass fraction */
// double zmin; /* Minimum redshift */
// 
// double Omega_m; /* Total Matter density */
// double H0;/*s^-1*/
// double fnu;
// double zeq;
// double Omega_Lambda;
// double Omega_r;/*Omega_m/(1+zeq);*/
// 
// double Rho_cr;
// 
// 
// double zmax; /* Maximum redshift */
// 
// 
// 
// /* Dark matter annihilation example */
// 
// double mdm; /* Dark matter mass - !!!!! This will be eventually a function taking value from micOMEGAs */
// double f_Gamma; /* Fraction of energy that goes into particles with electromagnetic interaction and is deposited in the plasma */
// 
// /* Chemical potential calculation parameters from arXiv: 1203.2601v2 */
// double C;
// double B;
// double z_dC;
// double z_br;
// double z_eps;
// double eps;
// double z_dC_prime;
// double z_br_prime;

/* Integration functions */
extern double integrate_trapezoid(double (*integrand)(double),double min, double max, int subdivisions);

/*Romberg method integration */
// extern double integrate_romberg(double *function,double min, double max, double subdivisions);


/*Cosmological parameters functions */
extern double ndm_z(double z);
extern double H(double z);
extern int photons_per_annihilation();
/*Energy injection calculation from eq. (5.5) in arXiv: 1203.2601v2 */
extern double dQdz(double z); /* Energy injection rate */
extern double dNdz(double z); /* Photons injection rate */
extern double dummy_energy_injection(double z); /* Energy injection for a fixed WIMP candidate as in arXiv: 1203.2601v2 */
/*Blackbody optical depth calculation from eq. (4.5) in arXiv: 1203.2601v2 */
extern double tau(double z);
/*Chemical potential calculation from eq. (3.6) in arXiv: 1203.2601v2 */
extern double mu(double z);
extern double mu_integrand(double z);

extern double mu_0(double z_i, double z_min, int subdivisions);

/* Interface to micOMEGAs functions */
extern double mdm_calc(double mass_in_GeV);
extern double Sigma_v(double sigma_v);


/* Test functions */
extern double test_function(double z);



#endif
