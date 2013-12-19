/*======  Calculation of mu-type and y-type distortion for a given SUSY model  ========= 
   Omer Tzuk, 12 November 2013
   cliffon@gmail.com
========================================================================================*/ 

/*======  Spectrum calculator  ========= 
   Choose RGE from the list below. SuSpect is included 
   in micrOMEGAs, to use another code define the path 
   to the corresponding package in lib/Makefile
=====================================*/ 
#define RGE  suspect
     /* choose 'suspect','isajet','softSusy','spheno'*/

/*=========   SUSY scenario  ==========
  One can define SUGRA, AMSB, EWSB (for low scale input). 
  By default the program reads SLHA data file 
=======================================*/
#define SUGRA 
//#define SUGRANUH
//#define AMSB 
//#define EWSB 

/*====== Modules ===============
   Keys to switch on 
   various modules of micrOMEGAs
================================*/

#define OBTAIN_LSP
#define OBTAIN_CROSS_SECTION
#define CALCULATION_OF_MU
#define TAKE_VALUES_FROM_LSP_OF_MICROMEGAS

#include "mucalc.h"
#include "../sources/micromegas.h"
#include "../sources/micromegas_aux.h"
#include "lib/pmodel.h"


#define SUGRAMODEL_(A) A ## SUGRA
#define SUGRAMODEL(A) SUGRAMODEL_(A)

#define SUGRANUHMODEL_(A) A ## SUGRAnuh
#define SUGRANUHMODEL(A) SUGRANUHMODEL_(A)

#define AMSBMODEL_(A) A ## AMSB
#define AMSBMODEL(A) AMSBMODEL_(A)

#define EWSBMODEL_(A) A ## EwsbMSSM
#define EWSBMODEL(A) EWSBMODEL_(A)

#define PRINTRGE_(A) printf(" Spectrum calculator is %s\n", #A)
#define PRINTRGE(A)  PRINTRGE_(A)


/*Physical constants */
double PI = acos(-1.);
double Sigma_T = 6.65246 * pow(10,-25); /* Compton scattering crossections */
// double Sigma_x = Sigma_T * 3/4 //1.+x//x^3 //2 x /1+x////1+2 x/-Log[1+2 x]/+1//2 x/ Log[1+2 x]-/1+3 x///1+2 x/^2/;
double G = 6.6742 * pow(10,-8);/*cgs*/
double a = 7.5658 * pow(10,-15);/*radiation constant cgs*/
double c = 2.99792458 * pow(10,10);/*cm/s*/
double parsec = 3.0856776 * pow(10,18);/*cm*/
double mparsec = parsec * pow(10,6);/*cm*/
double Alpha_fs= 1.0/137.036;/* Fine structure constant */
double kb=1.3806505 * pow(10,-16);/*Boltzmann cgs*/
double hb=1.05457168* pow(10,-27);/*Planck/2\[Pi] cgs*/
double b = (2.404)/(pow(PI,2))*pow((kb/(hb*c)),3); /*radiation density constant cgs*/
// double h=2 * \[Pi]*hb;
double me= 9.1093826 * pow(10,-28);/*mass of electron g*/
double ev=1.60217646* pow(10,-12);
double GeV= ev * pow(10,9);
/*Helium fraction Y=0.24*/
double Nnu = 3.046; /* Number of neutrinos */

/*Cosmological  Parameters */

double TCMB = 2.725; /* CMB Temperature in K */
double Omega_b = 4.8999 * pow(10,-2); /* Baryon density */
double Omega_cdm = 2.6709 * pow(10,-1); /* Cold Dark matter density */
double h0 = 0.6711; /* Hubble parameter */
double fmHe = 0.24; /*Helium mass fraction */
double zmin = 200.0; /* Minimum redshift */

double Omega_m = Omega_cdm + Omega_b; /* Total Matter density */
double H0 = h0*100*3.24077649*pow(10,-20) ;/*s^-1*/
double fnu = Nnu * (7.0/8) *  pow((4/11),(4/3));
double zeq = (3 *pow((H0 * c),2) * Omega_m)/(8 * PI * G * a * pow(TCMB,4) * (1+fnu))-1;
double Omega_Lambda = 0.6825;
double Omega_r = 1 - (Omega_Lambda + Omega_m);/*Omega_m/(1+zeq);*/

double Rho_cr = (3 * pow(H0,2)) / (8 * PI * G);


double zmax = 5.0 * pow(10,6); /* Maximum redshift */



/* Dark matter annihilation example */

double mdm = (10 * GeV)/pow(c,2); /* Dark matter mass - !!!!! This will be eventually a function taking value from micOMEGAs */
double f_Gamma = 1; /* Fraction of energy that goes into particles with electromagnetic interaction and is deposited in the plasma */
double sigma_v = (3 * pow(10,-27))/(Omega_cdm * pow(h0,2));

/* Chemical potential calculation parameters from arXiv: 1203.2601v2 */
double C = 0.7768;
double B = 1.803;
double z_dC = 1.96 * pow(10,6);
double z_br = 1.05 * pow(10,7);
double z_eps = 3.67 * pow(10,5);
double eps = 0.0151;
double z_dC_prime = 7.11 * pow(10,6);
double z_br_prime = 5.41 * pow(10,11);

int main(int argc,char** argv)
{
	int err;
	char cdmName[10];
	int spin2, charge3,cdim;

	ForceUG=0;  /* to Force Unitary Gauge assign 1 */
// sysTimeLim=1000; 
/*
	if you would like to work with superIso
	setenv("superIso","./superiso_v3.1",1);  
*/
#ifdef SUGRA
	{
		double m0,mhf,a0,tb;
		double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  sgn, gMHu,  gMHd,
  gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3;
         
  printf("\n========= mSUGRA scenario =====\n");
  PRINTRGE(RGE);

  if(argc<5) 
  { 
	  printf(" This program needs 4 parameters:\n"
			  "   m0      common scalar mass at GUT scale\n"
			  "   mhf     common gaugino mass at GUT scale\n"
			  "   a0      trilinear soft breaking parameter at GUT scale\n"
			  "   tb      tan(beta) \n");
	  printf(" Auxiliary parameters are:\n"
			  "   sgn     +/-1,  sign of Higgsino mass term (default 1)\n"    
			  "   Mtp     top quark pole mass\n"
			  "   MbMb    Mb(Mb) scale independent b-quark mass\n"
			  "   alfSMZ  strong coupling at MZ\n");
	  /*    printf("Example: ./main 70 250 -300 10\n");  */
	  printf("Example: ./main 120 500 -350 10 1 173.1 \n");
	  exit(1); 
  } else  
  {  double Mtp,MbMb,alfSMZ;
  sscanf(argv[1],"%lf",&m0);
  sscanf(argv[2],"%lf",&mhf);
  sscanf(argv[3],"%lf",&a0);
  sscanf(argv[4],"%lf",&tb);
  if(argc>5)sscanf(argv[5],"%lf",&sgn); else sgn=1;
  if(argc>6){ sscanf(argv[6],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
  if(argc>7){ sscanf(argv[7],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
  if(argc>8){ sscanf(argv[8],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
  }

  /*==== simulation of mSUGRA =====*/
  gMG1=mhf, gMG2=mhf,gMG3=mhf;
  gAl=a0,   gAt=a0,  gAb=a0;  gMHu=m0,  gMHd=m0;
  gMl2=m0,  gMl3=m0, gMr2=m0, gMr3=m0;
  gMq2=m0,  gMq3=m0, gMu2=m0, gMd2=m0, gMu3=m0, gMd3=m0;

  err= SUGRAMODEL(RGE) (tb,  
		  gMG1, gMG2, gMG3,  gAl,  gAt, gAb,  sgn, gMHu, gMHd,
    gMl2, gMl3, gMr2, gMr3, gMq2,  gMq3, gMu2, gMu3, gMd2, gMd3); 
	}
#elif defined(SUGRANUH)
	{
		double m0,mhf,a0,tb;
		double gMG1, gMG2, gMG3,  gAl, gAt, gAb,  gMl2, gMl3, gMr2, gMr3, gMq2, gMq3, gMu2, gMu3, gMd2, gMd3,mu,MA;
         
		printf("\n========= mSUGRA non-universal Higgs scenario =====\n");
		PRINTRGE(RGE);

		if(argc<7) 
		{ 
			printf(" This program needs 6 parameters:\n"
					"   m0      common scalar mass at GUT scale\n"
					"   mhf     common gaugino mass at GUT scale\n"
					"   a0      trilinear soft breaking parameter at GUT scale\n"
					"   tb      tan(beta) \n" 
					"   mu      mu(EWSB)\n"
					"   MA      mass of pseudoscalar Higgs\n");     
			printf(" Auxiliary parameters are:\n"
					"   Mtp     top quark pole mass\n"
					"   MbMb    Mb(Mb) scale independent b-quark mass\n"
					"   alfSMZ  strong coupling at MZ\n");
			/*    printf("Example: ./main 70 250 -300 10\n");  */
			printf("Example: ./main 120 500 -350 10 680 760  \n");
			exit(1); 
		} else  
		{  double Mtp,MbMb,alfSMZ;
		sscanf(argv[1],"%lf",&m0);
		sscanf(argv[2],"%lf",&mhf);
		sscanf(argv[3],"%lf",&a0);
		sscanf(argv[4],"%lf",&tb);
		sscanf(argv[5],"%lf",&mu);
		sscanf(argv[6],"%lf",&MA); 
		if(argc>7){ sscanf(argv[7],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
		if(argc>8){ sscanf(argv[8],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
		if(argc>9){ sscanf(argv[9],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
		}

		/*==== simulation of mSUGRA =====*/
		gMG1=mhf, gMG2=mhf,gMG3=mhf;
		gAl=a0,   gAt=a0,  gAb=a0;
		gMl2=m0,  gMl3=m0, gMr2=m0, gMr3=m0;
		gMq2=m0,  gMq3=m0, gMu2=m0, gMd2=m0, gMu3=m0, gMd3=m0;

		err= SUGRANUHMODEL(RGE) (tb,gMG1,gMG2,gMG3,gAl,gAt,gAb,gMl2,gMl3,gMr2,gMr3,gMq2,gMq3,gMu2,gMu3,gMd2,gMd3,mu,MA); 
	}
#elif defined(AMSB)
	{
		double m0,m32,sgn,tb;

		printf("\n========= AMSB scenario =====\n");
		PRINTRGE(RGE);
		if(argc<4) 
		{ 
			printf(" This program needs 3 parameters:\n"
					"   m0      common scalar mass at GUT scale\n"
					"   m3/2    gravitino mass\n"
					"   tb      tan(beta) \n");
			printf(" Auxiliary parameters are:\n"
					"   sgn     +/-1,  sign of Higgsino mass term (default 1)\n"    
					"   Mtp     top quark pole mass\n"
					"   MbMb    Mb(Mb) scale independent b-quark mass\n"
					"   alfSMZ  strong coupling at MZ\n");
			printf("Example: ./main 450  60000 10\n");                                                                          
			exit(1); 
		} else  
		{  double Mtp,MbMb,alfSMZ;
		sscanf(argv[1],"%lf",&m0);
		sscanf(argv[2],"%lf",&m32);
		sscanf(argv[3],"%lf",&tb);
		if(argc>4)sscanf(argv[4],"%lf",&sgn); else sgn=1;
		if(argc>5){ sscanf(argv[5],"%lf",&Mtp);    assignValW("Mtp",Mtp);      }
		if(argc>6){ sscanf(argv[6],"%lf",&MbMb);   assignValW("MbMb",MbMb);    }
		if(argc>7){ sscanf(argv[7],"%lf",&alfSMZ); assignValW("alfSMZ",alfSMZ);}
		}

		err= AMSBMODEL(RGE)(m0,m32,tb,sgn);
 
	}
#elif defined(EWSB)
	{ 
		printf("\n========= EWSB scale input =========\n");
		PRINTRGE(RGE);

		if(argc <2) 
		{  printf("The program needs one argument:the name of file with MSSM parameters.\n"
				"Example: ./main mssm1.par \n");
				exit(1);
		}  
   
		printf("Initial file  \"%s\"\n",argv[1]);
     
		err=readVarMSSM(argv[1]);
          
		if(err==-1)     { printf("Can not open the file\n"); exit(2);}
		else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(3);}

		err=EWSBMODEL(RGE)();
	}
#else
	{
		printf("\n========= SLHA file input =========\n");

		if(argc <2) 
		{  printf("The program needs one argument:the name of SLHA input file.\n"
				"Example: ./main suspect2_lha.out \n");
				exit(1);
		}  
   
		printf("Initial file  \"%s\"\n",argv[1]);
		err=lesHinput(argv[1]);
		if(err) exit(2);
	}
#endif
#ifdef OBTAIN_LSP          
	if(err==-1)     { printf("Can not open the file\n"); exit(2);}
	else if(err>0)  { printf("Wrong file contents at line %d\n",err);exit(3);}
  
	{ int nw;
	printf("Warnings from spectrum calculator:\n");
	nw=slhaWarnings(stdout);
	if(nw==0) printf(" .....none\n");
	} 

	if(err) exit(1);
	err=sortOddParticles(cdmName);
	if(err) { printf("Can't calculate %s\n",cdmName); return 1;}

	qNumbers(cdmName,&spin2, &charge3, &cdim);
	printf("\nDark matter candidate is '%s' with spin=%d/2  mass=%.2E\n",
	       cdmName,       spin2, Mcdm); 
  
	if(charge3) { printf("Dark Matter has electric charge %d/3\n",charge3); exit(1);}
	if(cdim!=1) { printf("Dark Matter is a color particle\n"); exit(1);}
	if(strcmp(cdmName,"~o1")) printf(" ~o1 is not CDM\n"); 
	else o1Contents(stdout);
#endif
#ifdef OBTAIN_CROSS_SECTION
{ 
  int err,i;
  double Emin=1,SMmev=320;/*Energy cut in GeV and solar potential in MV*/
  double  sigmaV;
  double vcs_gz,vcs_gg;
  char txt[100];
  double SpA[NZ],SpE[NZ],SpP[NZ];
  double FluxA[NZ],FluxE[NZ],FluxP[NZ];
  double SpNe[NZ],SpNm[NZ],SpNl[NZ];  
//  double * SpNe=NULL,*SpNm=NULL,*SpNl=NULL;
  double Etest=Mcdm/2;
 
/* default DarkSUSY parameters */

/*
    K_dif=0.036;
    L_dif=4;  
    Delta_dif=0.6; 
    Vc_dif=10;
    Rdisk=30;
    SMmev=320;
*/                        
  
printf("\n==== Indirect detection =======\n");  

  sigmaV=calcSpectrum( 2+4,SpA,SpE,SpP,SpNe,SpNm,SpNl ,&err);
    /* Returns sigma*v in cm^3/sec.     SpX - calculated spectra of annihilation.
       Use SpectdNdE(E, SpX) to calculate energy distribution in  1/GeV units.
       
       First parameter 1-includes W/Z polarization
                       2-includes gammas for 2->2+gamma
                       4-print cross sections             
    */
  printf("sigmav=%.2E[cm^3/s]\n",sigmaV);
  //sigma_v = Sigma_v(sigmaV);

  if(SpA)
  { 
     double fi=0.,dfi=M_PI/180.; /* angle of sight and 1/2 of cone angle in [rad] */ 
                                                   /* dfi corresponds to solid angle 1.E-3sr */                                             
     printf("Photon flux  for angle of sight f=%.2f[rad]\n"
     "and spherical region described by cone with angle %.4f[rad]\n",fi,2*dfi);
     gammaFluxTab(fi,dfi, sigmaV, SpA, FluxA);

#ifdef SHOWPLOTS
     sprintf(txt,"Photon flux for angle of sight %.2f[rad] and cone angle %.2f[rad]",fi,2*dfi);
     displaySpectrum(FluxA,txt,Emin,Mcdm,1);
#endif
     printf("Photon flux = %.2E[cm^2 s GeV]^{-1} for E=%.1f[GeV]\n",SpectdNdE(Etest, FluxA), Etest);       
#ifdef LoopGAMMA
     if(loopGamma(&vcs_gz,&vcs_gg)==0)
     {
         printf("Gamma  ray lines:\n");
         printf("E=%.2E[GeV]  vcs(Z,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm-91.19*91.19/4/Mcdm,vcs_gz,
                               gammaFlux(fi,dfi,vcs_gz));  
         printf("E=%.2E[GeV]  vcs(A,A)= %.2E[cm^3/s], flux=%.2E[cm^2 s]^{-1}\n",Mcdm,vcs_gg, 
                             2*gammaFlux(fi,dfi,vcs_gg));
     }
#endif     
  }

  if(SpE)
  { 
    posiFluxTab(Emin, sigmaV, SpE, FluxE);
    if(SMmev>0)  solarModulation(SMmev,0.0005,FluxE,FluxE);    
#ifdef SHOWPLOTS     
    displaySpectrum(FluxE,"positron flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("Positron flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxE),  Etest); 
  }
  
  if(SpP)
  {
    pbarFluxTab(Emin, sigmaV, SpP,  FluxP); 
    
    if(SMmev>0)  solarModulation(SMmev,1,FluxP,FluxP);     
#ifdef SHOWPLOTS    
     displaySpectrum(FluxP,"antiproton flux [cm^2 s sr GeV]^{-1}" ,Emin,Mcdm,1);
#endif
    printf("Antiproton flux  =  %.2E[cm^2 sr s GeV]^{-1} for E=%.1f[GeV] \n",
    SpectdNdE(Etest, FluxP),  Etest);     
  }
}  
#endif
                             
	
#ifdef CALCULATION_OF_MU
	{
#ifdef 	TAKE_VALUES_FROM_LSP_OF_MICROMEGAS
	{
	mdm = mdm_calc(Mcdm);
	
	}
#endif
	double z = pow(10,5);
	printf("***************** Values used for energy calculation ********************\n");
	printf("H0 - %.2e \n", H0);
	printf("zeq - %.2e \n", zeq);
	printf("Omega_m - %.2e \n", Omega_m);
	printf("Omega_r - %.2e \n", Omega_r);
	printf("Omega_Lambda - %.2e \n", Omega_Lambda);
	printf("Rho_cr - %.2e \n", Rho_cr);
	printf("mdm - %.2e \n", mdm);
	printf("The H_z - is %.2e \n",H(z));
	printf("The ndm_z is %.2e \n",ndm_z(z));
	printf("The sigma_v is %.2e \n",sigma_v);
	printf("The tau is %.2e \n",tau(z));
	printf("*************************************************************************\n");
	double value = dQdz(z) * ( exp(-tau(z))/ H(z));
	printf("For redshift %.2e ", z);
	printf("the energy injection is %.2e \n",value);
	double value_paper = dummy_energy_injection(z)* ( exp (-tau(z))/ H(z));
	printf("While the calculation as in the paper is %.2e \n", value_paper);
	double z_min = 5 * pow(10,4);
	double z_i = 6 * pow(10,6);
	int subdivisions = 1000;
	value = mu_0(z_i, z_min, subdivisions);
	printf("Mu at our time is  %.2e \n", value);
	}
#endif 
	return 0;
}


/* Integration functions */
extern double integrate_trapezoid(double (*integrand)(double),double min, double max, int subdivisions)
{
	double h = (max-min)/subdivisions;
	double s=0;
	for (int k=0; k < (subdivisions) ; k++)
	{
		s = s + (h/2) * ((*integrand)(min + k*h) + (*integrand)(min + (k+1)*h));
	}
	return s;
}

/*Romberg method integration */
// extern double integrate_romberg(double *function,double min, double max, double subdivisions)
// {
// 	double s[subdivisions];
// 	int i,k;
// 	double var ;
// 	for (i = 1; i< subdivisions; i++)
// 		s[i] = 1;
//  
// 	for (k=1; k< subdivisions ; k++)
// 	{
// 		for (i=1; i <=k; i++)
// 		{
// 			if (i==1)
// 			{
// 				var = s[i];
// 				s[i] = integrate_trapezoid(0, 1, pow(2, k-1));     // sub-routine trapeze 
// 			}                                       // integrated from 0 and 1
// 			/* pow() is the number of subdivisions*/
// 			else
// 			{
// 				s[k]= ( pow(4 , i-1)*s[i-1]-var )/(pow(4, i-1) - 1); 
//  
// 				var = s[i];
// 				s[i]= s[k];  
// 			}
// 		}
//  
// 		for (i=1; i <=k; i++)
// 			printf ("  %f  ", s[i]);
//  
// 		printf ("\n");
// 	}
//  
// 	return 0;
// }





/*Cosmological parameters functions */
extern double ndm_z(double z)
{
	double ndm = (Rho_cr * Omega_cdm / mdm) * pow((1+z),3);
	return ndm;
}

extern double H(double z)
{
	double H_z = H0 * pow((pow((1+z),3)*Omega_m + pow((1+z),4)*Omega_r + (1-(Omega_m+Omega_r))),0.5);
	return H_z;
}

/*WIMP cross section*/ 
extern double Sigma_v(double sigma_v)
{
	//double sigma_v = (3 * pow(10,-27))/(Omega_cdm * pow(h0,2));
	return sigma_v;// /(Omega_cdm * pow(h0,2));
}

extern int photons_per_annihilation()
{
	return 2;
}
/*Energy injection calculation from eq. (5.5) in arXiv: 1203.2601v2 */
extern double dQdz(double z) /* Energy injection rate */
{
	return f_Gamma * (mdm * pow(c,2) * pow(ndm_z(z),2) * sigma_v)/(a * pow((TCMB * (1+z)),4));
}

extern double dNdz(double z)
{
	return photons_per_annihilation() * (pow(ndm_z(z),2) * sigma_v)/(b * pow((TCMB * (1+z)),3) );
}
extern double dummy_energy_injection(double z) /* Energy injection for a fixed WIMP candidate as in arXiv: 1203.2601v2 */
{
	double epsilon_dot = 1.4 * pow(10,-29) * pow((1+z),2) * f_Gamma *(10/*GeV*/ / 10 /* WIMP mass in GeV*/) * ((Omega_cdm * pow(h0,2))/0.105 );
	return epsilon_dot;
}


/*Blackbody optical depth calculation from eq. (4.5) in arXiv: 1203.2601v2 */
extern double tau(double z)
{
	double first_part_tau_1, second_part_tau_1, tau_1,tau_2,tau_precise;
	first_part_tau_1 = pow(((1 + z)/(1 + z_dC)),5);
	second_part_tau_1 = pow(((1 + z)/(1 + z_br)),(5/2));
	tau_1 = pow((pow(((1 + z)/(1 + z_dC)),5) + pow(((1 + z)/(1 + z_br)),(5/2))),0.5);
	tau_2 = eps * log((pow(((1 + z)/(1 + z_eps)),(5/4)))+(pow((1 + pow(((1 + z)/(1 + z_eps)),(5/2))),(1/2))));
	tau_precise = (pow(((1+z)/(1+z_dC_prime)),3) + pow(((1+z)/(1+z_br_prime)),0.5));
	return 1.007 * (tau_1+tau_2) + tau_precise;
}

/*Chemical potential calculation from eq. (3.6) in arXiv: 1203.2601v2 */
extern double mu(double z)
{
	double first_part = 7.43 * pow(10,-5) * ((1 + z)/(2 * pow(10,6)));
	double second_part = 1.07 * pow(10,-6) * pow(((1 + z)/(2 * pow(10,6))), (-3/2));
	double value = pow((first_part + second_part),(1/2));
	
	return value;
}

extern double mu_integrand(double z)
{
	return (1/((1+z)*H(z)))* (dQdz(z)-(4/3)*(dNdz(z))) * exp(-tau(z));
}

extern double mu_0(double z_i, double z_min, int subdivisions)
{
	double first_part = mu(z_i) * exp(-tau(z_i));
	double second_part = C*B*integrate_trapezoid((&mu_integrand),z_min, z_i, subdivisions);
	return first_part+second_part;
}

extern double mdm_calc(double mass_in_GeV)
{
	return ( mass_in_GeV * GeV)/pow(c,2);
}

