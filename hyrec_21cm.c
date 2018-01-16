/***********************************************************************
 High-redshift 21 cm module for Hyrec.

Computes the coefficients of the expansion (in Kelvins)
T_b = <T_b> + a \delta_H + b \delta_Tgas 
      + c \delta_H^2 + d \delta_Tgas^2 + e \delta_H \delta_Tgas,
where \delta_H = \delta n_H/n_H, \delta_Tgas = \delta(Tgas)/Tgas.
 
Does not yet account for Wouthuysen-Field coupling.
Written December 2014.
***********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "history.h"
#include "hyrec_params.h"

#define A10 2.85e-15  /* Spontaneous decay rate of the 21 cm transtion */
#define T10 0.068   /* Transition energy in Kelvins */
#define lambda21 21.1   /* 21cm wavelength */

/* Collisional couplings */
double kappa10_H(double Tgas) {
    return 3.1e-11 * exp(-32./Tgas + 0.357*log(Tgas));
}

double kappa10_e(double Tgas) {
  return 0.;
}

double kappa10_p(double Tgas) {
  return 0.;
}


/* Spin temperature minus Tcmb */
double Tspin_Tcmb(double xe, double Tgas, double Tcmb, double nh) {
  double n_kappa10; 
  n_kappa10 = nh *((1.-xe) *kappa10_H(Tgas) 
		   + xe *kappa10_e(Tgas) + xe *kappa10_p(Tgas)); 

  return (Tgas - Tcmb)/(1. + A10/n_kappa10 *Tgas/T10);
}

/* Spin temperature */
double Tspin(double xe, double Tgas, double Tcmb, double nh) {
  return Tspin_Tcmb(xe, Tgas, Tcmb, nh) + Tcmb;   
}

/* Optical depth */
double tau21(double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
  return 3.*T10 /(32.* M_PI*(Tcmb + Tspin_Tcmb(xe, Tgas, Tcmb, nh)))
          *lambda21*lambda21*lambda21 *(1.-xe) *nh *A10 /H_invsec;
}

/* 21 cm mean brightness temperature: tau_21*(Ts - Tcmb)/(1+z) */
double T21cm(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
    return  tau21(xe, Tgas, Tcmb, nh, H_invsec) *Tspin_Tcmb(xe, Tgas, Tcmb, nh)/(1.+z);
}

/* d T21/d\delta_H */
double T21cm_dh(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return ( T21cm(z, xe, Tgas, Tcmb, nh*(1.+eps), H_invsec) 
	    -T21cm(z, xe, Tgas, Tcmb, nh*(1.-eps), H_invsec))/(2.*eps);
}

/* 1/2 d^2 T21/(d\delta_H)^2 */
double T21cm_dh2(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return (   T21cm(z, xe, Tgas, Tcmb, nh*(1.+eps), H_invsec) 
	    + T21cm(z, xe, Tgas, Tcmb, nh*(1.-eps), H_invsec)
            - 2.*T21cm(z, xe, Tgas, Tcmb, nh, H_invsec) )/(2.*eps*eps);
}

/* d T21/d\delta_Tgas */
double T21cm_dtgas(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return ( T21cm(z, xe, Tgas*(1.+eps), Tcmb, nh, H_invsec) 
	    -T21cm(z, xe, Tgas*(1.-eps), Tcmb, nh, H_invsec))/(2.*eps);
}

/* 1/2 d^2 T21/(d\delta_Tgas)^2 */
double T21cm_dtgas2(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return (   T21cm(z, xe, Tgas*(1.+eps), Tcmb, nh, H_invsec) 
	    + T21cm(z, xe, Tgas*(1.-eps), Tcmb, nh, H_invsec)
            - 2.*T21cm(z, xe, Tgas, Tcmb, nh, H_invsec) )/(2.*eps*eps);
}

/* d^2 T21/ (d delta_H  d delta_Tgas ) */
double T21cm_dhdtgas(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec) {
   double eps = 1e-3;
   return (   T21cm(z, xe, Tgas*(1.+eps), Tcmb, nh*(1.+eps), H_invsec) 
	    + T21cm(z, xe, Tgas*(1.-eps), Tcmb, nh*(1.-eps), H_invsec)
            + 2.*T21cm(z, xe, Tgas, Tcmb, nh, H_invsec) 
            - T21cm(z, xe, Tgas, Tcmb, nh*(1.+eps), H_invsec) 
	    - T21cm(z, xe, Tgas, Tcmb, nh*(1.-eps), H_invsec)
            - T21cm(z, xe, Tgas*(1.+eps), Tcmb, nh, H_invsec) 
	    - T21cm(z, xe, Tgas*(1.-eps), Tcmb, nh, H_invsec))/(2.*eps*eps);
}

/* */


/* Print out, from zstart to zend and at every redshift, 
  z, the 21 cm brightness temperature, and the coefficients of its perturbative expansion */
void print_21cm_history(char *filename, REC_COSMOPARAMS *param, double *xe_output, double *Tm_output, 
			double zstart, double zend, int Nz) { 
  double z, xe, Tgas, dlogz;
  FILE *fp = fopen(filename, "w");
  int iz;

  dlogz = log((1.+zstart)/(1.+zend))/(Nz-1.);

  for (iz = 0; iz < Nz; iz++){
     z = (1.+zstart)*exp(-iz *dlogz);
     xe   = rec_interp1d(-log(1.+ZSTART), DLNA, xe_output, param->nz, -log(1.+z));
     Tgas = rec_interp1d(-log(1.+ZSTART), DLNA, Tm_output, param->nz, -log(1.+z));
     fprintf(fp, "%7.2lf %15.15lf %15.15lf  %15.15lf %15.15lf %15.15lf %15.15lf\n", z, 
             T21cm(z, xe, Tgas, param->T0*(1.+z), 1e-6*param->nH0*cube(1.+z), rec_HubbleConstant(param, z)), 
             T21cm_dh(z, xe, Tgas, param->T0*(1.+z), 1e-6*param->nH0*cube(1.+z), rec_HubbleConstant(param, z)), 
             T21cm_dtgas(z, xe, Tgas, param->T0*(1.+z), 1e-6*param->nH0*cube(1.+z), rec_HubbleConstant(param, z)),
	     T21cm_dh2(z, xe, Tgas, param->T0*(1.+z), 1e-6*param->nH0*cube(1.+z), rec_HubbleConstant(param, z)), 
	     T21cm_dtgas2(z, xe, Tgas, param->T0*(1.+z), 1e-6*param->nH0*cube(1.+z), rec_HubbleConstant(param, z)), 
             T21cm_dhdtgas(z, xe, Tgas, param->T0*(1.+z), 1e-6*param->nH0*cube(1.+z), rec_HubbleConstant(param, z)));
  }

  fclose(fp);
}
