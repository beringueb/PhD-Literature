/*********************************************************************************************************/
/*                          HYREC: Hydrogen and Helium Recombination Code                                */
/*                     Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                          */
/*                                                                                                       */
/*         hyrec.c: main module                                                                          */
/*                                                                                                       */
/*         Version: January 2015                                                                         */
/*                                                                                                       */
/*         Revision history:                                                                             */
/*            - written November 2010                                                                    */
/*            - January 2011: changed various switches (notably for post-Saha expansions)                */
/*                             so that they remain valid for arbitrary cosmologies                       */
/*            - November 2011: extended integration down to z = 0 with Peeble's model for z < 20         */
/*                             changed dTm/dlna so it can be called at all times                         */
/*            - May 2012: - included explicit dependence on fine-structure constant and electron mass    */
/*                          - the user can now extract the Lyman-lines distortion (up to Ly-gamma)       */
/*            - January 2015: - changed the way the spectrum parameters are entered.                     */
/*                            - possibility to print 21 cm temperatures                                  */
/*            - January 2016 (B.Beringue) - changed rec_build_history to output fraction of ionized H    */
/*                            - Implemented PCA 
/*********************************************************************************************************/ 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "hyrectools.h"
#include "helium.h"
#include "hydrogen.h"
#include "history.h"
#include "hyrec_params.h"



int main(void) {

   REC_COSMOPARAMS param;
   HRATEEFF rate_table;
   TWO_PHOTON_PARAMS twog_params;
   RAD_OUTPUTS rad;
   double *xe_output, *xrayleigh_output, *delta_u_output, *Tm_output, **Dfnu_hist, *Dfminus_Ly_hist[3];
   long iz;
   REC_EXTRAS extras;
   
   double dz = -1;
   unsigned nz = 8001;
   double z, xe, Tm;

   FILE *fp;
   double *z_spec;
   unsigned b;
   long izmax_spec;
   double prefact, Dfnu;

   /* Allocate and read tables of effective rates */
   read_rates(&rate_table);

   /* Read two-photon rate tables */
   read_twog_params(&twog_params);
  
   /* Get cosmological parameters */
   rec_get_cosmoparam(stdin, &param, &extras);   

   /* allocate memory for output (only nzrt for spectrum arrays) */
   xe_output          = create_1D_array(param.nz);
   xrayleigh_output   = create_1D_array(param.nz);
   Tm_output          = create_1D_array(param.nz);
   delta_u_output     = create_1D_array(param.nz);
   allocate_rad_output(&rad, param.nzrt);
 

   /* Compute the recombination history */
   rec_build_history(&param, &rate_table, &twog_params, xe_output, xrayleigh_output, delta_u_output, Tm_output, rad.Dfnu, rad.Dfminus_Ly);
    
   
   /* Interpolate at the desired output redshifts */
   for(iz=0; iz<nz; iz++) {
       z = ZSTART + dz * iz;    /* print output every dz */
       xe = rec_interp1d(-log(1.+ZSTART), DLNA, xe_output, param.nz, -log(1.+z));
       Tm = rec_interp1d(-log(1.+ZSTART), DLNA, Tm_output, param.nz, -log(1.+z));
       printf("%7.2lf %15.15lf %15.13lf\n", z, xe, Tm/param.T0/(1.+z));
    }

 
   /***** Printing out 21 cm global signal if desired (added December 2014) ******/
   if (extras.print_21cm == 1) {
      print_21cm_history(extras.file_21cm, &param, xe_output, Tm_output,
			 extras.zmax_21cm, extras.zmin_21cm, extras.Nz_21cm);
   }

    /***** Printing out the Lyman-lines spectral distortion *****/

   if (extras.print_spec == 1) {
        z_spec = create_1D_array(extras.Nz_spec);
        for (iz = 0; iz < extras.Nz_spec; iz++) {
           z_spec[iz] = extras.zmin_spec + (extras.zmax_spec - extras.zmin_spec)
	                                  /(extras.Nz_spec - 1.) *iz;             /* Redshifts at which spectrum is printed */
        }   
        izmax_spec = (long) floor(2+log((1.+param.zH0)/(1.+extras.zmin_spec))/DLNA);    /* Max index used for interpolation */
      
        fp = fopen(extras.file_spec, "w");
        
        fprintf(fp, "Number of spectral distortion photons per hydrogen atom per log-frequency interval, (8 pi nu^3)/(c^3 nH) Delta f_nu"); 
        fprintf(fp, "nu/nuLya  z=");
        for (iz = 0; iz < extras.Nz_spec; iz++) fprintf(fp, " %E", z_spec[iz]);
        fprintf(fp, "\n");
      
        for (b = 0; b < NSUBLYA; b++) { /* Sub-Lyman alpha bins */
           fprintf(fp, "%E", twog_params.Eb_tab[b]/E21);
           for (iz = 0; iz < extras.Nz_spec; iz++) {
	     Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, rad.Dfnu[b], izmax_spec, -log(1.+z_spec[iz]));
             prefact = 8.*M_PI*cube(twog_params.Eb_tab[b]/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
             fprintf(fp," %E", prefact*Dfnu);
           }
           fprintf(fp, "\n");
        }
        fprintf(fp, "%E", 1.);    /* Lyman-alpha */
        for (iz = 0; iz < extras.Nz_spec; iz++) {
    	    Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, rad.Dfminus_Ly[0], izmax_spec, -log(1.+z_spec[iz]));
            prefact = 8.*M_PI/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
            fprintf(fp," %E", prefact*Dfnu);
        }
        fprintf(fp, "\n");
        for (; b < NSUBLYB; b++) { /* Bins between Ly-alpha and Ly-beta */
           fprintf(fp, "%E", twog_params.Eb_tab[b]/E21);
           for (iz = 0; iz < extras.Nz_spec; iz++) {
              Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, rad.Dfnu[b], izmax_spec, -log(1.+z_spec[iz]));
              prefact = 8.*M_PI*cube(twog_params.Eb_tab[b]/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
              fprintf(fp," %E", prefact*Dfnu);                
           }
           fprintf(fp, "\n");
        }
        fprintf(fp, "%E", E31/E21);    /* Lyman-beta */
        for (iz = 0; iz < extras.Nz_spec; iz++) {
             Dfnu =  interp_Dfnu(-log(1.+param.zH0), DLNA, rad.Dfminus_Ly[1], izmax_spec, -log(1.+z_spec[iz])); 
             prefact = 8.*M_PI*cube(E31/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
             fprintf(fp," %E", prefact*Dfnu);
        }
        fprintf(fp, "\n");
        for (; b < NVIRT; b++) {   /* Bins between Ly-beta and Ly-gamma */
           fprintf(fp, "%E", twog_params.Eb_tab[b]/E21);
           for (iz = 0; iz < extras.Nz_spec; iz++) {
              Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, rad.Dfnu[b], izmax_spec, -log(1.+z_spec[iz]));
              prefact = 8.*M_PI*cube(twog_params.Eb_tab[b]/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));
              fprintf(fp," %E", prefact*Dfnu);
           }
           fprintf(fp, "\n");
        }
        fprintf(fp, "%E", E41/E21);    /* Lyman-gamma */
        for (iz = 0; iz < extras.Nz_spec; iz++) {
            Dfnu = interp_Dfnu(-log(1.+param.zH0), DLNA, rad.Dfminus_Ly[2], izmax_spec, -log(1.+z_spec[iz]));
            prefact = 8.*M_PI*cube(E41/E21)/(param.nH0 * cube((1.+z_spec[iz]) * 1216e-10));    
            fprintf(fp," %E", prefact*Dfnu);
        }
        fprintf(fp, "\n");

        fclose(fp);
        free(z_spec);
     }
    /*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***/  


    /* Cleanup */
    free_rates(&rate_table);
    free(xe_output);
    free(Tm_output);
    free_rad_output(&rad);

  return(0);
}

/****************************************************************************************************/
