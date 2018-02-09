/*************************************************************************************************/
/*                 HYREC: Hydrogen and Helium Recombination Code                                 */
/*         Written by Yacine Ali-Haimoud and Chris Hirata (Caltech)                              */
/*                                                                                               */
/*         history.h: functions for recombination history                                        */
/*                                                                                               */
/*************************************************************************************************/


/**** Structure for cosmological parameters. 
      Include information on starting and ending redshit and timestep  
      Added May 2012: the period of Hydrogen recombination is evaluated right away so
      tables of radiation field can be reduced to only what is needed and memory is saved. 
      Modified January 2018 (B. Beringue) modify rec_build_hisory to output ionized fraction of H and to ouput delta_u, gaussian bump in recombination history for a PCA analysis. And modify REC_COSMOPARAMS to include position (mean_x) and width (width_x) of the gaussian bump
      
****/

typedef struct {
   double T0;                          /* CMB temperature today in K*/
   double obh2, omh2, odeh2, okh2;     /* density parameters */
   double Y;                           /* primordial helium abundance */
   double Nnueff;                      /* effective number of neutrinos */

   double fsR, meR;              /* alpha_fs/alpha_fs(today) and me/me(today) (Added April 2012)*/
     
   double pann;                  /* Dark matter annihilation parameter (added January 2015 ) */
   /*Gaussian perturbation to Recombination history parameters (added Januray 2018)*/
   double mean_x;  /*position (in redshift space) of the bump*/
   double width_x;  /*width of the bump*/

   /* Secondary parameters, to avoid recalculating every time */
   double nH0;                  /* density of hydrogen today in m^{-3} */  
   double fHe;                  /* Helium fraction by number */

   long nz;                     /* total number of redshift steps */
   long izH0;                   /* index when H recombination starts to be considered */  
   double zH0;                  /* Redshift at which H recombination starts (zH0 = z[izH0]) */
   long nzrt;                   /* number of redshift steps while radiative transfer is computed */
} REC_COSMOPARAMS;

/**** Structure for switches for additional output ****/

typedef struct {
  int print_spec, print_21cm;
  char file_spec[100];
  char file_21cm[100];
  int Nz_spec, Nz_21cm;
  double zmin_spec, zmax_spec, zmin_21cm, zmax_21cm;
} REC_EXTRAS;

/**** Structure for radiation field outputs. Added January 2015 ****/

typedef struct {
  double **Dfnu;
  double *Dfminus_Ly[3];
} RAD_OUTPUTS;

void allocate_rad_output(RAD_OUTPUTS *rad, long int nzrt);
void free_rad_output(RAD_OUTPUTS *rad);





void rec_get_cosmoparam(FILE *fin, REC_COSMOPARAMS *param, REC_EXTRAS *extras);
double rec_HubbleConstant(REC_COSMOPARAMS *param, double z);
double rec_Tmss(double xe, double Tr, double H, double fHe, 
                double fsR, double meR, double dE_dtdV, double nh);
double rec_dTmdlna(double xe, double Tm, double Tr, double H, double fHe, 
                   double fsR, double meR, double dE_dtdV, double nh);
void rec_get_xe_next1_He(REC_COSMOPARAMS *param, double z_in, double *xHeII, 
                         double *dxHeIIdlna_prev, double *dxHeIIdlna_prev2, int *post_saha); 
double rec_xH1s_postSaha(REC_COSMOPARAMS *param, unsigned iz_out, double z_out, double xHeII_out, 
                         HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
		         double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist, int *post_saha);
void get_rec_next2_HHe(REC_COSMOPARAMS *param, unsigned iz_in, double z_in, double Tm_in, double *xH1s, double *xHeII,
                       HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], double **Dfnu_hist,
                       double *dxHIIdlna_prev,  double *dxHeIIdlna_prev, double *dxHIIdlna_prev2, double *dxHeIIdlna_prev2, int *post_saha);
void rec_get_xe_next1_H(REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out,
                        HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                        double **Dfnu_hist, double *dxedlna_prev, double *dxedlna_prev2, int *post_saha);
void rec_get_xe_next2_HTm(int func_select, REC_COSMOPARAMS *param, double z_in, double xe_in, double Tm_in, double *xe_out, double *Tm_out,
                          HRATEEFF *rate_table, unsigned iz, TWO_PHOTON_PARAMS *twog_params, double **Dfminus_hist, double *Dfminus_Ly_hist[], 
                          double **Dfnu_hist, double *dxedlna_prev, double *dTmdlna_prev, double *dxedlna_prev2, double *dTmdlna_prev2);
void rec_build_history(REC_COSMOPARAMS *param, HRATEEFF *rate_table, TWO_PHOTON_PARAMS *twog_params,
                       double *xe_output, double *xrayleigh_output, double *delta_u_output, double *Tm_output, double **Dfnu_hist, double *Dfminus_Ly_hist[3]);

