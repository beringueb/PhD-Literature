double kappa10_H(double Tgas);
double kappa10_e(double Tgas);
double kappa10_p(double Tgas);
double Tspin_Tcmb(double xe, double Tgas, double Tcmb, double nh);
double Tspin(double xe, double Tgas, double Tcmb, double nh);
double tau21(double xe, double Tgas, double Tcmb, double nh, double H_invsec);
double T21cm(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec);
double T21cm_dh(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec);
double T21cm_dh2(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec);
double T21cm_dtgas(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec);
double T21cm_dtgas2(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec);
double T21cm_dhdtgas(double z, double xe, double Tgas, double Tcmb, double nh, double H_invsec);
void print_21cm_history(char *filename, REC_COSMOPARAMS *param, double *xe_output, double *Tm_output, 
			double zstart, double zend, int Nz);
