# PhD
This repository contains some of the literature material used for my PhD at DAMTP.

16/01/2018 : Added HyRec folder wich contains a modified version of HyRec to account for Rayleigh scattering of the CMB off neutral species. Woking version including Hydrogen only. 

17/01/2018 : - Added python scripts used for Rayleigh Scattering. 
             - HyRec is still NOT working correctly, trying to find the correct e volution of x_rayleigh. 

18/01/2018 : - Hyrec is working properly. Outputs the correct neutral species fraction (~0.1% on the spectrum with Recfast)
             - Passed DM_Pann as a parameter in CAMB to constrain it.

26/01/2018 : - Incorporated new noise generator for SO. More accurate and distinction between Temperature and polarization. 
             - Modified noise_gen.py and fisher_simple.py to account for this difference in polarization. Now noise files are 2*nb_fr + 1 columns with : ell,freq1_T, ..., freqn_T,freq1_pol, ..., freqn_pol. 
             - HUGE difference with previous results on forecasts. Need to do some consistency checks -> not reliable for now ... 
            
 
