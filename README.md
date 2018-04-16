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

30/01/2018 : - After a couple of checks, implementation seems to be correct. Similar results when using the new implementation of noise curves and quite insensitive to Recfast/HyRec
             - The new noise curves are significantly higher at low l due to atmosphere. This results in overall poorer constrains but increases the effect of Rayleigh scattering. To be further investigated ...
             - That shows how dependent on the noise modelling our forecasts are. Seems irrelevant to use previsous model for CCAT-p and CMB-S4 as it's very likely to be false. Need to refine our model. 

05/02/2018 : - Added some Literature

07/02/2018 : - Back-up version of HyRec and python before implementing PCA. To come back to this state use git checkout f4b5c50c485ca85ccfe9452002ac276802573fac

09/02/2018 : - Working version of PCA. Modified HyRec (history and hyrec) as well as CAMB (params.ini and HyRec.F90) in order to pass the position as width of the perturbation as parameters. Then modified HyRec to actually perturb the recombination history.
             - Created fisher_simple_PCA.py and PCA.py to actually perform the forecasts. Based on fisher_simple.py, have some issues for now with derivative of the covariance matrix (I took the relative amplitude of the perturbation (1%) whereas I think I should have taken the absolute one)
             - Need to output for each perturbation the perturbed recombination history in order to then be able to combine them linearly and get an idea of the shape of the most constrained modes.

13/02/2018 : - Fixed errors from previous post. Amplitude is now correctly accounted for and each perturbation output correctly. However, it seems that there is an error in the output of the xe. so I recomputed them from scratch, which is faster anyway ... Need to further check that though ! 
             - Test an implementation of +/-1 amplitude to correctly implement the derivative of the covariance matrix. check if i get the results from paper.
 
16/04/2018 : - Fresh start. Re uploaded the repo with CAMB/HyRec/fishpy/quicklens/python_scripts with a different subrepos.            
 
