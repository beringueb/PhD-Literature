#!/usr/bin/python3

import numpy as np

def V2L_SO_sensitivity_calculator_r1(N_LF,N_MF,N_HF,N_UHF,f_sky):
    '''senstivity calculator for Simons Observatory V2 Large telescope'''
    '''vary the number of dicrhoic optics tubes in each band'''
    '''by default you should use 7 total tubes'''
    ## fixed varibles for paramters you shouldn't vary
    # sensitivity per band in uK-rt(s)
    S_LF_27 = 20.8
    S_LF_39 = 14.3
    S_MF_90 = 6.5
    S_MF_150 = 8.1
    S_HF_150 = 6.9
    S_HF_220 = 15.8
    S_UHF_220 = 13.8
    S_UHF_270 = 35.6
    #Observing duration
    n_years = 5.0
    efficiency = 0.4*.5*.85 ## the product of 40% observing efficiency
    ## 50% cuts retention
    ## and keeping only 85% of the map area
    ################
    ## check that the total nunber of tubes is correct
    total_tubes = N_LF+ N_MF+ N_HF+ N_UHF
    if (N_LF+ N_MF+ N_HF+ N_UHF) != 7:
        print("WARNING! You requested:",total_tubes, "optics tune while V2 includes budget for 7")
    ## calculate the sensitivity in each band
    S_27 = S_39 = S_90 = S_150 = S_220 = S_270 = 1e9 ## e.g., make the noise irrelvently high by default
    # include LF receiver contirbutions
    S_27 = 1./np.sqrt( N_LF * S_LF_27**-2. + S_27**-2.)
    S_39 = 1./np.sqrt( N_LF * S_LF_39**-2. + S_39**-2.)
    # include MF receiver contirbutions
    S_90 = 1./np.sqrt( N_MF * S_MF_90**-2. + S_90**-2.)
    S_150 = 1./np.sqrt( N_MF * S_MF_150**-2. + S_150**-2.)
    # include HF receiver contirbutions
    S_150 = 1./np.sqrt( N_HF * S_HF_150**-2. + S_150**-2.)
    S_220 = 1./np.sqrt( N_HF * S_HF_220**-2. + S_220**-2.)
    # include UHF receiver contirbutions
    S_220 = 1./np.sqrt( N_UHF * S_UHF_220**-2. + S_220**-2.)
    S_270 = 1./np.sqrt( N_UHF * S_UHF_270**-2. + S_270**-2.)
    ### calculate the Noise per arcminute, per band
    integration_time = n_years *365.* 24. * 3600. * efficiency
    sky_area = 4.*np.pi * (180/np.pi)**2. * 3600. * f_sky
    N_27 = S_27 * np.sqrt(sky_area / integration_time)
    N_39 = S_39 * np.sqrt(sky_area / integration_time)
    N_90 = S_90 * np.sqrt(sky_area / integration_time)
    N_150 = S_150 * np.sqrt(sky_area / integration_time)
    N_220 = S_220 * np.sqrt(sky_area / integration_time)
    N_270 = S_270 * np.sqrt(sky_area / integration_time)
    ## return these sensitivity and bands
    bands = np.array([27,39,90,150,220,270.]) ## in GHz
    noise_per_arcminute = np.array([N_27,N_39,N_90,N_150,N_220,N_270])
    return(bands,noise_per_arcminute)


## example of running this code
f_sky = 0.4 ## 16,500 square degrees
N_LF = 1 ## number of tubes in the LF band (30/40 GHz)
N_MF = 4 ## number of tubes in the MF band (90/150 GHz)
N_HF = 1 ## number of tubes in the HF band (150/220 GHz)
N_UHF = 1 ## number of tubes in the UHF band (220/270 GHz)
band_centers_GHZ, map_noise_uk_arcmin = V2L_SO_sensitivity_calculator_r1(N_LF,N_MF,N_HF,N_UHF,f_sky)
print("band centers: ", band_centers_GHZ, "[GHz]")
print("map noise per band:",map_noise_uk_arcmin, "[uK/arcmin_CMB]")
print("calcualted for f_sky = ",f_sky, "-or-", f_sky* 4.*np.pi * (180/np.pi)**2., "degrees^2")
