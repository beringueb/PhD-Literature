#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

def Simons_Observatroy_V3_LA_bands():
    ## returns the band centers in GHz for a CMB spectrum
    ## if your studies require color corrections ask and we can estimates these for you
    return(np.array([27,39,93,145,225,280]))

def Simons_Observatroy_V3_LA_beams():
    ## returns the LAC beams in arcminuts
    beam_LAC_27 = 7.4
    beam_LAC_39 = 5.1
    beam_LAC_93 = 2.2
    beam_LAC_145 = 1.4
    beam_LAC_225 =1.0
    beam_LAC_280 = 0.9
    return(np.array([beam_LAC_27,beam_LAC_39,beam_LAC_93,beam_LAC_145,beam_LAC_225,beam_LAC_280]))

def Simons_Observatroy_V3_LA_noise(sensitivity_mode,f_sky,ell_max,delta_ell):
    """retuns noise curves, including the impact of the beam for the SO large aperature telescopes 
    - sensitivity_mode
        - 0: threshold,
        - 1: baseline,
        - 2: goal
    - f_sky: number from 0-1
    - ell_max: the maximum value of ell used in the computation of N(ell)
    - delta_ell: the step size for computing N_ell """

    ### Internal variables ###
    # Large aperture configuration #
    NTubes_LF = 1
    NTubes_MF = 4.
    NTubes_UHF = 2.
    # sensitivities 
    S_LA_27 = np.array([61,48,35])
    S_LA_39 = np.array([30,24,18])
    S_LA_93 = np.array([6.5,5.4,3.9])
    S_LA_145 = np.array([8.1,6.7,4.2])
    S_LA_225 = np.array([17,15,10])
    S_LA_280 = np.array([42,36,25])
    # 1/f pol: see http://simonsobservatory.wikidot.com/review-of-hwp-large-aperture-2017-10-04 #
    f_knee_pol_LA_27 = 700.
    f_knee_pol_LA_39 = 700.
    f_knee_pol_LA_93 = 700.
    f_knee_pol_LA_145 = 700.
    f_knee_pol_LA_225 = 700.
    f_knee_pol_LA_280 = 700.
    alpha_pol =-1.4
    # atmospheric 1/f temp from matthew's model
    C_27 = 200.
    C_39 = 7.7
    C_93 = 1800.
    C_145 = 12000.
    C_225 = 68000.
    C_280 =124000.
    alpha_temp = -3.5

    ## calculate the survey area and time ##
    t = 5* 365. * 24. * 3600
    t = t * 0.2 ## retention after observing efficiency and cuts
    t = t* 0.85 ## a kluge for the noise non-uniformity of the map edges
    A_SR = 4 * np.pi * f_sky ## sky areas in Steridians
    A_deg = A_SR * (180/np.pi)**2 ## sky area in square degrees
    A_arcmin = A_deg * 3600.
    print("sky area: ", A_deg, "degrees^2")
    ## make the ell array for the output noise curves
    ell = np.arange(2,ell_max,delta_ell)
    ### CALCULATE N(ell) for Temperature ###
    # calculate the experimental weight #
    W_T_27 = S_LA_27[sensitivity_mode] / np.sqrt(t)
    W_T_39 = S_LA_39[sensitivity_mode] / np.sqrt(t)
    W_T_93 = S_LA_93[sensitivity_mode] / np.sqrt(t)
    W_T_145 = S_LA_145[sensitivity_mode] / np.sqrt(t)
    W_T_225 = S_LA_225[sensitivity_mode] / np.sqrt(t)
    W_T_280 = S_LA_280[sensitivity_mode] / np.sqrt(t)
    ## calculate the map noise level (white) for the survey in uK_arcmin for temperature ##
    MN_T_27 = W_T_27 * np.sqrt(A_arcmin)
    MN_T_39 = W_T_39 * np.sqrt(A_arcmin)
    MN_T_93 = W_T_93 * np.sqrt(A_arcmin)
    MN_T_145 = W_T_145 * np.sqrt(A_arcmin)
    MN_T_225 = W_T_225 * np.sqrt(A_arcmin)
    MN_T_280 = W_T_280 * np.sqrt(A_arcmin)
    print("white noise level: ", np.array([MN_T_27,MN_T_39,MN_T_93,MN_T_145,MN_T_225,MN_T_280]),"[uK-arcmin]")
    ## calculate the astmospheric contribution for hews model ##
    ell_pivot = 1000.
    AN_T_27 = C_27 * (ell/ell_pivot)**alpha_temp*A_SR / 18-01t / np.sqrt(NTubes_LF)
    AN_T_39 = C_39 * (ell/ell_pivot)**alpha_temp*A_SR / t / np.sqrt(NTubes_LF)
    AN_T_93 = C_93 * (ell/ell_pivot)**alpha_temp*A_SR / t / np.sqrt(NTubes_MF)
    AN_T_145 = C_145 * (ell/ell_pivot)**alpha_temp*A_SR / t /np.sqrt(NTubes_MF)
    AN_T_225 = C_225 * (ell/ell_pivot)**alpha_temp*A_SR / t /np.sqrt(NTubes_UHF)
    AN_T_280 = C_280 * (ell/ell_pivot)**alpha_temp*A_SR / t /np.sqrt(NTubes_UHF)
    ## calculate
    N_ell_T_27 = (W_T_27* A_SR)**2. + AN_T_27
    N_ell_T_39 = (W_T_39* A_SR)**2. + AN_T_39
    N_ell_T_93 = (W_T_93* A_SR)**2. + AN_T_93
    N_ell_T_145 = (W_T_145* A_SR)**2. + AN_T_145
    N_ell_T_225 = (W_T_225* A_SR)**2. + AN_T_225
    N_ell_T_280 = (W_T_280* A_SR)**2. + AN_T_280

    ## include the impact of the beam ##
    LA_beams = Simons_Observatroy_V3_LA_beams() / np.sqrt(8. * np.log(2)) /60. * np.pi/180.
    ## lac beams as a sigma expressed in radians ##
    N_ell_T_27 *= np.exp( ell*(ell+1)* LA_beams[0]**2 )
    N_ell_T_39 *= np.exp( ell*(ell+1)* LA_beams[1]**2 )
    N_ell_T_93 *= np.exp( ell*(ell+1)* LA_beams[2]**2 )
    N_ell_T_145 *= np.exp( ell*(ell+1)* LA_beams[3]**2 )
    N_ell_T_225 *= np.exp( ell*(ell+1)* LA_beams[4]**2 )
    N_ell_T_280 *= np.exp( ell*(ell+1)* LA_beams[5]**2 )
    ## make an array of nosie curves for T ##
    N_ell_T_LA = np.array([N_ell_T_27,N_ell_T_39,N_ell_T_93,N_ell_T_145,N_ell_T_225,N_ell_T_280])
    ## CALCULATE N(ell) for Polarization ##
    ## calculate the astmospheric contribution for P
    AN_P_27 = (ell / f_knee_pol_LA_27 )**alpha_pol + 1.
    AN_P_39 = (ell / f_knee_pol_LA_39 )**alpha_pol + 1.
    AN_P_93 = (ell / f_knee_pol_LA_93 )**alpha_pol + 1.
    AN_P_145 = (ell / f_knee_pol_LA_145)**alpha_pol + 1.
    AN_P_225 = (ell / f_knee_pol_LA_225)**alpha_pol + 1.
    AN_P_280 = (ell / f_knee_pol_LA_280)**alpha_pol + 1.
    ## calculate and include beam cotribution
    N_ell_P_27 = (W_T_27 * np.sqrt(2) * A_SR)**2 * AN_P_27 
    N_ell_P_27 *= np.exp(ell*(ell+1) * LA_beams[0]**2)
    N_ell_P_39 = (W_T_39 * np.sqrt(2) * A_SR)**2 * AN_P_39 
    N_ell_P_39 *= np.exp(ell*(ell+1) * LA_beams[1]**2)
    N_ell_P_93 = (W_T_93 * np.sqrt(2) * A_SR)**2 * AN_P_93 
    N_ell_P_93 *= np.exp(ell*(ell+1) * LA_beams[2]**2)
    N_ell_P_145 = (W_T_145 * np.sqrt(2) * A_SR)**2 * AN_P_145 
    N_ell_P_145 *= np.exp(ell*(ell+1) * LA_beams[3]**2)
    N_ell_P_225 = (W_T_225 * np.sqrt(2) * A_SR)**2 * AN_P_225 
    N_ell_P_225 *= np.exp(ell*(ell+1) * LA_beams[4]**2)
    N_ell_P_280 = (W_T_280 * np.sqrt(2) * A_SR)**2 * AN_P_280
    N_ell_P_280 *= np.exp(ell*(ell+1) * LA_beams[5]**2)

    ## make an array of nosie curves for T
    N_ell_P_LA = np.array([N_ell_P_27,N_ell_P_39,N_ell_P_93,N_ell_P_145,N_ell_P_225,N_ell_P_280])
    return(ell, N_ell_T_LA,N_ell_P_LA)

print("band centers: ", Simons_Observatroy_V3_LA_bands(), "[GHz]")
print("beam sizes: " , Simons_Observatroy_V3_LA_beams(), "[arcmin]")
## run the code to generate noise curves
ell, N_ell_LA_T,N_ell_LA_Pol = Simons_Observatroy_V3_LA_noise(1,0.4,4001,1)
## plot the temperature noise curves
l = np.reshape(ell,(3999,1))
noiseT = np.transpose(N_ell_LA_T)
noiseP = np.transpose(N_ell_LA_Pol)
data_write = np.zeros((np.shape(l)[0],11))
data_write[:,0] = ell
data_write[:,1] = ell*(ell+1)/(2*np.pi)*noiseT[:,3]
data_write[:,2] = ell*(ell+1)/(2*np.pi)*noiseT[:,2]
data_write[:,3] = ell*(ell+1)/(2*np.pi)*noiseT[:,3]
data_write[:,4] = ell*(ell+1)/(2*np.pi)*noiseT[:,4]
data_write[:,5] = ell*(ell+1)/(2*np.pi)*noiseT[:,5]
data_write[:,6] = ell*(ell+1)/(2*np.pi)*noiseP[:,3]
data_write[:,7] = ell*(ell+1)/(2*np.pi)*noiseP[:,2]
data_write[:,8] = ell*(ell+1)/(2*np.pi)*noiseP[:,3]
data_write[:,9] = ell*(ell+1)/(2*np.pi)*noiseP[:,4]
data_write[:,10] = ell*(ell+1)/(2*np.pi)*noiseP[:,5]


np.savetxt("../noise/noise_T_SO_v3.txt",data_write , delimiter = "   ", newline = "\n")
N_bands = np.size(Simons_Observatroy_V3_LA_bands())
i = 0

color = ('c','b','k','y','g','r')
freq = ( 27,  39,  93, 145, 225, 280)

plt.figure()

while (i < N_bands):
    plt.loglog(ell,N_ell_LA_T[i], color = color[i])
    i+=1

plt.title("N($\ell$) Temperature")
plt.ylabel("N($\ell$) [$\mu K^2$-SR] " )
plt.xlabel("$\ell$")
plt.ylim(1e-6,1)
plt.xlim(100,4000)

plt.figure()

## plot the polarization noise curves
i = 0
while (i < N_bands):
    plt.loglog(ell,N_ell_LA_Pol[i],color = color[i],label = freq[i])
    i+=1
plt.title("N($\ell$) Polairiation")
plt.ylabel("N($\ell$) [$\mu K^2$-SR] " )
#plt.legend(loc = 'lower left')
plt.xlabel("$\ell$")
plt.ylim(1e-6,1)
plt.xlim(100,4000)
plt.show()





