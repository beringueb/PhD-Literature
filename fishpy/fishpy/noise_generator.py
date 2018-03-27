#!/usr/bin/python3

"""Modules that contains functions to generate noise power spectra for different CMB experiments. SO is based on V3 calculator. """

import numpy as np
import matplotlib.pyplot as plt

def SO(freqs,sensi,fsky,lmax):
    """ SO V3 noise calculator, includes atmospheric contribution.
        - freqs : frequency bands to keep (usually discard the two lower frequency bands as they are used for foreground removal)
        - sensi : senisibility 0 : threshold, 1 : baseline, 2 : goal
        - fsky : sky fraction observed
        returns : NlTT (lmax+1,nfreqs+2) Temperature noise power spectra (ell, primary (combined later), freqs) 
                  NlEE (lmax+1,nfreqs+2) Polarization noise power spectra (ell, primary (combined later), freqs)
    """
    
    LA_bands = np.array([27,39,93,145,225,280]) #Large aperture telescope bands in GHz
    LA_beam_widths = np.array([7.4,5.1,2.2,1.4,1.0,0.9]) #Large aperture telescope beam FWHM in arcmin  
    
    ### Internal variables ###
    # Large aperture configuration #
    NTubes = np.array([1.,1.,4.,4.,2.,2.]) #Number of dichroic tubes (LF : 27,39, HF : 93,145, UHF : 225,280)
    # sensitivities 
    sensitivities = np.array([[61,30,6.5,8.1,17,42],[48,24,5.4,6.7,15,36],[35,18,3.9,4.2,10,25]])
 #Sensitivies array (3*nfreqs) row 0 : threshold, row 1 : baseline, row 2 : goal
    # 1/f pol: see http://simonsobservatory.wikidot.com/review-of-hwp-large-aperture-2017-10-04 #
    f_knee_pol = np.array([ 700.,700.,700.,700.,700.,700.]) # Array of polarization f_knee 
    alpha_pol =-1.4
    # atmospheric 1/f temp from matthew's model
    C_atmo_temp = np.array([200.,7.7,1800.,12000.,68000.,124000.]) # Array of atmospheric 1/f temp
    alpha_temp = -3.5
    ell_pivot = 1000.
    
    ## calculate the survey area and time ##
    t = 5* 365. * 24. * 3600 #in secs, assuming a 5 years survey
    t = t * 0.2 ## retention after observing efficiency and cuts
    t = t* 0.85 ## a kluge for the noise non-uniformity of the map edges
    A_SR = 4 * np.pi * fsky ## sky areas in Steridians
    A_deg = A_SR * (180/np.pi)**2 ## sky area in square degrees
    A_arcmin = A_deg * 3600.
    
    ## calculate noise temperature spectra ##
    ell = np.linspace(0,lmax,lmax+1)
    weight_temp = sensitivities[sensi,:] / np.sqrt(t) # experimental weight
    map_white_noise = weight_temp * np.sqrt(A_arcmin) # map level white noise level in muK-arcmin
    AN_T = np.zeros((lmax+1,6)) # atmospheric noise contribution    
    for i in range(len(LA_bands)):
        AN_T[:,i] = C_atmo_temp[i] * (ell/ell_pivot)**alpha_temp * A_SR / t / (2. * NTubes[i]) 
    N_ell_T = np.zeros((lmax+1,len(freqs)))
    i = 0
    for fr in freqs:
        if fr in LA_bands:
            index = LA_bands.tolist().index(fr)
            N_ell_T[:,i] = ell * (ell + 1.) / (2.*np.pi) *( (weight_temp[index]) ** 2. * A_SR +0.* AN_T[:,index])
            N_ell_T[:,i] *= np.exp(ell * (ell + 1.) * (LA_beam_widths[index] * np.pi/10800.) ** 2 / (8. * np.log(2)))
            i+=1
        else:
            print("No {:d} GHz channel in SO, please check your ini file ...".format(fr))
            
    ## calculate polarization noise power spectra ##
    AN_P = np.zeros((lmax+1,6))
    for i in range(len(LA_bands)):
        AN_P[:,i] = (ell / f_knee_pol[i]) ** alpha_pol + 1.
    N_ell_P = np.zeros((lmax+1,len(freqs)))
    i = 0
    for fr in freqs:
        if fr in LA_bands:
            index = LA_bands.tolist().index(fr)
            N_ell_P[:,i] = ell * (ell + 1.) / (2.*np.pi) *  (weight_temp[index] * np.sqrt(2)) ** 2. * A_SR *(1. +0.* AN_P[:,index])
            N_ell_P[:,i] *= np.exp(ell * (ell + 1.) * (LA_beam_widths[index] * np.pi/10800.) ** 2 / (8. * np.log(2)) ) 
            i+=1
        else:
            print("No {:d} GHz channel in SO, please check your ini file ...".format(fr))
    
    NlTT = np.concatenate([ell.reshape(lmax+1,1),np.zeros((lmax+1,1)),N_ell_T],axis = 1)
    NlEE = np.concatenate([ell.reshape(lmax+1,1),np.zeros((lmax+1,1)),N_ell_P],axis = 1)
    
    return NlTT,NlEE 
    
def no_atmospheric_noise(freqs,noise_pix_T,beam_FWHM,lmax,noise_pix_P = None):
    """ Routine to compute noise power spectra assuming no atmospheric contribution. In practice, used for PLANCK.
        - freqs : list of frequency bands of the experiment
        - noise_pix_T : list of temperature noise per pixel for each frequency band in muK-arcmin
        - noise_pix_P (optional) : list of polarization noise per pixel for each frequency band in muK-arcmin. If not specified, assumed = sqrt(2)*noise_T
        - beam_FWHM : list of beam FWHM for each frequency bands in arcmin
        - lmax : ell  up to which noise is evaluated
        return : NlTT (lmax+1,nfreqs+2) Temperature noise power spectra (ell, primary (combined later), freqs) 
                 NlEE (lmax+1,nfreqs+2) Polarization noise power spectra (ell, primary (combined later), freqs)
    """
    
    try : 
        assert len(freqs) == len(noise_pix_T)
        assert len(freqs) == len(beam_FWHM)
    except AssertionError: 
        print("Noise per pixel or beam FWHM do not have the correct number of bands, check your ini file!")
    
    if noise_pix_P is None:
        noise_pix_P = np.asarray(noise_pix_T) * np.sqrt(2)
    ell = np.linspace(0,lmax,lmax+1)
    N_ell_T = np.zeros((lmax+1,len(freqs)))
    for i in range(len(freqs)):
        N_ell_T[:,i] = ell * (ell + 1.) / (2.*np.pi) * (noise_pix_T[i] * np.pi / 10800.) ** 2 
        N_ell_T[:,i] *= np.exp(ell * (ell + 1.) * (beam_FWHM[i] * np.pi / 10800.) ** 2 /(8. * np.log(2)) )
    N_ell_P = np.zeros((lmax+1,len(freqs)))
    for i in range(len(freqs)):
        N_ell_P[:,i] = ell * (ell + 1.) / (2.*np.pi) * (noise_pix_P[i] * np.pi / 10800.) ** 2 
        N_ell_P[:,i] *= np.exp(ell * (ell + 1.) * (beam_FWHM[i] * np.pi / 10800.) ** 2 /(8. * np.log(2)) )
        
    NlTT = np.concatenate([ell.reshape(lmax+1,1),np.zeros((lmax+1,1)),N_ell_T],axis = 1)
    NlEE = np.concatenate([ell.reshape(lmax+1,1),np.zeros((lmax+1,1)),N_ell_P],axis = 1)
    
    return NlTT, NlEE
    
def atmospheric_noise(freqs,noise_pix_T,beam_FWHM,lmax,alpha_temp, alpha_pol, ell_pivot_temp, ell_pivot_pol, c_atmo_temp, fsky, noise_pix_P = None):
    """ Routine to compute noise power spectra assuming an atmospehric contribution.
        - freqs : list of frequency bands of the experiment
        - noise_pix_T : list of temperature noise per pixel for each frequency band in muK-arcmin.
        - noise_pix_P (optional) : list of polarization noise per pixel for each frequency band in muK-arcmin. If not specified, assume = sqrt(2)*noise_T
        - beam_FWHM : list of beam FWHM for each frequency band in arcmin
        - lmax : ell up to which noise is evaluated.
        - alpha_temp : temperature atmospheric contribution coefficient
        - alpha_pol : polarization atmospheric contribution coefficient
        - ell_pivot_temp : list of pivot spherical harmonic for temperature atmospheric noise
        - ell_pivot_pol : list of pivot spherical harmonic for polarization atmospheric noise
        - c_atmo_temp : atmospheric noise per frequency band
        - fsky : fraction of the sky observed
        return : NlTT (lmax+1,nfreqs+2) Temperature noise power spectra (ell, primary (combined later), freqs) 
                 NlEE (lmax+1,nfreqs+2) Polarization noise power spectra (ell, primary (combined later), freqs)
    """
    try : 
        assert len(freqs) == len(noise_pix_T)
        assert len(freqs) == len(beam_FWHM)
    except AssertionError: 
        print("Noise per pixel or beam FWHM do not have the correct number of bands, check your ini file!")   
    try : 
        assert len(freqs) == len(c_atmo_temp)
    except AssertionError: 
        print("c_atmo_temp do not have the correct number of bands, check your ini file!")  
    try : 
        assert isinstance(ell_pivot_temp,list)
        assert isinstance(ell_pivot_pol,list)
    except AssertionError:
        print("ell_pivot_temp and ell_pivot_pol must be a list ! It can have only one element if it's the same for all frequencies")
    try:
        assert len(freqs) == len(ell_pivot_temp) or len(ell_pivot_temp) == 1
        assert len(freqs) == len(ell_pivot_pol) or len(ell_pivot_pol) == 1
    except AssertionError:
        print("Pivot scale must be either specified for all frequencies or only a one element list.")
        
    if len(ell_pivot_temp) == 1:
        ell_pivot_temp = [ell_pivot_temp[0] for i in range(len(freqs))]
    if len(ell_pivot_pol) == 1:
        ell_pivot_pol = [ell_pivot_pol[0] for i in range(len(freqs))]
    if noise_pix_P is None:
        noise_pix_P = np.asarray(noise_pix_T) * np.sqrt(2)
    ell = np.linspace(0,lmax,lmax+1)
    print(ell_pivot_temp)
    print(ell_pivot_pol)
    
    t = 5* 365. * 24. * 3600 #in secs, assuming a 5 years survey
    t = t * 0.2 ## retention after observing efficiency and cuts
    t = t* 0.85 ## a kluge for the noise non-uniformity of the map edges
    A_SR = 4 * np.pi * fsky ## sky areas in Steridians
    A_deg = A_SR * (180/np.pi)**2 ## sky area in square degrees
    A_arcmin = A_deg * 3600.
    
    AN_T = np.zeros((lmax+1,len(freqs))) # atmospheric noise contribution    
    for i in range(len(freqs)):
        AN_T[:,i] = c_atmo_temp[i] * (ell/ell_pivot_temp[i])**alpha_temp * A_SR / t 
    N_ell_T = np.zeros((lmax+1,len(freqs)))
    for i in range(len(freqs)):
        N_ell_T[:,i] = ell * (ell + 1.) / (2.*np.pi) *((noise_pix_T[i] * np.pi / 10800.) ** 2  + AN_T[:,i])
        N_ell_T[:,i] *= np.exp(ell * (ell + 1.) * (beam_FWHM[i] * np.pi / 10800.) ** 2 /(8. * np.log(2)) )
    
    AN_P = np.zeros((lmax+1,len(freqs)))
    for i in range(len(freqs)):
        AN_P[:,i] = (ell / ell_pivot_pol[i]) ** alpha_pol + 1.
    N_ell_P = np.zeros((lmax+1,len(freqs)))
    for i in range(len(freqs)):
        N_ell_P[:,i] = ell * (ell + 1.) / (2.*np.pi) * (noise_pix_P[i] * np.pi / 10800.) ** 2 * AN_P[:,i]
        N_ell_P[:,i] *= np.exp(ell * (ell + 1.) * (beam_FWHM[i] * np.pi / 10800.) ** 2 /(8. * np.log(2)) ) 
        
    NlTT = np.concatenate([ell.reshape(lmax+1,1),np.zeros((lmax+1,1)),N_ell_T],axis = 1)
    NlEE = np.concatenate([ell.reshape(lmax+1,1),np.zeros((lmax+1,1)),N_ell_P],axis = 1)
    
    return NlTT, NlEE
        
         
    
def plot_noise(freqs,NlTT,NlEE):
    """ Routine to plot the noise power spectra (temperature and polarization)
        - NlTT (lmax+1,nfreqs+2) Temperature noise power spectra (ell, primary (combined later), freqs) 
        - NlEE (lmax+1,nfreqs+2) Polarization noise power spectra (ell, primary (combined later), freqs)
    """
    
    n_freqs = len(freqs) 
    color = ('c','b','k','y','g','r','gold','sienna','coral','navy')
    ell = NlTT[:,0]
    plt.figure()
    lmax = int(NlTT[-1,0])
    for i in range(n_freqs):
        plt.loglog(ell,NlTT[:,i+2]*2.*np.pi/ell/(ell+1.), color = color[i], label = '{:d} GHz'.format(freqs[i]))
    plt.title("N($\ell$) Temperature")
    plt.ylabel("N($\ell$) [$\mu K^2$-SR] " )
    plt.xlabel("$\ell$")
    plt.ylim(1e-6,1)
    plt.legend(loc = 'lower right')
    plt.xlim(100,lmax)
    
    plt.figure()
    for i in range(n_freqs):
        plt.loglog(ell,NlEE[:,i+2]*2.*np.pi/ell/(ell+1.), color = color[i], label = '{:d} GHz'.format(freqs[i]))
    plt.title("N($\ell$) Polarization")
    plt.ylabel("N($\ell$) [$\mu K^2$-SR] " )
    plt.xlabel("$\ell$")
    plt.ylim(1e-6,1)
    plt.legend(loc = 'lower right')
    plt.xlim(100,lmax)
    
    plt.show()
    
    

        
        
        
    
if __name__ == "__main__":
    #NlTT,NlEE = no_atmospheric_noise([150,226,273,350,405,862],[6.,5.,6.,30.,70.,70000.],[1.4,1.0,0.8,0.6,0.5,0.3],10000)
    #NlTT,NlEE = atmospheric_noise([150,226,273,350,405,862],[6.,5.,6.,30.,70.,70000.],[1.4,1.0,0.8,0.6,0.5,0.3],10000, alpha_pol = -1.4, alpha_temp = -3.5, ell_pivot_temp = [1000.], ell_pivot_pol = [700.], c_atmo_temp = [1800.,12000.,68000.,124000.,6e7,7e8], fsky = 0.24)
    #NlTT, NlEE = no_atmospheric_noise([30,44,70,100,143,217,353],[145,149,137,65,43,66,200],[3,23,14,10,7,5,5,5],5000,[np.inf,np.inf,450,103,81,134,406])
    freqs = [93,145,225,280]
    lmax = 5000
    for sensi in [0,1,2]:
        for fsky in [0.1,0.2,0.4]:
            NlTT,NlEE = SO(freqs,sensi,fsky,lmax)
            np.savetxt("/home/bb510/Code/Rayleigh/noise/SO/NlTT_{:d}_{:3.1f}.dat".format(sensi,fsky),NlTT)
            np.savetxt("/home/bb510/Code/Rayleigh/noise/SO/NlEE_{:d}_{:3.1f}.dat".format(sensi,fsky),NlEE)
            print("Done for {} {}".format(sensi,fsky))
    #plot_noise([150,226,273,350,405,862],NlTT,NlEE)
    #np.savetxt("/home/bb510/Code/Rayleigh/noise/SO/NlTT_PLANCK.dat",NlTT)
    #np.savetxt("/home/bb510/Code/Rayleigh/noise/SO/NlEE_PLANCK.dat",NlEE)
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
