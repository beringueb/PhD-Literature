#!/usr/bin/python3

"""Modules that contains classes and methods to read a .ini file (in ../input/) and initialize several classes. """

import os
import configparser
import numpy as np
from . import noise_generator as noise  

QUICKLENS_LOC = "/home/bb510/code/quicklens/" # Location of quicklens, hard coded

class Experiment() : 
    """ Class that contains information about a given CMB experiment : 
        - name : name of the experiment
        - include : whether to include the experiemnt in the total analysis
        - include_P : whether to include polarization data from this experiment
        - include_lens : whether to include lensing potential information from this experiment
        - include_rayleigh : whether to include to Rayleigh scattering signal
        - freqs : list of frequency bands for this experiment (does not include 0 GHz channel as in CAMB)
        - noise_pix_T : list of temperature noise per pixel (in muK-arcmin)
        - noise_pix_P : list of polarization noise per pixel (in muK-arcmin)
        - beam_fwhm : list of beam FWHM (in arcmin)
        - f_sky : fraction of the sky covered by the experiment
        - lmax_T : ell up to which include temperature data
        - lmax_P : ell up to which include polarization data
        - lmin : ell from which to include data
        - NlTT : (lmax+1,n_freqs+2) temperature noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 0 to max lmax, units are ell*(ell+1.)/2pi
        - NlEE : (lmax+1,n_freqs+2) polarization noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 0 to max lmax, units are ell*(ell+1.)/2pi
        - NlPP : (lmax+1,2) lensing noise power spectra. Columns are ell, primary channel (combined freqs), ell from 0 to max lmax, units are (ell*(ell+1.))**2/2pi
    """
    def __init__(self,name,include=True,include_P=True,include_lens=False,include_rayleigh=False,freqs,noise_pix_T,noise_pix_P,beam_fwhm,f_sky=1.,lmax_T=3000,lmax_P=5000,lmin=2):
        """Initialization of an experiment"""
        lmax = max(lmax_T,lmax_P)
        self.name = name
        self.include = include
        self.include_P = include_P
        self.include_lens = include_lens
        self.include_rayleigh = include_rayleigh
        self.freqs = freqs
        self.noise_pix_T = noise_pix_T
        self.noise_pix_P = noise_pix_P
        self.beam_fwhm = beam_fwhm
        self.f_sky = f_sky
        self.lmax_T = lmax_T
        self.lmax_P = lmax_P
        self.lmin = lmin
        self.NlTT = np.zeros((lmax-1,len(freqs)+2)) + np.inf
        self.NlEE = np.zeros((lmax-1,len(freqs)+2)) + np.inf
        self.NlPP = np.zeros((lmax-1,2)) + np.inf
            
    def write_noise_file(self,noise_root):
        """ Writing noise power spectra into a file located at noise_root """
        
        lmax = max(self.lmax_T, self.lmax_P)
        ell = np.linspace(0,lmax,lmax+1).reshape(lmax+1,1)
        if self.include_lens:
            noise_tot = np.concatenate([ell,NlTT[:,1:],NlEE[:,1:],NlPP[:,1]], axis = 1)
        else : 
            noise_tot = np.concatenate([ell,NlTT[:,1:],NlEE[:,1:]], axis = 1)
        header = "{} : ell, Primary, {} GHz + Polarization + Lensing".format(self.name,str(self.freqs))
        file_name = os.path.join(noise_root, "noise_{}_lensing.dat".format(self.name)
        print("Writting noise power spectra to {} ... ".format(file_name), end="")
        np.savetxt(file_name, noise_tot, header = header, fmt = "%d    " + 2*(len(self.freqs)+1)*"%.10e    " + "%.10e", newline = "\n")
        print("Done !")
        
    def read_noise_file(self,noise_root):
        """ Reading noise power spectra from a file located at noise_root. Assumes file already exists."""
        
        file_name = os.path.join(noise_root,"noise_{}_lensing.dat".format(self.name))
        print("Reading noise power spectra from {} ... ".format(file_name), end = "")
        noise_read = np.loadtxt(file_name)
        self.NlTT[:,:] = np.concatenate([noise_read[:,0],noise_read[:,1:len(self.freqs)+2]],axis = 1)
        self.NlEE[:,:] = np.concatenate([noise_read[:,0],noise_read[:,len(self.freqs)+2:-1]],axis = 1)
        self.NlPP[:,:] = np.concatenate([noise_read[:,0],noise_read[:,-1]],axis = 1)
        print("Done ! ")

    def max_min(self):
        """Set noise to infinity for ell outside range"""
        
        #Deleting contributions from ell < lmin
        
        self.NlTT[0:self.l_min+1,1:] = np.inf
        self.NlEE[0:self.l_min+1,1:] = np.inf
        self.NlPP[0:self.l_min+1,1] = np.inf
        
        #Deleting contributions from ell > lmax
        
        self.NlTT[self.lmax_T:, 1:] = np.inf
        self.NlEE[self.lmax_P:, 1:] = np.inf
        if self.include_lens:
            self.NlPP[self.lmax_P:, 1] = np.inf
        else:
            self.NlPP[:,1] = np.inf
        

    def combine_primary(self, list_freqs = 'all'):
        """ Combine frequency channels to get the the primary CMB noise
            if list_freqs is all then all frequncy bands are used
            otherwise uses frequencies in list_freqs
        """
        
        try : 
            assert isinstance(list_freqs,list) or list_freqs == 'all'
        except AssertionError : 
            print("list_freqs must be either a list (not a tuple) or 'all' ... ")
        lmax = max(self.lmax_T,self.lmax_P)    
        temp_T = np.zeros((lmax+1,1))
        temp_E = np.zeros((lmax+1,1))
        if list_freqs == 'all':
            list_freqs = self.freqs
        print("Combining {} GHz channels into primary CMB noise ...".format(str(list_freqs)), end = "    ")
        for freq in list_freqs:
            index = list_freqs.index(freq)
            temp_T += 1./self.NlTT[:,index+2]
            temp_E += 1./self.NlPP[:,index + 2 + len(self.freqs)]
        self.NlTT[:,1] = 1./temp_T
        self.NlEE[:,1+len(self.freqs)] = 1./temp_E
        print("Done !")
        
        
    def get_noise(self):
        """ Method to get the noise power spectra using function in noise_genrator
        """
        lmax = max(self.lmax_T,self.lmax_P)
        if self.name == 'SO':
            try:
                assert hasattr(self,'sensibility')
            except AssertionError:
                print("For So, sensibility parameter should be specified, check .ini file")
            else :
                NlTT, NlEE = noise.SO(self.freqs,self.sensibility,self.fsky,lmax)
                print("Using SO V3 noise genrator")
        elif hasattr(self,'alpha_temp'):
            NlTT, NlEE = noise.atmospheric_noise(self.freqs,self.noise_pix_T,self.beam_FWHM,lmax,self.alpha_temp, self.alpha_pol, self.ell_pivot_temp, self.ell_pivot_pol, self.c_atmo_temp, self.fsky, self.noise_pix_P)
            print("Using atmospheric noise model for {}".format(self.name))
        else:
            NlTT, NlEE = noise.noatmospheric_noise(self.freqs,self.noise_pix_T,self.beam_FWHM,lmax,self.noise_pix_P) 
            print("Using simple noise model (no atmosphere) for {}".format(self.name))
        
        self.NlTT = NlTT
        if self.include_P:
            self.NlEE = NlEE
        if self.include_lens:                
            self.get_lensing(QUICKLENS_LOC,self.lensing_estimators)   
                      

        
        
    def get_lensing(self,lensing_estimators = ["TT","EB"]):
        """ Method to get the lensing noise using quicklens module. Unfortunately for now, quicklens only works with python2 and we'd like to use python3 instead. 
            - lensing_estimators = ["TT","EE","TE","TB","EB"] : list of quadratic estimators to combine (inverse variance weigthing) for the lensing noise estimator, default is all of them but there is some non diagnoal contributions that should be looked for.
        """
        lmax = max(self.lmax_T,self.lmax_P)
        #First save NlTT and NlEE to two files in quicklens/temp
        if not os.path.exists(os.path.join(QUICKLENS_LOC,"temp")):
            os.mkdir(os.path.join(QUICKLENS_LOC,"temp"))
        if os.path.exists(os.path.join(QUICKLENS_LOC,"temp","NlTT.temp")):
            os.remove(os.path.join(QUICKLENS_LOC,"temp","NlTT.temp"))
        if os.path.exists(os.path.join(QUICKLENS_LOC,"temp","NlEE.temp")):
            os.remove(os.path.join(QUICKLENS_LOC,"temp","NlEE.temp"))
        if os.path.exists(os.path.join(QUICKLENS_LOC,"temp","NlPP.temp")):
            os.remove(os.path.join(QUICKLENS_LOC,"temp","NlPP.temp"))
        
        np.savetxt(os.path.join(QUICKLENS_LOC,"temp","NlTT.temp"), self.NlTT)
        np.savetxt(os.path.join(QUICKLENS_LOC,"temp","NlTT.temp"), self.NlTT) # NEED TO BE SURE THE FIRST LINES ARE NOT NAN OR INF
        
        newtext = ""
        with open(os.path.join(QUICKLENS_LOC,"get_lensing.py"),"rw") as f:
            for line in f:
                if "lensing_list =" in line:
                    line1,line2 = line.split('=')
                    line2 = str(lensing_estimators)
                    line = line1 + " = " + line2 + "\n"
                newtext += line
            f.write(newtext)
        
        try:
            import subprocess 
        except ImportError:
            print("Module suprocess needs to be installed")
        subprocess.call("python get_lensing.py"), shell = True, cwd = QUICKLENS_LOC)
        
        self.NlPP = np.loadtxt(os.path.join(QUICKLENS_LOC,"temp","NlPP.temp"))
        
        
        

    def plot(self,list_freqs = None):
        """ Method to plot TT and EE noise curves.
            - list_freqs : if None plot primary channels, else plot frequencies from list_freqs
        """
        
        import matplotlib.pyplot as plt
        ell = self.NlTT[:,0]
        f, axs = plt.subplots(1,2, sharex = True, sharey = True)
        color = ['c','b','k','y','g','r','gold','sienna','coral','navy']
        if list_freqs is None:
            axs[0].loglog(ell,NlTT[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
            axs[1].loglog(ell,NlEE[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
        else :
            for fr in list_freqs:
                i = self.freqs.index(fr)
                axs[0].loglog(ell,NlTT[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
                axs[1].loglog(ell,NlEE[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
        axs[0].set_title("N($\ell$) Temperature")
        axs[1].set_title("N($\ell$) Polarization")
        axs[0].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[1].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[0].set_xlabel("$\ell$")
        axs[1].set_xlabel("$\ell$")
        axs[0].set_ylim(1e-6,10)
        axs[0].set_legend(loc = 'lower right')
        axs[1].set_legend(loc = 'lower right')
        axs[0].set_xlim(10, int(ell[-1]))
        plt.suptitle("Noise power spectra for {}".format(self.name))
        plt.show()
                
 
        
class CombinedExperiment():
    """ Class that contains conbined data for experiments observing the same sky fracion
        - list_experiments : list of the experiments to combine
        - name : name given to the combination of these experiments ('name1' + 'name2' + ...)
        - include_P : whether to include polarization data from these experiments (include_P1 or include_P2 or ...)
        - include_lens : whether to include lensing potential information from these experiments (include_lens1 or include_lens2 or ...)
        - include_rayleigh : whether to include to Rayleigh scattering signal from these experiments (include_rayleigh1 or include_rayleigh2 or ...)
        - freqs : list of frequency bands for these combined experiments, concatenation of frequency channels for experiments that include rayleigh
        - f_sky : fraction of the sky covered by the experiments, min of sky fractions of considered experiments
        - lmax_T : ell up to which include temperature data, max of lmax_T of considered experiments
        - lmax_P : ell up to which include polarization data, max of lmax_P of considered experiments
        - lmin : ell from which to include data, max of lmin of considered experiments
        - NlTT : (lmax+1,n_freqs+2) temperature noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 0 to max lmax, units are ell*(ell+1.)/2pi. Inverse variance sum of NlTT for combined experiments.
        - NlEE : (lmax+1,n_freqs+2) polarization noise power spectra. Columns are ell, primary channel (combined freqs), freqs bands, ell from 0 to max lmax, units are ell*(ell+1.)/2pi.
Inverse variance sum of NlEE for combined experiments that include polarization.
        - NlPP : (lmax+1,2) lensing noise power spectra. Columns are ell, primary channel (combined freqs), ell from 0 to max lmax, units are (ell*(ell+1.))**2/2pi. Inverse variance sum of NlPP for combined experiments that include lensing.
    """
    def __init__(self, list_experiments):
        """ Initialization of a combined experiment
            - list_experiments : list of the experiments to combine.
        """
        self.list_experiments = list_experiments
        
        name = ''
        for expe in list_experiments:
            name += '{} + '.format(expe.name)
        self.name = name[0:-3] # name of the experiment
        
        include_P = False
        include_lens = False
        include_rayleigh = False
        for expe in list_experiments:
            include_P = include_P or expe.include_P
            include_lens = include_lens or expe.include_lens
            include_rayleigh = include_rayleigh or exper.include_rayleigh
        self.include_P = include_P
        self.include_lens = include_lens
        self.include_rayleigh = include_rayleigh
        
        freqs = []
        for expe in list_experiments:
            if expe.include_rayleigh:
                freqs += expe.freqs  # include the freqs of an experiment only if include scattering from that scattering
        self.freqs = freqs.sort()
        
        self.lmax_T = max([expe.lmax_T for expe in list_experiments])
        self.lmax_P = max([expe.lmax_P for expe in list_experiments])
        self.lmin = min([expe.lmin for expe in list_experiments])
        self.fsky = min([expe.fsky for expe in list_experiments])
        
        lmax = max(self.lmax_T,self.lmax_P)
        NlTT = np.zeros((lmax-1,2 + len(self.freqs))) 
        NlEE = np.zeros((lmax-1,2 + len(self.freqs))) 
        noiseTT_tmp = np.zeros((lmax-1,1)) + 1E-15
        noiseEE_tmp = np.zeros((lmax-1,1)) + 1E-15
        NlPP = np.zeros((lmax-1,2)) + np.inf
        i = 0
        for freq in self.freqs:
            noiseTT_tmp = np.zeros((lmax-1,1)) + 1E-15
            noiseEE_tmp = np.zeros((lmax-1,1)) + 1E-15
        
            #look in each experiment if there is a channel at that frequency +/- 10 GHz (arbitrary for now, maybe some reasons for that ?)
            for expe in list_experiments:
                try:
                    ind = expe.freqs.index(freq)
                except ValueError:
                    pass #freq not in this experiment ... 
                else:
                    noiseTT_tmp[0:expe.lmax_T] += 1./(self.fsky/expe.fsky * expe.NlTT[:,2+ind]) #scaling with fsky
                    noiseEE_tmp[0:expe.lmax_P] += 1./(self.fsky/expe.fsky * expe.NlEE[:,2+ind])
            NlTT[:,2+i] += 1./noiseTT_tmp
            NlEE[:,2+i] += 1./noiseEE_tmp
        noiseTT_tmp = np.zeros((lmax-1,1)) + 1E-15
        noiseEE_tmp = np.zeros((lmax-1,1)) + 1E-15
        noisePP_tmp = np.zeros((lmax-1,1)) + 1E-15    
        for expe in list_experiments:
            noiseTT_tmp[0:expe.lmax_T-1] += 1./(self.fsky/expe.fsky * expe.NlTT[:,1])
            noiseEE_tmp[0:expe.lmax_P-1] += 1./(self.fsky/expe.fsky * expe.NlEE[:,1])
            noisePP_tmp[0:max(exp.lmax_T,expe.lmax_P)-1] += 1./(self.fsky/expe.fsky * expe.NlPP[:,1])
        NlTT[:,1] = 1./noiseTT_tmp
        NlEE[:,1] = 1./noiseEE_tmp
        NlPP[:,1] = 1./noisePP_tmp
        self.NlTT = NlTT
        self.NlEE = NlEE
        self.NlPP = NlPP
        
        
    def plot(self,list_freqs = None):
        """ Method to plot TT and EE noise curves.
            - list_freqs : if None plot primary channels, else plot frequencues from list_freqs
        """
        
        import matplotlib.pyplot as plt
        ell = self.NlTT[:,0]
        f, axs = plt.subplots(1,2, sharex = True, sharey = True)
        color = ['c','b','k','y','g','r','gold','sienna','coral','navy']
        if list_freqs is None:
            axs[0].loglog(ell,NlTT[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
            axs[1].loglog(ell,NlEE[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
        else :
            for fr in list_freqs:
                i = self.freqs.index(fr)
                axs[0].loglog(ell,NlTT[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
                axs[1].loglog(ell,NlEE[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
        axs[0].set_title("N($\ell$) Temperature")
        axs[1].set_title("N($\ell$) Polarization")
        axs[0].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[1].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[0].set_xlabel("$\ell$")
        axs[1].set_xlabel("$\ell$")
        axs[0].set_ylim(1e-6,10)
        axs[0].set_legend(loc = 'lower right')
        axs[1].set_legend(loc = 'lower right')
        axs[0].set_xlim(10, int(ell[-1]))
        plt.suptitle("Noise power spectra for {}".format(self.name))
        plt.show()


class Setup():
    """ Class that summarizes the setup used for the Fisher Forecast
        - param_list : lsit of paramters used for Fisher analysis
        - noise_root : directory where noise data are written
        - data_root : directory where power spectra are written
        - fiducial : dict with the fiducial values for the parameters used
        - step : dict of step sizes taken for the derivatives
        - list_experiments : list of (combined) experiments used in the analysis
        - use_BBN : 1 if we want to use BBN consistency for Y_He. If fixed to a given value, assume this value for y_He. Overiden if y_He is in param_list
        - mass_neutrinos : mass of the neutinos (in eV)
        - lmax : ell max used in the analysis, used by CAMB but doesn't override lmax for experiments
        - l_boost : CAMB parameter, keep more terms in the hierarchy evolution
        - boost : CAMB parameter, increase accuracy_boost to decrease time steps, use more k values
    """
    def __init__(self,param_list,noise_root,data_root,list_experiments, use_BBN, mass_neutrinos = 60., l_boost = 1, boost = 1):
        """ Initialization """
        self.param_list = param_list
        self.noise_root = noise_root
        self.data_root = data_root
        self.list_experiments = list_experiments 
        self.use_BBN = use_BBN
        self.mass_neutrinos = mass_neutrinos
        self.l_boost = l_boost
        self.boost = boost
        self.lmax = max([max(expe.lmax_T,expe.lmax_P) for expe in list_experiments])
    
    def get_fiducial(self, fid_file):
        """ Read fiducial values and step sizes for parameters from an external file 
            - fid_file : file from which values are read. Lines must be : name | fiducial value | step size with one header line starting with #
        """  
        fid = {}
        step = {}
        with open(fid_file,'w') as f:
            for line in f:
                 if '#' in line:
                    pass
                 else:  
                    name,fiducial_value,step_size = line.split('|')
                    fid[name.strip()] = float(fiducial_value)
                    step[name.strip()] = float(step_size)
        try:
            assert set(self.param_list).issubset(fid.keys())
        except AssertionError:
            print("Couldn't load fiducial values and step sizes for some parameters. Check {} !".format(fid_file))
        
        self.fiducial = fid
        self.step = step
        

def parser(inifile):
    """ Function that parses arguments from .ini file
        - inifile : loaction of the .ini file
    """
    config = configparser.SafeConfigParser()
    config.read(inifile)        
    sections = config.sections()    
    
    ## READING EXPERIMENTS CONFIGURATIONS ##
    
    for experiment in sections[1:]:
        name = config.get(experiment,'name')
        fsky = config.get(experiment,'fsky')
        include = config.getboolean(experiment,'include')
        include_P = config.getboolean(experiment,'include_polarization')
        include_lens = config.getboolean(experiment,'include_lensing')
        lmin = config.getint(experiment,'lmin')
        lmax_T = config.getint(experiment,'lmax_T')
        lmax_P = config.getint(experiment,'lmax_P')
        freqs = [ int(fr) for fr in config.get(experiment,'freqs').split(',')]
        if congig.has_option(experiment,'noise_pixel_T'):
            noise_pix_T = [float(noise) for noise in config.get(experiment,'noise_pixel_T').split(',')]
        if congig.has_option(experiment,'noise_pixel_P'):
            noise_pix_P = [float(noise) for noise in config.get(experiment,'noise_pixel_P').split(',')]
        if congig.has_option(experiment,'beam_FWHM'):
            beam_FWHM = [float(beam) for noise in config.get(experiment,'beam_FWHM').split(',')]
        if congig.has_option(experiment,'alpha_temp'):
            alpha_temp = [float(alpha) for noise in config.get(experiment,'alpha_temp').split(',')]
         if congig.has_option(experiment,'alpha_pol'):
            alpha_pol = [float(alpha) for noise in config.get(experiment,'alpha_pol').split(',')]
            
            
            
            
self.freqs,self.noise_pix_T,self.beam_FWHM,lmax,self.alpha_temp, self.alpha_pol, self.ell_pivot_temp, self.ell_pivot_pol, self.c_atmo_temp, self.fsky, self.noise_pix_P
        
        
        
            
                  
        
            
    
    
        
        
       

            
            
        
        
    
    








sections = config.sections()

print(sections)
