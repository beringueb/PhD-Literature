#!/usr/bin/python3

"""Modules that contains classes and methods to read a .ini file (in ../input/) and initialize several classes. """

import os
import configparser
import numpy as np 

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
        
        self.name = name
        self.include = include
        self.include_P = True
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
        self.NlTT = np.zeros((lmax+1,len(freqs)+2))
        self.NlEE = np.zeros((lmax+1,len(freqs)+2))
        self.NlPP = np.zeros((lmax+1,2))
            
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
        print("Writting noise power spectra to {} ...".format(file_name), end="   ")
        np.savetxt(file_name, noise_tot, header = header, fmt = "%d    " + 2*(len(self.freqs)+1)*"%.10e    " + "%.10e", newline = "\n")
        print("Done !")
        
    def read_noise_file(self,noise_root):
        """ Reading noise power spectra from a file located at noise_root. Assumes file already exists."""
        
        file_name = os.path.join(noise_root,"noise_{}_lensing.dat".format(self.name))
        print("Reading noise power spectra from {} ...".format(file_name), end = "    ")
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
        
    def get_lensing(self, quicklens_loc,lensing_estimators = ["TT","EE","TE","TB","EB"]):
        """ Method to get the lensing noise using quicklens module. Unfortunately for now, quicklens only works with python2 and we'd like to use python3 instead. 
            - quicklens_loc : location of the quicklens main repository.
            - lensing_estimators = ["TT","EE","TE","TB","EB"] : list of quadratic estimators to combine (inverse variance weigthing) for the lensing noise estimator, default is all of them but there is some non diagnoal contributions that should be looked for.
        """
        lmax = max(self.lmax_T,self.lmax_P)
        #First save NlTT and NlEE to two files in quicklens/temp
        if not os.path.exists(os.path.join(quicklens_loc,"temp")):
            os.mkdir(os.path.join(quicklens_loc,"temp"))
        if os.path.exists(os.path.join(quicklens_loc,"temp","NlTT.temp")):
            os.remove(os.path.join(quicklens_loc,"temp","NlTT.temp"))
        if os.path.exists(os.path.join(quicklens_loc,"temp","NlEE.temp")):
            os.remove(os.path.join(quicklens_loc,"temp","NlEE.temp"))
        if os.path.exists(os.path.join(quicklens_loc,"temp","NlPP.temp")):
            os.remove(os.path.join(quicklens_loc,"temp","NlPP.temp"))
        
        np.savetxt(os.path.join(quicklens_loc,"temp","NlTT.temp"), self.NlTT)
        np.savetxt(os.path.join(quicklens_loc,"temp","NlTT.temp"), self.NlTT)
        
        newtext = ""
        with open(os.path.join(quicklens_loc,"get_lensing.py"),"rw") as f:
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
        subprocess.call("python get_lensing.py"), shell = True, cwd = quicklens_loc)
        
        self.NlPP = np.loadtxt(os.path.join(quicklens_loc,"temp","NlPP.temp"))

    def plot(self,list_freqs = None):
        """ Method to plot TT and EE noise curves.
            - list_freqs : if None plot primary channels, else plot frequencues fron list-freqs
        """
        
        import matplotlib.pyplot as plt
        ell = self.NlTT[:,0]
        f, axs = plt.subplots(1,2, sharex = True, sharey = True)
        color = ['c','b','k','y','g','r','gold','sienna','coral','navy']
        if list_freqs is None:
            axs[0,0].loglog(ell,NlTT[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
            axs[0,1].loglog(ell,NlEE[:,1]*2*np.pi/ell/(ell+1), label = "Primary channel")
        else :
            for fr in list_freqs:
                i = self.freqs.index(fr)
                axs[0,0].loglog(ell,NlTT[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
                axs[0,1].loglog(ell,NlEE[:,2+i]*2*np.pi/ell/(ell+1), c = color[k], label = "{:d} GHz".format(fr))
        axs[0,0].set_title("N($\ell$) Temperature")
        axs[0,1].set_title("N($\ell$) Polarization")
        axs[0,0].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[0,1].set_ylabel("N($\ell$) [$\mu K^2$-SR] " )
        axs[0,0].set_xlabel("$\ell$")
        axs[0,1].set_xlabel("$\ell$")
        axs[0,0].set_ylim(1e-6,10)
        axs[0,0].set_legend(loc = 'lower right')
        axs[0,1].set_legend(loc = 'lower right')
        axs[0,0].set_xlim(100, int(ell[-1]))
        plt.suptitle("Noise power spectra for {}".format(self.name))
        plt.show()
                
        
        
        
       

            
            
        
        
    
    









config = configparser.SafeConfigParser()
config.read('../input/test_in.ini')





sections = config.sections()

print(sections)
