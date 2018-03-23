""" Module CAMC_wrap that provides a wrapper for using CAMB. It will generate one .ini file for the fiducial model and two additional .ini files, one for each direction of a step. 
    It requires the python module camb only for setting background cosmology. The reason camb is not called via pycamb is that pycamb does not include rayleigh scattering.
"""

import numpy as np
import os 

try:
    import subprocess 
except ImportError:
    print("Module suprocess needs to be installed in order to call CAMB")
try:
    import camb
except ImportError:
    print("pycamb (module camb) needs to be installed in order to set background cosmology")
    
CAMB_ROOT = "/home/bb510/Code/CAMB/"  # Location of the CAMB folder.

def init_file(setup):
    """ Function to generate the .ini file for the fiducial cosmology.
        - setup : setup object defined in 
    """ ### TO COMPLETE
    parameter_list = setup.param_list
    pars_fid = camb.CAMBparams()

    if 'H0' in parameter_list:
        pars_fid.set_cosmology(H0=setup.fiducial['H0'],cosmomc_theta=None,ombh2=setup.fiducial['ombh2'],omch2=setup.fiducial['omch2'],mnu = setup.mass_neutrinos/1000.,nnu=setup.fiducial['N_eff'],tau=setup.fiducial['tau'],YHe=None)
    else : 
        pars_fid.set_cosmology(H0=None,cosmomc_theta=setup.fiducial['theta_MC'],ombh2=setup.fiducial['ombh2'],omch2=setup.fiducial['omch2'],mnu = setup.mass_neutrinos/1000.,nnu=setup.fiducial['N_eff'],tau=setup.fiducial['tau'],YHe=None)

    
    param_file = os.path.join(CAMB_ROOT,"params.ini")
    with open(param_file,'r') as f :
        for line in f :
            for param in parameter_list:
                if param == ''
    



