#!/usr/bin/python3
"""Main module that drives the analysis"""

import numpy as np
import matplotlib.pyplot as plt
import ini_driver as ini
import CAMB_wrap as cw
import fisher as fi

if __name__ == "__main__":
    setup, list_experiments_included  = ini.parser('../input/test_in.ini')
    ### GET NOISE POWER SPECTRA FOR EVERY EXPERIMENT ###
    for experiment in list_experiments_included:
        if experiment.update:
            experiment.get_noise()
            experiment.max_min()
            experiment.write_noise_file(setup.noise_root)
            #experiment.plot()
            #experiment.plot(experiment.freqs)
        else:
            experiment.read_noise_file(setup.noise_root)
            #experiment.plot()
            #experiment.plot(experiment.freqs)
            
    ### COMBINE EXPERIMENT ON OVERLAPPING SKY AREAS ###
    fskycurrent = 0.
    list_combined_experiments = [] 
    while len(list_experiments_included) != 0:
        fsky = min([expe.fsky - fskycurrent for expe in list_experiments_included])
        combined_experiment = ini.CombinedExperiment(list_experiments_included,fsky)
        list_combined_experiments.append(combined_experiment)
        for experiment in list_experiments_included:
            if experiment.fsky - fskycurrent == fsky:
                list_experiments_included.remove(experiment)
        fskycurrent = fsky
        
    ### RUN CAMB FOR EACH COMBINED EXPERIMENT ###
    setup.get_list_experiments(list_combined_experiments)
    setup.get_fiducial(fid_file = '/home/bb510/Code/fishpy/input/fiducial.txt')
    cw.init_file(setup)
    for experiment in list_combined_experiments : # loop on every combined experiments
        print(experiment.fsky)       



        

    
    
    
            



