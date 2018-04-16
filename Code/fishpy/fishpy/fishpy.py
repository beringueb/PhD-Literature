#!/usr/bin/python3
"""Main module that drives the analysis"""

import numpy as np
import matplotlib.pyplot as plt
import ini_driver as ini
import CAMB_wrap as cw
import fisher as fi
import pickle
import os

def loop(ini_file):
    setup, list_experiments_included  = ini.parser(ini_file)
    ### GET NOISE POWER SPECTRA FOR EVERY EXPERIMENT ###
    for experiment in list_experiments_included:
        if experiment.name == 'SO':
            sensi = experiment.sensitivity
            fsky_write = experiment.fsky    
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
        #combined_experiment.plot()
        list_combined_experiments.append(combined_experiment)
        for experiment in list_experiments_included:
            if experiment.fsky - fskycurrent == fsky:
                list_experiments_included.remove(experiment)
        fskycurrent = fsky
        
    ### RUN CAMB FOR EACH COMBINED EXPERIMENT ###
    setup.get_list_experiments(list_combined_experiments)
    setup.get_fiducial(fid_file = '/home/bb510/Code/fishpy/input/fiducial.txt')
    cw.init_file(setup)
    fisher_list = []
    root = '/home/bb510/Code/Rayleigh/fisher_matrices/SO_{:d}_{:3.1f}'.format(sensi,fsky_write)
    if os.path.exists(root):
        pass 
    else:
        os.mkdir(root)
        os.mkdir(os.path.join(root,'pickles'))
        
    for experiment in list_combined_experiments : # loop on every combined experiments
        if update:
            cw.parameter_files(setup, experiment)      
            cw.compile_CAMB(experiment)
            cw.run_CAMB(experiment)
            
        fisher = fi.FisherMatrix(setup.param_list,experiment)
        fisher.get_fisher(setup)
        fisher.write_to_txt(root)
        fisher_list.append(fisher)    

    
if __name__ == "__main__": 
    update = True
    loop('../input/SO.ini')
    update = False
    ini_file_list = os.listdir('../input/SO_ini/')
    for ini_file in ini_file_list:
        filename = os.path.join('../input/SO_ini/',ini_file)
        loop(filename)
           
    
    
    




        

    
    
    
            



