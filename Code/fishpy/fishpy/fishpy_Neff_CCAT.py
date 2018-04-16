#!/usr/bin/python3
"""Main module that drives the analysis"""

import numpy as np
import matplotlib.pyplot as plt
import ini_driver as ini
import CAMB_wrap as cw
import fisher as fi
from scipy.optimize import curve_fit

def noise_model(ell,noise_pix,beam):
    noise = ell*(ell+1.)/(2*np.pi) * (noise_pix*np.pi/60/180) ** 2 * np.exp(ell * (ell + 1.) * (beam * np.pi / 10800.) ** 2 /(8. * np.log(2)) )
    return noise

iteration = 10
noise_pix_init = [5.8,6.3, 16.98, 21.47, 86.58, 249.41, 242487.09]
div = np.linspace(0.3,3,iteration)
noise_pix_list = [list(np.asarray(noise_pix_init)/i) for i in div]
error = np.zeros(iteration)
equiv_noise = np.zeros(iteration)
noise_lowest = np.zeros(iteration)

for k in range(iteration):
    setup, list_experiments_included  = ini.parser('../input/test_in.ini')
    ### GET NOISE POWER SPECTRA FOR EVERY EXPERIMENT ###
    for experiment in list_experiments_included:
        if experiment.name == 'CCAT':
            experiment.noise_pix_T = noise_pix_list[k]
            experiment.get_noise()            
            popt, pcov = curve_fit(noise_model, experiment.NlTT[:,0], experiment.NlTT[:,1]) 
            equiv_noise[k] = popt[0]
            noise_lowest[k] = experiment.noise_pix_T[1]
            experiment.max_min()
            experiment.write_noise_file(setup.noise_root)
            
            #experiment.plot()
            #experiment.plot(experiment.freqs)
        else:
            experiment.get_noise()
            experiment.max_min()
            experiment.write_noise_file(setup.noise_root)
            
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
    for expe in list_combined_experiments:
        np.savetxt(expe.name + 'TT',expe.NlTT)
        
    ### RUN CAMB FOR EACH COMBINED EXPERIMENT ###
    setup.get_list_experiments(list_combined_experiments)
    setup.get_fiducial(fid_file = '/home/bb510/Code/fishpy/input/fiducial.txt')
    cw.init_file(setup)
    fisher_list = []
    for experiment in list_combined_experiments : # loop on every combined experiments
        #cw.parameter_files(setup, experiment)      
        #cw.compile_CAMB(experiment)
        #cw.run_CAMB(experiment)
        fisher = fi.FisherMatrix(setup.param_list,experiment)
        fisher.get_fisher(setup)
        fisher.write_to_txt('/home/bb510/Code/Rayleigh/fisher_matrices/tests')
        fisher_list.append(fisher)    
    error[k] = fi.get_error(fisher_list[0] + fisher_list[1],'N_eff')
    
print(np.shape(noise_lowest))
print(np.shape(equiv_noise))
print(np.shape(error))

tot_write = np.concatenate([[noise_lowest],[equiv_noise],[error]]).transpose()
np.savetxt('error_N_eff.dat',tot_write)   

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.plot(equiv_noise, error)

ax1.set_xlim(1, 20)
ax1.set_xlabel("Equivalent sensitivity when combining every frequency channel [muK-arcmin]")
ax1.set_ylabel("1-sigma error on N_eff")

ax2 = ax1.twiny()
ax2.set_xlabel("Sensitivity at 150 Ghz [muK-arcmin]")
ax2.set_xlim(1, 20)
ax2.set_xticks(list(equiv_noise))
ax2.set_xticklabels(['{:3.1f}'.format(i) for i in noise_lowest])

plt.show()

print(tot_write) 

    




        

    
    
    
            



