#!/usr/bin/python3

import numpy as np 

mean_Y_p = 0.2477055
mean_N_eff = 3.046

cov_prism_nor = np.loadtxt('out_nor_p_full_CCAT_S4.txt')
cov_prism_r = np.loadtxt('out_r_p_full_CCAT_S4.txt')

error_Y_p_nor = np.sqrt(cov_prism_nor[0,0])
error_N_eff_nor = np.sqrt(cov_prism_nor[1,1])

error_Y_p_r = np.sqrt(cov_prism_r[0,0])
error_N_eff_r = np.sqrt(cov_prism_r[1,1])

print('Relative errors : \n - Y_p nor = {:6.4f}  N_eff nor = {:6.4f} \n - Y_p r = {:6.4f}  N_eff r = {:6.4f}'.format(error_Y_p_nor/mean_Y_p*100,error_N_eff_nor/mean_N_eff*100,error_Y_p_r/mean_Y_p*100,error_N_eff_r/mean_N_eff*100))

experiment = 'full_PRISM_S4'

fisher = np.linalg.inv(np.loadtxt('out_r_p_{}.txt'.format(experiment)))

fisher_y_p = np.delete(np.delete(fisher,1,0),1,1)
fisher_n_eff = np.delete(np.delete(fisher,0,0),0,1)

print('Experiment : {} \n - error on Y_p : {:6.4e} \n - error on N_eff : {:6.4e}'.format(experiment,np.sqrt(np.linalg.inv(fisher_y_p)[0,0]),np.sqrt(np.linalg.inv(fisher_n_eff)[0,0])))
