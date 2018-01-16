#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


file_out_vis_no_back = '/home/bb510/Code/CAMB/visibilities_rayleigh_only_no_back_approx.txt' 
file_out_vis_back = '/home/bb510/Code/CAMB/visibilities_rayleigh_only_back_approx.txt'
file_out_taudot_back = '/home/bb510/Code/CAMB/taudot_back_approx.txt' 
file_out_taudot_no_back = '/home/bb510/Code/CAMB/taudot_back_approx.txt' 

col1_vis_nob,col2_vis_nob,col3_vis_nob,col4_vis_nob,col5_vis_nob,col6_vis_nob,col7_vis_nob,col8_vis_nob,col9_vis_nob = np.loadtxt(file_out_vis_no_back, usecols = (0,1,2,3,4,5,6,7,8), unpack = True, skiprows = 1)

col1_vis_b,col2_vis_b,col3_vis_b,col4_vis_b,col5_vis_b,col6_vis_b,col7_vis_b,col8_vis_b,col9_vis_b = np.loadtxt(file_out_vis_back, usecols = (0,1,2,3,4,5,6,7,8), unpack = True, skiprows = 1)

col1_taudot_nob,col2_taudot_nob,col3_taudot_nob,col4_taudot_nob,col5_taudot_nob,col6_taudot_nob,col7_taudot_nob,col8_taudot_nob,col9_taudot_nob = np.loadtxt(file_out_taudot_no_back, usecols = (0,1,2,3,4,5,6,7,8), unpack = True, skiprows = 1)

col1_taudot_b,col2_taudot_b,col3_taudot_b,col4_taudot_b,col5_taudot_b,col6_taudot_b,col7_taudot_b,col8_taudot_b,col9_taudot_b = np.loadtxt(file_out_taudot_back, usecols = (0,1,2,3,4,5,6,7,8), unpack = True, skiprows = 1)

int_col5 =  integrate.quad(lambda x: col5_vis_b[x], 0, np.shape(col8_vis_b)[0])
int_col6 =  integrate.quad(lambda x: col6_vis_b[x], 0, np.shape(col8_vis_b)[0])
int_col7 =  integrate.quad(lambda x: col7_vis_b[x], 0, np.shape(col8_vis_b)[0])
int_col8 =  integrate.quad(lambda x: col8_vis_b[x], 0, np.shape(col8_vis_b)[0])
int_col9 =  integrate.quad(lambda x: col9_vis_b[x], 0, np.shape(col8_vis_b)[0])


plt.figure('Visibilities',figsize = (15,10))

plt.plot(col1_vis_nob,col4_vis_nob/np.amax(col4_vis_nob), label = 'Thomson')

plt.plot(col1_vis_nob,col5_vis_nob/np.amax(col4_vis_nob)/int_col5[0]+col6_vis_nob/np.amax(col4_vis_nob)/int_col6[0]+col7_vis_nob/np.amax(col4_vis_nob)/int_col7[0]+ col8_vis_nob/np.amax(col4_vis_nob)/int_col8[0]+col9_vis_nob/np.amax(col4_vis_nob)/int_col9[0], label = 'Rayleigh  (normalized)')
plt.xlabel('Conformal time')
plt.ylabel('Visibility')

plt.plot(col1_vis_nob,col9_vis_nob/np.amax(col4_vis_nob), label = 'Rayleigh freq 857')


plt.xlim(200,500)
plt.legend()
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/visibilities.pdf',format = 'pdf')
plt.figure('Optical depth', figsize = (15,10))

plt.plot(col1_taudot_nob,col4_taudot_nob, label = 'Thomson')
plt.plot(col1_taudot_nob,col7_taudot_nob, label = 'Rayleigh freq 353')
plt.plot(col1_taudot_nob,col5_taudot_nob, label = 'Rayleigh freq 545')
plt.plot(col1_taudot_nob,col9_taudot_nob, label = 'Rayleigh freq 857')
plt.xlabel('Conformal time')
plt.ylabel('Optical depth')


plt.xlim(100,600)
plt.yscale('log')
plt.ylim(1e-8,10)
plt.legend()
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/optical_depth.pdf',format = 'pdf')


plt.figure('Diff')
plt.plot(col1_vis_b,(col4_vis_nob-col4_vis_b)/col4_vis_nob*100, label = 'CMB vis')
plt.plot(col1_vis_b,(col9_vis_nob-col9_vis_b)/col9_vis_nob*100, label = '853 vis')

plt.plot(col1_taudot_b,(col4_taudot_nob-col4_taudot_b)/col4_taudot_nob*100, label = 'CMB taudot')
plt.plot(col1_taudot_b,(col9_taudot_nob-col9_taudot_b)/col9_taudot_nob*100, label = '853 taudot')
plt.legend()
plt.xlim(200,600)

plt.show()



