#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

file_name_diff = '/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_diff_TTT__'  
file_name_nodiff = '/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_diff_TTT__'  

file_name_diff_nor = '/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_diff_FFF__'  
file_name_nodiff_nor = '/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_nodiff_FFF__' 

boolean = False

num_freq = 6
#order == l, TT, EE, BB, TE freq : 0, 143, 217, 353, 545, 857
freq = ('Primary',143, 217, 353, 545, 857)
cl_diff = np.zeros((3999, 4, num_freq,num_freq))

for i in range(num_freq) :
    for j in range(num_freq) :
        l, cl_diff[:,0,i,j],  cl_diff[:,1,i,j], cl_diff[:,2,i,j], cl_diff[:,3,i,j] = np.loadtxt(file_name_diff + '%d_%d' %(i+1,j+1), usecols = (0,1,2,3,4), unpack = True)

cl_nodiff = np.zeros((3999, 4, num_freq,num_freq))

for i in range(num_freq) :
    for j in range(num_freq) :
        l, cl_nodiff[:,0,i,j],  cl_nodiff[:,1,i,j], cl_nodiff[:,2,i,j], cl_nodiff[:,3,i,j] = np.loadtxt(file_name_nodiff + '%d_%d' %(i+1,j+1), usecols = (0,1,2,3,4), unpack = True)

if boolean :
    cl_diff_nor = np.zeros((3999, 4, num_freq,num_freq))

    for i in range(num_freq) :
        for j in range(num_freq) :
            l, cl_diff_nor[:,0,i,j],  cl_diff_nor[:,1,i,j], cl_diff_nor[:,2,i,j], cl_diff_nor[:,3,i,j] = np.loadtxt(file_name_diff_nor + '%d_%d' %(i+1,j+1), usecols = (0,1,2,3,4), unpack = True)

    cl_nodiff_nor = np.zeros((3999, 4, num_freq,num_freq))

    for i in range(num_freq) :
        for j in range(num_freq) :
            l, cl_nodiff_nor[:,0,i,j],  cl_nodiff_nor[:,1,i,j], cl_nodiff_nor[:,2,i,j], cl_nodiff_nor[:,3,i,j] = np.loadtxt(file_name_nodiff_nor + '%d_%d' %(i+1,j+1), usecols = (0,1,2,3,4), unpack = True)


f, axs = plt.subplots(2, 2,figsize = (20,15))
colour = ('y', 'r', 'm', 'g', 'c', 'k','b')
for i in [0,1,2,3,4,5] : 
   axs[0,0].plot(l, np.abs(cl_diff[:,0,i,i]), color = colour[i], linestyle = '-.')
   axs[0,0].plot(l, np.abs(cl_diff[:,0,0,i]), color = colour[i],  label = '%s' %(freq[i]))

axs[0,0].legend(frameon = True, loc = 'upper left',prop={'size': 10})
axs[0,0].set_xscale('log')
axs[0,0].set_yscale('log')
axs[0,0].set_xlim(2,4000)
axs[0,0].set_ylim(1e-6,1e4)
axs[0,0].set_xlabel('l')
axs[0,0].set_ylabel('l(l+1)Cl/2pi')
axs[0,0].set_title('TT power spectrum')


for i in [0,1,2,3,4,5] : 
   axs[1,1].plot(l, cl_diff[:,1,i,i], color = colour[i], linestyle = '-.')
   axs[1,1].plot(l, np.abs(cl_diff[:,1,0,i]), color = colour[i],  label = '%s' %(freq[i]))

axs[1,1].legend(frameon = True, loc = 'upper left',prop={'size': 10})
axs[1,1].set_xscale('log')
axs[1,1].set_yscale('log')
axs[1,1].set_xlim(2,4000)
axs[1,1].set_ylim(1e-6,1e4)
axs[1,1].set_xlabel('l')
axs[1,1].set_ylabel('l(l+1)Cl/2pi')
axs[1,1].set_title('EE power spectrum')


for i in [0,1,2,3,4,5] : 
   axs[0,1].plot(l, np.abs(cl_diff[:,3,i,i]), color = colour[i], linestyle = '-.')
   axs[0,1].plot(l, np.abs(cl_diff[:,3,0,i]), color = colour[i],  label = '%s' %(freq[i]))
 
axs[0,1].legend(frameon = True, loc = 'upper left',prop={'size': 10})
axs[0,1].set_xscale('log')
axs[0,1].set_yscale('log')
axs[0,1].set_xlim(2,4000)
axs[0,1].set_ylim(1e-6,1e4)
axs[0,1].set_xlabel('l')
axs[0,1].set_ylabel('l(l+1)Cl/2pi')
axs[0,1].set_title('TE power spectrum')
   
 
for i in [0,1,2,3,4,5] : 
   axs[1,0].plot(l, np.abs(cl_diff[:,3,i,i]), color = colour[i], linestyle = '-.')
   axs[1,0].plot(l, np.abs(cl_diff[:,3,i,0]), color = colour[i],  label = '%s' %(freq[i]))

axs[1,0].legend(frameon = True, loc = 'upper left',prop={'size': 10})
axs[1,0].set_xscale('log')
axs[1,0].set_yscale('log')
axs[1,0].set_xlim(2,4000)
axs[1,0].set_ylim(1e-6,1e4)
axs[1,0].set_xlabel('l')
axs[1,0].set_ylabel('l(l+1)Cl/2pi')
axs[1,0].set_title('ET power spectrum')




plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/Power_spectra.pdf',format = 'pdf')
plt.show()
