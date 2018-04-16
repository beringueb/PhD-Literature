#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

data_diff = "/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_diff_TTT__"
data_nodiff = "/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_nodiff_TTT__"
freq =  ('Primary', 143, 217, 353, 545, 857) #in GHz


#Reading power spectra
cl_diff = np.zeros((3999, 4, len(freq), len(freq)))


for i in range(len(freq)) :
    for j in range(len(freq)) :
        l, cl_diff[:,0,i,j],  cl_diff[:,1,i,j], cl_diff[:,2,i,j], cl_diff[:,3,i,j] = np.loadtxt(data_diff + '%d_%d' %(i+1,j+1), usecols = (0,1,2,3,4), unpack = True)


colour = ('y', 'r', 'm', 'g', 'c', 'k','b')
f, axs = plt.subplots(1, 1)

for i in [0,1,2,3,4,5] : 
   axs.plot(l, np.abs(cl_diff[:,0,i,i]), color = colour[i], linestyle = '-.')
   axs.plot(l, np.abs(cl_diff[:,0,0,i]), color = colour[i],  label = '%s' %(freq[i]))

axs.legend(frameon = True, loc = 'upper left',prop={'size': 10})
axs.set_xscale('log')
axs.set_yscale('log')
axs.set_xlim(2,4000)
axs.set_ylim(1e-6,1e4)
axs.set_xlabel('l')
axs.set_ylabel('l(l+1)Cl/2pi')
axs.set_title('TT power spectrum')
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/TT.pdf', format = 'pdf')
plt.show()
