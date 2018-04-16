#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

num_freq = 6
num_component = 3 #TT,EE,TE
freq_name = ('CMB', '143', '217', '353', '545', '847') 

Diff_file_freq = '/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_diff_TTT__'
CL_file_freq =  '/home/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_nodiff_TTT__'


diff_freq = np.zeros((3999,num_component,num_freq,num_freq))

for i in range(num_freq) :
    for j in range(num_freq) :
        l,diff_freq[:,0,i,j],diff_freq[:,1,i,j],diff_freq[:,2,i,j] = np.loadtxt(Diff_file_freq+'%d_%d' %(i+1,j+1), usecols = (0,1,2,4), unpack = True) 

cl_freq = np.zeros((l.shape[0],num_component,num_freq,num_freq))

for i in range(num_freq) :
    for j in range(num_freq) :
        cl_freq[:,0,i,j],cl_freq[:,1,i,j],cl_freq[:,2,i,j] = np.loadtxt(CL_file_freq+'%d_%d' %(i+1,j+1), usecols = (1,2,4), unpack = True)


f, axs = plt.subplots(3, 3,figsize = (20,15))
axs = axs.ravel()
graph = 0
for i in (2,3,4) :
    for j in (2,3,4) :
        axs[graph].plot(l, -(cl_freq[:,0,0,0]-cl_freq[:,0,i,j])/cl_freq[:,0,0,0], 'r-', label = 'TT') 
        axs[graph].plot(l, -(cl_freq[:,1,0,0]-cl_freq[:,1,i,j])/cl_freq[:,1,0,0], 'b-', label = 'EE') 
        axs[graph].plot(l, -(cl_freq[:,2,0,0]-cl_freq[:,2,i,j])/np.sqrt(cl_freq[:,1,0,0]*cl_freq[:,0,0,0]), 'g-', label = 'TE') 
        axs[graph].set_title(freq_name[i]+'x'+freq_name[j])
        axs[graph].legend(loc = 'upper right',prop = {'size':10}) 
        axs[graph].set_xlim(2,3000)
        axs[graph].grid()
        graph=graph+1 


#f_diff, axs_diff = plt.subplots(5, 1, sharex = True, sharey = True)
#axs_diff = axs_diff.ravel()
#for i in range(num_freq-1) : 
#    axs_diff[i].plot(l,cl_freq_no[:,0,i+1,i+1]-cl_freq[:,0,i+1,i+1]) 
#    axs_diff[i].set_title(freq_name[i+1]+'x'+freq_name[i+1]+'-CMBxCMB') 
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/Effect_rayleigh.pdf',format = 'pdf')
plt.show()
