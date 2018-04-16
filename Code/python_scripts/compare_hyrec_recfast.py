#!/usr/bin/python3 

import matplotlib.pyplot as plt
import numpy as np

root_recfast = "/home/bb510/Code/CAMB/tests/Rayleigh/test/recfast_"
root_hyrec = "/home/bb510/Code/CAMB/tests/Rayleigh/test/hyrec_"

#order = TT,EE,TE
freq = ('$Primary$', '$143 GHz$', '$217 GHz$', '$353 GHz$', '$545 GHz$', '$857 GHz$')

data_recfast = np.zeros((3999,3,6,6))
data_hyrec = np.zeros((3999,3,6,6))

for i in range(6):
    for j in range(6) : 
        l,data_recfast[:,0,i,j],data_recfast[:,1,i,j],data_recfast[:,2,i,j] = np.loadtxt(root_recfast + '_{:d}_{:d}'.format(i+1,j+1),usecols = (0,1,2,3),unpack = True)

for i in range(6):
    for j in range(6) : 
        l,data_hyrec[:,0,i,j],data_hyrec[:,1,i,j],data_hyrec[:,2,i,j] = np.loadtxt(root_hyrec + '_{:d}_{:d}'.format(i+1,j+1),usecols = (0,1,2,3),unpack = True)


plt.rc('text', usetex = True)  
fig,axs = plt.subplots(6,6,figsize = (30,15))
plt.subplots_adjust(wspace=0.15, hspace=0.15,left=0.05,right=0.98,top = 0.90,bottom=0.05)
for i in range(6) : 
    for j in range(6) :
        axs[i,j].plot(l[2:-1],(data_recfast[2:-1:,0,i,j] - data_hyrec[2:-1,0,i,j])/data_recfast[2:-1,0,i,j]*100.)
        
        if i == 5 :
            axs[i,j].set_xlabel('{}'.format(freq[j]))
        if j == 0 :
            axs[i,j].set_ylabel('{}'.format(freq[i]))
            axs[i,j].set_ylim(-0.5,0.5)
        if i == 0 :
            axs[i,j].set_ylim(-0.5,0.5)

        axs[i,j].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
plt.suptitle(r'Difference between Rayleigh signal calculated from Hyrec and Recfast',fontsize =  26)
#plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/diff_hyrec_recfast.pdf', format = 'pdf')
plt.show()

