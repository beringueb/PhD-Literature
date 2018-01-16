#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt 

data_diff = {}
freq_def = {}
noise_def = {}
cl_diff = {}
f_sky_def = {'PLANCK':0.6,'CCAT':0.24,'SO':0.4,'S4':0.4}
#data location
data_diff['PLANCK'] = '/mhome/damtp/r/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_diff_TTT__'
data_diff['CCAT'] = '/mhome/damtp/r/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_CCAT_diff_TTT__'
data_diff['SO'] = '/mhome/damtp/r/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_SO_diff_TTT__'
data_diff['S4'] = '/mhome/damtp/r/bb510/Code/CAMB/tests/Rayleigh/Rayleigh_S4_diff_TTT__'

colour = ('c', 'y', 'g', 'b', 'r','m','k',(1,0.55,0)) 


freq_def['PLANCK'] =  (0, 143, 217, 353, 545, 857) #in GHz
freq_def['CCAT'] = (0,150,226,273,350,405) #in GHz
freq_def['SO'] = (0,90,150,220,270)
freq_def['S4'] = (0,95,145,220,270)

#Reading power spectra



def reading_cl(experiment) :
    freq = freq_def[experiment]
    cl_diff = np.zeros((3999, 4, len(freq), len(freq)))
    for i in range(len(freq)) :
        for j in range(len(freq)) :
            l, cl_diff[:,0,i,j],  cl_diff[:,1,i,j], cl_diff[:,2,i,j], cl_diff[:,3,i,j] = np.loadtxt(data_diff[experiment] + '%d_%d' %(i+1,j+1), usecols = (0,1,2,3,4), unpack = True)
    return l, cl_diff


noise_def['CCAT'] = np.loadtxt('noise_CCAT.txt')
noise_def['PLANCK'] = np.loadtxt('noise_PLANCK.txt')
noise_def['SO'] = np.loadtxt('noise_SO.txt')
noise_def['S4'] = np.loadtxt('noise_S4.txt')


def compute_error(experiment):
    noise = noise_def[experiment]
    l,cl_diff = reading_cl(experiment)
    freq = freq_def[experiment]
    l_max = l[-1]
    f_sky = f_sky_def[experiment]
    sigma_min2 = np.zeros((np.shape(l)[0],len(freq)))
    for i in range(len(freq)) : 
        sigma_min2[:,i] = ((2*l+1)*f_sky*(cl_diff[:,0,0,0]*(noise[0:l_max,i+1]+noise[0:l_max,1]) + noise[0:l_max,1]*(noise[0:l_max,i+1]+2*noise[0:l_max,1])))/((cl_diff[:,0,0,0]*(noise[0:l_max,i+1]+noise[0:l_max,+1]) + noise[0:l_max,1]*noise[0:l_max,i+1])**2) 
    return l,sigma_min2


list_experiment = ('PLANCK','CCAT','SO','S4')


f, axs = plt.subplots(2,2,figsize=(15, 10))
k=0
for i in range(2) :
    for j in range(2) :
        experiment = list_experiment[k]
        k+=1
        freq = freq_def[experiment]
        l,sigma_min2 = compute_error(experiment)
        for n in range(len(freq)-1): 
            
            axs[i,j].plot(l, (1/np.sqrt(sigma_min2[:,n+1])), color = colour[n] , linewidth = 2,linestyle = '-', label = ' {} GHz '.format(freq[n+1])) 
        axs[i,j].legend(loc = 'upper left',prop = {'size':10})
        axs[i,j].set_yscale('log')
        axs[i,j].set_xscale('log')
        if k == 1 :
            axs[i,j].set_xlim(2,3000)
        else :
            axs[i,j].set_xlim(30,3000)
        axs[i,j].set_ylim(1e-3,1e5)
        axs[i,j].set_title(experiment + ' errors per frequency channels',fontdict = {'fontsize':14})
        axs[i,j].set_xlabel('l')
        axs[i,j].set_ylabel('Errors [muK^2]')
        axs[i,j].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
plt.suptitle('Errors on the detection of Rayleigh x Primary signal per frequency channels',fontdict = {'fontsize':18})
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/Errors.pdf',format = 'pdf')




f, axs = plt.subplots(2,2,figsize=(15, 10.))
k=0
for i in range(2) :
    for j in range(2) :
        experiment = list_experiment[k]
        k+=1
        freq = freq_def[experiment]
        l,cl_diff = reading_cl(experiment)
        l,sigma_min2 = compute_error(experiment)
        for n in range(len(freq)-1): 
            
            axs[i,j].plot(l, -(cl_diff[:,0,0,n+1]*np.sqrt(sigma_min2[:,n+1])), color = colour[n] , linewidth = 1.5,linestyle = '-', label = ' {} GHz '.format(freq[n+1])) 
        axs[i,j].legend(loc = 'upper right',prop = {'size':10})
        axs[i,j].set_xlim(2,3000)
        axs[i,j].set_ylim(bottom = 0)
        axs[i,j].set_title(experiment + ' S/N per frequency channels',fontdict = {'fontsize':14})
        axs[i,j].set_xlabel('l')
        axs[i,j].set_ylabel('S/N')
        axs[i,j].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
plt.suptitle('Signal-to-noise of the Primary x Rayleigh cross correlation per frequency channels',fontdict = {'fontsize':18})
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/SN.pdf',format = 'pdf')
plt.show()







