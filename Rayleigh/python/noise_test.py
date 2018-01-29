#!/usr/bin/python3 

import numpy as np 
import matplotlib.pyplot as plt 

freqs = ('Primary', '93GHz', '145GHz', '225GHz', '280GHz')
clr = ('k', 'b', 'c', 'r', 'y')

file_new = "../noise/noise_SO.txt"
#file_old = "../noise/noise_SO_old.txt"

data_new = np.loadtxt(file_new)
#data_old = np.loadtxt(file_old)



plt.figure('Temperature')
for i in range(len(freqs)) :
    plt.plot(data_new[:,0], data_new[:,i+1], label = '{}'.format(freqs[i]), color = clr[i])
    #plt.plot(data_old[:,0], data_old[:,i+1], color = clr[i], ls = 'dashed')
plt.xscale('log')
plt.yscale('log')
plt.xlim(100,4000)
#plt.ylim(1e-6,1000)
plt.legend()
plt.title('Temperature power spectra')

plt.figure('Polarization')
for i in range(len(freqs)) :
    plt.plot(data_new[:,0], data_new[:,6+i], label = '{}'.format(freqs[i]), color = clr[i])
    #plt.plot(data_old[:,0], 2*data_old[:,i+1], color = clr[i], ls = 'dashed')
plt.xscale('log')
plt.yscale('log')
plt.xlim(100,4000)
#plt.ylim(1e-6,1000)
plt.legend()
plt.title('Polarization power spectra')

plt.show()


f, axs = plt.subplots(1,2, sharex = True, sharey = True)
#for i in range(len(freqs)) :
    #axs[0].plot(data_new[:,0],(data_old[:,i+1]-data_new[:,1+i]) / data_old[:,i+1] * 100, label = '{}'.format(freqs[i]), color = clr[i])
    #axs[1].plot(data_new[:,0],(2*data_old[:,i+1]-data_new[:,6+i]) / (2*data_old[:,i+1]) * 100, label = '{}'.format(freqs[i]), color = clr[i])
axs[0].set_xscale('log')
axs[1].set_xscale('log')
axs[0].legend(loc = 'lower right')
axs[1].legend(loc = 'lower right')
axs[0].set_xlim(100,4000)

axs[0].set_title('Relative diff in temperature')
axs[1].set_title('Relative diff in polarization')






