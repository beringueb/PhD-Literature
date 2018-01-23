#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt


file_out_hyrec_vis = '/home/bb510/Code/Rayleigh/visibilities/vis_hyrec.txt'
file_out_recfast_vis = '/home/bb510/Code/Rayleigh/visibilities/vis_recfast.txt' 

file_out_xrayleigh_hyrec = '/home/bb510/Code/Rayleigh/visibilities/out_ray_hyrec.txt'
file_out_xrayleigh_recfast = '/home/bb510/Code/Rayleigh/visibilities/out_ray_recfast.txt'

file_out_hyrec_taudot = '/home/bb510/Code/Rayleigh/visibilities/taudot_hyrec.txt'
file_out_recfast_taudot = '/home/bb510/Code/Rayleigh/visibilities/taudot_recfast.txt' 

recfast_vis1,recfast_vis2,recfast_vis3,recfast_vis4,recfast_vis5,recfast_vis6,recfast_vis7,recfast_vis8,recfast_vis9 = np.loadtxt(file_out_recfast_vis,usecols = (0,1,2,3,4,5,6,7,8), unpack = True)

recfast_taudot1,recfast_taudot2,recfast_taudot3,recfast_taudot4,recfast_taudot5,recfast_taudot6,recfast_taudot7,recfast_taudot8,recfast_taudot9 = np.loadtxt(file_out_recfast_taudot,usecols = (0,1,2,3,4,5,6,7,8), unpack = True)

hyrec_vis1,hyrec_vis2,hyrec_vis3,hyrec_vis4,hyrec_vis5,hyrec_vis6,hyrec_vis7,hyrec_vis8,hyrec_vis9 = np.loadtxt(file_out_hyrec_vis,usecols = (0,1,2,3,4,5,6,7,8), unpack = True)

hyrec_taudot1,hyrec_taudot2,hyrec_taudot3,hyrec_taudot4,hyrec_taudot5,hyrec_taudot6,hyrec_taudot7,hyrec_taudot8,hyrec_taudot9 = np.loadtxt(file_out_hyrec_taudot,usecols = (0,1,2,3,4,5,6,7,8), unpack = True)

hyrec_ray1, hyrec_ray2 = np.loadtxt(file_out_xrayleigh_hyrec, usecols = (0,1), unpack = True, skiprows = 17, comments = 'Reion')
recfast_ray1, recfast_ray2 = np.loadtxt(file_out_xrayleigh_recfast, usecols = (0,1), unpack = True)
print(np.shape(hyrec_ray1), np.shape(recfast_ray1))

fig, axs = plt.subplots(1,1)
axs.plot(hyrec_ray1,hyrec_ray2,markersize = 2, label = 'HyRec',marker = '+',ls = '')
axs.plot(recfast_ray1,recfast_ray2,lw = 2,c = 'r', label = 'Recfast')
axs.set_ylim(-0.1,1.1)
axs.set_xlim(0,0.002)
axs.legend(loc = 'lower right')


fig, axs = plt.subplots(1,2)

axs[0].plot(recfast_vis1,recfast_vis4, label = 'Primary')
axs[0].plot(recfast_vis1,recfast_vis5, label = '143 GHz')
axs[0].plot(recfast_vis1,recfast_vis7, label = '353 GHz')
axs[0].plot(recfast_vis1,recfast_vis9, label = '857 GHz')
axs[0].set_title('Recfast')
axs[0].set_xlim(200,500)
#axs[0].set_xscale('log')
axs[0].legend()

axs[1].plot(hyrec_vis1,hyrec_vis4, label = 'Primary')
axs[1].plot(hyrec_vis1,hyrec_vis5, label = '143 GHz')
axs[1].plot(hyrec_vis1,hyrec_vis7, label = '353 GHz')
axs[1].plot(hyrec_vis1,hyrec_vis9, label = '857 GHz')
axs[1].set_xlim(200,500)
axs[1].set_title('HyRec')
axs[1].legend()
#axs[1].set_xscale('log')
plt.suptitle('Visibility functions')


fig, axs = plt.subplots(1,2)

axs[0].plot(recfast_taudot1,recfast_taudot4, label = 'Primary')
axs[0].plot(recfast_taudot1,recfast_taudot5, label = '143 GHz')
axs[0].plot(recfast_taudot1,recfast_taudot7, label = '353 GHz')
axs[0].plot(recfast_taudot1,recfast_taudot9, label = '857 GHz')
axs[0].set_xlim(100,600)
axs[0].set_yscale('log')
axs[0].set_title('Recfast')
axs[0].set_ylim(1e-7,1e1)
axs[0].legend()
axs[0].set_xscale('log')

axs[1].plot(hyrec_taudot1,hyrec_taudot4, label = 'Primary')
axs[1].plot(hyrec_taudot1,hyrec_taudot5, label = '143 GHz')
axs[1].plot(hyrec_taudot1,hyrec_taudot7, label = '353 GHz')
axs[1].plot(hyrec_taudot1,hyrec_taudot9, label = '857 GHz')
axs[1].set_xlim(100,600)
axs[1].set_title('HyRec')
axs[1].set_yscale('log')
axs[1].set_ylim(1e-7,1e1)
axs[1].legend()
axs[1].set_xscale('log')
plt.suptitle('Optical depths')


fig, axs = plt.subplots(1,2)

axs[0].plot(recfast_taudot1,(recfast_taudot4-hyrec_taudot4)/recfast_taudot4*100, label = 'Primary')
axs[0].plot(recfast_taudot1,(recfast_taudot5-hyrec_taudot5)/recfast_taudot5*100, label = '143 GHz')
axs[0].plot(recfast_taudot1,(recfast_taudot7-hyrec_taudot7)/recfast_taudot7*100, label = '353 GHz')
axs[0].plot(recfast_taudot1,(recfast_taudot9-hyrec_taudot9)/recfast_taudot9*100, label = '857 GHz')
axs[0].set_xlim(0,1000)
axs[0].set_ylim(-5,5)
axs[0].set_xscale('log')
axs[0].set_title('Relative difference in the optical depth')
axs[0].legend()

axs[1].plot(hyrec_vis1,(recfast_vis4-hyrec_vis4)/recfast_vis4*100, label = 'Primary')
axs[1].plot(hyrec_vis1,(recfast_vis5-hyrec_vis5)/recfast_vis5*100, label = '143 GHz')
axs[1].plot(hyrec_vis1,(recfast_vis7-hyrec_vis7)/recfast_vis7*100, label = '353 GHz')
axs[1].plot(hyrec_vis1,(recfast_vis9-hyrec_vis9)/recfast_vis9*100, label = '857 GHz')
axs[1].set_xlim(0,1000)
axs[1].set_ylim(-5,5)
axs[1].set_xscale('log')
axs[1].set_title('Relative difference in the visibility function')
axs[1].legend()

plt.show()








