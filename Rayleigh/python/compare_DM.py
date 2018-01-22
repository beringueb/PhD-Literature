#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

out_DM10 = "../visibilities/out_ray_tau_DM.txt"
out_noDM = "../visibilities/out_ray_tau_noDM.txt"

file_cl_noDM = "/home/bb510/Code/CAMB/tests/Rayleigh/test/hyrec_noDM_"
file_cl_DM10 = "/home/bb510/Code/CAMB/tests/Rayleigh/test/hyrec_10DM_"

#order = TT,EE,TE
freq = ('$Primary$', '$143 GHz$', '$217 GHz$', '$353 GHz$', '$545 GHz$', '$857 GHz$')
lmax = 4000

noDM_x, noDM_y = np.loadtxt(out_noDM, usecols = (0,1), unpack = True, skiprows = 18, comments = 'Rei')
DM10_x, DM10_y = np.loadtxt(out_DM10, usecols = (0,1), unpack = True, skiprows = 18, comments = 'Rei')

plt.rc('text', usetex = True)  

fig, axs = plt.subplots(1,2, sharex = True)
axs[0].plot(noDM_x, noDM_y, label = r'No DM annihilation')
axs[0].plot(DM10_x, DM10_y, label = r'DM ammihilation (sigma = 10)')
axs[1].plot(DM10_x, (DM10_y-noDM_y)/noDM_y*100)
axs[0].legend()
#axs[0].set_xlim(0,0.005)

plt.show()

