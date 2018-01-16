#!/usr/bin/python3

import numpy as np 
import matplotlib.pyplot as plt

data_file = "/home/bb510/Code/CAMB/tests/Rayleigh/Fisher_simple/helium_fraction_mean_scalCovCls.dat"

data = np.loadtxt(data_file)
print(np.shape(data)[0])

plt.figure()
plt.scatter(np.arange(np.shape(data)[1]),data[10,:])
plt.ylim(1e-7,1e4)
plt.yscale('log')
plt.show()
