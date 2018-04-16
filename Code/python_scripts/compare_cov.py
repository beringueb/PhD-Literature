#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

suffix = "_10_p_False_r_True.txt"

cov_new = np.loadtxt('cov_new' + suffix)
cov_old = np.loadtxt('cov' + suffix)

print('New cov {} \n Old cov {}'.format(cov_new, cov_old))
f,(ax1,ax2) = plt.subplots(1,2)
im = ax1.matshow(np.log10(cov_new))
clim=im.properties()['clim']
ax2.matshow(np.log10(cov_old),clim = clim)
f.colorbar(im)

plt.matshow(np.log10(cov_old))
plt.colorbar()
plt.show()


