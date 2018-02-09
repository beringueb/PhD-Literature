import numpy as np
import matplotlib.pyplot as plt

no_PCA = "../tests_PCA/xe_noPCA.txt"
PCA = "../tests_PCA/xe_PCA.txt"


tau, xe_fid = np.loadtxt(no_PCA, usecols = (0,1), unpack = True)
tau, xe_per = np.loadtxt(PCA, usecols = (0,1), unpack = True)

f, axs = plt.subplots(1,2, sharex = True)
axs[0].plot(tau, xe_fid, label = 'fiducial xe', color = 'k')
axs[0].plot(tau, xe_per, label = 'perturbed xe', color = 'r')
axs[1].plot(tau, (xe_fid-xe_per)/xe_fid*100)
axs[1].grid()
axs[0].set_xlim(100,750)
axs[0].legend()
plt.show()

