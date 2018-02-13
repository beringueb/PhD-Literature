import numpy as np
import matplotlib.pyplot as plt

no_PCA = "../visibilities/PCA/xe_pos_0000.00_1.txt"
PCA_plus = "../visibilities/PCA/xe_pos_1100.00_1.txt"
PCA_minus = "../visibilities/PCA/xe_pos_1100.00_-1.txt"

tau, xe_fid = np.loadtxt(no_PCA, usecols = (0,1), unpack = True)
tau, xe_plus = np.loadtxt(PCA_plus, usecols = (0,1), unpack = True)
tau, xe_minus = np.loadtxt(PCA_minus, usecols = (0,1), unpack = True)

f, axs = plt.subplots(1,2, sharex = True)
axs[0].plot(tau, xe_fid, label = 'fiducial xe', color = 'k')
axs[0].plot(tau, xe_plus, label = 'perturbed xe', color = 'r')
axs[0].plot(tau, xe_minus, label = 'perturbed xe minus', color = 'b')
axs[1].plot(tau, (xe_plus - xe_fid)/xe_fid*100, label = 'plus perturbation', color = 'r')
axs[1].plot(tau, (xe_minus - xe_fid)/xe_fid*100, label = 'minus perturbation', color = 'b')
axs[1].grid()
axs[0].set_xlim(750,1500)
axs[0].legend()
axs[1].legend()
plt.show()

