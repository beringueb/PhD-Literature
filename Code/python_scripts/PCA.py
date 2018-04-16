import numpy as np
import matplotlib.pyplot as plt

matrix_file = "../fisher_matrices/PCA/fisher_nor_p_PLANCK_PCA.txt"
xe_fid_file = "../visibilities/PCA/xe_pos_0000.00.txt"
xe_pert_root = "../visibilities/PCA/xe_pos_"

N=160
position_list_tmp = np.linspace(200,3000,N)
position_list = np.append([0],position_list_tmp)
z_new = np.linspace(200,3000,5000)
width = (3000-200)/(2*(N+1.))
fish = np.loadtxt(matrix_file)
print(np.shape(fish))
z,xe_fid = np.loadtxt(xe_fid_file,usecols = (0,1), unpack = True)

pert_modes_gauss = np.zeros((160,np.shape(z_new)[0]))

for i in range(160) :
    pert_modes_gauss[i,:] = np.exp(-(z_new-position_list[i+1])**2/(2*width**2))

value,vector = np.linalg.eig(fish)

fig, axs = plt.subplots(2,3, sharex = True)
k=0
for i in range(2) : 
    for j in range(3) : 
        axs[i,j].plot(z_new,vector[:,k].dot(pert_modes_gauss))
        k += 1
        axs[i,j].set_xlim(200,3000)
plt.show()

