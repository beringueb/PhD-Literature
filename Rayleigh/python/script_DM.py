#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm, chi2
from matplotlib import rc

f_sky = {'PLANCK':0.6,'CCAT':0.24,'SO':0.4,'S4':0.4,'250_GRID' : 0.5,'250_500_GRID' : 0.5,'500_750_GRID' : 0.5,'750_900_GRID' : 0.5,'full_GRID' : 0.5}
mean = (67.27,0.2225E-01 ,0.1198,2.21,0.9645,0.079)
delta = np.asarray((1.5 ,4E-4 ,4e-3,0.12 ,0.01,0.04))
param_names = ['DM_Pann','helium_fraction','massless_neutrinos','hubble','ombh2','omch2', 'scalar_amp','scalar_spectral_index','re_optical_depth']
name = ('$H_0$','$\Omega_bh^2$','$\Omega_ch^2$', '$10^9 A_s$','$n_s$','$/\tau$')
param_list = ['DM_Pann','hubble','ombh2','omch2', 'scalar_amp','scalar_spectral_index','re_optical_depth']

def plot_cov_ellipse(cov, pos, nstd=2, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    a = np.sqrt((cov[0,0]+cov[1,1])/2. + np.sqrt((cov[0,0]-cov[1,1])**2/4. + cov[0,1]**2))
    b = np.sqrt((cov[0,0]+cov[1,1])/2. - np.sqrt((cov[0,0]-cov[1,1])**2/4. + cov[0,1]**2))
    theta = np.degrees(np.arctan2(2*cov[0,1],cov[0,0]-cov[1,1]))/2
    alphas = np.asarray((2.3,6.17,11.8))
    alpha = np.sqrt(alphas[nstd-1]) 
    # Width and height are "full" widths, not radius
    ellip = Ellipse(xy=pos, width=2*alpha*a, height=2*alpha*b, angle=theta, **kwargs)
    ax.add_artist(ellip)
    return ellip

def covariance(experiment_list,param_list) :
    root_data = '/home/bb510/Code/Rayleigh/fisher_matrices/'
    fish = np.zeros((len(param_names),len(param_names)))
    index = [param_names.index(param_list[i]) for i in range(len(param_list))]
    for exp in experiment_list : 
        fish += np.loadtxt(root_data+'fisher_{}.txt'.format(exp))
    cov = np.linalg.inv(fish[np.ix_(index,index)])
    return cov

def get_cov(fisher,param_list) : 
    index = [param_names.index(param_list[i]) for i in range(len(param_list))]
    cov = np.linalg.inv(fisher[np.ix_(index,index)])
    return cov

def combine_experiment(experiment_list, param_list,rayleigh,pola) :
    root_data = '/home/bb510/Code/Rayleigh/fisher_matrices/'
    f_sky_list = []
    for exp in experiment_list :
        f_sky_list.append(f_sky[exp])
    f_sky_list.sort()
    f_sky_list.insert(0,0)
    if rayleigh : 
        str1 = 'r'
    else :
        str1 = 'nor'
    if pola :
        str2 = 'p'
    else :
        str2 = 'nop'
    title = root_data + 'fisher_{}_{}_'.format(str1,str2) + '_'.join(experiment_list) + '_full_params_DM.txt'
    fish = np.zeros((9,9))
    for i in range(len(f_sky_list)-1) : 
        if experiment_list[i] == 'PLANCK' :
            title = title.replace('fisher_r','fisher_nor')
        fish += (f_sky_list[i+1]-f_sky_list[i])*np.loadtxt(title)
        title = title.replace(experiment_list[i]+'_','')
    cov = get_cov(fish,param_list)
    return cov



cov_nor_planck = combine_experiment(['PLANCK'], param_list, rayleigh=False,pola=True)
cov_r_planck = combine_experiment(['PLANCK'], param_list, rayleigh=True,pola=True)

cov_r_SO_planck = combine_experiment(['SO','PLANCK'], param_list,rayleigh=True,pola=True)
cov_nor_SO_planck = combine_experiment(['SO','PLANCK'], param_list,rayleigh=False,pola=True)

cov_r_S4_planck = combine_experiment(['S4','PLANCK'], param_list,rayleigh=True,pola=True)
cov_nor_S4_planck = combine_experiment(['S4','PLANCK'], param_list,rayleigh=False,pola=True)

cov_r_ccat_planck = combine_experiment(['CCAT','PLANCK'], param_list,rayleigh=True,pola=True)
cov_nor_ccat_planck = combine_experiment(['CCAT','PLANCK'], param_list,rayleigh=False,pola=True)

cov_r_ccat_SO_planck = combine_experiment(['CCAT','SO','PLANCK'], param_list,rayleigh=True,pola=True)
cov_r_ccat_S4_planck = combine_experiment(['CCAT','S4','PLANCK'], param_list,rayleigh=True,pola=True)



print(np.shape(cov_nor_planck),np.shape(cov_r_planck))


def display_errors(cov1,cov2 = 'None', param_list= param_list) :
    if cov2 != 'None' :
        for i in range(len(param_list)) :

            print('Error on {} : {:5.3e} then {:5.3e} , improvement {:5.3f} %'.format(param_list[i],np.sqrt(cov1[i,i]),np.sqrt(cov2[i,i]),(np.sqrt(cov1[i,i])-np.sqrt(cov2[i,i]))/np.sqrt(cov1[i,i])*100   ))
    else :
        for i in range(len(param_list)) :
            print('Error on {} : {:5.3e}'.format(param_list[i], np.sqrt(cov1[i,i])))
    return


print('PLANCK \n')
display_errors(cov_nor_planck, cov2 = 'None', param_list = param_list )

print('\n CCAT + PLANCK \n')
display_errors(cov_nor_ccat_planck, cov_r_ccat_planck, param_list )

print('\n SO + PLANCK \n')
display_errors(cov_nor_SO_planck, cov_r_SO_planck, param_list )

print('\n S4 + PLANCK \n')
display_errors(cov_nor_S4_planck, cov_r_S4_planck, param_list )

print('\n CCAT + SO + PLANCK \n')
display_errors(cov_r_SO_planck, cov_r_ccat_SO_planck, param_list )

print('\n CCAT + S4 + PLANCK \n')
display_errors(cov_r_S4_planck, cov_r_ccat_S4_planck, param_list )



        
