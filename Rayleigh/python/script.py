#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm, chi2
from matplotlib import rc

f_sky = {'PLANCK':0.6,'CCAT':0.24,'SO':0.4,'CMBS4':0.4,'250_GRID' : 0.5,'250_500_GRID' : 0.5,'500_750_GRID' : 0.5,'750_900_GRID' : 0.5,'full_GRID' : 0.5}
mean = (0.25,0.99,67.27,0.2225E-01 ,0.1198,2.21,0.9645,0.079)
delta = np.asarray((0.07,1.0,1.5 ,4E-4 ,4e-3,0.12 ,0.01,0.04))
param_names = ['helium_fraction','massless_neutrinos','hubble','ombh2','omch2', 'scalar_amp','scalar_spectral_index','re_optical_depth']
name = ('$Y_p$','$N_{eff}$','$H_0$','$\Omega_bh^2$','$\Omega_ch^2$', '$10^9 A_s$','$n_s$','$/\tau$')
param_list = ['helium_fraction','massless_neutrinos','hubble','ombh2','omch2', 'scalar_amp','scalar_spectral_index','re_optical_depth']


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
    root_data = '/home/bb510/Code/'
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
    title = 'fisher_{}_{}_'.format(str1,str2) + '_'.join(experiment_list) + '_full_params.txt'
    fish = np.zeros((8,8))
    for i in range(len(f_sky_list)-1) : 
        if experiment_list[i] == 'PLANCK' :
            title = title.replace('fisher_r','fisher_nor')
        fish += (f_sky_list[i+1]-f_sky_list[i])*np.loadtxt(title)
        title = title.replace(experiment_list[i]+'_','')
    cov = get_cov(fish,param_list)
    return cov



cov_nor_planck = combine_experiment(['PLANCK'], param_list,rayleigh=False,pola=True)

cov_r_SO_planck = combine_experiment(['SO','PLANCK'], param_list,rayleigh=True,pola=True)
cov_nor_SO_planck = combine_experiment(['SO','PLANCK'], param_list,rayleigh=False,pola=True)

cov_r_S4_planck = combine_experiment(['CMBS4','PLANCK'], param_list,rayleigh=True,pola=True)
cov_nor_S4_planck = combine_experiment(['CMBS4','PLANCK'], param_list,rayleigh=False,pola=True)

cov_r_ccat_planck = combine_experiment(['CCAT','PLANCK'], param_list,rayleigh=True,pola=True)
cov_nor_ccat_planck = combine_experiment(['CCAT','PLANCK'], param_list,rayleigh=False,pola=True)

cov_r_ccat_SO_planck = combine_experiment(['CCAT','SO','PLANCK'], param_list,rayleigh=True,pola=True)
cov_r_ccat_S4_planck = combine_experiment(['CCAT','CMBS4','PLANCK'], param_list,rayleigh=True,pola=True)

cov_r_grid_250 = combine_experiment(['250_GRID','PLANCK'], param_list,rayleigh=True,pola=True)
cov_r_grid_250_500 = combine_experiment(['250_500_GRID','PLANCK'], param_list,rayleigh=True,pola=True)
cov_r_grid_500_750 = combine_experiment(['500_750_GRID','PLANCK'], param_list,rayleigh=True,pola=True)
cov_r_grid_750_900 = combine_experiment(['750_900_GRID','PLANCK'], param_list,rayleigh=True,pola=True)
cov_r_grid_full = combine_experiment(['full_GRID','PLANCK'], param_list,rayleigh=True,pola=True)
cov_nor_grid = combine_experiment(['250_500_GRID','PLANCK'], param_list,rayleigh=False,pola=True)

def display_errors(cov1,cov2 = 'None', param_list= param_list) :
    if cov2 != 'None' :
        for i in range(len(param_list)) :

            print('Error on {} : {:5.3e} then {:5.3e} , improvement {:5.3f} %'.format(param_list[i],np.sqrt(cov1[i,i]),np.sqrt(cov2[i,i]),(np.sqrt(cov1[i,i])-np.sqrt(cov2[i,i]))/np.sqrt(cov1[i,i])*100   ))
    else :
        for i in range(len(param_list)) :
            print('Error on {} : {:5.3e}'.format(param_list[i], np.sqrt(cov1[i,i])))
    return

print('PLANCK \n')
display_errors(cov_nor_planck)
print('\n SO + planck \n')
display_errors(cov_r_SO_planck,cov_r_ccat_SO_planck)
print('\n CMB-S4 + planck \n')
display_errors(cov_r_S4_planck,cov_r_ccat_S4_planck)

        
plt.rc('text', usetex = True)  
fig,axs = plt.subplots(2,2,figsize = (30,15))
plt.subplots_adjust(wspace=0, hspace=0,left=0.05,right=0.95,top = 0.95,bottom=0.05)
for i in range(2) : 
    for j in range(i+1) :
        cov1 = cov_nor_planck
        cov2 = cov_nor_SO_planck
        cov3 = cov_r_SO_planck
        cov4 = cov_nor_ccat_planck
        cov5 = cov_r_ccat_planck

        if i != j : 
            e = plot_cov_ellipse(cov1[np.ix_([j,i],[j,i])],[mean[j],mean[i]],ax=axs[i,j],nstd = 1,lw = 2, ec = 'k', fc = 'k',alpha = 0.2)
            e = plot_cov_ellipse(cov2[np.ix_([j,i],[j,i])],[mean[j],mean[i]],ax=axs[i,j],nstd = 1,lw = 2, ec = 'b', fc = 'None',alpha = 1, ls = 'dashed')
            e = plot_cov_ellipse(cov3[np.ix_([j,i],[j,i])],[mean[j],mean[i]],ax=axs[i,j],nstd = 1,lw = 2, ec = 'b', fc = 'None',alpha = 1)
            e = plot_cov_ellipse(cov4[np.ix_([j,i],[j,i])],[mean[j],mean[i]],ax=axs[i,j],nstd = 1,lw = 2, ec = 'r', fc = 'None',alpha = 1, ls = 'dashed')
            e = plot_cov_ellipse(cov5[np.ix_([j,i],[j,i])],[mean[j],mean[i]],ax=axs[i,j],nstd = 1,lw = 2, ec = 'r', fc = 'None',alpha = 1)

            axs[i,j].scatter(mean[j],mean[i], c = 'r', marker = '+')
            axs[i,j].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
            axs[i,j].set_xlim(mean[j]-delta[j],mean[j]+delta[j])
            axs[i,j].set_ylim(mean[i]-delta[i],mean[i]+delta[i])
            axs[i,j].xaxis.set_ticks(np.linspace(mean[j]-delta[j],mean[j]+delta[j],9))
            axs[i,j].yaxis.set_ticks(np.linspace(mean[i]-delta[i],mean[i]+delta[i],7))
            for label in axs[i,j].xaxis.get_ticklabels()[::2] :
                label.set_visible(False)
            for label in axs[i,j].yaxis.get_ticklabels()[::2] :
                label.set_visible(False)
            axs[i,j].tick_params('both',labelsize = 14)
            fig.delaxes(axs[j,i])
        else :
            xx = np.linspace(mean[i]-2*delta[i],mean[i]+2*delta[i],2000)
            yy1 = np.exp(-(xx-mean[i])**2/(2*cov1[i,i]))
            yy2 = np.exp(-(xx-mean[i])**2/(2*cov2[i,i]))
            yy3 = np.exp(-(xx-mean[i])**2/(2*cov3[i,i]))
            yy4 = np.exp(-(xx-mean[i])**2/(2*cov4[i,i]))
            yy5 = np.exp(-(xx-mean[i])**2/(2*cov5[i,i]))

            axs[i,j].plot(xx,yy1,'k',lw = 1.5, label = r'$TT+TE+EE \quad Planck \quad only$')
            axs[i,j].plot(xx,yy2,'b',lw = 1, ls = 'dashed', label = r'$TT+TE+EE \quad SO + Planck$')
            axs[i,j].plot(xx,yy3,'b',lw = 1, label = r'$TT+TE+EE+Rayleigh \quad SO + Planck$')
            axs[i,j].plot(xx,yy4,'r',lw = 1, ls = 'dashed', label = r'$TT+TE+EE \quad CCAT + Planck$')
            axs[i,j].plot(xx,yy5,'r',lw = 1, label = r'$TT+TE+EE+Rayleigh \quad CCAT + Planck$')

            axs[i,j].xaxis.set_ticks(np.linspace(mean[j]-delta[j],mean[j]+delta[j],9))
            for label in axs[i,j].xaxis.get_ticklabels()[::2] :
                label.set_visible(False)
            axs[i,j].set_xlim(mean[j]-delta[j],mean[j]+delta[j])
            axs[i,j].tick_params('both',labelsize = 14)
            axs[i,j].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
        if i != 0 :
             axs[i,0].set_ylabel(r'{}'.format(name[i]),fontsize = 22)
             
        axs[0,0].set_yticklabels([])
        axs[1,i].set_xlabel(r'{}'.format(name[i]),fontsize = 22) 
        axs[0,0].legend(loc = 'center',prop = {'size':22}, fancybox = True, shadow = True,bbox_to_anchor = (1.5,0.5))
        if i != 1 :
            axs[i,j].set_xticklabels([])
        if j != 0 :
            axs[i,j].set_yticklabels([])

            

plt.suptitle(r'\textbf{Constraints on }$\mathbf{Y_P/N_{eff}}$\textbf{ from TT+TE+EE for fictitous a CMB experiment}',fontsize =  24)
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/combine_CCAT_SO_PLANCK_Y_p.pdf', format = 'pdf')
plt.show()



