#!/usr/bin/python3


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess


def KroneckerDelta(i,j) :
    if i==j :
        res = 1
    else:
        res = 0
    return res

def modify_param(experiment,param_name, value, which) :
    param_file_init = param_file_root + "params_init.ini"
    with open(param_file_init) as f:
        newText = ''
        for line in f :
            if param_name in line and not('#' in line) :
                line1, line2 = line.split('=')
                if param_name == 'scalar_amp' :
                    line2 = '   {} \n'.format(value*1E-09)
                else :
                    line2 = '   {} \n'.format(value)     
                line = line1 + '=' + line2
            if "output_root" in line and not('#' in line) :
                line = "output_root  = /store/DAMTP/bb510/Rayleigh/Data/camb_out/{}_{}_{} \n".format(experiment,param_name,which)
            newText += line 
    out_file = param_file_root + "params_{}_{}.ini".format(param_name, which) 
    with open(out_file, "w") as f:
        f.write(newText)

def modify_param_init(param_names,mean,l_max) :
    param_file_init = param_file_root + "params.ini"
    #Writing param file for mean value
    with open(param_file_init) as f:
        newText = ''
        for line in f :
            for param_name in param_names :
                if param_name in line and not('#' in line) :
                    line1, line2 = line.split('=')
                    if param_name == 'scalar_amp' : 
                        line2 = '   {} \n'.format(mean[param_name]*1E-09)
                    else : 
                        line2 = '   {} \n'.format(mean[param_name])
                    line = line1 + '=' + line2
            if "l_max_scalar" in line and not('#' in line) :
                    line1, line2 = line.split('=')
                    line2 = '   {} \n'.format(l_max+100)
                    line = line1 + '=' + line2
            newText += line 
    out_file = param_file_root + "params_init.ini"
    with open(out_file, "w") as f:
        f.write(newText)

def run_camb(param_name,which) :
    subprocess.call("./camb {}".format(param_file_root + "params_{}_{}.ini".format(param_name, which)), shell = True, cwd = "/home/bb510/Code/CAMB/", stdout = subprocess.DEVNULL)

def reading_data(experiment,param_name, freqs, freq_tot, which) :
    data_file_root =  "/store/DAMTP/bb510/Rayleigh/Data/camb_out/"
    data = np.zeros((l_max-l_min+1, len(freqs), len(freqs), 3))
    for i in range(len(freqs)) : 
        for j in range(len(freqs)) :
            data_read = pd.read_csv(data_file_root+'{}_{}_{}__{:d}_{:d}'.format(experiment,param_name, which, freq_tot.index(freqs[i])+1, freq_tot.index(freqs[j])+1), header = None, skiprows = 0, sep = '\s+' )
            data[:,i,j,0] = data_read.values[l_min-2:l_max-1,1] #TT
            data[:,i,j,1] = data_read.values[l_min-2:l_max-1,2] #EE
            data[:,i,j,2] = data_read.values[l_min-2:l_max-1,4] #TE
    return data


def def_cov_matrix(data, do_pol, do_rayleigh, freqs) :  
    noise =  pd.read_csv(noise_file, header = None, sep = '\s+' ,skiprows = 1  )
    if do_rayleigh :
        if  do_pol : 
            cov = np.zeros((l_max-l_min+1,2*len(freqs),2*len(freqs)))
            for i in range(len(freqs)) : 
                for j in range(len(freqs)) :
                    cov[:,i,j] = data[:,0,0,0] + (1-KroneckerDelta(i,0))* data[:,i,0,0] + (1-KroneckerDelta(j,0))*data[:,0,j,0] + (1-KroneckerDelta(i,0))*(1-KroneckerDelta(j,0))*data[:,i,j,0] + KroneckerDelta(i,j)*noise.values[l_min-2:l_max-1,i+1]
            for i in range(len(freqs)) :
                for j in range(len(freqs)) : 
                    cov[:,i+len(freqs),j+len(freqs)] = data[:,0,0,1] + (1-KroneckerDelta(i,0))* data[:,i,0,1] + (1-KroneckerDelta(j,0))*data[:,0,j,1] + (1-KroneckerDelta(i,0))*(1-KroneckerDelta(j,0))*data[:,i,j,1] + KroneckerDelta(i,j)*noise.values[l_min-2:l_max-1,i + len(freqs)+1] 
            for i in range(len(freqs)) :
               for j in range(len(freqs)) : 
                   cov[:,i,j+len(freqs)] = cov[:,j+len(freqs),i] =  data[:,0,0,2] + (1-KroneckerDelta(i,0))* data[:,i,0,2] + (1-KroneckerDelta(j,0))*data[:,0,j,2] + (1-KroneckerDelta(i,0))*(1-KroneckerDelta(j,0))*data[:,i,j,2]           
        else :
            cov = np.zeros((l_max-l_min+1,len(freqs),len(freqs))) 
            for i in range(len(freqs)) : 
                for j in range(len(freqs)) :
                    cov[:,i,j] = data[:,0,0,0] + (1-KroneckerDelta(i,0))* data[:,i,0,0] + (1-KroneckerDelta(j,0))*data[:,0,j,0] + KroneckerDelta(i,j)*noise.values[l_min-2:l_max-1,i+1] + (1-KroneckerDelta(i,0))*(1-KroneckerDelta(j,0))*data[:,i,j,0]
                
    else :
        if do_pol :
            cov = np.zeros((l_max-l_min+1,2,2))
            cov[:,0,0] = data[:,0,0,0] + noise.values[l_min-2:l_max-1,1] 
            cov[:,0,1] = data[:,0,0,2]
            cov[:,1,0] = data[:,0,0,2]
            cov[:,1,1] = data[:,0,0,1] + noise.values[l_min-2:l_max-1,1+len(freqs)] 
        else :
            cov = np.zeros((l_max-l_min+1,1,1))
            cov[:,0,0] = data[:,0,0,0] + noise.values[l_min-2:l_max-1,1] 
    return cov


def derivative_cov(experiment,param_name, freqs, delta) : 
    cov = {}
    for which in ['min', 'max', 'mean'] :
        data_cl = reading_data(experiment,param_name, freqs, freq_tot, which)
        cov[which] = def_cov_matrix(data_cl, do_pol, do_rayleigh, freqs)
    if param_name == 'DM_Pann' :
        deriv_cov = (cov['max'] - cov['min']) / (0.5*delta) 
    else :
        deriv_cov = (cov['max'] - cov['min']) / delta 
    return cov['mean'], deriv_cov

def compute_error(cov_mean, deriv_cov_i, deriv_cov_j) :
    fisher_error = np.zeros((l_max-l_min+1))
    for ll in range(l_max-l_min+1) :
        inv_cov = np.linalg.inv(cov_mean[ll,:,:])
        fisher_error[ll] = np.trace(np.dot(inv_cov,np.dot(deriv_cov_i[ll,:,:],np.dot(inv_cov,deriv_cov_j[ll,:,:]))))
    fisher = np.sum((2*l+1.)/2*fisher_error)

    return fisher

def compute_fisher(param_names,mean,delta,do_pol,do_rayleigh,update) : 
    cov_mean = {}
    deriv_cov = {}
    for param_name in param_names :
        for which in ['mean','min','max'] :
            if which == 'mean' :
                value = mean[param_name]
            elif which == 'min' :
                value = max(mean[param_name] - delta[param_name]/2, 0)
            else : 
                value = mean[param_name] + delta[param_name]/2
            if update :
                modify_param(experiment,param_name,value,which)
                run_camb(param_name,which)
            print('Done {} for {}'.format(which,param_name))
        cov_mean[param_name], deriv_cov[param_name] = derivative_cov(experiment,param_name,freqs,delta[param_name])
    fisher_mat = np.zeros((len(param_names),len(param_names)))
    i=0
    for param_i in param_names :
        j = 0 
        for param_j in param_names :
            fisher_mat[i,j] = compute_error(cov_mean[param_i], deriv_cov[param_i], deriv_cov[param_j])
            j+=1
        i+=1
    return fisher_mat


#MAIN PROGRAM

param_names = ['DM_Pann', 'helium_fraction','massless_neutrinos','hubble','ombh2','omch2', 'scalar_amp','scalar_spectral_index','re_optical_depth']
mean = {'DM_Pann':0.,'helium_fraction' : 0.25, 'massless_neutrinos' : 2.99, 'hubble' : 67.27,'ombh2' : 0.2225E-01 ,'omch2' : 0.1198, 'scalar_amp' : 2.21 ,'scalar_spectral_index' : 0.9645,'re_optical_depth' : 0.058}
delta = {'DM_Pann':1.,'helium_fraction' : 0.07, 'massless_neutrinos' : 0.8, 'hubble' : 1.5 ,'ombh2' : 4E-4 ,'omch2' : 4e-3, 'scalar_amp' : 0.1 ,'scalar_spectral_index' : 0.01,'re_optical_depth' : 0.02}
freq_tot = [0,93,145,225,280,350,405,862]
freqs = [0,93,145,225,280]
experiment = 'CCAT-SO'
param_file_root = "/home/bb510/Code/CAMB/"
noise_file = '/home/bb510/Code/Rayleigh/noise/noise_SO_PLANCK.txt'
l_min = 2
l_max = 3000 
l = np.linspace(l_min,l_max, l_max-l_min+1)
plt.figure()

update = False

modify_param_init(param_names,mean,l_max)

for do_pol in [True,False] :
    for do_rayleigh in [True,False] :
        fisher = compute_fisher(param_names,mean,delta,do_pol,do_rayleigh,update)
        if do_pol :
           str1 = 'p'
        else :
           str1 = 'nop'
        if do_rayleigh :
           str2 = 'r'
        else : 
           str2 = 'nor'
        np.savetxt('../fisher_matrices/fisher_{}_{}_SO_PLANCK_full_params_DM.txt'.format(str2,str1), fisher, delimiter = '   ', newline = '\n')
        update = False
print("Done !")

