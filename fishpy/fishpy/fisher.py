#!/usr/bin/python3 
""" Module that contains the definition of the fisher_matrix class and methods in it"""

import numpy as np
import os
import time


class FisherMatrix():  
    """ Class that defines a fisher_matrix and methods useful to generate it, add two of them, reshuffle them, etc ... """
    
    def __init__(self,param_list,experiment):
        """ Initilaization of the fisher matrices with useful parameters
            - param_list :
            - experiment: 
        """ #TO COMPLETE
        self.experiment = experiment # experiment corresponding to the fisher matrix
        self.param_list = setup.param_list # list of the parameters to include in the fisher matrix
        #Give a name to the fisher matrix : TT + TE + EE + PP + Rayleigh
        name = ': TT ({:d}) + '.format(experiment.lmax_T)
        if experiment.include_P:
            name += 'TE + EE ({:d}) + '.format(experiment.lmax_P)
        if experiment.include_lensing:
            name += 'PP + '
        if experiment.include_rayleigh:
            name += 'Rayleigh + {}  '.format(str(experiment.freqs))
        self.name = experiment.name + name[0:-2]
        self.fisher = np.zeros(( len(self.param_list), len(self.param_list) )) # fisher matrix per se
        
    def write_to_txt(self, root):
        """ Method to write the fisher matrix to a txt file in root
            - root : directory to which the txt file is written
        """
        str0 = self.name[0:self.name.index(':')]
        str1 = 'p' if self.experiment.include_P else 'nop'
        str2 = 'l' if self.experiment.include_lensing else 'nol'
        str3 = 'r' if self.experiment.include_rayleigh else 'nor'
        file_name = os.path.join(root,"fisher_{}_{}_{}_{}.fish".format(str0.strip(),str1,str2,str3)) 
        header = '{} \n {}'.format(self.name,str(self.param_list))
        print("Saving fisher matrix to {} ... ".format(file_name),end = '')
        np.savetxt(file_name,self.fisher,header=header)
        print("Done !")
        
    def get_fisher(self,setup):
        """ Method to compute the fisher matrix :
            - setup:
            -experiment:
        """
        # TO COMPLETE
        
        data_root = os.path.join(setup.data_root,self.experiment.name)
        lmax = max(self.experiment.lmax_T,self.experiment.lmax_P)
        l = np.linspace(2,lmax,lmax-1)
        #Get Cl for the fiducial cosmology
        cls_fid = read_cl(data_root, lmax, len(self.experiment.freqs), 'fiducial')
        cov_fid = compute_covariance(self.experiment,cls_fid)
        deriv = {}
        for parameter in self.param_list:
            deriv[parameter] = derivative(setup,self.experiment,data_root,parameter)
        i = 0
        time_start = time.time()
        print("Computing fisher matrix of {} for {:d} parameters ... ".format(self.name,len(self.param_list)), end = '')
        for param_i in self.param_list:
            j = 0
            for param_j in self.param_list:
                self.fisher[i,j] = self.experiment.fsky * np.sum((2*l+1.)/2.*fish_trace(lmax,deriv[param_i],deriv[param_j],cov_fid)) 
                j += 1
            time_tmp = time.time()
            ETA = (len(self.param_list) - i)*(time_tmp - time_start) / i
            print("{:3.1f}% done, ETA : {:2.0f} min {:2.0f} secs".format(i/n * 100, ETA // 60, ETA % 60), end = "\r" )
            i += 1
        print("Done in {:2.0f} min {:2.0f} secs".format(ETA // 60, ETA % 60))
        
    def reshuffle(self,new_param_list):
        """ Method to reshuffle the fisher matrix witha a new parameter list 
            - new_param_list : New order of the parameter list, has to be of the same length than old one, separate function to fix parameter or marginalize.
        """
        old_param_list = self.param_list
        try : 
            assert len(old_param_list) == len(new_param_list)
        except AssertionError:
            print("Reshuffling a fisher matrix requires a new parameter list of the same length than initial one. Use fix method to shorten the param list")
        try :
            assert set(new_param_list).issubset(old_param_list)
        except AssertionError:
            print("The new parameter list needs to contain the same parameters than the old one !")
        indices = [old_param_list.index(param) for param in new_param_list]
        tmp = self.fisher
        self.fisher = temp[np.ix_(indices, indices)]
        self.param_list = new_param_list
                

        



def read_cl(data_root, lmax, n_freqs, parameter, direction=None):
    """ Function that gets the power spectra for a given parameter and a given direction. Uses pandas to speed up things a bit.
        - data_root : directory containing the power pectra.
        - n_freqs : number of frequencies used by the experiment
        - parameter : parameter that we are interested in or 'fiducial' if we want to get fiducial power spectra
        - direction : left or right. If None, then read fiducial power spectra.
        returns : - cls (lmax-1 * *n_freqs + 1 * n_freqs + 1 * 5 (TT+TE+EE+TP+PP) ) TP and PP data are only taken for primary channels (index 0)
    """
    pandas = True
    try:
        import pandas as pd
    except ImportError:
        print("pandas not found, using numpy instead")
        pandas = False
    
    if direrction is not None:
        file_name = os.path.join(data_root,"{}_{}_".format(parameter,direction)) 
    else:
        file_name = os.path.join(data_root,"{}_".format("fiducial"))
    cls = np.zeros((lmax-1,n_freqs+1,n_freqs+1,5))
    for i in range(n_freqs+1):
        for j in range(n_freqs+1):
            if pandas:
                data_read = pd.read_csv(filename + 'scalCls_lensed_{:d}_{:d}'.format(i+1,j+1), header = None, skiprows = 0, sep = '\s+').values
            else:
                data_read = np.loadtxt(filename + 'scalCls_lensed_{:d}_{:d}'.format(i+1,j+1))
            cls[:,i,j,0] = data_read[0:lmax-1,1] # TT all spectra starts at ell = 2
            cls[:,i,j,1] = data_read[0:lmax-1,2] # EE
            cls[:,i,j,2] = data_read[0:lmax-1,4] # TE
    if pandas:
        data_read_phi = pd.read_csv(filename + 'scalCls_lenspot.dat', header = None, skiprows = 1, sep = '\s+').values # reading lensing potential values, files have a header
    else:
        data_read_phi = np.loadtxt(filename + 'scalCls_lenspot.dat')
    cls[:,0,0,3] = data_read_phi[0:lmax-1,5] # PP only the primary channel is completed (for now)
    cls[:,0,0,4] = data_read_phi[0:lmax-1,6] # TP
    return data

def kron_delta(i,j):
    """ Kronecker delta symbol, returns 1 if i==j, 0 otherwise """
    if i == j:
        return 1
    else:
        return 0
    
def compute_covariance(experiment,cls):
    """ Function to compute the covariance matrix given set of cls and an experiment
        - experiment : corresponding CMB experiment, used to know whether lensing, rayleigh or polarization should be included
        - cls : cls used to populate the covariance matrix
        return : - cov : covariance matrix
    """
    lmax = max(experiment.lmax_T,experiment.lmax_P)
    n_freqs = len(experiment.freqs)
    # Matrix is constructed with lensing and lensing is removed at the end if unwanted
    if experiment.include_rayleigh:
        if experiment.include_P:
            cov = np.zeros((lmax-1,2*(n_freqs+1)+1,2*(n_freqs+1)+1))
            for i in range(n_freqs+1):
                for j in range(n_freqs+1):
                    cov[:,i,j] = cls[:,0,0,0] + kron_delta(i,j)*experiment.NlTT[:,i+1] + (1.-kron_delta(i,0))*cls[:,i,0,0] + (1.-kron_delta(j,0))*cls[:,0,j,0] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,0] # TT part
                    cov[:,i+n_freq+1,j+n_freq+1] = cls[:,0,0,1] + kron_delta(i,j)*experiment.NlTEE[:,i+1] + (1.-kron_delta(i,0))*cls[:,i,0,1] + (1.-kron_delta(j,0))*cls[:,0,j,1] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,1] # EE part
                    cov[:,i+n_freq+1,j] = cls[:,0,0,2] + (1.-kron_delta(i,0))*cls[:,i,0,2] + (1.-kron_delta(j,0))*cls[:,0,j,2] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,2] # TE part
                    cov[:,i,j+n_freq+1] = cls[:,0,0,2] + (1.-kron_delta(i,0))*cls[:,i,0,2] + (1.-kron_delta(j,0))*cls[:,0,j,2] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,2] # TE part
        else:
            cov = np.zeros((l_max-1,2*(n_freqs+1)+1))
            for i in range(n_freqs+1):
                for j in range(n_freqs +1):
                    cov[:,i,j] = cls[:,0,0,0] + kron_delta(i,j)*experiment.NlTT[:,i+1] + (1.-kron_delta(i,0))*cls[:,i,0,0] + (1.-kron_delta(j,0))*cls[:,0,j,0] + (1-kron_delta(i,0))*(1-kron_delta(j,0))*cls[:,i,j,0] # TT part
    else:
        if experiment.incude_P:
            cov = np.zeros((3,3))
            cov[:,0,0] = cls[:,0,0,0] + experiment.NlTT[:,1] 
            cov[:,0,1] = cls[:,0,0,2]
            cov[:,1,0] = cls[:,0,0,2]
            cov[:,1,1] = cls[:,0,0,1] + experiment.NlEE[:,1]
        else:
            cov = np.zeros((2,2))
            cov[:,0,0] = cls[:,0,0,0] + experiment.NlTT[:,1]
    if experiment.include_lensing:
        cov[:,-1,0] = cls[:,0,0,4] # TP part
        cov[:,0,-1] = cls[:,0,0,4] # TP part
        cov[:,-1,-1] = cls[:,0,0,3] #PP part
        return cov
    else:
        return cov[:,0:-1,0:-1]
        

def derivative(setup,experiment,data_root,parameter):
    """ Function to compute the derivative of the the covariance matrix with respect to a given parameter.
        - setup :  
        - experiment : corresponding CMB experiment, used to get the number of frequency channels, lmax
        - data_root : directory containing the powers pectra
        - parameter : parameter with respect to which take the derivative
        return : deriv
    """ #TO COMPLETE
    lmax = max(experiment.lmax_T,experiment.lmax_P)
    n_freqs = len(experiment.freqs)
    #Reading data in every direction
    cls_right = read_cl(data_root, lmax, n_freqs, parameter, direction='right')
    cls_left = read_cl(data_root, lmax, n_freqs, parameter, direction='left')
    #Compute the corresponding covariance matrices
    cov_right = compute_covariance(experiment,cls_right)
    cov_left = compute_covariance(experiment,cls_left)
    #Compute the derivative (2pts for now, maybe update at some point ?)
    deriv = (cov_right - cov_left)/(2.*setup.step[parameter])
    return deriv
    
def fish_trace(lmax,deriv_i,deriv_j,cov_fid):
    """ Function that comput the trace part of the fisher formula.
        - lmax : ell max
        - deriv_i : derivative of the covariance matrix with respect to parameter_i
        - deriv_j : derivative of the covariance matrix with respect to parameter_j
        - cov_fid : fiducial covariance matrix
        return : tr(inv(cov_fid)*deriv_i*inv(cov_fid)*deriv_j) array of lmax-1 entries
    """
    
    tr = np.zeros(lmax-1)
    for ell in range(lmax-1):
        inv = np.linalg.inv(cov_fid[ell,:,:])
        tr = np.trace(np.dot(inv,np.dot(deriv_i,np.dot(inv,deriv_j))))
    return tr
    
        
    
             
   
   
        
        
        
        
        
        
        
        
        
        

