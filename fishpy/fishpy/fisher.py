#!/usr/bin/python3 
""" Module that contains the definition of the fisher_matrix class and methods in it"""

import numpy as np
import os


class fisher_matrix():  
    """ Class that defines a fisher_matrix and methods useful to generate it, add two of them, reshuffle them, etc ... """
    
    def __init__(self,setup,experiment):
        """ Initilaization of the fisher matrices with useful parameters
            - setup :
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
        
    def write_to_txt(self, root)
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
        
        data_root = os.path.join(setup.data_root,experiment.name)
        
        
        
        
        
        
        
        
        
        

