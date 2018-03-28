import time
import pandas as pd
import numpy as np
import configparser


inifile = "/home/bb510/Code/fishpy/input/test_in.ini"

config = configparser.SafeConfigParser()
config.read(inifile)        
sections = config.sections()
print(sections)

lmax = config.getint('general_parameters','l_max')
include = config.getboolean('experiment_1','include')
list1 = config.get('general_parameters','parameters')
print(lmax,include)
print(isinstance(list1,str))

print(list1[1:-1])

def get_list_from_str(string):
    n_elem = string.count(',') + 1
    
