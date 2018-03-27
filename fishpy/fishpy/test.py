import time
import pandas as pd
import numpy as np

filename = "/home/bb510/Code/CAMB/tests/test_4_"

time_start = time.time()
for i in range(5):
    for j in range(5):
        data_read = pd.read_csv(filename + '{:d}_{:d}'.format(i+1,j+1), header = None, skiprows = 0, sep = '\s+').values

time_stop = time.time()
#print("pandas took {:3.1f}".format((time_stop-time_start)))

time_start = time.time()                

for i in range(5):
    for j in range(5):
        data_read = np.loadtxt(filename + '{:d}_{:d}'.format(i+1,j+1)) 

time_stop = time.time()

#print("numpy took {:3.1f}".format((time_stop-time_start)))

import time
from datetime import datetime
n = 100000
spaces = ''
for i in range(n):
    date_format = "%m/%d/%Y %H:%M"
    today = datetime.now()
    a = datetime.strptime('3/31/2018 09:00', date_format)
    delta = a - today
    #print('Encore {} avant de retrouver ma beaute celeste =( !'.format(delta), end = '\r')
    #print("A quel point je t'aime : |{}|".format(spaces), end = '\r')
    #print("Pourcentage d'amour : Elsa 100!%, Benjamin {:d}%" .format(i+100),end = '\r')
    print('Qui est la meilleure ?')
    a = input()
    if a.lower() == 'elsa':
        print("Oui tout a fait, elle est vraiment trop forte, c'est pour ca que je l'aime !!")
        break
    else:
        print("Non, debile, essaye encore !!")
    if i % 10 == 0:
        spaces += ' '
    time.sleep(0.05) 

