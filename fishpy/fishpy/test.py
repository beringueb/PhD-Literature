import numpy as np

ell = np.linspace(0,10,11)
col1 = np.linspace(51,60,11)
col2 = np.linspace(23,25,11).reshape(11,1)

tot = np.concatenate([ell.reshape(11,1),col1.reshape(11,1),col2], axis = 1)
numcol  = 2
filename = "test.txt"
print(np.shape(tot))

tot[2:5,1:] = np.inf

np.savetxt(filename,tot, fmt = "%d  " + 2*"%3.1f  ", newline = "\n")


