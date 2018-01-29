#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
from scipy import stats

#PRISM data are form : https://arxiv.org/pdf/1306.2259.pdf
# PLANCK data are from : https://arxiv.org/pdf/astro-ph/0604069.pdf

T_CMB = 2.7255
freq = {}
angle_beam = {}
noise_pix = {}
ndet = {}

freq['PLANCK'] =  (0, 143, 217, 353, 545, 857) #in GHz
freq['PRISM'] = (0,295,300,320,395,460,555,660,800) #in GHz
freq['CCAT'] = (0,150,226,273,350,405,862) #in GHz
freq['CMBS4'] = (0,95,145,220,270) #in GHz
freq['CCAT-SO'] = (0,90,150,220,270,350,405,862) # in GHz
freq['CCAT-S4'] = (0,95,145,220,270,350,405,862) # in GHz
freq['SO'] = (0,90,150,220,270) # in GHz
freq['CCAT-PLANCK'] = freq['CCAT']

ndet['PRISM'] =  np.asarray((350,350,350,350,350,350,300,300,200)) #prism


angle_beam['PLANCK'] = np.asarray((4.89,7.1,5.0,5.0,5.0,5.0))  #in arcmin
angle_beam['PRISM'] = np.asarray((2.5, 1.7, 1.7, 1.6, 1.3, 1.1, 55./60, 46./60, 38./60)) #in arcmin 
angle_beam['CCAT'] = np.asarray((1.4,1.4,1.0,0.8,0.6,0.5,0.3)) #in arcmin 
angle_beam['CMBS4'] = np.asarray((1.4,2.2,1.4,1,0.8)) #in arcmin
angle_beam['SO'] = np.asarray((1.8,2.8,1.8,1.2,1))
angle_beam['CCAT-S4'] = np.asarray((1.4,2.2,1.4,1,0.8,0.6,0.5,0.3))
angle_beam['CCAT-SO'] = np.asarray((1.8,2.8,1.8,1.2,1,0.6,0.5,0.3))
angle_beam['CCAT-PLANCK'] = np.asarray((4.89,1.4,1.0,0.8,0.6,0.5,0.3)) #in arcmin

l_min = {'PLANCK' : 2, 'CCAT' : 30, 'SO' : 30, 'CMBS4' : 30, 'GRID' : 30}


noise_pix['PLANCK'] = np.asarray((5.148,2.2,4.8,14.7,500,6700))*T_CMB*angle_beam['PLANCK']*np.pi/10800 #in muKrad
noise_pix['PRISM'] = np.asarray((67.1, 67.1, 67.1, 73.2, 107.,156., 297., 700., 20000.))/np.sqrt(ndet['PRISM'])*np.pi/10800  #in muKrad
noise_pix['CCAT'] = np.asarray((5.,6.,5.,6.,30.,70.,70000.))*np.pi/10800  #in muKrad
noise_pix['CMBS4'] = np.asarray((1.,0.7,1.,5.,5.))*np.pi/10800 #in muKrad
noise_pix['SO'] = np.asarray((5.2,4.84,5.20,15.47,53.0))*np.pi/10800
noise_pix['CCAT-SO'] = np.asarray((5.2,4.84,5.20,15.47,53.0,30.,70.,70000.))*np.pi/10800 #in muKrad
noise_pix['CCAT-S4'] = np.asarray((1.,0.7,1.,5.,5.,30.,70.,70000.))*np.pi/10800 #in muKrad
noise_pix['CCAT-PLANCK'] = np.asarray((68.6,6.,5.,6.,30.,70.,70000.))*np.pi/10800  #in muKrad

def noise(experiment,l_max) :
    n = len(freq[experiment])
    l = np.linspace(2,l_max,l_max-1)
    var = np.zeros((l_max-1, n+1))
    for i in range(n) : 
        var[:,i+1] = l*(l+1.)/(2*np.pi) * noise_pix[experiment][i] **2 * np.exp(l*(l+1.)*(angle_beam[experiment][i]*np.pi/10800)**2 /(8*np.log(2)))

    for ll in range(l_min[experiment]-2) : 
        var[ll,1:n+1] = 1E12
    var[:,0] = l
    return var

np.savetxt('../noise/noise_PLANCK.txt',noise('PLANCK',4000), delimiter = '   ', newline = '\n')
np.savetxt('../noise/noise_SO.txt',noise('SO',4000), delimiter = '   ', newline = '\n')
np.savetxt('../noise/noise_CMBS4.txt',noise('CMBS4',4000), delimiter = '   ', newline = '\n')
np.savetxt('../noise/noise_CCAT.txt',noise('CCAT',4000), delimiter = '   ', newline = '\n')

def co_add(list_experiment,list_freq_to_keep,header = 'None') : 
    l = noise(list_experiment[0],4000)[:,0]
    n = len(list_experiment) 
    freq_list = []
    for i in range(n) :
        freq_list += list_freq_to_keep[i]
    freq_list_s = list(sorted(freq_list))
    freq_list_f = freq_list_s.copy()
    index = []
    for i in range(len(freq_list_s)-1) :
        if (freq_list_s[i+1]- freq_list_s[i]) <= 10 :
             index.append(i)
             freq_list_f.pop(freq_list_f.index(freq_list_s[i]))
    
    noise_tot = np.zeros((np.shape(l)[0],len(freq_list_f)+1))
    for i in range(len(freq_list_f)) :
        noise_min1 = np.zeros((np.shape(l)[0]))
        for j in range(n) :
            noi = noise(list_experiment[j],4000)
            fr = list_freq_to_keep[j]
            arg = min(fr,key = lambda x:abs(x-freq_list_f[i]))
            if abs(arg-freq_list_f[i]) < 10 :
                ind = freq[list_experiment[j]].index(arg)
                noise_min1[:] += 1./noi[:,ind+1]
        noise_tot[:,i+1] = 1./noise_min1
    noise_tot[:,0] = l
    if header == 'None' :
        str1 = ''
    else :
        str1 = '_' + header
    file_name = '../noise/noise{}_'.format(str1) + '_'.join(list_experiment) +'.txt'
    header = '+'.join(list_experiment) + ' : ' + str(freq_list_f)    
    np.savetxt(file_name, noise_tot, header = header, delimiter = '  ', newline = '\n')   
    return

co_add(('CCAT','PLANCK'),((0,150,226,273,350,405,862),[0]))
co_add(('CCAT','SO'),((150,226,273,350,405,862),(0,90,150,220,270)))
co_add(('SO','PLANCK'),((0,90,150,220,270),[0]))
co_add(('CMBS4','PLANCK'),((0,95,145,220,270),[0]))
co_add(('CCAT','SO','PLANCK'),((150,226,273,350,405,862),(0,90,150,220,270),[0]))
co_add(('CCAT','CMBS4','PLANCK'),((150,226,273,350,405,862),(0,95,145,220,270),[0]))


def interpolate(experiment) : 
    fr = freq[experiment]
    noise = noise_pix[experiment]*10800/np.pi
    tck = interp.splrep(fr,np.log10(noise),s=0,k=min(len(fr)-1,2))
    fr_new = np.linspace(fr[0],900,1000)
    noise_new = np.power(10,interp.splev(fr_new,tck,der=0))
    return fr_new,noise_new

def linear_reg(experiment) :
    fr = freq[experiment]
    noise = noise_pix[experiment]*10800/np.pi
    slope, intercept, r_value, p_value, std_err = stats.linregress(fr[1:-1],np.log10(noise[1:-1]))
    fr_new = np.linspace(fr[1],1000,1000)
    noise_new = np.power(10,slope*fr_new + intercept)
    pol = np.polyfit(fr[1:-1],np.log(angle_beam[experiment][1:-1]),deg = 1)
    angle_new = np.exp(np.polyval(pol,fr_new))  
    print('Slope for {} is : {} \n Intercept for {} is : {}'.format(experiment,slope,experiment,intercept))
    print('fir for angle : {}'.format(pol))
    return fr_new,noise_new,angle_new

def gen_noise(intercept) :
    fr = np.asarray((0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900))
    slope = 0.0055
    noise_gr = np.power(10,slope*fr+intercept)
    pol = [-0.0040535,0.92009514]
    angle = np.exp(np.polyval(pol,fr))
    noise_pix['GRID'] = noise_gr*np.pi/10800
    noise_pix['GRID'][0] = 3.*np.pi/10800
    noise_pix['GRID'][1] = noise_pix['GRID'][2]
    angle_beam['GRID'] = angle
    angle_beam['GRID'][0] = angle_beam['GRID'][3]
    angle_beam['GRID'][1] = angle_beam['GRID'][2]
    freq['GRID'] = fr
    np.savetxt('../noise/noise_grid_3muK.txt',noise('GRID',4000),delimiter = '   ', newline = '\n')
    return list(fr),noise_pix['GRID'],angle_beam['GRID']

freq['GRID'],noise_pix['GRID'],angle_beam['GRID'] = gen_noise(-0.5)

co_add(('GRID','PLANCK'),((0,50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900),[0]), header = 'full')
co_add(('GRID','PLANCK'),((0,50,100,150,200,250),[0]), header = '250')
co_add(('GRID','PLANCK'),((0,250,300,350,400,450,500),[0]), header = '250_500')
co_add(('GRID','PLANCK'),((0,500,550,600,650,700,750),[0]), header = '500_750')
co_add(('GRID','PLANCK'),((0,750,800,850,900),[0]), header = '750_900')


color = {'SO' : 'r', 'CMBS4' : 'b', 'CCAT' : 'g','CCAT' : 'y','PLANCK' : 'k'}
f,axs = plt.subplots(1,2,figsize = (15,10))
for experiment in ['SO','CMBS4','CCAT','PLANCK'] :
    axs[0].scatter(freq[experiment],noise_pix[experiment]*10800/np.pi,c=color[experiment],marker = '*', s = 150,label = experiment)
    axs[1].scatter(freq[experiment], angle_beam[experiment],c = color[experiment],marker = '*',s = 150,label = experiment)

fr_grid, noise_grid, angle_grid = gen_noise(-0.5)
axs[0].scatter(fr_grid,noise_grid*10800/np.pi, c=(1,0.55,0), marker = 'o',s = 150, label = 'Grid test')
axs[1].scatter(fr_grid,angle_grid, c=(1,0.55,0), marker = 'o',s = 150, label = 'Grid test')
axs[0].set_title('Noise per pixel')
axs[1].set_title('Beam width')
axs[0].set_xlabel('Frequency [GHz]')
axs[0].set_ylabel('Noise per pixel [muKarcmin]')
axs[1].set_xlabel('Frequency [GHz]')
axs[1].set_ylabel('FWHM [arcmin]')
axs[0].set_yscale('log')
axs[1].set_ylim(0.1,6)
axs[0].set_ylim(1e-1,1e5)
axs[0].set_xlim(0,900)
axs[1].set_xlim(0,900)
axs[0].legend(loc = 'upper left') 
axs[1].legend(loc = 'upper left')
axs[0].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
axs[1].grid(color='k', linestyle = '-.', linewidth = 0.5, alpha = 0.5)
plt.savefig('/home/bb510/Documents/PhD/Rayleigh/fig/noise_grid.pdf', format = 'pdf')
plt.show()

