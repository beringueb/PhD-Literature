#/usr/bin/env python
# --
# get_nl_phi.py
# --
# calculates and plots the lensing reconstruction noise levels
# for a given experiment. Noise levels are given as an input. For now fiducial theory cl's are used. Do we need to get the cl's for the experiement ? 

import numpy as np
import pylab as pl
import os
import quicklens as ql

QUICKLENS_LOC = "/home/bb510/Code/quicklens/"


#List of quadratic estimators
lensing_list  = ['TT', 'EB']

print("Estimating lensing noise from {} ...".format(str(lensing_list)) )

watch = ql.util.stopwatch()
def calc_nlqq(qest, clXX, clXY, clYY, flX, flY):
    errs = np.geterr(); np.seterr(divide='ignore', invalid='ignore')
    
    print( "[%s]"%watch.elapsed(), "calculating flat-sky noise level for estimator of type", type(qest))
    clqq_flatsky = qest.fill_clqq(ql.maps.cfft(nx,dx), clXX*flX*flX, clXY*flX*flY, clYY*flY*flY)
    resp_flatsky = qest.fill_resp(qest, ql.maps.cfft(nx, dx), flX, flY)
    nlqq_flatsky = clqq_flatsky / resp_flatsky**2
    
    print("[%s]"%watch.elapsed(), "calculating full-sky noise level for estimator of type", type(qest))
    clqq_fullsky = qest.fill_clqq(np.zeros(lmax+1, dtype=np.complex), clXX*flX*flX, clXY*flX*flY, clYY*flY*flY)
    resp_fullsky = qest.fill_resp(qest, np.zeros(lmax+1, dtype=np.complex), flX, flY)
    nlqq_fullsky = clqq_fullsky / resp_fullsky**2

    np.seterr(**errs)
    return nlqq_flatsky, nlqq_fullsky

NlTT = np.loadtxt(os.path.join(QUICKLENS_LOC,"temp_fishpy","NlTT.temp"))
NlEE = np.loadtxt(os.path.join(QUICKLENS_LOC,"temp_fishpy","NlEE.temp"))

lmax = int(NlTT[-1,0])
nx   = 512  # number of pixels for flat-sky calc.
dx   = 1./60./180.*np.pi # pixel width in radians
cl_unl = ql.spec.get_camb_scalcl(lmax=lmax) # unlensed theory spectra.
cl_len = ql.spec.get_camb_lensedcl(lmax=lmax) # lensed theory spectra.
pix    = ql.maps.pix(nx,dx)
n_freq = np.shape(NlTT)[1]-2

ell = np.linspace(2,lmax,lmax-1)
nltt = np.concatenate([np.zeros(2),NlTT[:,1]*2*np.pi/ell/(ell+1.)])
nlee = nlbb =  np.concatenate([np.zeros(2),NlEE[:,1]*2*np.pi/ell/(ell+1.)])

# signal spectra
sltt       = cl_len.cltt
slte       = cl_len.clte
slee       = cl_len.clee
slbb       = cl_len.clbb
zero       = np.zeros(lmax+1)

# signal+noise spectra
cltt       = sltt + nltt
clee       = slee + nlee
clbb       = slbb + nlbb

# filter functions
flt        = np.zeros( lmax+1 ); flt[2:] = 1./cltt[2:]
fle        = np.zeros( lmax+1 ); fle[2:] = 1./clee[2:]
flb        = np.zeros( lmax+1 ); flb[2:] = 1./clbb[2:]


# intialize quadratic estimators
qest_TT    = ql.qest.lens.phi_TT(sltt)
qest_EE    = ql.qest.lens.phi_EE(slee)
qest_TE    = ql.qest.lens.phi_TE(slte)
qest_TB    = ql.qest.lens.phi_TB(slte)
qest_EB    = ql.qest.lens.phi_EB(slee)

t = lambda l: (l*(l+1.))**2/(2.*np.pi)  # scaling to apply to cl_phiphi when plotting or saving 
noise_temp = np.zeros(lmax+1)
for lensing_channel in lensing_list:
    if lensing_channel == "TT":
        nlpp_TT_flatsky, nlpp_TT_fullsky = calc_nlqq( qest_TT, cltt, cltt, cltt, flt, flt )
        noise_temp += 1./np.real(nlpp_TT_fullsky)
    elif lensing_channel == "EE":
        nlpp_EE_flatsky, nlpp_EE_fullsky = calc_nlqq( qest_EE, clee, clee, clee, fle, fle )
        noise_temp += 1./np.real(nlpp_EE_fullsky)
    elif lensing_channel == "TE":
        nlpp_TE_flatsky, nlpp_TE_fullsky = calc_nlqq( qest_TE, cltt, slte, clee, flt, fle )
        noise_temp += 1./np.real(nlpp_TE_fullsky)
    elif lensing_channel == "TB":
        nlpp_TB_flatsky, nlpp_TB_fullsky = calc_nlqq( qest_TB, cltt, zero, clbb, flt, flb )
        noise_temp += 1./np.real(nlpp_TB_fullsky)
    elif lensing_channel == "EB":
        nlpp_EB_flatsky, nlpp_EB_fullsky = calc_nlqq( qest_EB, clee, zero, clbb, fle, flb )
        noise_temp += 1./np.real(nlpp_EB_fullsky)
    
noise_tot = t(ell) * 1./noise_temp[2:]
nlpp = np.concatenate([ell.reshape(lmax-1,1),noise_tot.reshape(lmax-1,1)], axis = 1)
np.savetxt(os.path.join(QUICKLENS_LOC,"temp_fishpy","NlPP.temp"), nlpp)

print("Estimating lensing noise done !")
















