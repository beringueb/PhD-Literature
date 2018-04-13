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


def include_delensing_noise(experiment) : 
    """Read noise file, compute delensing noise for TT,TE,EE,TB,EB and append that at the end of the file"""
    file_name = "/home/bb510/Code/Rayleigh/noise/SO/noise_%s.dat" % experiment
    lmax = 5000
    print(file_name)
    try : 
        noise_in = np.loadtxt(file_name)[0:lmax+1,:]
    except : 
        print("File doesn't exist, need to generate it first !")
    
    cl_unl = ql.spec.get_camb_scalcl(lmax=lmax) # unlensed theory spectra.
    cl_len = ql.spec.get_camb_lensedcl(lmax=lmax) # lensed theory spectra.
    pix    = ql.maps.pix(nx,dx)
    n_freq = (np.shape(noise_in)[1]-1)/2
    print(n_freq)
    l = np.linspace(2,lmax,lmax-1)
    ell = np.linspace(0,lmax,lmax+1)
    nltt = np.concatenate([np.zeros(2),noise_in[2:,1]*2*np.pi/l/(l+1.)])
    nlee = nlbb =  np.concatenate([np.zeros(2),noise_in[2:,1+n_freq]*2*np.pi/l/(l+1.)])

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
    nlpp_TT_flatsky, nlpp_TT_fullsky = calc_nlqq( qest_TT, cltt, cltt, cltt, flt, flt )
    #nlpp_EE_flatsky, nlpp_EE_fullsky = calc_nlqq( qest_EE, clee, clee, clee, fle, fle )
    #nlpp_TE_flatsky, nlpp_TE_fullsky = calc_nlqq( qest_TE, cltt, slte, clee, flt, fle )
    #nlpp_TB_flatsky, nlpp_TB_fullsky = calc_nlqq( qest_TB, cltt, zero, clbb, flt, flb )
    nlpp_EB_flatsky, nlpp_EB_fullsky = calc_nlqq( qest_EB, clee, zero, clbb, fle, flb )
    
    noise_tot = t(ell) * 1./(1./np.real(nlpp_TT_fullsky) + 1./np.real(nlpp_EB_fullsky))
    
    print(np.shape(noise_in))
    print(np.shape(noise_tot))
    print(np.shape(l))
    noise_write = np.concatenate([ell.reshape(lmax+1,1),noise_in[:,1:],noise_tot.reshape(lmax+1,1)], axis = 1)
    print(np.shape(noise_write))
    np.savetxt(file_name.replace('.dat','_lens.dat'), noise_write)

    return  nlpp_TT_fullsky, nlpp_EB_fullsky

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




if __name__ == "__main__" :
    lmax = 5000
    # make plot
    nx   = 512  # number of pixels for flat-sky calc.
    dx   = 1./60./180.*np.pi # pixel width in radians
    ls         = np.arange(0,lmax+1)
    t = lambda l: (l*(l+1.))**2/(2.*np.pi)  # scaling to apply to cl_phiphi when plotting or saving 
    lbins      = np.linspace(10, lmax, 100)       # multipole bins.
    
    #nlpp_TT_fullsky, nlpp_EB_fullsky = include_delensing_noise('CCAT_atmo_primary') 
    #nlpp_TT_fullsky, nlpp_EB_fullsky = include_delensing_noise('CCAT_noatmo_primary') 
    for fsky in [0.1,0.2,0.4] : 
        for sensi in [0,1,2] :
            nlpp_TT_fullsky, nlpp_EB_fullsky = include_delensing_noise('SO_%d_%3.1f' %(sensi,fsky))
            print('Done fksy : %3.1f, sensi %d: ' %(fsky, sensi))
            
    #cl_unl = ql.spec.get_camb_scalcl(lmax=lmax) # unlensed theory spectra.
    #cl_len = ql.spec.get_camb_lensedcl(lmax=lmax) # lensed theory spectra.

    #cl_unl.plot( 'clpp', t=t, color='k' )
    #ql.spec.cl2cfft(cl_unl.clpp, ql.maps.cfft(nx,dx)).get_ml(lbins, t=t).plot(color='gray', ls='--')

    #pl.plot( ls, t(ls) * nlpp_TT_fullsky, color='r', label=r'TT' )
    #nlpp_TT_flatsky.get_ml(lbins, t=t).plot(color='k', ls='--')

    #pl.plot( ls, t(ls) * nlpp_EE_fullsky, color='g', label=r'EE' )
    #nlpp_EE_flatsky.get_ml(lbins, t=t).plot(color='k', ls='--')

    #pl.plot( ls, t(ls) * nlpp_TE_fullsky, color='b', label=r'TE' )
    #nlpp_TE_flatsky.get_ml(lbins, t=t).plot(color='k', ls='--')

    #pl.plot( ls, t(ls) * nlpp_TB_fullsky, color='y', label=r'TB' )
    #nlpp_TB_flatsky.get_ml(lbins, t=t).plot(color='k', ls='--')

    #pl.plot( ls, t(ls) * nlpp_EB_fullsky, color='m', label=r'EB' )
    #nlpp_EB_flatsky.get_ml(lbins, t=t).plot(color='k', ls='--', label='Flatsky')

    pl.xlim(1,3900)
    pl.semilogy()
    pl.xscale('log')

    pl.legend(loc='upper left', ncol=2)
    pl.setp( pl.gca().get_legend().get_frame(), visible=False)

    pl.xlabel(r'$L$')
    pl.ylabel(r'$[L(L+1)]^2 C_L^{\phi\phi} / 2\pi$')

    #pl.ion()
    #pl.show()
