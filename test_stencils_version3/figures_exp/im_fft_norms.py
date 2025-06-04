#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
from   numpy import linalg      as la
import sys
#==============================================================================

#==============================================================================
# Plot Set
#==============================================================================
import matplotlib.pyplot       as plt
import matplotlib.ticker       as mticker    
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib              import ticker
from   matplotlib              import cm
#==============================================================================

#==============================================================================
# Signal Packages
#==============================================================================
import pywt
from numpy.fft import rfft, rfftfreq, fft, fftfreq
import scipy.signal
from scipy import integrate
from scipy.signal import hilbert,envelope
#==============================================================================

#==============================================================================
# Auxilary Functions
#==============================================================================
def distf4(x,y):
    
    locres = x-y
    
    while(locres>np.pi or locres<-np.pi):
    
        if(locres>np.pi): 
    
            locres = locres - 2*np.pi
                        
        elif(locres<-np.pi): 
                        
            locres = locres + 2*np.pi
                        
    return locres
#==============================================================================

#==============================================================================
# FFT Norm 1
#==============================================================================
def fftnorm1(recref,recnum,teste,wcut,cutf):

    normtype        = 2         
    dt              = (teste.tn/(recref.shape[0]-1))/1000         
    vtime           = np.linspace(0,dt*recref.shape[0]-1,recref.shape[0])
    sampling_period = np.diff(vtime).mean()

    xfref    = fftfreq(len(recref),sampling_period)
    yfref    = fft(recref)
    xfnum    = fftfreq(len(recnum),sampling_period)
    yfnum    = fft(recnum)
    nsamples = yfref.shape[0]

    xfref    = xfref[:nsamples//2]
    yfref    = yfref[:nsamples//2]       
    xfnum    = xfnum[:nsamples//2]
    yfnum    = yfnum[:nsamples//2]
    nsamples = yfref.shape[0]

    if(wcut==1):
                            
        yfref = yfref[xfref<cutf]
        xfnum = xfnum[xfref<cutf]
        yfnum = yfnum[xfref<cutf] 
        xfref = xfref[xfref<cutf]

    nsamples = yfref.shape[0]

    angleref = np.angle(yfref)
    anglenum = np.angle(yfnum)
    absref   = (2/nsamples)*np.abs(yfref)
    absnum   = (2/nsamples)*np.abs(yfnum) 
    difangle = np.zeros(nsamples)

    for m1 in range(0,nsamples):
                        
        alpha1 = anglenum[m1]
        alpha2 = angleref[m1]
        difangle[m1] = distf4(alpha1,alpha2)
        
    nmeasureangle1 = (difangle/np.pi)*(np.abs(absnum)/np.abs(np.amax(absnum)))
    nnorm1         = (1/np.sqrt(nsamples))*la.norm(nmeasureangle1,normtype)
    
    return nnorm1
#==============================================================================

#==============================================================================
# FFT Norm 2
#==============================================================================
def fftnorm2(recref,recnum,teste,wcut,cutf):

    dt              = (teste.tn/(recref.shape[0]-1))/1000         
    vtime           = np.linspace(0,dt*recref.shape[0]-1,recref.shape[0])
    sampling_period = np.diff(vtime).mean()
    
    widths     = np.geomspace(1, 1024, num=100)
    wavelet    = "cmor1.5-1.0"
    cwt_method = 'conv'
    normtype   = 2

    cwtmatr_ref, freqs_ref = pywt.cwt(recref, widths, wavelet, sampling_period=sampling_period,method=cwt_method)
    cwtmatr_num, freqs_num = pywt.cwt(recnum, widths, wavelet, sampling_period=sampling_period,method=cwt_method)

    if(wcut==1):

        cwtmatr_ref = cwtmatr_ref[freqs_ref<cutf]
        cwtmatr_num = cwtmatr_num[freqs_ref<cutf]
        freqs_num   = freqs_num[freqs_ref<cutf]
        freqs_ref   = freqs_ref[freqs_ref<cutf]

    cwtmatr_ref_norm       = np.abs(cwtmatr_ref)
    cwtmatr_num_norm       = np.abs(cwtmatr_num)  
    cwtmatr_ref_angle      = np.angle(cwtmatr_ref)
    cwtmatr_num_angle      = np.angle(cwtmatr_num)
    cwtmatr_dif_angle      = np.zeros((cwtmatr_ref_angle.shape[0],cwtmatr_ref_angle.shape[1]))
    
    for m1 in range(0,cwtmatr_ref_angle.shape[0]):

        for m2 in range(0,cwtmatr_ref_angle.shape[1]):
                            
            alpha1 = cwtmatr_num_angle[m1,m2]
            alpha2 = cwtmatr_ref_angle[m1,m2]
            cwtmatr_dif_angle[m1,m2] = distf4(alpha1,alpha2)
        
    nmeasureangle1 = (cwtmatr_dif_angle/np.pi)*(np.abs(cwtmatr_num_norm)/np.abs(np.amax(cwtmatr_num_norm)))
    nnorm1 = (1/np.sqrt(cwtmatr_ref_angle.shape[0]))*(1/np.sqrt(cwtmatr_ref_angle.shape[1]))*la.norm(nmeasureangle1,normtype)
    
    return nnorm1
#==============================================================================

#==============================================================================
# FFT Norm 3
#==============================================================================
def fftnorm3(recref,recnum,teste):

    dt       = (teste.tn/(recref.shape[0]-1))/1000         
    vtime    = np.linspace(0,dt*recref.shape[0]-1,recref.shape[0])
    normtype = 2

    analytic_signal_ref     = hilbert(recref)
    amplitude_envelope_ref  = np.abs(analytic_signal_ref)
    instantaneous_phase_ref = np.angle(analytic_signal_ref)
                   
    analytic_signal_num     = hilbert(recnum)
    amplitude_envelope_num  = np.abs(analytic_signal_num)
    instantaneous_phase_num = np.angle(analytic_signal_num)
    
    nsamples1               = instantaneous_phase_num.shape[0]                
    instantaneous_phase_dif = np.zeros(nsamples1)

    for m1 in range(0,nsamples1):
                        
        alpha1 = instantaneous_phase_num[m1]
        alpha2 = instantaneous_phase_ref[m1]
        instantaneous_phase_dif[m1] = distf4(alpha1,alpha2)                    

    nmeasureangle1 = (instantaneous_phase_dif/np.pi)*(np.abs(amplitude_envelope_num)/np.abs(np.amax(amplitude_envelope_num)))
    nnorm1 = (1/np.sqrt(nsamples1))*la.norm(nmeasureangle1,normtype)
    
    return nnorm1
#==============================================================================