#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
from   numpy import linalg      as la
import sys
#==============================================================================

#==============================================================================
# My Modules
#==============================================================================
import testes_opt              as ttopt
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
from numpy.fft import fft, fftfreq
from scipy import integrate
from scipy.signal import hilbert,envelope
#==============================================================================

#==============================================================================
# Signal Comparison
#==============================================================================
from im_fft_norms import fftnorm1
from im_fft_norms import fftnorm2
from im_fft_norms import fftnorm3
#==============================================================================

#==============================================================================
# Range of Parameters
#==============================================================================
vptype    = [1] 
vdxref    = [1,2,4,8]
vdtref    = [1,2,4,6]
vfreqref  = [1,2,3] 

nvptype      = len(vptype) 
nvdxref      = len(vdxref)
nvdtref      = len(vdtref)
nvfreqref    = len(vfreqref) 
cont_me      = 0
lconfig_me   = []
cont_glob    = 0

for k0 in range(0,nvptype):
    
    for k1 in range(0,nvdxref):
        
        for k2 in range(0,nvdtref):
            
            for k3 in range(0,nvfreqref):

                ptype      = vptype[k0] 
                dx_ref     = vdxref[k1]
                dt_ref     = vdtref[k2]
                freq_ref   = vfreqref[k3]
                factor_ref = 16
                
                print('')
                print('Number    = %d'%cont_me)
                print('Test Type = %d'%ptype)
                print('dx_ref    = %d'%dx_ref)
                print('dt_ref    = %d'%dt_ref)
                print('freq_ref  = %d'%freq_ref)
                print('')
                
                lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
                cont_me = cont_me + 1
#==============================================================================

#==============================================================================
# Teste Select
#==============================================================================
                if(ptype==1): 
                    
                    teste_ref  = ttopt.teste1_ref1(freq_ref,factor_ref)
                    teste      = ttopt.teste1(dx_ref,dt_ref,freq_ref)
                    testname   = 'Homogeneos Velocity Model'

                    xivx1      = 1000
                    xfvx1      = 5000
                    nvx1       = 21
                    vx1        = np.linspace(xivx1,xfvx1,nvx1)  
                    yhigh      = 3000
                    vy1        = yhigh*np.ones(vx1.shape[0])        
                    xpositionv = vx1
                    ypositionv = vy1

                    xpvplot = xpositionv
                    ypvplot = ypositionv
  
                if(ptype==2): 
                    
                    teste_ref  = ttopt.teste2_ref1(freq_ref,factor_ref)
                    teste      = ttopt.teste2(dx_ref,dt_ref,freq_ref)
                    testname   = 'Heterogeneos Velocity Model'

                    xivx1      = 1000
                    xfvx1      = 3000
                    nvx1       = 21
                    vx1        = np.linspace(xivx1,xfvx1,nvx1)  
                    yhigh      = 2000
                    vy1        = yhigh*np.ones(vx1.shape[0])        
                    xpositionv = vx1
                    ypositionv = vy1

                    xpvplot = xpositionv
                    ypvplot = ypositionv
  
                if(ptype==3): 
                
                    teste_ref  = ttopt.teste3_ref1(freq_ref,factor_ref)
                    teste      = ttopt.teste3(dx_ref,dt_ref,freq_ref)
                    testname   = 'SEG/EAGE 2D Salt Velocity Model'

                    xivx1      = 2000
                    xfvx1      = 10000
                    nvx1       = 21
                    vx1        = np.linspace(xivx1,xfvx1,nvx1)  
                    yhigh      = 50
                    vy1        = yhigh*np.ones(vx1.shape[0])        
                    xpositionv = vx1
                    ypositionv = vy1

                    xpvplot = xpositionv
                    ypvplot = ypositionv

                if(ptype==4): 
                    
                    teste_ref   = ttopt.teste4_ref1(freq_ref,factor_ref)
                    teste       = ttopt.teste4(dx_ref,dt_ref,freq_ref)
                    testname    = 'Marmousi Velocity Model'

                    xivx1      = 6000
                    xfvx1      = 11000
                    nvx1       = 21
                    vx1        = np.linspace(xivx1,xfvx1,nvx1)  
                    yhigh      = 20
                    vy1        = yhigh*np.ones(vx1.shape[0])        
                    xpositionv = vx1
                    ypositionv = vy1

                    xpvplot = xpositionv
                    ypvplot = ypositionv

                locsave = 'signals/boards/teste%d/dx%ddt%dfreq%d/'%(ptype,dx_ref,dt_ref,freq_ref)
#==============================================================================

#==============================================================================
# Vetores de Configurações
#==============================================================================
                print('Starting Configs!')
                vmethod  = ['spatte','spectetheta','dispte','specls','displs'] 
                nvmethod = len(vmethod)
                vmvalue  = [1,2,3,4,5,6,7,8]
                nvmvalue = len(vmvalue)
                
                total_configs = 0
                list_config   = []
                    
                for m1 in range(0,nvmethod):
                        
                    method = vmethod[m1]
                        
                    if(method == 'spatte'):
                            
                        vshape = ['cl']
                            
                    elif(method == 'spectetheta'):
                            
                        vshape = ['cl']
                        
                    elif(method == 'dispte' or method == 'specls' or method == 'displs'):
                            
                        #vshape = ['crb']
                        vshape = ['crb','csq'] 
                        #vshape = ['rb',crb','csq','sq']
                            
                    nvshape = len(vshape)
                        
                    for m2 in range(0,nvshape):
                        
                        shape = vshape[m2]
                    
                        for m3 in range(0,nvmvalue):
                                
                            mvalue = vmvalue[m3]
                                
                            if(shape == 'cl'):
                                    
                                vnvalue = [1]
                                
                            elif(shape == 'rb'):
                                    
                                vnvalue = [mvalue]
                                
                            elif(shape == 'crb'):
                            
                                vnvalue  = np.arange(1,mvalue+1)

                            elif(shape == 'csq'):
                            
                                vnvalue  = np.arange(0,mvalue+1)

                            # else:
                            
                            #     vnvalue  = np.arange(1,mvalue+1)
                                
                            nvnvalue = len(vnvalue)
                                
                            for m4 in range(0,nvnvalue):
                                    
                                nvalue = vnvalue[m4]
                                config  = (shape,method,mvalue,nvalue,total_configs)                            
                                total_configs = total_configs + 1
                                list_config.append(config)
                    
                nconfig     = len(list_config)
                print('Finishing Configs!')
#==============================================================================

#==============================================================================
# Parameters Relation
#==============================================================================
                config        = list_config[0]
                mshape        = config[0]        
                method        = config[1]        
                sou           = int(2*config[2]) 
                mvalue        = int(config[2])  
                nvalue        = int(config[3])
#==============================================================================                    

#==============================================================================
# Obtenção de Parâmetros para Open
#==============================================================================
                lglobal           = []
                plotfig           = 0
                plotnormsanalysis = 1

                for k in range(0,nconfig):
                    
                    print('')
                    print('Test with Stencil: %d'%(k))
                    
                    nptx    = teste.nptx    
                    npty    = teste.npty    
                    x0      = teste.x0      
                    y0      = teste.y0      
                    x1      = teste.x1   
                    y1      = teste.y1   
                    hxv     = teste.hx      
                    hyv     = teste.hy      
                    t0      = teste.t0      
                    tn      = teste.tn      
                    f0      = teste.f0      
                    CFL     = teste.CFL
                    jump    = teste.jump
                    dt0     = 0.5*dt_ref
                    nt      = int((tn-t0)*dt0)
                
                    parameters = []
                    parameters = [nptx,npty,x0,y0,x1,y1,hxv,hyv,t0,tn,f0,CFL,jump,dt0,nt] 
                    
                    config        = list_config[k]
                    shape         = config[0]        
                    teste.shape   = shape
                    method        = config[1]        
                    teste.method  = method
                    sou           = int(2*config[2]) 
                    teste.sou     = sou
                    mvalue        = int(config[2])  
                    teste.mvalue  = mvalue
                    nvalue        = int(config[3])  
                    teste.nvalue  = nvalue
                    print('shape: %s - method: %s - sou: %d - mvalue: %d - nvalue: %d'%(shape,method,sou,mvalue,nvalue))
                    
                    if(shape=='cl'):
                        
                        npt = 2*mvalue + 1
                        npe = 0 
                        
                    elif(shape=='crb'):
                        
                        if(mvalue%2==0):
                        
                            npt  = int((nvalue**2+4*mvalue+1))
                        
                        else:
                            
                            npt  = int((nvalue**2-1+4*mvalue+1))
                            
                        npe = npt - (4*mvalue+1) 
                    
                    elif(shape=='csq'):
                        
                        npt     = int(4*((2*mvalue-nvalue+1)*nvalue-nvalue)+4*mvalue+1)  
                        npe     = npt - (4*mvalue+1) 
                            
                    sou         = teste.sou    
                    mvalue      = teste.mvalue  
                    nvalue      = teste.nvalue  
                    mshape      = teste.shape   
                    method      = teste.method

                    # Cutoff Frequency

                    cutf = 150
                    wcut = 1

                    if(nvalue==1):

                        print('Case N=1')
                    
                        locopenref = '../data_save/teste%d/reffreq%d/'%(ptype,freq_ref)

                        sol_ref        = np.load("%ssolplotcut_%d_%d_%d_%d.npy"%(locopenref,vptype[k0],vdxref[k1],vdtref[k2],vfreqref[k3]))
                        rec_ref        = np.load("%srecrefcut_%d_%d_%d_%d.npy"%(locopenref,vptype[k0],vdxref[k1],vdtref[k2],vfreqref[k3]))
                        rec_select_ref = np.load("%srecselectcut_%d_%d_%d_%d.npy"%(locopenref,vptype[k0],vdxref[k1],vdtref[k2],vfreqref[k3]))

                        locopensol    = '../data_save/teste%d/dx%ddt%dfreq%d/'%(ptype,dx_ref,dt_ref,freq_ref)

                        sol            = np.load("%ssolplot_%s_%s_%d_%d.npy"%(locopensol,mshape,method,mvalue,nvalue))
                        rec            = np.load("%srec_%s_%s_%d_%d.npy"%(locopensol,mshape,method,mvalue,nvalue))
                        rec_select     = np.load("%srec_select_%s_%s_%d_%d.npy"%(locopensol,mshape,method,mvalue,nvalue))

                        nts = rec_select.shape[1]

                        nfftrecselect1     = 0.0
                        ncwtrecselect1     = 0.0
                        nhilbertrecselect1 = 0.0
                        
                        nfftrecselect2     = 0.0
                        ncwtrecselect2     = 0.0
                        nhilbertrecselect2 = 0.0

                        n1recselect1       = 0.0
                        n2recselect1       = 0.0
                        nmaxrecselect1     = 0.0

                        n1recselect2       = 0.0
                        n2recselect2       = 0.0
                        nmaxrecselect2     = 0.0

                        llist1 = []
                        llist2 = []
                        llist3 = []
                        llist4 = []
                        llist5 = []
                        llist6 = []
                
                        for m1 in range(0,nts):

                            recnumloc = rec_select[:,m1]
                            recrefloc = rec_select_ref[:,m1]

                            a1 = fftnorm1(recrefloc,recnumloc,teste,wcut,cutf)
                            a2 = fftnorm2(recrefloc,recnumloc,teste,wcut,cutf)
                            a3 = fftnorm3(recrefloc,recnumloc,teste)
                            a4 = la.norm(recnumloc-recrefloc,1)/la.norm(recrefloc,1)
                            a5 = la.norm(recnumloc-recrefloc,2)/la.norm(recrefloc,2)
                            a6 = la.norm(recnumloc-recrefloc,np.inf)/la.norm(recrefloc,np.inf)
                            
                            nfftrecselect1     = nfftrecselect1     + a1
                            ncwtrecselect1     = ncwtrecselect1     + a2
                            nhilbertrecselect1 = nhilbertrecselect1 + a3
                            n1recselect1       = n1recselect1       + a4
                            n2recselect1       = n2recselect1       + a5
                            nmaxrecselect1     = nmaxrecselect1     + a6
                        
                            llist1.append(a1)
                            llist2.append(a2)
                            llist3.append(a3)
                            llist4.append(a4)
                            llist5.append(a5)
                            llist6.append(a6)
                    
                        nfftrecselect1     = nfftrecselect1/nts
                        ncwtrecselect1     = ncwtrecselect1/nts
                        nhilbertrecselect1 = nhilbertrecselect1/nts
                        n1recselect1       = n1recselect1/nts
                        n2recselect1       = n2recselect1/nts
                        nmaxrecselect1     = nmaxrecselect1/nts
            
                        nfftrecselect2     = np.amax(np.array(llist1))
                        ncwtrecselect2     = np.amax(np.array(llist2))
                        nhilbertrecselect2 = np.amax(np.array(llist3))
                        n1recselect2       = np.amax(np.array(llist4))
                        n2recselect2       = np.amax(np.array(llist5))
                        nmaxrecselect2     = np.amax(np.array(llist6))

                        rec_selectnorm1   = la.norm(rec_select-rec_select_ref,1)/la.norm(rec_select_ref,1)
                        rec_selectnorm2   = la.norm(rec_select-rec_select_ref,2)/la.norm(rec_select_ref,2)
                        rec_selectnormmax = la.norm(rec_select-rec_select_ref,np.inf)/la.norm(rec_select_ref,np.inf)

                        llocal = []
                        llocal.append(mshape)
                        llocal.append(method)
                        llocal.append(mvalue)
                        llocal.append(nvalue)
                        llocal.append(nfftrecselect1)
                        llocal.append(ncwtrecselect1)
                        llocal.append(nhilbertrecselect1)
                        llocal.append(nfftrecselect2)
                        llocal.append(ncwtrecselect2)
                        llocal.append(nhilbertrecselect2)
                        llocal.append(n1recselect1)
                        llocal.append(n2recselect1)
                        llocal.append(nmaxrecselect1)
                        llocal.append(n1recselect2)
                        llocal.append(n2recselect2)
                        llocal.append(nmaxrecselect2)
                        llocal.append(rec_selectnorm1)
                        llocal.append(rec_selectnorm2)
                        llocal.append(rec_selectnormmax)
                        lglobal.append(llocal)

                        if(plotfig==1):

                            vtimejump  = np.linspace(0,teste.tn,sol.shape[0])
                            itime      = 10
                            timeselect = vtimejump[itime]

                            sol     = sol[itime,:,:]
                            sol_ref = sol_ref[itime,:,:]

                            rec     = rec[:,:]
                            rec_ref = rec_ref[:,:]

                            solnorm1   = la.norm(sol-sol_ref,1)/la.norm(sol_ref,1)
                            solnorm2   = la.norm(sol-sol_ref,2)/la.norm(sol_ref,2)
                            solnormmax = la.norm(sol-sol_ref,np.inf)/la.norm(sol_ref,np.inf)

                            recnorm1   = la.norm(rec-rec_ref,1)/la.norm(rec_ref,1)
                            recnorm2   = la.norm(rec-rec_ref,2)/la.norm(rec_ref,2)
                            recnormmax = la.norm(rec-rec_ref,np.inf)/la.norm(rec_ref,np.inf)

                            ts     = int(nts/2)
                            recref = rec_select_ref[:,ts]
                            recnum = rec_select[:,ts]

                            rec_selectnorm1   = la.norm(recnum-recref,1)/la.norm(recref,1)
                            rec_selectnorm2   = la.norm(recnum-recref,2)/la.norm(recref,2)
                            rec_selectnormmax = la.norm(recnum-recref,np.inf)/la.norm(recref,np.inf)

                            tns    = teste.tn/1000
                            nt     = recref.shape[0]
                            dt     = (tns/(nt-1))
                            vtime = np.linspace(0,tns,nt)

                            def distf4(x,y):

                                locres = x-y

                                while(locres>np.pi or locres<-np.pi):
                                
                                    if(locres>np.pi): 
                                
                                        locres = locres - 2*np.pi
                                
                                    elif(locres<-np.pi): 
                                
                                        locres = locres + 2*np.pi
                                
                                return locres

                            # FFT Analysis
                            
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

                            if(wcut==1):
                                
                                yfref = yfref[xfref<cutf]
                                xfnum = xfnum[xfref<cutf]
                                yfnum = yfnum[xfref<cutf] 
                                xfref = xfref[xfref<cutf]

                            nsamples = yfref.shape[0]

                            absref   = (2/nsamples)*np.abs(yfref)
                            absnum   = (2/nsamples)*np.abs(yfnum)
                            difabs   = absnum - absref
                            angleref = np.angle(yfref)
                            anglenum = np.angle(yfnum)

                            difangle = np.zeros(nsamples)

                            for m1 in range(0,nsamples):
                                
                                alpha1 = anglenum[m1]
                                alpha2 = angleref[m1]
                                difangle[m1] = distf4(alpha1,alpha2)

                            ncfft  = (difangle/np.pi)*(np.abs(absnum)/np.abs(np.amax(absnum)))
                            nncfft = (1/np.sqrt(nsamples))*la.norm(ncfft,2)

                            # CWT Analysis

                            widths     = np.geomspace(1,1024,num=100)
                            wavelet    = "cmor1.5-1.0"
                            cwt_method = 'conv'

                            cwtmatr_ref, freqs_ref = pywt.cwt(recref, widths, wavelet, sampling_period=sampling_period,method=cwt_method)
                            cwtmatr_num, freqs_num = pywt.cwt(recnum, widths, wavelet, sampling_period=sampling_period,method=cwt_method)

                            if(wcut==1):

                                cwtmatr_ref = cwtmatr_ref[freqs_ref<cutf]
                                cwtmatr_num = cwtmatr_num[freqs_ref<cutf]
                                freqs_num   = freqs_num[freqs_ref<cutf]
                                freqs_ref   = freqs_ref[freqs_ref<cutf]

                            cwtmatr_ref_norm  = np.abs(cwtmatr_ref)
                            cwtmatr_num_norm  = np.abs(cwtmatr_num)
                            cwtmatr_dif_norm  = cwtmatr_num_norm - cwtmatr_ref_norm
                            cwtmatr_ref_angle = np.angle(cwtmatr_ref)
                            cwtmatr_num_angle = np.angle(cwtmatr_num)
                            cwtmatr_dif_angle = np.zeros((cwtmatr_ref_angle.shape[0],cwtmatr_ref_angle.shape[1]))

                            for m1 in range(0,cwtmatr_ref_angle.shape[0]):

                                for m2 in range(0,cwtmatr_ref_angle.shape[1]):
                                    
                                    alpha1 = cwtmatr_num_angle[m1,m2]
                                    alpha2 = cwtmatr_ref_angle[m1,m2]
                                    cwtmatr_dif_angle[m1,m2] = distf4(alpha1,alpha2)
                    
                            nccwt  = (cwtmatr_dif_angle/np.pi)*(np.abs(cwtmatr_num_norm)/np.abs(np.amax(cwtmatr_num_norm)))
                            nnccwt = (1/np.sqrt(cwtmatr_ref_angle.shape[0]))*(1/np.sqrt(cwtmatr_ref_angle.shape[1]))*la.norm(nccwt,2)
                            
                            # Hilbert Analysis

                            analytic_signal_ref         = hilbert(recref)
                            analytic_signal_num         = hilbert(recnum)
                            amplitude_envelope_ref      = np.abs(analytic_signal_ref)
                            instantaneous_phase_ref     = np.angle(analytic_signal_ref)
                            instantaneous_frequency_ref = -(10**3*teste.f0/(2*np.pi))*np.diff(instantaneous_phase_ref)
                            amplitude_envelope_num      = np.abs(analytic_signal_num)
                            instantaneous_phase_num     = np.angle(analytic_signal_num)
                            instantaneous_frequency_num = -(10**3*teste.f0/(2*np.pi))*np.diff(instantaneous_phase_ref)
                            amplitude_envelope_dif      = amplitude_envelope_num-amplitude_envelope_ref
                            nsamples1                   = instantaneous_phase_num.shape[0]                
                            instantaneous_phase_dif     = np.zeros(nsamples1)

                            for m1 in range(0,nsamples1):
                                
                                alpha1 = instantaneous_phase_num[m1]
                                alpha2 = instantaneous_phase_ref[m1]
                                instantaneous_phase_dif[m1] = distf4(alpha1,alpha2)                    

                            nchilbert  = (instantaneous_phase_dif/np.pi)*(np.abs(amplitude_envelope_num)/np.abs(np.amax(amplitude_envelope_num)))
                            nnchilbert = (1/np.sqrt(nsamples1))*la.norm(nchilbert,2)

                            columns1   = ['Norm Type','Values']
                            cell_text1 = [['1 Norm Relative',np.around(solnorm1,5)],['2 Norm Relative',np.around(solnorm2,5)],['Max Norm Relative',np.around(solnormmax,5)]]
                            
                            columns2   =['Norm Type','Values','Norm Type','Values']
                            cell_text2 = [['1 Norm Relative',np.around(recnorm1,5)],['2 Norm Relative',np.around(recnorm2,5)],['Max Norm Relative',np.around(recnormmax,5)]]
                            
                            columns3   = ['Norm Type','Values']
                            cell_text3 = [['1 Norm Relative',np.around(rec_selectnorm1,5)],['2 Norm Relative',np.around(rec_selectnorm2,5)],['Max Norm Relative',np.around(rec_selectnormmax,5)]]
                            
                            columns4   = ['Norm New Curve','Values']  
                            cell_text4 = [['FFT RMSE Norm',np.around(nncfft,5)],['CWT RMSE Norm',np.around(nnccwt,5)],['Hilbert RMSE Norm',np.around(nnchilbert,5)]]

                            columns5   = ['Norms','Values']  
                            cell_text5 = [['FFT RMSE - AVG',np.around(nfftrecselect1,5)],['CWT RMSE - AVG',np.around(ncwtrecselect1,5)],['Hilbert RMSE - AVG',np.around(nhilbertrecselect1,5)]]

                            columns6   = ['Norms','Values']  
                            cell_text6 = [['FFT RMSE - MAX',np.around(nfftrecselect2,5)],['CWT RMSE - MAX',np.around(ncwtrecselect2,5)],['Hilbert RMSE - MAX',np.around(nhilbertrecselect2,5)]]

                            columns7   = ['Norm Type','Values']
                            cell_text7 = [['1 Norm - AVG',np.around(n1recselect1,5)],['2 Norm - AVG',np.around(n2recselect1,5)],['Max Relative - AVG',np.around(nmaxrecselect1,5)]]
                                                
                            columns8   = ['Norm Type','Values']
                            cell_text8 = [['1 Norm - MAX',np.around(n1recselect2,5)],['2 Norm - MAX',np.around(n2recselect2,5)],['Max Relative - MAX',np.around(nmaxrecselect2,5)]]
                      
                            plt.figure(figsize = (50,20))
                            plt.suptitle('Signal Analysis - %s \n Shape = %s - Method = %s - M = %d - N = %d \n hx = %.2f m - hy = %.2f m - dt = %.4f s'%(testname,mshape,method,mvalue,nvalue,teste.hx,teste.hy,dt),fontsize=20)
                            grid  = plt.GridSpec(6,7,wspace=0.15,hspace=0.6)   

                            # Set of Receivers

                            plt.subplot(grid[0,5])
                            plt.box(on=None)
                            plt.title('New Receivers Set Metric - NOR = %d'%nts,y=1.05,fontsize=15)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table5 = plt.table(cellText=cell_text5,colLabels=columns5,loc='center',cellLoc='center')
                            table5.set_fontsize(20)
                            table5.scale(1.0, 2.5)

                            plt.subplot(grid[0,6])
                            plt.box(on=None)
                            plt.title('Clasic Receivers Set Metric - NOR = %d'%nts,y=1.05,fontsize=15)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table7 = plt.table(cellText=cell_text7,colLabels=columns7,loc='center',cellLoc='center')
                            table7.set_fontsize(20)
                            table7.scale(1.0, 2.5)
                
                            plt.subplot(grid[1,5])
                            plt.box(on=None)
                            plt.title('New Receivers Set Metric - NOR = %d'%nts,y=1.05,fontsize=15)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table6 = plt.table(cellText=cell_text6,colLabels=columns6,loc='center',cellLoc='center')
                            table6.set_fontsize(20)
                            table6.scale(1.0, 2.5)

                            plt.subplot(grid[1,6])
                            plt.box(on=None)
                            plt.title('Classic Receivers Set Metric - NOR = %d'%nts,y=1.05,fontsize=15)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table8 = plt.table(cellText=cell_text8,colLabels=columns8,loc='center',cellLoc='center')
                            table8.set_fontsize(20)
                            table8.scale(1.0, 2.5)

                            plt.subplot(grid[2,6])
                            fig1 = plt.plot(xpvplot/1000,ypvplot/1000)   
                            plt.grid()
                            plt.title('Receivers Distribution - NOR = %d'%nts,y=1.05,fontsize=15)
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            
                            # Displacement

                            vmin = min(np.amin(sol_ref),np.amin(sol),np.amin(sol_ref-sol))
                            vmax = max(np.amax(sol_ref),np.amax(sol),np.amax(sol_ref-sol))
                            scale1 = 10**(-3)
                            extent = [scale1*teste.x0,scale1*teste.x1,scale1*teste.y1,scale1*teste.y0]
                            factor = 50
                            scale  = vmax/factor

                            plt.subplot(grid[0,0])
                            fig1 = plt.imshow(np.transpose(sol_ref),cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                            plt.grid()
                            plt.title('Displacement - Reference - T = %.3f s'%(timeselect/1000))
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            ax = plt.gca()
                            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                            cbar.ax.locator_params(nbins=5)
                            tick_locator = ticker.MaxNLocator(nbins=5)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                            plt.subplot(grid[0,1])
                            fig1 = plt.imshow(np.transpose(sol),cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale) 
                            plt.grid()
                            plt.title('Displacement - Numerical - T = %.3f s'%(timeselect/1000))
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            ax = plt.gca()
                            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                            cbar.ax.locator_params(nbins=5)
                            tick_locator = ticker.MaxNLocator(nbins=5)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                            plt.subplot(grid[0,2])
                            fig1 = plt.imshow(np.transpose(sol_ref-sol),cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                            plt.grid()
                            plt.title('Displacement - Difference - T = %.3f s'%(timeselect/1000))
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            ax = plt.gca()
                            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                            cbar.ax.locator_params(nbins=5)
                            tick_locator = ticker.MaxNLocator(nbins=5)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                            plt.subplot(grid[0,3])
                            plt.box(on=None)
                            plt.title('Displacement - Classic Metrics - T = %.3f s'%(timeselect/1000),y=1.05,fontsize=15)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table1 = plt.table(cellText=cell_text1,colLabels=columns1,loc='center',cellLoc='center')
                            table1.set_fontsize(20)
                            table1.scale(1.0, 2.5)

                            # Seismogram

                            vmin = min(np.amin(rec),np.amin(rec_ref),np.amin(rec_ref-rec))
                            vmax = max(np.amax(rec),np.amax(rec_ref),np.amax(rec_ref-rec))
                            scale1 = 10**(-3)
                            extent = [scale1*teste.x0,scale1*teste.x1,scale1*teste.tn,scale1*teste.t0]
                            factor = 50
                            scale  = vmax/factor

                            plt.subplot(grid[1,0])
                            fig1 = plt.imshow(rec_ref,cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                            plt.grid()
                            plt.title('Seismogram - Reference')
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                            ax = plt.gca()
                            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                            cbar.ax.locator_params(nbins=5)
                            tick_locator = ticker.MaxNLocator(nbins=5)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                            plt.subplot(grid[1,1])
                            fig1 = plt.imshow(rec,cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                            plt.grid()
                            plt.title('Seismogram - Numerical')
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                            ax = plt.gca()
                            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                            cbar.ax.locator_params(nbins=5)
                            tick_locator = ticker.MaxNLocator(nbins=5)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                            plt.subplot(grid[1,2])
                            fig1 = plt.imshow(rec_ref-rec,cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                            plt.grid()
                            plt.title('Seismogram - Difference')
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                            ax = plt.gca()
                            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                            cbar.ax.locator_params(nbins=5)
                            tick_locator = ticker.MaxNLocator(nbins=5)
                            cbar.locator = tick_locator
                            cbar.update_ticks()

                            plt.subplot(grid[1,3])
                            plt.box(on=None)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table2 = plt.table(cellText=cell_text2,colLabels=columns2,loc='center',cellLoc='center')
                            plt.title('Seismogram - Classic Metrics',y=1.05,fontsize=15)
                            table2.set_fontsize(20)
                            table2.scale(1.0, 2.5)

                            # Seismic Trace

                            vmin  = 1.2*min(np.amin(recref),np.amin(recnum),np.amin(recnum-recref))
                            vmax  = 1.2*max(np.amax(recref),np.amax(recnum),np.amax(recnum-recref))

                            plt.subplot(grid[2,0])
                            plt.plot(vtime,recref)
                            plt.grid()
                            plt.title('Reference - Seismic Trace \n xposition = %.1f m - yposition = %.1f m'%(xpositionv[ts],ypositionv[ts]))
                            plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                            plt.ylim((vmin,vmax))
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f s'))

                            plt.subplot(grid[2,1])
                            plt.plot(vtime,recnum)
                            plt.grid()
                            plt.title('Numerical - Seismic Trace \n xposition = %.1f m - yposition = %.1f m'%(xpositionv[ts],ypositionv[ts]))
                            plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                            plt.ylim((vmin,vmax))
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f s'))
        
                            plt.subplot(grid[2,2])
                            plt.plot(vtime,recnum-recref)
                            plt.grid()
                            plt.title('Difference - Seismic Trace \n xposition = %.1f m - yposition = %.1f m'%(xpositionv[ts],ypositionv[ts]))
                            plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                            plt.ylim((vmin,vmax))
                            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f s'))
        
                            plt.subplot(grid[2,3])
                            plt.box(on=None)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table3 = plt.table(cellText=cell_text3,colLabels=columns3,loc='center',cellLoc='center')
                            plt.title('Seismic Trace - Classic Metrics \n xposition = %.1f m - yposition = %.1f m'%(xpositionv[ts],ypositionv[ts]),y=1.05,fontsize=15)
                            table3.set_fontsize(20)
                            table3.scale(1.0, 2.5)

                            plt.subplot(grid[2,4])
                            plt.box(on=None)
                            ax = plt.gca()
                            ax.get_xaxis().set_visible(False)
                            ax.get_yaxis().set_visible(False)
                            table4 = plt.table(cellText=cell_text4,colLabels=columns4,loc='center',cellLoc='center')
                            table4.set_fontsize(20)
                            plt.title('Seismic Trace - New Metrics \n xposition = %.1f m - yposition = %.1f m'%(xpositionv[ts],ypositionv[ts]),y=1.05,fontsize=15)
                            table4.scale(1.0, 2.5)

                            # FFT Results

                            plt.subplot(grid[3,0])
                            plt.semilogx(xfref, absref)
                            plt.xlabel("Frequency (Hz)")
                            plt.title('FFT - Amplitude - Reference')
                            plt.grid()

                            plt.subplot(grid[3,1])
                            plt.semilogx(xfnum, absnum)                        
                            plt.xlabel("Frequency (Hz)")
                            plt.title('FFT - Amplitude - Numerical')
                            plt.grid()

                            plt.subplot(grid[3,2])
                            plt.semilogx(xfnum, difabs)
                            plt.xlabel("Frequency (Hz)")
                            plt.title('FFT - Difference of Amplitude')
                            plt.grid()

                            plt.subplot(grid[3,3])
                            plt.semilogx(xfref, angleref)
                            plt.xlabel("Frequency (Hz)")
                            plt.title('FFT - Phase - Reference')
                            plt.grid()

                            plt.subplot(grid[3,4])
                            plt.semilogx(xfnum, anglenum)
                            plt.xlabel("Frequency (Hz)")
                            plt.title('FFT - Phase - Numerical')
                            plt.grid()

                            plt.subplot(grid[3,5])
                            plt.semilogx(xfnum, difangle)
                            plt.xlabel("Frequency (Hz)")
                            plt.title('FFT - Difference of Phases')
                            plt.grid()

                            plt.subplot(grid[3,6])
                            plt.semilogx(xfnum, ncfft)
                            plt.xlabel("Frequency (Hz)")
                            plt.title('New Curve - FFT')
                            plt.grid()

                            # CWT Results
                    
                            plt.subplot(grid[4,0])
                            plt.pcolormesh(vtime, freqs_ref, cwtmatr_ref_norm)
                            plt.yscale('log')
                            plt.xlabel("Time (s)")
                            plt.ylabel("Frequency (Hz)")
                            plt.title('CWT - Power - Reference')
                            plt.colorbar()

                            plt.subplot(grid[4,1])
                            plt.pcolormesh(vtime, freqs_num, cwtmatr_num_norm)
                            plt.yscale('log')
                            plt.xlabel("Time (s)")
                            plt.ylabel("Frequency (Hz)")
                            plt.title('CWT - Power - Numerical')
                            
                            plt.subplot(grid[4,2])
                            plt.pcolormesh(vtime, freqs_num, cwtmatr_dif_norm)
                            plt.yscale('log')
                            plt.xlabel("Time (s)")
                            plt.ylabel("Frequency (Hz)")
                            plt.title('CWT - Difference of Powers')
                            plt.colorbar()

                            plt.subplot(grid[4,3])
                            plt.pcolormesh(vtime, freqs_ref, cwtmatr_ref_angle)
                            plt.yscale('log')
                            plt.xlabel("Time (s)")
                            plt.ylabel("Frequency (Hz)")
                            plt.title('CWT - Phase - Reference')
                            plt.colorbar()

                            plt.subplot(grid[4,4])
                            plt.pcolormesh(vtime, freqs_num, cwtmatr_num_angle)
                            plt.yscale('log')
                            plt.xlabel("Time (s)")
                            plt.ylabel("Frequency (Hz)")
                            plt.title('CWT - Phase - Numerical')
                            plt.colorbar()

                            plt.subplot(grid[4,5])
                            plt.pcolormesh(vtime, freqs_num, cwtmatr_dif_angle)
                            plt.yscale('log')
                            plt.xlabel("Time (s)")
                            plt.ylabel("Frequency (Hz)")
                            plt.title('CWT - Difference of Phases')
                            plt.colorbar()

                            plt.subplot(grid[4,6])
                            plt.pcolormesh(vtime, freqs_num, nccwt)
                            plt.yscale('log')
                            plt.xlabel("Time (s)")
                            plt.ylabel("Frequency (Hz)")
                            plt.title('New Curve - CWT')
                            plt.colorbar()

                            # Hilber Results

                            plt.subplot(grid[5,0])
                            plt.plot(vtime,amplitude_envelope_ref)
                            plt.title('Hilbert - Amplitude Envelope - Reference')
                            plt.xlabel("Time (s)")
                            plt.grid()

                            plt.subplot(grid[5,1])
                            plt.plot(vtime,amplitude_envelope_num)
                            plt.title('Hilbert - Amplitude Envelope')
                            plt.xlabel("Time (s)")
                            plt.grid()

                            plt.subplot(grid[5,2])
                            plt.plot(vtime,amplitude_envelope_dif)
                            plt.title('Hilbert - Difference of Amplitude Envelope')
                            plt.xlabel("Time (s)")
                            plt.grid()

                            plt.subplot(grid[5,3])
                            plt.plot(vtime,instantaneous_phase_ref)
                            plt.xlabel("Time (s)")
                            plt.title('Hilbert Instantaneous Phase')
                            plt.ylabel('Phase in rad')
                            plt.grid()

                            plt.subplot(grid[5,4])
                            plt.plot(vtime,instantaneous_phase_num)
                            plt.title('Hilbert - Instantaneous Phase - Numerical')
                            plt.xlabel("Time (s)")
                            plt.ylabel('Phase in rad')
                            plt.grid()

                            plt.subplot(grid[5,5])
                            plt.plot(vtime,instantaneous_phase_dif)
                            plt.xlabel("Time (s)")
                            plt.ylabel('Phase in rad')
                            plt.title('Hilbert - Difference of Instantaneous Phase')
                            plt.grid()
                            
                            plt.subplot(grid[5,6])
                            plt.plot(vtime, nchilbert)
                            plt.xlabel('Time (s)')
                            plt.title('New Curve - Hilbert')
                            plt.grid()

                            plt.savefig('%s/shape_%s_method_%s_m_%d_n_%d_teste_%d_board.jpeg'%(locsave,mshape,method,mvalue,nvalue,ptype),dpi=200,bbox_inches='tight')
                            plt.close()
#==============================================================================

#==============================================================================
                if(plotnormsanalysis==1):

                    nplot  = int(len(lglobal)/8)
                    vorder = np.linspace(2,14,7)
                            
                    for b0 in range(0,6):

                        if(b0==0): normname = 'FFT RMSE - AVG'
                        elif(b0==1): normname = 'CWT RMSE - AVG'
                        elif(b0==2): normname = 'Hilbert RMSE - AVG'
                        elif(b0==3): normname = 'FFT RMSE - MAX'
                        elif(b0==4): normname = 'CWT RMSE - MAX'
                        elif(b0==5): normname = 'Hilbert RMSE - MAX'
                        
                        ci      = 0
                        cf      = 8
                        vticks  = ['s', '+', '+', '+',  '+',  '^',   '^',   'D',      'D','s']
                        vline   = ['-', '-', '-', '--', '-.', '--',  '-.',  '--',     '-.','-']
                        vcolors = ['b', 'g', 'r', 'c',  'm',  'y',   'b',   'purple', 'teal','lime']
                        linep   = [vticks,vline,vcolors]

                        fontsizevalue1="20"
                        fontsizevalue2=20
                        
                        plt.figure(figsize = (32,20)) # largura x altura
                        plt.suptitle('Signal Analysis - Norm Comp \n %s \n dx = %.3fm - dy = %.3fm \n T= %.3fs - dt = %.3fms - freq = %.3fHz'%(testname,teste.hx,teste.hy,teste.tn/1000,dt0,teste.f0),fontsize=fontsizevalue2)
                        grid  = plt.GridSpec(3,3,wspace=0.4,hspace=0.5) # largura x altura
    
                        for b1 in range(0,nplot):
                            
                            lplot        = []
                            lplot1       = []
                            lplot2       = []
                            lplotmax     = []
                            lplot1full   = []
                            lplot2full   = []
                            lplotmaxfull = []

                            for b2 in range(ci,cf-1):

                                lplot.append(lglobal[b2][4+b0])

                                lplot1full.append(lglobal[b2][16])  
                                lplot2full.append(lglobal[b2][17])  
                                lplotmaxfull.append(lglobal[b2][18])
                                    
                                l1namefull   = 'Full Receiver - L1 Norm'
                                l2namefull   = 'Full Receiver - L2 Norm'
                                maxnamefull  = 'Full Receiver - Max Norm'
                                
                                if(b0==0 or b0==1 or b0==2):

                                    lplot1.append(lglobal[b2][10])  
                                    lplot2.append(lglobal[b2][11])  
                                    lplotmax.append(lglobal[b2][12])

                                    l1name   = 'AVG L1 Norm - Seismic Traces'
                                    l2name   = 'AVG L2 Norm - Seismic Traces'
                                    maxname  = 'AVG Max Norm - Seismic Traces'

                                if(b0==3 or b0==4 or b0==5):
    
                                    lplot1.append(lglobal[b2][13])  
                                    lplot2.append(lglobal[b2][14])  
                                    lplotmax.append(lglobal[b2][15])

                                    l1name   = 'MAX L1 Norm - Seismic Traces'
                                    l2name   = 'MAX L2 Norm - Seismic Traces'
                                    maxname  = 'MAX Max Norm - Seismic Traces'

                            # print(b1,lglobal[ci][0],lglobal[ci][1])

                            if(b1!=3 and b1!=4 and b1!=7):
                                
                                plt.rc('xtick',labelsize=fontsizevalue2)
                                plt.rc('ytick',labelsize=fontsizevalue2)
                                labelpadvalue = 14

                                plt.subplot(grid[0,1])       
                                #plt.semilogy(vorder,lplot,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.semilogy(vorder,lplot,label='%s'%(lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.title('%s'%normname,fontsize=fontsizevalue2)
                                plt.grid(True)
                                plt.xlabel('[Order]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.ylabel('[Norm]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.legend(loc=1,bbox_to_anchor=(-0.2, 1.0),fontsize=fontsizevalue1)
        
                                plt.subplot(grid[1,0])       
                                #plt.semilogy(vorder,lplot1,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.semilogy(vorder,lplot1,label='%s'%(lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.title('%s'%l1name,fontsize=fontsizevalue2)
                                plt.grid(True)
                                plt.xlabel('[Order]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.ylabel('[Norm]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                #plt.legend(loc=1,bbox_to_anchor=(-0.2, 1.0),fontsize=fontsizevalue1)  

                                plt.subplot(grid[1,1])       
                                #plt.semilogy(vorder,lplot2,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.semilogy(vorder,lplot2,label='%s'%(lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.title('%s'%l2name,fontsize=fontsizevalue2)
                                plt.grid(True)
                                plt.xlabel('[Order]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.ylabel('[Norm]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                #plt.legend(loc=1,bbox_to_anchor=(-0.2, 1.0),fontsize=fontsizevalue1)          

                                plt.subplot(grid[1,2])       
                                #plt.semilogy(vorder,lplotmax,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.semilogy(vorder,lplotmax,label='%s'%(lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.title('%s'%maxname,fontsize=fontsizevalue2)
                                plt.grid(True)
                                plt.xlabel('[Order]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.ylabel('[Norm]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                #plt.legend(loc=1,bbox_to_anchor=(-0.2, 1.0),fontsize=fontsizevalue1)

                                plt.subplot(grid[2,0])       
                                #plt.semilogy(vorder,lplot1,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.semilogy(vorder,lplot1full,label='%s'%(lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.title('%s'%l1namefull,fontsize=fontsizevalue2)
                                plt.grid(True)
                                plt.xlabel('[Order]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.ylabel('[Norm]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                #plt.legend(loc=1,bbox_to_anchor=(-0.2, 1.0),fontsize=fontsizevalue1)  

                                plt.subplot(grid[2,1])       
                                #plt.semilogy(vorder,lplot2,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.semilogy(vorder,lplot2full,label='%s'%(lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.title('%s'%l2namefull,fontsize=fontsizevalue2)
                                plt.grid(True)
                                plt.xlabel('[Order]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.ylabel('[Norm]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                #plt.legend(loc=1,bbox_to_anchor=(-0.2, 1.0),fontsize=fontsizevalue1)          

                                plt.subplot(grid[2,2])       
                                #plt.semilogy(vorder,lplotmax,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.semilogy(vorder,lplotmaxfull,label='%s'%(lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                                plt.title('%s'%maxnamefull,fontsize=fontsizevalue2)
                                plt.grid(True)
                                plt.xlabel('[Order]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                plt.ylabel('[Norm]',fontsize=fontsizevalue2,labelpad=labelpadvalue)
                                #plt.legend(loc=1,bbox_to_anchor=(-0.2, 1.0),fontsize=fontsizevalue1)  

                            ci = ci + 8
                            cf = cf + 8

                        plt.savefig('%s/norm%s.png'%(locsave,b0),dpi=400,bbox_inches='tight')
                        plt.close()
#==============================================================================