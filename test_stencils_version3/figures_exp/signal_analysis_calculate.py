#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
from   numpy import linalg      as la
import sys
import pickle
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
vptype    = [4] 
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
    
    testresults = []
    
    for k1 in range(0,nvdxref):
        
        for k2 in range(0,nvdtref):
            
            for k3 in range(0,nvfreqref):

                ptype      = vptype[k0] 
                dx_ref     = vdxref[k1]
                dt_ref     = vdtref[k2]
                freq_ref   = vfreqref[k3]
                factor_ref = 16
                lglobal    = []

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
                                
                            else:
                                
                                vnvalue  = np.arange(1,mvalue+1)
                                
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
# OOpen Reference
#============================================================================== 
                locopenref = '../data_save/teste%d/reffreq%d/'%(ptype,freq_ref)
                sol_ref        = np.load("%ssolplotcut_%d_%d_%d_%d.npy"%(locopenref,ptype,dx_ref,dt_ref,freq_ref))
                rec_ref        = np.load("%srecrefcut_%d_%d_%d_%d.npy"%(locopenref,ptype,dx_ref,dt_ref,freq_ref))
                rec_select_ref = np.load("%srecselectcut_%d_%d_%d_%d.npy"%(locopenref,ptype,dx_ref,dt_ref,freq_ref))
                   
                if(dt_ref==6 and ptype==1):

                    rec_ref        = rec_ref[1:,:]
                    rec_select_ref = rec_select_ref[1:,:]
#==============================================================================

#==============================================================================
# Obtenção de Parâmetros para Open
#==============================================================================                
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
                    
                    try:
                        
                        rec_selectnorm1 = la.norm(rec_select-rec_select_ref,1)/la.norm(rec_select_ref,1)

                    except:

                        rec_selectnorm1 = np.nan

                    try:
                        
                        rec_selectnorm2 = la.norm(rec_select-rec_select_ref,2)/la.norm(rec_select_ref,2)
                    
                    except:

                        rec_selectnorm2 = np.nan 

                    try:
                        
                        rec_selectnormmax = la.norm(rec_select-rec_select_ref,np.inf)/la.norm(rec_select_ref,np.inf)

                    except:

                        rec_selectnormmax = np.nan

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
            
                    testresults.append([cont_me,ptype,dx_ref,dt_ref,freq_ref,mshape,method,mvalue,nvalue,npt,npe,llocal,parameters,cont_glob])

                    lglobal.append(llocal)
                    cont_glob = cont_glob + 1
#==============================================================================

#==============================================================================
# Save Results1
#==============================================================================
                locname = '../testresults/test%d_results_norms_fch_%d_%d_%d'%(ptype,dx_ref,dt_ref,freq_ref)
                with open(locname, 'wb') as f: 
                    pickle.dump(lglobal, f) 
#==============================================================================

#==============================================================================
# Save Results2
#==============================================================================
    locname = '../testresults/test%d_results_norms_fch'%(ptype)
    with open(locname, 'wb') as f: 
        pickle.dump(testresults, f) 
#==============================================================================