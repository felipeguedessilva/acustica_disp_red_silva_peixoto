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
vptype    = [1,2,3,4] 
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
                locopen = '../testresults/test%d_results_norms_fch_%d_%d_%d'%(ptype,dx_ref,dt_ref,freq_ref)
                with open(locopen, 'rb') as f: 
                    test_results  = pickle.load(f) 

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
# Obtenção de Parâmetros para Open
#==============================================================================
                lglobal = []
                ntr     = len(test_results)

                for l0 in range(0,ntr):

                    nvalueloc = test_results[l0][3]

                    if(nvalueloc==1): lglobal.append(test_results[l0])

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
                    plt.suptitle('Signal Analysis - Norm Comp \n %s \n dx = %.3fm - dy = %.3fm \n T= %.3fs - dt = %.3fms - freq = %.3fHz'%(testname,teste.hx,teste.hy,teste.tn/1000,0.5*dt_ref,teste.f0),fontsize=fontsizevalue2)
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

                            ncom      = len(lplot)
                            factorloc = 10000
                            min1      = np.amin(np.array(lplot))
                            min2      = np.amin(np.array(lplot1))
                            min3      = np.amin(np.array(lplot2))
                            min4      = np.amin(np.array(lplotmax))
                            min5      = np.amin(np.array(lplot1full))
                            min6      = np.amin(np.array(lplot2full))
                            min7      = np.amin(np.array(lplotmaxfull))

                            for c1 in range(0,ncom):

                                # if(lplot[c1]<factorloc*min1 or lplot[c1]==np.inf):        lplot[c1] = np.nan
                                # if(lplot1[c1]<factorloc*min2 or lplot1[c1]==np.inf):       lplot1[c1] = np.nan
                                # if(lplot2[c1]<factorloc*min3 or lplot2[c1]==np.inf):       lplot2[c1] = np.nan
                                # if(lplotmax[c1]<factorloc*min4 or lplotmax[c1]==np.inf):     lplotmax[c1] = np.nan
                                # if(lplot1full[c1]<factorloc*min5 or lplot1full[c1]==np.inf):   lplot1full[c1] = np.nan
                                # if(lplot2full[c1]<factorloc*min6 or lplot2full[c1]==np.inf):   lplot2full[c1] = np.nan
                                # if(lplotmaxfull[c1]<factorloc*min7 or lplotmaxfull[c1]==np.inf): lplotmaxfull[c1] = np.nan

                                if(lplot[c1]>factorloc or lplot[c1]==np.inf):        lplot[c1] = np.nan
                                if(lplot1[c1]>factorloc or lplot1[c1]==np.inf):       lplot1[c1] = np.nan
                                if(lplot2[c1]>factorloc or lplot2[c1]==np.inf):       lplot2[c1] = np.nan
                                if(lplotmax[c1]>factorloc or lplotmax[c1]==np.inf):     lplotmax[c1] = np.nan
                                if(lplot1full[c1]>factorloc or lplot1full[c1]==np.inf):   lplot1full[c1] = np.nan
                                if(lplot2full[c1]>factorloc or lplot2full[c1]==np.inf):   lplot2full[c1] = np.nan
                                if(lplotmaxfull[c1]>factorloc or lplotmaxfull[c1]==np.inf): lplotmaxfull[c1] = np.nan

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
                    print('Finish!')
#==============================================================================