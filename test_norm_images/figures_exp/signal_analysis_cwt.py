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
from numpy.fft import rfft, rfftfreq, fft, fftfreq
#==============================================================================

#==============================================================================
# Range of Parameters
#==============================================================================
vptype    = [1] 
vdxref    = [1]
vdtref    = [1]
vfreqref  = [1] 

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
                    xpositionv = np.array([750.0,2250.0, 750.0,2250.0])
                    ypositionv = np.array([750.0, 750.0,2250.0,2250.0])
                    tf_test    = 900

                if(ptype==2): 
                    
                    teste_ref  = ttopt.teste2_ref1(freq_ref,factor_ref)
                    teste      = ttopt.teste2(dx_ref,dt_ref,freq_ref)
                    testname   = 'Heterogeneos Velocity Model'
                    xpositionv = np.array([500.0,1500.0, 500.0,1500.0])
                    ypositionv = np.array([500.0, 500.0,1500.0,1500.0])
                    tf_test    = 3000

                if(ptype==3): 
                
                    teste_ref  = ttopt.teste3_ref1(freq_ref,factor_ref)
                    teste      = ttopt.teste3(dx_ref,dt_ref,freq_ref)
                    testname   = 'SEG/EAGE 2D Salt Velocity Model'
                    xpositionv = np.array([4000.0,4000.0,4000.0,6000.0,6000.0,6000.0,8000.0,8000.0,8000.0])   
                    ypositionv = np.array([2000.0,2500.0,3000.0,2000.0,2500.0,3000.0,2000.0,2500.0,3000.0]) 
                    tf_test    = 3000

                if(ptype==4): 
                    
                    teste_ref   = ttopt.teste4_ref1(freq_ref,factor_ref)
                    teste       = ttopt.teste4(dx_ref,dt_ref,freq_ref)
                    testname    = 'Marmousi Velocity Model'
                    xpositionv  = np.array([6000.0,6000.0,6000.0,8000.0,8000.0,8000.0,10000.0,10000.0,10000.0,12000.0,12000.0,12000.0])
                    ypositionv  = np.array([1000.0,2000.0,3000.0,1000.0,2000.0,3000.0,1000.0,2000.0,3000.0,1000.0,2000.0,3000.0])
                    tf_test     = 3000

                locopen    = 'signals/signal_files/teste%d/'%(ptype)
                locsave    = 'signals/figs/teste%d/'%(ptype)
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
                lglobal = []
    
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
                    
                    rec_select     = np.load("%srec_select_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                    rec_select_ref = np.load("%srec_select_ref.npy"%(locopen))
                    
                    pi = 0
                    pf = rec_select_ref.shape[0]  

                    # pi = 650
                    # pf = 900  

                    dt     = (tf_test/(rec_select_ref.shape[0]-1))/1000
                    recref = rec_select_ref[pi:pf,0]
                    recnum = rec_select[pi:pf,0]

                    vmin  = 1.2*min(np.amin(recref),np.amin(recnum))
                    vmax  = 1.2*max(np.amax(recref),np.amax(recnum))
                    vtime = np.linspace(0,tf_test,recref.shape[0])/1000

                    vmindif = 1.2*np.amin(recref-recnum)
                    vmaxdif = 1.2*np.amax(recref-recnum)

                    sampling_period = np.diff(vtime).mean()
                    widths = np.geomspace(1, 1024, num=100)
                    wavelet = "cmor1.5-1.0"
                    normtype = 1
     
                    cwtmatr_ref, freqs_ref = pywt.cwt(recref, widths, wavelet, sampling_period=sampling_period,method='fft')
                    cwtmatr_num, freqs_num = pywt.cwt(recnum, widths, wavelet, sampling_period=sampling_period,method='fft')

                    cwtmatr_ref_norm = np.abs(cwtmatr_ref)
                    cwtmatr_num_norm = np.abs(cwtmatr_num)
                    cwtmatr_dif_norm = np.abs(cwtmatr_ref_norm - cwtmatr_num_norm) 

                    cwtmatr_ref_angle = np.abs(np.angle(cwtmatr_ref))
                    cwtmatr_num_angle = np.abs(np.angle(cwtmatr_num))
                    cwtmatr_dif_angle = np.abs(cwtmatr_ref_angle - cwtmatr_num_angle)

                    nf_abs_ref = la.norm(cwtmatr_ref_norm,normtype)
                    nf_abs_num = la.norm(cwtmatr_num_norm,normtype)
                    nf_abs_dif = la.norm(cwtmatr_dif_norm,normtype)/nf_abs_ref

                    nf_ang_ref = la.norm(cwtmatr_ref_angle,normtype)
                    nf_ang_num = la.norm(cwtmatr_num_angle,normtype)
                    nf_ang_dif = la.norm(cwtmatr_dif_angle,normtype)/nf_ang_ref

                    maxval   = max(np.amax(cwtmatr_num),np.amax(cwtmatr_ref))
                    minval   = min(np.amin(cwtmatr_num),np.amin(cwtmatr_ref))
                    valratio = 10**(-4)

                    cwtmatr_num_mod = cwtmatr_num
                    cwtmatr_ref_mod = cwtmatr_ref

                    cwtmatr_ref_mod[cwtmatr_ref_mod<valratio*maxval] = np.nan
                    cwtmatr_num_mod[cwtmatr_num_mod<valratio*maxval] = np.nan
                    
                    cwtmatr_ratio = cwtmatr_num_mod/cwtmatr_ref_mod

                    cwtmatr_ratio_norm  = np.abs(cwtmatr_ratio)
                    cwtmatr_ratio_angle = np.abs(np.angle(cwtmatr_ratio))
   
                    cwtmatr_ratio_norm_mod = np.nan_to_num(cwtmatr_ratio)
                    cwtmatr_ratio_angle_mod = np.nan_to_num(cwtmatr_ratio_angle)

                    nf_attuation  = la.norm(cwtmatr_ratio_norm_mod,normtype)
                    nf_wavenumber = la.norm(cwtmatr_ratio_angle_mod,normtype)
                    

                    if(nvalue==1):

                        llocal = []
                        llocal.append(mshape)
                        llocal.append(method)
                        llocal.append(mvalue)
                        llocal.append(nvalue)
                        llocal.append(nf_abs_num)
                        llocal.append(nf_abs_dif)
                        llocal.append(nf_ang_num)
                        llocal.append(nf_ang_dif)
                        llocal.append(nf_attuation)
                        llocal.append(nf_wavenumber)
                        lglobal.append(llocal)

                    plotfig = 0

                    if(plotfig==1):

                        plt.figure(figsize = (24,14))
                        plt.suptitle('Signal Analysis \n %s \n xposition = %.1f m - yposition = %.1f m'%(testname,xpositionv[0],ypositionv[0]))
                        grid  = plt.GridSpec(4,4,wspace=0.15,hspace=0.5)  

                        plt.subplot(grid[0,0])
                        plt.plot(vtime,recref)
                        plt.grid()
                        plt.title('Reference')
                        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                        plt.ylim((vmin,vmax))
                        plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                        plt.ylabel('[Error]')

                        plt.subplot(grid[0,1])
                        plt.plot(vtime,recnum)
                        plt.grid()
                        plt.title('Numerical \n Shape = %s - Method = %s - M = %d - N = %d'%(mshape,method,mvalue,nvalue))
                        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                        plt.ylim((vmin,vmax))
                        plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                        plt.ylabel('[Error]')

                        plt.subplot(grid[0,2])
                        plt.plot(vtime,recref,label='Reference')
                        plt.plot(vtime,recnum,label='Numerical')
                        plt.grid()
                        plt.title('Comparision Reference and Numerical \n Shape = %s - Method = %s - M = %d - N = %d'%(mshape,method,mvalue,nvalue))
                        plt.legend()
                        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                        plt.ylim(vmin,vmax)
                        plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                        plt.ylabel('[Error]')
                
                        plt.subplot(grid[0,3])
                        plt.plot(vtime,recref-recnum)
                        plt.grid()
                        plt.title('Difference Reference and Numerical \n Shape = %s - Method = %s - M = %d - N = %d'%(mshape,method,mvalue,nvalue))
                        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                        plt.ylim(vmindif,vmaxdif)
                        plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                        plt.ylabel('[Error]')

                        plt.subplot(grid[1,0])
                        plt.pcolormesh(vtime, freqs_ref, cwtmatr_ref_norm[:,:])
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Norm - Reference - NF = %.3f"%nf_abs_ref)
                        plt.colorbar()

                        plt.subplot(grid[1,1])
                        plt.pcolormesh(vtime, freqs_num, cwtmatr_num_norm)
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Norm - Numerical - NF = %.3f"%nf_abs_num)
                        plt.colorbar()

                        plt.subplot(grid[1,2])
                        plt.pcolormesh(vtime, freqs_num, cwtmatr_dif_norm)
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Norm - Difference - NF = %.3f"%nf_abs_dif)
                        plt.colorbar()

                        plt.subplot(grid[2,0])
                        plt.pcolormesh(vtime, freqs_ref, cwtmatr_ref_angle)
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Angle - Reference - NF = %.3f"%nf_ang_ref)
                        plt.colorbar()

                        plt.subplot(grid[2,1])
                        plt.pcolormesh(vtime, freqs_num, cwtmatr_num_angle)
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Angle - Numerical- NF = %.3f"%nf_ang_num)
                        plt.colorbar()

                        plt.subplot(grid[2,2])
                        plt.pcolormesh(vtime, freqs_num, cwtmatr_dif_angle)
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Angle - Difference- NF = %.3f"%nf_ang_dif)
                        plt.colorbar()

                        plt.subplot(grid[3,0])
                        plt.pcolormesh(vtime, freqs_num, cwtmatr_ratio_norm)
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Attenuation A(f) \n NF = %.3f"%nf_attuation)
                        plt.colorbar()
                        
                        plt.subplot(grid[3,1])
                        plt.pcolormesh(vtime, freqs_num, cwtmatr_ratio_angle)
                        plt.yscale('log')
                        plt.xlabel("Time (s)")
                        plt.ylabel("Frequency (Hz)")
                        plt.title("Wavenumber Relation: k(f) -fk'(f) \n NF = %.3f"%nf_wavenumber)
                        plt.colorbar()
                        plt.savefig('%s/sa_shape_%s_method_%s_m_%d_n_%d_teste_%d_cwt.jpeg'%(locsave,mshape,method,mvalue,nvalue,ptype),dpi=200,bbox_inches='tight')
                        plt.close()
#==============================================================================

#==============================================================================
# CWT Group Norm Analysis
#==============================================================================

#==============================================================================
                plotnormsanalysis = 0

                if(plotnormsanalysis==1):

                    nplot  = int(len(lglobal)/8)
                    vorder = np.linspace(2,16,8)
        
                    for b0 in range(0,6):

                        if(b0==0): normname = 'CWT - Norm Abs   - Numerical'
                        elif(b0==1): normname = 'CWT - Norm Abs   - Difference'
                        elif(b0==2): normname = 'CWT - Norm Angle - Numerical'
                        elif(b0==3): normname = 'CWT - Norm Angle - Difference'
                        elif(b0==4): normname = 'CWT - Norm Attenuation'
                        elif(b0==5): normname = 'CWT - Norm Wavenumber Relation'

                        ci      = 0
                        cf      = 8
                        vticks  = ['s', '+', '+', '+',  '+',  '^',   '^',   'D',      'D','s']
                        vline   = ['-', '-', '-', '--', '-.', '--',  '-.',  '--',     '-.','-']
                        vcolors = ['b', 'g', 'r', 'c',  'm',  'y',   'b',   'purple', 'teal','lime']
                        linep   = [vticks,vline,vcolors]

                        for b1 in range(0,nplot):
                            
                            lplot = []

                            for b2 in range(ci,cf):
                            
                                lplot.append(lglobal[b2][4+b0])  

                            plt.plot(vorder,lplot,label='%s-%s'%(lglobal[ci][0],lglobal[ci][1]),color=linep[2][b1],linestyle=linep[1][b1],marker=linep[0][b1])
                            plt.title('%s'%normname)
                            plt.grid(True)
                            plt.xlabel('[Order]')
                            plt.ylabel('[Norm]')
                            plt.legend(loc=1,bbox_to_anchor=(-0.15, 1.0))
                            ci = ci + 8
                            cf = cf + 8
                        
                        plt.savefig('%s/norm%s_cwt.png'%(locsave,b0),dpi=400,bbox_inches='tight')
                        plt.close()
#==============================================================================