#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
from   numpy import linalg      as la
import sys
import pickle
import gc
#==============================================================================

#==============================================================================
# My Modules
#==============================================================================
import testes_opt              as ttopt
#==============================================================================

#==============================================================================
# Signal Comparison
#==============================================================================
from im_fft_norms import fftnorm1
from im_fft_norms import fftnorm2
from im_fft_norms import fftnorm3
#==============================================================================

#==============================================================================
# Plot Set
#==============================================================================
import matplotlib.pyplot       as plt
import matplotlib.ticker       as mticker    
#from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib              import ticker
from   matplotlib              import cm
#==============================================================================

#==============================================================================
# Range of Parameters
#==============================================================================
save_fields   = 0 

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
    
    testresults = []

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
                    
                    teste_ref = ttopt.teste1_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste1(dx_ref,dt_ref,freq_ref)
                    
                if(ptype==2): 
                    
                    teste_ref = ttopt.teste2_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste2(dx_ref,dt_ref,freq_ref)
                
                if(ptype==3): 
                
                    teste_ref = ttopt.teste3_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste3(dx_ref,dt_ref,freq_ref)
                
                if(ptype==4): 
                    
                    teste_ref = ttopt.teste4_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste4(dx_ref,dt_ref,freq_ref)
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
# Open and Close Referencia
#==============================================================================    
                print('')
                print('Open Ref Files!')
                locopenref     = 'teste%d/reffreq%d'%(ptype,freq_ref)
                rec_ref        = np.load("../data_save/%s/rec_ref.npy"%(locopenref))   
                solplot_ref    = np.load("../data_save/%s/solplot_ref.npy"%(locopenref))            
                rec_select_ref = np.load("../data_save/%s/rec_select_ref.npy"%(locopenref))
                print('')
                print('Close Ref Files!')
                print('')
#==============================================================================

#==============================================================================
# Open and Close Tests
#==============================================================================    
                locopen = 'teste%d/dx%ddt%dfreq%d'%(ptype,dx_ref,dt_ref,freq_ref)
                print('Open Test Files!')
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
                    
                rec        = np.load("../data_save/%s/rec_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))    
                solplot    = np.load("../data_save/%s/solplot_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                rec_select = np.load("../data_save/%s/rec_select_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                    
                a              = rec_ref.shape[0]-1
                b              = rec.shape[0]-1
                c              = rec_ref.shape[1]-1
                d              = rec.shape[1]-1
                irecref0       = int(a/b)
                irecref1       = int(c/d)
                                        
                a              = solplot_ref.shape[0]-1
                b              = solplot.shape[0]-1
                c              = solplot_ref.shape[1]-1
                d              = solplot.shape[1]-1
                e              = solplot_ref.shape[2]-1
                f              = solplot.shape[2]-1
                isolplotref0   = int(a/b)
                isolplotref1   = int(c/d)
                isolplotref2   = int(e/f)
                                        
                a              = rec_select_ref.shape[0]-1
                b              = rec_select.shape[0]-1
                c              = rec_select_ref.shape[1]-1
                d              = rec_select.shape[1]-1
                irecselectref0 = int(a/b)
                irecselectref1 = int(c/d)
                  
                recrefcut    = rec_ref[0::irecref0,0::irecref1]
                solplotcut   = solplot_ref[0::isolplotref0,0::isolplotref1,0::isolplotref2]
                recselectcut = rec_select_ref[0::irecselectref0,0::irecselectref1]

                locosaverefcut = '../data_save/teste%d/reffreq%d/'%(ptype,freq_ref)
                np.save("%ssolplotcut_%d_%d_%d_%d.npy"%(locosaverefcut,vptype[k0],vdxref[k1],vdtref[k2],vfreqref[k3]),solplotcut[:,:,:])
                np.save("%srecrefcut_%d_%d_%d_%d.npy"%(locosaverefcut,vptype[k0],vdxref[k1],vdtref[k2],vfreqref[k3]),recrefcut)
                np.save("%srecselectcut_%d_%d_%d_%d.npy"%(locosaverefcut,vptype[k0],vdxref[k1],vdtref[k2],vfreqref[k3]),recselectcut)
                
                del rec_ref, solplot_ref, rec_select_ref
                gc.collect()
#==============================================================================