#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy                   as np
import sys
import matplotlib.pyplot       as plot
import time                    as tm
#==============================================================================

#==============================================================================
# Shell Comands
#==============================================================================
import os
os.system('clear')
#==============================================================================

#==============================================================================
# Print Configurations
#==============================================================================
np.set_printoptions(formatter={'float': '{: 0.3e}'.format})
#==============================================================================

#==============================================================================
# Basic Coefs
#==============================================================================
import ger_coef1 as coefdispred1
import ger_coef2 as coefdispred2
import ger_coef3 as coefdispred3
#==============================================================================

#==============================================================================
plot.close("all")
#==============================================================================

#============================================================================================================
vdx      = [30,20,40,40]
vdt      = [0.5,0.5,0.5,0.5]
vvmax    = [3.0,3.0,5.0,5.0]
nparam   = len(vdx)
vmethod  = ['spatte','spectetheta','dispte','specls','displs'] 
nvmethod = len(vmethod)
vmvalue  = [1,2,3,4,5,6,7,8,9,10]
nvmvalue = len(vmvalue)

for k0 in range(0,nparam):

    ntestes  = 0
    dx       = vdx[k0]
    dt       = vdt[k0]
    vmax     = vvmax[k0]
    cur      = (vmax*dt)/(dx)
    pteste   = 'Teste {}'.format(k0+1)

    for k1 in range(0,nvmethod):
    
        method = vmethod[k1]
    
        if(method == 'spatte'):
        
            vshape = ['cl']
        
        elif(method == 'spectetheta'):
        
            vshape = ['cl']
    
        elif(method == 'dispte' or method == 'specls' or method == 'displs'):
        
            vshape = ['crb','csq']

        nvshape = len(vshape)
    
        for k2 in range(0,nvshape):
    
            shape = vshape[k2]

            for k3 in range(0,nvmvalue):
            
                mvalue = vmvalue[k3]
            
                if(shape == 'cl'):
                
                    vnvalue = [1]
            
                elif(shape == 'rb'):
                
                    vnvalue = [mvalue]
            
                else:
            
                    vnvalue  = np.arange(1,mvalue+1)
            
                nvnvalue = len(vnvalue)
            
                for k4 in range(0,nvnvalue):
                
                    nvalue = vnvalue[k4]
                    T0     = coefdispred3.calccoef(method,shape,mvalue,nvalue,cur)

                    print('')
                    print('=========================================================================')
                    print('%s - Config - %d'%(pteste,ntestes))
                    print('Teste : %d'%ntestes)
                    print('Shape : %s'%shape)
                    print('Method: %s'%method)
                    print('mvalue: %d'%mvalue)
                    print('nvalue: %d'%nvalue)
                    print('dx    : %.2e m'%dx)
                    print('dt    : %.2e ms'%dt)
                    print('vmax  : %.2e KM/s'%vmax)
                    print('cur   : %.2e'%cur)
                    print('Status: OK!')
                    print('=========================================================================')
                    print('')
                  
                    ntestes = ntestes + 1
#============================================================================================================

#============================================================================================================                
print('O total de testes é: %d'%ntestes)
#============================================================================================================