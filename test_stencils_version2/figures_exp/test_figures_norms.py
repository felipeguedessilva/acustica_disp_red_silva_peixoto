#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy          as np
from   rp_files_norms import *
#==============================================================================

#==============================================================================
# Multiple Executations - Comum Tests
#==============================================================================
# vptype    = [1,2,3,4] 
# vdxref    = [1]
# vdtref    = [1]
# vfreqref  = [1,2,3] 

vptype    = [1] 
vdxref    = [2]
vdtref    = [1]
vfreqref  = [1] 

nvptype      = len(vptype) 
nvdxref      = len(vdxref)
nvdtref      = len(vdtref)
nvfreqref    = len(vfreqref) 
cont_me      = 0
lconfig_me   = []

for k0 in range(0,nvptype):
    
    for k1 in range(0,nvdxref):
        
        for k2 in range(0,nvdtref):
            
            for k3 in range(0,nvfreqref):

                ptype    = vptype[k0] 
                dx_ref   = vdxref[k1]
                dt_ref   = vdtref[k2]
                freq_ref = vfreqref[k3]

                print('')
                print('Number    = %d'%cont_me)
                print('Test Type = %d'%ptype)
                print('dx_ref    = %d'%dx_ref)
                print('dt_ref    = %d'%dt_ref)
                print('freq_ref  = %d'%freq_ref)
                print('')
                
                gera_norms(ptype,dx_ref,dt_ref,freq_ref)
                
                lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
                
                cont_me = cont_me + 1
#==============================================================================