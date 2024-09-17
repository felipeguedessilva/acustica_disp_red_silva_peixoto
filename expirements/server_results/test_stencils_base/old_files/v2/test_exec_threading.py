#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy                    as np
import threading
from   devito_acustica_one_test import *
from   devito_acustica_gera     import *
from   devito_acustica_ref      import *
#==============================================================================

#==============================================================================
# Multiple Executations - Comum Tests
#==============================================================================
if __name__ == "__main__": 

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
                    
                    print("Starting Job!") 
                    p1 = threading.Thread(target=acoustic_operator_multiple_exec, args=(dx_ref,dt_ref,freq_ref,ptype, ))                  
                    p1.start()
                    p1.join() 
                    print("Job Done!") 
                    
                    lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
                    
                    cont_me = cont_me + 1
#==============================================================================