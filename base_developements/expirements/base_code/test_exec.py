#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy                    as np
from   devito_acustica_one_test import *
from   devito_acustica_gera     import *
from   devito_acustica_ref      import *
#==============================================================================

#==============================================================================
# One Executation - Comum Test
#==============================================================================
# vptype    = [1,2,3,4] 
# vdxref    = [1,2,4,8]
# vdtref    = [1,2,4,6]
# vfreqref  = [1,2,3] 

# nvptype      = len(vdxref) 
# nvdxref      = len(vptype)
# nvdtref      = len(vdtref)
# nvfreqref    = len(vfreqref) 
# cont_ct      = 0
# lconfig_ct   = []

# for k0 in range(0,nvptype):
    
#     for k1 in range(0,nvdxref):
        
#         for k2 in range(0,nvdtref):
            
#             for k3 in range(0,nvfreqref):

#                 ptype    = vptype[k0] 
#                 dx_ref   = vdxref[k1]
#                 dt_ref   = vdtref[k2]
#                 freq_ref = vfreqref[k3]

#                 print('')
#                 print('Number    = %d'%cont_ct)
#                 print('Test Type = %d'%ptype)
#                 print('dx_ref    = %d'%dx_ref)
#                 print('dt_ref    = %d'%dt_ref)
#                 print('freq_ref  = %d'%freq_ref)
#                 print('')
                
#                 acoustic_operator_one_exec(dx_ref,dt_ref,freq_ref,ptype)
                
#                 lconfig_ct.append((cont_ct,ptype,dx_ref,dt_ref,freq_ref))
                
#                 cont_ct = cont_ct + 1
#==============================================================================

#==============================================================================
# Multiple Executation - Reference Tests
#==============================================================================
# vptype     = [1,2,3,4] 
# vfreqref   = [1,2,3] 
# nvptype    = len(vptype) 
# nvfreqref  = len(vfreqref) 
# cont_rt    = 0
# lconfig_rt = []

# for k0 in range(0,nvptype):

#     for k1 in range(0,nvfreqref):

#         ptype      = vptype[k0] 
#         freq_ref   = vfreqref[k1]
#         factor_ref = 10

#         print('')
#         print('Number    = %d'%cont_rt)
#         print('Test Type = %d'%ptype)
#         print('freq_ref  = %d'%freq_ref)
#         print('')
                
#         acoustic_operator_one_exec_ref(freq_ref,factor_ref,ptype)
                
#         lconfig_rt.append((cont_rt,ptype,freq_ref))
                
#         cont_rt = cont_rt + 1
#==============================================================================

#==============================================================================
# Multiple Executations - Comum Tests
#==============================================================================
# vptype    = [1,2,3,4] 
# vdxref    = [1,2,4,8]
# vdtref    = [1,2,4,6]
# vfreqref  = [1,2,3] 

# nvptype      = len(vdxref) 
# nvdxref      = len(vptype)
# nvdtref      = len(vdtref)
# nvfreqref    = len(vfreqref) 
# cont_me      = 0
# lconfig_me   = []

# for k0 in range(0,nvptype):
    
#     for k1 in range(0,nvdxref):
        
#         for k2 in range(0,nvdtref):
            
#             for k3 in range(0,nvfreqref):

#                 ptype    = vptype[k0] 
#                 dx_ref   = vdxref[k1]
#                 dt_ref   = vdtref[k2]
#                 freq_ref = vfreqref[k3]

#                 print('')
#                 print('Number    = %d'%cont_me)
#                 print('Test Type = %d'%ptype)
#                 print('dx_ref    = %d'%dx_ref)
#                 print('dt_ref    = %d'%dt_ref)
#                 print('freq_ref  = %d'%freq_ref)
#                 print('')
                
#                 acoustic_operator_multiple_exec(dx_ref,dt_ref,freq_ref,ptype)
                
#                 lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
                
#                 cont_me = cont_me + 1
#==============================================================================