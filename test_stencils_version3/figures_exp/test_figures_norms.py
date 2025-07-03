#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy          as np
from   rp_files_norms import *
#==============================================================================

#==============================================================================
# Multiple Executations - Comum Tests
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

for k0 in range(0,nvptype):
    
    for k1 in range(0,nvdxref):
        
        for k2 in range(0,nvdtref):
            
            for k3 in range(0,nvfreqref):

                ptype    = vptype[k0] 
                dx_ref   = vdxref[k1]
                dt_ref   = vdtref[k2]
                freq_ref = vfreqref[k3]

                if(ptype==1 and dx_ref==2 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==1 and dx_ref==4 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==1 and dx_ref==4 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==1 and dx_ref==8 and dt_ref==2):

                    print('Jump Test!')

                elif(ptype==1 and dx_ref==8 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==1 and dx_ref==8 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==2 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==2 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==4 and dt_ref==2):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==4 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==4 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==8 and dt_ref==1):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==8 and dt_ref==2):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==8 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==2 and dx_ref==8 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==3 and dx_ref==2 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==3 and dx_ref==4 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==3 and dx_ref==4 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==3 and dx_ref==8 and dt_ref==2):

                    print('Jump Test!')

                elif(ptype==3 and dx_ref==8 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==3 and dx_ref==8 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==4 and dx_ref==2 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==4 and dx_ref==4 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==4 and dx_ref==4 and dt_ref==6):

                    print('Jump Test!')

                elif(ptype==4 and dx_ref==8 and dt_ref==2):

                    print('Jump Test!')

                elif(ptype==4 and dx_ref==8 and dt_ref==4):

                    print('Jump Test!')

                elif(ptype==4 and dx_ref==8 and dt_ref==6):

                    print('Jump Test!')

                else:

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

# #==============================================================================
# # Multiple Executations - Comum Tests
# #==============================================================================
# vptype    = [1] 
# vdxref    = [1,2,4,8]
# vdtref    = [1,2,4,6]
# vfreqref  = [1,2,3] 

# nvptype      = len(vptype) 
# nvdxref      = len(vdxref)
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
                
#                 gera_norms(ptype,dx_ref,dt_ref,freq_ref)
                
#                 lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
                
#                 cont_me = cont_me + 1
# #==============================================================================