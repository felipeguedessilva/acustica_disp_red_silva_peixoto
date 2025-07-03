#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
#import numpy                    as np
from   devito_acustica_gera     import *
from   devito_acustica_ref      import *
#==============================================================================

#==============================================================================
execsol = 1
execref = 0
#==============================================================================

#==============================================================================
# Multiple Executations - Comum Tests
#==============================================================================
if(execsol==1):

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
    totaltests   = nvptype*nvdxref*nvdtref*nvfreqref

    for k0 in range(0,nvptype):
        
        for k1 in range(0,nvdxref):
            
            for k2 in range(0,nvdtref):
                
                for k3 in range(0,nvfreqref):

                    ptype        = vptype[k0] 
                    dx_ref       = vdxref[k1]
                    dt_ref       = vdtref[k2]
                    freq_ref     = vfreqref[k3]
                    percent_glob = (100*(cont_me+1))/(totaltests)
            
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
                    
                    acoustic_operator_multiple_exec(dx_ref,dt_ref,freq_ref,ptype,percent_glob)
                    lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
                    cont_me = cont_me + 1
#==============================================================================

# #==============================================================================
# # Multiple Executations - Comum Tests
# #==============================================================================
# if(execsol==1):

#     vptype    = [1] 
#     vdxref    = [1,2,4,8]
#     vdtref    = [1,2,4,6]
#     vfreqref  = [1,2,3] 

#     nvptype      = len(vptype) 
#     nvdxref      = len(vdxref)
#     nvdtref      = len(vdtref)
#     nvfreqref    = len(vfreqref) 
#     cont_me      = 0
#     lconfig_me   = []
#     totaltests   = nvptype*nvdxref*nvdtref*nvfreqref

#     for k0 in range(0,nvptype):
        
#         for k1 in range(0,nvdxref):
            
#             for k2 in range(0,nvdtref):
                
#                 for k3 in range(0,nvfreqref):

#                     ptype        = vptype[k0] 
#                     dx_ref       = vdxref[k1]
#                     dt_ref       = vdtref[k2]
#                     freq_ref     = vfreqref[k3]
#                     percent_glob = (100*(cont_me+1))/(totaltests)
            
#                     print('')
#                     print('Number    = %d'%cont_me)
#                     print('Test Type = %d'%ptype)
#                     print('dx_ref    = %d'%dx_ref)
#                     print('dt_ref    = %d'%dt_ref)
#                     print('freq_ref  = %d'%freq_ref)
#                     print('')
                    
#                     acoustic_operator_multiple_exec(dx_ref,dt_ref,freq_ref,ptype,percent_glob)
#                     lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
#                     cont_me = cont_me + 1
# #==============================================================================

#==============================================================================
# Multiple Executation - Reference Tests
#==============================================================================
if(execref==1):

    vptype     = [4] 
    vfreqref   = [1,2,3]
    
    nvptype    = len(vptype) 
    nvfreqref  = len(vfreqref) 
    cont_rt    = 0
    lconfig_rt = []

    for k0 in range(0,nvptype):

        for k1 in range(0,nvfreqref):

            ptype      = vptype[k0] 
            freq_ref   = vfreqref[k1]
            factor_ref = 16

            print('')
            print('Number    = %d'%cont_rt)
            print('Test Type = %d'%ptype)
            print('freq_ref  = %d'%freq_ref)
            print('')
                    
            acoustic_operator_one_exec_ref(freq_ref,factor_ref,ptype)
                    
            lconfig_rt.append((cont_rt,ptype,freq_ref))
                    
            cont_rt = cont_rt + 1
#==============================================================================