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
# Range of Parameters
#==============================================================================
vptype        = [1] 
vdxref        = [1,2,4,8]
vdtref        = [1,2,4,6]
vfreqref      = [1,2,3] 

# vptype        = [1] 
# vdxref        = [1]
# vdtref        = [1]
# vfreqref      = [1] 

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
                            
                        vshape = ['crb']
                            
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
                rec_ref        = np.load("data_save/%s/rec_ref.npy"%(locopenref))   
                solplot_ref    = np.load("data_save/%s/solplot_ref.npy"%(locopenref))            
                rec_select_ref = np.load("data_save/%s/rec_select_ref.npy"%(locopenref))
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
                    
                rec        = np.load("data_save/%s/rec_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))    
                solplot    = np.load("data_save/%s/solplot_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                rec_select = np.load("data_save/%s/rec_select_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                    
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
                        
                if(ptype==1 and dt_ref==6): 
                            
                    recrefcut    = recrefcut[1:,:]
                    recselectcut = recselectcut[1:,:]
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
                    dt0     = 0.5/dt_ref
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
                            
                    sou        = teste.sou    
                    mvalue     = teste.mvalue  
                    nvalue     = teste.nvalue  
                    mshape     = teste.shape   
                    method     = teste.method
                    rec        = np.load("data_save/%s/rec_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))    
                    solplot    = np.load("data_save/%s/solplot_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                    rec_select = np.load("data_save/%s/rec_select_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                    
                    normrec = []
                    
                    try:

                        normrec1      = la.norm(recrefcut-rec,1)
                        normrec2      = la.norm(recrefcut-rec,2)
                        normrecmax    = la.norm(recrefcut-rec,np.inf)                        
                        normrecrel1   = la.norm(recrefcut-rec,1)/la.norm(recrefcut,1)            
                        normrecrel2   = la.norm(recrefcut-rec,2)/la.norm(recrefcut,2)
                        normrecrelmax = la.norm(recrefcut-rec,np.inf)/la.norm(recrefcut,np.inf)
                        
                    except:

                        normrec1      = 'NC'
                        normrec2      = 'NC'
                        normrecmax    = 'NC'                        
                        normrecrel1   = 'NC'
                        normrecrel2   = 'NC'
                        normrecrelmax = 'NC'

                    normrec.append(normrec1)
                    normrec.append(normrec2)
                    normrec.append(normrecmax)                        
                    normrec.append(normrecrel1)
                    normrec.append(normrecrel2)
                    normrec.append(normrecrelmax)

                    normrecselect = [] 
                    n1            = []
                    n2            = []
                    nmax          = []
                    n1rel         = []
                    n2rel         = []
                    nmaxrel       = []
                    
                    for i in range(0,rec_select.shape[1]):

                        try:
                            
                            normloc1      = la.norm(recselectcut[:,i]-rec_select[:,i],1)
                            normloc2      = la.norm(recselectcut[:,i]-rec_select[:,i],2)  
                            normlocmax    = la.norm(recselectcut[:,i]-rec_select[:,i],np.inf)  
                            normlocrel1   = la.norm(recselectcut[:,i]-rec_select[:,i],1)/la.norm(recselectcut[:,i],1)  
                            normlocrel2   = la.norm(recselectcut[:,i]-rec_select[:,i],2)/la.norm(recselectcut[:,i],2)  
                            normlocrelmax = la.norm(recselectcut[:,i]-rec_select[:,i],np.inf)/la.norm(recselectcut[:,i],np.inf)  
                            
                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                            
                        except:
                            
                            normloc1      = 'NC'
                            normloc2      = 'NC'
                            normlocmax    = 'NC'
                            normlocrel1   = 'NC'
                            normlocrel2   = 'NC'
                            normlocrelmax = 'NC'
                            
                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                    
                    normrecselect.append(n1)
                    normrecselect.append(n2)
                    normrecselect.append(nmax)
                    normrecselect.append(n1rel)
                    normrecselect.append(n2rel)
                    normrecselect.append(nmaxrel)
                
                    normsolplot   = []
                    n1            = []
                    n2            = []
                    nmax          = []
                    n1rel         = []
                    n2rel         = []
                    nmaxrel       = []

                    timesolplot = [i*dt0*jump for i in range(0,solplot.shape[0])]
                    
                    for i in range(0,solplot.shape[0]):
                        
                        try:
                            
                            normloc1      = la.norm(solplotcut[i,:]-solplot[i,:],1)
                            normloc2      = la.norm(solplotcut[i,:]-solplot[i,:],2)
                            normlocmax    = la.norm(solplotcut[i,:]-solplot[i,:],np.inf)
                            normlocrel1   = la.norm(solplotcut[i,:]-solplot[i,:],1)/la.norm(solplotcut[i,:],1)
                            normlocrel2   = la.norm(solplotcut[i,:]-solplot[i,:],2)/la.norm(solplotcut[i,:],2)
                            normlocrelmax = la.norm(solplotcut[i,:]-solplot[i,:],np.inf)/la.norm(solplotcut[i,:],np.inf)
                            
                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                            
                        except:
                            
                            normloc1      = 'NC'
                            normloc2      = 'NC'
                            normlocmax    = 'NC'
                            normlocrel1   = 'NC'
                            normlocrel2   = 'NC'
                            normlocrelmax = 'NC'
                            
                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                    
                    normsolplot.append(n1)
                    normsolplot.append(n2)
                    normsolplot.append(nmax)
                    normsolplot.append(n1rel)
                    normsolplot.append(n2rel)
                    normsolplot.append(nmaxrel)
                    
                    testresults.append([cont_me,ptype,dx_ref,dt_ref,freq_ref,mshape,method,mvalue,nvalue,npt,npe,normrec,normrecselect,normsolplot,timesolplot,parameters,cont_glob])
                    
                    cont_glob = cont_glob + 1
#==============================================================================

#==============================================================================
                print('')
                print('Close Test Files!')
                print('')
#==============================================================================
 
#==============================================================================
# Save Results
#==============================================================================
    locname     = 'testresults/test%d_results'%(ptype)
    
    with open(locname, 'wb') as f: 
        pickle.dump(testresults, f) 
#==============================================================================