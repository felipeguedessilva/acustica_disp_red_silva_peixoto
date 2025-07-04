#==============================================================================
# -*- encoding: utf-8 -*-
#==============================================================================

#==============================================================================
# Módulos Importados do Python / Devito / Examples
#==============================================================================

#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy                   as np
import math                    as mt
import sys
from   numpy                   import linalg as la
#==============================================================================

#==============================================================================
# Testes de Leitura de Dados
#==============================================================================
ptype    = 5
normtype = 2
#==============================================================================

#==============================================================================
# Vetores de Configurações
#==============================================================================
vmethod0 = np.array([0,1,1,1,1,1,1,1,1,1])
vmethod1 = np.array([1,1,2,3,3,3,3,3,3,4])
vmethod2 = np.array([1,4,2,1,2,1,2,1,2,1])
vmethod3 = np.array([2,4,6,8,10,12,14,16,18,20])
nteste   = vmethod3.shape[0]

l1    = np.zeros(nteste)
l2    = np.zeros(nteste)
l3    = np.zeros(nteste)
l4    = np.zeros(nteste)
l5    = np.zeros(nteste)
l6    = np.zeros(nteste)
l7    = np.zeros(nteste)
l8    = np.zeros(nteste)
l9    = np.zeros(nteste)
l10   = np.zeros(nteste)

for i in range(0,nteste):
    
    l1[i]  = 1
    l2[i]  = 1
    l3[i]  = 1
    l4[i]  = 1
    l5[i]  = 1
    l6[i]  = int(vmethod3[i]/2)
    l7[i]  = int(vmethod3[i]/2)
    l8[i]  = int(vmethod3[i]/4 + 1)
    l9[i]  = int(vmethod3[i]/4 + 1)
    l10[i] = 1

vmethod4 = [l1,l2,l3,l4,l5,l6,l7,l8,l9,l10]

s1    = np.zeros(nteste)
s2    = np.zeros(nteste)
s3    = np.zeros(nteste)
s4    = np.zeros(nteste)
s5    = np.zeros(nteste)
s6    = np.zeros(nteste)
s7    = np.zeros(nteste)
s8    = np.zeros(nteste)
s9    = np.zeros(nteste)
s10   = np.zeros(nteste)

for i in range(0,nteste):
    
    s1[i]  = 1
    s2[i]  = 1
    s3[i]  = 1
    s4[i]  = 1
    s5[i]  = 1
    s6[i]  = 1
    s7[i]  = 1
    s8[i]  = 1
    s9[i]  = 1
    s10[i] = int(0.5*(int(0.5*(int(vmethod3[i]/2)**2+int(vmethod3[i]/2)))))

vmethod5 = [s1,s2,s3,s4,s5,s6,s7,s8,s9,s10]

total_configs = 0
list_config   = []

for i in range(0,vmethod0.shape[0]):
    
    scheme   = i
    peso     = vmethod0[i]
    wauthor  = vmethod1[i]
    wtype    = vmethod2[i]
    vnvalue  = vmethod4[i]
    vexttrap = vmethod5[i]
    
    for l in range(0,vmethod3.shape[0]):
        
        mvalue  = vmethod3[l]
        nvalue  = int(vnvalue[l])
        exttrap = int(vexttrap[l])
        config  = (peso,wauthor,wtype,mvalue,nvalue,exttrap,scheme)                            
        total_configs = total_configs + 1
        list_config.append(config)

list_config = list(set(list_config))
nconfig     = len(list_config)
vtime_exec  = np.zeros(nconfig)
#==============================================================================

#==============================================================================
# Looping Sobre os dt's
#==============================================================================
ndts   = 6
dttype = 0

for root in range(0,ndts):

    dttype = root + 1
#==============================================================================

#==============================================================================        
    if(ptype==1): 
        
        nptx        = 201
        npty        = 201
        dx          = 15
        dz          = 15    
        t0          = 0      
        tn          = 3000
        nrefesp     = 5
        xpositionv  = np.array([750.0,2250.0, 750.0,2250.0])
        ypositionv  = np.array([750.0, 750.0,2250.0,2250.0])
        ratiojump   = 150

        if(dttype==1): 
            
            datasave    = 'data_save/teste1/dt1/'
            nt          = 6000
            nreftemp    = 5
        
        if(dttype==2): 
            
            datasave    = 'data_save/teste1/dt2/'
            nt          = 3000
            nreftemp    = 10

        if(dttype==3): 
            
            datasave    = 'data_save/teste1/dt3/'
            nt          = 2000
            nreftemp    = 15

        if(dttype==4): 
            
            datasave    = 'data_save/teste1/dt4/'
            nt          = 1500
            nreftemp    = 20

        if(dttype==5): 
            
            datasave    = 'data_save/teste1/dt5/'
            nt          = 1200
            nreftemp    = 25
        
        if(dttype==6): 
            
            datasave    = 'data_save/teste1/dt6/'
            nt          = 1000
            nreftemp    = 30
#==============================================================================

#==============================================================================        
    if(ptype==2): 
        
        nptx        = 201
        npty        = 201
        dx          = 10
        dz          = 10    
        t0          = 0      
        tn          = 3000
        nrefesp     = 5
        xpositionv  = np.array([500.0,1500.0,500.0,1500.0])
        ypositionv  = np.array([500.0,500.0,1500.0,1500.0])
        ratiojump   = 180

        if(dttype==1): 
            
            datasave    = 'data_save/teste2/dt1/'
            nt          = 6000
            nreftemp    = 5
        
        if(dttype==2): 
            
            datasave    = 'data_save/teste2/dt2/'
            nt          = 3000
            nreftemp    = 10

        if(dttype==3): 
            
            datasave    = 'data_save/teste2/dt3/'
            nt          = 2000
            nreftemp    = 15

        if(dttype==4): 
            
            datasave    = 'data_save/teste2/dt4/'
            nt          = 1500
            nreftemp    = 20

        if(dttype==5): 
            
            datasave    = 'data_save/teste2/dt5/'
            nt          = 1200
            nreftemp    = 25
        
        if(dttype==6): 
            
            datasave    = 'data_save/teste2/dt6/'
            nt          = 1000
            nreftemp    = 30
#==============================================================================

#==============================================================================
    if(ptype==3): 
        
        nptx        = 601 
        npty        = 201
        dx          = 20
        dz          = 20    
        t0          = 0      
        tn          = 3000
        nrefesp     = 5
        xpositionv  = np.array([4000.0,4000.0,4000.0,6000.0,6000.0,6000.0,8000.0,8000.0,8000.0])   
        ypositionv  = np.array([2000.0,2500.0,3000.0,2000.0,2500.0,3000.0,2000.0,2500.0,3000.0])    
        ratiojump   = 150

        if(dttype==1): 
            
            datasave    = 'data_save/teste3/dt1/'
            nt          = 6000
            nreftemp    = 5

        if(dttype==2): 
            
            datasave    = 'data_save/teste3/dt2/'
            nt          = 3000
            nreftemp    = 10

        if(dttype==3): 
            
            datasave    = 'data_save/teste3/dt3/'
            nt          = 2000
            nreftemp    = 15

        if(dttype==4): 
            
            datasave    = 'data_save/teste3/dt4/'
            nt          = 1500
            nreftemp    = 20

        if(dttype==5): 
            
            datasave    = 'data_save/teste3/dt5/'
            nt          = 1200
            nreftemp    = 25
        
        if(dttype==6): 
            
            datasave    = 'data_save/teste3/dt6/'
            nt          = 1000
            nreftemp    = 30
#==============================================================================        

#==============================================================================        
    if(ptype==4): 
        
        nptx        = 626 
        npty        = 311
        dx          = 32  
        dz          = 32    
        t0          = 0      
        tn          = 3000
        nrefesp     = 5
        xpositionv  = np.array([30000.0,30000.0,30000.0,40000.0,40000.0,40000.0])
        ypositionv  = np.array([2500.0,5000.0,7500.0,2500.0,5000.0,7500.0])
        ratiojump   = 150

        if(dttype==1): 
            
            datasave    = 'data_save/teste4/dt1/'
            nt          = 6000
            nreftemp    = 5

        if(dttype==2): 
            
            datasave    = 'data_save/teste4/dt2/'
            nt          = 3000
            nreftemp    = 10

        if(dttype==3): 
            
            datasave    = 'data_save/teste4/dt3/'
            nt          = 2000
            nreftemp    = 15

        if(dttype==4): 
            
            datasave    = 'data_save/teste4/dt4/'
            nt          = 1500
            nreftemp    = 20

        if(dttype==5): 
            
            datasave    = 'data_save/teste4/dt5/'
            nt          = 1200
            nreftemp    = 25
        
        if(dttype==6): 
            
            datasave    = 'data_save/teste4/dt6/'
            nt          = 1000
            nreftemp    = 30
#==============================================================================

#==============================================================================        
    if(ptype==5): 
        
        nptx        = 901 
        npty        = 321
        dx          = 10  
        dz          = 10    
        t0          = 0      
        tn          = 3000
        nrefesp     = 5
        xpositionv  = np.array([6000.0,6000.0,6000.0,8000.0,8000.0,8000.0,10000.0,10000.0,10000.0,12000.0,12000.0,12000.0])
        ypositionv  = np.array([1000.0,2000.0,3000.0,1000.0,2000.0,3000.0,1000.0,2000.0,3000.0,1000.0,2000.0,3000.0])
        ratiojump   = 150

        if(dttype==1): 
            
            datasave    = 'data_save/teste5/dt1/'
            nt          = 6000
            nreftemp    = 5

        if(dttype==2): 
            
            datasave    = 'data_save/teste5/dt2/'
            nt          = 3000
            nreftemp    = 10

        if(dttype==3): 
            
            datasave    = 'data_save/teste5/dt3/'
            nt          = 2000
            nreftemp    = 15

        if(dttype==4): 
            
            datasave    = 'data_save/teste5/dt4/'
            nt          = 1500
            nreftemp    = 20

        if(dttype==5): 
            
            datasave    = 'data_save/teste5/dt5/'
            nt          = 1200
            nreftemp    = 25
        
        if(dttype==6): 
            
            datasave    = 'data_save/teste5/dt6/'
            nt          = 1000
            nreftemp    = 30
#==============================================================================

#==============================================================================

#==============================================================================
# Carregando Soluções
#==============================================================================
    for m in range(0,nconfig):
#==============================================================================

#==============================================================================
        config  = list_config[m]
        npesos  = int(config[0])          
        wauthor = int(config[1]) 
        wtype   = int(config[2])        
        mvalue  = int(config[3]/2)           
        nvalue  = int(config[4])
        exttrap = int(config[5])
#==============================================================================

#==============================================================================
        if(wauthor!=4):

            solplot_num2      = np.load('%ssolplot_%d%d%d%d%d.npy'%(datasave,npesos,wauthor,wtype,mvalue,nvalue))    
            rec_num2          = np.load('%srec_%d%d%d%d%d.npy'%(datasave,npesos,wauthor,wtype,mvalue,nvalue))    
            rec_select_num2   = np.load('%srec_select_%d%d%d%d%d.npy'%(datasave,npesos,wauthor,wtype,mvalue,nvalue))    

        else:
            
            solplot_num2      = np.load('%ssolplot_%d%d%d%d%d.npy'%(datasave,npesos,wauthor,wtype,mvalue,exttrap))    
            rec_num2          = np.load('%srec_%d%d%d%d%d.npy'%(datasave,npesos,wauthor,wtype,mvalue,exttrap))    
            rec_select_num2   = np.load('%srec_select_%d%d%d%d%d.npy'%(datasave,npesos,wauthor,wtype,mvalue,exttrap))    
#==============================================================================

#==============================================================================
        if(m==0):

            solplot_ref_read     = np.load('%ssolplot_ref.npy'%datasave)
            rec_ref_read         = np.load('%srec_ref.npy'%datasave)
            rec_select_ref_read  = np.load('%srec_select_ref.npy'%datasave)
        
            solplot = np.zeros((nconfig,solplot_num2.shape[0],solplot_num2.shape[1],solplot_num2.shape[2]))
            rec_num = np.zeros((nconfig,rec_num2.shape[0],rec_num2.shape[1]))
            rec_select_num = np.zeros((nconfig,rec_select_num2.shape[0],rec_select_num2.shape[1]))

            solplot_ref    = np.zeros((solplot_num2.shape[0],solplot_num2.shape[1],solplot_num2.shape[2]))
            rec_ref        = np.zeros((rec_num2.shape[0],rec_num2.shape[1]))
            rec_select_ref = np.zeros((rec_select_num2.shape[0],rec_select_num2.shape[1]))
    
            nconty = 0

            for j in range(0,npty):
  
                ncontx = 0
    
                for i in range(0,nptx):
        
                    solplot_ref[:,i,j] = solplot_ref_read[:,ncontx,nconty]

                    ncontx = ncontx + nrefesp
            
                nconty = nconty + nrefesp
        
            ntimes     = rec_num2.shape[0]
            nconttempo = 0

            for j in range(0,ntimes):

                ncontx = 0
        
                if(j==0):    
        
                    for i in range(0,nptx):
                
                        rec_ref[0,i] = rec_ref_read[0,ncontx]
            
                        ncontx = ncontx + nrefesp
        
                    rec_select_ref[0,:] = rec_select_ref_read[0,:]
        
                    nconttempo = nconttempo + 0

                else:
                      
                    for i in range(0,nptx):
                
                        rec_ref[j,i] = rec_ref_read[nconttempo,ncontx]
            
                        ncontx = ncontx + nrefesp
                
                    rec_select_ref[j,:] = rec_select_ref_read[nconttempo,:]

                    nconttempo = nconttempo + nreftemp
            
            ntimes_disp = solplot_num2.shape[0]
            timev_disp  = np.append(np.arange(t0,tn,ratiojump),tn) 
            #timev_disp  = np.linspace(t0,tn,ntimes_disp)
            ntimes_rec  = rec_num2.shape[0]
            timev_rec   = np.linspace(t0,tn,ntimes_rec)
            ordersv     = np.zeros(nteste)

            for i in range(0,nteste):
    
                ordersv[i] = 2*(i+1)
        
            mnormas_disp = np.zeros((nconfig,solplot_num2.shape[0]))
            mnormas_rec  = np.zeros((nconfig,rec_num2.shape[0]))
#==============================================================================

#==============================================================================
        solplot[m,:,:,:]      = solplot_num2[:,:,:]
        rec_num[m,:,:]        = rec_num2[:,:]
        rec_select_num[m,:,:] = rec_select_num2[:,:]
#==============================================================================

#==============================================================================
    for i in range(0,nconfig):
    
        for j in range(0,ntimes_disp):
        
            v1                = solplot_ref[j,:,:]
            v2                = solplot[i,j,:,:]
            nt_val            = np.sqrt(nt)*max(dx,dz)           
            a                 = np.linalg.norm(np.reshape(v1-v2,-1),normtype)/nt_val
            mnormas_disp[i,j] = a
        
    for i in range(0,nconfig):
    
        for j in range(0,ntimes_rec):
        
            v1               = rec_ref[j,:]
            v2               = rec_num[i,j,:]
            nt_val           = np.sqrt(nt)*max(dx,dz)
            a                = np.linalg.norm(np.reshape(v1-v2,-1),normtype)/nt_val
            mnormas_rec[i,j] = a
#==============================================================================

#==============================================================================
    nposition      = xpositionv.shape[0]
    solplot_select = np.zeros((rec_select_num.shape[0],rec_select_num.shape[2],rec_select_num.shape[1],))

    for i in range(0,nconfig):

        for k in range(0,nposition):
            
            solplot_select[i,k,:] = rec_select_num[i,:,k]
        
    mnormas_disp_select = np.zeros((nconfig,nposition))

    for i in range(0,nconfig):

        for j in range(0,nposition):
    
            dtv         = np.array([0.5,1.0,1.5,2.0,2.5,3.0])
            if(ptype==1): time_select = 900
            if(ptype==2): time_select = 540
            if(ptype==3): time_select = 1800
            if(ptype==4): time_select = 2550
            if(ptype==5): time_select = 1800
            
            posfinal    = int(time_select/dtv[root])+1
            v1          = rec_select_ref[0:posfinal,j]
            v2          = solplot_select[i,j,0:posfinal]
            nt_val      = np.sqrt(nt)*max(dx,dz) 
            a           = np.linalg.norm(np.reshape(v1-v2,-1),normtype)/nt_val
            
            mnormas_disp_select[i,j] = a
#==============================================================================

#==============================================================================
# Save Data 
#==============================================================================
    np.save("%smnormas_disp_%d"%(datasave,dttype),mnormas_disp)
    np.save("%smnormas_rec_%d"%(datasave,dttype),mnormas_rec)
    np.save("%smnormas_disp_select_%d"%(datasave,dttype),mnormas_disp_select)
    np.save("%srec_num_%d"%(datasave,dttype),rec_num)
    np.save("%ssolplot_%d"%(datasave,dttype),solplot)
    np.save("%srec_select_num_%d"%(datasave,dttype),rec_select_num)
    np.save("%ssolplot_select_%d"%(datasave,dttype),solplot_select)
    np.save("%ssolplot_ref_%d"%(datasave,dttype),solplot_ref)
    np.save("%srec_ref_%d"%(datasave,dttype),rec_ref)
    np.save("%srec_select_ref_%d"%(datasave,dttype),rec_select_ref)
    np.save("%stimev_disp_%d"%(datasave,dttype),timev_disp)
    np.save("%stimev_rec_%d"%(datasave,dttype),timev_rec)
    np.save("%sordersv_%d"%(datasave,dttype),ordersv)
#==============================================================================