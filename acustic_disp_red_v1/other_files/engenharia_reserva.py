#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy as np
import math  as mt 
#==============================================================================

#==============================================================================
# Test Parameters
#==============================================================================
ptype  = 1
ref    = 0

if(ptype==1):
    
    dx         = 30.0
    dz         = 30.0
    vmax       = 3.0
    t0         = 0.0
    tn         = 3000.0
    CFLv       = np.array([0.05])
    vdt        = np.array([0.5])
    jumpv      = np.array([300])
    ref_factor = 5
    
    if(ref==1): 
        
        dx    = dx/ref_factor
        dz    = dz/ref_factor
        CFLv  = np.array([0.05])
        vdt   = np.array([0.1])
        jumpv = np.array([1500])

if(ptype==2):
    
    dx         = 20.0
    dz         = 20.0
    vmax       = 3.0
    t0         = 0.0
    tn         = 3000.0
    CFLv       = np.array([0.075])
    vdt        = np.array([0.5])
    jumpv      = np.array([360])
    ref_factor = 5
    
    if(ref==1): 
        
        dx    = dx/ref_factor
        dz    = dz/ref_factor
        CFLv  = np.array([0.075])
        vdt   = np.array([0.1])
        jumpv = np.array([1800])

if(ptype==3):
    
    dx         = 40.0
    dz         = 40.0
    vmax       = 5.0
    t0         = 0.0
    tn         = 3000.0
    CFLv       = np.array([0.0625])
    vdt        = np.array([0.5])
    jumpv      = np.array([300])
    ref_factor = 5
    
    if(ref==1): 
        
        dx    = dx/ref_factor
        dz    = dz/ref_factor
        CFLv  = np.array([0.0625])
        vdt   = np.array([0.1])
        jumpv = np.array([1500])
        
if(ptype==4):
    
    dx         = 40.0
    dz         = 40.0
    vmax       = 5.0
    t0         = 0.0
    tn         = 3000.0
    CFLv       = np.array([0.0625])
    vdt        = np.array([0.5])
    jumpv      = np.array([300])
    ref_factor = 5
    
    if(ref==1): 
        
        dx    = dx/ref_factor
        dz    = dz/ref_factor
        CFLv  = np.array([0.0625])
        vdt   = np.array([0.1])
        jumpv = np.array([1500])
#==============================================================================

#==============================================================================
# Cacl Delta t
#==============================================================================
ntimes = vdt.shape[0]

for i in range(0,ntimes):
    
    dtmax  = np.round(min(dx,dz)*CFLv[i]/vmax,8)
    ntmax  = int((tn-t0)/dtmax)
    dt0    = (tn-t0)/(ntmax)
    nplot  = mt.ceil(ntmax/jumpv[i]) + 1
    ntjump = jumpv[i]*dt0 
    print(vdt[i],dtmax,ntmax,dt0,nplot,ntjump)
#==============================================================================