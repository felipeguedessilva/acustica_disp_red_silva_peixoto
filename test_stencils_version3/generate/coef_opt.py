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
import testes_opt              as ttopt
import macustica               as mc
import sys
#==============================================================================

#==============================================================================
# FD Coefficients Generator
#==============================================================================
sys.path.insert(0, './coef_opt/')
import ger_coef_opt as coefdispred
#==============================================================================

#==============================================================================
class coefopt1:
#==============================================================================

#==============================================================================    
    def __init__(self,teste,MV,ptype):
        self.teste = teste
        self.ptype = ptype
        self.MV    = mc.acusdevito(self.teste,ptype)
#==============================================================================

#==============================================================================
    def calccoef(self,method,shape,mvalue,nvalue,cur):

        mcoef = coefdispred.calccoef(method,shape,mvalue,nvalue,cur)
            
        return mcoef
#==============================================================================

#==============================================================================
    def eqconstuct1(self,mcoef,u,t,x,y):
        
        npx      = mcoef.shape[0]
        npy      = mcoef.shape[1]
        npxm     = int(npx/2)
        npym     = int(npy/2)
        initialx = -npxm
        initialy =  npym
        pdeaux   = 0
        contcoef = 0 
        
        for i in range(0,npx):
            
            for j in range(0,npy):
                                
                a   = int(initialx)
                b   = int(initialy)
                pxs = x + a
                pys = y + b                
                                
                if(mcoef[i,j]!=0): contcoef = contcoef + 1
                
                pdeaux = pdeaux + u[t,pxs,pys]*mcoef[i,j]
                                
                initialx = initialx + 1

            initialx = -npxm
            initialy =  initialy - 1

        return pdeaux, contcoef
#==============================================================================