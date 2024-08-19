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
import matplotlib.pyplot       as plot
import matplotlib.ticker       as mticker    
from   scipy.interpolate       import CubicSpline
from   mpl_toolkits.mplot3d    import Axes3D
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib              import cm
from   matplotlib.animation    import FFMpegWriter
from   matplotlib              import ticker
from   matplotlib.colors       import LogNorm
import matplotlib.tri          as tri
#==============================================================================

#==============================================================================
# Data Save
#==============================================================================        
def datasave(teste,rec,solplot,rec_select,ref,ptype,locsave):
    
    sou    = teste.sou    
    mvalue = teste.mvalue  
    nvalue = teste.nvalue  
    mshape = teste.shape   
    method = teste.method
    
    if(ref==0):
          
        np.save("data_save/%s/rec_%s_%s_%d_%d"%(locsave,mshape,method,mvalue,nvalue),rec)    
        np.save("data_save/%s/solplot_%s_%s_%d_%d"%(locsave,mshape,method,mvalue,nvalue),solplot)
        np.save("data_save/%s/rec_select_%s_%s_%d_%d"%(locsave,mshape,method,mvalue,nvalue),rec_select)
        
    if(ref==1):
        
        lsave = ['teste%d/dx1dt1freq%d'%(ptype,locsave),
                 'teste%d/dx1dt2freq%d'%(ptype,locsave),
                 'teste%d/dx1dt4freq%d'%(ptype,locsave),
                 'teste%d/dx1dt6freq%d'%(ptype,locsave),
                 'teste%d/dx2dt1freq%d'%(ptype,locsave),
                 'teste%d/dx2dt2freq%d'%(ptype,locsave),
                 'teste%d/dx2dt4freq%d'%(ptype,locsave),
                 'teste%d/dx2dt6freq%d'%(ptype,locsave),
                 'teste%d/dx4dt1freq%d'%(ptype,locsave),
                 'teste%d/dx4dt2freq%d'%(ptype,locsave),
                 'teste%d/dx4dt4freq%d'%(ptype,locsave),
                 'teste%d/dx4dt6freq%d'%(ptype,locsave),
                 'teste%d/dx8dt1freq%d'%(ptype,locsave),
                 'teste%d/dx8dt2freq%d'%(ptype,locsave),
                 'teste%d/dx8dt4freq%d'%(ptype,locsave),
                 'teste%d/dx8dt6freq%d'%(ptype,locsave)
                 ]
        
        nsaves = len(lsave)
        
        for i in range(0,nsaves):
        
            locsavemod = lsave[i]
        
            np.save("data_save/%s/rec_ref"%(locsavemod),rec)   
            np.save("data_save/%s/solplot_ref"%(locsavemod),solplot)            
            np.save("data_save/%s/rec_select_ref"%(locsavemod),rec_select)   
#==============================================================================