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
          
        np.save("../data_save/%s/rec_%s_%s_%d_%d"%(locsave,mshape,method,mvalue,nvalue),rec)    
        np.save("../data_save/%s/solplot_%s_%s_%d_%d"%(locsave,mshape,method,mvalue,nvalue),solplot)
        np.save("../data_save/%s/rec_select_%s_%s_%d_%d"%(locsave,mshape,method,mvalue,nvalue),rec_select)
        
    if(ref==1):
        
        locsavemod = 'teste%d/reffreq%d'%(ptype,locsave)
        np.save("../data_save/%s/rec_ref"%(locsavemod),rec)   
        np.save("../data_save/%s/solplot_ref"%(locsavemod),solplot)            
        np.save("../data_save/%s/rec_select_ref"%(locsavemod),rec_select)   
#==============================================================================