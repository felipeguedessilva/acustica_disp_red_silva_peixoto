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
def datasave(teste,rec,solplot,rec_select,ref,ptype,dttype):
    
    sou    = teste.sou    
    mvalue = teste.mvalue  
    nvalue = teste.nvalue  
    mshape = teste.shape   
    method = teste.method
    
    loc_save1_list = ['teste1/','teste2/','teste3/','teste4/']
    loc_save2_list = ['dt1/','dt2/','dt3/','dt4/']
    
    nteste = len(loc_save1_list)
    ntimes = len(loc_save2_list)
    
    if(ptype==1): local_save1 = loc_save1_list[0]
    if(ptype==2): local_save1 = loc_save1_list[1]
    if(ptype==3): local_save1 = loc_save1_list[2]
    if(ptype==4): local_save1 = loc_save1_list[3]

    if(ref==0):
        
        if(dttype==0): local_save2 = loc_save2_list[0]
        if(dttype==1): local_save2 = loc_save2_list[1]
        if(dttype==2): local_save2 = loc_save2_list[2]
        if(dttype==3): local_save2 = loc_save2_list[3]
     
        np.save("data_save/%s%srec_%s_%s_%d_%d"%(local_save1,local_save2,mshape,method,mvalue,nvalue),rec)    
        np.save("data_save/%s%ssolplot_%s_%s_%d_%d"%(local_save1,local_save2,mshape,method,mvalue,nvalue),solplot)
        np.save("data_save/%s%srec_select_%s_%s_%d_%d"%(local_save1,local_save2,mshape,method,mvalue,nvalue),rec_select)
        
    if(ref==1):
        
        for i in range(0,ntimes):
            
            local_save2 = loc_save2_list[i]
            np.save("data_save/%s%srec_ref"%(local_save1,local_save2,),rec)   
            np.save("data_save/%s%ssolplot_ref"%(local_save1,local_save2,),solplot)            
            np.save("data_save/%s%srec_select_ref"%(local_save1,local_save2,),rec_select)   
#==============================================================================

#==============================================================================
# Plot do Deslocamento
#==============================================================================
def graph2d(U,teste,i):

    x0     = teste.x0
    x1     = teste.x1
    y0     = teste.y0
    y1     = teste.y1
    sou    = teste.sou    
    mvalue = teste.mvalue  
    nvalue = teste.nvalue  
    mshape = teste.shape   
    method = teste.method  
    
    plot.figure(figsize = (14,10))
    fscale =  10**(-3)

    scale = np.amax(U)/10
    extent = [fscale*x0,fscale*x1, fscale*y1, fscale*y0]
    fig = plot.imshow(np.transpose(U), vmin=-scale, vmax=scale, cmap=cm.gray, extent=extent,interpolation='bicubic', aspect='auto')

    plot.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
    plot.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
    plot.axis('equal')
    
    if(i==0): plot.title('Shape = %s - Method = %s - mvalue = %d - nvalue = %d'%(mshape,method,mvalue,nvalue),fontsize=10)
    if(i==1): plot.title('Rererence Solution',fontsize=10)
    
    plot.grid()
    ax = plot.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.025)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar = plot.colorbar(fig, cax=cax, format='%.2e')
    cbar.locator = tick_locator
    cbar.update_ticks()
    plot.draw()
    
    if(i==0): plot.savefig('figures/temporario/disp_%s_%s_%d_%d.png'%(mshape,method,mvalue,nvalue),dpi=100)
    if(i==1): plot.savefig('figures/temporario/disp_ref.png',dpi=100)
    
    plot.show()
    plot.close()
#==============================================================================

#==============================================================================
# Plot dos Receivers
#==============================================================================
def graph2drec(rec,teste,i):
    
    x0 = teste.x0
    x1 = teste.x1
    y0 = teste.y0
    y1 = teste.y1
    t0 = teste.t0
    tn = teste.tn
 
    sou    = teste.sou    
    mvalue = teste.mvalue  
    nvalue = teste.nvalue  
    mshape = teste.shape   
    method = teste.method  

    plot.figure(figsize = (14,10))
    fscale =  10**(-3)    
    scale = np.amax(rec)/50.
    extent = [fscale*x0,fscale*x1, fscale*tn, fscale*t0]
    fig = plot.imshow(rec, vmin=-scale, vmax=scale, cmap=cm.gray, extent=extent)
        
    plot.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
    plot.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
    plot.axis('equal')
    if(i==0): plot.title('Shape = %s - Method = %s - mvalue = %d - nvalue = %d'%(mshape,method,mvalue,nvalue),fontsize=10)
    if(i==1): plot.title('Rererence Solution',fontsize=10)
    plot.grid()
    ax = plot.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.025)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar = plot.colorbar(fig, cax=cax, format='%.2e')
    cbar.locator = tick_locator
    cbar.update_ticks()
    if(i==0): plot.savefig('figures/temporario/rec_%s_%s_%d_%d.png'%(mshape,method,mvalue,nvalue),dpi=100)
    if(i==1): plot.savefig('figures/temporario/rec_ref.png')
    plot.show()
    plot.close()
#==============================================================================

#==============================================================================
# Plot Velocidades
#==============================================================================
def graph2dvel(vel,teste):
    
    x0 = teste.x0
    x1 = teste.x1
    y0 = teste.y0
    y1 = teste.y1        
  
    plot.figure(figsize = (14,10))
    fscale =  10**(-3)
           
    scale  = np.amax(vel)
    extent = [fscale*x0,fscale*x1, fscale*y1, fscale*y0]
    fig = plot.imshow(np.transpose(vel), vmin=0,vmax=scale, cmap=cm.jet, extent=extent)
          
    plot.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
    plot.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
    plot.title('Velocity Profile',fontsize=10)
    plot.grid()
    ax = plot.gca()
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.025)
    tick_locator = ticker.MaxNLocator(nbins=5)
    cbar = plot.colorbar(fig, cax=cax, format='%.2e')
    cbar.locator = tick_locator
    cbar.update_ticks()
    cbar.set_label('Velocity [km/s]')
    plot.savefig('figures/temporario/vel_map.png')
    plot.show()
    plot.close()
#==============================================================================