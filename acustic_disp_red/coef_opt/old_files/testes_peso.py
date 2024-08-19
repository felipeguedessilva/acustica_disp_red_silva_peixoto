#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy                   as np
import sys
import matplotlib.pyplot       as plot
import time                    as tm
#==============================================================================

#==============================================================================
# Shell Comands
#==============================================================================
import os
os.system('clear')
#==============================================================================

#==============================================================================
# Print Configurations
#==============================================================================
np.set_printoptions(formatter={'float': '{: 0.3e}'.format})
#==============================================================================

#==============================================================================
# Basic Coefs
#==============================================================================
import ger_coef  as coefdispred
import ger_coef1 as coefdispred1
import ger_coef2 as coefdispred2
import ger_coef3 as coefdispred3
#==============================================================================

#==============================================================================
plot.close("all")
print_allow = 1
plot_allow  = 0
#==============================================================================

#==============================================================================
# Coef Dispersion Reduction
#==============================================================================
print('')
print('=========================================================================')
print('Dispersion Reduction Scheme: ')
print('=========================================================================')
print('')

mvalue    = 2
nvalue    = 1#np.random.randint(low=1,high=mvalue+1,size=1)[0]
dx        = 30
dt        = 0.5
vmax      = 5.0
cur       = (vmax*dt)/(dx)
vshape    = ['cl','rb','crb','sq','csq']
vmethod   = ['spatte','specte','specls','dispte','spectetheta','displs'] 
shape     = vshape[0]
method    = vmethod[0]

if(print_allow==1):
        
    print('')
    print('=========================================================================')
    print('Shape : %s'%shape)
    print('Method: %s'%method)
    print('mvalue: %d'%mvalue)
    print('nvalue: %d'%nvalue)
    print('=========================================================================')

start0 = tm.time()
T0     = coefdispred.calccoef(method,shape,mvalue,nvalue,cur)
end0   = tm.time()
tt0    = end0 - start0
print('Tempo Total T0: %f'%tt0)

# start0 = tm.time()
# T0     = coefdispred1.calccoef(method,shape,mvalue,nvalue,cur)
# end0   = tm.time()
# tt0    = end0 - start0
# print('Tempo Total T0: %f'%tt0)

# start1 = tm.time()
# T1     = coefdispred2.calccoef(method,shape,mvalue,nvalue,cur)
# end1   = tm.time()
# tt1    = end1 - start1
# print('Tempo Total T1: %f'%tt1)

#start2 = tm.time()
#T2     = coefdispred3.calccoef(method,shape,mvalue,nvalue,cur)
#end2   = tm.time()
#tt2    = end2 - start2
#print('Tempo Total T2: %f'%tt2)
#==============================================================================

#==============================================================================
if(print_allow==1):
            
    print('')
    print('=========================================================================')
    print(T0)
    print('=========================================================================')
    print('')
#==============================================================================

#==============================================================================
if(plot_allow==1):

    fig1 = plot.figure(figsize = (12,12))
    grid = plot.GridSpec(1,1,wspace=0.4,hspace=0.2)
    plot.suptitle('Stencil Geometry')
    position_plot_listx = np.array([0,0])
    position_plot_listy = np.array([0,1])
    xpos = int(position_plot_listx[0])
    ypos = int(position_plot_listy[0])    
    plot.subplot(grid[xpos,ypos])
    plot.spy(T0,aspect='equal',markersize=8)
    plot.title('Method = %s - Shape = %s - mvalue = %d - nvalue = %d'%(shape,method,mvalue,nvalue))
    ax = plot.gca()
    plot.grid()
    vtikx = np.arange(-mvalue,mvalue+1) 
    ax.set_xticks(np.arange(0,len(vtikx)))    
    ax.axes.xaxis.set_ticklabels(vtikx)
    ax.set_yticks(np.arange(0,len(vtikx)))
    ax.axes.yaxis.set_ticklabels(np.flip(vtikx))
    ax.xaxis.set_ticks_position('bottom')
    plot.xlabel('[X Direction]')
    plot.ylabel('[Z Direction]')
    plot.show()
#==============================================================================