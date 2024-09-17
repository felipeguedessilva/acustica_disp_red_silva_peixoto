#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
import pandas                   as pd
from   numpy import linalg      as la
import sys
import pickle 
import math                    as mt
import sys
import time                    as tm
#==============================================================================

#==============================================================================
# Plot Set
#==============================================================================
import matplotlib.pyplot       as plt
import plotly.express          as px
import seaborn                 as sns
import matplotlib.ticker       as mticker    
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib              import ticker
from   matplotlib              import cm
pd.set_option("future.no_silent_downcasting", True)
pd.options.mode.chained_assignment = None
from plotly.offline import plot
#==============================================================================

#==============================================================================
plt.close("all")
#==============================================================================

#=============================================================================
# List Content
#=============================================================================
#cont_me 0 
#ptype 1
#dx_ref 2
#dt_ref 3
#freq_ref 4
#mshape 5
#method 6
#mvalue 7
#nvalue 8
#npt 9 
#npe 10

#normrec 11
#    normrec1 - scalar 110 
#    normrec2  - scalar  111   
#    normrecmax - scalar  112                          
#    normrecrel1 - scalar 113  
#    normrecrel2 - scalar  114
#    normrecrelmax - scalar 115

#normrecselect 12
#    normrec1 - vector with 4 positons 120
#    normrec2  - vector with 4 positons   121 
#    normrecmax - vector with 4 positons  122                          
#    normrecrel1 - vector with 4 positons  123  
#    normrecrel2 - vector with 4 positons  124
#    normrecrelmax - vector with 4 positons 125

#normsolplot 13
#    normrec1 - vector with nsave positons 130
#    normrec2  - vector with nsave positons 131   
#    normrecmax - vector with nsave positons 132                           
#    normrecrel1 - vector with nsave positons 133   
#    normrecrel2 - vector with nsave positons 134 
#    normrecrelmax - vector with nsave positons 135 

#timesolplot14

#parameters 15
#    nptx 150
#    npty 151
#    x0 152
#    y0 153
#    compx 154
#    compy 155
#    hxv 156
#    hyv 157
#    t0 158
#    tn 159
#    f0 1510
#    CFL 1511

#cont_glob 16

# method spatte                 - shape cl
# method spectetheta            - shape cl
# method dispte, specls, displs - shape crb, cl, rb
#=============================================================================

#==============================================================================
# Read Txt Files
#==============================================================================
ptype    = 1
locopen  = '../testresults/test%d_results_norms_fields'%(ptype)

with open(locopen, 'rb') as f: 

    if(ptype==1): 
            
        test1_results = pickle.load(f) 
        ntr1          = len(test1_results)
        testname      = 'Homogeneos Velocity Model'
        xpositionv    = np.array([750.0,2250.0, 750.0,2250.0])
        ypositionv    = np.array([750.0, 750.0,2250.0,2250.0])
            
    elif(ptype==2): 
            
        test2_results = pickle.load(f) 
        ntr2          = len(test2_results)
        testname      = 'Heterogeneos Velocity Model'
        xpositionv    = np.array([500.0,1500.0, 500.0,1500.0])
        ypositionv    = np.array([500.0, 500.0,1500.0,1500.0])
        
    elif(ptype==3): 
            
        test3_results = pickle.load(f) 
        ntr3          = len(test3_results)
        testname      = 'Marmousi Velocity Model'
        xpositionv    = np.array([4000.0,4000.0,4000.0,6000.0,6000.0,6000.0,8000.0,8000.0,8000.0])   
        ypositionv    = np.array([2000.0,2500.0,3000.0,2000.0,2500.0,3000.0,2000.0,2500.0,3000.0]) 

    elif(ptype==4): 
            
        test4_results = pickle.load(f)
        ntr4          = len(test4_results)
        testname      = 'SEG/EAGE 2D Salt Velocity Model'
        xpositionv    = np.array([30000.0,30000.0,30000.0,40000.0,40000.0,40000.0])
        ypositionv    = np.array([2500.0,5000.0,7500.0,2500.0,5000.0,7500.0])
#==============================================================================

#==============================================================================
# Data Frame
#==============================================================================
ntest1_results = len(test1_results)

vdx      = []
vdt      = []
vfreq    = []
vshape   = []
vmethod  = []
vmvalue  = []
vnvalue  = []
vnpt     = []
vnpe     = []
vnormrec = []
vnormsol = []
vdxv     = []
vdtv     = []
vfreqv   = []

for k1 in range(0,ntest1_results):
    
    vdx.append(test1_results[k1][2])
    vdt.append(test1_results[k1][3])
    vfreq.append(test1_results[k1][4])
    vshape.append(test1_results[k1][5])
    vmethod.append(test1_results[k1][6])
    vmvalue.append(test1_results[k1][7])
    vnvalue.append(test1_results[k1][8])
    vnpt.append(test1_results[k1][9])
    vnpe.append(test1_results[k1][10])
    vnormrec.append(test1_results[k1][11][4])
    vnormsol.append(test1_results[k1][13][4][-1])
    vdxv.append(test1_results[k1][15][6])
    vdtv.append(test1_results[k1][15][13])
    vfreqv.append(test1_results[k1][15][10])
    

dfrec = pd.DataFrame()
dfsol = pd.DataFrame()

dfrec['dx']      = vdx
dfrec['dt']      = vdt
dfrec['freq']    = vfreq
dfrec['shape']   = vshape
dfrec['method']  = vmethod
dfrec['mvalue']  = vmvalue
dfrec['nvalue']  = vnvalue
dfrec['npt']     = vnpt
dfrec['npe']     = vnpe
dfrec['normrec'] = vnormrec
dfrec['dxv'] = vdxv
dfrec['dtv'] = vdtv
dfrec['freqv'] = vfreqv

dfsol['dx']      = vdx
dfsol['dt']      = vdt
dfsol['freq']    = vfreq
dfsol['shape']   = vshape
dfsol['method']  = vmethod
dfsol['mvalue']  = vmvalue
dfsol['nvalue']  = vnvalue
dfsol['npt']     = vnpt
dfsol['npe']     = vnpe
dfsol['normsol'] = vnormsol

#catcol = ['dx', 'dt', 'freq', 'shape','method','mvalue','nvalue','npt','npe']
catcol = ['shape','method']

for col in catcol:
    
    dfrec[col] = dfrec[col].astype('category')
    dfsol[col] = dfsol[col].astype('category')

dfrec = dfrec.replace('NC', np.nan)
dfsol = dfsol.replace('NC', np.nan)

dfrec['normrec'] = dfrec['normrec'].astype('float')
dfsol['normsol'] = dfsol['normsol'].astype('float')

vallimit_rec = 3
mask         = (dfrec['normrec'] < vallimit_rec)
dfrec_short  = dfrec[mask]
a = 100*dfrec_short.shape[0]/dfrec.shape[0]
print('Percent Valid in Recnorm: %.3f'%a)

vallimit_sol = 3
mask         = dfsol['normsol'] < vallimit_sol
dfsol_short  = dfsol[mask]
a = 100*dfsol_short.shape[0]/dfsol.shape[0]
print('Percent Valid in Solnorm: %.3f'%a)

conversion_dict1 = {1: 30.0, 2: 15.0, 4: 7.5,8: 3.75}
conversion_dict2 = {1: 0.5, 2: 1.0, 4: 2.0, 6: 3.0}

dfrec_short['dx'] = dfrec_short['dx'].replace(conversion_dict1)
dfrec_short['dt'] = dfrec_short['dt'].replace(conversion_dict2)

dfsol_short['dx'] = dfsol_short['dx'].replace(conversion_dict1)
dfsol_short['dt'] = dfsol_short['dt'].replace(conversion_dict2)
#==============================================================================

#==============================================================================
# Plot Data Sets
#==============================================================================
mask = (dfrec_short['method'] == 'dispte') | (dfrec_short['shape'] == 'crb') | (dfrec_short['freq'] == 1) 

x = 


# Define dimensions
Nx, Ny, Nz = 100, 300, 500
X, Y, Z = np.meshgrid(np.arange(Nx), np.arange(Ny), -np.arange(Nz))

# Create fake data
data = (((X+100)**2 + (Y-20)**2 + 2*Z)/1000+1)

kw = {
    'vmin': data.min(),
    'vmax': data.max(),
    'levels': np.linspace(data.min(), data.max(), 10),
}

# Create a figure with 3D ax
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(111, projection='3d')

# Plot contour surfaces
_ = ax.contourf(
    X[:, :, 0], Y[:, :, 0], data[:, :, 0],
    zdir='z', offset=0, **kw
)
_ = ax.contourf(
    X[0, :, :], data[0, :, :], Z[0, :, :],
    zdir='y', offset=0, **kw
)
C = ax.contourf(
    data[:, -1, :], Y[:, -1, :], Z[:, -1, :],
    zdir='x', offset=X.max(), **kw
)
# --


# Set limits of the plot from coord limits
xmin, xmax = X.min(), X.max()
ymin, ymax = Y.min(), Y.max()
zmin, zmax = Z.min(), Z.max()
ax.set(xlim=[xmin, xmax], ylim=[ymin, ymax], zlim=[zmin, zmax])

# Plot edges
edges_kw = dict(color='0.4', linewidth=1, zorder=1e3)
ax.plot([xmax, xmax], [ymin, ymax], 0, **edges_kw)
ax.plot([xmin, xmax], [ymin, ymin], 0, **edges_kw)
ax.plot([xmax, xmax], [ymin, ymin], [zmin, zmax], **edges_kw)

# Set labels and zticks
ax.set(
    xlabel='X [km]',
    ylabel='Y [km]',
    zlabel='Z [m]',
    zticks=[0, -150, -300, -450],
)

# Set zoom and angle view
ax.view_init(40, -30, 0)
ax.set_box_aspect(None, zoom=0.9)

# Colorbar
fig.colorbar(C, ax=ax, fraction=0.02, pad=0.1, label='Name [units]')

# Show Figure
plt.show()

#==============================================================================