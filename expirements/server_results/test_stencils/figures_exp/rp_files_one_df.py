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
            
        test_results  = pickle.load(f) 
        ntr           = len(test_results)
        testname      = 'Homogeneos Velocity Model'
        xpositionv    = np.array([750.0,2250.0, 750.0,2250.0])
        ypositionv    = np.array([750.0, 750.0,2250.0,2250.0])
            
    elif(ptype==2): 
            
        test_results  = pickle.load(f) 
        ntr           = len(test_results)
        testname      = 'Heterogeneos Velocity Model'
        xpositionv    = np.array([500.0,1500.0, 500.0,1500.0])
        ypositionv    = np.array([500.0, 500.0,1500.0,1500.0])
        
    elif(ptype==3): 
            
        test_results  = pickle.load(f) 
        ntr           = len(test_results)
        testname      = 'SEG/EAGE 2D Salt Velocity Model'
        xpositionv    = np.array([4000.0,4000.0,4000.0,6000.0,6000.0,6000.0,8000.0,8000.0,8000.0])   
        ypositionv    = np.array([2000.0,2500.0,3000.0,2000.0,2500.0,3000.0,2000.0,2500.0,3000.0]) 

    elif(ptype==4): 
            
        test_results  = pickle.load(f)
        ntr           = len(test_results)
        testname      = 'Marmousi Velocity Model'
        xpositionv  = np.array([6000.0,6000.0,6000.0,8000.0,8000.0,8000.0,10000.0,10000.0,10000.0,12000.0,12000.0,12000.0])
        ypositionv  = np.array([1000.0,2000.0,3000.0,1000.0,2000.0,3000.0,1000.0,2000.0,3000.0,1000.0,2000.0,3000.0])
#==============================================================================

#==============================================================================
# Data Frame Construction
#==============================================================================
ntest_results = len(test_results)

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
vcflv    = []

for k1 in range(0,ntest_results):
    
    vdx.append(test_results[k1][2])
    vdt.append(test_results[k1][3])
    vfreq.append(test_results[k1][4])
    vshape.append(test_results[k1][5])
    vmethod.append(test_results[k1][6])
    vmvalue.append(test_results[k1][7])
    vnvalue.append(test_results[k1][8])
    vnpt.append(test_results[k1][9])
    vnpe.append(test_results[k1][10])
    vnormrec.append(test_results[k1][11][4])
    vnormsol.append(test_results[k1][13][4][-1])
    vdxv.append(test_results[k1][15][6])
    vdtv.append(test_results[k1][15][13])
    vfreqv.append(1000*test_results[k1][15][10])
    vcflv.append(test_results[k1][15][11])

dfrec = pd.DataFrame()
dfsol = pd.DataFrame()

dfrec['dx']      = vdx
dfrec['dt']      = vdt
dfrec['freq']    = vfreq
dfrec['freqv']   = vfreqv
dfrec['cfl']     = vcflv
dfrec['shape']   = vshape
dfrec['method']  = vmethod
dfrec['mvalue']  = vmvalue
dfrec['nvalue']  = vnvalue
dfrec['npt']     = vnpt
dfrec['npe']     = vnpe
dfrec['normrec'] = vnormrec

dfsol['dx']      = vdx
dfsol['dt']      = vdt
dfsol['freq']    = vfreq
dfsol['freqv']   = vfreqv
dfsol['cfl']     = vcflv
dfsol['shape']   = vshape
dfsol['method']  = vmethod
dfsol['mvalue']  = vmvalue
dfsol['nvalue']  = vnvalue
dfsol['npt']     = vnpt
dfsol['npe']     = vnpe
dfsol['normsol'] = vnormsol

if(ptype==1):

    conversion_dict1  = {1: 30.0, 2: 15.0, 4: 7.5,8: 3.75}
    conversion_dict2  = {1: 0.5, 2: 1.0, 4: 2.0, 6: 3.0}

dfrec['dx'] = dfrec['dx'].replace(conversion_dict1)
dfrec['dt'] = dfrec['dt'].replace(conversion_dict2)
dfsol['dx'] = dfsol['dx'].replace(conversion_dict1)
dfsol['dt'] = dfsol['dt'].replace(conversion_dict2)

catcol = ['shape','method']

for col in catcol:
    
    dfrec[col] = dfrec[col].astype('category')
    dfsol[col] = dfsol[col].astype('category')

dfrec = dfrec.replace('NC', np.nan)
dfsol = dfsol.replace('NC', np.nan)

dfrec['normrec'] = dfrec['normrec'].astype('float')
dfsol['normsol'] = dfsol['normsol'].astype('float')
#==============================================================================

#==============================================================================
# Norm Check
#==============================================================================
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
#==============================================================================

#==============================================================================
# Plot Data Sets
#==============================================================================
# mask =   (dfrec_short['method'] == 'displs') \
#         & (dfrec_short['shape']  == 'crb') \
#         & (dfrec_short['dt']     == 0.5) \
#         & (dfrec_short['freq']   == 1)  

mask =   (dfrec_short['method'] == 'displs') \
        & (dfrec_short['shape']  == 'crb') \
        & (dfrec_short['freq']   == 1) 

dfrec_short_select  = dfrec_short[mask].reset_index(drop=True)
xv = dfrec_short_select['mvalue'].values
yv = dfrec_short_select['nvalue'].values
zv = dfrec_short_select['dx'].values
cv = dfrec_short_select['normrec'].values
sv = 400*dfrec_short_select['dt'].values

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
data_plot = ax.scatter(xv, yv, zv, c=cv, s=sv,cmap="jet",alpha=0.5)
cbar = fig.colorbar(data_plot,pad=0.2)
cbar.ax.set_ylabel('[Error]')
ax.set_xlabel('mvalue')
ax.set_ylabel('nvalue')
ax.set_zlabel('dx[m]')
ax.set_title('ssssss')

# mask =   (dfrec_short['method'] == 'displs') \
#         & (dfrec_short['shape']  == 'crb') \
#         & (dfrec_short['freq']   == 1) 

# dfrec_short_select  = dfrec_short[mask].reset_index(drop=True)
# xv = dfrec_short_select['mvalue'].values
# yv = dfrec_short_select['nvalue'].values
# zv = dfrec_short_select['dx'].values
# cv = dfrec_short_select['normrec'].values
# sv = 500*dfrec_short_select['dt'].values

# mask2 =   (dfrec_short['method'] == 'displs') \
#         & (dfrec_short['shape']  == 'crb') \
#         & (dfrec_short['freq']   == 2) 

# dfrec_short_select2  = dfrec_short[mask2].reset_index(drop=True)
# xv2 = dfrec_short_select2['mvalue'].values
# yv2 = dfrec_short_select2['nvalue'].values
# zv2 = dfrec_short_select2['dx'].values
# cv2 = dfrec_short_select2['normrec'].values
# sv2 = 500*dfrec_short_select2['dt'].values

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# data_plot = ax.scatter(xv, yv, zv, c=cv, s=sv,cmap="jet")
# data_plot = ax.scatter(xv2, yv2, zv2, c=cv2, s=sv2,cmap="jet",marker='d')

# cbar = fig.colorbar(data_plot,pad=0.2)
# cbar.ax.set_ylabel('[Error]')
# ax.set_xlabel('mvalue')
# ax.set_ylabel('nvalue')
# ax.set_zlabel('dx[m]')
# ax.set_title('ssssss')
#==============================================================================