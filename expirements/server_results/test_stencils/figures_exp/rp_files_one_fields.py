#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
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
import matplotlib.ticker       as mticker    
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   matplotlib              import ticker
from   matplotlib              import cm
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

#parameters 17
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


#=============================================================================

#==============================================================================
# Read Txt Files
#==============================================================================
ptype    = 1
dx_ref   = 1
dt_ref   = 1 
freq_ref = 1

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
# Selecting Data
#==============================================================================
lf_select = []

for k1 in range(0,ntr1):
    
    if(test1_results[k1][4]==freq_ref):
    
        lf_select.append(test1_results[k1])   

nlf = len(lf_select)

lfdxdt_select = []

for k1 in range(0,nlf):
    
    if(lf_select[k1][2]==1 and lf_select[k1][3]==1):
        
        lfdxdt_select.append(lf_select[k1])   

nlfdxdt = len(lfdxdt_select)

lfdxdt_select_spatte_cl       = []
lfdxdt_select_spectetheta_cl  = []
lfdxdt_select_dispte_crb      = []
lfdxdt_select_displs_crb      = []
lfdxdt_select_specls_crb      = []

for k1 in range(0,nlfdxdt):
    
    if(lfdxdt_select[k1][5]=='cl'  and lfdxdt_select[k1][6]=='spatte'):      lfdxdt_select_spatte_cl.append(lfdxdt_select[k1])
    if(lfdxdt_select[k1][5]=='cl'  and lfdxdt_select[k1][6]=='spectetheta'): lfdxdt_select_spectetheta_cl.append(lfdxdt_select[k1])    
    if(lfdxdt_select[k1][5]=='crb' and lfdxdt_select[k1][6]=='dispte'):      lfdxdt_select_dispte_crb.append(lfdxdt_select[k1])
    if(lfdxdt_select[k1][5]=='crb' and lfdxdt_select[k1][6]=='displs'):      lfdxdt_select_displs_crb.append(lfdxdt_select[k1])
    if(lfdxdt_select[k1][5]=='crb' and lfdxdt_select[k1][6]=='specls'):      lfdxdt_select_specls_crb.append(lfdxdt_select[k1])
    
ncl         = len(lfdxdt_select_spatte_cl)
ncrb        = len(lfdxdt_select_specls_crb)
nprecselect = len(lfdxdt_select_spatte_cl[0][12][0])
npsolselect = len(lfdxdt_select_spatte_cl[0][13][0])
timesolplot = lfdxdt_select_dispte_crb[0][14] 
#==============================================================================

#==============================================================================
# Find Values
#==============================================================================
mv1 = 1
nv1 = 1

mv2 = 2
nv2 = 1

mv3 = 4
nv3 = 1

mv4 = 6
nv4 = 1

mv5 = 8
nv5 = 1

nlfdxdt_select_spatte_cl = len(lfdxdt_select_spatte_cl)

for k1 in range(0,nlfdxdt_select_spatte_cl):
    
    if(lfdxdt_select_spatte_cl[k1][7]==mv1 and lfdxdt_select_spatte_cl[k1][8]==nv1):
        
        i11 = k1
        i21 = k1
    
    elif(lfdxdt_select_spatte_cl[k1][7]==mv2 and lfdxdt_select_spatte_cl[k1][8]==nv2):
        
        i12 = k1
        i22 = k1

    if(lfdxdt_select_spatte_cl[k1][7]==mv3 and lfdxdt_select_spatte_cl[k1][8]==nv3):
        
        i13 = k1
        i23 = k1

    if(lfdxdt_select_spatte_cl[k1][7]==mv4 and lfdxdt_select_spatte_cl[k1][8]==nv4):
        
        i14 = k1
        i24 = k1

    if(lfdxdt_select_spatte_cl[k1][7]==mv5 and lfdxdt_select_spatte_cl[k1][8]==nv5):
        
        i15 = k1
        i25 = k1

nlfdxdt_select_dispte_crb = len(lfdxdt_select_dispte_crb)

mv1 = 1
nv1 = 1

mv2 = 4
nv2 = 1

mv3 = 4
nv3 = 4

mv4 = 8
nv4 = 1

mv5 = 8
nv5 = 8

for k1 in range(0,nlfdxdt_select_dispte_crb):
    
    if(lfdxdt_select_dispte_crb[k1][7]==mv1 and lfdxdt_select_dispte_crb[k1][8]==nv1):
        
        i31 = k1
        i31 = k1
        
        i41 = k1
        i41 = k1
        
        i51 = k1
        i51 = k1
    
    elif(lfdxdt_select_dispte_crb[k1][7]==mv2 and lfdxdt_select_dispte_crb[k1][8]==nv2):
        
        i32 = k1
        i32 = k1
        
        i42 = k1
        i42 = k1
        
        i52 = k1
        i52 = k1

    if(lfdxdt_select_dispte_crb[k1][7]==mv3 and lfdxdt_select_dispte_crb[k1][8]==nv3):
        
        i33 = k1
        i33 = k1
        
        i43 = k1
        i43 = k1
        
        i53 = k1
        i53 = k1

    if(lfdxdt_select_dispte_crb[k1][7]==mv4 and lfdxdt_select_dispte_crb[k1][8]==nv4):
        
        i34 = k1
        i34 = k1
        
        i44 = k1
        i44 = k1
        
        i54 = k1
        i54 = k1

    if(lfdxdt_select_dispte_crb[k1][7]==mv5 and lfdxdt_select_dispte_crb[k1][8]==nv5):
        
        i35 = k1
        i35 = k1
        
        i45 = k1
        i45 = k1
        
        i55 = k1
        i55 = k1
#==============================================================================

#==============================================================================
# Plot Routines 1
#==============================================================================
def plot1(vsols,vnames,xpos,ypos,extent,vparameters):
    
    plt.figure(figsize = (24,16))
    plt.suptitle('Difference of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
    grid  = plt.GridSpec(6,5,wspace=0.5,hspace=0.5)    
    nfigs = len(xpos)
    vmin  = 0
    vmax  = 0
    
    for k1 in range(0,nfigs):
    
        vmin = min(vmin,np.amin(vsols[k1]))
        vmax = max(vmax,np.amax(vsols[k1]))
    
    factor = 50
    scale  = vmax/factor
        
    for k1 in range(0,nfigs):
    
        plt.subplot(grid[xpos[k1],ypos[k1]])
        
        if(k1==0):
            
            fig1 = plt.imshow(vsols[k1],cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
            plt.grid()
            plt.title('%s'%(vnames[k1]))
            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
            ax = plt.gca()
            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
            cbar.ax.locator_params(nbins=5)
            tick_locator = ticker.MaxNLocator(nbins=5)
            cbar.locator = tick_locator
            cbar.update_ticks()
        
        else:

            fig1 = plt.imshow(vsols[0]-vsols[k1],cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
            plt.grid()
            plt.title('%s'%(vnames[k1]))
            plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
            plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
            ax = plt.gca()
            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
            cbar.ax.locator_params(nbins=5)
            cbar.ax.set_ylabel('[Error]')
            tick_locator = ticker.MaxNLocator(nbins=5)
            cbar.locator = tick_locator
            cbar.update_ticks()
            
    plt.savefig('%s/%s.jpeg'%(locsave,figname),dpi=200,bbox_inches='tight')
    plt.close()

    return
#==============================================================================

#==============================================================================
# Plot Routines 2
#==============================================================================
def plot2(vsols,vnames,xpos,ypos,extent,vparameters):
    
    nfigs  = len(xpos)
    ntimes = vsols[0].shape[0]
    vmax   = 0
    vmin   = 0
    
    for k3 in range(0,ntimes):
        
        for k1 in range(0,nfigs):
        
            vmin = min(vmin,np.amin(vsols[k1][k3]))
            vmax = max(vmax,np.amax(vsols[k1][k3]))
    
    factor = 50
    scale  = vmax/factor
    
    for k3 in range(0,ntimes):
    
        plt.figure(figsize = (24,16))
        plt.suptitle('Difference of Full Displacement at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3][k3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
        grid  = plt.GridSpec(6,5,wspace=0.5,hspace=0.5)    
            
        for k1 in range(0,nfigs):
        
            plt.subplot(grid[xpos[k1],ypos[k1]])
            
            if(k1==0):
                
                fig1 = plt.imshow(vsols[k1][k3],cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                plt.grid()
                plt.title('%s'%(vnames[k1]))
                plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                ax = plt.gca()
                cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                cbar.ax.locator_params(nbins=5)
                tick_locator = ticker.MaxNLocator(nbins=5)
                cbar.locator = tick_locator
                cbar.update_ticks()
            
            else:
    
                fig1 = plt.imshow(vsols[0][k3]-vsols[k1][k3],cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                plt.grid()
                plt.title('%s'%(vnames[k1]))
                plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                ax = plt.gca()
                cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                cbar.ax.locator_params(nbins=5)
                cbar.ax.set_ylabel('[Error]')
                tick_locator = ticker.MaxNLocator(nbins=5)
                cbar.locator = tick_locator
                cbar.update_ticks()
                
        plt.savefig('%s/%s_ntime%d.jpeg'%(locsave,figname,k3),dpi=200,bbox_inches='tight')
        plt.close()

    return
#==============================================================================

#==============================================================================
# Plot Routines 3
#==============================================================================
def plot3(vsols,vnames,xpos,ypos,extent,vparameters,xpositionv,ypositionv,vrectime):
    
    nfigs  = len(xpos)
    ntimes = vsols[0].shape[1]
    vmax   = 0
    vmin   = 0
    
    for k3 in range(0,ntimes):
        
        for k1 in range(0,nfigs):
        
            vmin = min(vmin,np.amin(vsols[k1][:,k3]))
            vmax = max(vmax,np.amax(vsols[k1][:,k3]))
    
    factor = 50
    scale  = vmax/factor
    vmin   = 1.2*vmin
    vmax   = 1.2*vmax
    
    for k3 in range(0,ntimes):
    
        plt.figure(figsize = (24,16))
        plt.suptitle('Diference of Full Receivers at time  %.3f s \n dx = %.4f m - dt = %.4f s - freq = %.3f Hz \n %s \n xpos = %.2f m - ypos = %.2f m'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5],xpositionv[k3],ypositionv[k3]))
        grid  = plt.GridSpec(6,5,wspace=0.5,hspace=0.5)    
            
        for k1 in range(0,nfigs):
        
            plt.subplot(grid[xpos[k1],ypos[k1]])
            
            if(k1==0):
                
                plt.plot(vrectime,vsols[k1][:,k3])
                plt.grid()
                plt.title('Reference')
                plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                plt.ylim((vmin,vmax))
                plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                plt.ylabel('[Error]')
            
            else:
                
                plt.plot(vrectime,vsols[0][:,k3]-vsols[k1][:,k3])
                plt.grid()
                plt.title('%s'%(vnames[k1]))
                plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                plt.ylim((vmin,vmax))
                plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                plt.ylabel('[Error]')
                
        plt.savefig('%s/%s_ntime%d.jpeg'%(locsave,figname,k3),dpi=200,bbox_inches='tight')
        plt.close()

    return
#==============================================================================


#==============================================================================
# Loc Save
#==============================================================================
locsave = 'comp_fig/teste%d/dx%ddt%dfreq%d/'%(ptype,dx_ref,dt_ref,freq_ref) 
#==============================================================================

#==============================================================================
# Plot Infos
#==============================================================================
xpos        = [0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5]
ypos        = [2,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4]
fscale      = 10**(-3)
extent1     = [fscale*lfdxdt_select_spatte_cl[0][15][2],fscale*lfdxdt_select_spatte_cl[0][15][4], fscale*lfdxdt_select_spatte_cl[0][15][9], fscale*lfdxdt_select_spatte_cl[0][15][8]]
extent2     = [fscale*lfdxdt_select_spatte_cl[0][15][2],fscale*lfdxdt_select_spatte_cl[0][15][4], fscale*lfdxdt_select_spatte_cl[0][15][3], fscale*lfdxdt_select_spatte_cl[0][15][5]]
vnames      = ['Reference',
              'spatte-cl M=1','spatte-cl M=2','spatte-cl M=4','spatte-cl M=6','spatte-cl M=8',
              'spectetheta-cl M=1','spectetheta-cl M=2','spectetheta-cl M=4','spectetheta-cl M=6','spectetheta-cl M=8',
              'dispte-crb M=1 and N=1','dispte-crb M=4 and N=1','dispte-crb M=4 and N=4','dispte-crb M=8 and N=1','dispte-crb M=8 and N=8',
              'displs-crb M=1 and N=1','displs-crb M=4 and N=1','displs-crb M=4 and N=4','displs-crb M=8 and N=1','displs-crb M=8 and N=8',
              'specls-crb M=1 and N=1','specls-crb M=4 and N=1','specls-crb M=4 and N=4','specls-crb M=8 and N=1','specls-crb M=8 and N=8']

vrectime    = np.linspace(fscale*lfdxdt_select_spatte_cl[0][15][8],fscale*lfdxdt_select_spatte_cl[0][15][9],lfdxdt_select_spatte_cl[0][17][3].shape[0])
#==============================================================================

#==============================================================================
# Plot1 Execute
#==============================================================================
vsols = [lfdxdt_select_spatte_cl[0][17][1],
          lfdxdt_select_spatte_cl[i11][17][0],lfdxdt_select_spatte_cl[i12][17][0],lfdxdt_select_spatte_cl[i13][17][0],lfdxdt_select_spatte_cl[i14][17][0],lfdxdt_select_spatte_cl[i15][17][0],
          lfdxdt_select_spectetheta_cl[i21][17][0],lfdxdt_select_spectetheta_cl[i22][17][0],lfdxdt_select_spectetheta_cl[i23][17][0],lfdxdt_select_spectetheta_cl[i24][17][0],lfdxdt_select_spectetheta_cl[i25][17][0],
          lfdxdt_select_dispte_crb[i31][17][0],lfdxdt_select_dispte_crb[i32][17][0],lfdxdt_select_dispte_crb[i33][17][0],lfdxdt_select_dispte_crb[i34][17][0],lfdxdt_select_dispte_crb[i35][17][0],
          lfdxdt_select_displs_crb[i41][17][0],lfdxdt_select_displs_crb[i42][17][0],lfdxdt_select_displs_crb[i43][17][0],lfdxdt_select_displs_crb[i44][17][0],lfdxdt_select_displs_crb[i45][17][0],
          lfdxdt_select_specls_crb[i51][17][0],lfdxdt_select_specls_crb[i52][17][0],lfdxdt_select_specls_crb[i53][17][0],lfdxdt_select_specls_crb[i54][17][0],lfdxdt_select_specls_crb[i55][17][0]] 

vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

figname     = 'fieldsrec_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

P1          = plot1(vsols,vnames,xpos,ypos,extent1,vparameters)
#==============================================================================

#==============================================================================
# Plot2 Execute
#==============================================================================
vsols = [lfdxdt_select_spatte_cl[0][17][5],
          lfdxdt_select_spatte_cl[i11][17][4],lfdxdt_select_spatte_cl[i12][17][4],lfdxdt_select_spatte_cl[i13][17][4],lfdxdt_select_spatte_cl[i14][17][4],lfdxdt_select_spatte_cl[i15][17][4],
          lfdxdt_select_spectetheta_cl[i21][17][4],lfdxdt_select_spectetheta_cl[i22][17][4],lfdxdt_select_spectetheta_cl[i23][17][4],lfdxdt_select_spectetheta_cl[i24][17][4],lfdxdt_select_spectetheta_cl[i25][17][4],
          lfdxdt_select_dispte_crb[i31][17][4],lfdxdt_select_dispte_crb[i32][17][4],lfdxdt_select_dispte_crb[i33][17][4],lfdxdt_select_dispte_crb[i34][17][4],lfdxdt_select_dispte_crb[i35][17][4],
          lfdxdt_select_displs_crb[i41][17][4],lfdxdt_select_displs_crb[i42][17][4],lfdxdt_select_displs_crb[i43][17][4],lfdxdt_select_displs_crb[i44][17][4],lfdxdt_select_displs_crb[i45][17][4],
          lfdxdt_select_specls_crb[i51][17][4],lfdxdt_select_specls_crb[i52][17][4],lfdxdt_select_specls_crb[i53][17][4],lfdxdt_select_specls_crb[i54][17][4],lfdxdt_select_specls_crb[i55][17][4]] 


vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                np.array(timesolplot)/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

figname     = 'fieldssolplot_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

P2          = plot2(vsols,vnames,xpos,ypos,extent1,vparameters)
#==============================================================================

#==============================================================================
# Plot3 Execute
#==============================================================================
vsols = [lfdxdt_select_spatte_cl[0][17][3],
          lfdxdt_select_spatte_cl[i11][17][2],lfdxdt_select_spatte_cl[i12][17][2],lfdxdt_select_spatte_cl[i13][17][2],lfdxdt_select_spatte_cl[i14][17][2],lfdxdt_select_spatte_cl[i15][17][2],
          lfdxdt_select_spectetheta_cl[i21][17][2],lfdxdt_select_spectetheta_cl[i22][17][2],lfdxdt_select_spectetheta_cl[i23][17][2],lfdxdt_select_spectetheta_cl[i24][17][2],lfdxdt_select_spectetheta_cl[i25][17][2],
          lfdxdt_select_dispte_crb[i31][17][2],lfdxdt_select_dispte_crb[i32][17][2],lfdxdt_select_dispte_crb[i33][17][2],lfdxdt_select_dispte_crb[i34][17][2],lfdxdt_select_dispte_crb[i35][17][2],
          lfdxdt_select_displs_crb[i41][17][2],lfdxdt_select_displs_crb[i42][17][2],lfdxdt_select_displs_crb[i43][17][2],lfdxdt_select_displs_crb[i44][17][2],lfdxdt_select_displs_crb[i45][17][2],
          lfdxdt_select_specls_crb[i51][17][2],lfdxdt_select_specls_crb[i52][17][2],lfdxdt_select_specls_crb[i53][17][2],lfdxdt_select_specls_crb[i54][17][2],lfdxdt_select_specls_crb[i55][17][2]] 

vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

figname     = 'fieldsrecselect_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

P3          = plot3(vsols,vnames,xpos,ypos,extent1,vparameters,xpositionv,ypositionv,vrectime)
#==============================================================================