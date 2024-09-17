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

lfdxdt_normrecrel1_spatte_cl                = np.zeros((ncl,1))
lfdxdt_normrecrel1_spectetheta_cl           = np.zeros((ncl,1))
lfdxdt_normrecrel1_dispte_crb               = np.zeros((ncl,ncl))
lfdxdt_normrecrel1_displs_crb               = np.zeros((ncl,ncl))
lfdxdt_normrecrel1_specls_crb               = np.zeros((ncl,ncl))

lfdxdt_normrecrel2_spatte_cl                = np.zeros((ncl,1))
lfdxdt_normrecrel2_spectetheta_cl           = np.zeros((ncl,1))
lfdxdt_normrecrel2_dispte_crb               = np.zeros((ncl,ncl))
lfdxdt_normrecrel2_displs_crb               = np.zeros((ncl,ncl))
lfdxdt_normrecrel2_specls_crb               = np.zeros((ncl,ncl))

lfdxdt_normrecselectrel1_spatte_cl          = np.zeros((nprecselect,ncl,1))
lfdxdt_normrecselectrel1_spectetheta_cl     = np.zeros((nprecselect,ncl,1))
lfdxdt_normrecselectrel1_dispte_crb         = np.zeros((nprecselect,ncl,ncl))
lfdxdt_normrecselectrel1_displs_crb         = np.zeros((nprecselect,ncl,ncl))
lfdxdt_normrecselectrel1_specls_crb         = np.zeros((nprecselect,ncl,ncl))

lfdxdt_normrecselectrel2_spatte_cl          = np.zeros((nprecselect,ncl,1))
lfdxdt_normrecselectrel2_spectetheta_cl     = np.zeros((nprecselect,ncl,1))
lfdxdt_normrecselectrel2_dispte_crb         = np.zeros((nprecselect,ncl,ncl))
lfdxdt_normrecselectrel2_displs_crb         = np.zeros((nprecselect,ncl,ncl))
lfdxdt_normrecselectrel2_specls_crb         = np.zeros((nprecselect,ncl,ncl))

lfdxdt_normsolplotrel1_spatte_cl            = np.zeros((npsolselect,ncl,1))
lfdxdt_normsolplotrel1_spectetheta_cl       = np.zeros((npsolselect,ncl,1))
lfdxdt_normsolplotrel1_dispte_crb           = np.zeros((npsolselect,ncl,ncl))
lfdxdt_normsolplotrel1_displs_crb           = np.zeros((npsolselect,ncl,ncl))
lfdxdt_normsolplotrel1_specls_crb           = np.zeros((npsolselect,ncl,ncl))

lfdxdt_normsolplotrel2_spatte_cl            = np.zeros((npsolselect,ncl,1))
lfdxdt_normsolplotrel2_spectetheta_cl       = np.zeros((npsolselect,ncl,1))
lfdxdt_normsolplotrel2_dispte_crb           = np.zeros((npsolselect,ncl,ncl))
lfdxdt_normsolplotrel2_displs_crb           = np.zeros((npsolselect,ncl,ncl))
lfdxdt_normsolplotrel2_specls_crb           = np.zeros((npsolselect,ncl,ncl))

for k1 in range(0,ncl):
    
    lfdxdt_normrecrel1_spatte_cl[k1,0] = lfdxdt_select_spatte_cl[k1][11][3]
    lfdxdt_normrecrel2_spatte_cl[k1,0] = lfdxdt_select_spatte_cl[k1][11][4]
    
    lfdxdt_normrecrel1_spectetheta_cl[k1,0] = lfdxdt_select_spectetheta_cl[k1][11][3]
    lfdxdt_normrecrel2_spectetheta_cl[k1,0] = lfdxdt_select_spectetheta_cl[k1][11][4]
    
    for k2 in range(0,nprecselect):
        
        lfdxdt_normrecselectrel1_spatte_cl[k2,k1,0] = lfdxdt_select_spatte_cl[k1][12][3][k2]
        lfdxdt_normrecselectrel2_spatte_cl[k2,k1,0] = lfdxdt_select_spatte_cl[k1][12][4][k2]
        
        lfdxdt_normrecselectrel1_spectetheta_cl[k2,k1,0] = lfdxdt_select_spectetheta_cl[k1][12][3][k2]
        lfdxdt_normrecselectrel2_spectetheta_cl[k2,k1,0] = lfdxdt_select_spectetheta_cl[k1][12][4][k2]
        
    for k3 in range(0,npsolselect):
        
        lfdxdt_normsolplotrel1_spatte_cl[k3,k1,0] = lfdxdt_select_spatte_cl[k1][13][3][k3]
        lfdxdt_normsolplotrel2_spatte_cl[k3,k1,0] = lfdxdt_select_spatte_cl[k1][13][4][k3]
    
        lfdxdt_normsolplotrel1_spectetheta_cl[k3,k1,0] = lfdxdt_select_spectetheta_cl[k1][13][3][k3]
        lfdxdt_normsolplotrel2_spectetheta_cl[k3,k1,0] = lfdxdt_select_spectetheta_cl[k1][13][4][k3]

contglob        = 0
contloc         = 0
npte_dispte_crb = []
npte_displs_crb = []
npte_specls_crb = []

for k1 in range(0,ncl):

    for k2 in range(0,contloc+1):

        npte_dispte_crb.append([lfdxdt_select_dispte_crb[contglob][7],lfdxdt_select_dispte_crb[contglob][8],lfdxdt_select_dispte_crb[contglob][9],lfdxdt_select_dispte_crb[contglob][10]])
        npte_displs_crb.append([lfdxdt_select_displs_crb[contglob][7],lfdxdt_select_displs_crb[contglob][8],lfdxdt_select_displs_crb[contglob][9],lfdxdt_select_displs_crb[contglob][10]])
        npte_specls_crb.append([lfdxdt_select_specls_crb[contglob][7],lfdxdt_select_specls_crb[contglob][8],lfdxdt_select_specls_crb[contglob][9],lfdxdt_select_specls_crb[contglob][10]])
        
        lfdxdt_normrecrel1_dispte_crb[k2,k1] = lfdxdt_select_dispte_crb[contglob][11][3]
        lfdxdt_normrecrel2_dispte_crb[k2,k1] = lfdxdt_select_dispte_crb[contglob][11][4]
        
        lfdxdt_normrecrel1_displs_crb[k2,k1] = lfdxdt_select_displs_crb[contglob][11][3]
        lfdxdt_normrecrel2_displs_crb[k2,k1] = lfdxdt_select_displs_crb[contglob][11][4]
        
        lfdxdt_normrecrel1_specls_crb[k2,k1] = lfdxdt_select_specls_crb[contglob][11][3]
        lfdxdt_normrecrel2_specls_crb[k2,k1] = lfdxdt_select_specls_crb[contglob][11][4]
                
        for k3 in range(0,nprecselect):
            
            lfdxdt_normrecselectrel1_dispte_crb[k3,k2,k1] = lfdxdt_select_dispte_crb[contglob][12][3][k3]
            lfdxdt_normrecselectrel2_dispte_crb[k3,k2,k1] = lfdxdt_select_dispte_crb[contglob][12][4][k3]
            
            lfdxdt_normrecselectrel1_displs_crb[k3,k2,k1] = lfdxdt_select_displs_crb[contglob][12][3][k3]
            lfdxdt_normrecselectrel2_displs_crb[k3,k2,k1] = lfdxdt_select_displs_crb[contglob][12][4][k3]
            
            lfdxdt_normrecselectrel1_specls_crb[k3,k2,k1] = lfdxdt_select_specls_crb[contglob][12][3][k3]
            lfdxdt_normrecselectrel2_specls_crb[k3,k2,k1] = lfdxdt_select_specls_crb[contglob][12][4][k3]
            
        for k4 in range(0,npsolselect):
          
            lfdxdt_normsolplotrel1_dispte_crb[k4,k2,k1] = lfdxdt_select_dispte_crb[contglob][13][3][k4]
            lfdxdt_normsolplotrel2_dispte_crb[k4,k2,k1] = lfdxdt_select_dispte_crb[contglob][13][4][k4]

            lfdxdt_normsolplotrel1_displs_crb[k4,k2,k1] = lfdxdt_select_displs_crb[contglob][13][3][k4]
            lfdxdt_normsolplotrel2_displs_crb[k4,k2,k1] = lfdxdt_select_displs_crb[contglob][13][4][k4]

            lfdxdt_normsolplotrel1_specls_crb[k4,k2,k1] = lfdxdt_select_specls_crb[contglob][13][3][k4]
            lfdxdt_normsolplotrel2_specls_crb[k4,k2,k1] = lfdxdt_select_specls_crb[contglob][13][4][k4]
        
        contglob = contglob + 1

    contloc = contloc+1
    
nnpte_dispte_crb = len(npte_dispte_crb)
vnptt            = []
vnpte            = []

for k1 in range(0,nnpte_dispte_crb):
 
    vnptt.append(npte_dispte_crb[k1][2])
    vnpte.append(npte_dispte_crb[k1][3])
    
vnptt = sorted(set(vnptt))
vnpte = sorted(set(vnpte))

timesolplot  = lfdxdt_select_dispte_crb[0][14]
ordersv      = [2*i for i in range(0,9)]
vticks       = ['s', '+', '+', '+',  '+',  '^',   '^',   'D',      'D','s']
vline        = ['-', '-', '-', '--', '-.', '--',  '-.',  '--',     '-.','-']
vcolors      = ['b', 'g', 'r', 'c',  'm',  'y',   'b',   'purple', 'teal','lime']
linep        = [vticks,vline,vcolors]
#==============================================================================

#==============================================================================
# Plot Routines 1
#==============================================================================
def plot1(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep):

    nvcl     = len(vcl)
    nordersv = len(ordersv)
    nvnpte   = len(vnpte)
    
    plt.figure(figsize = (16,12))
    plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
    grid = plt.GridSpec(2,2,wspace=0.4,hspace=0.2)    
    plt.subplot(grid[xpos[0],ypos[0]])
    min_value = 0.8*min(min(vcl[0]),min(vcl[1]),min(vcl[2]),min(vcl[3]),min(vcl[4]))
    max_value = 1.2*max(max(vcl[0]),max(vcl[1]),max(vcl[2]),max(vcl[3]),max(vcl[4])) 
    
    for k1 in range(0,nvcl):
        
        plt.plot(ordersv[1::],vcl[k1],label=clnames[k1],color=linep[2][k1],linestyle=linep[1][k1],marker=linep[0][k1])
    
    plt.grid()
    plt.title('cross-line')
    plt.xticks(ordersv)
    plt.legend()
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    plt.ylim((min_value,max_value))
    plt.xlabel('[Order]')
    plt.ylabel('[Error]')
    plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))
    ax = plt.gca()
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter(ordersv))
    ax.xaxis.set_major_locator(plt.MaxNLocator(nordersv+1))
        
    for k1 in range(1,4):
    
        plt.subplot(grid[xpos[k1],ypos[k1]])
        
        vcrb[k1-1][vcrb[k1-1]==0] = np.nan 
        
        fig1 = plt.imshow(np.transpose(vcrb[k1-1]),cmap='jet',interpolation='kaiser')  
        plt.grid()
        plt.title('%s'%(crbnames[k1-1]))
        plt.xlabel('[Number of Extra Points]')
        plt.ylabel('[Order]')
        ax = plt.gca()
        ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte)]))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte))
        ax.xaxis.set_major_locator(plt.MaxNLocator(len(vnpte)))
        ax.yaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
        ax.yaxis.set_major_formatter(ticker.FixedFormatter(ordersv[1::]))
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
def plot2(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep):

    nvcl     = len(vcl)
    nordersv = len(ordersv)
    nvnpte   = len(vnpte)
    ntimes   = vcl[0].shape[0]
    
    for k3 in range(0,ntimes):
    
        plt.figure(figsize = (16,12))
        plt.suptitle('Relative Quadratic Error of Full Displacements at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3][k3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
        grid = plt.GridSpec(2,2,wspace=0.4,hspace=0.2)    
        plt.subplot(grid[xpos[0],ypos[0]])
        min_value = 0.8*min(min(vcl[0][k3,:]),min(vcl[1][k3,:]),min(vcl[2][k3,:]),min(vcl[3][k3,:]),min(vcl[4][k3,:]))
        max_value = 1.2*max(max(vcl[0][k3,:]),max(vcl[1][k3,:]),max(vcl[2][k3,:]),max(vcl[3][k3,:]),max(vcl[4][k3,:])) 
        
        for k1 in range(0,nvcl):
            
            plt.plot(ordersv[1::],vcl[k1][k3,:],label=clnames[k1],color=linep[2][k1],linestyle=linep[1][k1],marker=linep[0][k1])
        
        plt.grid()
        plt.title('cross-line')
        plt.xticks(ordersv)
        plt.legend()
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.ylim((min_value,max_value))
        plt.xlabel('[Order]')
        plt.ylabel('[Error]')
        plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))
        ax = plt.gca()
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(ordersv))
        ax.xaxis.set_major_locator(plt.MaxNLocator(nordersv+1))
            
        for k1 in range(1,4):
        
            plt.subplot(grid[xpos[k1],ypos[k1]])
            
            vcrb[k1-1][vcrb[k1-1]==0] = np.nan 
            
            fig1 = plt.imshow(np.transpose(vcrb[k1-1][k3,:]),cmap='jet',interpolation='kaiser')  
            plt.grid()
            plt.title('%s'%(crbnames[k1-1]))
            plt.xlabel('[Number of Extra Points]')
            plt.ylabel('[Order]')
            ax = plt.gca()
            ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte)]))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte))
            ax.xaxis.set_major_locator(plt.MaxNLocator(len(vnpte)))
            ax.yaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
            ax.yaxis.set_major_formatter(ticker.FixedFormatter(ordersv[1::]))
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
# Plot Routines 2
#==============================================================================
def plot3(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,xpositionv,ypositionv):

    nvcl     = len(vcl)
    nordersv = len(ordersv)
    nvnpte   = len(vnpte)
    ntimes   = vcl[0].shape[0]
    
    for k3 in range(0,ntimes):
    
        plt.figure(figsize = (16,12))
        plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4f m - dt = %.4f s - freq = %.3f Hz \n %s \n xpos = %.2f m - ypos = %.2f m'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5],xpositionv[k3],ypositionv[k3]))

        grid = plt.GridSpec(2,2,wspace=0.4,hspace=0.2)    
        plt.subplot(grid[xpos[0],ypos[0]])
        min_value = 0.8*min(min(vcl[0][k3,:]),min(vcl[1][k3,:]),min(vcl[2][k3,:]),min(vcl[3][k3,:]),min(vcl[4][k3,:]))
        max_value = 1.2*max(max(vcl[0][k3,:]),max(vcl[1][k3,:]),max(vcl[2][k3,:]),max(vcl[3][k3,:]),max(vcl[4][k3,:])) 
        
        for k1 in range(0,nvcl):
            
            plt.plot(ordersv[1::],vcl[k1][k3,:],label=clnames[k1],color=linep[2][k1],linestyle=linep[1][k1],marker=linep[0][k1])
        
        plt.grid()
        plt.title('cross-line')
        plt.xticks(ordersv)
        plt.legend()
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        plt.ylim((min_value,max_value))
        plt.xlabel('[Order]')
        plt.ylabel('[Error]')
        plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.2e'))
        ax = plt.gca()
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
        ax.xaxis.set_major_formatter(ticker.FixedFormatter(ordersv))
        ax.xaxis.set_major_locator(plt.MaxNLocator(nordersv+1))
            
        for k1 in range(1,4):
        
            plt.subplot(grid[xpos[k1],ypos[k1]])
            
            vcrb[k1-1][vcrb[k1-1]==0] = np.nan 
            
            fig1 = plt.imshow(np.transpose(vcrb[k1-1][k3,:]),cmap='jet',interpolation='kaiser')  
            plt.grid()
            plt.title('%s'%(crbnames[k1-1]))
            plt.xlabel('[Number of Extra Points]')
            plt.ylabel('[Order]')
            ax = plt.gca()
            ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte)]))
            ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte))
            ax.xaxis.set_major_locator(plt.MaxNLocator(len(vnpte)))
            ax.yaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
            ax.yaxis.set_major_formatter(ticker.FixedFormatter(ordersv[1::]))
            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
            cbar.ax.locator_params(nbins=5)
            cbar.ax.set_ylabel('[Error]')
            tick_locator = ticker.MaxNLocator(nbins=5)
            cbar.locator = tick_locator
            cbar.update_ticks()
    
        plt.savefig('%s/%s_npos%d.jpeg'%(locsave,figname,k3),dpi=200,bbox_inches='tight')
        
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
xpos        = [0,0,1,1]
ypos        = [0,1,0,1]
clnames     = ['spatte','spectetheta','dispte-crb N=1','displs-crb N=1','displs-crb N=1','specls-crb N=1']
crbnames    = ['dispte-crb','displs-crb','specls-crb']
#==============================================================================

#==============================================================================
# Plot1 Execute
#==============================================================================
vcl         = [lfdxdt_normrecrel2_spatte_cl,lfdxdt_normrecrel2_spectetheta_cl,lfdxdt_normrecrel2_dispte_crb[0,:],
                lfdxdt_normrecrel2_displs_crb[0,:],lfdxdt_normrecrel2_specls_crb[0,:]]

vcrb        = [lfdxdt_normrecrel2_dispte_crb,lfdxdt_normrecrel2_displs_crb,lfdxdt_normrecrel2_specls_crb]

vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

figname     = 'normrec_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

P1          = plot1(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep)
#==============================================================================

#==============================================================================
# Plot2 Execute
#==============================================================================
vcl         = [lfdxdt_normsolplotrel2_spatte_cl,lfdxdt_normsolplotrel2_spectetheta_cl,lfdxdt_normsolplotrel2_dispte_crb[:,0,:],
                lfdxdt_normsolplotrel2_displs_crb[:,0,:],lfdxdt_normsolplotrel2_specls_crb[:,0,:]]

vcrb        = [lfdxdt_normsolplotrel2_dispte_crb,lfdxdt_normsolplotrel2_displs_crb,lfdxdt_normsolplotrel2_specls_crb]

vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                np.array(timesolplot)/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

figname     = 'normsolplot_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

P2          = plot2(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep)
#==============================================================================

#==============================================================================
# Plot3 Execute
#==============================================================================
vcl         = [lfdxdt_normrecselectrel2_spatte_cl,lfdxdt_normrecselectrel2_spectetheta_cl,lfdxdt_normrecselectrel2_dispte_crb[:,0,:],
                lfdxdt_normrecselectrel2_displs_crb[:,0,:],lfdxdt_normrecselectrel2_specls_crb[:,0,:]]

vcrb        = [lfdxdt_normrecselectrel2_dispte_crb,lfdxdt_normrecselectrel2_displs_crb,lfdxdt_normrecselectrel2_specls_crb]

vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

figname     = 'normrecselect_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

P3          = plot3(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,xpositionv,ypositionv)
#==============================================================================