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

#==============================================================================
# Read Txt Files
#==============================================================================
<<<<<<< HEAD
ptype        = 1
locopen      = '../testresults/test%d_results_norms_fields'%(ptype)
conversion_dict2  = {1: 0.5, 2: 1.0, 4: 2.0, 6: 3.0}

=======
ptype        = 4
locopen      = '../testresults/test%d_results_norms_fields'%(ptype)
conversion_dict2  = {1: 0.5, 2: 1.0, 4: 2.0, 6: 3.0}



>>>>>>> origin/main
with open(locopen, 'rb') as f: 

    if(ptype==1): 
            
        test_results     = pickle.load(f) 
        ntr              = len(test_results)
        testname         = 'Homogeneos Velocity Model'
        solpos           = -1
        conversion_dict1 = {1: 30.0, 2: 15.0, 4: 7.5,8: 3.75}
        dxselect         = 30.0
        dtselect         = 0.5
        nvalueselect     = 1 
        vparameters      = [0.9,dxselect,dtselect]
        
        if(dxselect==30.0): 
            
            dx_ref = 1
        
        elif(dxselect==15.0): 

            dx_ref = 2
            
        elif(dxselect==7.5): 

            dx_ref = 4
        
        elif(dxselect==3.75): 

            dx_ref = 8

    elif(ptype==2): 
            
        test_results   = pickle.load(f) 
        ntr              = len(test_results)
        testname         = 'Heterogeneos Velocity Model'
        solpos           = 3
        conversion_dict1 = {1: 20.0, 2: 10.0, 4: 5.0,8: 2.5}
        dxselect         = 20.0
        dtselect         = 0.5
        nvalueselect     = 1 
        vparameters      = [0.540,dxselect,dtselect]
        
        if(dxselect==20.0): 
            
            dx_ref = 1
        
        elif(dxselect==10.0): 

            dx_ref = 2
            
        elif(dxselect==5.0): 

            dx_ref = 4
        
        elif(dxselect==2.5): 

            dx_ref = 8

    elif(ptype==3): 
            
        test_results     = pickle.load(f) 
        ntr              = len(test_results)
        testname         = 'SEG/EAGE 2D Salt Velocity Model'
        solpos           = 2
        conversion_dict1 = {1: 40.0, 2: 20.0, 4: 10.0,8: 5.0}
        dxselect         = 40.0
        dtselect         = 0.5
        nvalueselect     = 1 
        vparameters      = [0.300,dxselect,dtselect]
        
        if(dxselect==40.0): 
            
            dx_ref = 1
        
        elif(dxselect==20.0): 

            dx_ref = 2
            
        elif(dxselect==10.0): 

            dx_ref = 4
        
        elif(dxselect==5.0): 
            
            dx_ref = 8

    elif(ptype==4): 
            
        test_results     = pickle.load(f)
        ntr              = len(test_results)
        testname         = 'Marmousi Velocity Model'
        solpos           = 2
        conversion_dict1 = {1: 40.0, 2: 20.0, 4: 10.0,8: 5.0}
        dxselect         = 40.0
        dtselect         = 0.5
        nvalueselect     = 1 
        vparameters      = [0.300,dxselect,dtselect]
        
        if(dxselect==40.0): 
            
            dx_ref = 1
        
        elif(dxselect==20.0): 

            dx_ref = 2
            
        elif(dxselect==10.0): 

            dx_ref = 4
        
        elif(dxselect==5.0): 
            
            dx_ref = 8

if(dtselect==0.5): 
            
    dt_ref = 1
        
elif(dtselect==1.0): 

    dt_ref = 2
            
elif(dtselect==2.0): 

    dt_ref = 4
        
elif(dtselect==3.0): 
            
    dt_ref = 6
#==============================================================================

#==============================================================================
# Data Frame Construction
#==============================================================================
ntest_results = len(test_results)

vdx         = []
vdt         = []
vfreq       = []
vshape      = []
vmethod     = []
vmvalue     = []
vnvalue     = []
vnpt        = []
vnpe        = []
vnormrecl2  = []
vnormrecmax = []
vnormsoll2  = []
vnormsolmax = []
vdxv        = []
vdtv        = []
vfreqv      = []
vcflv       = []

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
    vnormrecl2.append(test_results[k1][11][4])
    vnormrecmax.append(test_results[k1][11][5])
    vnormsoll2.append(test_results[k1][13][4][solpos])
    vnormsolmax.append(test_results[k1][13][5][solpos])
    vdxv.append(test_results[k1][15][6])
    vdtv.append(test_results[k1][15][13])
    vfreqv.append(1000*test_results[k1][15][10])
    vcflv.append(test_results[k1][15][11])

dfcomp = pd.DataFrame()

dfcomp['dx']         = vdx
dfcomp['dt']         = vdt
dfcomp['freq']       = vfreq
dfcomp['freqv']      = vfreqv
dfcomp['cfl']        = vcflv
dfcomp['shape']      = vshape
dfcomp['method']     = vmethod
dfcomp['mvalue']     = vmvalue
dfcomp['nvalue']     = vnvalue
dfcomp['npt']        = vnpt
dfcomp['npe']        = vnpe
dfcomp['normrecl2']  = vnormrecl2
dfcomp['normrecmax'] = vnormrecmax
dfcomp['normsoll2']  = vnormsoll2
dfcomp['normsolmax'] = vnormsolmax

dfcomp['dx'] = dfcomp['dx'].replace(conversion_dict1)
dfcomp['dt'] = dfcomp['dt'].replace(conversion_dict2)

catcol = ['shape','method']

for col in catcol:
    
    dfcomp[col] = dfcomp[col].astype('category')

dfcomp = dfcomp.replace('NC', np.nan)
dfcomp = dfcomp.replace('NC', np.nan)

dfcomp['normrecl2']  = dfcomp['normrecl2'].astype('float')
dfcomp['normrecmax'] = dfcomp['normrecmax'].astype('float')

dfcomp['normsoll2']  = dfcomp['normsoll2'].astype('float')
dfcomp['normsolmax'] = dfcomp['normsolmax'].astype('float')
#==============================================================================

#==============================================================================
# Norm Check
#==============================================================================
vallimit_rec  = 3

mask          = (dfcomp['normrecl2'] < vallimit_rec)
dfcomp_short  = dfcomp[mask]
a = 100*dfcomp_short.shape[0]/dfcomp.shape[0]
print('Percent Valid in Recnorml2: %.3f'%a)

mask          = (dfcomp['normrecmax'] < vallimit_rec)
dfcomp_short  = dfcomp[mask]
a = 100*dfcomp_short.shape[0]/dfcomp.shape[0]
print('Percent Valid in Recnormmax: %.3f'%a)

vallimit_sol  = 3

mask          = dfcomp['normsoll2'] < vallimit_sol
dfcomp_short  = dfcomp[mask]
a = 100*dfcomp_short.shape[0]/dfcomp.shape[0]
print('Percent Valid in Solnorml2: %.3f'%a)

mask          = dfcomp['normsolmax'] < vallimit_sol
dfcomp_short  = dfcomp[mask]
a = 100*dfcomp_short.shape[0]/dfcomp.shape[0]
print('Percent Valid in Solnormmax: %.3f'%a)
<<<<<<< HEAD

dfcomp_short  = dfcomp.copy()

=======
>>>>>>> origin/main
#==============================================================================

#==============================================================================
# Select Info
#==============================================================================
mask =   (dfcomp_short['dx']     == dxselect) \
       & (dfcomp_short['dt']     == dtselect) \
       & (dfcomp_short['nvalue'] == nvalueselect)
<<<<<<< HEAD

       
dfselect  = dfcomp_short[mask]
lmethod   = ['spatte','spectetheta','dispte','displs','specls']
lshape    = ['cl','crb','csq']
lfreq     = dfselect['freqv'].unique().tolist()
nlmethod  = len(lmethod)
nlshape   = len(lshape)
=======
       
dfselect  = dfcomp_short[mask]
lmethod   = ['spatte','spectetheta','dispte','displs','specls']
lfreq     = dfselect['freqv'].unique().tolist()
nlmethod  = len(lmethod)
>>>>>>> origin/main
nlfreq    = len(lfreq)

lnormsoll2  = []
lnormsolmax = []

lnormrecl2  = []
lnormrecmax = []

lv1min = []
lv1max = []

lv2min = []
lv2max = []

lv3min = []
lv3max = []

lv4min = []
lv4max = []

for k0 in range(0,nlfreq):

    locfreq   = lfreq[k0]

    locsoll2  = []
    locsolmax = []

    locrecl2  = []
    locrecmax = []

    for k1 in range(0,nlmethod):

        locmethod = lmethod[k1]

<<<<<<< HEAD
        for k2 in range(0,nlshape):

            locshape = lshape[k2]

            mask = (dfcomp_short['method']  == locmethod) & (dfcomp_short['freqv'] == locfreq) & (dfcomp_short['shape'] == locshape) 
            
            if(locshape=='cl'and locmethod=='dispte'): continue
            elif(locshape=='cl' and locmethod=='displs'): continue
            elif(locshape=='cl' and locmethod=='specls'): continue
            elif(locshape=='crb' and locmethod=='spatte'): continue
            elif(locshape=='crb' and locmethod=='spectetheta'): continue
            elif(locshape=='csq' and locmethod=='spatte'): continue
            elif(locshape=='csq' and locmethod=='spectetheta'): continue

            v1 = dfselect[mask]['normsoll2'].to_numpy()
            v2 = dfselect[mask]['normsolmax'].to_numpy()
            v3 = dfselect[mask]['normrecl2'].to_numpy()
            v4 = dfselect[mask]['normrecmax'].to_numpy()

            print(locshape,locmethod)
                

            locsoll2.append(v1)
            locsolmax.append(v2)
            locrecl2.append(v3)
            locrecmax.append(v4)
            
            lv1min.append(min(v1))
            lv1max.append(max(v1))
            
            lv2min.append(min(v2))
            lv2max.append(max(v2))
            
            lv3min.append(min(v3))
            lv3max.append(max(v3))
            
            lv4min.append(min(v4))
            lv4max.append(max(v4))
=======
        mask = (dfcomp_short['method']  == locmethod) & (dfcomp_short['freqv'] == locfreq) 

        v1 = dfselect[mask]['normsoll2'].to_numpy()
        v2 = dfselect[mask]['normsolmax'].to_numpy()
        v3 = dfselect[mask]['normrecl2'].to_numpy()
        v4 = dfselect[mask]['normrecmax'].to_numpy()
        
        locsoll2.append(v1)
        locsolmax.append(v2)
        locrecl2.append(v3)
        locrecmax.append(v4)
        
        lv1min.append(min(v1))
        lv1max.append(max(v1))
        
        lv2min.append(min(v2))
        lv2max.append(max(v2))
        
        lv3min.append(min(v3))
        lv3max.append(max(v3))
        
        lv4min.append(min(v4))
        lv4max.append(max(v4))
>>>>>>> origin/main
    
    lnormsoll2.append(locsoll2)
    lnormsolmax.append(locsolmax)
    lnormrecl2.append(locrecl2)
    lnormrecmax.append(locrecmax)

vsolminl2  = min(np.array(lv1min))
vsolmaxl2  = max(np.array(lv1max))

vsolminmax = min(np.array(lv2min))
vsolmaxmax = max(np.array(lv2max))

vrecminl2  = min(np.array(lv3min))
vrecmaxl2  = max(np.array(lv3max))

vrecminmax = min(np.array(lv4min))
vrecmaxmax = max(np.array(lv4max))
#==============================================================================

#==============================================================================
# Comum Infos
#==============================================================================
ordersv      = [2*i for i in range(0,9)]
vticks       = ['s', '+', '+', '+',  '+',  '^',   '^',   'D',      'D','s']
vline        = ['-', '-', '-', '--', '-.', '--',  '-.',  '--',     '-.','-']
vcolors      = ['b', 'g', 'r', 'c',  'm',  'y',   'b',   'purple', 'teal','lime']
linep        = [vticks,vline,vcolors]
#==============================================================================

#==============================================================================
# Plot1
#==============================================================================
<<<<<<< HEAD
def plot1(lnorm,xpos,ypos,clnames,normtype,vmin,vmax,lfreq,ordersv,figname,figtype):
=======
def plot1(lnorm,xpos,ypos,clnames,normtype,vmin,vmax,lfreq,ordersv,figname):
>>>>>>> origin/main

    npos  = len(xpos)
    nfreq = len(lfreq)
    nordersv = len(ordersv)
        
<<<<<<< HEAD
    plt.figure(figsize = (20,18))

    if(figtype=='rec'):

        if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
        if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))

    if(figtype=='sol'):

        if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Displacement at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
        if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Displacement at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
        
    grid = plt.GridSpec(3,3,wspace=0.4,hspace=0.2)
=======
    plt.figure(figsize = (16,12))
    if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
    if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
    grid = plt.GridSpec(3,2,wspace=0.2,hspace=0.4)
>>>>>>> origin/main
        
    min_value = -0.2*vmin
    max_value = 1.2*vmax
    
    for k0 in range(0,npos):
        
        plt.subplot(grid[xpos[k0],ypos[k0]])

        for k1 in range(0,nfreq):
                
            plt.plot(ordersv[1::],lnorm[k1][k0],label=str(lfreq[k1])+'Hz',color=linep[2][k1],linestyle=linep[1][k1],marker=linep[0][k1])
            
        plt.grid()
        plt.title(clnames[k0])
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
                
    plt.savefig('%s/%s_%s.jpeg'%(locsave,figname,normtype),dpi=200,bbox_inches='tight')
        
    plt.close()
    
    return
#==============================================================================

#==============================================================================
# Plot2
#==============================================================================
<<<<<<< HEAD
def plot2(lnorm,xpos,ypos,clnames,normtype,vmin,vmax,lfreq,ordersv,figname,figtype):
=======
def plot2(lnorm,xpos,ypos,clnames,normtype,vmin,vmax,lfreq,ordersv,figname):
>>>>>>> origin/main

    npos     = len(xpos)
    nfreq    = len(lfreq)
    nordersv = len(ordersv)
    nmethods = len(clnames)
        
<<<<<<< HEAD
    plt.figure(figsize = (16,6))
    
    if(figtype=='rec'):

        if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
        if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))

    if(figtype=='sol'):

        if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Displacement at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
        if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Displacement at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))

    grid = plt.GridSpec(1,3,wspace=0.4,hspace=0.2)
=======
    plt.figure(figsize = (16,12))
    if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
    if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs'%(vparameters[0],vparameters[1],vparameters[2]))
    grid = plt.GridSpec(2,2,wspace=0.2,hspace=0.4)
>>>>>>> origin/main
        
    min_value = -0.2*vmin
    max_value = 1.2*vmax
    
    for k0 in range(0,nfreq):

        plt.subplot(grid[xpos[k0],ypos[k0]])

        for k1 in range(0,nmethods):
                
            plt.plot(ordersv[1::],lnorm[k0][k1],label=clnames[k1],color=linep[2][k1],linestyle=linep[1][k1],marker=linep[0][k1])
            
        plt.grid()
        plt.title(str(lfreq[k0])+'Hz')
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
                
    plt.savefig('%s/%s_%s.jpeg'%(locsave,figname,normtype),dpi=200,bbox_inches='tight')
        
    plt.close()
    
    return
#==============================================================================

#==============================================================================
# Plot Infos
#==============================================================================
<<<<<<< HEAD
xpos        = [0,0,1,1,1,2,2,2]
ypos        = [0,1,0,1,2,0,1,2]
clnames     = ['spatte','spectetheta',
               'dispte-crb N=1','dispte-csq N=0',
               'displs-crb N=1','displs-csq N=0',
               'specls-crb N=1','specls-csq N=0']

=======
xpos        = [0,0,1,1,2]
ypos        = [0,1,0,1,0]
clnames     = ['spatte','spectetheta','dispte-crb N=1','displs-crb N=1','specls-crb N=1']
>>>>>>> origin/main
locsave     = 'comp_fig/teste%d/compfreq/'%(ptype) 
#==============================================================================

#==============================================================================
# Plot1 Execute
#==============================================================================
normtype = 'n2'
<<<<<<< HEAD

figtype  = 'sol'
figname  = '1normsol_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,normtype)
P1       = plot1(lnormsoll2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)

figtype  = 'rec'
figname  = '1normrec_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,normtype)
P2       = plot1(lnormrecl2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)

normtype = 'nmax'

figtype  = 'sol'
figname  = '1normsol_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,normtype)
P3       = plot1(lnormsolmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)

figtype  = 'sol'
figname  = '1normrec_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,normtype)
P4       = plot1(lnormrecmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)
#==============================================================================

#==============================================================================
# Plot2 Execute
#==============================================================================
xpos     = [0,0,0]
ypos     = [0,1,2]

normtype = 'n2'

figtype  = 'sol'
figname  = '2normsol_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,figtype)
P1       = plot2(lnormsoll2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)

figtype  = 'rec'
figname  = '2normrec_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,figtype)
P2       = plot2(lnormrecl2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)

normtype = 'nmax'

figtype  = 'sol'
figname  = '2normsol_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,figtype)
P3       = plot2(lnormsolmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)

figtype  = 'rec'
figname  = '2normrec_p%ddx%ddt%d%s'%(ptype,dx_ref,dt_ref,figtype)
P4       = plot2(lnormrecmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname,figtype)
=======
figname  = '1normsol_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P1       = plot1(lnormsoll2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
figname  = '1normrec_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P2       = plot1(lnormrecl2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
normtype = 'nmax'
figname  = '1normsol_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P3       = plot1(lnormsolmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
figname  = '1normrec_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P4       = plot1(lnormrecmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
#==============================================================================

#==============================================================================
# Plot1 Execute - L2
#==============================================================================
normtype = 'n2'
figname  = '2normsol_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P1       = plot2(lnormsoll2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
figname  = '2normrec_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P2       = plot2(lnormrecl2,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
normtype = 'nmax'
figname  = '2normsol_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P3       = plot2(lnormsolmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
figname  = '2normrec_p%ddx%ddt%d'%(ptype,dx_ref,dt_ref)
P4       = plot2(lnormrecmax,xpos,ypos,clnames,normtype,vsolminl2,vsolmaxl2,lfreq,ordersv,figname)
>>>>>>> origin/main
#==============================================================================