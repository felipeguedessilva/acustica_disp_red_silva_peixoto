#==============================================================================
def gera_fields(ptype,dx_ref,dt_ref,freq_ref):
#==============================================================================

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
    # Image Comparisson
    #==============================================================================
    from skimage.metrics import structural_similarity as ssim
    from skimage.metrics import hausdorff_distance as hd
    from skimage.metrics import mean_squared_error as mse
    from skimage.metrics import normalized_root_mse as nrmse
    from skimage.metrics import adapted_rand_error as are
    from skimage.metrics import normalized_mutual_information as nmi
    from skimage.metrics import peak_signal_noise_ratio as psnr
    from skimage.metrics import variation_of_information as voi
    from skimage.io import imsave
    import cv2
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
    
    #==============================================================================
    # Read Txt Files
    #==============================================================================
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
    # Selecting Data
    #==============================================================================
    lf_select = []

    for k1 in range(0,ntr):
        
        if(test_results[k1][4]==freq_ref):
        
            lf_select.append(test_results[k1])   

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
    lfdxdt_select_dispte_csq      = []
    lfdxdt_select_displs_csq      = []
    lfdxdt_select_specls_csq      = []

    for k1 in range(0,nlfdxdt):
        
        if(lfdxdt_select[k1][5]=='cl'  and lfdxdt_select[k1][6]=='spatte'):      lfdxdt_select_spatte_cl.append(lfdxdt_select[k1])
        if(lfdxdt_select[k1][5]=='cl'  and lfdxdt_select[k1][6]=='spectetheta'): lfdxdt_select_spectetheta_cl.append(lfdxdt_select[k1])    
        if(lfdxdt_select[k1][5]=='crb' and lfdxdt_select[k1][6]=='dispte'):      lfdxdt_select_dispte_crb.append(lfdxdt_select[k1])
        if(lfdxdt_select[k1][5]=='crb' and lfdxdt_select[k1][6]=='displs'):      lfdxdt_select_displs_crb.append(lfdxdt_select[k1])
        if(lfdxdt_select[k1][5]=='crb' and lfdxdt_select[k1][6]=='specls'):      lfdxdt_select_specls_crb.append(lfdxdt_select[k1])
        if(lfdxdt_select[k1][5]=='csq' and lfdxdt_select[k1][6]=='dispte'):      lfdxdt_select_dispte_csq.append(lfdxdt_select[k1])
        if(lfdxdt_select[k1][5]=='csq' and lfdxdt_select[k1][6]=='displs'):      lfdxdt_select_displs_csq.append(lfdxdt_select[k1])
        if(lfdxdt_select[k1][5]=='csq' and lfdxdt_select[k1][6]=='specls'):      lfdxdt_select_specls_csq.append(lfdxdt_select[k1])

    ncl         = len(lfdxdt_select_spatte_cl)
    ncrb        = len(lfdxdt_select_specls_crb)
    ncsq        = len(lfdxdt_select_specls_csq)
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
            i41 = k1            
            i51 = k1
        
        elif(lfdxdt_select_dispte_crb[k1][7]==mv2 and lfdxdt_select_dispte_crb[k1][8]==nv2):
            
            i32 = k1
            i42 = k1
            i52 = k1
            
        if(lfdxdt_select_dispte_crb[k1][7]==mv3 and lfdxdt_select_dispte_crb[k1][8]==nv3):
            
            i33 = k1
            i43 = k1            
            i53 = k1

        if(lfdxdt_select_dispte_crb[k1][7]==mv4 and lfdxdt_select_dispte_crb[k1][8]==nv4):
            
            i34 = k1
            i44 = k1            
            i54 = k1

        if(lfdxdt_select_dispte_crb[k1][7]==mv5 and lfdxdt_select_dispte_crb[k1][8]==nv5):
            
            i35 = k1
            i45 = k1            
            i55 = k1
        
    nlfdxdt_select_dispte_csq = len(lfdxdt_select_dispte_csq)

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

    for k1 in range(0,nlfdxdt_select_dispte_csq):
        
        if(lfdxdt_select_dispte_csq[k1][7]==mv1 and lfdxdt_select_dispte_csq[k1][8]==nv1):
            
            i61 = k1            
            i71 = k1
            i81 = k1
        
        elif(lfdxdt_select_dispte_csq[k1][7]==mv2 and lfdxdt_select_dispte_csq[k1][8]==nv2):
            
            i62 = k1
            i72 = k1
            i82 = k1

        if(lfdxdt_select_dispte_csq[k1][7]==mv3 and lfdxdt_select_dispte_csq[k1][8]==nv3):
            
            i63 = k1
            i73 = k1            
            i83 = k1

        if(lfdxdt_select_dispte_csq[k1][7]==mv4 and lfdxdt_select_dispte_csq[k1][8]==nv4):
            
            i64 = k1
            i74 = k1            
            i84 = k1

        if(lfdxdt_select_dispte_csq[k1][7]==mv5 and lfdxdt_select_dispte_csq[k1][8]==nv5):
            
            i65 = k1
            i75 = k1            
            i85 = k1
    #==============================================================================

    #==============================================================================
    # Plot Routines 1
    #==============================================================================
    def plot1(vsols,vnames,xpos,ypos,extent,vparameters):
        
        plt.figure(figsize = (26,16))
        plt.suptitle('Difference of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
        grid  = plt.GridSpec(9,5,wspace=0.5,hspace=0.5) 
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
        
            plt.figure(figsize = (26,16))
            plt.suptitle('Difference of Full Displacement at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3][k3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
            grid  = plt.GridSpec(9,5,wspace=0.5,hspace=0.5)    
                
            for k1 in range(0,nfigs):
            
                plt.subplot(grid[xpos[k1],ypos[k1]])
                
                if(k1==0):
                    
                    fig1 = plt.imshow(np.transpose(vsols[k1][k3]),cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                    plt.grid()
                    plt.title('%s'%(vnames[k1]))
                    plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                    plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                    ax = plt.gca()
                    cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                    cbar.ax.locator_params(nbins=5)
                    tick_locator = ticker.MaxNLocator(nbins=5)
                    cbar.locator = tick_locator
                    cbar.update_ticks()
                
                else:
                    
                    fig1 = plt.imshow(np.transpose(vsols[0][k3]-vsols[k1][k3]),cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                    plt.grid()
                    plt.title('%s'%(vnames[k1]))
                    plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                    plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
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
        
            plt.figure(figsize = (26,16))
            plt.suptitle('Diference of Full Receivers at time  %.3f s \n dx = %.4f m - dt = %.4f s - freq = %.3f Hz \n %s \n xpos = %.2f m - ypos = %.2f m'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5],xpositionv[k3],ypositionv[k3]))
            grid  = plt.GridSpec(9,5,wspace=0.5,hspace=0.5)    
                
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
    # Plot Routines 4
    #==============================================================================
    def plot4(vsols,vnames,xpos,ypos,extent,vparameters):
                
        nfigs = len(xpos)
        vmin  = 0
        vmax  = 0
        
        for k1 in range(0,nfigs):
        
            vmin = min(vmin,np.amin(vsols[k1]))
            vmax = max(vmax,np.amax(vsols[k1]))
        
        factor = 50
        scale  = vmax/factor
        
        plt.figure(figsize = (26,16))
        plt.suptitle('Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
        grid  = plt.GridSpec(9,5,wspace=0.5,hspace=0.5)    
     
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

                fig1 = plt.imshow(vsols[k1],cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                plt.grid()
                plt.title('%s'%(vnames[k1]),fontsize=5)
                plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                ax = plt.gca()
                cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                cbar.ax.locator_params(nbins=5)
                tick_locator = ticker.MaxNLocator(nbins=5)
                cbar.locator = tick_locator
                cbar.update_ticks()
                
        plt.savefig('%s/%s.jpeg'%(locsave,figname),dpi=200,bbox_inches='tight')
        plt.close()

        return
    #==============================================================================

    #==============================================================================
    # Plot Routines 5
    #==============================================================================
    def plot5(vsols,vnames,xpos,ypos,extent,vparameters):
        
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
        
            list_ssim = []
            list_l2 = []
            
            for k1 in range(0,nfigs):
                
                arrayim1 = vsols[0][k3]
                plt.imshow(arrayim1,cmap='binary',interpolation='kaiser',aspect='auto',vmin=-scale,vmax=scale)
                plt.axis('off')
                plt.savefig('%s/%s.jpeg'%(locsave,'im_ref'),transparent = True, bbox_inches = 'tight', pad_inches = 0)
                plt.close()
                
                arrayim2 = vsols[k1][k3]
                plt.imshow(arrayim2,cmap='binary',interpolation='kaiser',aspect='auto',vmin=-scale,vmax=scale)
                plt.axis('off')
                plt.savefig('%s/%s.jpeg'%(locsave,'im_num'),transparent = True, bbox_inches = 'tight', pad_inches = 0)
                plt.close()
                
                l2norm = la.norm(arrayim1-arrayim2,2)/la.norm(arrayim1,2)
                
                # print(arrayim1.shape,arrayim2.shape)
                
                im1 = cv2.imread('%s/%s.jpeg'%(locsave,'im_ref'))
                im2 = cv2.imread('%s/%s.jpeg'%(locsave,'im_num'))
                
                # print(im1.shape,im2.shape)

                im1_gray = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
                im2_gray = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
                
                # print(im1_gray.shape,im2_gray.shape)

                if(k1>0):
                    (score, diff) = ssim(im1_gray, im2_gray, full=True)
                    #print("Image similarity with Image %f"%k1, score)
                    list_l2.append(l2norm)
                    list_ssim.append(score)        
            
            plt.clf()
            plt.figure(figsize = (26,16))
            plt.suptitle('Full Displacement at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3][k3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
            grid  = plt.GridSpec(9,5,wspace=0.5,hspace=0.5)    
                
            for k1 in range(0,nfigs):
            
                plt.subplot(grid[xpos[k1],ypos[k1]])
                
                if(k1==0):
                    
                    fig1 = plt.imshow(np.transpose(vsols[k1][k3]),cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                    plt.grid()
                    plt.title('%s'%(vnames[k1]))
                    plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                    plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                    ax = plt.gca()
                    cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                    cbar.ax.locator_params(nbins=5)
                    tick_locator = ticker.MaxNLocator(nbins=5)
                    cbar.locator = tick_locator
                    cbar.update_ticks()
                
                else:
        
                    fig1 = plt.imshow(np.transpose(vsols[k1][k3]),cmap='binary',interpolation='kaiser',extent=extent,aspect='auto',vmin=-scale,vmax=scale)   
                    plt.grid()
                    plt.title('%s - L2_NORM = %.4f - SSIM_NORM = %.4f'%(vnames[k1],list_l2[k1-1],list_ssim[k1-1]),fontsize=5)
                    plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                    plt.gca().yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f km'))
                    ax = plt.gca()
                    cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                    cbar.ax.locator_params(nbins=5)
                    tick_locator = ticker.MaxNLocator(nbins=5)
                    cbar.locator = tick_locator
                    cbar.update_ticks()
                    
            plt.savefig('%s/%s_ntime%d.jpeg'%(locsave,figname,k3),dpi=200,bbox_inches='tight')
            plt.close()

        return
    #==============================================================================

    #==============================================================================
    # Plot Routines 6
    #==============================================================================
    def plot6(vsols,vnames,xpos,ypos,extent,vparameters,xpositionv,ypositionv,vrectime):
        
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
        
            plt.figure(figsize = (26,16))
            plt.suptitle('Full Receivers at time  %.3f s \n dx = %.4f m - dt = %.4f s - freq = %.3f Hz \n %s \n xpos = %.2f m - ypos = %.2f m'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5],xpositionv[k3],ypositionv[k3]))
            grid  = plt.GridSpec(9,5,wspace=0.5,hspace=0.5)    
                
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
                    
                    plt.plot(vrectime,vsols[k1][:,k3])
                    plt.grid()
                    plt.title('%s'%(vnames[k1]))
                    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
                    plt.ylim((vmin,vmax))
                    plt.gca().xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f s'))
                    
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
    xpos        = [0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8,8,8,8,8]
    ypos        = [2,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4,0,1,2,3,4]
    fscale      = 10**(-3)
    extent1     = [fscale*lfdxdt_select_spatte_cl[0][15][2],fscale*lfdxdt_select_spatte_cl[0][15][4], fscale*lfdxdt_select_spatte_cl[0][15][9], fscale*lfdxdt_select_spatte_cl[0][15][8]]
    extent2     = [fscale*lfdxdt_select_spatte_cl[0][15][2],fscale*lfdxdt_select_spatte_cl[0][15][4], fscale*lfdxdt_select_spatte_cl[0][15][3], fscale*lfdxdt_select_spatte_cl[0][15][5]]
    vnames      = ['Reference',
                  'spatte-cl M=1','spatte-cl M=2','spatte-cl M=4','spatte-cl M=6','spatte-cl M=8',
                  'spectetheta-cl M=1','spectetheta-cl M=2','spectetheta-cl M=4','spectetheta-cl M=6','spectetheta-cl M=8',
                  'dispte-crb M=1 and N=1','dispte-crb M=4 and N=1','dispte-crb M=4 and N=4','dispte-crb M=8 and N=1','dispte-crb M=8 and N=8',
                  'displs-crb M=1 and N=1','displs-crb M=4 and N=1','displs-crb M=4 and N=4','displs-crb M=8 and N=1','displs-crb M=8 and N=8',
                  'specls-crb M=1 and N=1','specls-crb M=4 and N=1','specls-crb M=4 and N=4','specls-crb M=8 and N=1','specls-crb M=8 and N=8',
                  'dispte-csq M=1 and N=0','dispte-csq M=4 and N=0','dispte-csq M=4 and N=4','dispte-csq M=8 and N=0','dispte-csq M=8 and N=8',
                  'displs-csq M=1 and N=0','displs-csq M=4 and N=0','displs-csq M=4 and N=4','displs-csq M=8 and N=0','displs-csq M=8 and N=8',
                  'specls-csq M=1 and N=0','specls-csq M=4 and N=0','specls-csq M=4 and N=4','specls-csq M=8 and N=0','specls-csq M=8 and N=8']

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
              lfdxdt_select_specls_crb[i51][17][0],lfdxdt_select_specls_crb[i52][17][0],lfdxdt_select_specls_crb[i53][17][0],lfdxdt_select_specls_crb[i54][17][0],lfdxdt_select_specls_crb[i55][17][0],
              lfdxdt_select_dispte_csq[i61][17][0],lfdxdt_select_dispte_csq[i62][17][0],lfdxdt_select_dispte_csq[i63][17][0],lfdxdt_select_dispte_csq[i64][17][0],lfdxdt_select_dispte_csq[i65][17][0],
              lfdxdt_select_displs_csq[i71][17][0],lfdxdt_select_displs_csq[i72][17][0],lfdxdt_select_displs_csq[i73][17][0],lfdxdt_select_displs_csq[i74][17][0],lfdxdt_select_displs_csq[i75][17][0],
              lfdxdt_select_specls_csq[i81][17][0],lfdxdt_select_specls_csq[i82][17][0],lfdxdt_select_specls_csq[i83][17][0],lfdxdt_select_specls_csq[i84][17][0],lfdxdt_select_specls_csq[i85][17][0]] 

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'diffieldsrec_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P1          = plot1(vsols,vnames,xpos,ypos,extent1,vparameters)
    #==============================================================================

    #==============================================================================
    # Plot2 Execute
    #==============================================================================
    vsols = [lfdxdt_select_spatte_cl[0][17][5],
              lfdxdt_select_spatte_cl[i11][17][4],lfdxdt_select_spatte_cl[i12][17][4],lfdxdt_select_spatte_cl[i13][17][4],lfdxdt_select_spatte_cl[i14][17][4],lfdxdt_select_spatte_cl[i15][17][4],
              lfdxdt_select_spectetheta_cl[i21][17][4],lfdxdt_select_spectetheta_cl[i22][17][4],lfdxdt_select_spectetheta_cl[i23][17][4],lfdxdt_select_spectetheta_cl[i24][17][4],lfdxdt_select_spectetheta_cl[i25][17][4],
              lfdxdt_select_dispte_crb[i31][17][4],lfdxdt_select_dispte_crb[i32][17][4],lfdxdt_select_dispte_crb[i33][17][4],lfdxdt_select_dispte_crb[i34][17][4],lfdxdt_select_dispte_crb[i35][17][4],
              lfdxdt_select_displs_crb[i41][17][4],lfdxdt_select_displs_crb[i42][17][4],lfdxdt_select_displs_crb[i43][17][4],lfdxdt_select_displs_crb[i44][17][4],lfdxdt_select_displs_crb[i45][17][4],
              lfdxdt_select_specls_crb[i51][17][4],lfdxdt_select_specls_crb[i52][17][4],lfdxdt_select_specls_crb[i53][17][4],lfdxdt_select_specls_crb[i54][17][4],lfdxdt_select_specls_crb[i55][17][4],
              lfdxdt_select_dispte_csq[i61][17][4],lfdxdt_select_dispte_csq[i62][17][4],lfdxdt_select_dispte_csq[i63][17][4],lfdxdt_select_dispte_csq[i64][17][4],lfdxdt_select_dispte_csq[i65][17][4],
              lfdxdt_select_displs_csq[i71][17][4],lfdxdt_select_displs_csq[i72][17][4],lfdxdt_select_displs_csq[i73][17][4],lfdxdt_select_displs_csq[i74][17][4],lfdxdt_select_displs_csq[i75][17][4],
              lfdxdt_select_specls_csq[i81][17][4],lfdxdt_select_specls_csq[i82][17][4],lfdxdt_select_specls_csq[i83][17][4],lfdxdt_select_specls_csq[i84][17][4],lfdxdt_select_specls_csq[i85][17][4]] 

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    np.array(timesolplot)/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'diffieldssolplot_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P2          = plot2(vsols,vnames,xpos,ypos,extent2,vparameters)
    #==============================================================================

    #==============================================================================
    # Plot3 Execute
    #==============================================================================
    vsols = [lfdxdt_select_spatte_cl[0][17][3],
              lfdxdt_select_spatte_cl[i11][17][2],lfdxdt_select_spatte_cl[i12][17][2],lfdxdt_select_spatte_cl[i13][17][2],lfdxdt_select_spatte_cl[i14][17][2],lfdxdt_select_spatte_cl[i15][17][2],
              lfdxdt_select_spectetheta_cl[i21][17][2],lfdxdt_select_spectetheta_cl[i22][17][2],lfdxdt_select_spectetheta_cl[i23][17][2],lfdxdt_select_spectetheta_cl[i24][17][2],lfdxdt_select_spectetheta_cl[i25][17][2],
              lfdxdt_select_dispte_crb[i31][17][2],lfdxdt_select_dispte_crb[i32][17][2],lfdxdt_select_dispte_crb[i33][17][2],lfdxdt_select_dispte_crb[i34][17][2],lfdxdt_select_dispte_crb[i35][17][2],
              lfdxdt_select_displs_crb[i41][17][2],lfdxdt_select_displs_crb[i42][17][2],lfdxdt_select_displs_crb[i43][17][2],lfdxdt_select_displs_crb[i44][17][2],lfdxdt_select_displs_crb[i45][17][2],
              lfdxdt_select_specls_crb[i51][17][2],lfdxdt_select_specls_crb[i52][17][2],lfdxdt_select_specls_crb[i53][17][2],lfdxdt_select_specls_crb[i54][17][2],lfdxdt_select_specls_crb[i55][17][2],
              lfdxdt_select_dispte_csq[i61][17][2],lfdxdt_select_dispte_csq[i62][17][2],lfdxdt_select_dispte_csq[i63][17][2],lfdxdt_select_dispte_csq[i64][17][2],lfdxdt_select_dispte_csq[i65][17][2],
              lfdxdt_select_displs_csq[i71][17][2],lfdxdt_select_displs_csq[i72][17][2],lfdxdt_select_displs_csq[i73][17][2],lfdxdt_select_displs_csq[i74][17][2],lfdxdt_select_displs_csq[i75][17][2],
              lfdxdt_select_specls_csq[i81][17][2],lfdxdt_select_specls_csq[i82][17][2],lfdxdt_select_specls_csq[i83][17][2],lfdxdt_select_specls_csq[i84][17][2],lfdxdt_select_specls_csq[i85][17][2]] 

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'diffieldsrecselect_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P3          = plot3(vsols,vnames,xpos,ypos,extent1,vparameters,xpositionv,ypositionv,vrectime)
    #==============================================================================
    
    #==============================================================================
    # Plot4 Execute
    #==============================================================================
    vsols = [lfdxdt_select_spatte_cl[0][17][1],
              lfdxdt_select_spatte_cl[i11][17][0],lfdxdt_select_spatte_cl[i12][17][0],lfdxdt_select_spatte_cl[i13][17][0],lfdxdt_select_spatte_cl[i14][17][0],lfdxdt_select_spatte_cl[i15][17][0],
              lfdxdt_select_spectetheta_cl[i21][17][0],lfdxdt_select_spectetheta_cl[i22][17][0],lfdxdt_select_spectetheta_cl[i23][17][0],lfdxdt_select_spectetheta_cl[i24][17][0],lfdxdt_select_spectetheta_cl[i25][17][0],
              lfdxdt_select_dispte_crb[i31][17][0],lfdxdt_select_dispte_crb[i32][17][0],lfdxdt_select_dispte_crb[i33][17][0],lfdxdt_select_dispte_crb[i34][17][0],lfdxdt_select_dispte_crb[i35][17][0],
              lfdxdt_select_displs_crb[i41][17][0],lfdxdt_select_displs_crb[i42][17][0],lfdxdt_select_displs_crb[i43][17][0],lfdxdt_select_displs_crb[i44][17][0],lfdxdt_select_displs_crb[i45][17][0],
              lfdxdt_select_specls_crb[i51][17][0],lfdxdt_select_specls_crb[i52][17][0],lfdxdt_select_specls_crb[i53][17][0],lfdxdt_select_specls_crb[i54][17][0],lfdxdt_select_specls_crb[i55][17][0],
              lfdxdt_select_dispte_csq[i61][17][0],lfdxdt_select_dispte_csq[i62][17][0],lfdxdt_select_dispte_csq[i63][17][0],lfdxdt_select_dispte_csq[i64][17][0],lfdxdt_select_dispte_csq[i65][17][0],
              lfdxdt_select_displs_csq[i71][17][0],lfdxdt_select_displs_csq[i72][17][0],lfdxdt_select_displs_csq[i73][17][0],lfdxdt_select_displs_csq[i74][17][0],lfdxdt_select_displs_csq[i75][17][0],
              lfdxdt_select_specls_csq[i81][17][0],lfdxdt_select_specls_csq[i82][17][0],lfdxdt_select_specls_csq[i83][17][0],lfdxdt_select_specls_csq[i84][17][0],lfdxdt_select_specls_csq[i85][17][0]] 

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'fieldsrec_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    P4          = plot4(vsols,vnames,xpos,ypos,extent1,vparameters)
    #==============================================================================

    #==============================================================================
    # Plot5 Execute
    #==============================================================================
    vsols = [lfdxdt_select_spatte_cl[0][17][5],
              lfdxdt_select_spatte_cl[i11][17][4],lfdxdt_select_spatte_cl[i12][17][4],lfdxdt_select_spatte_cl[i13][17][4],lfdxdt_select_spatte_cl[i14][17][4],lfdxdt_select_spatte_cl[i15][17][4],
              lfdxdt_select_spectetheta_cl[i21][17][4],lfdxdt_select_spectetheta_cl[i22][17][4],lfdxdt_select_spectetheta_cl[i23][17][4],lfdxdt_select_spectetheta_cl[i24][17][4],lfdxdt_select_spectetheta_cl[i25][17][4],
              lfdxdt_select_dispte_crb[i31][17][4],lfdxdt_select_dispte_crb[i32][17][4],lfdxdt_select_dispte_crb[i33][17][4],lfdxdt_select_dispte_crb[i34][17][4],lfdxdt_select_dispte_crb[i35][17][4],
              lfdxdt_select_displs_crb[i41][17][4],lfdxdt_select_displs_crb[i42][17][4],lfdxdt_select_displs_crb[i43][17][4],lfdxdt_select_displs_crb[i44][17][4],lfdxdt_select_displs_crb[i45][17][4],
              lfdxdt_select_specls_crb[i51][17][4],lfdxdt_select_specls_crb[i52][17][4],lfdxdt_select_specls_crb[i53][17][4],lfdxdt_select_specls_crb[i54][17][4],lfdxdt_select_specls_crb[i55][17][4],
              lfdxdt_select_dispte_csq[i61][17][4],lfdxdt_select_dispte_csq[i62][17][4],lfdxdt_select_dispte_csq[i63][17][4],lfdxdt_select_dispte_csq[i64][17][4],lfdxdt_select_dispte_csq[i65][17][4],
              lfdxdt_select_displs_csq[i71][17][4],lfdxdt_select_displs_csq[i72][17][4],lfdxdt_select_displs_csq[i73][17][4],lfdxdt_select_displs_csq[i74][17][4],lfdxdt_select_displs_csq[i75][17][4],
              lfdxdt_select_specls_csq[i81][17][4],lfdxdt_select_specls_csq[i82][17][4],lfdxdt_select_specls_csq[i83][17][4],lfdxdt_select_specls_csq[i84][17][4],lfdxdt_select_specls_csq[i85][17][4]] 

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    np.array(timesolplot)/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'fieldssolplot_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P5          = plot5(vsols,vnames,xpos,ypos,extent2,vparameters)
    #==============================================================================

    #==============================================================================
    # Plot6 Execute
    #==============================================================================
    vsols = [lfdxdt_select_spatte_cl[0][17][3],
              lfdxdt_select_spatte_cl[i11][17][2],lfdxdt_select_spatte_cl[i12][17][2],lfdxdt_select_spatte_cl[i13][17][2],lfdxdt_select_spatte_cl[i14][17][2],lfdxdt_select_spatte_cl[i15][17][2],
              lfdxdt_select_spectetheta_cl[i21][17][2],lfdxdt_select_spectetheta_cl[i22][17][2],lfdxdt_select_spectetheta_cl[i23][17][2],lfdxdt_select_spectetheta_cl[i24][17][2],lfdxdt_select_spectetheta_cl[i25][17][2],
              lfdxdt_select_dispte_crb[i31][17][2],lfdxdt_select_dispte_crb[i32][17][2],lfdxdt_select_dispte_crb[i33][17][2],lfdxdt_select_dispte_crb[i34][17][2],lfdxdt_select_dispte_crb[i35][17][2],
              lfdxdt_select_displs_crb[i41][17][2],lfdxdt_select_displs_crb[i42][17][2],lfdxdt_select_displs_crb[i43][17][2],lfdxdt_select_displs_crb[i44][17][2],lfdxdt_select_displs_crb[i45][17][2],
              lfdxdt_select_specls_crb[i51][17][2],lfdxdt_select_specls_crb[i52][17][2],lfdxdt_select_specls_crb[i53][17][2],lfdxdt_select_specls_crb[i54][17][2],lfdxdt_select_specls_crb[i55][17][2],
              lfdxdt_select_dispte_csq[i61][17][2],lfdxdt_select_dispte_csq[i62][17][2],lfdxdt_select_dispte_csq[i63][17][2],lfdxdt_select_dispte_csq[i64][17][2],lfdxdt_select_dispte_csq[i65][17][2],
              lfdxdt_select_displs_csq[i71][17][2],lfdxdt_select_displs_csq[i72][17][2],lfdxdt_select_displs_csq[i73][17][2],lfdxdt_select_displs_csq[i74][17][2],lfdxdt_select_displs_csq[i75][17][2],
              lfdxdt_select_specls_csq[i81][17][2],lfdxdt_select_specls_csq[i82][17][2],lfdxdt_select_specls_csq[i83][17][2],lfdxdt_select_specls_csq[i84][17][2],lfdxdt_select_specls_csq[i85][17][2]] 

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'fieldsrecselect_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P6          = plot6(vsols,vnames,xpos,ypos,extent1,vparameters,xpositionv,ypositionv,vrectime)
    #==============================================================================

    #==============================================================================
    return
    #==============================================================================