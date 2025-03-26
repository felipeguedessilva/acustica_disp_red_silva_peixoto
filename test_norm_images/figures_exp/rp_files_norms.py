#==============================================================================
def gera_norms(ptype,dx_ref,dt_ref,freq_ref):
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
                
            test_results = pickle.load(f) 
            ntr          = len(test_results)
            testname      = 'Heterogeneos Velocity Model'
            xpositionv    = np.array([500.0,1500.0, 500.0,1500.0])
            ypositionv    = np.array([500.0, 500.0,1500.0,1500.0])
            
        elif(ptype==3): 
                
            test_results = pickle.load(f) 
            ntr          = len(test_results)
            testname      = 'SEG/EAGE 2D Salt Velocity Model'
            xpositionv    = np.array([4000.0,4000.0,4000.0,6000.0,6000.0,6000.0,8000.0,8000.0,8000.0])   
            ypositionv    = np.array([2000.0,2500.0,3000.0,2000.0,2500.0,3000.0,2000.0,2500.0,3000.0]) 
    
        elif(ptype==4): 
                
            test_results = pickle.load(f)
            ntr          = len(test_results)
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
    
    lfdxdt_normrecrel1_spatte_cl                = np.zeros((ncl,1))
    lfdxdt_normrecrel1_spectetheta_cl           = np.zeros((ncl,1))
    lfdxdt_normrecrel1_dispte_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel1_displs_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel1_specls_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel1_dispte_csq               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel1_displs_csq               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel1_specls_csq               = np.zeros((ncl,ncl))
    
    lfdxdt_normrecrel2_spatte_cl                = np.zeros((ncl,1))
    lfdxdt_normrecrel2_spectetheta_cl           = np.zeros((ncl,1))
    lfdxdt_normrecrel2_dispte_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel2_displs_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel2_specls_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel2_dispte_csq               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel2_displs_csq               = np.zeros((ncl,ncl))
    lfdxdt_normrecrel2_specls_csq               = np.zeros((ncl,ncl))
    
    lfdxdt_normrecremax_spatte_cl               = np.zeros((ncl,1))
    lfdxdt_normrecremax_spectetheta_cl          = np.zeros((ncl,1))
    lfdxdt_normrecremax_dispte_crb              = np.zeros((ncl,ncl))
    lfdxdt_normrecremax_displs_crb              = np.zeros((ncl,ncl))
    lfdxdt_normrecremax_specls_crb              = np.zeros((ncl,ncl))
    lfdxdt_normrecremax_dispte_csq              = np.zeros((ncl,ncl))
    lfdxdt_normrecremax_displs_csq              = np.zeros((ncl,ncl))
    lfdxdt_normrecremax_specls_csq              = np.zeros((ncl,ncl))
    
    lfdxdt_normrecreim_spatte_cl                = np.zeros((ncl,1))
    lfdxdt_normrecreim_spectetheta_cl           = np.zeros((ncl,1))
    lfdxdt_normrecreim_dispte_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecreim_displs_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecreim_specls_crb               = np.zeros((ncl,ncl))
    lfdxdt_normrecreim_dispte_csq               = np.zeros((ncl,ncl))
    lfdxdt_normrecreim_displs_csq               = np.zeros((ncl,ncl))
    lfdxdt_normrecreim_specls_csq               = np.zeros((ncl,ncl))
    
    lfdxdt_normrecselectrel1_spatte_cl          = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectrel1_spectetheta_cl     = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectrel1_dispte_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel1_displs_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel1_specls_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel1_dispte_csq         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel1_displs_csq         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel1_specls_csq         = np.zeros((nprecselect,ncl,ncl))
    
    lfdxdt_normrecselectrel2_spatte_cl          = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectrel2_spectetheta_cl     = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectrel2_dispte_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel2_displs_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel2_specls_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel2_dispte_csq         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel2_displs_csq         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectrel2_specls_csq         = np.zeros((nprecselect,ncl,ncl))
    
    lfdxdt_normrecselectremax_spatte_cl         = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectremax_spectetheta_cl    = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectremax_dispte_crb        = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectremax_displs_crb        = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectremax_specls_crb        = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectremax_dispte_csq        = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectremax_displs_csq        = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectremax_specls_csq        = np.zeros((nprecselect,ncl,ncl))

    lfdxdt_normrecselectreim_spatte_cl          = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectreim_spectetheta_cl     = np.zeros((nprecselect,ncl,1))
    lfdxdt_normrecselectreim_dispte_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectreim_displs_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectreim_specls_crb         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectreim_dispte_csq         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectreim_displs_csq         = np.zeros((nprecselect,ncl,ncl))
    lfdxdt_normrecselectreim_specls_csq         = np.zeros((nprecselect,ncl,ncl))

    lfdxdt_normsolplotrel1_spatte_cl            = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotrel1_spectetheta_cl       = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotrel1_dispte_crb           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel1_displs_crb           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel1_specls_crb           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel1_dispte_csq           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel1_displs_csq           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel1_specls_csq           = np.zeros((npsolselect,ncl,ncl))
    
    lfdxdt_normsolplotrel2_spatte_cl            = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotrel2_spectetheta_cl       = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotrel2_dispte_crb           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel2_displs_crb           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel2_specls_crb           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel2_dispte_csq           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel2_displs_csq           = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotrel2_specls_csq           = np.zeros((npsolselect,ncl,ncl))
    
    lfdxdt_normsolplotremax_spatte_cl           = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotremax_spectetheta_cl      = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotremax_dispte_crb          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotremax_displs_crb          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotremax_specls_crb          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotremax_dispte_csq          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotremax_displs_csq          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotremax_specls_csq          = np.zeros((npsolselect,ncl,ncl))
    
    lfdxdt_normsolplotreim_spatte_cl           = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotreim_spectetheta_cl      = np.zeros((npsolselect,ncl,1))
    lfdxdt_normsolplotreim_dispte_crb          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotreim_displs_crb          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotreim_specls_crb          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotreim_dispte_csq          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotreim_displs_csq          = np.zeros((npsolselect,ncl,ncl))
    lfdxdt_normsolplotreim_specls_csq          = np.zeros((npsolselect,ncl,ncl))
    
    for k1 in range(0,ncl):
        
        lfdxdt_normrecrel1_spatte_cl[k1,0]  = lfdxdt_select_spatte_cl[k1][11][3]
        lfdxdt_normrecrel2_spatte_cl[k1,0]  = lfdxdt_select_spatte_cl[k1][11][4]
        lfdxdt_normrecremax_spatte_cl[k1,0] = lfdxdt_select_spatte_cl[k1][11][5]
        lfdxdt_normrecreim_spatte_cl[k1,0]  = lfdxdt_select_spatte_cl[k1][11][6]

        lfdxdt_normrecrel1_spectetheta_cl[k1,0]  = lfdxdt_select_spectetheta_cl[k1][11][3]
        lfdxdt_normrecrel2_spectetheta_cl[k1,0]  = lfdxdt_select_spectetheta_cl[k1][11][4]
        lfdxdt_normrecremax_spectetheta_cl[k1,0] = lfdxdt_select_spectetheta_cl[k1][11][5]
        lfdxdt_normrecreim_spectetheta_cl[k1,0]  = lfdxdt_select_spectetheta_cl[k1][11][6]
    
        for k2 in range(0,nprecselect):
            
            lfdxdt_normrecselectrel1_spatte_cl[k2,k1,0]  = lfdxdt_select_spatte_cl[k1][12][3][k2]
            lfdxdt_normrecselectrel2_spatte_cl[k2,k1,0]  = lfdxdt_select_spatte_cl[k1][12][4][k2]
            lfdxdt_normrecselectremax_spatte_cl[k2,k1,0] = lfdxdt_select_spatte_cl[k1][12][5][k2]
            lfdxdt_normrecselectreim_spatte_cl[k2,k1,0]  = lfdxdt_select_spatte_cl[k1][12][6][k2]
            
            lfdxdt_normrecselectrel1_spectetheta_cl[k2,k1,0]  = lfdxdt_select_spectetheta_cl[k1][12][3][k2]
            lfdxdt_normrecselectrel2_spectetheta_cl[k2,k1,0]  = lfdxdt_select_spectetheta_cl[k1][12][4][k2]
            lfdxdt_normrecselectremax_spectetheta_cl[k2,k1,0] = lfdxdt_select_spectetheta_cl[k1][12][5][k2]
            lfdxdt_normrecselectreim_spectetheta_cl[k2,k1,0]  = lfdxdt_select_spectetheta_cl[k1][12][6][k2]
        
        for k3 in range(0,npsolselect):
            
            lfdxdt_normsolplotrel1_spatte_cl[k3,k1,0]  = lfdxdt_select_spatte_cl[k1][13][3][k3]
            lfdxdt_normsolplotrel2_spatte_cl[k3,k1,0]  = lfdxdt_select_spatte_cl[k1][13][4][k3]
            lfdxdt_normsolplotremax_spatte_cl[k3,k1,0] = lfdxdt_select_spatte_cl[k1][13][5][k3]
            lfdxdt_normsolplotreim_spatte_cl[k3,k1,0]  = lfdxdt_select_spatte_cl[k1][13][6][k3]
    
            lfdxdt_normsolplotrel1_spectetheta_cl[k3,k1,0]  = lfdxdt_select_spectetheta_cl[k1][13][3][k3]
            lfdxdt_normsolplotrel2_spectetheta_cl[k3,k1,0]  = lfdxdt_select_spectetheta_cl[k1][13][4][k3]
            lfdxdt_normsolplotremax_spectetheta_cl[k3,k1,0] = lfdxdt_select_spectetheta_cl[k1][13][5][k3]
            lfdxdt_normsolplotreim_spectetheta_cl[k3,k1,0]  = lfdxdt_select_spectetheta_cl[k1][13][6][k3]

    contglob        = 0
    contloc         = 0
    npte_dispte_crb = []
    npte_displs_crb = []
    npte_specls_crb = []
    npte_dispte_csq = []
    npte_displs_csq = []
    npte_specls_csq = []
    
    for k1 in range(0,ncl):
    
        for k2 in range(0,contloc+1):
    
            npte_dispte_crb.append([lfdxdt_select_dispte_crb[contglob][7],lfdxdt_select_dispte_crb[contglob][8],lfdxdt_select_dispte_crb[contglob][9],lfdxdt_select_dispte_crb[contglob][10]])
            npte_displs_crb.append([lfdxdt_select_displs_crb[contglob][7],lfdxdt_select_displs_crb[contglob][8],lfdxdt_select_displs_crb[contglob][9],lfdxdt_select_displs_crb[contglob][10]])
            npte_specls_crb.append([lfdxdt_select_specls_crb[contglob][7],lfdxdt_select_specls_crb[contglob][8],lfdxdt_select_specls_crb[contglob][9],lfdxdt_select_specls_crb[contglob][10]])
            
            npte_dispte_csq.append([lfdxdt_select_dispte_csq[contglob][7],lfdxdt_select_dispte_csq[contglob][8],lfdxdt_select_dispte_csq[contglob][9],lfdxdt_select_dispte_csq[contglob][10]])
            npte_displs_csq.append([lfdxdt_select_displs_csq[contglob][7],lfdxdt_select_displs_csq[contglob][8],lfdxdt_select_displs_csq[contglob][9],lfdxdt_select_displs_csq[contglob][10]])
            npte_specls_csq.append([lfdxdt_select_specls_csq[contglob][7],lfdxdt_select_specls_csq[contglob][8],lfdxdt_select_specls_csq[contglob][9],lfdxdt_select_specls_csq[contglob][10]])
            
            lfdxdt_normrecrel1_dispte_crb[k2,k1]  = lfdxdt_select_dispte_crb[contglob][11][3]
            lfdxdt_normrecrel2_dispte_crb[k2,k1]  = lfdxdt_select_dispte_crb[contglob][11][4]
            lfdxdt_normrecremax_dispte_crb[k2,k1] = lfdxdt_select_dispte_crb[contglob][11][5]
            lfdxdt_normrecreim_dispte_crb[k2,k1]  = lfdxdt_select_dispte_crb[contglob][11][6]
    
            lfdxdt_normrecrel1_displs_crb[k2,k1]  = lfdxdt_select_displs_crb[contglob][11][3]
            lfdxdt_normrecrel2_displs_crb[k2,k1]  = lfdxdt_select_displs_crb[contglob][11][4]
            lfdxdt_normrecremax_displs_crb[k2,k1] = lfdxdt_select_displs_crb[contglob][11][5]
            lfdxdt_normrecreim_displs_crb[k2,k1]  = lfdxdt_select_displs_crb[contglob][11][6]
            
            lfdxdt_normrecrel1_specls_crb[k2,k1]  = lfdxdt_select_specls_crb[contglob][11][3]
            lfdxdt_normrecrel2_specls_crb[k2,k1]  = lfdxdt_select_specls_crb[contglob][11][4]
            lfdxdt_normrecremax_specls_crb[k2,k1] = lfdxdt_select_specls_crb[contglob][11][5]
            lfdxdt_normrecreim_specls_crb[k2,k1]  = lfdxdt_select_specls_crb[contglob][11][6]

            lfdxdt_normrecrel1_dispte_csq[k2,k1]  = lfdxdt_select_dispte_csq[contglob][11][3]
            lfdxdt_normrecrel2_dispte_csq[k2,k1]  = lfdxdt_select_dispte_csq[contglob][11][4]
            lfdxdt_normrecremax_dispte_csq[k2,k1] = lfdxdt_select_dispte_csq[contglob][11][5]
            lfdxdt_normrecreim_dispte_csq[k2,k1]  = lfdxdt_select_dispte_csq[contglob][11][6]
            
            lfdxdt_normrecrel1_displs_csq[k2,k1]  = lfdxdt_select_displs_csq[contglob][11][3]
            lfdxdt_normrecrel2_displs_csq[k2,k1]  = lfdxdt_select_displs_csq[contglob][11][4]
            lfdxdt_normrecremax_displs_csq[k2,k1] = lfdxdt_select_displs_csq[contglob][11][5]
            lfdxdt_normrecreim_displs_csq[k2,k1]  = lfdxdt_select_displs_csq[contglob][11][6]
            
            lfdxdt_normrecrel1_specls_csq[k2,k1]  = lfdxdt_select_specls_csq[contglob][11][3]
            lfdxdt_normrecrel2_specls_csq[k2,k1]  = lfdxdt_select_specls_csq[contglob][11][4]
            lfdxdt_normrecremax_specls_csq[k2,k1] = lfdxdt_select_specls_csq[contglob][11][5]
            lfdxdt_normrecreim_specls_csq[k2,k1]  = lfdxdt_select_specls_csq[contglob][11][6]
            
            for k3 in range(0,nprecselect):
                
                lfdxdt_normrecselectrel1_dispte_crb[k3,k2,k1]  = lfdxdt_select_dispte_crb[contglob][12][3][k3]
                lfdxdt_normrecselectrel2_dispte_crb[k3,k2,k1]  = lfdxdt_select_dispte_crb[contglob][12][4][k3]
                lfdxdt_normrecselectremax_dispte_crb[k3,k2,k1] = lfdxdt_select_dispte_crb[contglob][12][5][k3]
                lfdxdt_normrecselectreim_dispte_crb[k3,k2,k1]  = lfdxdt_select_dispte_crb[contglob][12][6][k3]

                lfdxdt_normrecselectrel1_displs_crb[k3,k2,k1]  = lfdxdt_select_displs_crb[contglob][12][3][k3]
                lfdxdt_normrecselectrel2_displs_crb[k3,k2,k1]  = lfdxdt_select_displs_crb[contglob][12][4][k3]
                lfdxdt_normrecselectremax_displs_crb[k3,k2,k1] = lfdxdt_select_displs_crb[contglob][12][5][k3]
                lfdxdt_normrecselectreim_displs_crb[k3,k2,k1]  = lfdxdt_select_displs_crb[contglob][12][6][k3]

                lfdxdt_normrecselectrel1_specls_crb[k3,k2,k1]  = lfdxdt_select_specls_crb[contglob][12][3][k3]
                lfdxdt_normrecselectrel2_specls_crb[k3,k2,k1]  = lfdxdt_select_specls_crb[contglob][12][4][k3]
                lfdxdt_normrecselectremax_specls_crb[k3,k2,k1] = lfdxdt_select_specls_crb[contglob][12][5][k3]
                lfdxdt_normrecselectreim_specls_crb[k3,k2,k1]  = lfdxdt_select_specls_crb[contglob][12][6][k3]

                lfdxdt_normrecselectrel1_dispte_csq[k3,k2,k1]  = lfdxdt_select_dispte_csq[contglob][12][3][k3]
                lfdxdt_normrecselectrel2_dispte_csq[k3,k2,k1]  = lfdxdt_select_dispte_csq[contglob][12][4][k3]
                lfdxdt_normrecselectremax_dispte_csq[k3,k2,k1] = lfdxdt_select_dispte_csq[contglob][12][5][k3]
                lfdxdt_normrecselectreim_dispte_csq[k3,k2,k1]  = lfdxdt_select_dispte_csq[contglob][12][6][k3]

                lfdxdt_normrecselectrel1_displs_csq[k3,k2,k1]  = lfdxdt_select_displs_csq[contglob][12][3][k3]
                lfdxdt_normrecselectrel2_displs_csq[k3,k2,k1]  = lfdxdt_select_displs_csq[contglob][12][4][k3]
                lfdxdt_normrecselectremax_displs_csq[k3,k2,k1] = lfdxdt_select_displs_csq[contglob][12][5][k3]
                lfdxdt_normrecselectreim_displs_csq[k3,k2,k1]  = lfdxdt_select_displs_csq[contglob][12][6][k3]

                lfdxdt_normrecselectrel1_specls_csq[k3,k2,k1]  = lfdxdt_select_specls_csq[contglob][12][3][k3]
                lfdxdt_normrecselectrel2_specls_csq[k3,k2,k1]  = lfdxdt_select_specls_csq[contglob][12][4][k3]
                lfdxdt_normrecselectremax_specls_csq[k3,k2,k1] = lfdxdt_select_specls_csq[contglob][12][5][k3]
                lfdxdt_normrecselectreim_specls_csq[k3,k2,k1]  = lfdxdt_select_specls_csq[contglob][12][6][k3]

            for k4 in range(0,npsolselect):
              
                lfdxdt_normsolplotrel1_dispte_crb[k4,k2,k1]  = lfdxdt_select_dispte_crb[contglob][13][3][k4]
                lfdxdt_normsolplotrel2_dispte_crb[k4,k2,k1]  = lfdxdt_select_dispte_crb[contglob][13][4][k4]
                lfdxdt_normsolplotremax_dispte_crb[k4,k2,k1] = lfdxdt_select_dispte_crb[contglob][13][5][k4]
                lfdxdt_normsolplotreim_dispte_crb[k4,k2,k1]  = lfdxdt_select_dispte_crb[contglob][13][6][k4]

                lfdxdt_normsolplotrel1_displs_crb[k4,k2,k1]  = lfdxdt_select_displs_crb[contglob][13][3][k4]
                lfdxdt_normsolplotrel2_displs_crb[k4,k2,k1]  = lfdxdt_select_displs_crb[contglob][13][4][k4]
                lfdxdt_normsolplotremax_displs_crb[k4,k2,k1] = lfdxdt_select_displs_crb[contglob][13][5][k4]
                lfdxdt_normsolplotreim_displs_crb[k4,k2,k1]  = lfdxdt_select_displs_crb[contglob][13][6][k4]

                lfdxdt_normsolplotrel1_specls_crb[k4,k2,k1]  = lfdxdt_select_specls_crb[contglob][13][3][k4]
                lfdxdt_normsolplotrel2_specls_crb[k4,k2,k1]  = lfdxdt_select_specls_crb[contglob][13][4][k4]
                lfdxdt_normsolplotremax_specls_crb[k4,k2,k1] = lfdxdt_select_specls_crb[contglob][13][5][k4]
                lfdxdt_normsolplotreim_specls_crb[k4,k2,k1]  = lfdxdt_select_specls_crb[contglob][13][6][k4]

                lfdxdt_normsolplotrel1_dispte_csq[k4,k2,k1]  = lfdxdt_select_dispte_csq[contglob][13][3][k4]
                lfdxdt_normsolplotrel2_dispte_csq[k4,k2,k1]  = lfdxdt_select_dispte_csq[contglob][13][4][k4]
                lfdxdt_normsolplotremax_dispte_csq[k4,k2,k1] = lfdxdt_select_dispte_csq[contglob][13][5][k4]
                lfdxdt_normsolplotreim_dispte_csq[k4,k2,k1]  = lfdxdt_select_dispte_csq[contglob][13][6][k4]

                lfdxdt_normsolplotrel1_displs_csq[k4,k2,k1]  = lfdxdt_select_displs_csq[contglob][13][3][k4]
                lfdxdt_normsolplotrel2_displs_csq[k4,k2,k1]  = lfdxdt_select_displs_csq[contglob][13][4][k4]
                lfdxdt_normsolplotremax_displs_csq[k4,k2,k1] = lfdxdt_select_displs_csq[contglob][13][5][k4]
                lfdxdt_normsolplotreim_displs_csq[k4,k2,k1]  = lfdxdt_select_displs_csq[contglob][13][6][k4]

                lfdxdt_normsolplotrel1_specls_csq[k4,k2,k1]  = lfdxdt_select_specls_csq[contglob][13][3][k4]
                lfdxdt_normsolplotrel2_specls_csq[k4,k2,k1]  = lfdxdt_select_specls_csq[contglob][13][4][k4]
                lfdxdt_normsolplotremax_specls_csq[k4,k2,k1] = lfdxdt_select_specls_csq[contglob][13][5][k4]
                lfdxdt_normsolplotreim_specls_csq[k4,k2,k1]  = lfdxdt_select_specls_csq[contglob][13][6][k4]

            contglob = contglob + 1
    
        contloc = contloc+1
        
    nnpte_dispte_crb = len(npte_dispte_crb)
    vnptt_crb        = []
    vnpte_crb        = []
    
    for k1 in range(0,nnpte_dispte_crb):
     
        vnptt_crb.append(npte_dispte_crb[k1][2])
        vnpte_crb.append(npte_dispte_crb[k1][3])

    vnptt_crb = sorted(set(vnptt_crb))
    vnpte_crb = sorted(set(vnpte_crb))
    
    nnpte_dispte_csq = len(npte_dispte_csq)
    vnptt_csq        = []
    vnpte_csq        = []
    
    for k1 in range(0,nnpte_dispte_csq):
     
        vnptt_csq.append(npte_dispte_csq[k1][2])
        vnpte_csq.append(npte_dispte_csq[k1][3])

    vnptt_csq = sorted(set(vnptt_csq))
    vnpte_csq = sorted(set(vnpte_csq))

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
    def plot1(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq):
    
        nvcl       = len(vcl)
        nordersv   = len(ordersv)
        nvnpte_crb = len(vnpte_crb)
        vnpte_csq  = [vnpte_csq[3*i] for i in range(0,int(len(vnpte_csq)/3))]  
        nvnpte_csq = len(vnpte_csq)
        
        plt.figure(figsize = (20,16))
        if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
        if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
        if(normtype=='nim'): plt.suptitle('IM Norm of Full Receivers at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
        grid = plt.GridSpec(3,3,wspace=0.4,hspace=0.2)    
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
            
        for k1 in range(1,7):
        
            plt.subplot(grid[xpos[k1],ypos[k1]])
            
            vcrb[k1-1][vcrb[k1-1]==0] = np.nan 
            
            fig1 = plt.imshow(np.transpose(vcrb[k1-1]),cmap='jet',interpolation='kaiser')  
            plt.grid()
            plt.title('%s'%(crbnames[k1-1]))
            plt.xlabel('[Number of Extra Points]')
            plt.ylabel('[Order]')
            ax = plt.gca()
            
            if(xpos[k1]==1):

                ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte_crb)]))
                ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte_crb))
                ax.xaxis.set_major_locator(plt.MaxNLocator(nvnpte_crb))

            if(xpos[k1]==2):

                ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte_csq)]))
                ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte_csq))
                ax.xaxis.set_major_locator(plt.MaxNLocator(nvnpte_csq))

            ax.yaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
            ax.yaxis.set_major_formatter(ticker.FixedFormatter(ordersv[1::]))
            cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
            cbar.ax.locator_params(nbins=5)
            cbar.ax.set_ylabel('[Error]')
            tick_locator = ticker.MaxNLocator(nbins=5)
            cbar.locator = tick_locator
            cbar.update_ticks()
    
        plt.savefig('%s/%s_%s.jpeg'%(locsave,figname,normtype),dpi=200,bbox_inches='tight')
        
        plt.close()
    
        return
    #==============================================================================
    
    #==============================================================================
    # Plot Routines 2
    #==============================================================================
    def plot2(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq):
    
        nvcl       = len(vcl)
        nordersv   = len(ordersv)
        nvnpte_crb = len(vnpte_crb)
        vnpte_csq  = [vnpte_csq[3*i] for i in range(0,int(len(vnpte_csq)/3))]  
        nvnpte_csq = len(vnpte_csq)        
        ntimes     = vcl[0].shape[0]
        
        for k3 in range(0,ntimes):
        
            plt.figure(figsize = (20,16))
            if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Displacements at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3][k3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
            if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Displacements at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3][k3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
            if(normtype=='nim'): plt.suptitle('IM Norm of Full Displacements at time  %.3f s \n dx = %.4fm - dt = %.4fs - freq = %.3f Hz \n %s'%(vparameters[3][k3],vparameters[1],vparameters[2],vparameters[6],vparameters[5]))
            grid = plt.GridSpec(3,3,wspace=0.4,hspace=0.2)    
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
                
            for k1 in range(1,7):
            
                plt.subplot(grid[xpos[k1],ypos[k1]])
                
                vcrb[k1-1][vcrb[k1-1]==0] = np.nan 
                
                fig1 = plt.imshow(np.transpose(vcrb[k1-1][k3,:]),cmap='jet',interpolation='kaiser')  
                plt.grid()
                plt.title('%s'%(crbnames[k1-1]))
                plt.xlabel('[Number of Extra Points]')
                plt.ylabel('[Order]')
                ax = plt.gca()
                
                if(xpos[k1]==1):

                    ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte_crb)]))
                    ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte_crb))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(nvnpte_crb))

                if(xpos[k1]==2):

                    ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte_csq)]))
                    ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte_csq))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(nvnpte_csq))

                ax.yaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(ordersv[1::]))
                cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                cbar.ax.locator_params(nbins=5)
                cbar.ax.set_ylabel('[Error]')
                tick_locator = ticker.MaxNLocator(nbins=5)
                cbar.locator = tick_locator
                cbar.update_ticks()
        
            plt.savefig('%s/%s_ntime%d_%s.jpeg'%(locsave,figname,k3,normtype),dpi=200,bbox_inches='tight')
            
            plt.close()
    
        return
    #==============================================================================
    
    #==============================================================================
    # Plot Routines 3
    #==============================================================================
    def plot3(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,xpositionv,ypositionv,normtype,vnpte_crb,vnpte_csq):
    
        nvcl       = len(vcl)
        nordersv   = len(ordersv)
        nvnpte_crb = len(vnpte_crb)
        vnpte_csq  = [vnpte_csq[3*i] for i in range(0,int(len(vnpte_csq)/3))]  
        nvnpte_csq = len(vnpte_csq)        
        ntimes     = vcl[0].shape[0]
        
        for k3 in range(0,ntimes):
        
            plt.figure(figsize = (20,16))
            if(normtype=='n2'): plt.suptitle('Relative Quadratic Error of Full Receivers at time  %.3f s \n dx = %.4f m - dt = %.4f s - freq = %.3f Hz \n %s \n xpos = %.2f m - ypos = %.2f m'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5],xpositionv[k3],ypositionv[k3]))
            if(normtype=='nmax'): plt.suptitle('Relative Maximum Error of Full Receivers at time  %.3f s \n dx = %.4f m - dt = %.4f s - freq = %.3f Hz \n %s \n xpos = %.2f m - ypos = %.2f m'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5],xpositionv[k3],ypositionv[k3]))
            if(normtype=='nim'): plt.suptitle('IM Norm of Full Receivers at time  %.3f s \n dx = %.4f m - dt = %.4f s - freq = %.3f Hz \n %s \n xpos = %.2f m - ypos = %.2f m'%(vparameters[3],vparameters[1],vparameters[2],vparameters[6],vparameters[5],xpositionv[k3],ypositionv[k3]))
            grid = plt.GridSpec(3,3,wspace=0.4,hspace=0.2)    
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
                
            for k1 in range(1,7):
            
                plt.subplot(grid[xpos[k1],ypos[k1]])
                
                vcrb[k1-1][vcrb[k1-1]==0] = np.nan 
                
                fig1 = plt.imshow(np.transpose(vcrb[k1-1][k3,:]),cmap='jet',interpolation='kaiser')  
                plt.grid()
                plt.title('%s'%(crbnames[k1-1]))
                plt.xlabel('[Number of Extra Points]')
                plt.ylabel('[Order]')
                ax = plt.gca()
                
                if(xpos[k1]==1):

                    ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte_crb)]))
                    ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte_crb))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(nvnpte_crb))
               
                if(xpos[k1]==2):

                    ax.xaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nvnpte_csq)]))
                    ax.xaxis.set_major_formatter(ticker.FixedFormatter(vnpte_csq))
                    ax.xaxis.set_major_locator(plt.MaxNLocator(nvnpte_csq))
                
                ax.yaxis.set_major_locator(ticker.FixedLocator([i for i in range(0,nordersv)]))
                ax.yaxis.set_major_formatter(ticker.FixedFormatter(ordersv[1::]))
                cbar = plt.colorbar(fig1,format='%.2e',orientation='vertical') 
                cbar.ax.locator_params(nbins=5)
                cbar.ax.set_ylabel('[Error]')
                tick_locator = ticker.MaxNLocator(nbins=5)
                cbar.locator = tick_locator
                cbar.update_ticks()
        
            plt.savefig('%s/%s_npos%d_%s.jpeg'%(locsave,figname,k3,normtype),dpi=200,bbox_inches='tight')
            
            plt.close()
    
        return
    #==============================================================================
    
    #==============================================================================
    # Loc Save
    #==============================================================================
    locsave = 'comp_fig/teste%d/dx%ddt%dfreq%d/'%(ptype,dx_ref,dt_ref,freq_ref) 
    #==============================================================================
    
    #==============================================================================
    normtype = 'n2'
    #==============================================================================
  
    #==============================================================================
    # Plot Infos
    #==============================================================================
    xpos        = [0,1,1,1,2,2,2]
    ypos        = [1,0,1,2,0,1,2]

    clnames     = ['spatte','spectetheta',
                   'dispte-crb N=1','displs-crb N=1','specls-crb N=1',
                   'dispte-csq N=0','displs-csq N=0','specls-csq N=0']
    
    crbnames    = ['dispte-crb','displs-crb','specls-crb',
                   'dispte-csq','displs-csq','specls-csq']
    #==============================================================================
    
    #==============================================================================
    # Plot1 Execute - L2
    #==============================================================================
    vcl         = [lfdxdt_normrecrel2_spatte_cl,lfdxdt_normrecrel2_spectetheta_cl,
                   lfdxdt_normrecrel2_dispte_crb[0,:],lfdxdt_normrecrel2_displs_crb[0,:],lfdxdt_normrecrel2_specls_crb[0,:],
                   lfdxdt_normrecrel2_dispte_csq[0,:],lfdxdt_normrecrel2_displs_csq[0,:],lfdxdt_normrecrel2_specls_csq[0,:]]
    
    vcrb        = [lfdxdt_normrecrel2_dispte_crb,lfdxdt_normrecrel2_displs_crb,lfdxdt_normrecrel2_specls_crb,
                   lfdxdt_normrecrel2_dispte_csq,lfdxdt_normrecrel2_displs_csq,lfdxdt_normrecrel2_specls_csq]
    
    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]
    
    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]
    
    figname     = 'normrec_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])
    
    #P1          = plot1(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================
    
    #==============================================================================
    # Plot2 Execute - L2
    #==============================================================================
    vcl         = [lfdxdt_normsolplotrel2_spatte_cl,lfdxdt_normsolplotrel2_spectetheta_cl,
                   lfdxdt_normsolplotrel2_dispte_crb[:,0,:],lfdxdt_normsolplotrel2_displs_crb[:,0,:],lfdxdt_normsolplotrel2_specls_crb[:,0,:],
                   lfdxdt_normsolplotrel2_dispte_csq[:,0,:],lfdxdt_normsolplotrel2_displs_csq[:,0,:],lfdxdt_normsolplotrel2_specls_csq[:,0,:]]
    
    vcrb        = [lfdxdt_normsolplotrel2_dispte_crb,lfdxdt_normsolplotrel2_displs_crb,lfdxdt_normsolplotrel2_specls_crb,
                   lfdxdt_normsolplotrel2_dispte_csq,lfdxdt_normsolplotrel2_displs_csq,lfdxdt_normsolplotrel2_specls_csq]
    
    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    np.array(timesolplot)/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]
    
    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]
    
    figname     = 'normsolplot_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])
    
    #P2          = plot2(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================
    
    #==============================================================================
    # Plot3 Execute - L2
    #==============================================================================
    vcl         = [lfdxdt_normrecselectrel2_spatte_cl,lfdxdt_normrecselectrel2_spectetheta_cl,
                   lfdxdt_normrecselectrel2_dispte_crb[:,0,:],lfdxdt_normrecselectrel2_displs_crb[:,0,:],lfdxdt_normrecselectrel2_specls_crb[:,0,:],
                   lfdxdt_normrecselectrel2_dispte_csq[:,0,:],lfdxdt_normrecselectrel2_displs_csq[:,0,:],lfdxdt_normrecselectrel2_specls_csq[:,0,:]]
    
    vcrb        = [lfdxdt_normrecselectrel2_dispte_crb,lfdxdt_normrecselectrel2_displs_crb,lfdxdt_normrecselectrel2_specls_crb,
                   lfdxdt_normrecselectrel2_dispte_csq,lfdxdt_normrecselectrel2_displs_csq,lfdxdt_normrecselectrel2_specls_csq]
    
    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]
    
    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]
    
    figname     = 'normrecselect_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])
    
    #P3          = plot3(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,xpositionv,ypositionv,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================
    
    #==============================================================================
    normtype = 'nmax'
    #==============================================================================
  
    #==============================================================================
    # Plot1 Execute - MAX
    #==============================================================================
    vcl         = [lfdxdt_normrecremax_spatte_cl,lfdxdt_normrecremax_spectetheta_cl,
                   lfdxdt_normrecremax_dispte_crb[0,:],lfdxdt_normrecremax_displs_crb[0,:],lfdxdt_normrecremax_specls_crb[0,:],
                   lfdxdt_normrecremax_dispte_csq[0,:],lfdxdt_normrecremax_displs_csq[0,:],lfdxdt_normrecremax_specls_csq[0,:]]

    vcrb        = [lfdxdt_normrecremax_dispte_crb,lfdxdt_normrecremax_displs_crb,lfdxdt_normrecremax_specls_crb,
                   lfdxdt_normrecremax_dispte_csq,lfdxdt_normrecremax_displs_csq,lfdxdt_normrecremax_specls_csq]

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'normrec_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P4          = plot1(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================

    #==============================================================================
    # Plot2 Execute - MAX
    #==============================================================================
    vcl         = [lfdxdt_normsolplotremax_spatte_cl,lfdxdt_normsolplotremax_spectetheta_cl,
                   lfdxdt_normsolplotremax_dispte_crb[:,0,:],lfdxdt_normsolplotremax_displs_crb[:,0,:],lfdxdt_normsolplotremax_specls_crb[:,0,:],
                   lfdxdt_normsolplotremax_dispte_csq[:,0,:],lfdxdt_normsolplotremax_displs_csq[:,0,:],lfdxdt_normsolplotremax_specls_csq[:,0,:]]

    vcrb        = [lfdxdt_normsolplotremax_dispte_crb,lfdxdt_normsolplotremax_displs_crb,lfdxdt_normsolplotremax_specls_crb,
                   lfdxdt_normsolplotremax_dispte_csq,lfdxdt_normsolplotremax_displs_csq,lfdxdt_normsolplotremax_specls_csq]

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    np.array(timesolplot)/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'normsolplot_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P5          = plot2(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================

    #==============================================================================
    # Plot3 Execute - MAX
    #==============================================================================
    vcl         = [lfdxdt_normrecselectremax_spatte_cl,lfdxdt_normrecselectremax_spectetheta_cl,
                   lfdxdt_normrecselectremax_dispte_crb[:,0,:],lfdxdt_normrecselectremax_displs_crb[:,0,:],lfdxdt_normrecselectremax_specls_crb[:,0,:],
                   lfdxdt_normrecselectremax_dispte_csq[:,0,:],lfdxdt_normrecselectremax_displs_csq[:,0,:],lfdxdt_normrecselectremax_specls_csq[:,0,:]]

    vcrb        = [lfdxdt_normrecselectremax_dispte_crb,lfdxdt_normrecselectremax_displs_crb,lfdxdt_normrecselectremax_specls_crb,
                   lfdxdt_normrecselectremax_dispte_csq,lfdxdt_normrecselectremax_displs_csq,lfdxdt_normrecselectremax_specls_csq]

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'normrecselect_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    #P6          = plot3(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,xpositionv,ypositionv,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================

    #==============================================================================
    normtype = 'nim'
    #==============================================================================
  
    #==============================================================================
    # Plot1 Execute - SSIM
    #==============================================================================
    vcl         = [lfdxdt_normrecreim_spatte_cl,lfdxdt_normrecreim_spectetheta_cl,
                   lfdxdt_normrecreim_dispte_crb[0,:],lfdxdt_normrecreim_displs_crb[0,:],lfdxdt_normrecreim_specls_crb[0,:],
                   lfdxdt_normrecreim_dispte_csq[0,:],lfdxdt_normrecreim_displs_csq[0,:],lfdxdt_normrecreim_specls_csq[0,:]]

    vcrb        = [lfdxdt_normrecreim_dispte_crb,lfdxdt_normrecreim_displs_crb,lfdxdt_normrecreim_specls_crb,
                   lfdxdt_normrecreim_dispte_csq,lfdxdt_normrecreim_displs_csq,lfdxdt_normrecreim_specls_csq]

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'normrec_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    P7          = plot1(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================

    #==============================================================================
    # Plot2 Execute - SSIM
    #==============================================================================
    vcl         = [lfdxdt_normsolplotreim_spatte_cl,lfdxdt_normsolplotreim_spectetheta_cl,
                   lfdxdt_normsolplotreim_dispte_crb[:,0,:],lfdxdt_normsolplotreim_displs_crb[:,0,:],lfdxdt_normsolplotreim_specls_crb[:,0,:],
                   lfdxdt_normsolplotreim_dispte_csq[:,0,:],lfdxdt_normsolplotreim_displs_csq[:,0,:],lfdxdt_normsolplotreim_specls_csq[:,0,:]]

    vcrb        = [lfdxdt_normsolplotreim_dispte_crb,lfdxdt_normsolplotreim_displs_crb,lfdxdt_normsolplotreim_specls_crb,
                   lfdxdt_normsolplotreim_dispte_csq,lfdxdt_normsolplotreim_displs_csq,lfdxdt_normsolplotreim_specls_csq]

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    np.array(timesolplot)/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'normsolplot_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    P8          = plot2(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================

    #==============================================================================
    # Plot3 Execute - SSIM
    #==============================================================================
    vcl         = [lfdxdt_normrecselectreim_spatte_cl,lfdxdt_normrecselectreim_spectetheta_cl,
                   lfdxdt_normrecselectreim_dispte_crb[:,0,:],lfdxdt_normrecselectreim_displs_crb[:,0,:],lfdxdt_normrecselectreim_specls_crb[:,0,:],
                   lfdxdt_normrecselectreim_dispte_csq[:,0,:],lfdxdt_normrecselectreim_displs_csq[:,0,:],lfdxdt_normrecselectreim_specls_csq[:,0,:]]

    vcrb        = [lfdxdt_normrecselectreim_dispte_crb,lfdxdt_normrecselectreim_displs_crb,lfdxdt_normrecselectreim_specls_crb,
                   lfdxdt_normrecselectreim_dispte_csq,lfdxdt_normrecselectreim_displs_csq,lfdxdt_normrecselectreim_specls_csq]

    vparameters = [lfdxdt_select_spatte_cl[0][15][6],lfdxdt_select_spatte_cl[0][15][7],lfdxdt_select_spatte_cl[0][15][13]/1000,
                    lfdxdt_select_spatte_cl[0][15][9]/1000,lfdxdt_select_spatte_cl[0][15][10],testname,lfdxdt_select_spatte_cl[0][15][10]*1000]

    vnamefig    = [lfdxdt_select[0][1],lfdxdt_select[0][2],lfdxdt_select[0][3],lfdxdt_select[0][4]]

    figname     = 'normrecselect_p%ddx%ddt%dfreq%d'%(vnamefig[0],vnamefig[1],vnamefig[2],vnamefig[3])

    P9          = plot3(vcl,clnames,xpos,ypos,vcrb,crbnames,vparameters,vnamefig,figname,locsave,ordersv,linep,xpositionv,ypositionv,normtype,vnpte_crb,vnpte_csq)
    #==============================================================================    
    
    #==============================================================================
    return
    #==============================================================================