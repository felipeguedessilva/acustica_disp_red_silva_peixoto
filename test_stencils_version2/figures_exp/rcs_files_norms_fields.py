#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
from   numpy import linalg      as la
import sys
import pickle 
#==============================================================================

#==============================================================================
# My Modules
#==============================================================================
import testes_opt              as ttopt
# import im_fuzzy_norms          as ifn
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
# Range of Parameters
#==============================================================================
save_fields   = 1
fuzzymember   = 0
from im_fuzzy_norms import fuzzynorms1 as fuzzynormfunction

# vptype    = [1,2,3,4] 
# vdxref    = [1]
# vdtref    = [1]
# vfreqref  = [1,2,3] 

vptype    = [1] 
vdxref    = [1]
vdtref    = [1]
vfreqref  = [1] 

nvptype      = len(vptype) 
nvdxref      = len(vdxref)
nvdtref      = len(vdtref)
nvfreqref    = len(vfreqref) 
cont_me      = 0
lconfig_me   = []
cont_glob    = 0

for k0 in range(0,nvptype):
    
    testresults = []

    for k1 in range(0,nvdxref):
        
        for k2 in range(0,nvdtref):
            
            for k3 in range(0,nvfreqref):

                ptype      = vptype[k0] 
                dx_ref     = vdxref[k1]
                dt_ref     = vdtref[k2]
                freq_ref   = vfreqref[k3]
                factor_ref = 16
                
                print('')
                print('Number    = %d'%cont_me)
                print('Test Type = %d'%ptype)
                print('dx_ref    = %d'%dx_ref)
                print('dt_ref    = %d'%dt_ref)
                print('freq_ref  = %d'%freq_ref)
                print('')
                
                lconfig_me.append((cont_me,ptype,dx_ref,dt_ref,freq_ref))
                cont_me = cont_me + 1
#==============================================================================

#==============================================================================
# Teste Select
#==============================================================================
                if(ptype==1): 
                    
                    teste_ref = ttopt.teste1_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste1(dx_ref,dt_ref,freq_ref)
                    
                if(ptype==2): 
                    
                    teste_ref = ttopt.teste2_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste2(dx_ref,dt_ref,freq_ref)
                
                if(ptype==3): 
                
                    teste_ref = ttopt.teste3_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste3(dx_ref,dt_ref,freq_ref)
                
                if(ptype==4): 
                    
                    teste_ref = ttopt.teste4_ref1(freq_ref,factor_ref)
                    teste     = ttopt.teste4(dx_ref,dt_ref,freq_ref)
#==============================================================================

#==============================================================================
# Vetores de Configurações
#==============================================================================
                print('Starting Configs!')
                vmethod  = ['spatte','spectetheta','dispte','specls','displs'] 
                nvmethod = len(vmethod)
                vmvalue  = [1,2,3,4,5,6,7,8]
                nvmvalue = len(vmvalue)
                
                total_configs = 0
                list_config   = []
                    
                for m1 in range(0,nvmethod):
                        
                    method = vmethod[m1]
                        
                    if(method == 'spatte'):
                            
                        vshape = ['cl']
                            
                    elif(method == 'spectetheta'):
                            
                        vshape = ['cl']
                        
                    elif(method == 'dispte' or method == 'specls' or method == 'displs'):
                            
                        #vshape = ['crb']
                        vshape = ['crb','csq'] 
                        #vshape = ['rb',crb','csq','sq']
                            
                    nvshape = len(vshape)
                        
                    for m2 in range(0,nvshape):
                        
                        shape = vshape[m2]
                    
                        for m3 in range(0,nvmvalue):
                                
                            mvalue = vmvalue[m3]
                                
                            if(shape == 'cl'):
                                    
                                vnvalue = [1]
                                
                            elif(shape == 'rb'):
                                    
                                vnvalue = [mvalue]
                                
                            else:
                                
                                vnvalue  = np.arange(1,mvalue+1)
                                
                            nvnvalue = len(vnvalue)
                                
                            for m4 in range(0,nvnvalue):
                                    
                                nvalue = vnvalue[m4]
                                config  = (shape,method,mvalue,nvalue,total_configs)                            
                                total_configs = total_configs + 1
                                list_config.append(config)
                    
                nconfig     = len(list_config)
                print('Finishing Configs!')
#==============================================================================

#==============================================================================
# Open and Close Referencia
#==============================================================================    
                print('')
                print('Open Ref Files!')
                locopenref     = 'teste%d/reffreq%d'%(ptype,freq_ref)
                rec_ref        = np.load("../data_save/%s/rec_ref.npy"%(locopenref))   
                solplot_ref    = np.load("../data_save/%s/solplot_ref.npy"%(locopenref))            
                rec_select_ref = np.load("../data_save/%s/rec_select_ref.npy"%(locopenref))
                print('')
                print('Close Ref Files!')
                print('')
#==============================================================================

#==============================================================================
# Open and Close Tests
#==============================================================================    
                locopen = 'teste%d/dx%ddt%dfreq%d'%(ptype,dx_ref,dt_ref,freq_ref)
                print('Open Test Files!')
#==============================================================================   

#==============================================================================
# Parameters Relation
#==============================================================================
                config        = list_config[0]
                mshape        = config[0]        
                method        = config[1]        
                sou           = int(2*config[2]) 
                mvalue        = int(config[2])  
                nvalue        = int(config[3])  
                    
                rec        = np.load("../data_save/%s/rec_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))    
                solplot    = np.load("../data_save/%s/solplot_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                rec_select = np.load("../data_save/%s/rec_select_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                    
                a              = rec_ref.shape[0]-1
                b              = rec.shape[0]-1
                c              = rec_ref.shape[1]-1
                d              = rec.shape[1]-1
                irecref0       = int(a/b)
                irecref1       = int(c/d)
                                        
                a              = solplot_ref.shape[0]-1
                b              = solplot.shape[0]-1
                c              = solplot_ref.shape[1]-1
                d              = solplot.shape[1]-1
                e              = solplot_ref.shape[2]-1
                f              = solplot.shape[2]-1
                isolplotref0   = int(a/b)
                isolplotref1   = int(c/d)
                isolplotref2   = int(e/f)
                                        
                a              = rec_select_ref.shape[0]-1
                b              = rec_select.shape[0]-1
                c              = rec_select_ref.shape[1]-1
                d              = rec_select.shape[1]-1
                irecselectref0 = int(a/b)
                irecselectref1 = int(c/d)
                  
                recrefcut    = rec_ref[0::irecref0,0::irecref1]
                solplotcut   = solplot_ref[0::isolplotref0,0::isolplotref1,0::isolplotref2]
                recselectcut = rec_select_ref[0::irecselectref0,0::irecselectref1]
                        
                if(ptype==1 and dt_ref==6): 
                            
                    recrefcut    = recrefcut[1:,:]
                    recselectcut = recselectcut[1:,:]

                locsave = 'signals/signal_files/teste%d/'%(ptype)
                np.save("%srec_select_ref.npy"%(locsave),recselectcut)
#==============================================================================                    

#==============================================================================
# Obtenção de Parâmetros para Open
#==============================================================================
                for k in range(0,nconfig):
                    
                    print('')
                    print('Test with Stencil: %d'%(k))
                    
                    nptx    = teste.nptx    
                    npty    = teste.npty    
                    x0      = teste.x0      
                    y0      = teste.y0      
                    x1      = teste.x1   
                    y1      = teste.y1   
                    hxv     = teste.hx      
                    hyv     = teste.hy      
                    t0      = teste.t0      
                    tn      = teste.tn      
                    f0      = teste.f0      
                    CFL     = teste.CFL
                    jump    = teste.jump
                    dt0     = 0.5*dt_ref
                    nt      = int((tn-t0)*dt0)
                
                    parameters = []
                    parameters = [nptx,npty,x0,y0,x1,y1,hxv,hyv,t0,tn,f0,CFL,jump,dt0,nt] 
                    
                    config        = list_config[k]
                    shape         = config[0]        
                    teste.shape   = shape
                    method        = config[1]        
                    teste.method  = method
                    sou           = int(2*config[2]) 
                    teste.sou     = sou
                    mvalue        = int(config[2])  
                    teste.mvalue  = mvalue
                    nvalue        = int(config[3])  
                    teste.nvalue  = nvalue
                    print('shape: %s - method: %s - sou: %d - mvalue: %d - nvalue: %d'%(shape,method,sou,mvalue,nvalue))
                    
                    if(shape=='cl'):
                        
                        npt = 2*mvalue + 1
                        npe = 0 
                        
                    elif(shape=='crb'):
                        
                        if(mvalue%2==0):
                        
                            npt  = int((nvalue**2+4*mvalue+1))
                        
                        else:
                            
                            npt  = int((nvalue**2-1+4*mvalue+1))
                            
                        npe = npt - (4*mvalue+1) 
                    
                    elif(shape=='csq'):
                        
                        npt     = int(4*((2*mvalue-nvalue+1)*nvalue-nvalue)+4*mvalue+1)  
                        npe     = npt - (4*mvalue+1) 
                            
                    sou         = teste.sou    
                    mvalue      = teste.mvalue  
                    nvalue      = teste.nvalue  
                    mshape      = teste.shape   
                    method      = teste.method
                    
                    try:
                    
                        rec         = np.load("../data_save/%s/rec_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))    
                        solplot     = np.load("../data_save/%s/solplot_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                        rec_select  = np.load("../data_save/%s/rec_select_%s_%s_%d_%d.npy"%(locopen,mshape,method,mvalue,nvalue))
                        fmark       = 0
                        
                        locsave = 'signals/signal_files/teste%d/'%(ptype)
                        np.save("%srec_select_%s_%s_%d_%d.npy"%(locsave,mshape,method,mvalue,nvalue),rec_select)

                    except:
                         
                        rec        = recrefcut.copy()
                        solplot    = solplotcut.copy()
                        rec_select = recselectcut.copy()
                        fmark      = 1

                    fields_save = []
                    normrec     = []
                                        
                    try:

                        if(fmark==1):
                            
                            normrec1      = np.nan
                            normrec2      = np.nan
                            normrecmax    = np.nan                        
                            normrecrel1   = np.nan
                            normrecrel2   = np.nan
                            normrecrelmax = np.nan
                            normrecim     = np.nan
                            
                        else:

                            normrec1      = la.norm(recrefcut-rec,1)
                            normrec2      = la.norm(recrefcut-rec,2)
                            normrecmax    = la.norm(recrefcut-rec,np.inf)                        
                            normrecrel1   = la.norm(recrefcut-rec,1)/la.norm(recrefcut,1)            
                            normrecrel2   = la.norm(recrefcut-rec,2)/la.norm(recrefcut,2)
                            normrecrelmax = la.norm(recrefcut-rec,np.inf)/la.norm(recrefcut,np.inf)
                            
                            plt.imshow(recrefcut,cmap='binary',interpolation='kaiser',aspect='auto',vmin=-10,vmax=10)
                            plt.axis('off')
                            plt.savefig('../testresults/im_ref1.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
                            plt.close()
                            plt.imshow(rec,cmap='binary',interpolation='kaiser',aspect='auto',vmin=-10,vmax=10)
                            plt.axis('off')
                            plt.savefig('../testresults/im_num1.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
                            plt.close()
                            im1       = cv2.imread('../testresults/im_ref1.png')
                            im2       = cv2.imread('../testresults/im_num1.png')
                            im1_gray  = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
                            im2_gray  = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
                            locmax    = max(np.amax(im1_gray),np.amax(im2_gray))
                            locmin    = min(np.amin(im1_gray),np.amin(im2_gray))
                            vssim     = ssim(im1_gray, im2_gray,data_range=locmax-locmin)
                            vhd       = hd(im1_gray, im2_gray)
                            vmse      = mse(im1_gray, im2_gray)
                            vnrmse    = nrmse(im1_gray, im2_gray)
                            vare      = are(im1_gray, im2_gray)
                            vnme      = nmi(im1_gray, im2_gray)
                            vpsnr     = psnr(im1_gray, im2_gray)
                            vvoi      = voi(im1_gray, im2_gray)
                            vfuzzy    = fuzzynormfunction(im1_gray,im2_gray,fuzzymember)
                            normrecim = vfuzzy

                    except:

                        normrec1      = np.nan
                        normrec2      = np.nan
                        normrecmax    = np.nan                        
                        normrecrel1   = np.nan
                        normrecrel2   = np.nan
                        normrecrelmax = np.nan
                        normrecim     = np.nan

                    normrec.append(normrec1)
                    normrec.append(normrec2)
                    normrec.append(normrecmax)                        
                    normrec.append(normrecrel1)
                    normrec.append(normrecrel2)
                    normrec.append(normrecrelmax)
                    normrec.append(normrecim)

                    normrecselect = [] 
                    n1            = []
                    n2            = []
                    nmax          = []
                    n1rel         = []
                    n2rel         = []
                    nmaxrel       = []
                    nim           = []
                                        
                    for i in range(0,rec_select.shape[1]):
                        
                        try:
                            
                            if(fmark==1):
                            
                                normloc1      = np.nan
                                normloc2      = np.nan
                                normlocmax    = np.nan
                                normlocrel1   = np.nan
                                normlocrel2   = np.nan
                                normlocrelmax = np.nan
                                normlocim     = np.nan

                            else:
                            
                                normloc1      = la.norm(recselectcut[:,i]-rec_select[:,i],1)
                                normloc2      = la.norm(recselectcut[:,i]-rec_select[:,i],2)  
                                normlocmax    = la.norm(recselectcut[:,i]-rec_select[:,i],np.inf)  
                                normlocrel1   = la.norm(recselectcut[:,i]-rec_select[:,i],1)/la.norm(recselectcut[:,i],1)  
                                normlocrel2   = la.norm(recselectcut[:,i]-rec_select[:,i],2)/la.norm(recselectcut[:,i],2)  
                                normlocrelmax = la.norm(recselectcut[:,i]-rec_select[:,i],np.inf)/la.norm(recselectcut[:,i],np.inf)  
                                
                                plt.plot(recselectcut[:,i])
                                plt.axis('off')
                                plt.savefig('../testresults/im_ref2.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
                                plt.close()
                                plt.plot(rec_select[:,i])
                                plt.axis('off')
                                plt.savefig('../testresults/im_num2.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
                                plt.close()
                                im1       = cv2.imread('../testresults/im_ref2.png')
                                im2       = cv2.imread('../testresults/im_num2.png')
                                im1_gray  = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
                                im2_gray  = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
                                locmax    = max(np.amax(im1_gray),np.amax(im2_gray))
                                locmin    = min(np.amin(im1_gray),np.amin(im2_gray))
                                vssim     = ssim(im1_gray, im2_gray,data_range=locmax-locmin)      
                                vhd       = hd(im1_gray, im2_gray)
                                vmse      = mse(im1_gray, im2_gray)
                                vnrmse    = nrmse(im1_gray, im2_gray)
                                vare      = are(im1_gray, im2_gray)
                                vnme      = nmi(im1_gray, im2_gray)
                                vpsnr     = psnr(im1_gray, im2_gray)
                                vvoi      = voi(im1_gray, im2_gray)
                                vfuzzy    = fuzzynormfunction(im1_gray,im2_gray,fuzzymember)
                                normlocim = vfuzzy

                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                            nim.append(normlocim)
                            
                        except:
                            
                            normloc1      = np.nan
                            normloc2      = np.nan
                            normlocmax    = np.nan
                            normlocrel1   = np.nan
                            normlocrel2   = np.nan
                            normlocrelmax = np.nan
                            normlocim     = np.nan
                            
                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                            nim.append(normlocim)

                    normrecselect.append(n1)
                    normrecselect.append(n2)
                    normrecselect.append(nmax)
                    normrecselect.append(n1rel)
                    normrecselect.append(n2rel)
                    normrecselect.append(nmaxrel)
                    normrecselect.append(nim)

                    normsolplot   = []
                    n1            = []
                    n2            = []
                    nmax          = []
                    n1rel         = []
                    n2rel         = []
                    nmaxrel       = []
                    nim           = []
                    timesolplot   = [i*dt0*jump for i in range(0,solplot.shape[0])]
                                        
                    for i in range(0,solplot.shape[0]):
                        
                        try:
                            
                            if(fmark==1):
                            
                                normloc1      = np.nan
                                normloc2      = np.nan
                                normlocmax    = np.nan
                                normlocrel1   = np.nan
                                normlocrel2   = np.nan
                                normlocrelmax = np.nan
                                normlocim     = np.nan

                            else:
                            
                                normloc1      = la.norm(solplotcut[i,:]-solplot[i,:],1)
                                normloc2      = la.norm(solplotcut[i,:]-solplot[i,:],2)
                                normlocmax    = la.norm(solplotcut[i,:]-solplot[i,:],np.inf)
                                normlocrel1   = la.norm(solplotcut[i,:]-solplot[i,:],1)/la.norm(solplotcut[i,:],1)
                                normlocrel2   = la.norm(solplotcut[i,:]-solplot[i,:],2)/la.norm(solplotcut[i,:],2)
                                normlocrelmax = la.norm(solplotcut[i,:]-solplot[i,:],np.inf)/la.norm(solplotcut[i,:],np.inf)
                          
                                plt.imshow(solplotcut[i,:],cmap='binary',interpolation='kaiser',aspect='auto',vmin=-10,vmax=10)
                                plt.axis('off')
                                plt.savefig('../testresults/im_ref3.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
                                plt.close()
                                plt.imshow(solplot[i,:],cmap='binary',interpolation='kaiser',aspect='auto',vmin=-10,vmax=10)
                                plt.axis('off')
                                plt.savefig('../testresults/im_num3.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
                                plt.close()
                                im1       = cv2.imread('../testresults/im_ref3.png')
                                im2       = cv2.imread('../testresults/im_num3.png')
                                im1_gray  = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
                                im2_gray  = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
                                locmax    = max(np.amax(im1_gray),np.amax(im2_gray))
                                locmin    = min(np.amin(im1_gray),np.amin(im2_gray))
                                vssim     = ssim(im1_gray, im2_gray,data_range=locmax-locmin)
                                vhd       = hd(im1_gray, im2_gray)
                                vmse      = mse(im1_gray, im2_gray)
                                vnrmse    = nrmse(im1_gray, im2_gray)
                                vare      = are(im1_gray, im2_gray)
                                vnme      = nmi(im1_gray, im2_gray)
                                vpsnr     = psnr(im1_gray, im2_gray)
                                vvoi      = voi(im1_gray, im2_gray)
                                vfuzzy    = fuzzynormfunction(im1_gray,im2_gray,fuzzymember)
                                normlocim = vfuzzy
                                
                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                            nim.append(normlocim)
                            
                        except:
                            
                            normloc1      = np.nan
                            normloc2      = np.nan
                            normlocmax    = np.nan
                            normlocrel1   = np.nan
                            normlocrel2   = np.nan
                            normlocrelmax = np.nan
                            normlocim     = np.nan
                            
                            n1.append(normloc1)
                            n2.append(normloc2)
                            nmax.append(normlocmax)
                            n1rel.append(normlocrel1)
                            n2rel.append(normlocrel2)
                            nmaxrel.append(normlocrelmax)
                            nim.append(normlocim)
                    
                    normsolplot.append(n1)
                    normsolplot.append(n2)
                    normsolplot.append(nmax)
                    normsolplot.append(n1rel)
                    normsolplot.append(n2rel)
                    normsolplot.append(nmaxrel)
                    normsolplot.append(nim)

                    if(save_fields==1):

                        if(fmark==1):

                            rec[:]        = np.nan
                            solplot[:]    = np.nan
                            rec_select[:] = np.nan
                
                        fields_save.append(rec)
                        fields_save.append(recrefcut)
                        fields_save.append(rec_select)
                        fields_save.append(recselectcut)
                        fields_save.append(solplot)
                        fields_save.append(solplotcut)
                        testresults.append([cont_me,ptype,dx_ref,dt_ref,freq_ref,mshape,method,mvalue,nvalue,npt,npe,normrec,normrecselect,normsolplot,timesolplot,parameters,cont_glob,fields_save])
                    
                    else:
                        
                        testresults.append([cont_me,ptype,dx_ref,dt_ref,freq_ref,mshape,method,mvalue,nvalue,npt,npe,normrec,normrecselect,normsolplot,timesolplot,parameters,cont_glob])
                    
                    cont_glob = cont_glob + 1
#==============================================================================

#==============================================================================
                print('')
                print('Close Test Files!')
                print('')
#==============================================================================
 
#==============================================================================
# Save Results
#==============================================================================
    if(save_fields==0): locname = '../testresults/test%d_results_norms'%(ptype)
    if(save_fields==1): locname = '../testresults/test%d_results_norms_fields'%(ptype)
    
    with open(locname, 'wb') as f: 
        pickle.dump(testresults, f) 
#==============================================================================