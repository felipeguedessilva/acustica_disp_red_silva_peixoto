#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
from   numpy import linalg      as la
import sys
import os
import pickle
from scipy.interpolate       import interp1d,RectBivariateSpline
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
# Similarity Function
#==============================================================================
def simfun(a,b,param):
    
    if(np.abs(a-b)<param):

        sval = 1 - ((np.abs(a-b)/param))

    else:

        sval = 0.0

    return sval
#==============================================================================

#==============================================================================
# Norm 1
#==============================================================================
def normnew1(im1_gray,im2_gray,sqsize):

    locmax_full = max(np.amax(im1_gray),np.amax(im2_gray))
    locmin_full = min(np.amin(im1_gray),np.amin(im2_gray))
    vssim_full  = ssim(im1_gray, im2_gray,data_range=locmax_full-locmin_full)
    print(vssim_full)

    nlin = im1_gray.shape[0]
    ncol = im1_gray.shape[1]

    nlinsqz = int(np.floor(nlin/sqsize) + 1)
    ncolsqz = int(np.floor(ncol/sqsize) + 1)

    nsquare = nlinsqz*ncolsqz
    
    nlinnew = sqsize*nlinsqz
    ncolnew = sqsize*ncolsqz

    nlindif = nlinnew - nlin
    ncoldif = ncolnew - ncol

    im_gray_extend1  = np.pad(im1_gray,((0,nlindif),(0,ncoldif)),'edge')
    im_gray_extend2  = np.pad(im2_gray,((0,nlindif),(0,ncoldif)),'edge')
 
    # xnew1 = np.linspace(0,nlin-1,nlin)
    # ynew1 = np.linspace(0,ncol-1,ncol)

    # xnew2 = np.linspace(0,nlinnew-1,nlinnew)
    # ynew2 = np.linspace(0,ncolnew-1,ncolnew)

    # im_gray1_interp1 = RectBivariateSpline(xnew1,ynew1,im1_gray.astype(float))
    # im_gray_extend1  = im_gray1_interp1(xnew2,ynew2)

    # im_gray2_interp2 = RectBivariateSpline(xnew1,ynew1,im2_gray.astype(float))
    # im_gray_extend2  = im_gray2_interp2(xnew2,ynew2)

    imnorm = 0.0

    contx1 = 0
    contx2 = sqsize

    nsquarecont = nsquare

    for k1 in range(0,nlinsqz):

        conty1 = 0
        conty2 = sqsize
        
        for k2 in range(0,ncolsqz):

            lval1 = []
            lval2 = []

            for k3 in range(contx1,contx2):

                for k4 in range(conty1,conty2):

                    lval1.append(im_gray_extend1[k3,k4])
                    lval2.append(im_gray_extend2[k3,k4])
               
            lval1max = max(lval1)
            lval1min = min(lval1)
            
            lval2max = max(lval2)
            lval2min = min(lval2)

            param    = 1/2
            ha       = simfun(lval1max,lval1min,param)
            hb       = simfun(lval2max,lval2min,param)
            weight   = simfun(ha,hb,param)

            locfig1 = np.array(lval1).reshape(sqsize,sqsize)
            locfig2 = np.array(lval2).reshape(sqsize,sqsize)

            plt.imshow(locfig1,cmap='binary',interpolation='kaiser',aspect='auto',vmin=-10,vmax=10)
            plt.axis('off')
            plt.savefig('../testresults/im_locfig1.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
            plt.close()
            plt.imshow(locfig2,cmap='binary',interpolation='kaiser',aspect='auto',vmin=-10,vmax=10)
            plt.axis('off')
            plt.savefig('../testresults/im_locfig2.png',transparent = True, bbox_inches = 'tight', pad_inches = 0)
            plt.close()
            imloc1       = cv2.imread('../testresults/im_locfig1.png')
            imloc2       = cv2.imread('../testresults/im_locfig2.png')

            imloc1_gray  = cv2.cvtColor(imloc1, cv2.COLOR_BGR2GRAY)
            imloc2_gray  = cv2.cvtColor(imloc2, cv2.COLOR_BGR2GRAY)

            locmax    = max(np.amax(imloc1_gray),np.amax(imloc2_gray))
            locmin    = min(np.amin(imloc1_gray),np.amin(imloc2_gray))
            vssim     = ssim(imloc1_gray,imloc2_gray,data_range=locmax-locmin)

            if(np.isnan(vssim)==True): 
                
                vssim = 0.0
                nsquarecont = nsquarecont - 1
           
            # vhd       = hd(im1_gray, im2_gray)
            # vmse      = mse(im1_gray, im2_gray)
            # vnrmse    = nrmse(im1_gray, im2_gray)
            # vare      = are(im1_gray, im2_gray)
            # vnme      = nmi(im1_gray, im2_gray)
            # vpsnr     = psnr(im1_gray, im2_gray)
            # vvoi      = voi(im1_gray, im2_gray)

            conty1 = conty1 + sqsize
            conty2 = conty2 + sqsize            

            imnorm = imnorm + vssim

        contx1 = contx1 + sqsize
        contx2 = contx2 + sqsize

    imnorm = (1/nsquarecont)*imnorm

    return imnorm 
#==============================================================================

#==============================================================================
# IM Test
#==============================================================================
im1       = cv2.imread('../testresults/im_ref1.png')
im2       = cv2.imread('../testresults/im_num1.png')
im1_gray  = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
im2_gray  = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
#==============================================================================

#==============================================================================
# IM New Norm
#==============================================================================
sqsize   = 20
imnorms1 = normnew1(im1_gray,im2_gray,sqsize)
print('')
print('IM New Norm 1 - Square Size %d x %d'%(sqsize,sqsize))
print(imnorms1)
print('')
#==============================================================================