#==============================================================================
# Python Modules and Imports
#==============================================================================
import numpy                    as np
from   numpy import linalg      as la
import sys
import os
import pickle
from scipy.interpolate       import interp1d,RectBivariateSpline
# os.system('clear')
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
# Fuzzy Operations
#==============================================================================
def fuzzyunion(aset,bset):

    cset = np.fmax(aset,bset)
    #cset = aset+bset-aset*bset
    #cset = np.fmin(np.ones((aset.shape[0],aset.shape[1])),aset+bset)

    return cset

def fuzzyintersect(aset,bset):

    cset = np.fmin(aset,bset)
    #cset = aset*bset
    #cset = np.fmax(np.zeros((aset.shape[0],aset.shape[1])),aset+bset-1) 

    return cset

def fuzzycomplement(aset):

    cset = 1 - aset

    return cset

def fuzzycardn(aset):

    cardn = np.sum(aset)

    return cardn

def fuzzyabarb(aset,bset):

    bcomplementset = fuzzycomplement(bset)

    cset = fuzzyintersect(aset,bcomplementset)

    return cset

def fuzzyatrib(aset,bset):

    if((aset==bset).all()==True):
        
        cset = np.zeros((aset.shape[0],aset.shape[1]))
    
    else:

        acomplementset = fuzzycomplement(aset)
        bcomplementset = fuzzycomplement(bset)

        abarb = fuzzyabarb(aset,bset)
        bbara = fuzzyabarb(bset,aset)
        cset  = fuzzyunion(abarb,bbara)

    # acomplementset = fuzzycomplement(aset)
    # bcomplementset = fuzzycomplement(bset)

    # abarb = fuzzyabarb(aset,bset)
    # bbara = fuzzyabarb(bset,aset)
    # cset  = fuzzyunion(abarb,bbara)

    return cset
#==============================================================================

#==============================================================================
# Fuzzy Membership 0 - Gaussian - 2D
#==============================================================================
def fuzzymember0(im_gray):

    im_gray  = im_gray.astype(float)
    vmin_im  = np.amin(im_gray)
    vmax_im  = np.amax(im_gray)
    nlin     = im_gray.shape[0]
    ncol     = im_gray.shape[1]
    im_fuzzy = np.zeros((nlin,ncol))

    for k1 in range(0,nlin):

        for k2 in range(0,ncol):

            val = im_gray[k1,k2]

            if(val<vmin_im): 
            
                val_fuzzy = 0.0
            
            elif(vmin_im<=val<=vmax_im): 
                
                val_fuzzy = np.exp(-val**2)

            else:

                val_fuzzy = 0.0

            im_fuzzy[k1,k2] = val_fuzzy

    return im_fuzzy
#==============================================================================

#==============================================================================
# Fuzzy Membership 1 - Triangular - 2D
#==============================================================================
def fuzzymember1(im_gray):

    im_gray  = im_gray.astype(float)
    vmin_im  = np.amin(im_gray)
    vmax_im  = np.amax(im_gray)
    nlin     = im_gray.shape[0]
    ncol     = im_gray.shape[1]
    im_fuzzy = np.zeros((nlin,ncol))
    vbar     = 0.5*(vmax_im+vmin_im)

    for k1 in range(0,nlin):

        for k2 in range(0,ncol):

            val = im_gray[k1,k2]

            if(val<vmin_im): 
            
                val_fuzzy = 0.0
            
            elif(vmin_im<=val<vbar): 
                
                val_fuzzy = (val-vmin_im)/(vbar-vmin_im)

            elif(val==vbar): 
                
                val_fuzzy = 1.0

            elif(vbar<val<=vmax_im): 
                
                val_fuzzy = (vmax_im-val)/(vmax_im-vbar)

            else:

                val_fuzzy = 0.0

            im_fuzzy[k1,k2] = val_fuzzy

    return im_fuzzy
#==============================================================================

#==============================================================================
# Fuzzy Membership 2 - Trapezoidal - 2D
#==============================================================================
def fuzzymember2(im_gray):

    im_gray  = im_gray.astype(float)
    vmin_im  = np.amin(im_gray)
    vmax_im  = np.amax(im_gray)
    nlin     = im_gray.shape[0]
    ncol     = im_gray.shape[1]
    im_fuzzy = np.zeros((nlin,ncol))
    vbar     = 0.5*(vmax_im+vmin_im)
    vbar1    = 0.5*(vbar+vmin_im)
    vbar2    = 0.5*(vmax_im+vbar)

    for k1 in range(0,nlin):

        for k2 in range(0,ncol):

            val = im_gray[k1,k2]

            if(val<vmin_im): 
            
                val_fuzzy = 0.0
            
            elif(vmin_im<=val<vbar1): 
                
                val_fuzzy = (val-vmin_im)/(vbar-vmin_im)

            elif(vbar1<=val<=vbar2): 
                
                val_fuzzy = 1.0

            elif(vbar2<val<=vmax_im): 
                
                val_fuzzy = (vmax_im-val)/(vmax_im-vbar)

            else:

                val_fuzzy = 0.0

            im_fuzzy[k1,k2] = val_fuzzy

    return im_fuzzy
#==============================================================================

#==============================================================================
# Fuzzy Norm S1
#==============================================================================
def fuzzynorms1(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    rval      = 2
    mval      = 255
    nval      = mval    
    imnorm    = np.sum(np.abs(im_fuzzy1-im_fuzzy2)**(rval))
    imnorm    = 1 - ((1/(mval*nval))*imnorm)**(1/rval)

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S2
#==============================================================================
def fuzzynorms2(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    imnorm1   = np.sum(np.abs(im_fuzzy1-im_fuzzy2))
    imnorm2   = np.sum(np.abs(im_fuzzy1+im_fuzzy2))
   
    if(imnorm2==0):
        
        imnorm = 0.0

    else:

         imnorm = 1 - (imnorm1/imnorm2)

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S3
#==============================================================================
def fuzzynorms3(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    mval      = 255
    nval      = mval
    imnorm    = np.sum((im_fuzzy1-im_fuzzy2)*np.log((1+im_fuzzy1)/(1+im_fuzzy2))+(im_fuzzy2-im_fuzzy1)*np.log((2-im_fuzzy1)/(2-im_fuzzy2)))
    imnorm    = 1 - ((1/(mval*nval*np.log(2)))*imnorm)
 
    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S4
#==============================================================================
def fuzzynorms4(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    ainterb     = fuzzyintersect(im_fuzzy1,im_fuzzy2)
    aunionb     = fuzzyunion(im_fuzzy1,im_fuzzy2)
    cardainterb = fuzzycardn(ainterb)
    cardaunionb = fuzzycardn(aunionb)

    if(cardaunionb==0):

        imnorm = 0.0

    else:

        imnorm = cardainterb/cardaunionb

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S4c
#==============================================================================
def fuzzynorms4c(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    acomplement                     = fuzzycomplement(im1_gray)
    bcomplement                     = fuzzycomplement(im2_gray)
    acomplementinterbcomplement     = fuzzyintersect(acomplement,bcomplement)
    acomplementunionbcomplement     = fuzzyunion(acomplement,bcomplement)
    cardacomplementinterbcomplement = fuzzycardn(acomplementinterbcomplement)
    cardacomplementunionbcomplement = fuzzycardn(acomplementunionbcomplement)

    if(cardacomplementunionbcomplement==0):

        imnorm = 0.0

    else:

        imnorm = cardacomplementinterbcomplement/cardacomplementunionbcomplement

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S5
#==============================================================================
def fuzzynorms5(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    ainterb     = fuzzyintersect(im_fuzzy1,im_fuzzy2)
    cardainterb = fuzzycardn(ainterb)
    carda       = fuzzycardn(im_fuzzy1)
    cardb       = fuzzycardn(im_fuzzy2)

    if(max(carda,cardb)==0):

        imnorm = 0.0

    else:

        imnorm = cardainterb/max(carda,cardb)
            
    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S5c
#==============================================================================
def fuzzynorms5c(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    acomplement                     = fuzzycomplement(im_fuzzy1)
    bcomplement                     = fuzzycomplement(im_fuzzy2)
    acomplementinterbcomplement     = fuzzyintersect(acomplement,bcomplement)
    cardacomplementinterbcomplement = fuzzycardn(acomplementinterbcomplement)
    cardacomplement                 = fuzzycardn(acomplement)
    cardbcomplement                 = fuzzycardn(bcomplement)

    if(max(cardacomplement,cardbcomplement)==0):

        imnorm = 1.0

    else:

        imnorm = cardacomplementinterbcomplement/max(cardacomplement,cardbcomplement)

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S6
#==============================================================================
def fuzzynorms6(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    ainterb     = fuzzyintersect(im_fuzzy1,im_fuzzy2)
    cardainterb = fuzzycardn(ainterb)
    carda       = fuzzycardn(im_fuzzy1)
    cardb       = fuzzycardn(im_fuzzy2)

    if(min(carda,cardb)==0):

        imnorm = 0.0

    else:

        imnorm = cardainterb/min(carda,cardb)

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S6c
#==============================================================================
def fuzzynorms6c(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    acomplement                     = fuzzycomplement(im_fuzzy1)
    bcomplement                     = fuzzycomplement(im_fuzzy2)
    acomplementinterbcomplement     = fuzzyintersect(acomplement,bcomplement)
    cardacomplementinterbcomplement = fuzzycardn(acomplementinterbcomplement)
    cardacomplement                 = fuzzycardn(acomplement)
    cardbcomplement                 = fuzzycardn(bcomplement)

    if(min(cardacomplement,cardbcomplement)==0):

        imnorm = 1.0

    else:

        imnorm = cardacomplementinterbcomplement/min(cardacomplement,cardbcomplement)

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S7
#==============================================================================
def fuzzynorms7(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    atrib               = fuzzyatrib(im_fuzzy1,im_fuzzy2)
    atribcomplement     = fuzzycomplement(atrib)
    cardatribcomplement = fuzzycardn(atribcomplement) 

    abarb               = fuzzyabarb(im_fuzzy1,im_fuzzy2)
    bbara               = fuzzyabarb(im_fuzzy2,im_fuzzy1)
    abarbcomplement     = fuzzycomplement(abarb)
    bbaracomplement     = fuzzycomplement(bbara)
    cardabarbcomplement = fuzzycardn(abarbcomplement)
    cardbbaracomplement = fuzzycardn(bbaracomplement)

    if(max(cardabarbcomplement,cardbbaracomplement)==0):

        imnorm = 0.0

    else:

        imnorm = cardatribcomplement/max(cardabarbcomplement,cardbbaracomplement)

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S8
#==============================================================================
def fuzzynorms8(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    atrib               = fuzzyatrib(im_fuzzy1,im_fuzzy2)
    atribcomplement     = fuzzycomplement(atrib)
    cardatribcomplement = fuzzycardn(atribcomplement) 

    abarb               = fuzzyabarb(im_fuzzy1,im_fuzzy2)
    bbara               = fuzzyabarb(im_fuzzy2,im_fuzzy1)
    abarbcomplement     = fuzzycomplement(abarb)
    bbaracomplement     = fuzzycomplement(bbara)
    cardabarbcomplement = fuzzycardn(abarbcomplement)
    cardbbaracomplement = fuzzycardn(bbaracomplement)

    if(min(cardabarbcomplement,cardbbaracomplement)==0):

        imnorm = 0.0

    else:

        imnorm = cardatribcomplement/min(cardabarbcomplement,cardbbaracomplement)

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S9
#==============================================================================
def fuzzynorms9(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    mval   = im1_gray.shape[0]
    nval   = im1_gray.shape[1]
    atrib  = fuzzyatrib(im_fuzzy1,im_fuzzy2)
    
    imnorm = 1 - ((1/(mval*nval))*np.sum(atrib))

    return imnorm 
#==============================================================================

#==============================================================================
# Fuzzy Norm S10
#==============================================================================
def fuzzynorms10(im1_gray,im2_gray,fuzzymember):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    mval   = im1_gray.shape[0]
    nval   = im1_gray.shape[1]

    val1 = np.fmin(im_fuzzy1,im_fuzzy2)
    val2 = np.fmax(im_fuzzy1,im_fuzzy2)

    val3 = val1/val2
    val3[np.isnan(val3)==True] = 1.0

    imnorm = (1/(mval*nval))*np.sum(val3)

    return imnorm 
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
# Fuzzy Norm S11
#==============================================================================
def fuzzynorms11(im1_gray,im2_gray,fuzzymember,fuzzylocalnorm,sqsize):

    if(fuzzymember==0):

        im_fuzzy1 = fuzzymember0(im1_gray)
        im_fuzzy2 = fuzzymember0(im2_gray)

    elif(fuzzymember==1):

        im_fuzzy1 = fuzzymember1(im1_gray)
        im_fuzzy2 = fuzzymember1(im2_gray)

    elif(fuzzymember==2):

        im_fuzzy1 = fuzzymember2(im1_gray)
        im_fuzzy2 = fuzzymember2(im2_gray)

    nlin = im1_gray.shape[0]
    ncol = im1_gray.shape[1]

    nlinsqz = int(np.floor(nlin/sqsize) + 1)
    ncolsqz = int(np.floor(ncol/sqsize) + 1)

    nsquare = nlinsqz*ncolsqz
    
    nlinnew = sqsize*nlinsqz
    ncolnew = sqsize*ncolsqz

    nlindif = nlinnew - nlin
    ncoldif = ncolnew - ncol
    
    extendtype = 2

    if(extendtype==1):

        im_gray_extend1  = np.pad(im1_gray,((0,nlindif),(0,ncoldif)),'edge')
        im_gray_extend2  = np.pad(im2_gray,((0,nlindif),(0,ncoldif)),'edge')
        im_fuzzy_extend1 = np.pad(im_fuzzy1,((0,nlindif),(0,ncoldif)),'edge')
        im_fuzzy_extend2 = np.pad(im_fuzzy2,((0,nlindif),(0,ncoldif)),'edge')

    if(extendtype==2):

        xnew1 = np.linspace(0,nlin-1,nlin)
        ynew1 = np.linspace(0,ncol-1,ncol)

        xnew2 = np.linspace(0,nlinnew-1,nlinnew)
        ynew2 = np.linspace(0,ncolnew-1,ncolnew)

        im_gray1_interp1 = RectBivariateSpline(xnew1,ynew1,im1_gray.astype(float))
        im_gray_extend1  = im_gray1_interp1(xnew2,ynew2)

        im_gray2_interp2 = RectBivariateSpline(xnew1,ynew1,im2_gray.astype(float))
        im_gray_extend2  = im_gray2_interp2(xnew2,ynew2)

        im_fuzzy1_interp1 = RectBivariateSpline(xnew1,ynew1,im_fuzzy1)
        im_fuzzy_extend1  = im_fuzzy1_interp1(xnew2,ynew2)

        im_fuzzy2_interp2 = RectBivariateSpline(xnew1,ynew1,im_fuzzy2)
        im_fuzzy_extend2  = im_fuzzy2_interp2(xnew2,ynew2)

    imnorm = 0.0

    contx1 = 0
    contx2 = sqsize

    for k1 in range(0,nlinsqz):

        conty1 = 0
        conty2 = sqsize
        
        for k2 in range(0,ncolsqz):

            lval1 = []
            lval2 = []

            lval3 = []
            lval4 = []

            for k3 in range(contx1,contx2):

                for k4 in range(conty1,conty2):

                    lval1.append(im_fuzzy_extend1[k3,k4])
                    lval2.append(im_fuzzy_extend2[k3,k4])

                    lval3.append(im_gray_extend1[k3,k4])
                    lval4.append(im_gray_extend2[k3,k4])
               
            lval3max = max(lval3)
            lval3min = min(lval3)
            
            lval4max = max(lval4)
            lval4min = min(lval4)

            param    = 1/2
            ha       = simfun(lval3max,lval3min,param)
            hb       = simfun(lval4max,lval4min,param)
            weight   = simfun(ha,hb,param)

            locfuzzy1 = np.array(lval1).reshape(sqsize,sqsize)
            locfuzzy2 = np.array(lval2).reshape(sqsize,sqsize)

            if(fuzzylocalnorm==0):    fuzzynorm = fuzzynorms1(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==1):  fuzzynorm = fuzzynorms2(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==2):  fuzzynorm = fuzzynorms3(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==3):  fuzzynorm = fuzzynorms4(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==4):  fuzzynorm = fuzzynorms4c(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==5):  fuzzynorm = fuzzynorms5(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==6):  fuzzynorm = fuzzynorms5c(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==7):  fuzzynorm = fuzzynorms6(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==8):  fuzzynorm = fuzzynorms6c(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==9): fuzzynorm = fuzzynorms7(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==10): fuzzynorm = fuzzynorms8(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==11): fuzzynorm = fuzzynorms9(locfuzzy1,locfuzzy2,fuzzymember)
            elif(fuzzylocalnorm==12): fuzzynorm = fuzzynorms10(locfuzzy1,locfuzzy2,fuzzymember)

            conty1 = conty1 + sqsize
            conty2 = conty2 + sqsize

            imnorm = imnorm + weight*fuzzynorm

        contx1 = contx1 + sqsize
        contx2 = contx2 + sqsize

    imnorm = (1/nsquare)*imnorm

    return imnorm 
#==============================================================================

#==============================================================================
# IM Test
#==============================================================================
# im1       = cv2.imread('../testresults/im_ref1.png')
# im2       = cv2.imread('../testresults/im_num1.png')
# im1_gray  = cv2.cvtColor(im1, cv2.COLOR_BGR2GRAY)
# im2_gray  = cv2.cvtColor(im2, cv2.COLOR_BGR2GRAY)
#==============================================================================

#==============================================================================
# vnorm    = ['1','2','3','4','4c','5','5c','6','6c','7','8','9','10']
# nvnorm   = len(vnorm)
# vmember  = [0,1,2]
# nvmember = len(vmember)

# vmember0 = []
# vmember1 = []
# vmember2 = []

# im1_gray = im2_gray.copy()

# for k1 in range(0,nvnorm):

#     for k2 in range(0,nvmember):

#         fuzzylocalnorm_name = vnorm[k1]  
#         fuzzylocalnorm      = k1

#         fuzzymember_name    = vmember[k2]
#         fuzzymember         = k2

#         sqsize              = 4

#         imnorms11           = fuzzynorms11(im1_gray,im2_gray,fuzzymember,fuzzylocalnorm,sqsize)

#         if(fuzzymember==0): vmember0.append('NORM %s = %f'%(fuzzylocalnorm_name,imnorms11))
#         if(fuzzymember==1): vmember1.append('NORM %s = %f'%(fuzzylocalnorm_name,imnorms11))
#         if(fuzzymember==2): vmember2.append('NORM %s = %f'%(fuzzylocalnorm_name,imnorms11))

#         if(imnorms11<0.95):

#             print('')
#             print('IM Fuzzy Norm Segmented - Membership %s - %d x %d square size - Fuzzy Local Norm %s'%(fuzzymember_name,sqsize,sqsize,fuzzylocalnorm_name))       
#             print(imnorms11)
#             print('')

#         if(fuzzylocalnorm==0):    fuzzynormbase = fuzzynorms1(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==1):  fuzzynormbase = fuzzynorms2(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==2):  fuzzynormbase = fuzzynorms3(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==3):  fuzzynormbase = fuzzynorms4(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==4):  fuzzynormbase = fuzzynorms4c(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==5):  fuzzynormbase = fuzzynorms5(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==6):  fuzzynormbase = fuzzynorms5c(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==7):  fuzzynormbase = fuzzynorms6(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==8):  fuzzynormbase = fuzzynorms6c(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==9): fuzzynormbase = fuzzynorms7(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==10): fuzzynormbase = fuzzynorms8(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==11): fuzzynormbase = fuzzynorms9(im1_gray,im2_gray,fuzzymember)
#         elif(fuzzylocalnorm==12): fuzzynormbase = fuzzynorms10(im1_gray,im2_gray,fuzzymember)
#         # print('')
#         # print('IM Fuzzy Norm Base - Membership %s - Fuzzy Local Norm %s'%(fuzzymember_name,fuzzylocalnorm_name))
#         # print(fuzzynormbase)
#         # print('')

#         if(fuzzymember==0): vmember0.append('NORM FULL %s = %f'%(fuzzylocalnorm_name,fuzzynormbase))
#         if(fuzzymember==1): vmember1.append('NORM FULL %s = %f'%(fuzzylocalnorm_name,fuzzynormbase))
#         if(fuzzymember==2): vmember2.append('NORM FULL %s = %f'%(fuzzylocalnorm_name,fuzzynormbase))

# print('')
# print('Norms Membership 0')
# print(vmember0)
# print('')

# print('')
# print('Norms Membership 1')
# print(vmember1)
# print('')

# print('')
# print('Norms Membership 2')
# print(vmember2)
# print('')
#==============================================================================

#==============================================================================
# IM Fuzzy Nomr
#==============================================================================
# fuzzymember = 0
# sqsize      = 4

# im2_gray = im1_gray.copy()

# imnorms1  = fuzzynorms1(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,0,sqsize)
# print('')
# print('IM Fuzzy Norm S1')
# print(imnorms1)
# print('IM Fuzzy Norm S1 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms2 = fuzzynorms2(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,1,sqsize)
# print('')
# print('IM Fuzzy Norm S2')
# print(imnorms2)
# print('IM Fuzzy Norm S2 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms3 = fuzzynorms3(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,2,sqsize)
# print('')
# print('IM Fuzzy Norm S3')
# print(imnorms3)
# print('IM Fuzzy Norm S3 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms4 = fuzzynorms4(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,3,sqsize)
# print('')
# print('IM Fuzzy Norm S4')
# print(imnorms4)
# print('IM Fuzzy Norm S4 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms4c = fuzzynorms4c(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,4,sqsize)
# print('')
# print('IM Fuzzy Norm S4c')
# print(imnorms4c)
# print('IM Fuzzy Norm S4c - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms5 = fuzzynorms5(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,5,sqsize)
# print('')
# print('IM Fuzzy Norm S5')
# print(imnorms5)
# print('IM Fuzzy Norm S5 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms5c = fuzzynorms5c(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,6,sqsize)
# print('')
# print('IM Fuzzy Norm S5c')
# print(imnorms5c)
# print('IM Fuzzy Norm S5c - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms6 = fuzzynorms6(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,7,sqsize)
# print('')
# print('IM Fuzzy Norm S6')
# print(imnorms6)
# print('IM Fuzzy Norm S6 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms6c = fuzzynorms6c(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,8,sqsize)
# print('')
# print('IM Fuzzy Norm S6c')
# print(imnorms6c)
# print('IM Fuzzy Norm S6c - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms7 = fuzzynorms7(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,9,sqsize)
# print('')
# print('IM Fuzzy Norm S7')
# print(imnorms7)
# print('IM Fuzzy Norm S7 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms8 = fuzzynorms8(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,10,sqsize)
# print('')
# print('IM Fuzzy Norm S8')
# print(imnorms8)
# print('IM Fuzzy Norm S8 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms9 = fuzzynorms9(im1_gray,im1_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,11,sqsize)
# print('')
# print('IM Fuzzy Norm S9')
# print(imnorms9)
# print('IM Fuzzy Norm S9 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')

# imnorms10 = fuzzynorms10(im1_gray,im2_gray,fuzzymember)
# imnorms11 = fuzzynorms11(im1_gray,im2_gray,fuzzymember,12,sqsize)
# print('')
# print('IM Fuzzy Norm S10')
# print(imnorms10)
# print('IM Fuzzy Norm S10 - Square Size = %d'%sqsize)
# print(imnorms11)
# print('')
#==============================================================================