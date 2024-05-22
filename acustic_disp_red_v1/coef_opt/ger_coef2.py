#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy               as np
import math                as mt
import sys
import time
import scipy.special
from   scipy               import integrate
from   scipy.integrate     import nquad
from   scipy.integrate     import quad
from   scipy.integrate     import dblquad
from   numpy               import linalg as la
from   scipy.linalg        import lu_factor, lu_solve
from   scipy.linalg        import solve
#==============================================================================

#==============================================================================
# Auxiliary Functions
#==============================================================================

#==============================================================================
def fun1(rvalue,theta,p,q):

    phisum = 0
    
    for s in range(0,rvalue+1):
        
        a      = scipy.special.binom(2*rvalue,2*s)
        b      = p**(2*(rvalue-s))
        c      = q**(2*s)
        d      = np.cos(theta)**(2*(rvalue-s))
        e      = np.sin(theta)**(2*s)
        phisum = phisum + a*b*c*d*e

    return phisum
#==============================================================================

#==============================================================================
def fun2(beta,theta,mvalue,mloc):
    
    phisum = 0
    
    for r in range(0,mvalue+1):
        
        a      = (2*mloc**(2*r)*beta**(2*r))/(np.math.factorial(2*r))
        b      = (-1)**(r)
        c      = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
        phisum = phisum + a*b*c
    
    phisum = 4 - phisum
    
    return phisum
#==============================================================================

#==============================================================================
def fun3(beta,theta,mvalue,ploc,qloc):
    
    phisum = 0
    
    for r in range(0,mvalue+1):
        
        for s in range(0,r+1):
        
            a      = scipy.special.binom(2*r,2*s)
            b      = (4*(ploc**(2*r-2*s))*(qloc**(2*s))*(beta**(2*r)))/(np.math.factorial(2*r))
            c      = (-1)**(r)
            d      = np.cos(theta)**(2*r-2*s) + np.sin(theta)**(2*s)
            phisum = phisum + a*b*c*d
    
    phisum = 4 - phisum
        
    return phisum
#==============================================================================

#==============================================================================
def fun4(beta):
    
    a = beta**2
    
    return a
#==============================================================================

#==============================================================================
def fun5(beta,cur):
    
    a = (2-2*np.cos(cur*beta))/(cur**2)
    
    return a
#==============================================================================

#==============================================================================
def fun6(beta,theta,mvalue,mloc1,mloc2):
    
    a = fun2(beta,theta,mvalue,mloc1)
    b = fun2(beta,theta,mvalue,mloc2)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun7(beta,theta,mvalue,ploc,qloc,mloc):
    
    a = fun3(beta,theta,mvalue,ploc,qloc)
    b = fun2(beta,theta,mvalue,mloc)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun8(beta,theta,mvalue,mloc):
    
    a = fun4(beta)
    b = fun2(beta,theta,mvalue,mloc)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun9(beta,theta,mvalue,mloc,ploc,qloc):

    a = fun2(beta,theta,mvalue,mloc)    
    b = fun3(beta,theta,mvalue,ploc,qloc)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun10(beta,theta,mvalue,ploc1,qloc1,ploc2,qloc2):
    
    a = fun3(beta,theta,mvalue,ploc1,qloc1)
    b = fun3(beta,theta,mvalue,ploc2,qloc2)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun11(beta,theta,mvalue,ploc,qloc):
    
    a = fun4(beta)
    b = fun3(beta,theta,mvalue,ploc,qloc)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun12(beta,theta,mvalue,mloc,cur):
    
    a = fun5(beta,cur)
    b = fun2(beta,theta,mvalue,mloc)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun13(beta,theta,mvalue,ploc,qloc,cur):
    
    a = fun5(beta,cur)
    b = fun3(beta,theta,mvalue,ploc,qloc)
    c = a*b
    
    return c
#==============================================================================

#==============================================================================
def fun6new(bint,thetaint,mvalue,mloc1,mloc2):
    
    int2d   = 0
    phisum0 = 0 
    phisum1 = 0
    phisum2 = 0
    phisum3 = 0
    
    #f1         = lambda beta: beta**0
    #f2         = lambda theta: theta**0
    #i1, error1 = quad(f1,0,bint)
    #i2, error2 = quad(f2,0,thetaint)
    #phisum0    = 16*i1*i2 
    
    phisum0 = 16*bint*thetaint 
    
    for r1 in range(0,mvalue+1):
        
        f1 = lambda beta: ((-1)**(r1))*(2*mloc1**(2*r1)*beta**(2*r1))/(np.math.factorial(2*r1))
        f2 = lambda theta: np.cos(theta)**(2*r1) + np.sin(theta)**(2*r1)
        
        i1, error1 = quad(f1,0,bint)
        i2, error2 = quad(f2,0,thetaint)
        
        phisum1 = phisum1 + i1*i2
        
    phisum1 = -4*phisum1    
    
    for r2 in range(0,mvalue+1):
        
        f1 = lambda beta: ((-1)**(r2))*(2*mloc2**(2*r2)*beta**(2*r2))/(np.math.factorial(2*r2))
        f2 = lambda theta: np.cos(theta)**(2*r2) + np.sin(theta)**(2*r2)
        
        i1, error1 = quad(f1,0,bint)
        i2, error2 = quad(f2,0,thetaint)
        
        phisum2 = phisum2 + i1*i2
        
    phisum2 = -4*phisum2
    
    for r1 in range(0,mvalue+1):
        
        for r2 in range(0,mvalue+1):
            
            f1 = lambda beta: (((-1)**(r1))*(2*mloc1**(2*r1)*beta**(2*r1))/(np.math.factorial(2*r1)))*(((-1)**(r2))*(2*mloc2**(2*r2)*beta**(2*r2))/(np.math.factorial(2*r2)))                              
            f2 = lambda theta: (np.cos(theta)**(2*r1) + np.sin(theta)**(2*r1))*(np.cos(theta)**(2*r2) + np.sin(theta)**(2*r2))
        
            i1, error1 = quad(f1,0,bint)
            i2, error2 = quad(f2,0,thetaint)
            
            phisum3 = phisum3 + i1*i2

    int2d = phisum0 + phisum1 + phisum2 + phisum3
        
    return int2d
#==============================================================================

#==============================================================================
def fun7new(bint,thetaint,mvalue,ploc1,qloc1,mloc1):
    
    int2d   = 0
    phisum0 = 0 
    phisum1 = 0
    phisum2 = 0
    phisum3 = 0
    
    #f1         = lambda beta: beta**0
    #f2         = lambda theta: theta**0
    #i1, error1 = quad(f1,0,bint)
    #i2, error2 = quad(f2,0,thetaint)
    #phisum0    = 16*i1*i2 
    
    phisum0 = 16*bint*thetaint 
    
    for r1 in range(0,mvalue+1):
        
        f1 = lambda beta: ((-1)**(r1))*(2*mloc1**(2*r1)*beta**(2*r1))/(np.math.factorial(2*r1))
        f2 = lambda theta: np.cos(theta)**(2*r1) + np.sin(theta)**(2*r1)
        
        i1, error1 = quad(f1,0,bint)
        i2, error2 = quad(f2,0,thetaint)
        
        phisum1 = phisum1 + i1*i2
        
    phisum1 = -4*phisum1    
    
    for r1 in range(0,mvalue+1):
    
        for s1 in range(0,r1+1):
            
            f1 = lambda beta: (scipy.special.binom(2*r1,2*s1))*((4*(ploc1**(2*r1-2*s1))*(qloc1**(2*s1))*(beta**(2*r1)))/(np.math.factorial(2*r1)))*((-1)**(r1))
            f2 = lambda theta:np.cos(theta)**(2*r1-2*s1) + np.sin(theta)**(2*s1)
            
            i1, error1 = quad(f1,0,bint)
            i2, error2 = quad(f2,0,thetaint)
            
            phisum2 = phisum2 + i1*i2

    phisum2 = -4*phisum2

    for r1 in range(0,mvalue+1):
    
        for s1 in range(0,r1+1):
            
            for r2 in range(0,mvalue+1):

                f1 = lambda beta: ((scipy.special.binom(2*r1,2*s1))*((4*(ploc1**(2*r1-2*s1))*(qloc1**(2*s1))*(beta**(2*r1)))/(np.math.factorial(2*r1)))*((-1)**(r1)))*(((-1)**(r2))*(2*mloc1**(2*r2)*beta**(2*r2))/(np.math.factorial(2*r2)))      
                f2 = lambda theta: (np.cos(theta)**(2*r1-2*s1) + np.sin(theta)**(2*s1))*(np.cos(theta)**(2*r2) + np.sin(theta)**(2*r2))
            
                i1, error1 = quad(f1,0,bint)
                i2, error2 = quad(f2,0,thetaint)
            
                phisum3 = phisum3 + i1*i2

    int2d = phisum0 + phisum1 + phisum2 + phisum3
        
    return int2d
#==============================================================================

#==============================================================================
def fun8new(bint,thetaint,mvalue,mloc1):
    
    int2d   = 0
    phisum0 = 0 
    phisum1 = 0
    
    #f1         = lambda beta: beta**2
    #f2         = lambda theta: theta**0
    #i1, error1 = quad(f1,0,bint)
    #i2, error2 = quad(f2,0,thetaint)
    #phisum0    = 4*i1*i2 
    
    phisum0 = 4*(1/3)*(bint**3)*(thetaint) 
    
    for r1 in range(0,mvalue+1):
        
        f1 = lambda beta: ((beta**2))*(((-1)**(r1))*(2*mloc1**(2*r1)*beta**(2*r1))/(np.math.factorial(2*r1)))
        f2 = lambda theta: np.cos(theta)**(2*r1) + np.sin(theta)**(2*r1)
        
        i1, error1 = quad(f1,0,bint)
        i2, error2 = quad(f2,0,thetaint)
        
        phisum1 = phisum1 + i1*i2
    
    phisum1 = -1*phisum1
    
    int2d = phisum0 + phisum1 
    
    return int2d
#==============================================================================

#==============================================================================
def fun10new(bint,thetaint,mvalue,ploc1,qloc1,ploc2,qloc2):
    
    int2d   = 0
    phisum0 = 0 
    phisum1 = 0
    phisum2 = 0
    phisum3 = 0
    
    #f1         = lambda beta: beta**0
    #f2         = lambda theta: theta**0
    #i1, error1 = quad(f1,0,bint)
    #i2, error2 = quad(f2,0,thetaint)
    #phisum0    = 16*i1*i2 
    
    phisum0 = 16*bint*thetaint 
       
    for r1 in range(0,mvalue+1):
    
        for s1 in range(0,r1+1):
            
            f1 = lambda beta: (scipy.special.binom(2*r1,2*s1))*((4*(ploc1**(2*r1-2*s1))*(qloc1**(2*s1))*(beta**(2*r1)))/(np.math.factorial(2*r1)))*((-1)**(r1))
            f2 = lambda theta:np.cos(theta)**(2*r1-2*s1) + np.sin(theta)**(2*s1)
            
            i1, error1 = quad(f1,0,bint)
            i2, error2 = quad(f2,0,thetaint)
            
            phisum1 = phisum1 + i1*i2

    phisum1 = -4*phisum1
    
    for r2 in range(0,mvalue+1):
    
        for s2 in range(0,r2+1):
            
            f1 = lambda beta: (scipy.special.binom(2*r2,2*s2))*((4*(ploc2**(2*r2-2*s2))*(qloc2**(2*s2))*(beta**(2*r2)))/(np.math.factorial(2*r2)))*((-1)**(r2))
            f2 = lambda theta:np.cos(theta)**(2*r2-2*s2) + np.sin(theta)**(2*s2)
            
            i1, error1 = quad(f1,0,bint)
            i2, error2 = quad(f2,0,thetaint)
            
            phisum2 = phisum2 + i1*i2

    phisum2 = -4*phisum2
    
    for r1 in range(0,mvalue+1):
    
        for s1 in range(0,r1+1):
            
            for r2 in range(0,mvalue+1):
                
                for s2 in range(0,r2+1):
            
                    f1 = lambda beta: ((scipy.special.binom(2*r1,2*s1))*((4*(ploc1**(2*r1-2*s1))*(qloc1**(2*s1))*(beta**(2*r1)))/(np.math.factorial(2*r1)))*((-1)**(r1)))*((scipy.special.binom(2*r2,2*s2))*((4*(ploc2**(2*r2-2*s2))*(qloc2**(2*s2))*(beta**(2*r2)))/(np.math.factorial(2*r2)))*((-1)**(r2)))
                    f2 = lambda theta: (np.cos(theta)**(2*r1-2*s1) + np.sin(theta)**(2*s1))*(np.cos(theta)**(2*r2-2*s2) + np.sin(theta)**(2*s2))
                
                    i1, error1 = quad(f1,0,bint)
                    i2, error2 = quad(f2,0,thetaint)
            
                    phisum3 = phisum3 + i1*i2

    int2d = phisum0 + phisum1 + phisum2 + phisum3
        
    return int2d
#==============================================================================

#==============================================================================
def fun11new(bint,thetaint,mvalue,ploc1,qloc1):
    
    int2d   = 0
    phisum0 = 0 
    phisum1 = 0
    
    #f1         = lambda beta: beta**2
    #f2         = lambda theta: theta**0
    #i1, error1 = quad(f1,0,bint)
    #i2, error2 = quad(f2,0,thetaint)
    #phisum0    = 4*i1*i2 
    
    phisum0 = 4*(1/3)*(bint**3)*(thetaint) 
    
    for r1 in range(0,mvalue+1):
    
        for s1 in range(0,r1+1):
            
            f1 = lambda beta: (beta**2)*((scipy.special.binom(2*r1,2*s1))*((4*(ploc1**(2*r1-2*s1))*(qloc1**(2*s1))*(beta**(2*r1)))/(np.math.factorial(2*r1)))*((-1)**(r1)))
            f2 = lambda theta:np.cos(theta)**(2*r1-2*s1) + np.sin(theta)**(2*s1)
            
            i1, error1 = quad(f1,0,bint)
            i2, error2 = quad(f2,0,thetaint)
            
            phisum1 = phisum1 + i1*i2

    phisum1 = -1*phisum1
    
    int2d = phisum0 + phisum1 
    
    return int2d
#==============================================================================

#==============================================================================
def fun12new(bint,thetaint,mvalue,mloc1,cur):
    
    int2d   = 0
    phisum0 = 0 
    phisum1 = 0
    
    f1         = lambda beta: (2-2*np.cos(cur*beta))/(cur**2)
    f2         = lambda theta: theta**0
    i1, error1 = quad(f1,0,bint)
    i2, error2 = quad(f2,0,thetaint)
    phisum0    = 4*i1*i2 
     
    for r1 in range(0,mvalue+1):
        
        f1 = lambda beta: ((2-2*np.cos(cur*beta))/(cur**2))*(((-1)**(r1))*(2*mloc1**(2*r1)*beta**(2*r1))/(np.math.factorial(2*r1)))
        f2 = lambda theta: np.cos(theta)**(2*r1) + np.sin(theta)**(2*r1)
        
        i1, error1 = quad(f1,0,bint)
        i2, error2 = quad(f2,0,thetaint)
        
        phisum1 = phisum1 + i1*i2
    
    phisum1 = -1*phisum1
    
    int2d = phisum0 + phisum1 
    
    return int2d
#==============================================================================

#==============================================================================
def fun13new(bint,thetaint,mvalue,ploc1,qloc1,cur):
    
    int2d   = 0
    phisum0 = 0 
    phisum1 = 0
    
    f1         = lambda beta: (2-2*np.cos(cur*beta))/(cur**2)
    f2         = lambda theta: theta**0
    i1, error1 = quad(f1,0,bint)
    i2, error2 = quad(f2,0,thetaint)
    phisum0    = 4*i1*i2 
    
    for r1 in range(0,mvalue+1):
    
        for s1 in range(0,r1+1):
            
            f1 = lambda beta: ((2-2*np.cos(cur*beta))/(cur**2))*((scipy.special.binom(2*r1,2*s1))*((4*(ploc1**(2*r1-2*s1))*(qloc1**(2*s1))*(beta**(2*r1)))/(np.math.factorial(2*r1)))*((-1)**(r1)))
            f2 = lambda theta:np.cos(theta)**(2*r1-2*s1) + np.sin(theta)**(2*s1)
            
            i1, error1 = quad(f1,0,bint)
            i2, error2 = quad(f2,0,thetaint)
            
            phisum1 = phisum1 + i1*i2

    phisum1 = -1*phisum1
    
    int2d = phisum0 + phisum1 
    
    return int2d
#==============================================================================

#==============================================================================
def spatte_cl(mvalue,nvalue,nrow,ncol):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((ncol,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        if(r==1): bsis[prow,0] = 1.0
        
        pcol = 0
            
        for m in range(1,ncol+1):
                
            asis[pcol,prow] = (m)**(2*r)
            pcol            = pcol + 1
            
        prow = prow + 1
            
    try:
    
        csis = np.linalg.solve(asis,bsis)
    
    except:
            
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    return csis
#==============================================================================

#==============================================================================
def specte_cl(mvalue,nvalue,nrow,ncol,theta):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        if(r==1): bsis[prow,0] = 1.0
        
        pcol = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                         

    return csis
#==============================================================================

#==============================================================================
def specte_rb(mvalue,nvalue,nrow,ncol,theta):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        if(r==1): bsis[prow,0] = 1.0
        
        pcol = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(mvalue/2)+1):
            
            for p in range(q,mvalue-q+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                  

    return csis
#==============================================================================

#==============================================================================
def specte_crb(mvalue,nvalue,nrow,ncol,theta):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        if(r==1): bsis[prow,0] = 1.0
        
        pcol = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(nvalue/2)+1):
            
            for p in range(q,nvalue-q+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                  

    return csis
#==============================================================================

#==============================================================================
def specte_sq(mvalue,nvalue,nrow,ncol,theta):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        if(r==1): bsis[prow,0] = 1.0
        
        pcol = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,mvalue+1):
            
            for p in range(q,mvalue+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 
        
    return csis
#==============================================================================

#==============================================================================
def specte_csq(mvalue,nvalue,nrow,ncol,theta):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        if(r==1): bsis[prow,0] = 1.0
        
        pcol = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,nvalue+1):
            
            for p in range(q,mvalue+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                    

    return csis
#==============================================================================

#==============================================================================
def specls_cl(mvalue,nvalue,nrow,ncol,thetaint,bint):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun8,[[0,bint],[0,thetaint]],args=(mvalue,r,))
        res          = fun8new(bint,thetaint,mvalue,r)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
    
    for qline in range(1,0+1):
            
        for pline in range(qline,0+1):
                
            #res, err     = nquad(fun11,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,))
            res          = fun11new(bint,thetaint,mvalue,pline,qline)
            bsis[prow,0] = res 
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,0+1):
            
                for p in range(q,0+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                           

    return csis
#==============================================================================

#==============================================================================
def specls_rb(mvalue,nvalue,nrow,ncol,thetaint,bint):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun8,[[0,bint],[0,thetaint]],args=(mvalue,r,))
        res          = fun8new(bint,thetaint,mvalue,r)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(mvalue/2)+1):
            
            for p in range(q,mvalue-q+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,mt.floor(mvalue/2)+1):
            
        for pline in range(qline,mvalue-qline+1):
            
            #res, err     = nquad(fun11,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,))
            res          = fun11new(bint,thetaint,mvalue,pline,qline)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,mt.floor(mvalue/2)+1):
            
                for p in range(q,mvalue-q+1):
                    
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                            

    return csis
#==============================================================================

#==============================================================================
def specls_crb(mvalue,nvalue,nrow,ncol,thetaint,bint):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun8,[[0,bint],[0,thetaint]],args=(mvalue,r,))
        res          = fun8new(bint,thetaint,mvalue,r)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(nvalue/2)+1):
            
            for p in range(q,nvalue-q+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,mt.floor(nvalue/2)+1):
            
        for pline in range(qline,nvalue-qline+1):
                
            #res, err     = nquad(fun11,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,))
            res          = fun11new(bint,thetaint,mvalue,pline,qline)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,mt.floor(nvalue/2)+1):
            
                for p in range(q,nvalue-q+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                   

    return csis
#==============================================================================

#==============================================================================
def specls_sq(mvalue,nvalue,nrow,ncol,thetaint,bint):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun8,[[0,bint],[0,thetaint]],args=(mvalue,r,))
        res          = fun8new(bint,thetaint,mvalue,r)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,mvalue+1):
            
            for p in range(q,mvalue+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,mvalue+1):
            
        for pline in range(qline,mvalue+1):
                
            #res, err     = nquad(fun11,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,))
            res          = fun11new(bint,thetaint,mvalue,pline,qline)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,mvalue+1):
            
                for p in range(q,mvalue+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                  
        
    return csis
#==============================================================================

#==============================================================================
def specls_csq(mvalue,nvalue,nrow,ncol,thetaint,bint):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun8,[[0,bint],[0,thetaint]],args=(mvalue,r,))
        res          = fun8new(bint,thetaint,mvalue,r)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,nvalue+1):
            
            for p in range(q,mvalue+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,nvalue+1):
            
        for pline in range(qline,mvalue+1):
                
            #res, err     = nquad(fun11,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,))
            res          = fun11new(bint,thetaint,mvalue,pline,qline)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,nvalue+1):
            
                for p in range(q,mvalue+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    return csis
#==============================================================================

#==============================================================================
def dispte_cl(mvalue,nvalue,nrow,ncol,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            asis[prow,pcol] = a
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                asis[prow,pcol] = (p)**(2*r)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    for r in range(1,mvalue+1):
        
        for s in range(1,int(r/2)+1):
            
            a            = 0.5
            b            = (np.math.factorial(r))*(np.math.factorial(2*s))*(np.math.factorial(2*r-2*s))
            c            = (np.math.factorial(2*r))*(np.math.factorial(r-s))*(np.math.factorial(s))
            d            = cur**(2*r-2)
            bsis[prow,0] = a*(b/c)*d
            pcol         = mvalue
              
            for q in range(1,0+1):
            
                for p in range(q,0+1):
                    
                    a               = (p)**(2*r-2*s)
                    b               = (q)**(2*s)
                    asis[prow,pcol] = a*b
                    pcol            = pcol + 1
        
            prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                             

    return csis
#==============================================================================

#==============================================================================
def dispte_rb(mvalue,nvalue,nrow,ncol,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            asis[prow,pcol] = a
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                asis[prow,pcol] = (p)**(2*r)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    for r in range(1,mvalue+1):
        
        for s in range(1,int(r/2)+1):
            
            a            = 0.5
            b            = (np.math.factorial(r))*(np.math.factorial(2*s))*(np.math.factorial(2*r-2*s))
            c            = (np.math.factorial(2*r))*(np.math.factorial(r-s))*(np.math.factorial(s))
            d            = cur**(2*r-2)
            bsis[prow,0] = a*(b/c)*d
            pcol         = mvalue
              
            for q in range(1,mt.floor(mvalue/2)+1):
            
                for p in range(q,mvalue-q+1):
                    
                    a               = (p)**(2*r-2*s)
                    b               = (q)**(2*s)
                    asis[prow,pcol] = a*b
                    pcol            = pcol + 1
        
            prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                              
 

    return csis
#==============================================================================

#==============================================================================
def dispte_crb(mvalue,nvalue,nrow,ncol,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            asis[prow,pcol] = a
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                asis[prow,pcol] = (p)**(2*r)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    for r in range(1,mvalue+1):
        
        for s in range(1,int(r/2)+1):
            
            a            = 0.5
            b            = (np.math.factorial(r))*(np.math.factorial(2*s))*(np.math.factorial(2*r-2*s))
            c            = (np.math.factorial(2*r))*(np.math.factorial(r-s))*(np.math.factorial(s))
            d            = cur**(2*r-2)
            bsis[prow,0] = a*(b/c)*d
            pcol         = mvalue
              
            for q in range(1,mt.floor(nvalue/2)+1):
            
                for p in range(q,nvalue-q+1):
                    
                    a               = (p)**(2*r-2*s)
                    b               = (q)**(2*s)
                    asis[prow,pcol] = a*b
                    pcol            = pcol + 1
        
            prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                           
    
    return csis
#==============================================================================

#==============================================================================
def dispte_sq(mvalue,nvalue,nrow,ncol,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            asis[prow,pcol] = a
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                asis[prow,pcol] = (p)**(2*r)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    for r in range(1,mvalue+1):
        
        for s in range(1,int(r/2)+1):
            
            a            = 0.5
            b            = (np.math.factorial(r))*(np.math.factorial(2*s))*(np.math.factorial(2*r-2*s))
            c            = (np.math.factorial(2*r))*(np.math.factorial(r-s))*(np.math.factorial(s))
            d            = cur**(2*r-2)
            bsis[prow,0] = a*(b/c)*d
            pcol         = mvalue
              
            for q in range(1,mvalue+1):
            
                for p in range(q,mvalue+1):
                    
                    a               = (p)**(2*r-2*s)
                    b               = (q)**(2*s)
                    asis[prow,pcol] = a*b
                    pcol            = pcol + 1
        
            prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                        
        
    return csis
#==============================================================================

#==============================================================================
def dispte_csq(mvalue,nvalue,nrow,ncol,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            asis[prow,pcol] = a
            pcol            = pcol + 1
        
        for q in range(1,nvalue+1):
            
            for p in range(q,mvalue+1):
                
                asis[prow,pcol] = (p)**(2*r)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    for r in range(1,mvalue+1):
        
        for s in range(1,int(r/2)+1):
            
            a            = 0.5
            b            = (np.math.factorial(r))*(np.math.factorial(2*s))*(np.math.factorial(2*r-2*s))
            c            = (np.math.factorial(2*r))*(np.math.factorial(r-s))*(np.math.factorial(s))
            d            = cur**(2*r-2)
            bsis[prow,0] = a*(b/c)*d
            pcol         = mvalue

            for q in range(1,nvalue+1):
            
                for p in range(q,mvalue+1):
                    
                    a               = (p)**(2*r-2*s)
                    b               = (q)**(2*s)
                    asis[prow,pcol] = a*b
                    pcol            = pcol + 1
        
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    return csis
#==============================================================================

#==============================================================================
def spectetheta_cl(mvalue,nvalue,nrow,ncol,theta,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 
         

    return csis
#==============================================================================

#==============================================================================
def spectetheta_rb(mvalue,nvalue,nrow,ncol,theta,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(mvalue/2)+1):
            
            for p in range(q,mvalue-q+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                    

    return csis
#==============================================================================

#==============================================================================
def spectetheta_crb(mvalue,nvalue,nrow,ncol,theta,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(nvalue/2)+1):
            
            for p in range(q,nvalue-q+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                   

    return csis
#==============================================================================

#==============================================================================
def spectetheta_sq(mvalue,nvalue,nrow,ncol,theta,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,mvalue+1):
            
            for p in range(q,mvalue+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    return csis
#==============================================================================

#==============================================================================
def spectetheta_csq(mvalue,nvalue,nrow,ncol,theta,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,nrow+1):
            
        bsis[prow,0] = cur**(2*r-2)
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            a               = (m)**(2*r)
            b               = np.cos(theta)**(2*r) + np.sin(theta)**(2*r)
            asis[prow,pcol] = a*b
            pcol            = pcol + 1
        
        for q in range(1,nvalue+1):
            
            for p in range(q,mvalue+1):
                
                asis[prow,pcol] = 2*fun1(pcol,theta,p,q)
                pcol            = pcol + 1
        
        prow = prow + 1
    
    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                      

    return csis
#==============================================================================

#==============================================================================
def displs_cl(mvalue,nvalue,nrow,ncol,thetaint,bint,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun12,[[0,bint],[0,thetaint]],args=(mvalue,r,cur,))
        res          = fun12new(bint,thetaint,mvalue,r,cur)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,0+1):
            
            for p in range(q,0+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
    
    for qline in range(1,0+1):
            
        for pline in range(qline,0+1):
                
            #res, err     = nquad(fun13,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,cur))
            res          = fun13new(bint,thetaint,mvalue,pline,qline,cur)
            bsis[prow,0] = res 
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,0+1):
            
                for p in range(q,0+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                          

    return csis
#==============================================================================

#==============================================================================
def displs_rb(mvalue,nvalue,nrow,ncol,thetaint,bint,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0

    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun12,[[0,bint],[0,thetaint]],args=(mvalue,r,cur,))
        res          = fun12new(bint,thetaint,mvalue,r,cur)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(mvalue/2)+1):
            
            for p in range(q,mvalue-q+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,mt.floor(mvalue/2)+1):
            
        for pline in range(qline,mvalue-qline+1):
                
            #res, err     = nquad(fun13,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,cur,))
            res          = fun13new(bint,thetaint,mvalue,pline,qline,cur)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,mt.floor(mvalue/2)+1):
            
                for p in range(q,mvalue-q+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                             

    return csis
#==============================================================================

#==============================================================================
def displs_crb(mvalue,nvalue,nrow,ncol,thetaint,bint,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun12,[[0,bint],[0,thetaint]],args=(mvalue,r,cur,))
        res          = fun12new(bint,thetaint,mvalue,r,cur)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,mt.floor(nvalue/2)+1):
            
            for p in range(q,nvalue-q+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,mt.floor(nvalue/2)+1):
            
        for pline in range(qline,nvalue-qline+1):
                
            #res, err     = nquad(fun13,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,cur,))
            res          = fun13new(bint,thetaint,mvalue,pline,qline,cur)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,mt.floor(nvalue/2)+1):
            
                for p in range(q,nvalue-q+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                    

    return csis
#==============================================================================

#==============================================================================
def displs_sq(mvalue,nvalue,nrow,ncol,thetaint,bint,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun12,[[0,bint],[0,thetaint]],args=(mvalue,r,cur,))
        res          = fun12new(bint,thetaint,mvalue,r,cur)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,mvalue+1):
            
            for p in range(q,mvalue+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,mvalue+1):
            
        for pline in range(qline,mvalue+1):
                
            #res, err     = nquad(fun13,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,cur,))
            res          = fun13new(bint,thetaint,mvalue,pline,qline,cur)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,mvalue+1):
            
                for p in range(q,mvalue+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 
        
    return csis
#==============================================================================

#==============================================================================
def displs_csq(mvalue,nvalue,nrow,ncol,thetaint,bint,cur):
    
    asis   = np.zeros((nrow,ncol))
    bsis   = np.zeros((nrow,1))
    prow   = 0
    
    for r in range(1,mvalue+1):
        
        #res, err     = nquad(fun12,[[0,bint],[0,thetaint]],args=(mvalue,r,cur,))
        res          = fun12new(bint,thetaint,mvalue,r,cur)
        bsis[prow,0] = res
        pcol         = 0
        
        for m in range(1,mvalue+1):
            
            #res, err        = nquad(fun6,[[0,bint],[0,thetaint]],args=(mvalue,m,r,))
            res             = fun6new(bint,thetaint,mvalue,m,r)
            asis[prow,pcol] = res
            pcol            = pcol + 1
        
        for q in range(1,nvalue+1):
            
            for p in range(q,mvalue+1):
                
                #res, err        = nquad(fun7,[[0,bint],[0,thetaint]],args=(mvalue,p,q,r,))
                res             = fun7new(bint,thetaint,mvalue,p,q,r)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
        prow = prow + 1
            
    for qline in range(1,nvalue+1):
            
        for pline in range(qline,mvalue+1):
                
            #res, err     = nquad(fun13,[[0,bint],[0,thetaint]],args=(mvalue,pline,qline,cur))
            res          = fun13new(bint,thetaint,mvalue,pline,qline,cur)
            bsis[prow,0] = res
            pcol         = 0
        
            for m in range(1,mvalue+1):
            
                #res, err        = nquad(fun9,[[0,bint],[0,thetaint]],args=(mvalue,m,pline,qline,))
                res             = fun7new(bint,thetaint,mvalue,pline,qline,m)
                asis[prow,pcol] = res
                pcol            = pcol + 1
        
            for q in range(1,nvalue+1):
            
                for p in range(q,mvalue+1):
                
                    #res, err        = nquad(fun10,[[0,bint],[0,thetaint]],args=(mvalue,p,q,pline,qline,))
                    res             = fun10new(bint,thetaint,mvalue,p,q,pline,qline)
                    asis[prow,pcol] = res
                    pcol            = pcol + 1
                
            prow = prow + 1

    if(nrow==ncol):
    
        try:
    
            csis = np.linalg.solve(asis,bsis)
    
        except:
            
            csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    else:
        
        csis, rcsis, rankcsis, singcsis = la.lstsq(asis,bsis,rcond=None)                 

    return csis
#==============================================================================

#==============================================================================
def stencil_distribute(csis,mvalue,nvalue,norder,listcoef):

    mcoef                = np.zeros((norder+1,norder+1))
    mcoef[mvalue,mvalue] = -4*np.sum(csis)
    nlistcoef            = len(listcoef)
    
    for k1 in range(0,nlistcoef):
        
        p = listcoef[k1][0]
        q = listcoef[k1][1]
        
        mcoef[mvalue+q,mvalue+p] = csis[k1,0]
        mcoef[mvalue+p,mvalue+q] = mcoef[mvalue+q,mvalue+p]
        
        mcoef[mvalue+q,mvalue-p] = mcoef[mvalue+q,mvalue+p]
        mcoef[mvalue+p,mvalue-q] = mcoef[mvalue+q,mvalue-p]
        
        mcoef[mvalue-q,mvalue+p] = mcoef[mvalue+q,mvalue+p]
        mcoef[mvalue-p,mvalue+q] = mcoef[mvalue-q,mvalue+p]
        
        mcoef[mvalue-q,mvalue-p] = mcoef[mvalue+q,mvalue+p]
        mcoef[mvalue-p,mvalue-q] = mcoef[mvalue-q,mvalue-p]
    
    return mcoef
#==============================================================================

#==============================================================================
# System Construction
#==============================================================================
def calccoef(method,shape,mvalue,nvalue,cur):
#==============================================================================

#==============================================================================    
    norder    = 2*mvalue
    theta     = np.pi/8
    bint      = 2.74
    thetaint  = 2*np.pi
    
    if(shape=='cl'):
        
        nvalue = 1
    
    elif(shape=='rb'):
    
        nvalue = mvalue
#==============================================================================

#==============================================================================
    listcoef = []
    
    for p in range(1,mvalue+1):
        
        a = p
        b = 0
        pair = (a,b)
        listcoef.append(pair)
    
    if(shape=='cl'):
 
        for q in range(1,0+1):
        
            for p in range(q,0+1):
            
                a = p
                b = q
                pair = (a,b)
                listcoef.append(pair)

    if(shape=='rb'):
 
        for q in range(1,mt.floor(mvalue/2)+1):
        
            for p in range(q,mvalue-q+1):
            
                a = p
                b = q
                pair = (a,b)
                listcoef.append(pair)

    if(shape=='crb'):
 
        for q in range(1,mt.floor(nvalue/2)+1):
        
            for p in range(q,nvalue-q+1):
            
                a = p
                b = q
                pair = (a,b)
                listcoef.append(pair)
   
    if(shape=='sq'):

        for q in range(1,mvalue+1):
        
            for p in range(q,mvalue+1):
            
                a = p
                b = q
                pair = (a,b)
                listcoef.append(pair)
  
    if(shape=='csq'):

        for q in range(1,nvalue+1):
        
            for p in range(q,mvalue+1):
            
                a = p
                b = q
                pair = (a,b)
                listcoef.append(pair)
#==============================================================================

#==============================================================================        
    if(method=='spatte' and shape=='cl'):
        
        nrow   = mvalue 
        ncol   = len(listcoef)
        csis   = spatte_cl(mvalue,nvalue,nrow,ncol)
        mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
#==============================================================================

#==============================================================================        
    elif(method=='specte'):
        
        nrow = mvalue 
        ncol = len(listcoef)
        
        if(shape=='cl'):
            
            csis   = specte_cl(mvalue,nvalue,nrow,ncol,theta)      
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)

        elif(shape=='rb'):

            csis   = specte_rb(mvalue,nvalue,nrow,ncol,theta)
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='crb'):

            csis  = specte_crb(mvalue,nvalue,nrow,ncol,theta)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='sq'):
        
            csis  = specte_sq(mvalue,nvalue,nrow,ncol,theta)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='csq'):
        
            csis  = specte_csq(mvalue,nvalue,nrow,ncol,theta)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
#==============================================================================

#==============================================================================        
    elif(method=='specls'):
        
        nrow = len(listcoef) 
        ncol = len(listcoef)
        
        if(shape=='cl'):
        
            csis   = specls_cl(mvalue,nvalue,nrow,ncol,thetaint,bint)      
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)

        elif(shape=='rb'):

            csis   = specls_rb(mvalue,nvalue,nrow,ncol,thetaint,bint)
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='crb'):
             
            csis  = specls_crb(mvalue,nvalue,nrow,ncol,thetaint,bint)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='sq'):
        
            csis  = specls_sq(mvalue,nvalue,nrow,ncol,thetaint,bint)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='csq'):
        
            csis  = specls_csq(mvalue,nvalue,nrow,ncol,thetaint,bint)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
#==============================================================================

#==============================================================================        
    elif(method=='dispte'):
        
        nrow = mvalue
        
        for r in range(1,mvalue+1):
            
            for s in range(1,int(r/2)+1):
                
                nrow = nrow + 1
        
        ncol = len(listcoef)
                        
        if(shape=='cl'):
            
            csis   = dispte_cl(mvalue,nvalue,nrow,ncol,cur)      
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)

        elif(shape=='rb'):

            csis   = dispte_rb(mvalue,nvalue,nrow,ncol,cur)
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='crb'):

            csis  = dispte_crb(mvalue,nvalue,nrow,ncol,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='sq'):
        
            csis  = dispte_sq(mvalue,nvalue,nrow,ncol,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='csq'):
        
            csis  = dispte_csq(mvalue,nvalue,nrow,ncol,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
#==============================================================================

#==============================================================================        
    elif(method=='spectetheta'):
        
        nrow = mvalue 
        ncol = len(listcoef)
        
        if(shape=='cl'):
            
            csis   = spectetheta_cl(mvalue,nvalue,nrow,ncol,theta,cur)      
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)

        elif(shape=='rb'):

            csis   = spectetheta_rb(mvalue,nvalue,nrow,ncol,theta,cur)
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='crb'):

            csis  = spectetheta_crb(mvalue,nvalue,nrow,ncol,theta,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='sq'):
        
            csis  = spectetheta_sq(mvalue,nvalue,nrow,ncol,theta,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='csq'):
        
            csis  = spectetheta_csq(mvalue,nvalue,nrow,ncol,theta,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
#==============================================================================

#==============================================================================        
    elif(method=='displs'):
        
        nrow = len(listcoef) 
        ncol = len(listcoef)
        
        if(shape=='cl'):
        
            csis   = displs_cl(mvalue,nvalue,nrow,ncol,thetaint,bint,cur)      
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)

        elif(shape=='rb'):

            csis   = displs_rb(mvalue,nvalue,nrow,ncol,thetaint,bint,cur)
            mcoef  = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='crb'):
             
            csis  = displs_crb(mvalue,nvalue,nrow,ncol,thetaint,bint,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='sq'):
        
            csis  = displs_sq(mvalue,nvalue,nrow,ncol,thetaint,bint,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
        
        elif(shape=='csq'):
        
            csis  = displs_csq(mvalue,nvalue,nrow,ncol,thetaint,bint,cur)
            mcoef = stencil_distribute(csis,mvalue,nvalue,norder,listcoef)
#==============================================================================

#==============================================================================
    else: print('Configuration Not Possible!')
#==============================================================================

#==============================================================================
    return mcoef
#==============================================================================