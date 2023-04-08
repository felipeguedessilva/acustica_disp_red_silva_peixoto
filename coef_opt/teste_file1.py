#============================================================================================================
import numpy as np
import math  as mt
import time  as tm
import scipy.special
from   scipy.integrate     import quad, romberg, fixed_quad
#============================================================================================================

#============================================================================================================
bint     = 2.74
thetaint = 2*np.pi
mvalue   = 10
start1   = tm.time()
#============================================================================================================

#============================================================================================================
for qline in range(1,mt.floor(mvalue/2)+1):
    
    for pline in range(qline,mvalue-qline+1):
                
        for q in range(1,mt.floor(mvalue/2)+1):
            
            for p in range(q,mvalue-q+1):
            
                ploc1 = p
                qloc1 = q
                
                ploc2 = pline
                qloc2 = qline

                for r1 in range(0,mvalue+1):
    
                    for s1 in range(0,r1+1):
            
                        for r2 in range(0,mvalue+1):
                
                            for s2 in range(0,r2+1):
            
                                f1 = lambda beta: ((scipy.special.binom(2*r1,2*s1))*((4*(ploc1**(2*r1-2*s1))*(qloc1**(2*s1))*(beta**(2*r1)))/(np.math.factorial(2*r1)))*((-1)**(r1)))*((scipy.special.binom(2*r2,2*s2))*((4*(ploc2**(2*r2-2*s2))*(qloc2**(2*s2))*(beta**(2*r2)))/(np.math.factorial(2*r2)))*((-1)**(r2)))
                                f2 = lambda theta: (np.cos(theta)**(2*r1-2*s1) + np.sin(theta)**(2*s1))*(np.cos(theta)**(2*r2-2*s2) + np.sin(theta)**(2*s2))
                
                                # i1, error1 = quad(f1,0,bint)
                                # i2, error2 = quad(f2,0,thetaint)
    
                                # i1 = romberg(f1,0,bint)
                                # i2 = romberg(f2,0,thetaint)
        
                                i1, error1 = fixed_quad(f1,0,bint,n=5)
                                i2, error2 = fixed_quad(f2,0,thetaint,n=5)
                                
end1   = tm.time()
tt1    = end1 - start1
print('Tempo Total 1: %f'%tt1)
#============================================================================================================