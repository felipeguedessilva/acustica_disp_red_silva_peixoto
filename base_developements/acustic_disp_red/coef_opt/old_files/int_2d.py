#==============================================================================
import numpy as np
#==============================================================================

#==============================================================================
def pontopeso(npt):

    pontos = np.zeros(npt)
    pesos  = np.zeros(npt)

    if(npt==1):
        
        pontos[0] = 0.0;
        pesos[0]  = 2.0;

    elif(npt==2):
     
        pontos[0] = -0.577350269189626
        pontos[1] = -pontos[0]
        pesos[0]  =  1.0
        pesos[1]  =  pesos[0]
    
    elif(npt==3):
     
        pontos[0]  = -0.774596669241483
        pontos[1]  =  0.0
        pontos[2]  = -pontos[0]
        pesos[0]   = 0.555555555555556
        pesos[1]   = 0.888888888888889
        pesos[2]   = pesos[0]

    elif(npt==4):
     
        pontos[0]  = -0.861136311594053
        pontos[1]  = -0.339981043584856
        pontos[2]  = -pontos[1]
        pontos[3]  = -pontos[0]
        pesos[0]   = 0.347854845137454
        pesos[1]   = 0.652145154862546
        pesos[2]   = pesos[1]
        pesos[3]   = pesos[0]
 
    elif(npt==5):
    
        pontos[0]  = -0.906179845938664
        pontos[1]  = -0.538469310105683
        pontos[2]  =  0.0;
        pontos[3]  = -pontos[1];
        pontos[4]  = -pontos[0];
        pesos[0]   = 0.236926885056189
        pesos[1]   = 0.478628670499366
        pesos[2]   = 0.568888888888889
        pesos[3]   = pesos[1]
        pesos[4]   = pesos[0]
    
    return pontos, pesos
#==============================================================================

#==============================================================================
def int2d(x1,x2,y1,y2,npix,npiy,f):

    pontosx, pesosx = pontopeso(npix)
    pontosy, pesosy = pontopeso(npiy)

    S     = 0
    int2d = 0
    
    for k1 in range(0,npix):
        
        for k2 in range(0,npiy):
      
            Lx = (x2-x1)*0.5*pontosx[k1]+(x2+x1)*0.5
            Ly = (y2-y1)*0.5*pontosy[k2]+(y2+y1)*0.5
            S  = S + pesosx[k1]*pesosy[k2]*f(Lx,Ly);

    int2d = (x2-x1)*0.5*(y2-y1)*0.5*S

    return int2d
#==============================================================================

#==============================================================================
def int2dv(x1,x2,y1,y2,npix,npiy,f):

    pontosx, pesosx = pontopeso(npix)
    pontosy, pesosy = pontopeso(npiy)

    S     = 0
    int2d = 0
    cte   = (x2-x1)*0.5*(y2-y1)*0.5
    
    Lx = (x2-x1)*0.5*pontosx[:]+(x2+x1)*0.5
    Ly = (y2-y1)*0.5*pontosy[:]+(y2+y1)*0.5
    
    LX, LY = np.meshgrid(Lx,Ly)
    PX,PY  = np.meshgrid(pesosx,pesosy)
    F1     = f(LX,LY)
    int2d  = cte*np.sum(PX*PY*F1)
        
    return int2d
#==============================================================================

#==============================================================================
def f1(x,y):
    
    a = x**1
    
    return a
#==============================================================================

#==============================================================================
x1 = 0
x2 = 2

y1 = 0
y2 = 1

npix = 5
npiy = 5

res1 = int2d(x1,x2,y1,y2,npix,npiy,f1)
res2 = int2dv(x1,x2,y1,y2,npix,npiy,f1)
print(res1,res2)
#==============================================================================