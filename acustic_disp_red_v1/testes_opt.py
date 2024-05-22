#==============================================================================
# -*- encoding: utf-8 -*-
#==============================================================================

#==============================================================================
# Módulos Importados do Python / Devito / Examples
#==============================================================================

#==============================================================================
# Pyhton Modules and Imports
#==============================================================================
import numpy as np
import math  as mt
#==============================================================================

#==============================================================================
class teste1:
# Configuração Non-Absorbing para o Teste 1
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 101     # Número de Pontos Direção X
    npty   = 101     # Número de Pontos Direção Y
    x0     =    0.   # Ponto Inicial Direção X
    x1     = 3000.   # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     = 3000.   # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 900.     # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
    f0     = 0.03     # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 1500.    # Posição da Fonte em X
    yposf  = 1500.    # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  = 1500                      # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.05])      # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]          # Condição CFL
    jumpv      = np.array([300])       # Vetor com Diferentes Valores de Jump
    jump       = int(jumpv[dttype])               # Intervalo de Dados a Serem Salvos
    tou        = 2                                # Time Order Displacement 
    sou        = 20                               # Space Order Displacement
    mvalue     = int(sou/2)                       # First Parameter for Stencil
    nvalue     = 3                                # Second Parameter for Stencils
    shape      = 'csq'                            # Stencil Geometry
    method     = 'dispte'                         # FD Method
    btype      = 3                                # Boundary Type    
    ftype      = 0                                # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#==============================================================================

#==============================================================================
class teste1_ref1:
# Configuração Non-Absorbing para o Teste 1
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 501     # Número de Pontos Direção X
    npty   = 501     # Número de Pontos Direção Y
    x0     =    0.   # Ponto Inicial Direção X
    x1     = 3000.   # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     = 3000.   # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 900.     # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
    f0     = 0.03     # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 1500.    # Posição da Fonte em X
    yposf  = 1500.    # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  = 1500                      # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.05])      # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]                    # Condição CFL
    jumpv      = np.array([1500])                # Vetor Intervalo de Dados a Serem Salvos 
    jump       = int(jumpv[dttype])              # Intervalo de Dados a Serem Salvos
    tou        = 2                               # Time Order Displacement 
    sou        = 40                              # Space Order Displacement
    mvalue     = int(sou/2)                      # First Parameter for Stencil
    nvalue     = 1                               # Second Parameter for Stencils
    shape      = 'cl'                            # Stencil Geometry
    method     = 'spatte'                        # FD Method
    btype      = 3                               # Boundary Type    
    ftype      = 0                               # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#==============================================================================

#==============================================================================
class teste2:
# Configuração Non-Absorbing para o Teste 1
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 101     # Número de Pontos Direção X
    npty   = 101     # Número de Pontos Direção Y
    x0     =    0.   # Ponto Inicial Direção X
    x1     = 2000.   # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     = 2000.   # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 540.     # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
    f0     = 0.025    # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 1000.    # Posição da Fonte em X
    yposf  = 800.     # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  = 1000.                     # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.075])                # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]                     # Condição CFL
    jumpv      = np.array([360])                  # Vetor com Diferentes Valores de Jump
    jump       = int(jumpv[dttype])               # Intervalo de Dados a Serem Salvos
    tou        = 2                                # Time Order Displacement 
    sou        = 2                                # Space Order Displacement
    mvalue     = int(sou/2)                       # First Parameter for Stencil
    nvalue     = 1                                # Second Parameter for Stencils
    shape      = 'cl'                             # Stencil Geometry
    method     = 'spatte'                         # FD Method
    btype      = 1                                # Boundary Type    
    ftype      = 0                                # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#==============================================================================

#==============================================================================
class teste2_ref1:
# Configuração Non-Absorbing para o Teste 1
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 501     # Número de Pontos Direção X
    npty   = 501     # Número de Pontos Direção Y
    x0     =    0.   # Ponto Inicial Direção X
    x1     = 2000.   # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     = 2000.   # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
    f0     = 0.025    # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 1000.    # Posição da Fonte em X
    yposf  = 800.     # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  = 1000.                     # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.075])     # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]                    # Condição CFL
    jumpv      = np.array([1800])                # Vetor Intervalo de Dados a Serem Salvos 
    jump       = int(jumpv[dttype])              # Intervalo de Dados a Serem Salvos
    tou        = 2                               # Time Order Displacement 
    sou        = 40                              # Space Order Displacement
    mvalue     = int(sou/2)                      # First Parameter for Stencil
    nvalue     = 1                               # Second Parameter for Stencils
    shape      = 'cl'                            # Stencil Geometry
    method     = 'spatte'                        # FD Method
    btype      = 1                               # Boundary Type    
    ftype      = 0                               # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#==============================================================================

#==============================================================================
class teste3:
# Configuração Non-Absorbing para o Teste 2
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 301     # Número de Pontos Direção X
    npty   = 101     # Número de Pontos Direção Y
    x0     =    0.   # Ponto Inicial Direção X
    x1     = 12000.  # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     =  4000.  # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 1800.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms)      
    f0     = 0.019    # Frequência da Fonte (1 Hz = 0.001 kHz)
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 6000.    # Posição da Fonte em X
    yposf  =  100.    # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  =  20.                      # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.0625])    # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]                     # Condição CFL
    jumpv      = np.array([300])                  # Vetor com Diferentes Valores de Jump
    jump       = int(jumpv[dttype])               # Intervalo de Dados a Serem Salvos
    tou        = 2                                # Time Order Displacement 
    sou        = 2                                # Space Order Displacement
    mvalue     = int(sou/2)                       # First Parameter for Stencil
    nvalue     = 1                                # Second Parameter for Stencils
    shape      = 'cl'                             # Stencil Geometry
    method     = 'spatte'                         # FD Method
    btype      = 2                                # Boundary Type
    ftype      = 0                                # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#==============================================================================

#==============================================================================
class teste3_ref1:
# Configuração Non-Absorbing para o Teste 2
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 1501    # Número de Pontos Direção X
    npty   = 501     # Número de Pontos Direção Y
    x0     =    0.   # Ponto Inicial Direção X
    x1     = 12000.  # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     =  4000.  # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms)      
    f0     = 0.034    # Frequência da Fonte (1 Hz = 0.001 kHz)
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 6000.    # Posição da Fonte em X
    yposf  =  100.    # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  =  20.                      # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.0625])    # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]                    # Condição CFL
    jumpv      = np.array([1500])                # Vetor Intervalo de Dados a Serem Salvos 
    jump       = int(jumpv[dttype])              # Intervalo de Dados a Serem Salvos
    tou        = 2                               # Time Order Displacement 
    sou        = 40                              # Space Order Displacement
    mvalue     = int(sou/2)                      # First Parameter for Stencil
    nvalue     = 1                               # Second Parameter for Stencils
    shape      = 'cl'                            # Stencil Geometry
    method     = 'spatte'                        # FD Method
    btype      = 2                               # Boundary Type
    ftype      = 0                               # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#==============================================================================

#==============================================================================
class teste4:
# Configuração Non-Absorbing para o Teste 5
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 226     # Número de Pontos Direção X - Original 401
    npty   = 81      # Número de Pontos Direção Y
    x0     = 4000.   # Ponto Inicial Direção X
    x1     = 13000.  # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     =  3200.  # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 1800.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) # tF = 3000
    f0     = 0.015    # Frequência da Fonte (1 Hz = 0.001 kHz)
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 8500     # Posição da Fonte em X
    yposf  =   50.    # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  =  50.                      # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.0625])    # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]                    # Condição CFL
    jumpv      = np.array([300])                 # Vetor com Diferentes Valores de Jump
    jump       = int(jumpv[dttype])              # Intervalo de Dados a Serem Salvos
    tou        = 2                               # Time Order Displacement 
    sou        = 2                               # Space Order Displacement
    mvalue     = int(sou/2)                      # First Parameter for Stencil
    nvalue     = 1                               # Second Parameter for Stencils
    shape      = 'cl'                            # Stencil Geometry
    method     = 'spatte'                        # FD Method
    btype      = 2                               # Boundary Type
    ftype      = 0                               # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                   
#==============================================================================

#==============================================================================
class teste4_ref1:
# Configuração Non-Absorbing para o Teste 5
#==============================================================================

#==============================================================================
#Parâmetros de Malha e Tempo
#==============================================================================
    dttype = 0       # Choose dt Value
    nptx   = 1126    # Número de Pontos Direção X - Original 2001
    npty   = 401     # Número de Pontos Direção Y
    x0     = 4000.   # Ponto Inicial Direção X
    x1     = 13000.  # Ponto Final Direção X 
    compx  = x1-x0   # Comprimento do Domínio em X
    y0     =    0.   # Ponto Inicial Direção Y
    y1     =  3200.  # Ponto Final Direção Y
    compy  = y1-y0   # Comprimento do Domínio em Y
    hx     = (x1-x0)/(nptx-1)           # Delta x
    hy     = (y1-y0)/(npty-1)           # Delta y    
    X0     = np.linspace(x0,x1,nptx)    # Malha Direção X
    Y0     = np.linspace(y0,y1,npty)    # Malha Direção Y
    X0grid,Y0grid = np.meshgrid(X0,Y0)  # Grid Auxiliar X0 x Y0 
    t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
    tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms)      
    f0     = 0.015    # Frequência da Fonte (1 Hz = 0.001 kHz)
    nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
    xposf  = 8500.    # Posição da Fonte em X
    yposf  =   50.    # Posição da Fonte em Y
    nrec   = nptx     # Número de Receivers
    nxpos  = np.linspace(x0,x1,nrec)   # Posição dos Receivers em X
    nypos  =  50.                      # Posição dos Receivers em Y
    datainter  = 1                     # Interpolação de Dados de Velocidade
    dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
    CFLv       = np.array([0.0625])              # Vetor com Diferentes Condições CFL
    CFL        = CFLv[dttype]                    # Condição CFL
    jumpv      = np.array([1500])                # Vetor com Diferentes Valores de Jump
    jump       = int(jumpv[dttype])              # Intervalo de Dados a Serem Salvos
    tou        = 2                               # Time Order Displacement 
    sou        = 40                              # Space Order Displacement
    mvalue     = int(sou/2)                      # First Parameter for Stencil
    nvalue     = 1                               # Second Parameter for Stencils
    shape      = 'cl'                            # Stencil Geometry
    method     = 'spatte'                        # FD Method
    btype      = 2                               # Boundary Type
    ftype      = 0                               # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                
#==============================================================================