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
    def __init__(self,dx_ref,dt_ref,freq_ref):
        
        self.dx_ref   = dx_ref
        self.dt_ref   = dt_ref
        self.freq_ref = freq_ref
        self.nptx     = dx_ref*200+1 # Número de Pontos Direção X
        self.npty     = dx_ref*200+1 # Número de Pontos Direção Y
        self.x0       =    0.   # Ponto Inicial Direção X
        self.x1       = 6000.   # Ponto Final Direção X 
        self.compx    = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0       =    0.   # Ponto Inicial Direção Y
        self.y1       = 6000.   # Ponto Final Direção Y
        self.compy    = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx       = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy       = (self.y1-self.y0)/(self.npty-1)           # Delta y

        self.deltax = 0.1*self.compx  # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy  # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X
        
        self.X0       = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0       = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)    # Grid Auxiliar X0 x Y0 
        self.t0         = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn         = 1500.     # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
        self.f0         = np.around(0.030 - 0.010*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte     = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf      = 3000.    # Posição da Fonte em X
        self.yposf      = 3000.    # Posição da Fonte em Y
        self.nrec       = self.nptx     # Número de Receivers
        self.nxpos      = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos      = 3000                      # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.05*self.dx_ref*self.dt_ref     # Condição CFL
        self.jump       = int(300/self.dt_ref)             # Intervalo de Dados a Serem Salvos
        self.tou        = 2                                # Time Order Displacement 
        self.sou        = 2                                # Space Order Displacement
        self.mvalue     = int(self.sou/2)                  # First Parameter for Stencil
        self.nvalue     = 1                                # Second Parameter for Stencils
        self.shape      = 'cl'                             # Stencil Geometry
        self.method     = 'spatte'                         # FD Method
        self.btype      = 4                                # Boundary Type    
        self.ftype      = 0                                # Source type
#==============================================================================

#==============================================================================
class teste1_ref1:
# Configuração Non-Absorbing para o Teste 1
#==============================================================================

#==============================================================================
    def __init__(self,freq_ref,factor_ref):
        
        self.freq_ref   = freq_ref
        self.factor_ref = factor_ref
        self.nptx   = factor_ref*200+1 # Número de Pontos Direção X
        self.npty   = factor_ref*200+1 # Número de Pontos Direção Y
        self.x0     =    0.   # Ponto Inicial Direção X
        self.x1     = 6000.   # Ponto Final Direção X 
        self.compx  = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0     =    0.   # Ponto Inicial Direção Y
        self.y1     = 6000.   # Ponto Final Direção Y
        self.compy  = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx     = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy     = (self.y1-self.y0)/(self.npty-1)           # Delta y

        self.deltax = 0.1*self.compx   # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy   # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X

        self.X0     = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0     = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)  # Grid Auxiliar X0 x Y0 
        self.t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn     = 1500.     # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
        self.f0     = np.around(0.030 - 0.010*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf  = 3000.    # Posição da Fonte em X
        self.yposf  = 3000.    # Posição da Fonte em Y
        self.nrec   = self.nptx     # Número de Receivers
        self.nxpos  = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos  = 3000                      # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.05                  # Condição CFL
        self.jump       = int(300*factor_ref)             # Intervalo de Dados a Serem Salvos
        self.tou        = 2                               # Time Order Displacement 
        self.sou        = 20                              # Space Order Displacement
        self.mvalue     = int(self.sou/2)                 # First Parameter for Stencil
        self.nvalue     = 1                               # Second Parameter for Stencils
        self.shape      = 'cl'                            # Stencil Geometry
        self.method     = 'spatte'                        # FD Method
        self.btype      = 4                               # Boundary Type    
        self.ftype      = 0                               # Source type
#==============================================================================

#==============================================================================
class teste2:
# Configuração Non-Absorbing para o Teste 1
#==============================================================================

#==============================================================================
    def __init__(self,dx_ref,dt_ref,freq_ref):

        self.dx_ref   = dx_ref
        self.dt_ref   = dt_ref
        self.freq_ref = freq_ref
        self.nptx   = dx_ref*200+1   # Número de Pontos Direção X
        self.npty   = dx_ref*200+1   # Número de Pontos Direção Y
        self.x0     =    0.   # Ponto Inicial Direção X
        self.x1     = 4000.   # Ponto Final Direção X 
        self.compx  = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0     =    0.   # Ponto Inicial Direção Y
        self.y1     = 4000.   # Ponto Final Direção Y
        self.compy  = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx     = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy     = (self.y1-self.y0)/(self.npty-1)           # Delta y

        self.deltax = 0.1*self.compx  # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy  # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X

        self.X0     = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0     = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)  # Grid Auxiliar X0 x Y0 
        self.t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
        self.f0     = np.around(0.025 - 0.010*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf  = 2000.    # Posição da Fonte em X
        self.yposf  = 1600.     # Posição da Fonte em Y
        self.nrec   = self.nptx     # Número de Receivers
        self.nxpos  = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos  = 2000.                     # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.075*dx_ref*dt_ref   # Condição CFL
        self.jump       = int(360/dt_ref)       # Intervalo de Dados a Serem Salvos
        self.tou        = 2                                # Time Order Displacement 
        self.sou        = 2                                # Space Order Displacement
        self.mvalue     = int(self.sou/2)                  # First Parameter for Stencil
        self.nvalue     = 1                                # Second Parameter for Stencils
        self.shape      = 'cl'                             # Stencil Geometry
        self.method     = 'spatte'                         # FD Method
        self.btype      = 4                                # Boundary Type    
        self.ftype      = 0                                # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                  
#==============================================================================

#==============================================================================
class teste2_ref1:
# Configuração Non-Absorbing para o Teste 1
#==============================================================================

#==============================================================================
    def __init__(self,freq_ref,factor_ref):
        
        self.freq_ref   = freq_ref
        self.factor_ref = factor_ref
        self.nptx   = factor_ref*200+1   # Número de Pontos Direção X
        self.npty   = factor_ref*200+1   # Número de Pontos Direção Y
        self.x0     =    0.   # Ponto Inicial Direção X
        self.x1     = 4000.   # Ponto Final Direção X 
        self.compx  = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0     =    0.   # Ponto Inicial Direção Y
        self.y1     = 4000.   # Ponto Final Direção Y
        self.compy  = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx     = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy     = (self.y1-self.y0)/(self.npty-1)           # Delta y

        self.deltax = 0.1*self.compx  # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy  # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X

        self.X0     = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0     = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)  # Grid Auxiliar X0 x Y0 
        self.t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) #min de 0.2s até 2.0s      
        self.f0     = np.around(0.025 - 0.010*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf  = 2000.    # Posição da Fonte em X
        self.yposf  = 1600.    # Posição da Fonte em Y
        self.nrec   = self.nptx     # Número de Receivers
        self.nxpos  = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos  = 2000.                     # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.075                           # Condição CFL
        self.jump       = int(360*factor_ref)             # Intervalo de Dados a Serem Salvos
        self.tou        = 2                               # Time Order Displacement 
        self.sou        = 20                              # Space Order Displacement
        self.mvalue     = int(self.sou/2)                 # First Parameter for Stencil
        self.nvalue     = 1                               # Second Parameter for Stencils
        self.shape      = 'cl'                            # Stencil Geometry
        self.method     = 'spatte'                        # FD Method
        self.btype      = 4                               # Boundary Type    
        self.ftype      = 0                               # Source type  
#==============================================================================

#==============================================================================
class teste3:
# Configuração Non-Absorbing para o Teste 3
#==============================================================================

#==============================================================================
    def __init__(self,dx_ref,dt_ref,freq_ref):

        self.dx_ref   = dx_ref
        self.dt_ref   = dt_ref
        self.freq_ref = freq_ref
        self.nptx   = dx_ref*300+1   # Número de Pontos Direção X
        self.npty   = dx_ref*100+1   # Número de Pontos Direção Y
        self.x0     =    0.   # Ponto Inicial Direção X
        self.x1     = 12000.  # Ponto Final Direção X 
        self.compx  = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0     =    0.   # Ponto Inicial Direção Y
        self.y1     =  4000.  # Ponto Final Direção Y
        self.compy  = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx     = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy     = (self.y1-self.y0)/(self.npty-1)           # Delta y

        self.deltax = 0.1*self.compx  # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy  # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X

        self.X0     = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0     = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)  # Grid Auxiliar X0 x Y0 
        self.t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms)      
        self.f0     = np.around(0.020 - 0.005*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf  = 6000.    # Posição da Fonte em X
        self.yposf  =  100.    # Posição da Fonte em Y
        self.nrec   = self.nptx     # Número de Receivers
        self.nxpos  = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos  =  20.                      # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.0625*dx_ref*dt_ref  # Condição CFL
        self.jump       = int(300/dt_ref)       # Intervalo de Dados a Serem Salvos
        self.tou        = 2                                # Time Order Displacement 
        self.sou        = 2                                # Space Order Displacement
        self.mvalue     = int(self.sou/2)                  # First Parameter for Stencil
        self.nvalue     = 1                                # Second Parameter for Stencils
        self.shape      = 'cl'                             # Stencil Geometry
        self.method     = 'spatte'                         # FD Method
        self.btype      = 2                                # Boundary Type
        self.ftype      = 0                                # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#==============================================================================

#==============================================================================
class teste3_ref1:
# Configuração Non-Absorbing para o Teste 2
#==============================================================================

#==============================================================================
    def __init__(self,freq_ref,factor_ref):
        
        self.freq_ref   = freq_ref
        self.factor_ref = factor_ref
        self.nptx   = factor_ref*300+1   # Número de Pontos Direção X
        self.npty   = factor_ref*100+1   # Número de Pontos Direção Y
        self.x0     =    0.   # Ponto Inicial Direção X
        self.x1     = 12000.  # Ponto Final Direção X 
        self.compx  = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0     =    0.   # Ponto Inicial Direção Y
        self.y1     =  4000.  # Ponto Final Direção Y
        self.compy  = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx     = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy     = (self.y1-self.y0)/(self.npty-1)           # Delta y

        self.deltax = 0.1*self.compx  # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy  # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X

        self.X0     = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0     = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)  # Grid Auxiliar X0 x Y0 
        self.t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms)      
        self.f0     = np.around(0.020 - 0.005*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf  = 6000.    # Posição da Fonte em X
        self.yposf  =  100.    # Posição da Fonte em Y
        self.nrec   = self.nptx     # Número de Receivers
        self.nxpos  = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos  =  20.                      # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.0625                # Condição CFL
        self.jump       = int(300*factor_ref)             # Intervalo de Dados a Serem Salvos
        self.tou        = 2                               # Time Order Displacement 
        self.sou        = 20                              # Space Order Displacement
        self.mvalue     = int(self.sou/2)                 # First Parameter for Stencil
        self.nvalue     = 1                               # Second Parameter for Stencils
        self.shape      = 'cl'                            # Stencil Geometry
        self.method     = 'spatte'                        # FD Method
        self.btype      = 2                               # Boundary Type
        self.ftype      = 0                               # Source type
#==============================================================================

#==============================================================================
class teste4:
# Configuração Non-Absorbing para o Teste 5
#==============================================================================

#==============================================================================
    def __init__(self,dx_ref,dt_ref,freq_ref):

        self.dx_ref   = dx_ref
        self.dt_ref   = dt_ref
        self.freq_ref = freq_ref
        self.nptx   = dx_ref*225+1   # Número de Pontos Direção X - Original 401
        self.npty   = dx_ref*80+1    # Número de Pontos Direção Y
        self.x0     = 4000.   # Ponto Inicial Direção X
        self.x1     = 13000.  # Ponto Final Direção X 
        self.compx  = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0     =    0.   # Ponto Inicial Direção Y
        self.y1     =  3200.  # Ponto Final Direção Y
        self.compy  = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx     = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy     = (self.y1-self.y0)/(self.npty-1)           # Delta y

        self.deltax = 0.1*self.compx  # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy  # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X

        self.X0     = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0     = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)  # Grid Auxiliar X0 x Y0 
        self.t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms) # tF = 3000
        self.f0     = np.around(0.015 - 0.005*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf  = 8500     # Posição da Fonte em X
        self.yposf  =   50.    # Posição da Fonte em Y
        self.nrec   = self.nptx     # Número de Receivers
        self.nxpos  = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos  =  50.                      # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.0625*dx_ref*dt_ref  # Condição CFL
        self.jump       = int(300/dt_ref)       # Intervalo de Dados a Serem Salvos
        self.tou        = 2                               # Time Order Displacement 
        self.sou        = 2                               # Space Order Displacement
        self.mvalue     = int(self.sou/2)                 # First Parameter for Stencil
        self.nvalue     = 1                               # Second Parameter for Stencils
        self.shape      = 'cl'                            # Stencil Geometry
        self.method     = 'spatte'                        # FD Method
        self.btype      = 2                               # Boundary Type
        self.ftype      = 0                               # Source type                                                                                                                                                                                                                                                                                                                                                                                                                                                   
#==============================================================================

#==============================================================================
class teste4_ref1:
# Configuração Non-Absorbing para o Teste 5
#==============================================================================

#==============================================================================
    def __init__(self,freq_ref,factor_ref):
        
        self.freq_ref   = freq_ref
        self.factor_ref = factor_ref
        self.nptx   = factor_ref*225+1   # Número de Pontos Direção X - Original 401
        self.npty   = factor_ref*80+1    # Número de Pontos Direção Y
        self.x0     = 4000.   # Ponto Inicial Direção X
        self.x1     = 13000.  # Ponto Final Direção X 
        self.compx  = self.x1-self.x0   # Comprimento do Domínio em X
        self.y0     =    0.   # Ponto Inicial Direção Y
        self.y1     =  3200.  # Ponto Final Direção Y
        self.compy  = self.y1-self.y0   # Comprimento do Domínio em Y
        self.hx     = (self.x1-self.x0)/(self.nptx-1)           # Delta x
        self.hy     = (self.y1-self.y0)/(self.npty-1)           # Delta y    
        
        self.deltax = 0.1*self.compx  # Comprimento do Absorving Layer em X
        self.deltay = 0.1*self.compy  # Comprimento do Absorving Layer em Y
        self.npmlx  = int(self.deltax/self.hx)  # Número de Pontos no Absroving Layer em X
        self.npmly  = int(self.deltay/self.hy)  # Número de Pontos no Absroving Layer em Y
        
        self.x0    = self.x0 - self.deltax
        self.x1    = self.x1 + self.deltax
        self.compx = self.compx + 2*self.deltax
        self.nptx  = self.nptx  + 2*self.npmlx 
        self.y0    = self.y0 - self.deltay
        self.y1    = self.y1 + self.deltay
        self.compy = self.compy + 2*self.deltay
        self.npty  = self.npty  + 2*self.npmly
        
        self.x0pml  = self.x0 + self.npmlx*self.hx # Ponto Inicial do PML em X
        self.x1pml  = self.x1 - self.npmlx*self.hx # Ponto Final do PML em X
        self.y0pml  = self.y0 + self.npmly*self.hy # Ponto Inicial do PML em X
        self.y1pml  = self.y1 - self.npmly*self.hy # Ponto Final do PML em X

        self.X0     = np.linspace(self.x0,self.x1,self.nptx)    # Malha Direção X
        self.Y0     = np.linspace(self.y0,self.y1,self.npty)    # Malha Direção Y
        self.X0grid,self.Y0grid = np.meshgrid(self.X0,self.Y0)  # Grid Auxiliar X0 x Y0 
        self.t0     = 0.       # Tempo Inicial da Simulação em Milisegundos
        self.tn     = 3000.    # Tempo Final   da Simulação em Milisegundos (1 Segundo = 1000 Ms)      
        self.f0     = np.around(0.015 - 0.005*(self.freq_ref-1),4)   # Frequência da Fonte (1 Hz = 0.001 kHz) #min de 0.02 até 0.05
        self.nfonte = 1        # Número de Fontes  (Para nfonte>1 -> Vetor em xposf e yposf)
        self.xposf  = 8500.    # Posição da Fonte em X
        self.yposf  =   50.    # Posição da Fonte em Y
        self.nrec   = self.nptx     # Número de Receivers
        self.nxpos  = np.linspace(self.x0,self.x1,self.nrec)   # Posição dos Receivers em X
        self.nypos  =  50.                      # Posição dos Receivers em Y
        self.datainter  = 1                     # Interpolação de Dados de Velocidade
        self.dataintera = 0                     # Interpolação de Dados de Velocidade Artificial
        self.CFL        = 0.0625                # Condição CFL
        self.jump       = int(300*factor_ref)             # Intervalo de Dados a Serem Salvos
        self.tou        = 2                               # Time Order Displacement 
        self.sou        = 20                              # Space Order Displacement
        self.mvalue     = int(self.sou/2)                 # First Parameter for Stencil
        self.nvalue     = 1                               # Second Parameter for Stencils
        self.shape      = 'cl'                            # Stencil Geometry
        self.method     = 'spatte'                        # FD Method
        self.btype      = 2                               # Boundary Type
        self.ftype      = 0                               # Source type
#==============================================================================