import numpy as np
from math import gamma

# =========================== Morrison ====================================
PI     = np.pi
RV     = 461.5 #GAS CONSTANT FOR WATER VAPOR
CP     = 1005 #SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
R      = 287.15 #GAS CONSTANT DRY AIR (USE FOR PARTIAL PRESSURES) sometimes called R_d in code
RHOSU  = 85000/(287.15*273.15) #STANDARD AIR DENSITY AT 850 MB
RHOI_M   = 500
rhoi_I   = 920
RHOSN  = 100
RHOG   = 400
c      = RHOI_M * PI / 6 #factor in the mass dimensional relation
EP_2   = R/RV  #constant for specific humidity calculation (dimensionless) from wrf constants file

#cut off diameter for small ice
DCS     = 125e-6 

#FALLSPEED PARAMETERS (V=AD^B)
BS       = 0.41
AG       = 19.3
BG       = 0.37
F1S, F2S = 0.86, 0.28 #VENTILATION PARAMETER FOR SNOW

#Size distribution parameters
DI, DG, DS = 3,3,3
CI      = RHOI_M*np.pi/6
CS      = RHOSN*np.pi/6
CG      = RHOG*PI/6

CONS1   = gamma(1+DS)*CS
CONS2   = 6*400*np.pi/6
CONS10  = gamma(5/2 + BS/2)
CONS11  = gamma(5/2+BG/2)
CONS12  = gamma(1+DI)*CI
CONS21  = 4/(DCS*RHOI_M)
CONS22  = PI*RHOI_M*DCS**3/6
CONS35  = 2.5+BS/2
CONS36  = 5/2+BG/2

#LAM limits 
LAMMAXI = 1./1e-6
LAMMINI = 1./(2.*DCS+100e-6)
LAMMAXS = 1/10e-6
LAMMINS = 1/2000e-6
LAMMAXG = 1/20e-6
LAMMING = 1/2000e-6

#P_c in gamma function = alpha-1 (alpha =1)
P_c     = 0
P_cs    = 0
P_cg    = 0



#ISHMAEL


QSMALL       = 1e-12   # Smallest ice mass (kg kg^-1)
rho_i        = 920    #Bulk ice density (kg m^-3)
NU           = 4     #Ice distribution parameter
nu           = NU
fourthirdspi = 4*np.pi/3
gammnu       = 6
i_gammnu     = 1/6
G_HOME       = 9.8 #(gravity)
T0           = 273.15
QSMALL       = 1e-12   # Smallest ice mass (kg kg^-1)
QNSMALL      = 1.25e-7 # Smallest ice number (# kg^-1)
QASMALL      = 1e-24   # Smallest ice volume (m^-3 kg^-1)
ao           = 0.1e-6  # alphstr=ao^(1-deltastr)  
RCP          = 0.285856574
CP           = 1004      #Heat capacity Air (J kg^-1 K^-1)
i_cp         = 1/CP
RD           = 287.15     # Dry air gas constant (J kg^-1 K^-1)
