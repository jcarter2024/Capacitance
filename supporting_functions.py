#Functions from capacitance 

from constants import *
import numpy as np
from math import gamma

def get_igr(temp):
    
    topython = -1 #to_python is a conversion as fortran indexes the first element as 1, but python as 0
                  #therefore igrdata[1] --> igrdata[0]
    igrdata =(0.910547, 0.81807, 0.6874, 0.60127, 1.59767, 2.32423, 2.08818,  
            1.61921, 1.15865, 0.863071, 0.617586, 0.453917, 0.351975, 0.28794, 
            0.269298, 0.28794, 0.333623, 0.418883, 0.56992, 0.796458, 1.14325, 
            1.64103, 1.90138, 1.82653, 1.61921, 1.47436, 1.32463, 1.25556,     
            1.22239, 1.206, 1.11522, 1.10751, 1.10738, 1.11484, 1.12234,       
            1.12221, 1.14529, 1.16884, 1.20104, 1.22573, 1.25094, 1.27666,     
            1.31183, 1.3388, 1.35704, 1.37553, 1.38479, 1.39411, 1.40349,      
            1.41294, 1.42245, 1.43202, 1.44166, 1.45137, 1.46114, 1.47097, 
            1.48087, 1.50105, 1.50087, 1.51098 )
    
    """Inherent growth data from chen and lamb"""
    
    #Inherent growth ratio (IGR) from data See Lamb and Scott (1972) and Chen and Lamb (1994)
    if((temp-T0) >= -59 and (temp-T0) <= -1):
        dum     = (abs(int(temp-T0)) + 1) - abs(temp-T0)
        igr1    = igrdata[max((int(temp-T0)*(-1))+topython  ,  1+topython)]
        igr2    = igrdata[min(((int(temp-T0)*(-1))+1)+topython  ,  60+topython)]
        get_igr = dum*igr1 + (1-dum)*igr2
        
    elif((temp-T0) > -1 and (temp-T0) < 0):
        dum     = 1 - abs(temp-T0)
        igr1    = 1
        igr2    = igrdata[1+topython]
        get_igr = dum*igr1 + (1-dum)*igr2
        
    elif((temp-T0) > -60 and (temp-T0) < -59):
        dum     = 60 - abs(temp-T0)
        igr1    = igrdata[59+topython]
        igr2    = igrdata[60+topython]
        get_igr = dum*igr1 + (1-dum)*igr2
    
    elif((temp-T0) < -60):
        get_igr = igrdata[60+topython]
    
    else:
        get_igr = 1
        
    return get_igr

################################################################################################################
################################################################################################################
################################################################################################################


def capacitance_gamma(ani, dsdum, alphstr):
    
    NU = 4 #constant
    nu=NU
    i_gammnu = 1/(gamma(NU))
    
    #Oblate Spheroid
    if (dsdum <= 1): 
        a1 = 0.6369427      
        a2 = 0.57*a1
        b1 = 0.0
        b2 = 0.95   
        c1 = a1*alphstr**b1
        c2 = a2*alphstr**b2
        d1 = b1*(dsdum - 1.0) + 1.0
        d2 = b2*(dsdum - 1.0) + 1.0
        
    #!.. Prolate Spheroid
    elif (dsdum > 1):
        a1 = 0.5714285
        a2 = 0.75*a1
        b1 = -1.0
        b2 = -0.18 
        c1 = a1*alphstr**(b1+1.0)
        c2 = a2*alphstr**(b2+1.0)
        d1 = b1*(dsdum - 1.0) + dsdum
        d2 = b2*(dsdum - 1.0) + dsdum
      
    
    #get capacitance 
    if(dsdum <= 1):
        gammad1 = (gamma(nu+d1))
        gammad2 = (gamma(nu+d2))
        capacitance_gamma = c1*ani**d1 * gammad1*i_gammnu + c2*ani**d2 * gammad2*i_gammnu
    elif(dsdum > 1):
        gammad1 = (gamma(nu+d1))
        gammad2 = (gamma(nu+d2))
        capacitance_gamma = c1*ani**d1 * gammad1*i_gammnu  + c2*ani**d2 * gammad2*i_gammnu
    
    return capacitance_gamma


################################################################################################################
################################################################################################################
################################################################################################################

def polysvp_M(t, TYPE):
    
    alog10 = lambda x : np.log10(x)
    a0i,a1i,a2i = 6.11147274, 0.503160820, 0.188439774e-1
    a3i,a4i,a5i = 0.420895665e-3, 0.615021634e-5,0.602588177e-7
    a6i,a7i,a8i = 0.385852041e-9, 0.146898966e-11, 0.252751365e-14

    #liquid V7
    a0,a1,a2 = 6.11239921, 0.443987641, 0.142986287e-1
    a3,a4,a5 = 0.264847430e-3, 0.302950461e-5, 0.206739458e-7
    a6,a7,a8 = 0.640689451e-10,-0.952447341e-13,-0.976195544e-15
    
    #ice
    if TYPE == 1:
        if (t > 195.8):
            dt=t-273.15
            polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt))))))) 
            polysvp = polysvp*100
        else:
            polysvp = 10**(-9.09718*(273.16/t-1)-3.56654*alog10(273.16/t)+0.876793*(1-t/273.16)+alog10(6.1071))*100

    #LIQUID
    if TYPE == 0:
        if (t > 202.0):
            dt = t-273.15
            polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
            polysvp = polysvp*100.
        else:
         #uncertain below -70 C, but produces physical values (non-negative) unlike flatau
            polysvp = 10.**(-7.90298*(373.16/t-1.)+  5.02808*alog10(373.16/t)- \
                       1.3816e-7*(10**(11.344*(1-t/373.16))-1)+ 8.1328e-3*(10**(-3.49149*(373.16/t-1))-1)+ \
             alog10(1013.246))*100                      

    return polysvp


################################################################################################################
################################################################################################################
################################################################################################################


def polysvp_I(t,TYPE):
    
    #ice
    a0i,a1i,a2i = 6.11147274, 0.503160820, 0.188439774e-1
    a3i,a4i,a5i = 0.420895665e-3, 0.615021634e-5,0.602588177e-7
    a6i,a7i,a8i = 0.385852041e-9, 0.146898966e-11, 0.252751365e-14
    
    #liquid V7
    a0,a1,a2 = 6.11239921, 0.443987641, 0.142986287e-1
    a3,a4,a5 = 0.264847430e-3, 0.302950461e-5, 0.206739458e-7
    a6,a7,a8 = 0.640689451e-10,-0.952447341e-13,-0.976195544e-15
    
    #ICE
    if TYPE == 1:
        dt = max(-80, t-273.16)
        polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt)))))))
        polysvp = polysvp*100
       
    #Liquid
    if TYPE ==0:
        dt = max(-80,t-273.16)
        polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
        polysvp = polysvp*100
        
    return polysvp



################################################################################################################
################################################################################################################
################################################################################################################

def vapour_grow_ISHMAEL(dt, ani, cni, rni, igr, nidum, temp, rimesum, presdum,
       NU, alphstr, sui, sup, qvs, qvi, mu, iwci, rhodum, qidum, dv, kt, ao,
       nsch, npr, gammnu, i_gammnu, fourthirdspi, svpi, xxls, xxlv, xxlf,
       capgam, rbdum, dsdum, CST_SF='NO'):
            
    
    #lambdas to save time 
    log    = lambda x : np.log(x)
    gammln = lambda x : np.log(gamma(x))
    
    QASMALL= 1.e-19   # Smallest ice size for sublimation  (squared) --> this is different to prev
    
    #required for input
    #dt, ani, cni, rni, igr, nidum, temp, presdum, NU, alphstr, mu, iwci, rhodum, qidum
    #dv, kt, ao, nsch, npr, gammnu, i_gammnu, fourthirdspi, svpi, xxls, xxlv, xxlf, capgam, 
    #rimesum, sui, sup, qvs, qvi, dsdum
            
    #returned
    #vtbarb, vtbarbm, vtbarbz, anf, cnf, rnf, iwcf, rdout, fvdum, fhdum, dsdumout, rbdum
    
    #note gamma_arg replaced with gi and gamma_tab by gamma
    
    fvdum    = 1.0
    fhdum    = 1.0
    dsdumout = dsdum
    gi       = NU-1+dsdum
    phii     = cni/ani*gamma(gi)*i_gammnu
    
    if CST_SF=='NO':
        fs       = capgam/rni
    elif CST_SF=='YES':
        capgam2 = capacitance_gamma(rni,1,1)
        fs       = capgam2/rni
   
    alphanr  = ani/rni**(3./(2.+igr))
    
  #Determine  coefficients for the Best number for fall speeds
    if(phii < 1 ):
        bl = 1
        al = 2
        aa = PI
        ba = 2
        qe = (1-phii)*(rbdum/rhoi_I) + phii
    elif(phii > 1):
        al = 2
        bl = 1
        aa = PI*alphstr
        ba = dsdum + 1
        qe = 1
    elif(phii == 1):
        bl = 1
        al = 2
        aa = PI
        ba = 2
        qe = 1

    qe = min(qe, 1.0)
    
  #Fall speed and ventilation (Best number formulation)
    xn =  2/rbdum*(rbdum-rhodum)*G_HOME*rhodum/mu**2 *(fourthirdspi*rbdum)*alphstr*al**2/(aa*qe) * qe**(3/4)
    bx = dsdum+2+2*bl-ba
    
  #Number-average Best Number
    xm = xn*ani**bx * (gamma(nu+bx))*i_gammnu
    
  #The following fall speed coefficients are from  Mitchell and Heymsfield 2005
    f_c1 = 4.0 / (5.83 * 5.83 * np.sqrt(0.6))
    f_c2 = (5.83 * 5.83) / 4.0
    
    bm   = ((f_c1 * np.sqrt(xm)) /  (2 * (np.sqrt(1 + f_c1 * np.sqrt(xm)) - 1) *
            (np.sqrt(1 + f_c1 * np.sqrt(xm))))) - (1e-5 * xm) / (f_c2 * (np.sqrt(1 + f_c1 * np.sqrt(xm)) - 1)**2)
        
    am = ((f_c2 * (np.sqrt(1 + f_c1 * np.sqrt(xm)) - 1)**2) - (1e-5 * xm)) / (xm**bm)
    
    if(xm > 1e8):
        am = 1.0865
        bm = 0.499
    
  #Reynolds Number
    Nre = am*xm**bm          
    
  #Number-averaged ice fall speed
    vtbarb = mu/rhodum*0.5 * am*(xn)**bm*ani**(bx*bm-1) *(gamma(nu+bx*bm-1))*i_gammnu
        
  #Mass-averaged ice fall speed
    vtbarbm = mu/rhodum*0.5 * am*(xn)**bm*ani**(bx*bm-1) *(gamma(nu+bx*bm-1+2+dsdum))/(gamma(nu+2+dsdum))
    
  #.. Reflectivity-weighted fall speed (if needed)
    vtbarbz = ( mu/rhodum*0.5 * am*(xn)**bm*ani**(bx*bm-1) * np.exp(gammln(nu+bx*bm-1+4+2*dsdum))/
               np.exp(gammln(nu+4+2*dsdum)) )
        
  #.. Calculate Ventilation
    xvent = nsch**0.333333333*Nre**0.5
    ntherm = Nre**0.5*npr**0.333333333
    
    if(xvent < 1.0):
        bv1 = 1.0
        bv2 = 0.14
        gv = 2.
    else:
        bv1 = 0.86
        bv2 = 0.28
        gv = 1.

    if(ntherm < 1.4):
        bt1 = 1.0
        bt2 = 0.108
        gt = 2.0
    else:
        bt1 = 0.78
        bt2 = 0.308
        gt = 1.0

    fvdum = bv1 + bv2*xvent**gv
    fhdum = bt1 + bt2*ntherm**gt
    
    # If T < T0 do vapor growth/sublimation (otherwise assume water on surface)
    if(temp < T0):
        vtbranch = vtbarbm  #Fall speed needed to determine when branching occurs
       
        if(sup > 0):
            maxsui = 1
        elif(sui > 0 and qvi < qvs):
            maxsui = ((sui+1)*qvi)/(qvs-qvi) - qvi/(qvs-qvi)
            maxsui = min(maxsui,1)
            maxsui = max(maxsui,0)
        else:
            maxsui = 0.
    
  # Vapor growth density 
        if(igr < 1):  #--> Planar
            if(vtbranch > 0):
                if(ani > np.sqrt((dv*PI*2*cni)/(vtbranch*nu))):
                        rhodep = (rhoi_I*igr)*maxsui + rhoi_I*(1-maxsui)
                else:
                    rhodep = rhoi_I
            else:
                rhodep = rhoi_I
                
        else:# --> Columnar
            rhodep = (rhoi_I/igr)*maxsui + rhoi_I*(1-maxsui)
  
          # high limit on rhodep for vapor growth = 700 kg m^-3 Cotton et al. 2012 (QJRMS)
        rhodep = min(rhodep,700)
       
      #Vapor growth solution which includes heating from rime mass [J.P. Chen's thesis]
        alpha = (dv*fvdum*svpi*xxls)/(RV*kt*fhdum*temp)
        if (nidum > 0 and rimesum > 0):
            del1 = (xxlf*(rimesum/nidum)/(4*PI*kt*fhdum*capgam)) *((temp + alpha*((xxls/(RV*temp))-1))**(-1))
        else:
            del1 = 0
     
        del2 = sui * (((temp/alpha) + ((xxls/(RV*temp))-1))**(-1))
        Del = del1 + del2
        afn = ((dv*fvdum*polysvp_I(temp, 1))/(RV*temp)) * (sui - Del*((xxls/(RV*temp))-1))
 
  #!.. During sublimation using polynomial removal of density
        if(afn < 0): 
            rhodep = rbdum     
            videp  = rni**3
            vmin   = (10e-6)**3
            if(vmin < videp):
                betavol = np.log(rhoi_I/rbdum)*1/(np.log(vmin/videp))
                rhodep  = rbdum*(1+betavol)
            else:
                rhodep  = rbdum

        rhodep = max(rhodep,50)
        rhodep = min(rhodep,rhoi_I)
       
        gi        = NU+2+dsdum
        gammnubet = gamma(gi)

    #   !..  characteristic r-axis and a-axis after growth timestep
        rnf = (max((rni**2 + 2*afn*fs/rhodep*gammnu/gammnubet*dt),QASMALL))**(0.5) 
        anf = alphanr*rnf**(3/(2+igr))
       
    #   !.. Do not sublimation change the shape of ice from 
    #   !.. prolate to oblate or vice versa
        phif = phii*(rnf**3/rni**3)**((igr-1)/(igr+2))
        if(sui < 0 or afn < 0):
            
            if(phii > 1 and phif < 1):
                phif    = phii
                alphanr = ani/rni
                anf     = alphanr*rnf

            if(phii < 1 and phif > 1):
                phif    = phii
                alphanr = ani/rni
                anf     = alphanr*rnf

#   !.. Do not let sublimation create extreme shapes
            if(phii > 1):
                if(phif > phii):
                    phif    = phii
                    alphanr = ani/rni
                    anf     = alphanr*rnf
   
            else:
                 if(phif < phii):
                    phif    = phii
                    alphanr = ani/rni
                    anf     = alphanr*rnf
     
       #recalculate vi and final volume
        vi       = fourthirdspi*rni**3*gamma(gi)*i_gammnu 
        vf       = fourthirdspi*rnf**3*gamma(gi)*i_gammnu
        rdout    = rhodep

        rbdumtmp = rbdum*(vi/vf) + rhodep*(1.-vi/vf)
        rbdumtmp = min(rbdumtmp,rhoi_I)
        iwcf     = nidum*rbdumtmp*vf
       
    #   !.. Update delta* from vapor growth/sublimation
        if(igr != 1):
            if(anf > (1.1*ao)):
                dsdumout = (3*log(rnf)-2*log(anf)-log(ao))/ (log(anf)-log(ao))
            else:
                dsdumout=1

    #   !.. Do not let particles sublimate to sizes that are too small
        if(afn < 0 and rnf < 1e-6):
            rbdum    = rhoi_I
            rdout    = rhoi_I
            dsdumout = 1
            phif     = 1
            alphanr  = ani/rni
            anf      = alphanr*rnf
       
        if(afn < 0 and anf < 1e-6):
            rbdum    = rhoi_I
            rdout    = rhoi_I
            dsdumout = 1
            phif     = 1
            alphanr  = ani/rni
            anf      = alphanr*rnf
      
       
    #   !.. Sublimation check
        if(afn < 0 and dsdumout < 0):
            dsdumout = 1
            anf      = rnf
       
    #   !.. C-axis after vapor growth
        gi  = NU-1+dsdumout
        cnf = phif*anf*gammnu/gamma(gi)
        
    else:  #!.. T > T0
        print(temp, '>', T0)
        anf  = ani
        cnf  = cni
        rnf  = rni
        iwcf = iwci
        
    return  anf, cnf, rnf, iwcf, rdout,dsdumout, rbdum, capgam, fs, afn, rhodep


################################################################################################################
################################################################################################################
################################################################################################################



def ISHMAEL_var_check(qidum, dsdum, ani, cni, rni, rbdum, nidum, aidum, cidum):
    """ 
    Note that var check will possibly alter all of the input variables 
    Checks deltstring, recalculates ani, aidum, cidum
    Checks rbdum, recalculates ani, cni, aidum, cidum
    Checks rni, recalculates nidum, ani, cni, aidum, cidum
    Recalculates RNI
    
    
    """
    #in NU, ao, fourthirdspi, gammnu, qidum
    #out :: rni, alphstr, alphv, betam
    
    #check deltastring ----------
    if(dsdum < 0.55):
#         print("LOW DELTATRING CHECK TRIGGERED")
        dsdum = 0.55
        ani   =(cni/(ao**(1.-dsdum)))**(1./dsdum)
        aidum =ani**2*cni*nidum
        cidum =cni**2*ani*nidum
        
    elif (dsdum > 1.3):
#         print("HIGH DELTATRING CHECK TRIGGERED")
        dsdum =1.3
        cni   =ao**(1.-dsdum)*ani**dsdum
        aidum =ani**2*cni*nidum
        cidum =cni**2*ani*nidum

    alphstr = ao**(1-dsdum)
    alphv   = fourthirdspi*alphstr
    betam   = 2 + dsdum
    gi      = NU+2+dsdum

  #form ice density
    if(ani > 2e-6):
        rbdum = qidum*gamma(NU)/(nidum*alphv*ani**betam*gamma(gi))
#         print('New Density', rbdum)
    else:
        rbdum = rhoi_I
    
#     Check ice density --> Keep ice density between 50 and rhoi_I=920 kg m^-3
    if(rbdum > rhoi_I):
#         print("HIGH DENSITY CHECK TRIGGERED")
        rbdum = rhoi_I
        #as a consequence of resetting the density, recalculate the axes
        ani=((qidum*gamma(NU))/(rbdum*nidum*alphv*gamma(gi)))**(1/betam)
        cni=ao**(1-dsdum)*ani**dsdum
        aidum=ani**2*cni*nidum
        cidum=cni**2*ani*nidum
    elif(rbdum < 50):
#         print("LOW DENSITY CHECK TRIGGERED")
        rbdum=50.
        ani=((qidum*gamma(NU))/(rbdum*nidum*alphv*gamma(gi)))**(1/betam)
        cni=ao**(1.-dsdum)*ani**dsdum
        aidum=ani**2*cni*nidum
        cidum=cni**2*ani*nidum
    
    #recalculate rni
    rni= (qidum*3/(nidum*rbdum*4*PI*(gamma(gi)/gamma(NU))))**(1/3) 
          
    #   !..Small ice limit Keep rni > 2 microns
    if(rni < 2e-6):
#         print("LOW RNI CHECK TRIGGERED")
        rni   =2e-6
        nidum =3*qidum*gamma(NU)/(4*PI*rbdum*rni**3*(gamma(gi)))
        ani   =((qidum*gamma(NU))/(rbdum*nidum*alphv*gamma(gi)))**(1/betam)
        cni   =ao**(1-dsdum)*ani**dsdum
        aidum =ani**2*cni*nidum
        cidum =cni**2*ani*nidum

#   !.. Large ice limit
#   !.. This is a number weighted diameter of 8 mm and a spherical mass weighted diameter of 14 mm 
    maxsize = max(ani,cni)
    if (maxsize > 1e-3):
        if(ani > cni):
#             print("HIGH ANI/CNI CHECK TRIGGERED -->HIGH ANI")
            ani = 1e-3
            #recalculate ni cni based on new ani
            nidum=qidum*gamma(NU)/(fourthirdspi*rbdum*ao**(1-dsdum)*ani**(2+dsdum)*(gamma(gi)))
            cni=ao**(1.-dsdum)*ani**dsdum
            #recalculate ai ci based on changes
            aidum=ani**2*cni*nidum
            cidum=cni**2*ani*nidum
        else:
#             print("HIGH ANI/CNI CHECK TRIGGERED -->HIGH CNI")
            cni = 1e-3
            ani = (cni/(ao**(1-dsdum)))**(1/dsdum)
            nidum=qidum*gamma(NU)/(fourthirdspi*rbdum*ao**(1-dsdum)*ani**(2+dsdum)*(gamma(gi)))
            aidum=ani**2*cni*nidum
            cidum=cni**2*ani*nidum

       #recalculate rni based on a change of a or c 
        rni= (qidum/(nidum*rbdum*fourthirdspi*(gamma(gi)/gamma(NU))))**0.333333333333 
    
    return qidum, dsdum, ani, cni, rni, rbdum, nidum, aidum, cidum, alphstr


################################################################################################################
################################################################################################################
################################################################################################################

def PRD_Morrison(QV3D, QNI3D, QG3D, QVI, ABI, EPSS, EPSI, EPSG, LAMI):
    
    PRD, PRDS, PRDG =0, 0, 0
    
    #PRD calculation for tail-end ice (largest ice) first
    DUM = (1-np.exp(-LAMI*DCS)*(1+LAMI*DCS))
    PRD = EPSI*(QV3D-QVI)/ABI*DUM 

    #ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
    if QNI3D > 1e-14:
        PRDS = EPSS*(QV3D-QVI)/ABI + EPSI*(QV3D-QVI)/ABI*(1-DUM)
#         print("Deposition to ice     = 0")
#         print("Deposition to snow    =", "{:.2e}".format(PRDS))

    #! OTHERWISE ADD TO CLOUD ICE
    else:
        PRD += PRD+EPSI*(QV3D-QVI)/ABI*(1.-DUM)
#         print("Deposition to ice     =", "{:.2e}".format(PRD))
#         print("Deposition to snow    = 0, no snow present")

    if QG3D > 1e-14:
        #VAPOR DEPOSITION ON GRAUPEL
        PRDG = EPSG*(QV3D-QVI)/ABI
#         print("Deposition to graupel =", "{:.2e}".format(PRDG))
    else:
#         print("Deposition to graupel = 0, no graupel present")
        pass

    return PRD, PRDS, PRDG
