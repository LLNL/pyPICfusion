
import numpy as np
import pylab as pl
pl.ion()

#import scipy.io
#from scipy.interpolate import interp1d

#font = {'family' : 'normal',
#        'size'   : 8}
#pl.rc('font', **font)

#pl.rcParams['axes.formatter.limits']= (-3,4)


def xsec(E,react='DD_nHe',units='mbarn'):
    """ Uses Bosch-Hale to give cross-section
    # E in keV"""
    
    if react=='DD_nHe' :
        BG = 31.3970 # sqrt(keV)
    
        A1 =  5.3701e4
        A2 =  3.3027e2
        A3 = -1.2706e-1
        A4 =  2.9327e-5
        A5 = -2.5151e-9
        B1 =  0.0
        B2 =  0.0
        B3 =  0.0
        B4 =  0.0
        
    elif react=='DD_pt' :

        BG = 31.3970 # sqrt(keV)
        
        A1 =  5.5576e4
        A2 =  2.1054e2
        A3 = -3.2638e-2
        A4 =  1.4987e-6
        A5 =  1.8181e-10
        B1 =  0.0
        B2 =  0.0
        B3 =  0.0
        B4 =  0.0

    elif react=='DT' :

        BG =  34.3827 # sqrt(keV)
        
        A1 =  6.927e4
        A2 =  7.454e8
        A3 =  2.050e6
        A4 =  5.2002e4
        A5 =  0.0
        B1 =  6.38e1
        B2 = -9.95e-1
        B3 =  6.981e-5
        B4 =  1.728e-4       

    if type(E) is np.ndarray :
        index = E > 0.01 
        Ei = E[index]
        SE_top = A1 + Ei*(A2 + Ei*(A3 + Ei*(A4 + Ei*A5))) 
        SE_bot = 1. + Ei*(B1 + Ei*(B2 + Ei*(B3 + Ei*B4)))
        
        SE = SE_top / SE_bot
        
        sigma = np.zeros_like(E)
        sigma[index] = SE / ( Ei * np.exp(BG / np.sqrt(Ei) ) )

        
    else :
        
        
        if E > 0.01 :  # > 10 eV  
    
            SE_top = A1 + E*(A2 + E*(A3 + E*(A4 + E*A5))) 
            SE_bot = 1. + E*(B1 + E*(B2 + E*(B3 + E*B4)))
            
            SE = SE_top / SE_bot
            
            sigma = SE / ( E * np.exp(BG / np.sqrt(E) ) )
        else :
            sigma = 0.0
        
    if units == 'mbarn' :
        sigma = sigma * 1.0
    elif units == 'barn' :
        sigma = sigma * 1.0e-3        
    elif units == 'cm2' :
        sigma = sigma * 1.0e-27 
    else : 
        print 'unit ' + units + ' is unknown, try mbarn, barn,or cm2'   
        1/0
        
        
    return sigma
    
    
def sigmav(T,react='DD_nHe',units='cm3/ns'):
    # T in keV
    
    if react=='DD_nHe' :
        BG = 31.3970 # sqrt(keV)
        
        C1 = 5.43360e-12
        C2 = 5.85778e-3
        C3 = 7.68222e-3
        C4 = 0.0
        C5 = -2.96400e-6
        C6 = 0.0
        C7 = 0.0
        
        mrc2 = 937814. # keV    

    elif react=='DT' :
        BG = 34.3827 # sqrt(keV)
        
        C1 =  1.17302e-9
        C2 =  1.51361e-2
        C3 =  7.51886e-2
        C4 =  4.60643e-3
        C5 =  1.35000e-2
        C6 = -1.06750e-4
        C7 =  1.36600e-5
        
        mrc2 = 1124656. # keV            

    Theta_top =       T * (C2 + T * (C4 + T*C6 ))
    Theta_bot = 1. +  T * (C3 + T * (C5 + T*C7 ))
    Theta = T / (1. - (Theta_top/Theta_bot) )
    xchi = ( BG**2 / (4*Theta) ) ** (1./3.)
    sigma_v = C1 * Theta * np.sqrt(xchi / (mrc2*T**3)) * np.exp(-3*xchi)
    
    if units == 'cm3/ns' :
        sigma_v *= 1e-9
    elif units == 'cm3/s' :
        pass
    else : 
        print 'unit ' + units + ' is unknown, try cm3/s, or cm3/ns'   
        1/0    

        
    return sigma_v # in cm3/ns

    
def sv_multiplier(n1,n2,V,react='DD_nHe'):
    
    if react in ['DD_nHe','DD_pt'] :
        delta12 = 1.0
    elif react in ['DT'] :
        delta12 = 0.0
        
    mult = ( (n1*n2)/(1+delta12) ) *  V
    return mult # in cm3  

  
def reactrate(n1,n2,T,V,react='DD_nHe'):

    mult = sv_multiplier(n1,n2,V)
    reactrate = mult * sigmav(T,react=react)
    
    return reactrate # in n/ns
    
    
