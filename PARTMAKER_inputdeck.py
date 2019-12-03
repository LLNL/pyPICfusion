# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 16:18:04 2016
@author: higginson2
"""
from __future__ import division
import numpy as np
import pylab as pl
#import FUSION_source as FS
import PARTMAKER_source as PS
#import MAKEPLOTS_source as MP
import _BoschHale_library as BH
import os
import sys

font = {'family' : 'sans-serif',
        'size'   : 12}

pl.rc('font', **font)
pl.rcParams['axes.formatter.limits']= (-3,3)
#pl.ion()
pl.ioff()


        

#==============================================================================
# Main program
#==============================================================================
if __name__ == '__main__':
    # background...
    

    #savename = 'Test1' 
    savename = 'Test2' 

    if savename == 'Test1' :

        nA = 1e16 # 1/cm3
        nB = 1e16 # 1/cm3
        ppcA = 400 # D
        ppcB = 400 # D


        vxA =  0.0 # cm/s
        vyA =  0.0 # cm/s
        vzA =  0.0 # cm/s
        TiA =  19.0 # keV

        vxB =  0.0 # cm/s
        vyB =  0.0 # cm/s
        vzB =  0.0 # cm/s
        TiB =  TiA # keV

    elif savename == 'Test2' :

        nA = 1e16 # 1/cm3
        nB = 1e06 # 1/cm3
        ppcA = 400 # D
        ppcB = 400 # D


        vxA =  0.0 # cm/s
        vyA =  0.0 # cm/s
        vzA =  2e8 # cm/s
        TiA =  0.0 # keV

        vxB =  0.0 # cm/s
        vyB =  0.0 # cm/s
        vzB =  0.0 # cm/s
        TiB =  TiA # keV

    print savename

    SA = PS.make_species(nA,vxA,vyA,vzA,ppcA,Ti=TiA) 
    SB = PS.make_species(nB,vxB,vyB,vzB,ppcB,Ti=TiB) 

    SA.expand('D',Volume=1) 
    SB.expand('D',Volume=1) 
    
    S = PS.combine_particles(SA,SB)
    #print S.__dict__
    #print S.__dict__.keys()

    PS.save_particle_reactants(savename+'_parts.hd5',S)
    OUT = PS.load_particle_reactants(savename+'_parts.hd5')

    #print S.vz
    #print OUT.vz


    sys.exit()





