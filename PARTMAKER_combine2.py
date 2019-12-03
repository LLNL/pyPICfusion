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
    

    filename1 = 'ionM.h5'
    filename2 = 'ionP.h5'
    
    S1 = PS.load_particle_reactants(filename1)
    S2 = PS.load_particle_reactants(filename2)
    Sout = PS.combine_particles(S1,S2)
    
    PS.save_particle_reactants('ion_all.h5',Sout)

    sys.exit()





