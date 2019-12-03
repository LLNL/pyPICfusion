# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 16:18:04 2016
@author: higginson2
"""
from __future__ import division
import numpy as np
import pylab as pl
#import FUSION_source as FS
#import PARTMAKER_source as PS
import MAKEPLOTS_source as MP
#import _BoschHale_library as BH
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
    

    #productfilename = 'Test1_product3.hd5'
    productfilename = 'K/Anna_000070_product3.hd5'
    p3 = MP.load_particle_products(productfilename)
    MP.calculations(p3)
    fig1,ax1 = MP.plot_1D_Energy_hist(p3)
    pl.sca(ax1)
    pl.yscale('log')

    fig2,ax2 = MP.plot_2D_Energy_vs_acosvzvt(p3)


    pl.show()

    sys.exit()





