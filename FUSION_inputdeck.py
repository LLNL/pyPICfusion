# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 16:18:04 2016
@author: higginson2
"""
from __future__ import division
import numpy as np
import pylab as pl
import FUSION_source as FS
import PARTMAKER_source as PS
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
    #partfilename = 'Test1_parts.hd5' 
    savename = 'Anna_5e3iter_000070'
    savedir = 'K'
    partfilename = 'K/ion_all.h5' 
    print partfilename

    Fmult = 1e13

    dt = 1e-9
    nsteps = 5000.
    do_plot1 = True # plot as function of step

    S1 = PS.load_particle_reactants(partfilename)
    S2 = None

    # S1 is an object with S.W, S.vx, S.vy, S.vz
    # given arrays W, vx, vy, vx this can be made via the following :
    # 
    # class holder(object):
    #    pass
    # S = holder()
    # 
    # S.W  = W
    # S.vx = vx 
    # S.vy = vy
    # S.vz = vz
    #

    do_mpi = True # This is not implemented, there is no mpi at the moment
    Volume = 1. # cm3

    compare_with = 1.0

#   Initialize Simulation
    Run = FS.run_parameters(S1,S2,dt,Volume,fusiontype='DD',comparewith=compare_with,savename=savename,do_plot1=do_plot1,Fmult=Fmult,savedir=savedir)
    Sim = FS.sim_parameters(nsteps,do_mpi=do_mpi,subtime=1.)

#   Run Simulation
    done = FS.run_simulation(Run,Sim)
    #print done, ' finished running'
    #os._exit(1)
    sys.exit()





