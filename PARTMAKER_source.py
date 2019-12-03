# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 10:46:54 2016
@author: higginson2
"""
from __future__ import division
import numpy as np
import pylab as pl
#from scipy.interpolate import interp1d
import datetime
import time
#
import os
import h5py
#import itertools as itt
#from multiprocessing import Pool, freeze_support
#from multiprocessing import Process
#import multiprocessing
#from mpi4py import MPI
import _BoschHale_library as BH


# font = {'family' : 'sans-serif',
#         'size'   : 12}
# 
# pl.rc('font', **font)
pl.rcParams['axes.formatter.limits']= (-3,3)
#pl.ion()

c_cms = 2.9979e10

#==============================================================================
# Making species etc 
#==============================================================================

#my_cmap = pl.get_cmap('gnuplot2'); 
class make_species(object):    
    
    def __init__(self,density,vx,vy,vz,nparts,Ti=None):
        self.den = density # 1/cm3
        self.velocity_x = vx # cm/s
        self.velocity_y = vy # cm/s
        self.velocity_z = vz # cm/s
        self.nparts = int(nparts)
        self.Ti = Ti # keV 
        #if  self.nparts%2 :
        #    print 'particle number must be even!'
        #    1/0

    def make_title(self) :
        if (self.velocity_x**2 < 1.) and (self.velocity_y**2 < 1.) and (self.velocity_z**2 < 1.) :
            vel_string = 'v:0'
        else : 
            vel_string = 'vx:%.1e, vy:%.1e, vz:%.1e cm/s' % (self.velocity_x, self.velocity_y, self.velocity_z) 
            
        string = 'n:%.0e/cm3; T:%g keV, N:%.0f, ' % (self.den, self.Ti, self.nparts) 
        string += vel_string
        return string

    def expand(self,iontype,Volume) :
        # print 'self.den',self.den
        # print 'self.nparts',self.nparts
        # print 'Volume',Volume
        self.weight =  self.den * Volume / float(self.nparts)
        self.W = np.ones(self.nparts) * self.weight
        if self.Ti is not None :

            if iontype == 'D' :
                mc2 = 1875.613e3 #keV
            elif iontype == 'T' :
                mc2 = 2809.432e3 #keV
            Temp = self.Ti # keV
            c_cms = 2.9979e10 # cm/s
            
            sigma_vth = c_cms  * np.sqrt(Temp/mc2)
            
            self.vx = np.random.normal(0.0, sigma_vth, self.nparts) + self.velocity_x
            self.vy = np.random.normal(0.0, sigma_vth, self.nparts) + self.velocity_y
            self.vz = np.random.normal(0.0, sigma_vth, self.nparts) + self.velocity_z
        else :
            self.vx = np.ones(self.nparts) * self.velocity_x
            self.vy = np.ones(self.nparts) * self.velocity_y
            self.vz = np.ones(self.nparts) * self.velocity_z



class combine_particles(object):    
    def __init__(self,S1,S2) :
         self.W = np.concatenate([S1.W,S2.W])
         self.vx = np.concatenate([S1.vx,S2.vx])
         self.vy = np.concatenate([S1.vy,S2.vy])
         self.vz = np.concatenate([S1.vz,S2.vz])

def match_particle_number(S1,S2) :
    nS1 = len(S1.W)
    nS2 = len(S2.W)
    if nS1 == nS2 :
        print 'S1 = S2'
    
    elif nS1 > nS2 :
        print 'S1 is larger'
        if nS1%nS2 > 0 :
            error_exit('ppc1 must be evenly divisible by ppc2')    
        mult = int(nS1/nS2)
        S2.W = np.repeat(S2.W,mult) / float(mult)
        S2.vx = np.repeat(S2.vx,mult)
        S2.vy = np.repeat(S2.vy,mult)
        S2.vz = np.repeat(S2.vz,mult)

    elif nS2 > nS1 :
        print 'S2 is larger'
        print 'nS1',nS1
        print 'nS2',nS2
        if nS2%nS1 > 0 :
            error_exit('ppc2 must be evenly divisible by ppc1')    
        mult = int(nS2/nS1)
        S1.W = np.repeat(S1.W,mult) / float(mult)
        S1.vx = np.repeat(S1.vx,mult)
        S1.vy = np.repeat(S1.vy,mult)
        S1.vz = np.repeat(S1.vz,mult)

    return S1, S2

    


#==============================================================================
# Write and read HDF Reactant Particle Files 
#==============================================================================



def save_particle_reactants(hdfname,S) :

    #print 'saving dhf5 part file'
    h5f = h5py.File(hdfname, 'w')
    #compression = 'lzf' # This compresses, but not readable by h5dump
    compression = None # larger file, but is readable by h5dump, good for testing

    h5f.create_dataset('W', data=S.W,compression=compression,shuffle=True)
    h5f.create_dataset('vx', data=S.vx,compression=compression,shuffle=True)
    h5f.create_dataset('vy', data=S.vy,compression=compression,shuffle=True)
    h5f.create_dataset('vz', data=S.vz,compression=compression,shuffle=True)
    h5f.close()
    #print 'FINISHED:: saving dhf5 part file'

def load_particle_reactants(hdfname) :
    class holder(object):
       pass
    S = holder()

    f = h5py.File(hdfname,'r')
    S.W  = f['W'].value
    S.vx = f['vx'].value
    S.vy = f['vy'].value
    S.vz = f['vz'].value
    f.close()

    return S


