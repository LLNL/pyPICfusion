# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 10:46:54 2016
@author: higginson2
"""
from __future__ import division
import numpy as np
import pylab as pl
import os
import h5py


# font = {'family' : 'sans-serif',
#         'size'   : 12}
# 
# pl.rc('font', **font)
pl.rcParams['axes.formatter.limits']= (-3,3)
#pl.ion()


c_cms = 2.9979e10

#==============================================================================
# Read HDF Reactant Product Files 
#==============================================================================


def load_particle_products(hdfname) :
    class holder(object):
       pass
    S = holder()

    f = h5py.File(hdfname,'r')
    S.W  = f['W'].value
    S.vx = f['vx'].value
    S.vy = f['vy'].value
    S.vz = f['vz'].value
    S.E = f['E'].value
    f.close()

    vtotal = np.sqrt(S.vx**2 + S.vy**2 + S.vz**2)
    S.acosvxvt = np.arccos(S.vx/vtotal)
    S.acosvyvt = np.arccos(S.vy/vtotal)
    S.acosvzvt = np.arccos(S.vz/vtotal)

    return S

def plot_1D_Energy_hist(P) :
    fig,ax =  pl.subplots(1,1)
    #print 'P3.E',P3.E
    #print 'P3.W',P3.W
    pl.hist(P.E,weights=P.W,bins=200)
    #pl.title(title,fontsize=10)    
    pl.xlim(2000.0,3000.0)
    pl.xlabel('Energy [keV]')
    pl.ylabel('PDF [number/keV]')
    
    return fig, ax


def plot_2D_Energy_vs_acosvzvt(P) :
    yvar = P.E
    xvar = np.rad2deg(P.acosvzvt)
    weights = P.W / ( 2 * np.pi * np.sin(np.deg2rad(xvar)) )

    xbins = np.linspace(0,180,100)
    ybins = 100


    Z, yedges, xedges = np.histogram2d(yvar,xvar, bins=(ybins,xbins), weights=weights )
    X, Y = np.meshgrid(xedges, yedges)
    Z = Z #/ ( 2 * np.pi * np.sin(np.deg2rad(X)) )

    cmax = 0.5e-9
    cmin = 1e-13
    do_log = True
    if do_log :
        Z = np.log10(Z)
        cmin = np.log10(cmin)
        cmax = np.log10(cmax)


    fig, ax = pl.subplots(1,1)
    im = ax.pcolormesh(X,Y,Z,vmin=cmin,vmax=cmax)
    pl.ylabel('Energy [keV]')
    pl.xlabel('Theta = arccos(vz/vtotal) [degree]')
    pl.colorbar(im)
    return fig, ax



def calculations(P) :

    Emin = np.min(P.E)
    Emax = np.max(P.E)
    E0 = np.sum(P.E*P.W) / np.sum(P.W)
    EFWHM = np.sqrt( 8. * np.log(2) * np.sum(P.W*(P.E-E0)**2) / np.sum(P.W) )
    print 'P.E  keV',  (P.E)
    print 'Emin  keV',  (Emin)
    print 'Emax  keV',  (Emax)
    print 'E0    keV',  (E0)
    print 'EFWHM keV',  (EFWHM)

