# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 10:46:54 2016
@author: higginson2
"""
from __future__ import division
import numpy as np
import pylab as pl
#import PARTMAKER_source as PS
import MAKEPLOTS_source as MP
import datetime
import time
#
import os
import h5py
#import itertools as itt
#from multiprocessing import Pool, freeze_support
#from multiprocessing import Process
#import multiprocessing
from mpi4py import MPI
import _BoschHale_library as BH

# font = {'family' : 'sans-serif',
#         'size'   : 12}
# 
# pl.rc('font', **font)
pl.rcParams['axes.formatter.limits']= (-3,3)
#pl.ion()

c_cms = 2.9979e10

def do_DD_fusion_old(S,Run,nsteps) :

    nsteps = int(nsteps) 

    if len(S.W)%2 > 0 :
        error_exit('len(W) must be even')    

    nhalf = int(len(S.W)/2)

    qqyield_a   = np.zeros(nsteps)
    qqyield_b   = np.zeros(nsteps)
    randyield_a = np.zeros(nsteps)
    randyield_b = np.zeros(nsteps)


    reduc_Yt = 0.5*( np.sum(S.W) )**2 - 0.5 *np.sum(S.W**2) 

    for n in range(nsteps) :
        if n%10 == 0 :
            print 'n: %.0f of %.0f, %.0f percent' % (n,nsteps,float(n)/nsteps*100.)
    
    
        index = range(len(S.W))
        np.random.shuffle(index)
        
        index1 = index[0:nhalf]
        index2 = index[nhalf::]
            
        W1  = S.W[index1]
        v1x = S.vx[index1]
        v1y = S.vy[index1]
        v1z = S.vz[index1]
        W2  = S.W[index2]
        v2x = S.vx[index2]
        v2y = S.vy[index2]
        v2z = S.vz[index2]
    
        vij = np.sqrt((v1x-v2x)**2 + (v1y-v2y)**2 + (v1z-v2z)**2 )
    
        Er = 0.5 * Run.mr * ( (vij/c_cms)**2 )
        sigma = BH.xsec(Er,react='DD_nHe',units='cm2')
        #sigma = xsec(Er,react='DD_nHe',units='cm2')

        reduc_Ys = np.sum(W1*W2) 

        Pna = reduc_Yt / reduc_Ys 
        #Pnb = float(len(S.W)) -1.
        Pnb = float(len(S.W)) 

        qqyield_a[n] = np.sum(W1*W2*vij*sigma*Run.dt/Run.Volume) * Pna
        qqyield_b[n] = np.sum(W1*W2*vij*sigma*Run.dt/Run.Volume) * Pnb

    return qqyield_a, qqyield_b        



def do_DD_fusion(S,Run,nsteps) :

    nsteps = int(nsteps) 

    if len(S.W)%2 > 0 :
        error_exit('len(W) must be even')    

    nhalf = int(len(S.W)/2)
    Nfact = float(len(S.W)) -1.
    print 'nhalf',nhalf

    W3 = []
    W4 = []
    V3 = []
    V4 = []

    qqyield   = np.zeros(nsteps)
    wwyield   = np.zeros(nsteps)

    for n in range(nsteps) :
        if n%10 == 0 :
            print 'n: %.0f of %.0f, %.0f percent' % (n,nsteps,float(n)/nsteps*100.)
    
        Wyield = 0.
    
        index = range(len(S.W))
        np.random.shuffle(index)
        
        index1 = index[0:nhalf]
        index2 = index[nhalf::]
            
        W1  = S.W[index1]
        v1x = S.vx[index1]
        v1y = S.vy[index1]
        v1z = S.vz[index1]
        W2  = S.W[index2]
        v2x = S.vx[index2]
        v2y = S.vy[index2]
        v2z = S.vz[index2]
    
        vij = np.sqrt((v1x-v2x)**2 + (v1y-v2y)**2 + (v1z-v2z)**2 )
    
        Er = 0.5 * Run.mr * ( (vij/c_cms)**2 )
        sigma = BH.xsec(Er,react='DD_nHe',units='cm2')

        Wmin = np.min([W1,W2],0)
        Wmax = np.max([W1,W2],0)
        Fmult = np.ones(nhalf) * Run.Fmult

        Pfusion = Nfact * Fmult * Wmax * vij * sigma * Run.dt / Run.Volume
        #print 'Pfusion', Pfusion
        while np.max(Pfusion) > .99 :
            #print '***Pfusion >1'
            pindx = Pfusion>.99
            #print 'Pfusion[pindx]',Pfusion[pindx]
    
            Fmult[pindx] = Fmult[pindx] / 10.
            Pfusion[pindx] = Nfact * Fmult[pindx] * Wmax[pindx] * vij[pindx] * sigma[pindx] * Run.dt / Run.Volume
            #print 'Pfusion[pindx]',Pfusion[pindx]
            if np.min(Fmult) < 1. :
                print 'np.min(Fmult)',np.min(Fmult)
                error_exit('Fmult too low')    
                

        
        Rands = np.random.random(nhalf)
        RlP_index = np.where(Rands<Pfusion)[0]
        #print 'Pfusion',Pfusion
        #print 'Rands',Rands
        #print 'Rands<Pfusion',Rands<Pfusion
        #print 'RlP_index',RlP_index
        for i in RlP_index :
              #print 'i',i
              Wp = Wmin[i] / Fmult[i]
              v3,v4 = create_fusion_products(Run,[v1x[i],v1y[i],v1z[i]],[v2x[i],v2y[i],v2z[i]])
              W3.append(Wp)
              W4.append(Wp)
              V3.append(v3)
              V4.append(v4)

              Wyield = Wyield + Wp


        qqyield[n] = np.sum(W1*W2*vij*sigma*Run.dt/Run.Volume) * Nfact
        wwyield[n] = Wyield
 
    #print W3
    #print V3
    Y3a = np.sum(qqyield)
    Y3b = np.sum(wwyield)
    print 'Y3a',Y3a
    print 'Y3b',Y3b
    print 'T3a/Y3b',(Y3a/Y3b)

    PROD3 = fusion_product_class(W3,V3,Run.m3)
    PROD4 = fusion_product_class(W4,V4,Run.m4)


    return qqyield, wwyield, PROD3, PROD4



def do_DT_fusion(S1,S2,Run,nsteps) :

    nsteps = int(nsteps) 

    qqyield_a   = np.zeros(nsteps)
    qqyield_b   = np.zeros(nsteps)
    randyield_a = np.zeros(nsteps)
    randyield_b = np.zeros(nsteps)


    reduc_Yt = np.sum(S1.W) * np.sum(S2.W) 

    for n in range(nsteps) :
        #if n%10 == 0 :
        #    print 'n: %.0f of %.0f, %.0f percent' % (n,nsteps,float(n)/nsteps*100.)
    
        index1 = range(len(S1.W))
        index2 = range(len(S2.W))
        np.random.shuffle(index1)
        np.random.shuffle(index2)
        
        # W1 is D
        W1  = S1.W[index1]
        v1x = S1.vx[index1]
        v1y = S1.vy[index1]
        v1z = S1.vz[index1]

        # W2 is T
        W2  = S2.W[index2]
        v2x = S2.vx[index2]
        v2y = S2.vy[index2]
        v2z = S2.vz[index2]
    
        vij = np.sqrt((v1x-v2x)**2 + (v1y-v2y)**2 + (v1z-v2z)**2 )
    
        Er = 0.5 * Run.mr * ( (vij/c_cms)**2 )
        sigma = BH.xsec(Er,react='DD_nHe',units='cm2')
        #sigma = xsec(Er,react='DT',units='cm2')

        reduc_Ys = np.sum(W1*W2) 

        Pna = reduc_Yt / reduc_Ys 
        Pnb = float(len(S1.W)) # Npairs

        qqyield_a[n] = np.sum(W1*W2*vij*sigma*Run.dt/Run.Volume) * Pna
        qqyield_b[n] = np.sum(W1*W2*vij*sigma*Run.dt/Run.Volume) * Pnb

    return qqyield_a, qqyield_b        

#==============================================================================
# Fusion Kinetmatics
#==============================================================================

def create_fusion_products(Run,v1,v2) :
    # here masses are in mc2 = keV
    # Q is in keV
    # velocity is in cm/s
    va = np.array(v1)
    vb = np.array(v2)

    ##################
    # Step A1
    ##################
    u = va - vb
    vab = np.sqrt( u[0]**2 + u[1]**2 + u[2]**2 )

    ##################
    # Step A2
    ##################
    if np.sqrt(u[0]**2 + u[1]**2 ) > 0 :
        cosphi = u[0] / np.sqrt(u[0]**2 + u[1]**2 ) 
        sinphi = u[1] / np.sqrt(u[0]**2 + u[1]**2 ) 
    else :
        cosphi = 0.0
        sinphi = 1.0
    costheta = u[2] / np.sqrt(u[0]**2 + u[1]**2 + u[2]**2 )
    sintheta = np.sqrt( (u[0]**2 + u[1]**2)  / (u[0]**2 + u[1]**2 + u[2]**2) )

    M = np.array([[cosphi*costheta, sinphi*costheta, -1.*sintheta],\
                  [-1.*sinphi        , cosphi       , 0.0 ],\
                  [cosphi*sintheta, sinphi*sintheta, costheta]])
    V = np.matmul(M, u)

    ##################
    # Step A3
    ##################
    Vcm = V * (Run.m1 / (Run.m1+Run.m2) )
    Vprime = V - Vcm

    ##################
    # Step B
    ##################
    V3mag = c_cms * np.sqrt( ( 2./(Run.m3+Run.m4) ) * (Run.m4/Run.m3) * ( 0.5*Run.mr*((vab/c_cms)**2) + Run.Q ) )

    cosTHETA = 2. * np.random.random() - 1.0
    sinTHETA = np.sqrt( 1. - cosTHETA**2) 
    PHI = 2. * np.pi * np.random.random()
    cosPHI = np.cos(PHI)
    sinPHI = np.sin(PHI)

    V3prime = V3mag * np.array( [sinTHETA*cosPHI,sinTHETA*sinPHI,cosTHETA] )
    V4prime = -1.* np.sqrt(Run.m3/Run.m4) * V3prime

    ##################
    # Step C1
    ##################
    V3 = V3prime + Vcm
    V4 = V4prime + Vcm
  
    ##################
    # Step C2
    ##################
    u3 = np.matmul(M.T, V3)
    u4 = np.matmul(M.T, V4)

    ##################
    # Step C3
    ##################
    v3 = u3 + vb
    v4 = u4 + vb

    return v3,v4

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

    

class run_parameters(object):    
    
    def __init__(self,S1,S2,dt,Volume,fusiontype='DD',comparewith=1.0,savename='out',Fmult=1,do_plot1=True,savedir='.') :
        # constants
        qe = 1.60217662e-19 # coulombs
        amu2kg = 1.66054e-27 # amu-> kg
        eps0 = 8.854187817e-12 # F/m
        h = 6.62607e-34 # J*s, Planck constant
        kB = 1.60218e-19 # J/eV 
        c_ms = 2.9979e8 # m/s
        
        mc2D = 1875.613 * 1e3
        mc2T = 2809.432 * 1e3 # keV
        mc2n = 939.56536* 1e3
        #mc2n = 932.57 * 1e3
        mc2He3 = 2809.414 * 1e3
        mc2He4 = 3728.401 * 1e3
        
        if fusiontype=='DD' :
            m1 = mc2D 
            m2 = mc2D 
            m3 = mc2n
            m4 = mc2He3
            Q = 3.269 * 1e3 # keV 
            S1type = 'D'
            S2type = 'D'
        elif fusiontype=='DT' :
            m1 = mc2D 
            m2 = mc2T
            m3 = mc2n
            m4 = mc2He4
            Q = 17.589293 * 1e3 # keV From Qcalc 
            S1type = 'D'
            S2type = 'T'

        mr = (m1*m2) / (m1+m2)
        mstar = (m3/m4) * (m3+m4)
        #vn0 = c_cms * np.sqrt(2*Q/mstar)

        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4
        self.Q  = Q
        self.mr = mr
        self.mstar = mstar

        self.S1type = S1type
        self.S2type = S2type

        self.t_nsec  = 0
        self.t_step = 0
        self.t_t0  = time.time()
        self.t_run = 0

        if fusiontype == 'DD' :
            self.deltaij = 1.
        elif fusiontype == 'DT' :
            self.deltaij = 0.0


        self.Fmult = Fmult
        self.S1 = S1
        self.S2 = S2
        self.dt = dt
        self.fusiontype = fusiontype
        self.Volume = Volume
        self.comparewith = comparewith
        self.savename = savename
        self.savedir = savedir
        self.do_plot1 = do_plot1

    def increase_time(self,step) :
        self.t_nsec  += (step * self.dt)
        self.t_step += step
        self.t_run  = np.round( time.time() - self.t_t0 )

def makedirectory(mydir) :
    try:
        os.makedirs(mydir)
    except OSError, e:
        if e.errno != 17:
            raise   
        # time.sleep might help here
        pass

##def makedirectory(directory) :
##    try :
##        os.stat(directory)
##    except :
##        os.makedirs(directory)    

class sim_parameters(object):    
    
    def __init__(self,nsteps,do_mpi=False,subtime=1.) :

        self.do_mpi = do_mpi
        self.nsteps = int(nsteps)
        self.maxsubtime = subtime # in hours
        self.maxsubtime_sec = subtime * 3600. # in seconds
        self.t0 = time.time()

    def init_time(self) :
        self.t0 = time.time()
    
    def timerun_sofar(self) :
        return time.time() - self.t0 

class fusion_product_class(object):    
    
    def __init__(self,W,V,mass) :

        W = np.array(W)
        V = np.array(V)
       
        Vx = V[:,0]
        Vy = V[:,1]
        Vz = V[:,2]

        self.mc2 = mass # in keV
        self.W = W # weight in number
        self.vx = Vx # cm/2
        self.vy = Vy # cm/2
        self.vz = Vz # cm/2
        vtotal = np.sqrt(Vx**2 + Vy**2 + Vz**2)
        self.E = 0.5 * mass * (vtotal/c_cms)**2
        self.acosvxvt = np.arccos(Vx/vtotal)
        self.acosvyvt = np.arccos(Vy/vtotal)
        self.acosvzvt = np.arccos(Vz/vtotal)


def save_partcle_products(hdfname,S) :

    #print 'saving dhf5 part file'
    h5f = h5py.File(hdfname, 'w')
    #compression = 'lzf' # This compresses, but not readable by h5dump
    compression = None # larger file, but is readable by h5dump, good for testing

    h5f.create_dataset('W',   data=S.W,  compression=compression,shuffle=True)
    h5f.create_dataset('vx',  data=S.vx, compression=compression,shuffle=True)
    h5f.create_dataset('vy',  data=S.vy, compression=compression,shuffle=True)
    h5f.create_dataset('vz',  data=S.vz, compression=compression,shuffle=True)
    h5f.create_dataset('E',   data=S.E,  compression=compression,shuffle=True)
    h5f.create_dataset('mc2', data=S.mc2,compression=compression)
    h5f.close()
    #print 'FINISHED:: saving dhf5 part file'

#==============================================================================
# Plotting
#==============================================================================





def print_run_params(Run) :
   
    print 'running the sim'
    #myprint( 'nparts = %.2e' % Run.nparts )
    #myprint( 'L Debye = %.2e m' % Run.debye )
    #myprint( 'L R0 = %.2e m' % Run.r0 )
    #myprint( 'b star = %.2e m' % Run.bstar )
    #myprint( 'LnLam = %.4f ' % Run.lnLam0 )
    #
    #myprint( 'sigma star = %.2e barns' % Run.sigma_star_barns )
    #myprint( 'tau  = %.2e s' % Run.tau0 )
    #myprint( 'dt  = %.2e s' % Run.dt )
    #myprint( 'Rate  = %.2e num/sec' % Run.Rstar0 )
    #myprint( 'Rate*dt  = %.2e num/step' % Run.Pstar0 )
    #myprint( 'Rate*dt^2  = %.2e num/step' % (Run.Pstar0**2) )
    
    #myprint( 'ncpus = %.0f' % multiprocessing.cpu_count() )
    
def savematrix(x,y,filename='out.txt',x_name='x',y_name='y',x_units='[]',y_units='[]',header=''):
     f = open(filename,"w")
     dirname=os.getcwd() + '/';
     f.write('# '+ dirname+filename + "\n")
     f.write('# '+ header + "\n")
     f.write('# x='+x_name + ", y=" + y_name + "\n")
     f.write('# xunit='+x_units + ", yunit=" + y_units + "\n")
     for i in range (0,np.size(x)) :
          lineout = '%6.0f %12.6e \n' % (x[i], y[i])
          f.write(lineout)
     f.close()    

def myprint(string,fname='outputlog.txt') :
    print string
    f = open(fname,'a')
    print >>f, string
    f.close()

def error_exit(error_string) :
    print ('\n'+'*'*50)*2
    print error_string
    print ('*'*50+'\n')*2
    exit()

def warning_noexit(error_string) :
    print ('\n'+'*'*50)*2
    print error_string
    print ('*'*50+'\n')*2

#==============================================================================
# The Simulation
#==============================================================================
def run_simulation(Run,Sim) :

    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()

    #S1title = 'S1: '+Run.S1.make_title()
    #S2title = 'S2: '+Run.S2.make_title()
    #title = S1title + '\n' + S2title
    title = Run.savename 
    do_plot1 = Run.do_plot1

    if do_plot1 :
        fig,ax = pl.subplots(1,1)

    if rank == 0 :
        Sim.init_time()
        print_run_params(Run)    
        myprint( 'ncores = %.0f' % nprocs )
        
        print 'Making Particles'
        if Run.fusiontype == 'DD' :
            #W, vx, vy, vz = combine_particles(Run.S1,Run.S2)
            #S = combine_particles(Run.S1,Run.S2)
            S = Run.S1
        elif Run.fusiontype == 'DT' :
            #W, vx, vy, vz = combine_particles(Run.S1,Run.S2)
            #S1,S2 = match_particle_number(Run.S1,Run.S2)
            S1 = Run.S1
            S2 = Run.S2
           

    if Run.fusiontype == 'DD' :
        #f1,f2 = do_DD_fusion_old(S,Run,Sim.nsteps) 
        f1,f2,P3,P4 = do_DD_fusion(S,Run,Sim.nsteps) 
        savename3 = os.path.join(Run.savedir,Run.savename+'_product3.hd5')
        savename4 = os.path.join(Run.savedir,Run.savename+'_product4.hd5')
        save_partcle_products(savename3,P3)
        save_partcle_products(savename4,P4)
    elif Run.fusiontype == 'DT' :
        pass
        # need to implement particle based method
        # f1,f2 = do_DT_fusion(S1,S2,Run,Sim.nsteps) 
    
    steps = np.arange(len(f1)) + 1. 
    cf1 = np.cumsum(f1) / steps
    cf2 = np.cumsum(f2) / steps

    sv1 = cf1 / Run.comparewith
    sv2 = cf2 / Run.comparewith
    
    if do_plot1 :
        pl.sca(ax)
        pl.plot(steps,sv1,'x-',ms=2,lw=2,dashes=(20,4))
        pl.plot(steps,sv2,'s-',ms=2)
        pl.title(title,fontsize=10)    
        savefigname = os.path.join(Run.savedir,'p1_'+Run.savename+'.png')
        pl.savefig(savefigname)



    MP.plot_1D_Energy_hist(P3)
    MP.calculations(P3)



    if do_plot1 or do_plot2 :
        pl.show()





    return 1
