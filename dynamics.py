#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 17:21:48 2020

@author: stefanvet
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi



def solveODE(par,dynamics,x0=np.ones(2),t_0=0.,t_f=50.,Dt=0.1,write=False):
    """ Solve ODE for parameters par """
    if write: print(par)
        
    ode = spi.ode(dynamics)
    ode.set_f_params(par)
    # BDF method suited to stiff systems of ODEs
    ode.set_integrator('vode',nsteps=10000,method='bdf')
    
    ode.set_initial_value(x0,t_0)

    ts = []
    ys = []

    # Integrate
    while (ode.successful() and ode.t < t_f):
        ode.integrate(ode.t + Dt)
        
        ts.append(ode.t)
        ys.append(ode.y)
    t_grid = np.vstack(ts)
    R = np.vstack(ys)
    return [t_grid,R]

def dynamics_Allee(t,y,par):
    dyn = (par["r"]*(par["C"]-y)*(par["a"]+y) - par["d"])*y
    return dyn    

def perCapita_Allee(y,par):
    dyn = (par["r"]*(par["C"]-y)*(par["a"]+y) - par["d"])
    return dyn    

def growthrate(x,s,p,par):
    """ growth rate of the chemostat system """
        
    mu = par["mu"]
    K0 = par["K0"]
    Kp = par["Kp"]
    
    rate = mu*s/(K0 + s)*p/(Kp + p)
    
    return rate

def dynamics_chem(t,y,par):
    dim = len(y)
    
    """ dynamics of the chemostat system """    
    phi = par["phi"]
    
    St = par["St"]
    Pt = par["Pt"]
    
    v0 = par["v0"]
    vp = par["vp"]
    ap = par["ap"]
    
    x = y[:2]
    s = y[2:4]
    p = y[4:]
    growth = growthrate(x,s,p,par)

    dyn = np.zeros(dim)
    dyn[:2] = (growth - phi)*x
    dyn[2:4] = phi*(St - s) - v0*growth*x
    dyn[4:] = phi*(Pt - p) - vp*growth*x + ap*np.flip(growth*x)
    return dyn

def dynamics_simpl(t,y,par):
    r=par["r"]
    C = par["C"]
    b = par["b"]
    C0 = par["C0"]
    d = par["d"]
    
    prod = (C + b*np.flip(y) - y)
    cons = C0 - y    
    
    dyn = r*cons*prod*y - d*y
    
    return dyn

def nullclines(x1,x2,par):
    r=par["r"]
    C = par["C"]
    b = par["b"]
    C0 = par["C0"]
    d = par["d"]
        
    prod = np.array([C[0] + b[0]*x2 - x1,C[1] + b[1]*x1 - x2])
    cons = np.array([C0[0] - x1,C0[1] - x2])
    
    r1,r2 = r
    eq1 = r1*prod[0]*cons[0] - d
    eq2 = r2*prod[1]*cons[1] - d
    return [eq1,eq2]


def phasePlane(x1,x2,par,nullclineFunc,name="phasePlane.png",save=False):
    print(par)
    """ Figure of the phase plane """
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    
    X,Y = np.meshgrid(x1,x2)

    NC = nullclineFunc(X,Y,par)
    
    for NCi in NC:
        ax.contour(X,Y,NCi,[0],linewidths=1.5)
    
    ax.set_xlim(0,x1.max())
    ax.set_ylim(0,x2.max())
    ax.set_xlabel(r'$\rho_1$')
    ax.set_ylabel(r"$\rho_2$")
    ax.grid()
    
    plt.show()
    if save: fig.savefig(name)
    return