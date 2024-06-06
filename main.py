#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 17:08:19 2020

@author: stefanvet
"""

import numpy as np
import matplotlib.pyplot as plt

from bifurcation import *
from parameters import *
from dynamics import *

P = Parameters()

print(P.par_Allee)
print(P.param)
print(P.par_simpl)

""" Make phase plane of the extended LV system """
x = np.linspace(0,P.C0[0],1000)
y = np.linspace(0,P.C0[1],1000)
phasePlane(x,y,P.par_simpl,nullclines) # Steady states are given by the origin and the intersection of the nullclines.

# Simulate Allee effect with different initial densities
t,R_a1 = solveODE(P.par_Allee,dynamics_Allee,x0=.3,t_f=100.)
t,R_a2 = solveODE(P.par_Allee,dynamics_Allee,x0=.1,t_f=100.)

# Simulated the reduced mutualistic system
t,R1 = solveODE(P.par_simpl,dynamics_simpl,x0 = [.6,.6],t_f=100.)
t,R2 = solveODE(P.par_simpl,dynamics_simpl,x0 = [.1,.1],t_f=100.)

# Reduce the dilution rate. Not necessary for high Monod constants.
P = Parameters(phi=0.1)

# Simulate the mutualistic system with nutrients 
t,R_chem1 = solveODE(P.param,dynamics_chem,x0 = [.9,.9,1.,1.,1.,1.],t_f=100.)
t,R_chem2 = solveODE(P.param,dynamics_chem,x0 = [.1,.1,1.,1.,.01,.01],t_f=100.)


# Figure of the Allee effect for different initial densities.
fig,ax = plt.subplots(1,3,figsize=(8,3),tight_layout=True)
for i,R in enumerate([R_a1,R_a2]):
    ax[i].plot(t,R,c="blue")
    ax[i].set_xlabel("time t")
    ax[i].set_ylabel("Density")
    ax[i].set_ylim(0,1)
    
f_capita = perCapita_Allee(x,P.par_Allee)
ax[2].plot(x,f_capita,c="blue")
ax[2].axhline(c="red")
ax[2].set_xlabel("Density")
ax[2].set_ylabel("Per capita growth")
ax[2].set_xlim(0,1.)
plt.show()

# Figure of the reduced mutualistic system without nutrients
colors = ["blue","purple"]
fig,axes = plt.subplots(1,2,figsize=(8,3),tight_layout=True)
ax1,ax2 = axes
R = [R1,R2]
for i,ax in enumerate(axes):
    for k in range(2):
        ax.plot(t,R[i][:,k],c=colors[k])
        ax.set_xlabel("time t")
        ax.set_ylabel("Density")
        ax.set_ylim(0,1)
plt.show()

fig,axes = plt.subplots(2,2,figsize=(8,6),tight_layout=True)

# Figure of the complete mutualistic system with nutrients
colors=["blue","purple","red","green"]
for i,R in enumerate([R_chem1,R_chem2]):
    Rx = R[:,:2] # Species density
    Rs = R[:,2:] # Nutrient concentrations
    for j,Ri in enumerate([Rx,Rs]):
        ax = axes[j,i]
        ax.set_xlabel("time t")
        ax.set_ylim(0,1.5)
        for k in range(Ri.shape[1]):
            ax.plot(t,Ri[:,k],c=colors[k])
    axes[0,i].set_ylabel("Density")
    axes[1,i].set_ylabel("Concentration")
plt.title("Chemostat equations")
plt.show()
