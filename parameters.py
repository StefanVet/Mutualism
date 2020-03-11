#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 09:17:09 2019

@author: Stefan
"""

""" Chemostat parameters """
import numpy as np

class Parameters(object):
    
    def __init__(self,**params_opt):
        
        # Dilution rate
        self.phi=0.2
        
        # consumption 1st resource
        mu1=2.
        St1=1.
        K01=2.
        v01=1.
        
        # consumption 2st resource
        mu2=2.
        St2=1.
        K02=2.
        v02=1.
        
        # consumption 1st cross-feeding nutrient
        Pt1 = 0.1
        a1 = 2.
        
        # consumption 1st cross-feeding nutrient
        K1 = 1.
        v1 = 1.
        
        # consumption 2nd cross-feeding nutrient
        Pt2 = 0.1
        a2 = 2.
        
        # consumption 2nd cross-feeding nutrient
        K2 = 1.
        v2 = 1.
        
        self.mu = np.array([mu1,mu2])
        self.St = np.array([St1,St2])
        self.Pt = np.array([Pt1,Pt2])
        self.K0 = np.array([K01,K02])
        self.Kp = np.array([K1,K2])
        self.v0 = np.array([v01,v02])
        self.vp = np.array([v1,v2])
        self.ap = np.array([a1,a2])
        
        self.param = dict(mu=self.mu,St=self.St,Pt=self.Pt,
                                  K0=self.K0,Kp=self.Kp,v0=self.v0,vp=self.vp,
                                  ap=self.ap,phi=self.phi)
        
        for key,value in params_opt.items():
            self.param[key] = value
            
        self.mu = self.param["mu"]
        self.St = self.param["St"]
        self.Pt = self.param["Pt"]
        self.K0 = self.param["K0"]
        self.Kp = self.param["Kp"]
        self.v0 = self.param["v0"]
        self.vp = self.param["vp"]
        self.ap = self.param["ap"]
        self.phi = self.param["phi"]        
        """ Reduced parameters """    
    
        r1 = self.mu[0]*self.v0[0]*self.vp[0]/(self.K0[0]*self.Kp[0])
        C01 = self.St[0]/self.v0[0]
        C1 = self.Pt[0]/self.vp[0]
        b1 = self.ap[0]/self.vp[0]
        
        r2 = self.mu[1]*self.v0[1]*self.vp[1]/(self.K0[1]*self.Kp[1])
        C02 = self.St[1]/self.v0[1]
        C2 = self.Pt[1]/self.vp[1]
        b2 = self.ap[1]/self.vp[1]
        
        self.r = np.array([r1,r2])
        self.C0 = np.array([C01,C02])
        self.Cp = np.array([C1,C2])
        self.b = np.array([b1,b2])
        
        #par_simpl = dict(r1=r1,r2=r2,C0=C0,C1=C1,C2=C2,b1=b1,b2=b2,d=d)
        self.par_simpl = dict(r=self.r,C=self.Cp,
                              b=self.b,C0=self.C0,d=self.phi)
            
            
        """ Allee parameters """
        self.a = np.average(self.Cp)/(np.average(self.b) - 1)
        self.ra = np.average(self.r)*(np.average(self.b) - 1)
        self.Ca = np.average(self.C0)

        
        self.par_Allee = dict(r=self.ra,C=self.Ca,a=self.a,d=self.phi)
    
