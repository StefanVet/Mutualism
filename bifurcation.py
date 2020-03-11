#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 09:07:44 2019

@author: Stefan
"""




# Estimates the Jacobian of the system using h (time t).     
def jacobian(y,func,par):
    

    h = 10**(-3)
    dim = len(y)
    
    J = np.zeros((dim,dim))	# Initialization
    for j in range(dim):

        y_0 = y[j]	# Center of derivative evaluated between [y_0-h,y_0+h]
        y[j] += h

        dyn1 = func(0,y,par)	# !Array of size N! (determine rows of jacobian J)

        y[j] = y_0 - h
        dyn2 = func(0,y,par)	# !Array of size N!

        J[:,j] = (dyn1[:] - dyn2[:])/(2.*h)

    return J

def newtonRaphson(yN,dynamics,par,write=False):
    """ algorithm: y = y - J^(-1)*f(y) """
    count = 0
    converge = False 
    if write: print(par)
    f_y = dynamics(0,yN,par)

    while(not converge and count<1000):
        f_y = dynamics(0,yN,par)
        jacob = jacobian(yN.copy(),dynamics,par)
        try: 
            inv_jacob = la.inv(jacob)
        except:
            print("\nSingular matrix")
            print("\nLast evaluation: y: %r\nf_y: %r,jacob: %r"%(yN,f_y,jacob))
            break      
        
        dyN = np.dot(inv_jacob,f_y)
        yN = yN - dyN  
        yN[yN<0] = 10**(-14)
        # Check if condition is met
        if(np.all(np.abs(f_y)<10.**(-10))):
            converge = True
        else: count += 1      
    if(not converge):
        print("\nNo convergence after %r iterations"%count)
        print("\nLast evaluation: y: %r\nf_y: %r,jacob: %r"%(yN,f_y,jacob))
        sys.exit()
    elif(write): 
        print(yN,"yN")
        print("\nConvergence reached after %r iterations"%count)
        print("\nLast evaluation: y: %r\nf_y: %r,jacob: %r"%(yN,f_y,jacob))
    return [yN,converge,count]


def bifurcation(paramRange,y_st,par,func,paramName='phi',write=False):
        

    # empty lists for cross-section and corresponding eigenvalues
    y_cross_list = []
    eig_list = []
    
    
    for i,p in enumerate(paramRange):
        par[paramName] = p
        # DO newton-raphson
        
        try: 
            y_st,conv,cnts = newtonRaphson(y_st,func,par,write=False)
        except: 
            print("no convergence at "+paramName + "="+str(p))
            conv = False      	
        
        y_cross_list.append(y_st)          

        # Calculate eigenvalues
        if conv: 
            jacob = jacobian(y_st.copy(),func,par)
            eig = la.eig(jacob)[0]
            eig_list.append(eig.real)
        elif not conv: eig_list.append(eig_list[-1])
        
        if write: print("p",paramName,par[paramName], "y_st", y_st, "eig",eig)
        
    eig_list = np.array(eig_list)
    y_cross_list = np.array(y_cross_list)
    
    return [y_cross_list,eig_list]
