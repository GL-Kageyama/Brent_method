#===================================================================================
#--------------------------       Brent's method       -----------------------------
#===================================================================================
#Brent's method is a root-finding algorithm combining the bisection method, 
#the secant method and inverse quadratic interpolation.
#===================================================================================
# Pattern 1 : Basic implementation in Python
#===================================================================================

import copy
import numpy as np
import sys

def get_solution_by_bisection_method(b, c):
    return (c + b)/2

def get_solution_by_secant_method(a, b, fa, fb):
    return b - fb * (b - a) / (fb -fa)

def get_solution_by_inv_quad_interpolation(a, b, c, fa, fb, fc):
    d_pre = (fa - fb) / (a - b)
    d_blk = (fc - fb) / (c - b)
    
    return b - fb * (fc * d_blk - fa * d_pre) / (d_blk * d_pre * (fc - fa))

def brent_method(func, init_x1, init_x2, xtol, max_iter):
    
    f1 = func(init_x1)
    f2 = func(init_x2)
    
    if f1 * f2 > 0:
        print("The two specified points do not surround the root.")
        return
    
    if f1 == 0:
        return [init_x1, f1]
    if f2 == 0:
        return [init_x2, f2]
    
    a = init_x1
    b = init_x2
    c = copy.copy(a)
    fa = func(a)
    fb = func(b)
    fc = func(c)
    
    if np.abs(fc) < np.abs(fb):
        a = copy.copy(b)
        b = copy.copy(c)
        c = copy.copy(a)
        
        fa = copy.copy(fb)
        fb = copy.copy(fc)
        fc = copy.copy(fa)
        
    meps = sys.float_info.epsilon
    
    for n_iter in range(max_iter):
        delta = 0.5 * xtol + 2. * meps * np.abs(b)
        
        interval_hdist = (c - b) / 2
        improve_dist = b - a
        
        if (fb==0) | (np.abs(interval_hdist) < delta):
            return [b, fb]
        
        if (np.abs(improve_dist) > delta) & (np.abs(fb) < np.abs(fa)):
            
            if a==c:
                s = get_solution_by_secant_method(a, b, fa, fb)
                
            else:
                s = get_solution_by_inv_quad_interpolation(a, b, c, fa, fb, fc)
                
            tmp_iprove_dist = s - b
            
            if not 2 * np.abs(tmp_iprove_dist) < np.min([np.abs(improve_dist), 3 * np.abs(interval_hdist) - delta]):
                s = get_solution_by_bisection_method(b. c)
                
        else:
            s = get_solution_by_bisection_method(b, c)
            
        a = copy.copy(b)
        fa = copy.copy(fb)
        
        if np.abs(tmp_iprove_dist) > delta:
            b = s
            
        else:
            if tmp_iprove_dist > 0:
                b += delta
            else:
                b -= delta
                
        fb = func(b)
        print([n_iter, b, fb])
        
        if fa * fb < 0:
            c = copy.copy(a)
            fc = copy.copy(fa)
            
        if np.abs(fc) < np.abs(fb):
            a = copy.copy(b)
            b = copy.copy(c)
            c = copy.copy(a)
            
            fa = copy.copy(fb)
            fb = copy.copy(fc)
            fc = copy.copy(fa)
            
    return [b, fb]

def test_func(x):
    return x**2 - 100 * x - 10

init_xa = -100
init_xb = 100
max_iter = 100
print("===== Basic implementation in Python =====")
print(brent_method(test_func, init_xa, init_xb, xtol=10e-12, max_iter=max_iter))
print("")

#===================================================================================
# Pattern 2 : Brent's method using Scipy's API
#===================================================================================

from scipy.optimize import root_scalar 

sol = root_scalar(test_func, bracket=[init_xa, init_xb], method="brentq")
print("===== Brent's method using Scipy's API =====")
print("root: {0} interation: {1}".format(sol.root, sol.iterations))
print("")

