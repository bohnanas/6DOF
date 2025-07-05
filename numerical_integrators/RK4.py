import numpy as np

def rk4(f, t, x, h, vmod, amod):
    ######################################################################################
    # Fourth-order Runge-Kutta ODE solver.
    
    # Parameters:
    # f: governing equations
    # t: time
    # x: state
    # h: time step size
    # vmod: vehicle model
    # amod: atmospheric model
    
    # Returns:
    # t: time
    # x: state solution
    #######################################################################################

    for i in range(1, len(t)):
        k1 = h * f(t[i-1], x[:,i-1], vmod, amod)
        k2 = h * f(t[i-1] + 0.5 * h, x[:,i-1] + 0.5 * k1, vmod, amod)
        k3 = h * f(t[i-1] + 0.5 * h, x[:,i-1] + 0.5 * k2, vmod, amod)
        k4 = h * f(t[i-1] + h, x[:,i-1] + k3, vmod, amod) 
        x[:,i] = x[:,i-1] + 1.0/6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

    return t, x
