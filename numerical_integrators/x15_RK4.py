import numpy as np

def rk4(f, t, x, h, vmod, amod, aux):
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
        pk1, aux_data = f(t[i-1], x[:,i-1], vmod, amod)
        k1 = h * pk1
        pk2, _ = f(t[i-1] + 0.5 * h, x[:,i-1] + 0.5 * k1, vmod, amod)
        k2 = h * pk2
        pk3, _ = f(t[i-1] + 0.5 * h, x[:,i-1] + 0.5 * k2, vmod, amod)
        k3 = h * pk3
        pk4, _ = f(t[i-1] + h, x[:,i-1] + k3, vmod, amod) 
        k4 = h * pk4
        x[:,i] = x[:,i-1] + 1.0/6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        aux[:,i] = aux_data

    return t, x, aux
