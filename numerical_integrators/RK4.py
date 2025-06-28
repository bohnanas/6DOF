def rk4(f, y0, t0, t_end, h, vmod):
    """
    Fourth-order Runge-Kutta ODE solver.
    
    Parameters:
    - f: function f(t, y) representing dy/dt
    - y0: initial value of y (can be scalar or NumPy array)
    - t0: initial time
    - t_end: final time
    - h: time step
    
    Returns:
    - t_values: list of time points
    - y_values: list of y values at each time step
    """
    t = t0
    y = y0
    t_values = [t]
    y_values = [y]

    while t < t_end:
        if t + h > t_end:  # adjust final step to hit t_end exactly
            h = t_end - t

        k1 = f(t, y, vmod)
        k2 = f(t + h/2, y + h/2 * k1, vmod)
        k3 = f(t + h/2, y + h/2 * k2, vmod)
        k4 = f(t + h, y + h * k3, vmod)

        y = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        t = t + h

        t_values.append(t)
        y_values.append(y)

    return t_values, y_values
