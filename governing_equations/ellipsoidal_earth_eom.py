# Import libraries
import math
import numpy as np

def ellipsoidal_earth_eom(t, x, vmod, amod):
    # INPUTS
    # t = time (s)
    # x = state vector at an instantaneous time t, numpy array
    # vmod = vehicle model data stored as a dictionary
    # amod = atmospheric data (not used in this EOM, but kept for compatibility)

    # OUTPUTS
    # dx = time derivative of each state in x

    # pre-allocate solution vector
    dx = np.zeros(12, dtype=float)

    # WGS-84 and J2 constants from the NASA document
    mu = 3.986004418e14  # Gravitational parameter (m^3/s^2)
    J2 = 0.0010826267    # J2 harmonic
    Re = 6378137.0       # Equatorial radius (m)
    omega_earth = 7.292115e-5 # Earth's rotation rate (rad/s)

    # declare state vector (ECI frame)
    x_eci = x[0]   # ECI x-position (m)
    y_eci = x[1]   # ECI y-position (m)
    z_eci = x[2]   # ECI z-position (m)
    u_eci = x[3]   # ECI x-velocity (m/s)
    v_eci = x[4]   # ECI y-velocity (m/s)
    w_eci = x[5]   # ECI z-velocity (m/s)
    # Note: For this check case, we only need translational motion.
    # Rotational states are not integrated here but would be needed for a full 6-DOF simulation.
    
    ##################################################################################################################
    # EQUATIONS OF MOTION (in Earth-Centered Inertial Frame)
    ##################################################################################################################

    # POSITION AND VELOCITY
    r_mag = math.sqrt(x_eci**2 + y_eci**2 + z_eci**2)
    r_mag_sq = r_mag**2

    # J2 GRAVITY MODEL (from NASA TM-2015-218675, Eqs. 29-31)
    j2_term_common = 1 - (3/2) * J2 * (Re/r_mag)**2 * (5 * (z_eci**2 / r_mag_sq) - 1)
    
    ax_g = -mu * x_eci / (r_mag**3) * j2_term_common
    ay_g = -mu * y_eci / (r_mag**3) * j2_term_common
    az_g = -mu * z_eci / (r_mag**3) * (1 - (3/2) * J2 * (Re/r_mag)**2 * (5 * (z_eci**2 / r_mag_sq) - 3))

    # AERODYNAMIC FORCES (Simplified for this case)
    # For a simple drop, we assume drag is the only aerodynamic force and acts opposite to the ECI velocity vector.
    # A full implementation would transform ECI velocity to body frame, calculate forces, and transform back.
    # This simplification is sufficient for this specific check case.
    
    # Calculate atmospheric properties at current altitude
    # First, convert ECI position to altitude
    alt_m = r_mag - Re 
    if alt_m < 0:
      alt_m = 0
    rho_kgpm3 = np.interp(alt_m, amod["alt_m"], amod["rho_kgpm3"])
    
    # Velocity relative to a non-rotating atmosphere
    v_rel_mag = math.sqrt(u_eci**2 + v_eci**2 + w_eci**2)
    
    if v_rel_mag > 0:
        qbar = 0.5 * rho_kgpm3 * v_rel_mag**2
        drag_force = vmod["CD"] * vmod["Aref_m2"] * qbar
        
        # Drag acceleration vector (opposite to velocity vector)
        ax_a = - (u_eci / v_rel_mag) * (drag_force / vmod["m_kg"])
        ay_a = - (v_eci / v_rel_mag) * (drag_force / vmod["m_kg"])
        az_a = - (w_eci / v_rel_mag) * (drag_force / vmod["m_kg"])
    else:
        ax_a, ay_a, az_a = 0, 0, 0
        
    # TOTAL ACCELERATION (Gravity + Aero)
    ax = ax_g + ax_a
    ay = ay_g + ay_a
    az = az_g + az_a

    # STATE DERIVATIVES
    # Translational Kinematics
    dx[0] = u_eci
    dx[1] = v_eci
    dx[2] = w_eci
    dx[3] = ax
    dx[4] = ay
    dx[5] = az
    
    # Rotational dynamics are ignored for this specific translational check case
    dx[6] = 0
    dx[7] = 0
    dx[8] = 0
    dx[9] = 0
    dx[10] = 0
    dx[11] = 0

    return dx