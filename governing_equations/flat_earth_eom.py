# Import libraries
import math
import numpy as np

def flat_earth_eom(t, x, vmod, amod):
    # INPUTS    
    # t = time (s)
    # x = state vector at an instantaneous time t, numpy array
    # vmod = aircraft model data stored as a dictionary
    # amod = atmospheric data

    # OUTPUTS
    # dx = time derivative of each state in x

    # pre-allocate solution vector
    dx = np.zeros(12, dtype=float)

    # declare state vector
    u_b_mps = x[0]   # body x-axis velocity
    v_b_mps = x[1]   # body y-axis velocity
    w_b_mps = x[2]   # body z-axis velocity
    p_b_rps = x[3]   # roll angular velocity
    q_b_rps = x[4]   # pitch angular velocity
    r_b_rps = x[5]   # yaw angular velocity
    phi_rad = x[6]   # roll euler angle
    theta_rad = x[7] # pitch euler angle
    psi_rad = x[8]   # yaw euler angle
    p1_NED_m = x[9]  # x-axis position (NED Frame)
    p2_NED_m = x[10] # y-axis position (NED Frame)
    p3_NED_m = x[11] # z-axis position (NED Frame)

    # get mass and intertial parameters from vehicle model dictionary
    m_kg = vmod["m_kg"]
    Ixx_b_kgm2 = vmod["Ixx_b_kgm2"] # must be positive
    Iyy_b_kgm2 = vmod["Iyy_b_kgm2"] # must be positive
    Izz_b_kgm2 = vmod["Izz_b_kgm2"] # must be positive
    Ixz_b_kgm2 = vmod["Ixz_b_kgm2"]
    Ixy_b_kgm2 = vmod["Ixy_b_kgm2"] # symmetric about x-z plane? = 0
    Iyz_b_kgm2 = vmod["Iyz_b_kgm2"] # symmetric about x-z plane? = 0

    # Atmospheric model
    h_m = -p3_NED_m # current altitude
    rho_kgpm3 = np.interp(h_m, amod["alt_m"], amod["rho_kgpm3"]) # density as a function of altitude
    # rho_kgpm3 = 1.2 # density as a constant (sea-level)
    c_mps = np.interp(h_m, amod["alt_m"], amod["c_mps"], h_m)

    # AoA, Sideslip
    true_airspeed_mps = math.sqrt(u_b_mps**2 + v_b_mps**2 + w_b_mps**2)
    qbar_kgpms2 = 0.5 * rho_kgpm3 * true_airspeed_mps**2  

    alpha_rad = np.atan2(w_b_mps, u_b_mps) 
    beta_rad = np.asin(v_b_mps / true_airspeed_mps)

    # Lift, Drag, Sideforce, Thrust
    lift_N = 0 # (units = Newton)
    drag_N = vmod["CD"] * qbar_kgpms2 * vmod["Aref_m2"] # (units = Newton)
    side_N = 0 # (units = Newton)
    thrust_bx_N = 0 # thrust in body x axis
    thrust_by_N = 0 # thrust in body y axis
    thrust_bz_N = 0 # thrust in body z axis

    # gravity (assuming constant)
    # gz_n_mps2 = 9.81

    # gravity (calculating as a function of altitude)
    gz_n_mps2 = np.interp(h_m, amod["alt_m"], amod["g_mps2"])

    # external forces
    X_b_N = thrust_bx_N + lift_N * math.sin(alpha_rad) - drag_N * math.cos(alpha_rad) * math.cos(beta_rad) # (units = Newton)
    Y_b_N = thrust_by_N - drag_N * math.sin(beta_rad) + side_N * math.cos(beta_rad) # (units = Newton)
    Z_b_N = thrust_bz_N - lift_N * math.cos(alpha_rad) - drag_N * math.sin(alpha_rad) * math.cos(beta_rad) - side_N * math.sin(alpha_rad) * math.sin(beta_rad) # (units = Newton)

    # external moments
    L_b_Nm = 0 # (units = Newton Meter)
    M_b_Nm = 0 # (units = Newton Meter)
    N_b_Nm = 0 # (units = Newton Meter)

    ##################################################################################################################
    # EQUATIONS OF MOTION
    ##################################################################################################################

    # TRANSLATIONAL EQUATIONS
    # u dot
    dx[0] = (r_b_rps * v_b_mps) - (q_b_rps * w_b_mps) - (gz_n_mps2 * math.sin(theta_rad)) + (X_b_N / m_kg)

    # v dot
    dx[1] = (p_b_rps * w_b_mps) - (r_b_rps * u_b_mps) + (gz_n_mps2 * math.sin(phi_rad) * math.cos(theta_rad)) + (Y_b_N / m_kg)

    # w dot
    dx[2] = (q_b_rps * u_b_mps) - (p_b_rps * v_b_mps) + (gz_n_mps2 * math.cos(phi_rad) * math.cos(theta_rad)) + (Z_b_N / m_kg)

    # ROTATIONAL EQUATIONS
    denom = Ixx_b_kgm2 * Izz_b_kgm2 - Ixz_b_kgm2**2
    # p dot
    dx[3] = (1 / denom) * (Izz_b_kgm2 * L_b_Nm + Ixz_b_kgm2 * N_b_Nm + ((Ixz_b_kgm2 * p_b_rps * q_b_rps) * (Ixx_b_kgm2 - Iyy_b_kgm2 + Izz_b_kgm2)) - ((q_b_rps * r_b_rps) * (Izz_b_kgm2**2 - Izz_b_kgm2 * Iyy_b_kgm2 + Ixz_b_kgm2**2)))

    # q dot
    dx[4] = (1 / Iyy_b_kgm2) * (M_b_Nm - ((r_b_rps * p_b_rps) * (Ixx_b_kgm2 - Izz_b_kgm2)) - ((Ixz_b_kgm2) * (p_b_rps**2 - r_b_rps**2)))

    # r dot
    dx[5] = (1 / denom) * (Ixz_b_kgm2 * L_b_Nm + Ixx_b_kgm2 * N_b_Nm - ((Ixz_b_kgm2 * r_b_rps * q_b_rps) * (Ixx_b_kgm2 - Iyy_b_kgm2 + Izz_b_kgm2)) - ((p_b_rps * r_b_rps) * (Ixx_b_kgm2**2 - Ixx_b_kgm2 * Iyy_b_kgm2 + Ixz_b_kgm2**2)))

    # EULER KINEMATIC EQUATIONS
    # phi dot
    dx[6] = p_b_rps + q_b_rps * math.tan(theta_rad) * math.sin(phi_rad) + r_b_rps * math.tan(theta_rad) * math.cos(phi_rad)

    # theta dot
    dx[7] = q_b_rps * math.cos(phi_rad) - r_b_rps * math.sin(phi_rad)

    # psi dot    
    dx[8] = q_b_rps * (math.sin(phi_rad) / math.cos(theta_rad)) + r_b_rps * (math.cos(phi_rad) / math.cos(theta_rad))

    # NAVIGATION EQUATIONS
    # x dot
    dx[9] = u_b_mps * math.cos(theta_rad) * math.cos(psi_rad) + (v_b_mps * (math.sin(phi_rad) * math.sin(theta_rad) * math.cos(psi_rad) - math.cos(phi_rad) * math.sin(psi_rad))) + (w_b_mps * (math.cos(phi_rad) * math.sin(theta_rad) * math.cos(psi_rad) + math.sin(phi_rad) * math.sin(psi_rad)))

    # y dot
    dx[10] = u_b_mps * math.cos(theta_rad) * math.sin(psi_rad) + (v_b_mps * (math.sin(phi_rad) * math.sin(theta_rad) * math.sin(psi_rad) + math.cos(phi_rad) * math.cos(psi_rad))) + (w_b_mps * (math.cos(phi_rad) * math.sin(theta_rad) * math.sin(psi_rad) - math.sin(phi_rad) * math.cos(psi_rad)))

    # z dot
    dx[11] = -u_b_mps * math.sin(theta_rad) + v_b_mps * math.sin(phi_rad) * math.cos(theta_rad) + w_b_mps * math.cos(phi_rad) * math.cos(theta_rad)

    return dx