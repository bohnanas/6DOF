# Import libraries
import math
import numpy as np

# Coefficient tables for the X15
from vehicle_models.X15.aero import CDrag_X15
from vehicle_models.X15.aero import CLift_X15
from vehicle_models.X15.aero import CSide_X15
from vehicle_models.X15.aero import CRoll_X15
from vehicle_models.X15.aero import CPitch_X15
from vehicle_models.X15.aero import CYaw_X15

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
    aux_data = np.zeros(4, dtype=float)

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

    # Reference dimensions
    aref_m2 = vmod["Aref_m2"]
    b_m = vmod["b_m"]
    c_m = vmod["c_m"]

    # Drag Coefficients
    cd_wb = vmod["CD_wb"]
    cd_elev_pdeg = vmod["CD_elev_pdeg"]
    cd_sb_pdeg = vmod["CD_sb_pdeg"]

    # Lift coefficients
    cl_wb = vmod["CL_wb"]
    cl_alpha_pdeg = vmod["CL_alpha_pdeg"]
    cl_elev_pdeg = vmod["CL_elev_pdeg"]

    # Side force coefficients
    cy_beta_prad = vmod["CY_beta_prad"]
    cy_p_prps = vmod["CY_p_prps"]
    cy_r_prps = vmod["CY_r_prps"]
    cy_betadot_pdps = vmod["CY_betadot_pdps"]
    cy_ail_pdeg = vmod["CY_ail_pdeg"]
    cy_rud_pdeg = vmod["CY_rud_pdeg"]

    # Roll moment coefficients
    cl_beta_prad    = vmod["Cl_beta_prad"]
    cl_p_prps       = vmod["Cl_p_prps"]
    cl_r_prps       = vmod["Cl_r_prps"]
    cl_betadot_pdps = vmod["Cl_betadot_pdps"]
    cl_ail_pdeg    = vmod["Cl_ail_pdeg"]
    cl_rud_pdeg    = vmod["Cl_rud_pdeg"]
  
    # Pitch moment coefficients
    cm_wb            = vmod["Cm_wb"]
    cm_alpha_pdeg    = vmod["Cm_alpha_pdeg"]
    cm_q_prps        = vmod["Cm_q_prps"]
    cm_alphadot_pdps = vmod["Cm_alphadot_pdps"]
    cm_elev_pdeg     = vmod["Cm_elev_pdeg"]
    cm_sb_pdeg    = vmod["Cm_sb_pdeg"]

    # Yaw moment coefficients
    cn_beta_pdeg    = vmod["Cn_beta_pdeg"]
    cn_p_prps       = vmod["Cn_p_prps"]
    cn_r_prps       = vmod["Cn_r_prps"]
    cn_betadot_pdps = vmod["Cn_betadot_pdps"]
    cn_ail_pdeg    = vmod["Cn_ail_pdeg"]
    cn_rud_pdeg    = vmod["Cn_rud_pdeg"]

    # Control input data (aileron doublet)
    ail_deg = 0
    elev_deg = 0
    rud_deg = 0

    if t < 1:
        ail_deg = 0
    elif t < 3:
        ail_deg = -0.25*0
    elif t < 5:
        ail_deg = 0.25*0
    else:
        ail_deg = 0
  
    if t < 1:
        elev_deg = 0
    elif t < 3:
        elev_deg = -2*0
    elif t < 5:
        elev_deg = 2*0
    else:
        elev_deg = 0
  
    if t < 1:
        rud_deg = 0
    elif t < 3:
        rud_deg = -0.2
    elif t < 5:
        rud_deg = 0.2
    else:
        rud_deg = 0

    sb_deg = 0
    thrust_percent = 0

    # Atmospheric model
    h_m = -p3_NED_m # current altitude
    rho_kgpm3 = np.interp(h_m, amod["alt_m"], amod["rho_kgpm3"]) # density as a function of altitude
    c_mps = np.interp(h_m, amod["alt_m"], amod["c_mps"], h_m)

    # Air data
    true_airspeed_mps = math.sqrt(u_b_mps**2 + v_b_mps**2 + w_b_mps**2)
    qbar_kgpms2 = 0.5 * rho_kgpm3 * true_airspeed_mps**2  
    mach = true_airspeed_mps / c_mps
    alpha_rad = np.atan2(w_b_mps, u_b_mps)
    beta_rad = np.asin(v_b_mps / true_airspeed_mps)
    alpha_deg = alpha_rad * (180 / math.pi)
    beta_deg = beta_rad * (180 / math.pi)

    alphadot_dps = 0
    betadot_dps = 0

    # Lift, Drag, Sideforce, Thrust
    cl = CLift_X15(cl_wb, cl_alpha_pdeg, cl_elev_pdeg, alpha_deg, elev_deg)
    cd = CDrag_X15(cd_wb, cd_elev_pdeg, cd_sb_pdeg, elev_deg, sb_deg)
    cy = CSide_X15(cy_beta_prad, cy_p_prps, cy_r_prps, cy_betadot_pdps, cy_ail_pdeg, cy_rud_pdeg, beta_rad, p_b_rps, r_b_rps, betadot_dps, ail_deg, rud_deg, true_airspeed_mps, b_m)

    lift_N = cl * qbar_kgpms2 * aref_m2
    drag_N = cd * qbar_kgpms2 * aref_m2
    side_N = cy * qbar_kgpms2 * aref_m2
    thrust_bx_N = 0 # thrust in body x axis
    thrust_by_N = 0 # thrust in body y axis
    thrust_bz_N = 0 # thrust in body z axis

    # gravity (calculating as a function of altitude)
    gz_n_mps2 = np.interp(h_m, amod["alt_m"], amod["g_mps2"])

    # external forces
    X_b_N = thrust_bx_N + lift_N * math.sin(alpha_rad) - drag_N * math.cos(alpha_rad) * math.cos(beta_rad) 
    Y_b_N = thrust_by_N - drag_N * math.sin(beta_rad) + side_N * math.cos(beta_rad) 
    Z_b_N = thrust_bz_N - lift_N * math.cos(alpha_rad) - drag_N * math.sin(alpha_rad) * math.cos(beta_rad) - side_N * math.sin(alpha_rad) * math.sin(beta_rad)

    # external moments
    cL = CRoll_X15(cl_beta_prad, cl_p_prps, cl_r_prps, cl_betadot_pdps, cl_ail_pdeg, cl_rud_pdeg, beta_rad, betadot_dps, p_b_rps, r_b_rps, ail_deg, rud_deg, true_airspeed_mps, b_m)
    cM = CPitch_X15(cm_wb, cm_alpha_pdeg, cm_q_prps, cm_alphadot_pdps, cm_elev_pdeg, cm_sb_pdeg, alpha_deg, alphadot_dps, q_b_rps, elev_deg, sb_deg, true_airspeed_mps, c_m)
    cN = CYaw_X15(cn_beta_pdeg, cn_p_prps, cn_r_prps, cn_betadot_pdps, cn_ail_pdeg, cn_rud_pdeg, beta_rad, betadot_dps, p_b_rps, r_b_rps, ail_deg, rud_deg, true_airspeed_mps, b_m)

    L_b_Nm = cL * qbar_kgpms2 * aref_m2 * b_m
    M_b_Nm = cM * qbar_kgpms2 * aref_m2 * c_m
    N_b_Nm = cN * qbar_kgpms2 * aref_m2 * b_m

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
    dx[5] = (1 / denom) * (Ixz_b_kgm2 * L_b_Nm + Ixx_b_kgm2 * N_b_Nm - ((Ixz_b_kgm2 * r_b_rps * q_b_rps) * (Ixx_b_kgm2 - Iyy_b_kgm2 + Izz_b_kgm2)) + ((p_b_rps * q_b_rps) * (Ixx_b_kgm2**2 - Ixx_b_kgm2 * Iyy_b_kgm2 + Ixz_b_kgm2**2)))

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

    ##################################################################################################################
    # CONTROL SURFACE DATA
    ##################################################################################################################

    aux_data[0] = ail_deg
    aux_data[1] = elev_deg
    aux_data[2] = rud_deg
    aux_data[3] = thrust_percent

    return dx, aux_data