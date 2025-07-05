import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math 

"""validation_cases.py verifies the python 6DOF simulation against 
the various NASA code verification cases. The approach is to load 
data for each case and create plots comparing the results of the 
various simulations. 

The NASA verification data contains the following columns in CSV form
with the following variable names. The reader is referred to 

https://nescacademy.nasa.gov/flightsim for further documentation and all data verification files.

time
gePosition_ft_X
gePosition_ft_Y
gePosition_ft_Z
feVelocity_ft_s_X
feVelocity_ft_s_Y
feVelocity_ft_s_Z
altitudeMsl_ft
longitude_deg
latitude_deg
localGravity_ft_s2
eulerAngle_deg_Yaw
eulerAngle_deg_Pitch
eulerAngle_deg_Roll
bodyAngularRateWrtEi_deg_s_Roll
bodyAngularRateWrtEi_deg_s_Pitch
bodyAngularRateWrtEi_deg_s_Yaw
altitudeRateWrtMsl_ft_min
speedOfSound_ft_s
airDensity_slug_ft3
ambientPressure_lbf_ft2
ambientTemperature_dgR
aero_bodyForce_lbf_X
aero_bodyForce_lbf_Y
aero_bodyForce_lbf_Z
aero_bodyMoment_ftlbf_L
aero_bodyMoment_ftlbf_M
aero_bodyMoment_ftlbf_N
mach
dynamicPressure_lbf_ft2
trueAirspeed_nmi_h

# Python simulation column order
# 0.  t_s 
# 1.  u_b_mps, axial velocity of CM wrt inertial CS resolved in aircraft body fixed CS
# 2.  v_b_mps, lateral velocity of CM wrt inertial CS resolved in aircraft body fixed CS
# 3.  w_b_mps, vertical velocity of CM wrt inertial CS resolved in aircraft body fixed CS
# 4.  p_b_rps, roll angular velocity of body fixed CS with respect to inertial CS
# 5.  q_b_rps, pitch angular velocity of body fixed CS with respect to inertial CS
# 6.  r_b_rps, yaw angular velocity of body fixed CS with respect to inertial CS
# 7.  phi_rad, roll angle
# 8.  theta_rad, pitch angle
# 9.  psi_rad, yaw angle
# 10. p1_n_m, x-axis position of aircraft resolved in NED CS
# 11. p2_n_m, y-axis position of aircraft resolved in NED CS
# 12. p3_n_m, z-axis position of aircraft resolved in NED CS
# 13. Altitude_m, altitude (= -p3_n_m) (row vector)
# 14. Cs_mps, speed of sound interpolated from atmosphere model (row vector)
# 15. Rho_kgpm3, Air density interpolated from atmosphere model (row vector)
# 16. Mach, Mach number (row vector)
# 17. Alpha_rad, Angle of attack (row vector)
# 18. Beta_rad, Angle of sideslip (row vector)
# 19. True_Airspeed_mps, True airspeed (row vector)
# 20. C_phi, cosine of roll angle (row vector)
# 21. C_theta, cosine of pitch angle (row vector)
# 22. C_psi, cosine of yaw angle (row vector)
# 23. S_phi, sine of roll angle (row vector)
# 24. S_theta, sine of pitch angle (row vector)
# 25. S_psi, sine of yaw angle (row vector)
# 26. T_theta, tangent of pitch angle (row vector)
#

"""

# Atmospheric case 01: Dropped Sphere, dragless (WGS-84 earth)
# Read the CSV file, specifying the header row
data_02_sim_01 = pd.read_csv('./NASA_data/trajectory_data/Atmos_02_TumblingBrickNoDamping/Atmos_02_sim_01.csv', header=0)
data_02_sim_02 = pd.read_csv('./NASA_data/trajectory_data/Atmos_02_TumblingBrickNoDamping/Atmos_02_sim_02.csv', header=0)
data_02_sim_04 = pd.read_csv('./NASA_data/trajectory_data/Atmos_02_TumblingBrickNoDamping/Atmos_02_sim_04.csv', header=0)
data_02_sim_05 = pd.read_csv('./NASA_data/trajectory_data/Atmos_02_TumblingBrickNoDamping/Atmos_02_sim_05.csv', header=0)
data_02_sim_06 = pd.read_csv('./NASA_data/trajectory_data/Atmos_02_TumblingBrickNoDamping/Atmos_02_sim_06.csv', header=0)
data_02_pysim  = np.load('./my_saved_data/atmos02/atmos02.npy')

#print(data_02_sim_01.columns)
#print(data_02_sim_02.columns)
#print(data_02_sim_04.columns)
#print(data_02_sim_05.columns)
#print(data_02_sim_06.columns)

# A. Create figure of translational and rotational 6-DOF states
fig, axes = plt.subplots(2, 3, figsize=(10, 6))
fig.set_facecolor('black')  
title_suffix = 'NASA Verification Case 02: Tumbling Brick Without Drag or Damping'
title_prefix = "Six Degree of Freedom Translational and Rotational Equations\n"
fig.suptitle(title_prefix + title_suffix, fontsize=14, fontweight='normal', color='white')

# For the velocity plots we need to transform the velocity in body to velocity in NED
t_s = data_02_pysim[:,0]; nt_s = t_s.size
u_n_mps  = np.zeros((nt_s,1))
v_n_mps  = np.zeros((nt_s,1))
w_n_mps  = np.zeros((nt_s,1))
u_b_mps = data_02_pysim[:,1]
v_b_mps = data_02_pysim[:,2]
w_b_mps = data_02_pysim[:,3]
C_phi   = data_02_pysim[:,20]
C_theta = data_02_pysim[:,21]
C_psi   = data_02_pysim[:,22]
S_phi   = data_02_pysim[:,23]
S_theta = data_02_pysim[:,24]
S_psi   = data_02_pysim[:,25]

# The is the DCM of body to NED applied to the velocity in body. These equations were taken directly from the position equations in flat_earth_eom.py
for i, element in enumerate(t_s):
    C_b2n_11 =  C_theta[i]*C_psi[i]
    C_b2n_12 = -C_phi[i]*S_psi[i] + S_phi[i]*S_theta[i]*C_psi[i]
    C_b2n_13 =  S_phi[i]*S_psi[i] + C_phi[i]*S_theta[i]*C_psi[i]
    C_b2n_21 =  C_theta[i]*S_psi[i]
    C_b2n_22 =  C_phi[i]*C_psi[i] + S_phi[i]*S_theta[i]*S_psi[i]
    C_b2n_23 = -S_phi[i]*C_psi[i] + C_phi[i]*S_theta[i]*S_psi[i]
    C_b2n_31 = -S_theta[i]
    C_b2n_32 =  S_phi[i]*C_theta[i]
    C_b2n_33 =  C_phi[i]*C_theta[i]
    u_n_mps[i] = C_b2n_11*u_b_mps[i] + C_b2n_12*v_b_mps[i] + C_b2n_13*w_b_mps[i]
    v_n_mps[i] = C_b2n_21*u_b_mps[i] + C_b2n_22*v_b_mps[i] + C_b2n_23*w_b_mps[i]
    w_n_mps[i] = C_b2n_31*u_b_mps[i] + C_b2n_32*v_b_mps[i] + C_b2n_33*w_b_mps[i]

axes[0, 0].plot(data_02_sim_01['time'], data_02_sim_01['feVelocity_ft_s_X'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[0, 0].plot(data_02_sim_02['time'], data_02_sim_02['feVelocity_ft_s_X'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[0, 0].plot(data_02_sim_04['time'], data_02_sim_04['feVelocity_ft_s_X'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[0, 0].plot(data_02_sim_05['time'], data_02_sim_05['feVelocity_ft_s_X'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[0, 0].plot(data_02_sim_06['time'], data_02_sim_06['feVelocity_ft_s_X'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[0, 0].plot(data_02_pysim[:,0],                            u_n_mps*3.28, linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0, 0].set_xlabel('Time [s]', color='white')
axes[0, 0].set_ylabel('u Velocity [ft/s]', color='white')
axes[0, 0].set_facecolor('black')
axes[0, 0].grid(True)  
axes[0, 0].tick_params(colors = 'white')
axes[0, 0].legend(loc='lower left', fontsize=6, facecolor='none', labelcolor='white', frameon=False)

# Plot v velocity versus time
axes[0, 1].plot(data_02_sim_01['time'], data_02_sim_01['feVelocity_ft_s_Y'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[0, 1].plot(data_02_sim_02['time'], data_02_sim_02['feVelocity_ft_s_Y'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[0, 1].plot(data_02_sim_04['time'], data_02_sim_04['feVelocity_ft_s_Y'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[0, 1].plot(data_02_sim_05['time'], data_02_sim_05['feVelocity_ft_s_Y'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[0, 1].plot(data_02_sim_06['time'], data_02_sim_06['feVelocity_ft_s_Y'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[0, 1].plot(data_02_pysim[:,0],                            v_n_mps*3.28, linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0, 1].set_xlabel('Time [s]', color='white')
axes[0, 1].set_ylabel('v Velocity [ft/s]', color='white')
axes[0, 1].set_facecolor('black')
axes[0, 1].grid(True)  
axes[0, 1].tick_params(colors = 'white')

# Plot v velocity versus time
axes[0, 2].plot(data_02_sim_01['time'], data_02_sim_01['feVelocity_ft_s_Z'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[0, 2].plot(data_02_sim_02['time'], data_02_sim_02['feVelocity_ft_s_Z'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[0, 2].plot(data_02_sim_04['time'], data_02_sim_04['feVelocity_ft_s_Z'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[0, 2].plot(data_02_sim_05['time'], data_02_sim_05['feVelocity_ft_s_Z'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[0, 2].plot(data_02_sim_06['time'], data_02_sim_06['feVelocity_ft_s_Z'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[0, 2].plot(data_02_pysim[:,0],                            w_n_mps*3.28, linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0, 2].set_xlabel('Time [s]', color='white')
axes[0, 2].set_ylabel('w Velocity [ft/s]', color='white')
axes[0, 2].set_facecolor('black')
axes[0, 2].grid(True)  
axes[0, 2].tick_params(colors = 'white')

# Plot roll rate versus time
axes[1, 0].plot(data_02_sim_01['time'], data_02_sim_01['bodyAngularRateWrtEi_deg_s_Roll'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[1, 0].plot(data_02_sim_02['time'], data_02_sim_02['bodyAngularRateWrtEi_deg_s_Roll'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[1, 0].plot(data_02_sim_04['time'], data_02_sim_04['bodyAngularRateWrtEi_deg_s_Roll'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[1, 0].plot(data_02_sim_05['time'], data_02_sim_05['bodyAngularRateWrtEi_deg_s_Roll'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[1, 0].plot(data_02_sim_06['time'], data_02_sim_06['bodyAngularRateWrtEi_deg_s_Roll'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[1, 0].plot(    data_02_pysim[:,0],                    180/math.pi*data_02_pysim[:,4], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[1, 0].set_xlabel('Time [s]', color='white')
axes[1, 0].set_ylabel('Roll Rate [deg/s]', color='white')
axes[1, 0].set_facecolor('black')
axes[1, 0].grid(True)  
axes[1, 0].tick_params(colors = 'white')

# Plot pitch rate versus time
axes[1, 1].plot(data_02_sim_01['time'], data_02_sim_01['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[1, 1].plot(data_02_sim_02['time'], data_02_sim_02['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[1, 1].plot(data_02_sim_04['time'], data_02_sim_04['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[1, 1].plot(data_02_sim_05['time'], data_02_sim_05['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[1, 1].plot(data_02_sim_06['time'], data_02_sim_06['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[1, 1].plot(    data_02_pysim[:,0],                     180/math.pi*data_02_pysim[:,5], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[1, 1].set_xlabel('Time [s]', color='white')
axes[1, 1].set_ylabel('Pitch Rate [deg/s]', color='white')
axes[1, 1].set_facecolor('black')
axes[1, 1].grid(True)  
axes[1, 1].tick_params(colors = 'white')

# Plot yaw rate versus time
axes[1, 2].plot(data_02_sim_01['time'], data_02_sim_01['bodyAngularRateWrtEi_deg_s_Yaw'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[1, 2].plot(data_02_sim_02['time'], data_02_sim_02['bodyAngularRateWrtEi_deg_s_Yaw'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[1, 2].plot(data_02_sim_04['time'], data_02_sim_04['bodyAngularRateWrtEi_deg_s_Yaw'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[1, 2].plot(data_02_sim_05['time'], data_02_sim_05['bodyAngularRateWrtEi_deg_s_Yaw'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[1, 2].plot(data_02_sim_06['time'], data_02_sim_06['bodyAngularRateWrtEi_deg_s_Yaw'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[1, 2].plot(    data_02_pysim[:,0],                   180/math.pi*data_02_pysim[:,6], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[1, 2].set_xlabel('Time [s]', color='white')
axes[1, 2].set_ylabel('Yaw Rate [deg/s]', color='white')
axes[1, 2].set_facecolor('black')
axes[1, 2].grid(True)  
axes[1, 2].tick_params(colors = 'white')
plt.tight_layout()

# B. Create figure of Euler angles
fig, axes = plt.subplots(1, 3, figsize=(10, 3))
fig.set_facecolor('black')  
title_suffix = 'NASA Verification Case 02: Tumbling Brick Without Drag or Damping'
title_prefix = "Euler Angles\n"
fig.suptitle(title_prefix + title_suffix, fontsize=14, fontweight='normal', color='white')

axes[0].plot(data_02_sim_01['time'], data_02_sim_01['eulerAngle_deg_Roll'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[0].plot(data_02_sim_02['time'], data_02_sim_02['eulerAngle_deg_Roll'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[0].plot(data_02_sim_04['time'], data_02_sim_04['eulerAngle_deg_Roll'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[0].plot(data_02_sim_05['time'], data_02_sim_05['eulerAngle_deg_Roll'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[0].plot(data_02_sim_06['time'], data_02_sim_06['eulerAngle_deg_Roll'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[0].plot(data_02_pysim[:,0],            180/math.pi*data_02_pysim[:,7], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0].set_xlabel('Time [s]', color='white')
axes[0].set_ylabel('Roll Angle [deg]', color='white')
axes[0].set_facecolor('black')
axes[0].grid(True)  
axes[0].tick_params(colors = 'white')
axes[0].legend(loc='lower left', fontsize=6, facecolor='none', labelcolor='white', frameon=False)

axes[1].plot(data_02_sim_01['time'], data_02_sim_01['eulerAngle_deg_Pitch'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[1].plot(data_02_sim_02['time'], data_02_sim_02['eulerAngle_deg_Pitch'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[1].plot(data_02_sim_04['time'], data_02_sim_04['eulerAngle_deg_Pitch'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[1].plot(data_02_sim_05['time'], data_02_sim_05['eulerAngle_deg_Pitch'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[1].plot(data_02_sim_06['time'], data_02_sim_06['eulerAngle_deg_Pitch'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[1].plot(    data_02_pysim[:,0],         180/math.pi*data_02_pysim[:,8], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[1].set_xlabel('Time [s]', color='white')
axes[1].set_ylabel('Pitch Angle [deg]', color='white')
axes[1].set_facecolor('black')
axes[1].grid(True)  
axes[1].tick_params(colors = 'white')

axes[2].plot(data_02_sim_01['time'], data_02_sim_01['eulerAngle_deg_Yaw'], linewidth=1, color='blue',    label="Atmos_02_sim_01")  
axes[2].plot(data_02_sim_02['time'], data_02_sim_02['eulerAngle_deg_Yaw'], linewidth=1, color='green',   label="Atmos_02_sim_02")  
axes[2].plot(data_02_sim_04['time'], data_02_sim_04['eulerAngle_deg_Yaw'], linewidth=1, color='red',     label="Atmos_02_sim_04")  
axes[2].plot(data_02_sim_05['time'], data_02_sim_05['eulerAngle_deg_Yaw'], linewidth=1, color='white',   label="Atmos_02_sim_05")  
axes[2].plot(data_02_sim_06['time'], data_02_sim_06['eulerAngle_deg_Yaw'], linewidth=1, color='gray',    label="Atmos_02_sim_06")  
axes[2].plot(    data_02_pysim[:,0],       180/math.pi*data_02_pysim[:,9], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[2].set_xlabel('Time [s]', color='white')
axes[2].set_ylabel('Yaw Angle [deg]', color='white')
axes[2].set_facecolor('black')
axes[2].grid(True)  
axes[2].tick_params(colors = 'white')
plt.tight_layout()
plt.show()