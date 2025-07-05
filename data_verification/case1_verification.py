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

# x-combined from Python simulation column order
# 0: t_s 
# 1: u_b_mps, axial velocity of CM wrt inertial CS resolved in aircraft body fixed CS
# 2: v_b_mps, lateral velocity of CM wrt inertial CS resolved in aircraft body fixed CS
# 3: w_b_mps, vertical velocity of CM wrt inertial CS resolved in aircraft body fixed CS
# 4: p_b_rps, roll angular velocity of body fixed CS with respect to inertial CS
# 5: q_b_rps, pitch angular velocity of body fixed CS with respect to inertial CS
# 6: r_b_rps, yaw angular velocity of body fixed CS with respect to inertial CS
# 7: phi_rad, roll angle
# 8: theta_rad, pitch angle
# 9: psi_rad, yaw angle
# 10: p1_n_m, x-axis position of aircraft resolved in NED CS
# 11: p2_n_m, y-axis position of aircraft resolved in NED CS
# 12: p3_n_m, z-axis position of aircraft resolved in NED CS
# 13: Altitude_m, altitude (= -p3_n_m)
# 14: Cs_mps, speed of sound interpolated from atmosphere model
# 15: Rho_kgpm3, Air density interpolated from atmosphere model
# 16: Mach, Mach number
# 17: Alpha_rad, Angle of attack
# 18: Beta_rad, Angle of sideslip
# 19: True_Airspeed_mps, True airspeed
#

"""

# Atmospheric case 01: Dropped Sphere, dragless (WGS-84 earth)
# Read the CSV file, specifying the header row
data_01_sim_01 = pd.read_csv('./NASA_data/trajectory_data/Atmos_01_DroppedSphere/Atmos_01_sim_01.csv', header=0)
data_01_sim_02 = pd.read_csv('./NASA_data/trajectory_data/Atmos_01_DroppedSphere/Atmos_01_sim_02.csv', header=0)
data_01_sim_03 = pd.read_csv('./NASA_data/trajectory_data/Atmos_01_DroppedSphere/Atmos_01_sim_03.csv', header=0)
data_01_sim_04 = pd.read_csv('./NASA_data/trajectory_data/Atmos_01_DroppedSphere/Atmos_01_sim_04.csv', header=0)
data_01_sim_05 = pd.read_csv('./NASA_data/trajectory_data/Atmos_01_DroppedSphere/Atmos_01_sim_05.csv', header=0)
data_01_sim_06 = pd.read_csv('./NASA_data/trajectory_data/Atmos_01_DroppedSphere/Atmos_01_sim_06.csv', header=0)
data_01_pysim  = np.load('./my_saved_data/atmos01/atmos01.npy')

print(data_01_sim_01.columns)
print(data_01_sim_02.columns)
print(data_01_sim_03.columns)
print(data_01_sim_04.columns)
print(data_01_sim_05.columns)
print(data_01_sim_06.columns)

# A. Create figure of translational and rotation states
fig, axes = plt.subplots(2, 4, figsize=(10, 6))
fig.set_facecolor('black')  
fig.suptitle('NASA Verification Case 01: Dragless Sphere', fontsize=14, fontweight='bold', color='white')

# Plot altitude versus time
axes[0, 0].plot(data_01_sim_01['time'], data_01_sim_01['altitudeMsl_ft'], linewidth=1, color='blue',    label="Atmos_01_sim_01")  
axes[0, 0].plot(data_01_sim_02['time'], data_01_sim_02['altitudeMsl_ft'], linewidth=1, color='green',   label="Atmos_01_sim_02")  
axes[0, 0].plot(data_01_sim_03['time'], data_01_sim_03['altitudeMsl_ft'], linewidth=1, color='orange',  label="Atmos_01_sim_03")  
axes[0, 0].plot(data_01_sim_04['time'], data_01_sim_04['altitudeMsl_ft'], linewidth=1, color='red',     label="Atmos_01_sim_04")  
axes[0, 0].plot(data_01_sim_05['time'], data_01_sim_05['altitudeMsl_ft'], linewidth=1, color='white',   label="Atmos_01_sim_05")  
axes[0, 0].plot(data_01_sim_06['time'], data_01_sim_06['altitudeMsl_ft'], linewidth=1, color='gray',    label="Atmos_01_sim_06")  
axes[0, 0].plot(data_01_pysim[:,0],             data_01_pysim[:,13]*3.28, linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0, 0].set_xlabel('Time [s]', color='white')
axes[0, 0].set_ylabel('Altitude Mean Sea Level [ft]', color='white')
axes[0, 0].set_facecolor('black')
axes[0, 0].grid(True)  
axes[0, 0].tick_params(colors = 'white')
axes[0, 0].legend(loc='lower left', fontsize=6, facecolor='none', labelcolor='white', frameon=False)

# Plot Mach versus time
axes[0, 1].plot(data_01_sim_01['time'], data_01_sim_01['mach'], linewidth=1, color='blue',    label="Atmos_01_sim_01")  
axes[0, 1].plot(data_01_sim_02['time'], data_01_sim_02['mach'], linewidth=1, color='green',   label="Atmos_01_sim_02")  
#axes[0, 1].plot(data_01_sim_03['time'], data_01_sim_03['mach'], linewidth=1, color='orange',  label="Atmos_01_sim_03")  
axes[0, 1].plot(data_01_sim_04['time'], data_01_sim_04['mach'], linewidth=1, color='red',     label="Atmos_01_sim_04")  
axes[0, 1].plot(data_01_sim_05['time'], data_01_sim_05['mach'], linewidth=1, color='white',   label="Atmos_01_sim_05")  
axes[0, 1].plot(data_01_sim_06['time'], data_01_sim_06['mach'], linewidth=1, color='gray',    label="Atmos_01_sim_06")  
axes[0, 1].plot(data_01_pysim[:,0],        data_01_pysim[:,16], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0, 1].set_xlabel('Time [s]', color='white')
axes[0, 1].set_ylabel('Mach', color='white')
axes[0, 1].set_facecolor('black')
axes[0, 1].grid(True)  
axes[0, 1].tick_params(colors = 'white')

# Plot air density versus time
axes[0, 2].plot(data_01_sim_01['time'], data_01_sim_01['airDensity_slug_ft3'], linewidth=1, color='blue',    label="Atmos_01_sim_01")  
axes[0, 2].plot(data_01_sim_02['time'], data_01_sim_02['airDensity_slug_ft3'], linewidth=1, color='green',   label="Atmos_01_sim_02")  
axes[0, 2].plot(data_01_sim_03['time'], data_01_sim_03['airDensity_slug_ft3'], linewidth=1, color='orange',  label="Atmos_01_sim_03")  
axes[0, 2].plot(data_01_sim_04['time'], data_01_sim_04['airDensity_slug_ft3'], linewidth=1, color='red',     label="Atmos_01_sim_04")  
axes[0, 2].plot(data_01_sim_05['time'], data_01_sim_05['airDensity_slug_ft3'], linewidth=1, color='white',   label="Atmos_01_sim_05")  
axes[0, 2].plot(data_01_sim_06['time'], data_01_sim_06['airDensity_slug_ft3'], linewidth=1, color='gray',    label="Atmos_01_sim_06")  
axes[0, 2].plot(data_01_pysim[:,0],           0.001941811*data_01_pysim[:,15], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0, 2].set_xlabel('Time [s]', color='white')
axes[0, 2].set_ylabel('Air Density [slug/ft^3]', color='white')
axes[0, 2].set_facecolor('black')
axes[0, 2].grid(True)  
axes[0, 2].tick_params(colors = 'white')

# Plot Euler pitch angle versus time
axes[0, 3].plot(data_01_sim_01['time'], data_01_sim_01['eulerAngle_deg_Pitch'], linewidth=1, color='blue',    label="Atmos_01_sim_01")  
axes[0, 3].plot(data_01_sim_02['time'], data_01_sim_02['eulerAngle_deg_Pitch'], linewidth=1, color='green',   label="Atmos_01_sim_02")  
#axes[0, 3].plot(data_01_sim_03['time'], data_01_sim_03['eulerAngle_deg_Pitch'], linewidth=1, color='orange',  label="Atmos_01_sim_03")  
axes[0, 3].plot(data_01_sim_04['time'], data_01_sim_04['eulerAngle_deg_Pitch'], linewidth=1, color='red',     label="Atmos_01_sim_04")  
axes[0, 3].plot(data_01_sim_05['time'], data_01_sim_05['eulerAngle_deg_Pitch'], linewidth=1, color='white',   label="Atmos_01_sim_05")  
axes[0, 3].plot(data_01_sim_06['time'], data_01_sim_06['eulerAngle_deg_Pitch'], linewidth=1, color='gray',    label="Atmos_01_sim_06")  
axes[0, 3].plot(data_01_pysim[:,0],           180/math.pi*data_01_pysim[:,8], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[0, 3].set_xlabel('Time [s]', color='white')
axes[0, 3].set_ylabel('Pitch Angle [deg]', color='white')
axes[0, 3].set_facecolor('black')
axes[0, 3].grid(True)  
axes[0, 3].tick_params(colors = 'white')

# Plot speed of sound versus time
axes[1, 0].plot(data_01_sim_01['time'], data_01_sim_01['speedOfSound_ft_s'], linewidth=1, color='blue',    label="Atmos_01_sim_01")  
axes[1, 0].plot(data_01_sim_02['time'], data_01_sim_02['speedOfSound_ft_s'], linewidth=1, color='green',   label="Atmos_01_sim_02")  
axes[1, 0].plot(data_01_sim_03['time'], data_01_sim_03['speedOfSound_ft_s'], linewidth=1, color='orange',  label="Atmos_01_sim_03")  
axes[1, 0].plot(data_01_sim_04['time'], data_01_sim_04['speedOfSound_ft_s'], linewidth=1, color='red',     label="Atmos_01_sim_04")  
axes[1, 0].plot(data_01_sim_05['time'], data_01_sim_05['speedOfSound_ft_s'], linewidth=1, color='white',   label="Atmos_01_sim_05")  
axes[1, 0].plot(data_01_sim_06['time'], data_01_sim_06['speedOfSound_ft_s'], linewidth=1, color='gray',    label="Atmos_01_sim_06")  
axes[1, 0].plot(data_01_pysim[:,0],                3.28*data_01_pysim[:,14], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[1, 0].set_xlabel('Time [s]', color='white')
axes[1, 0].set_ylabel('Speed of Sound [ft/s]', color='white')
axes[1, 0].set_facecolor('black')
axes[1, 0].grid(True)  
axes[1, 0].tick_params(colors = 'white')

# Plot true airspeed versus time
axes[1, 1].plot(data_01_sim_01['time'], data_01_sim_01['trueAirspeed_nmi_h'], linewidth=1, color='blue',    label="Atmos_01_sim_01")  
axes[1, 1].plot(data_01_sim_02['time'], data_01_sim_02['trueAirspeed_nmi_h'], linewidth=1, color='green',   label="Atmos_01_sim_02")  
#axes[1, 1].plot(data_01_sim_03['time'], data_01_sim_03['trueAirspeed_nmi_h'], linewidth=1, color='orange',  label="Atmos_01_sim_03")  
#axes[1, 1].plot(data_01_sim_04['time'], data_01_sim_04['trueAirspeed_nmi_h'], linewidth=1, color='red',     label="Atmos_01_sim_04")  
axes[1, 1].plot(data_01_sim_05['time'], data_01_sim_05['trueAirspeed_nmi_h'], linewidth=1, color='white',   label="Atmos_01_sim_05")  
axes[1, 1].plot(data_01_sim_06['time'], data_01_sim_06['trueAirspeed_nmi_h'], linewidth=1, color='gray',    label="Atmos_01_sim_06")  
axes[1, 1].plot(data_01_pysim[:,0],              1.94384*data_01_pysim[:,19], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[1, 1].set_xlabel('Time [s]', color='white')
axes[1, 1].set_ylabel('True Airspeed [nmi_h]', color='white')
axes[1, 1].set_facecolor('black')
axes[1, 1].grid(True)  
axes[1, 1].tick_params(colors = 'white')

# Plot pitch rate versus time
axes[1, 2].plot(data_01_sim_01['time'], data_01_sim_01['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='blue',    label="Atmos_01_sim_01")  
axes[1, 2].plot(data_01_sim_02['time'], data_01_sim_02['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='green',   label="Atmos_01_sim_02")  
#axes[1, 2].plot(data_01_sim_03['time'], data_01_sim_03['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='orange',  label="Atmos_01_sim_03")  
axes[1, 2].plot(data_01_sim_04['time'], data_01_sim_04['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='red',     label="Atmos_01_sim_04")  
axes[1, 2].plot(data_01_sim_05['time'], data_01_sim_05['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='white',   label="Atmos_01_sim_05")  
axes[1, 2].plot(data_01_sim_06['time'], data_01_sim_06['bodyAngularRateWrtEi_deg_s_Pitch'], linewidth=1, color='gray',    label="Atmos_01_sim_06")  
axes[1, 2].plot(data_01_pysim[:,0],                         180/math.pi*data_01_pysim[:,5], linewidth=1, color='magenta', label="My Sim", linestyle='dashed')  
axes[1, 2].set_xlabel('Time [s]', color='white')
axes[1, 2].set_ylabel('Pitch rate [deg/s]', color='white')
axes[1, 2].set_facecolor('black')
axes[1, 2].grid(True)  
axes[1, 2].tick_params(colors = 'white')

plt.tight_layout()
plt.show()