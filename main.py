import matplotlib.pyplot as plt
import numpy as np
import math

from numerical_integrators import RK4
from governing_equations import flat_earth_eom

#=====================================================================================================
# INITIALIZE SIMULATION
#=====================================================================================================

# Define vehicle (sphere for now)
r_sphere_m = 0.08
m_sphere_kg = 5
I_sphere_kgm2 = 0.4 * m_sphere_kg * r_sphere_m**2

vmod = {
    "m_kg" : m_sphere_kg, \
    "Ixx_b_kgm2" : I_sphere_kgm2, \
    "Iyy_b_kgm2" : I_sphere_kgm2, \
    "Izz_b_kgm2" : I_sphere_kgm2, \
    "Ixz_b_kgm2" : 0, \
    "Ixy_b_kgm2" : 0, \
    "Iyz_b_kgm2" : 0 }

# Set initial state
u0_b_mps = 0
v0_b_mps = 0
w0_b_mps = 0
p0_b_mps = 0
q0_b_mps = 0
r0_b_mps = 0
phi0_rad = 0
theta0_rad = 90 * math.pi / 180
psi0_rad = 0
p10_NED_m = 0
p20_NED_m = 0
p30_NED_m = 0

# Put initial state into an array
x0 = np.array([
    u0_b_mps,   # x-axis body frame velocity (m/s)
    v0_b_mps,   # y-axis body frame velocity (m/s)
    w0_b_mps,   # z-axis body frame velocity (m/s)
    p0_b_mps,   # roll rate (rad/s)
    q0_b_mps,   # pitch rate (rad/s)
    r0_b_mps,   # yaw rate (rad/s)
    phi0_rad,   # roll angle (rad)
    theta0_rad, # pitch angle (rad)
    psi0_rad,   # yaw angle (rad)
    p10_NED_m,  # x-axis position in inertial frame (NED frame)
    p20_NED_m,  # y-axis position in inertial frame (NED frame)
    p30_NED_m,  # z-axis position in inertial frame (NED frame)
])

# Get number of elements in x0
nx0 = x0.size

# Set time conditions for fourth-order Runge Kutta
t0_s = 0.0  # initial time
tf_s = 10.0 # final time (simulation length)
h_s = 0.01  # time step

# Declare solution array
t_s = np.arange(t0_s, tf_s + h_s, h_s)      # array of time steps from initial to final, incremented by h
nt_s = t_s.size                             # number of time steps
x = np.empty((nx0, nt_s), dtype = float)    # solution array where state vector is stored

# Assign initial condition x0 the solution array x
x[:,0] = x0 # sets the first column of the solution array, x, to the initial state, x0, at time t_s[0]

#=====================================================================================================
# SIMULATE NUMERICALLY WITH FOURTH ORDER RUNGE KUTTA INTEGRATION
#=====================================================================================================

# Calls on RK4 function in RK4 file to do the integration over our defined time conditions
# Sends: governing equation function, initial state vector, initial time, final time, step size, and vehicle model
# Returns: array of time steps and the state at each time step
t_s, x = RK4.rk4(flat_earth_eom.flat_earth_eom, x0, t0_s, tf_s, h_s, vmod)

#=====================================================================================================
# PLOT DATA
#=====================================================================================================

x = np.array(x)
x = x.T

# Create subplots and layout
fig, axes = plt.subplots(3, 3, figsize=(10,6)) # 3 rows, 3 columns
fig.set_facecolor('black')

# Plot body x velocity, u
axes[0, 0].plot(t_s, x[0,:], color='yellow')
axes[0, 0].set_xlabel('Time (s)', color='white')
axes[0, 0].set_ylabel('u (m/s)', color='white')
axes[0, 0].set_facecolor('black')
axes[0, 0].grid(True)
axes[0, 0].tick_params(colors='white')

# Plot body y velocity, v
axes[0, 1].plot(t_s, x[1,:], color='yellow')
axes[0, 1].set_xlabel('Time (s)', color='white')
axes[0, 1].set_ylabel('v (m/s)', color='white')
axes[0, 1].set_facecolor('black')
axes[0, 1].grid(True)
axes[0, 1].tick_params(colors='white')

# Plot body z velocity, w
axes[0, 2].plot(t_s, x[2,:], color='yellow')
axes[0, 2].set_xlabel('Time (s)', color='white')
axes[0, 2].set_ylabel('w (m/s)', color='white')
axes[0, 2].set_facecolor('black')
axes[0, 2].grid(True)
axes[0, 2].tick_params(colors='white')

# Plot roll rate, p
axes[1, 0].plot(t_s, x[3,:], color='yellow')
axes[1, 0].set_xlabel('Time (s)', color='white')
axes[1, 0].set_ylabel('Roll rate (rad/s)', color='white')
axes[1, 0].set_facecolor('black')
axes[1, 0].grid(True)
axes[1, 0].tick_params(colors='white')

# Plot pitch rate, q
axes[1, 1].plot(t_s, x[4,:], color='yellow')
axes[1, 1].set_xlabel('Time (s)', color='white')
axes[1, 1].set_ylabel('Pitch rate (rad/s)', color='white')
axes[1, 1].set_facecolor('black')
axes[1, 1].grid(True)
axes[1, 1].tick_params(colors='white')

# Plot yaw rate, r
axes[1, 2].plot(t_s, x[5,:], color='yellow')
axes[1, 2].set_xlabel('Time (s)', color='white')
axes[1, 2].set_ylabel('Yaw rate (rad/s)', color='white')
axes[1, 2].set_facecolor('black')
axes[1, 2].grid(True)
axes[1, 2].tick_params(colors='white')

# Plot roll angle, phi
axes[2, 0].plot(t_s, x[6,:], color='yellow')
axes[2, 0].set_xlabel('Time (s)', color='white')
axes[2, 0].set_ylabel('Roll angle (rad)', color='white')
axes[2, 0].set_facecolor('black')
axes[2, 0].grid(True)
axes[2, 0].tick_params(colors='white')

# Plot pitch angle, theta
axes[2, 1].plot(t_s, x[7,:], color='yellow')
axes[2, 1].set_xlabel('Time (s)', color='white')
axes[2, 1].set_ylabel('Pitch angle (rad)', color='white')
axes[2, 1].set_facecolor('black')
axes[2, 1].grid(True)
axes[2, 1].tick_params(colors='white')

# Plot yaw angle, psi
axes[2, 2].plot(t_s, x[8,:], color='yellow')
axes[2, 2].set_xlabel('Time (s)', color='white')
axes[2, 2].set_ylabel('Yaw angle (rad)', color='white')
axes[2, 2].set_facecolor('black')
axes[2, 2].grid(True)
axes[2, 2].tick_params(colors='white')

plt.tight_layout()
plt.show()