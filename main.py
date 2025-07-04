# Import Libraries
import matplotlib.pyplot as plt
import numpy as np
import math
import ussa1976

# Import partner files
from numerical_integrators import RK4
from governing_equations import flat_earth_eom
from vehicle_models.sphere import sphere

#=====================================================================================================
# INITIALIZE SIMULATION
#=====================================================================================================

# Get atmospheric tables
atmosphere = ussa1976.compute()

# Store atmospheric vectors of interest into a dictionary
alt_m = atmosphere["z"].values # altitude
rho_kgpm3 = atmosphere["rho"].values # density
c_mps = atmosphere["cs"].values # speed of sound
g_mps2 = ussa1976.core.compute_gravity(alt_m) 

amod = {
    "alt_m"     : alt_m,
    "rho_kgpm3" : rho_kgpm3,
    "c_mps"     : c_mps,
    "g_mps2"    : g_mps2
}

# Get vehicle model
vmod = sphere()

print(f"Analytical terminal velocity: {vmod['Vterm_mps']:.2f} m/s") # print the terminal velocity

# Set initial condition
u0_b_mps = 0.00001 # set small to avoid dividing by zero in AoA and Sideslip
v0_b_mps = 0
w0_b_mps = 0
p0_b_rps = 0
q0_b_rps = 0
r0_b_rps = 0
phi0_rad = 0 * math.pi / 180
theta0_rad = 0 * math.pi / 180
psi0_rad = 0 * math.pi / 180
p10_NED_m = 0
p20_NED_m = 0
p30_NED_m = -30000

# Put initial conditions into an array
x0 = np.array([
    u0_b_mps,   # x-axis body frame velocity (m/s)
    v0_b_mps,   # y-axis body frame velocity (m/s)
    w0_b_mps,   # z-axis body frame velocity (m/s)
    p0_b_rps,   # roll rate (rad/s)
    q0_b_rps,   # pitch rate (rad/s)
    r0_b_rps,   # yaw rate (rad/s)
    phi0_rad,   # roll angle (rad)
    theta0_rad, # pitch angle (rad)
    psi0_rad,   # yaw angle (rad)
    p10_NED_m,  # x-axis position in inertial frame (NED frame)
    p20_NED_m,  # y-axis position in inertial frame (NED frame)
    p30_NED_m,  # z-axis position in inertial frame (NED frame)
])

# Get number of elements in x0
nx0 = x0.size

#=====================================================================================================
# SIMULATE NUMERICALLY WITH FOURTH ORDER RUNGE KUTTA INTEGRATION
#=====================================================================================================

# Set time conditions for integrator
t0_s = 0.0  # initial time
tf_s = 185.0 # final time (simulation length)
h_s = 0.01  # time step

# Initialize solution matrix
t_s = np.arange(t0_s, tf_s + h_s, h_s)      # array of time steps from initial to final, incremented by h
nt_s = t_s.size                             # number of time steps
x = np.empty((nx0, nt_s), dtype = float)    # solution array where state vector is stored

# Assign initial condition x0 the solution array x
x[:,0] = x0 # sets the first column of the solution array, x, to the initial state, x0, at time t_s[0]

# Calls on RK4 function in RK4 file to do the integration over our defined time conditions
# Sends: governing equation function, initial state vector, initial time, final time, step size, and vehicle model
# Returns: array of time steps and the state at each time step
t_s, x = RK4.rk4(flat_earth_eom.flat_earth_eom, x0, t0_s, tf_s, h_s, vmod, amod)

#=====================================================================================================
# DATA POST PROCESSING
#=====================================================================================================

x = np.array(x)
x = x.T

print(f"Numerical Terminal Velocity: {x[2,-1]:.2f} m/s")

#=====================================================================================================
# PLOT DATA
#=====================================================================================================

fig, axes = plt.subplots(3, 3, figsize=(20,10))
fig.set_facecolor('black')

# Variable names and labels for the y-axis
var_labels = [
    'u (m/s)', 'v (m/s)', 'w (m/s)',
    'Roll rate (rad/s)', 'Pitch rate (rad/s)', 'Yaw rate (rad/s)',
    'Roll angle (rad)', 'Pitch angle (rad)', 'Yaw angle (rad)'
]

for i in range(9):
    ax = axes[i // 3, i % 3]
    ax.plot(t_s, x[i, :], color='yellow')
    ax.set_xlabel('Time (s)', color='white')
    ax.set_ylabel(var_labels[i], color='white')
    ax.set_facecolor('black')
    ax.grid(True)
    ax.tick_params(colors='white')

fig2, axes2 = plt.subplots(2, 3, figsize=(17,10), constrained_layout=True)
fig2.set_facecolor('black')

pos_labels = ['North (m)', 'East (m)', 'Altitude (m)']

north = x[9,:]
east = x[10,:]
altitude = -x[11,:]  # invert down to altitude

# --- TOP ROW: position vs time ---
for i in range(3):
    ax = axes2[0, i]
    if i == 2:
        ax.plot(t_s, altitude, color='yellow')
    else:
        ax.plot(t_s, x[9+i,:], color='yellow')
    ax.set_xlabel('Time (s)', color='white')
    ax.set_ylabel(pos_labels[i], color='white')
    ax.set_facecolor('black')
    ax.grid(True)
    ax.tick_params(colors='white')
    plt.setp(ax.get_xticklabels(), rotation=30, ha='right', fontsize=8)

# --- BOTTOM ROW: position vs position ---
pos_pos_data = [
    (east, north),     # North vs East
    (east, altitude),  # Altitude vs East
    (north, altitude)  # Altitude vs North
]

pos_pos_labels = [
    ('East (m)', 'North (m)'),
    ('East (m)', 'Altitude (m)'),
    ('North (m)', 'Altitude (m)')
]

for i in range(3):
    ax = axes2[1, i]
    xdata, ydata = pos_pos_data[i]
    ax.plot(xdata, ydata, color='yellow')
    ax.set_xlabel(pos_pos_labels[i][0], color='white')
    ax.set_ylabel(pos_pos_labels[i][1], color='white')
    ax.set_facecolor('black')
    ax.grid(True)
    ax.tick_params(colors='white')

plt.tight_layout()
plt.show()