"""
Simple Flight Dynamics Model (FDM) example of the X-15

History:
Adapted from "Simple FDM Loop" example by Ben Dickinson
9/10/2024
"""
import time
import numpy as np
from flightgear_python.fg_if import FDMConnection

def fdm_callback(fdm_data, event_pipe):
    
    if event_pipe.child_poll():
        
        # Unpack tuple
        alt_m_child, phi_rad_child, theta_rad_child, psi_rad_child, \
        u_b_ftps_child, v_b_ftps_child, w_b_ftps_child, alpha_rad_child, \
        beta_rad_child, lat_rad_child, long_rad_child, dela_rad_child, \
        dele_rad_child, delr_rad_child, v_north_ft_per_s_parent, \
        v_east_ft_per_s_parent, v_down_ft_per_s_parent, \
        = event_pipe.child_recv()  
        
        # Set only the data that we need to (we can force our own values)
        fdm_data['alt_m']        = alt_m_child  
        fdm_data['phi_rad']      = phi_rad_child
        fdm_data['theta_rad']    = theta_rad_child
        fdm_data['psi_rad']      = psi_rad_child
        fdm_data['v_body_u']     = u_b_ftps_child
        fdm_data['v_body_v']     = v_b_ftps_child
        fdm_data['v_body_w']     = w_b_ftps_child
        fdm_data['alpha_rad']    = alpha_rad_child
        fdm_data['beta_rad']     = beta_rad_child
        fdm_data['lat_rad']      = lat_rad_child
        fdm_data['lon_rad']      = long_rad_child
        fdm_data['elevator']     = dele_rad_child
        fdm_data['v_north_ft_per_s_parent'] = v_north_ft_per_s_parent
        fdm_data['v_east_ft_per_s_parent']  = v_east_ft_per_s_parent
        fdm_data['v_down_ft_per_s_parent']  = v_down_ft_per_s_parent
        
    # Return the whole structure
    return fdm_data  

"""
Start FlightGear with: 
`--native-fdm=socket,out,30,localhost,5501,udp --native-fdm=socket,in,30,localhost,5502,udp`
(you probably also want `--fdm=null` and `--max-fps=30` to stop the simulation fighting with
these external commands)
"""
if __name__ == '__main__':  # NOTE: This is REQUIRED on Windows!
    
    fdm_conn = FDMConnection()
    fdm_event_pipe = fdm_conn.connect_rx('localhost', 5501, fdm_callback)
    fdm_conn.connect_tx('localhost', 5502)
    
    # Start the FDM RX/TX loop
    fdm_conn.start()  
    
    # Get Python 6-DOF simulation data
    data_pysim  = np.load('./my_saved_data/x15_test2/x15_test2.npy')
    
    # Get time and other variables
    t_s   = data_pysim[:,0]; nt_s = t_s.size

    i = 0
    while True:
        
        # Increment time step counter
        i += 1
        
        # Get present altitude
        v_body_u_parent         = data_pysim[i,1]*3.28 # (converts m/s to ft/s)
        v_body_v_parent         = data_pysim[i,2]*3.28
        v_body_w_parent         = data_pysim[i,3]*3.28
        phi_rad_parent          = data_pysim[i,7]
        theta_rad_parent        = data_pysim[i,8]
        psi_rad_parent          = data_pysim[i,9]
        alt_m_parent            = data_pysim[i,13]
        dela_rad_parent         = data_pysim[i,33]
        dele_rad_parent         = data_pysim[i,34]
        delr_rad_parent         = data_pysim[i,35] 
        alpha_rad_parent        = data_pysim[i,17]
        beta_rad_parent         = data_pysim[i,18]
        v_north_ft_per_s_parent = data_pysim[i,30]*3.28 # (converts m/s to ft/s)
        v_east_ft_per_s_parent  = data_pysim[i,31]*3.28
        v_down_ft_per_s_parent  = data_pysim[i,32]*3.28
        lat_rad_parent          = data_pysim[i,37]*3.14/180 # (convert degrees to radians)
        long_rad_parent         = data_pysim[i,38]*3.14/180
        
        # Send tuple (could also do `fdm_conn.event_pipe.parent_send` so you just need to pass around `fdm_conn`)
        fdm_event_pipe.parent_send((alt_m_parent, phi_rad_parent, theta_rad_parent, psi_rad_parent, \
            v_body_u_parent, v_body_v_parent, v_body_w_parent, alpha_rad_parent, beta_rad_parent, \
            lat_rad_parent, long_rad_parent, dela_rad_parent, dele_rad_parent, delr_rad_parent, \
            v_north_ft_per_s_parent, v_east_ft_per_s_parent, v_down_ft_per_s_parent))  
        
        if i == nt_s-1:
            print('\nVisualization completed.\n')
            break
        
        time.sleep(0.007)