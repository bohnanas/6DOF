import math
import numpy as np

def sphere():
    
    # import atmospheric data
    rho_kgpm3 = 1.225

    # Define vehicle (dragless sphere for now)
    m_sphere_kg = 14.5939
    Ixx = 4.8809446628
    Iyy = 4.8809446628
    Izz = 4.8809446628
    CD = 0.1
    Aref_m2 = 0.0182414654525
    G = 9.81

    Vterm_mps = math.sqrt((2 * m_sphere_kg * G) / (rho_kgpm3 * CD * Aref_m2))

    vmod = {
        "m_kg" : m_sphere_kg, \
        "CD"   : CD, \
        "Aref_m2" : Aref_m2, \
        "Ixx_b_kgm2" : Ixx, \
        "Iyy_b_kgm2" : Iyy, \
        "Izz_b_kgm2" : Izz, \
        "Ixz_b_kgm2" : 0, \
        "Ixy_b_kgm2" : 0, \
        "Iyz_b_kgm2" : 0,
        "Vterm_mps"  : Vterm_mps }

    return vmod