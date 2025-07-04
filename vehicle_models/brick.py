def brick():
    
    # Define vehicle (dragless sphere for now)
    m_kg = 2.2679619056149
    Ixx = 0.0025682175
    Iyy = 0.0084210111
    Izz = 0.009754656
    Aref_m2 = 0.02064491355
    b = 0.101598984
    cbar = 0.203201016
    CD = 0.01
    C_l_p = -1.0
    C_m_q = -1.0
    C_n_r = -1.0

    vmod = {
        "m_kg" : m_kg, \
        "Ixx_b_kgm2" : Ixx, \
        "Iyy_b_kgm2" : Iyy, \
        "Izz_b_kgm2" : Izz, \
        "Ixz_b_kgm2" : 0, \
        "Ixy_b_kgm2" : 0, \
        "Iyz_b_kgm2" : 0,
        "Aref_m2" : Aref_m2, \
        "b"    : b, \
        "cbar" : cbar, \
        "CD"   : CD, \
        "C_l_p" : C_l_p, \
        "C_m_q" : C_m_q, \
        "C_n_r" : C_n_r, }

    return vmod