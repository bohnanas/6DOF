def f16():
    
    # Define vehicle (dragless sphere for now)
    m_kg = 9300.11063
    Ixx = 12874.847366
    Iyy = 75673.623725
    Izz = 85552.113395
    Izx = 1331.4132386
    Aref_m2 = 27.8709
    b = 9.144
    cbar = 3.450336
    CD = 
    C_l_p = 
    C_m_q = 
    C_n_r = 

    vmod = {
        "m_kg" : m_kg, \
        "Ixx_b_kgm2" : Ixx, \
        "Iyy_b_kgm2" : Iyy, \
        "Izz_b_kgm2" : Izz, \
        "Ixz_b_kgm2" : Izx, \
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