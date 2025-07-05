def NASA_brick():
    # unit conversions
    slug_ft2_2_kg_m2 = 1.3558179619
    slug2kg = 14.594
    ft2_2_m2 = 0.09290304
    ft2m = 0.3048

    # Define vehicle (dragless sphere for now)
    m_kg = 0.155404754 * slug2kg
    Ixx = 0.001894220 * slug_ft2_2_kg_m2
    Iyy = 0.006211019 * slug_ft2_2_kg_m2
    Izz = 0.007194665 * slug_ft2_2_kg_m2
    Aref_m2 = 0.22222 * ft2_2_m2
    b = 0.33333 * ft2m
    cbar = 0.66667 * ft2m
    CD = 0.01
    C_l_p = -1.0
    C_m_q = -1.0
    C_n_r = -1.0

    vmod = {
        "V_Name": 'NASA Brick', \
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