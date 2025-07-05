import math
import numpy as np

def NASA_sphere():
    # unit conversions
    slug_ft2_2_kg_m2 = 1.3558179619
    slug2kg = 14.594
    ft2_2_m2 = 0.09290304

    # Define vehicle (dragless sphere for now)
    m_kg = 1 * slug2kg
    I = 3.6 * slug_ft2_2_kg_m2
    CD = 0.1
    Aref_m2 = 0.1963495 * ft2_2_m2

    vmod = {
        "V_Name": "NASA Sphere", \
        "m_kg" : m_kg, \
        "CD"   : CD, \
        "Aref_m2" : Aref_m2, \
        "Ixx_b_kgm2" : I, \
        "Iyy_b_kgm2" : I, \
        "Izz_b_kgm2" : I, \
        "Ixz_b_kgm2" : 0, \
        "Ixy_b_kgm2" : 0, \
        "Iyz_b_kgm2" : 0 }

    return vmod

def musket_ball():

    # set atmospheric data
    rho_kgpm3 = 1.2

    # inches to meters conversion
    conv = 0.0254

    # sphere density
    rho_lead_kgpm3 = 11300

    # sphere radius
    r_m = 0.495 * conv

    # Define vehicle
    m_kg, I, Aref_m2 = sphereCalc(r_m, rho_lead_kgpm3)
    CD = 0.5

    # Analytical terminal velocity
    G = 9.81

    Vterm_mps = math.sqrt((2 * m_kg * G) / (rho_kgpm3 * CD * Aref_m2))

    # return the model in a dictionary
    vmod = {
        "V_Name": 'Musket Ball', \
        "m_kg" : m_kg, \
        "CD"   : CD, \
        "Aref_m2" : Aref_m2, \
        "Ixx_b_kgm2" : I, \
        "Iyy_b_kgm2" : I, \
        "Izz_b_kgm2" : I, \
        "Ixz_b_kgm2" : 0, \
        "Ixy_b_kgm2" : 0, \
        "Iyz_b_kgm2" : 0,
        "Vterm_mps"  : Vterm_mps }

    return vmod

def bowling_ball():

    # set atmospheric data
    rho_kgpm3 = 1.2

    # inches to meters conversion
    conv = 0.0254

    # sphere density
    rho_bowling_kgpm3 = 1500

    # sphere radius
    r_m = 4.4 * conv

    # Define vehicle
    m_kg, I, Aref_m2 = sphereCalc(r_m, rho_bowling_kgpm3)
    CD = 0.5

    # Analytical terminal velocity
    G = 9.81

    Vterm_mps = math.sqrt((2 * m_kg * G) / (rho_kgpm3 * CD * Aref_m2))

    # return the model in a dictionary
    vmod = {
        "V_Name": 'Bowling Ball', \
        "m_kg" : m_kg, \
        "CD"   : CD, \
        "Aref_m2" : Aref_m2, \
        "Ixx_b_kgm2" : I, \
        "Iyy_b_kgm2" : I, \
        "Izz_b_kgm2" : I, \
        "Ixz_b_kgm2" : 0, \
        "Ixy_b_kgm2" : 0, \
        "Iyz_b_kgm2" : 0,
        "Vterm_mps"  : Vterm_mps }

    return vmod

def blueberry():

    # set atmospheric data
    rho_kgpm3 = 1.2

    # inches to meters conversion
    conv = 0.0254

    # sphere density
    rho_blueberry_kgpm3 = 786

    # sphere radius
    r_m = 0.3 * conv

    # Define vehicle
    m_kg, I, Aref_m2 = sphereCalc(r_m, rho_blueberry_kgpm3)
    CD = 0.5

    # Analytical terminal velocity
    G = 9.81

    Vterm_mps = math.sqrt((2 * m_kg * G) / (rho_kgpm3 * CD * Aref_m2))

    # return the model in a dictionary
    vmod = {
        "V_Name": 'Blueberry', \
        "m_kg" : m_kg, \
        "CD"   : CD, \
        "Aref_m2" : Aref_m2, \
        "Ixx_b_kgm2" : I, \
        "Iyy_b_kgm2" : I, \
        "Izz_b_kgm2" : I, \
        "Ixz_b_kgm2" : 0, \
        "Ixy_b_kgm2" : 0, \
        "Iyz_b_kgm2" : 0,
        "Vterm_mps"  : Vterm_mps }

    return vmod

def sphereCalc(r_m, rho_kgpm3):
    # calculates volume, mass, moment of inertia, and reference area of a sphere
    vol_m3 = 4/3*math.pi*r_m**3
    m_kg = rho_kgpm3*vol_m3
    MOI = 0.4*m_kg*r_m**2
    Aref_m2 = math.pi*r_m**2

    return m_kg, MOI, Aref_m2