import math

def X15_mach3():
        
        # Rocket motor burned out
        # Speed brake closed
        # Data taken from WT line on plot
        # Lower rudder on

        # Conversions
        ft2m = 0.3048

        # Name of vehicle 
        vehicle_name = "X15_Space_Plane_at_Mach_3"

        # Short name for file names
        short_name = "X15_at_Mach_3"
        
        # Mass
        m_kg = 11765

        # Reference wing span
        b_m = 22.36*ft2m

        # Reference mean aerodynamic chord
        c_m = 10.27*ft2m
        
        # Reference area
        A_ref_m2 = 18.6
        
        # Aerodynamic coefficients will be set to Mach 3 for an initial study
        #--------------------------------------------------------------------------------
        
        # Drag coefficients
        CDwb         = 0.062
        CDdele_pdeg  = 0
        CDdelsb_pdeg = 0
        
        # Lift coefficients
        CLwb         = 0
        CLalpha_pdeg = 0.03
        CLdele_pdeg  = 0.005
        
        # Side force coefficients
        CYbeta_prad    = -0.02
        CYp_prps       = 0.017
        CYr_prps       = -0.92
        CYbetadot_pdps = 0
        CYdela_pdeg    = -0.0003
        CYdelr_pdeg    = 0.0057
        
        # Roll moment coefficients
        Clbeta_prad    = 0
        Clp_prps       = -0.185
        Clr_prps       = -0.08
        Clbetadot_pdps = 0
        Cldela_pdeg    = 0.0007
        Cldelr_pdeg    = 0.00025
        
        # Pitch moment coefficients
        Cmwb            = 0
        Cmalpha_pdeg    = -0.01
        Cmq_prps        = -4.4
        Cmalphadot_pdps = -0.1
        Cmdele_pdeg     = -0.008 
        Cmdelsb_pdeg    = 0
        
        # Yaw moment coefficients
        Cnbeta_pdeg     = 0.0046
        Cnp_prps        = -0.012
        Cnr_prps        = -1.19
        Cnbetadot_pdps  = 0
        Cndela_pdeg     = 0.00045
        Cndelr_pdeg     = -0.0039
        
        # Moments and products of inertia fully burned out (Yancey64)
        Jyy_b_slugft2 = 86000
        Jyy_b_kgm2 = Jyy_b_slugft2*1.355
        Jzz_b_slugft2 = 88500
        Jzz_b_kgm2 = Jzz_b_slugft2*1.355
        Jxx_b_slugft2 = 3600
        Jxx_b_kgm2 = Jxx_b_slugft2*1.355
        Jxz_b_slugft2 = -700
        Jxz_b_kgm2 = Jxz_b_slugft2*1.355

        vmod = {"V_name"          : vehicle_name, 
                "m_kg"            : m_kg, 
                "Ixz_b_kgm2"      : Jxz_b_kgm2, 
                "Ixx_b_kgm2"      : Jxx_b_kgm2, 
                "Iyy_b_kgm2"      : Jyy_b_kgm2, 
                "Izz_b_kgm2"      : Jzz_b_kgm2, 
                "Ixy_b_kgm2"      : 0,
                "Iyz_b_kgm2"      : 0,
                "CD_wb"            : CDwb, 
                "CD_elev_pdeg"     : CDdele_pdeg, 
                "CD_sb_pdeg"    : CDdelsb_pdeg, 
                "CL_wb"            : CLwb, 
                "CL_alpha_pdeg"    : CLalpha_pdeg,
                "CL_elev_pdeg"     : CLdele_pdeg, 
                "CY_beta_prad"     : CYbeta_prad, 
                "CY_p_prps"        : CYp_prps, 
                "CY_r_prps"        : CYr_prps, 
                "CY_betadot_pdps"  : CYbetadot_pdps, 
                "CY_ail_pdeg"     : CYdela_pdeg, 
                "CY_rud_pdeg"     : CYdelr_pdeg, 
                "Cl_beta_prad"     : Clbeta_prad, 
                "Cl_p_prps"        : Clp_prps, 
                "Cl_r_prps"        : Clr_prps, 
                "Cl_betadot_pdps"  : Clbetadot_pdps, 
                "Cl_ail_pdeg"     : Cldela_pdeg, 
                "Cl_rud_pdeg"     : Cldelr_pdeg, 
                "Cm_wb"            : Cmwb,
                "Cm_alpha_pdeg"    : Cmalpha_pdeg, 
                "Cm_q_prps"        : Cmq_prps, 
                "Cm_alphadot_pdps" : Cmalphadot_pdps, 
                "Cm_elev_pdeg"     : Cmdele_pdeg, 
                "Cm_sb_pdeg"    : Cmdelsb_pdeg, 
                "Cn_beta_pdeg"     : Cnbeta_pdeg, 
                "Cn_p_prps"        : Cnp_prps,
                "Cn_r_prps"        : Cnr_prps, 
                "Cn_betadot_pdps"  : Cnbetadot_pdps, 
                "Cn_ail_pdeg"     : Cndela_pdeg,   
                "Cn_rud_pdeg"     : Cndelr_pdeg, 
                "Aref_m2"        : A_ref_m2, 
                "b_m"             : b_m, 
                "c_m"             : c_m, 
                "short_name"      : short_name} 

        return vmod