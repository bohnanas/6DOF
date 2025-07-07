# Drag Coefficient
def CDrag_X15(CDwb, CDdele_pdeg, CDdelsb_pdeg, dele_deg, delsb_deg):

    CD = CDwb + CDdele_pdeg*dele_deg + CDdelsb_pdeg*delsb_deg
    
    return CD

# Lift Coefficient
def CLift_X15(CLwb, CLAlpha_pdeg, CLdele_pdeg, Alpha_deg, dele_deg):

    CL = CLwb + CLAlpha_pdeg*Alpha_deg + CLdele_pdeg*dele_deg
    
    return CL

def CPitch_X15(Cmwb, Cmalpha_pdeg, Cmq_prps, Cmalphadot_pdps, Cmdele_pdeg, Cmdelsb_pdeg, alpha_deg, \
        alphadot_dps, q_rps, dele_deg, delsb_deg, true_airspeed_mps, c_m):

    Cm = Cmwb + \
        Cmalpha_pdeg*alpha_deg + \
        Cmq_prps*q_rps*c_m/(2*true_airspeed_mps) + \
        Cmalphadot_pdps*alphadot_dps*c_m/(2*true_airspeed_mps) + \
        Cmdele_pdeg*dele_deg + \
        Cmdelsb_pdeg*delsb_deg

    return Cm

def CRoll_X15(ClBeta_prad, ClP_prps, ClR_prps, ClBetadot_pdps, Cldela_pdeg, Cldelr_pdeg, Beta_rad, Betadot_dps, P_rps, r_rps, \
     dela_deg, delr_deg, true_airspeed_mps, b_m):

     Cl = ClBeta_prad*Beta_rad + ClP_prps*P_rps*b_m/(2*true_airspeed_mps) + ClR_prps*r_rps*b_m/(2*true_airspeed_mps) + \
     ClBetadot_pdps*Betadot_dps*b_m/(2*true_airspeed_mps) + Cldela_pdeg*dela_deg + Cldelr_pdeg*delr_deg

     return Cl

def CSide_X15(CYBeta_prad, CYP_prps, CYR_prps, CYBetadot_pdps, CYdela_pdeg, CYdelr_pdeg, Beta_rad, P_rps, R_rps, \
    Betadot_dps, dela_deg, delr_deg, true_airspeed_mps, b_m):

    CY = CYBeta_prad*Beta_rad + CYP_prps*P_rps*b_m/(2*true_airspeed_mps) + CYR_prps*R_rps*b_m/(2*true_airspeed_mps) + \
    CYBetadot_pdps*Betadot_dps*b_m/(2*true_airspeed_mps) + CYdela_pdeg*dela_deg + CYdelr_pdeg*delr_deg
    
    return CY

def CYaw_X15(CnBeta_pdeg, CnP_prps, CnR_prps, CnBetadot_pdps, Cndela_pdeg, Cndelr_pdeg, Beta_rad, Betadot_ps, P_rps, R_rps, \
    dela_deg, delr_deg, true_airspeed_mps, b_m):

    Cn = CnBeta_pdeg*Beta_rad + CnP_prps*P_rps*b_m/(2*true_airspeed_mps) + CnR_prps*R_rps*b_m/(2*true_airspeed_mps) + \
    CnBetadot_pdps*Betadot_ps*b_m/(2*true_airspeed_mps) + Cndela_pdeg*dela_deg + Cndelr_pdeg*delr_deg

    return Cn