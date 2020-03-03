import ROOT
from math import sqrt

PDG_MUON_MASS_T     = 0.1
PDG_MUON_MASS       = 0.1056583745
PDG_PION_MASS       = 0.13957018
PDG_KAON_MASS       = 0.493677
PDG_PROTON_MASS     = 0.9382720813
PDG_KSHORT_MASS     = 0.497611
PDG_KSHORT_DM       = 0.000013
PDG_KSHORT_TIME     = 0.8954 * 0.0000000001
PDG_KS_MASS         = PDG_KSHORT_MASS
PDG_LAMBDA_MASS     = 1.115683 
PDG_LAMBDA_DM       = 0.000006 
PDG_LAMBDA_TIME     = 2.632 * 0.0000000001
PDG_SIGMA0_MASS     = 1.192642 
PDG_XImunus_MASS    = 1.32171 
PDG_XImunus_DM      = 0.00007 
PDG_XImunus_TIME    = 1.639 * 0.0000000001
PDG_KSTAR_MASS      = 0.89581
PDG_KSTAR_GAMMA     = 0.0508
PDG_PHI_MASS        = 1.019461
PDG_PHI_GAMMA       = 0.004266
PDG_JPSI_MASS       = 3.096900
PDG_PSI2S_MASS      = 3.686097
PDG_X3872_MASS      = 3.87169
PDG_BU_MASS         = 5.27931
PDG_B0_MASS         = 5.27962
PDG_dm_BST_B        = 0.04518
PDG_BS_MASS         = 5.36682
PDG_BC_MASS         = 6.2751
PDG_LB_MASS         = 5.61951
PDG_XIB_MASS        = 5.7970
PDG_C               = 29979245800. ## in cm/c
PDG_LIFETIME_BU     = 1.638   * (10**-12) ## in s
PDG_LIFETIME_B0     = 1.52    * (10**-12) ## in s
PDG_LIFETIME_BS     = 1.511   * (10**-12) ## in s
PDG_LIFETIME_KS     = 0.8954  * (10**-10) ## in s
PDG_LIFETIME_LAMBDA = 2.632   * (10**-10) ## in s
PDG_LIFETIME_LB     = 1.466   * (10**-12) ## in s
PDG_BS2_MASS        = 5.83983
PDG_BS1_MASS        = 5.82878
PDG_B1_MASS         = 5.7249
PDG_B2ST_MASS       = 5.739
PDG_DMBstB          = 0.04538
PDG_OMEGA_MASS      = 1.67245

SAMEEVENT        = ROOT.RooRealVar('SAMEEVENT'      , 'SAMEEVENT', 0, 1)

Jpsi_Psi2S_m     = ROOT.RooRealVar('Jpsi_Psi2S_mass_Cmumu', 'M(#mu^{+} #mu^{-}), [GeV]', 2.9, 3.9)
Jpsi_Psi2S_pt    = ROOT.RooRealVar('Jpsi_Psi2S_pt_Cmumu'  , 'p_{T}(J\#psi), [GeV/c]'   , 0, 180)
Jpsi_Psi2S_pvdistsign = ROOT.RooRealVar('Jpsi_Psi2S_pvdistsignif2_Cmumu', '(J\#psi)_{pvdistsign}', 0, 800)

Lambda_m         = ROOT.RooRealVar('LAM_mass'       , 'M(#Lambda), [GeV]'        , 1.08, 1.17)
Lambda_pt        = ROOT.RooRealVar('LAM_pt'         , 'p_{T}(#Lambda), [GeV/c]'  , 0, 60)
Lambda_pr_charge = ROOT.RooRealVar('LAM_pr_charge'  , 'q(p)'                     , -1, 1)

Ks_m             = ROOT.RooRealVar('KS_mas1'        , 'M(#pi^{#pm} #pi^{#mp}), [GeV]', 0, 0.65)
B0_m             = ROOT.RooRealVar('B0_mass'        , 'M(J/#psi #pi^{#pm} #pi^{#mp}), [GeV]', 3, 20)
Omega_m          = ROOT.RooRealVar('XI_lamk_mass'   , 'M(#Lambda #K^{-}), [GeV]' , 1.2, 3.3)
Omegab_m         = ROOT.RooRealVar('OmegaB_mass'    , 'M(J/#psi #Lambda #K^{-}), [GeV]' , 5, 20)

Xi_m             = ROOT.RooRealVar('XI_mass_old'    , 'M(#Lambda #pi^{-}), [GeV]', 1.25, 1.39)
Xi_pt            = ROOT.RooRealVar('XI_pt'          , 'p_{T}(#Xi), [GeV/c]'      , 0, 80)
Xi_pi1_charge    = ROOT.RooRealVar('XI_PI_charge'   , 'q(#pi_{#Xi})'             , -1, 1)
Xi_pi_ips        = ROOT.RooRealVar('XI_PI_ips'      , '(#pi_{#Xi}_{ips})' , 0, 1000)

Xib_m            = ROOT.RooRealVar('XB_mass_Cjp_old', 'M(J/#psi #Xi^{-}), [GeV]' , 5.6, 6.0)
Xib_pvdistsign   = ROOT.RooRealVar('XB_pvdistsignif2_Cjp', '(#Xi_{b})_{pvdistsign}', 0, 1000)
Xib_pt           = ROOT.RooRealVar('XB_pt_Cjp'      , 'p_{T}(#Xi_{b}), [GeV/c]'    , 0, 200)
Xib_sum_m        = ROOT.RooRealVar('XB_mass_sum', 'M(J/#psi) + M(#Xi^{-}), [GeV]' , 0, 15)
Xib_delt_m       = ROOT.RooRealVar('XB_mass_delta', 'M(J/#psi #Xi^{-}) - M(J/#psi) - M(#Xi^{-}), [GeV]' , 0, 3)

Xib_star_m       = ROOT.RooRealVar('XB_star_mass_sum', 'M(#Xi_{b}^{-} #pi^{+}), [GeV]' , 5.7, 7.3)
Xib_star_pm_m    = ROOT.RooRealVar('XB_star_mass_sum_pm', 'M(#Xi_{b}^{-} #pi^{+}), [GeV]' , 5.7, 7.4)
Xib_star_delt_m  = ROOT.RooRealVar('XB_star_mass_delta', 'M(#Xi_{b}^{-} #pi^{+}) - M(#Xi_{b}^{-}) - M(#pi^{+}), [GeV]' , 0, 3)
Xib_star_delt_m2 = ROOT.RooRealVar('XB_star_mass_delta2', 'M(#Xi_{b}^{-} #pi^{+}) - M(#Xi_{b}^{-}) - M(#pi^{+}), [GeV]' , 0, 3)
Xib_pi2_pt        = ROOT.RooRealVar('XB_PI2_pt'       , 'p_{T}(#pi_{#Xi_{b}}), [GeV/c]'      , 0, 30)
Xib_pi2_charge    = ROOT.RooRealVar('XB_PI2_charge'   , 'q(#pi_{2#Xi_{b}})'                  , -1, 1)
Xib_pi2_eta       = ROOT.RooRealVar('XB_PI2_eta'      , 'eta_pi2'                            , -2.5, 2.5)

Xib_pipi_m       = ROOT.RooRealVar('XBPIPI_mass_sum', 'M(#Xi_{b}^{-} #pi^{+} #pi^{-}), [GeV]' , 6.077, 7.5)
Xib_pipi_delt_m  = ROOT.RooRealVar('XBPIPI_mass_delta', 'M(#Xi_{b}^{-} #pi^{+} #pi^{-}) - M(#Xi_{b}^{-}) - 2 M(#pi^{#pm}), [GeV]' , 0, 1.4)
Xib_pi3_charge   = ROOT.RooRealVar('XB_PI3_charge'   , 'q(#pi_{3#Xi_{b}})'                  , -1, 1)
Xib_pi3_pt       = ROOT.RooRealVar('XB_PI3_pt'       , 'p_{T}(#pi_{#Xi_{b}}), [GeV/c]'      , 0, 30)
Xib_pi3_eta      = ROOT.RooRealVar('XB_PI3_eta'      , 'eta_pi3'                            , -2.5, 2.5)

X_mass           = ROOT.RooRealVar('X_mass', 'M(#Xi_{b}^{-} #pi^{+}), [GeV]' , -0.1, 7.4)
Xib_star_pm_mX   = ROOT.RooRealVar('XB_star_mass_X', 'M(#Xi_{b}^{-} #pi^{+}), [GeV]' , 5.7, 7.4)
Xib_star_delt_mX = ROOT.RooRealVar('XB_star_mass_delta_X', 'M(#Xi_{b}^{-} #pi^{+}) - M(#Xi_{b}^{-}) - M(#pi^{+}), [GeV]' , 0, 1)

A_mass           = ROOT.RooRealVar('A_mass', 'M(#Xi_{b}^{-} #pi^{+}), [GeV]' , -1, 800)
Xib_star_pm_mA   = ROOT.RooRealVar('XB_star_mass_A', 'M(#Xi_{b}^{-} #pi^{+}), [GeV]' , 5.7, 7.4)
Xib_star_delt_mA = ROOT.RooRealVar('XB_star_mass_delta_A', 'M(#Xi_{b}^{-} #pi^{+}) - M(#Xi_{b}^{-}) - M(#pi^{+}), [GeV]' , 0, 1)

Xib_pipi_delt_mX = ROOT.RooRealVar('XBPIPI_mass_delta_X', 'M(#Xi_{b}^{-} #pi^{+} #pi^{-}) - M(#Xi_{b}^{-}) - 2 M(#pi^{#pm}), [GeV]' , 0, 2)
Xib_pipi_delt_mA = ROOT.RooRealVar('XBPIPI_mass_delta_A', 'M(#Xi_{b}^{-} #pi^{+} #pi^{-}) - M(#Xi_{b}^{-}) - 2 M(#pi^{#pm}), [GeV]' , 0, 2)

cos_PV_Xib       = ROOT.RooRealVar('XB_pvcos2_Cjp'  , 'cos(#Xi_{b}, PV)'         , -1, 1)
cos_Xib_Xi       = ROOT.RooRealVar('XI_XBcos2'      , 'cos(#Xi_{b}, #Xi)'        , -1, 1)
cos_Lam_Xi       = ROOT.RooRealVar('LAM_XIcos2_CL'  , 'cos(#Lambda, #Xi)'        , -1, 1)

# cos_Xib_Xi       = ROOT.RooRealVar('XI_XBcos2'      , 'cos(#Xi_{b}, #Xi)'        , -1, 1)
# LA_Bcos2         = ROOT.RooRealVar('LA_Bcos2'       , 'LA_Bcos2'                 , -1, 1)
# B_mass           = ROOT.RooRealVar('B_mass'         , 'M(J/#psi #Lambda #pi #pi), [GeV]', 4.5, 5.9)


lumi_13TeV = {2012: 19.6, 2016: 36.295, 2017: 42.14, 2018: 61.31}

### Tp = Tlab * sqrt(1-v2/c2) = d/v * sqrt(1-v2/c2)
### p = m*v / sqrt(1-v2/c2)
### sqrt(1-v2/c2) = m*v/p
### m2v2c2/p2c2 + v2/c2= 1
### v = c * sqrt(1/(1+m2c2/p2)) [=~ c * (1-m2c2/2p2)]
### Tp = d/v * sqrt(1-1/(1+m2c2/p2))
###


def DetachSignificance2(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2))


def DetachSignificance3(vtx, vtxE1, vtxE2):
    return sqrt( vtx.X()**2 / (vtxE1.X()**2 + vtxE2.X()**2) + vtx.Y()**2 / (vtxE1.Y()**2 + vtxE2.Y()**2) + vtx.Z()**2 / (vtxE1.Z()**2 + vtxE2.Z()**2))


def DirectionCos2 (v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() ) / (r1*r2 + 0.0000001)


def DirectionCos3(v1, v2):
    r1 = sqrt(v1.X()**2 + v1.Y()**2 + v1.Z()**2)
    r2 = sqrt(v2.X()**2 + v2.Y()**2 + v2.Z()**2)
    return ( v1.X() * v2.X() + v1.Y() * v2.Y() + v1.Z() * v2.Z() ) / (r1*r2 + 0.0000001)


def DirectionChi22(vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0
    dvtxE   = ROOT.TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) )
    Pscaled = P * (dvtx.Mag() / P.Mag())
    PscaledE= PE * (dvtx.Mag() / P.Mag())
    return DetachSignificance2 (Pscaled - dvtx, PscaledE, dvtxE)


def DirectionChi23(vtx0, vtx0E, vtx1, vtx1E, P, PE):
    dvtx    = vtx1 - vtx0 ## vertex difference
    dvtxE   = ROOT.TVector3( sqrt(vtx0E.X()**2 + vtx1E.X()**2), sqrt(vtx0E.Y()**2 + vtx1E.Y()**2), sqrt(vtx0E.Z()**2 + vtx1E.Z()**2) ) ## its error
    Pscaled = P * (dvtx.Mag() / P.Mag()) ## scaled momentum to be the same length as vertex difference
    PscaledE= PE * (dvtx.Mag() / P.Mag()) ## its error
    return DetachSignificance3 (Pscaled - dvtx, PscaledE, dvtxE)
