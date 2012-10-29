# produce the energy scale systematics plots for ggH130
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import ROOT
import glob
import math

tag = 'mc12a'
removeCrack = False

if (tag == 'mc11c'):
  input = glob.glob('/afs/cern.ch/user/v/vippolit/scratch0/testarea/AtlasProduction-17.5.0.1/tier2/output_test.17.root')
  outfile = ROOT.TFile('energy_scale_syst_2011.root', 'RECREATE')
elif (tag == 'mc12a'):
  input = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_v20/*ggH130*/*root*')
  outfile = ROOT.TFile('energy_scale_syst_2012b.root', 'RECREATE')

ROOT.gSystem.Load('libCint')
ROOT.gROOT.ProcessLine('.include egammaTRUNK')
ROOT.gROOT.ProcessLine('.L egammaTRUNK/Root/EnergyRescalerUpgrade.cxx+')
ROOT.gROOT.ProcessLine('.L egammaTRUNK/egammaAnalysisUtils/EnergyRescalerUpgrade.h+')

energyRescaler = ROOT.egRescaler.EnergyRescalerUpgrade()
if (tag == 'mc11c'):
  energyRescaler.Init('egammaTRUNK/share/EnergyRescalerData.root', '2011', 'es2010')
elif (tag == 'mc12a'):
  energyRescaler.Init('egammaTRUNK/share/EnergyRescalerData.root', '2012', 'es2010')

class syst_histos:
  def __init__(self, name):
    nbins = 2000
    minx = 0
    maxx = 1000

    self.name = name
    self.nominal = ROOT.TH1F('%s_syst_nominal' % (name), '', nbins, minx, maxx)
    self.PS1_up = ROOT.TH1F('%s_syst_PS1_up' % (name), '', nbins, minx, maxx)
    self.PS1_down = ROOT.TH1F('%s_syst_PS1_down' % (name), '', nbins, minx, maxx)
    self.PS2_up = ROOT.TH1F('%s_syst_PS2_up' % (name), '', nbins, minx, maxx)
    self.PS2_down = ROOT.TH1F('%s_syst_PS2_down' % (name), '', nbins, minx, maxx)
    self.MAT1_up = ROOT.TH1F('%s_syst_MAT1_up' % (name), '', nbins, minx, maxx)
    self.MAT1_down = ROOT.TH1F('%s_syst_MAT1_down' % (name), '', nbins, minx, maxx)
    self.MAT2_up = ROOT.TH1F('%s_syst_MAT2_up' % (name), '', nbins, minx, maxx)
    self.MAT2_down = ROOT.TH1F('%s_syst_MAT2_down' % (name), '', nbins, minx, maxx)
    self.ES_Z_up = ROOT.TH1F('%s_syst_ES_Z_up' % (name), '', nbins, minx, maxx)
    self.ES_Z_down = ROOT.TH1F('%s_syst_ES_Z_down' % (name), '', nbins, minx, maxx)
    self.ES_lowpt_up = ROOT.TH1F('%s_syst_ES_lowpt_up' % (name), '', nbins, minx, maxx)
    self.ES_lowpt_down = ROOT.TH1F('%s_syst_ES_lowpt_down' % (name), '', nbins, minx, maxx)
    self.MUONS_up = ROOT.TH1F('%s_syst_MUONS_up' % (name), '', nbins, minx, maxx)
    self.MUONS_down = ROOT.TH1F('%s_syst_MUONS_down' % (name), '', nbins, minx, maxx)
    self.lep_pt = ROOT.TH1F('%s_el_pt' % (name), '', 1000, 0, 1000);
    self.lep_PS1_up = ROOT.TH1F('%s_lep_deOe_syst_PS1_up' % (name), '', 1000, -2., 2.)
    self.lep_PS1_down = ROOT.TH1F('%s_lep_deOe_syst_PS1_down' % (name), '', 1000, -2., 2.)
    self.lep_PS2_up = ROOT.TH1F('%s_lep_deOe_syst_PS2_up' % (name), '', 1000, -2., 2.)
    self.lep_PS2_down = ROOT.TH1F('%s_lep_deOe_syst_PS2_down' % (name), '', 1000, -2., 2.)
    self.lep_MAT1_up = ROOT.TH1F('%s_lep_deOe_syst_MAT1_up' % (name), '', 1000, -2., 2.)
    self.lep_MAT1_down = ROOT.TH1F('%s_lep_deOe_syst_MAT1_down' % (name), '', 1000, -2., 2.)
    self.lep_MAT2_up = ROOT.TH1F('%s_lep_deOe_syst_MAT2_up' % (name), '', 1000, -2., 2.)
    self.lep_MAT2_down = ROOT.TH1F('%s_lep_deOe_syst_MAT2_down' % (name), '', 1000, -2., 2.)
    self.lep_ES_Z_up = ROOT.TH1F('%s_lep_deOe_syst_ES_Z_up' % (name), '', 1000, -2., 2.)
    self.lep_ES_Z_down = ROOT.TH1F('%s_lep_deOe_syst_ES_Z_down' % (name), '', 1000, -2., 2.)
    self.lep_ES_lowpt_up = ROOT.TH1F('%s_lep_deOe_syst_ES_lowpt_up' % (name), '', 1000, -2., 2.)
    self.lep_ES_lowpt_down = ROOT.TH1F('%s_lep_deOe_syst_ES_lowpt_down' % (name), '', 1000, -2., 2.)
    self.lep_MUONS_up = ROOT.TH1F('%s_lep_deOe_syst_MUONS_up' % (name), '', 1000, -2., 2.)
    self.lep_MUONS_down = ROOT.TH1F('%s_lep_deOe_syst_MUONS_down' % (name), '', 1000, -2., 2.)

  def set_global_ownership(self):
    self.nominal.SetDirectory(0)
    self.PS1_up.SetDirectory(0)
    self.PS1_down.SetDirectory(0)
    self.PS2_up.SetDirectory(0)
    self.PS2_down.SetDirectory(0)
    self.MAT1_up.SetDirectory(0)
    self.MAT1_down.SetDirectory(0)
    self.MAT2_up.SetDirectory(0)
    self.MAT2_down.SetDirectory(0)
    self.ES_Z_up.SetDirectory(0)
    self.ES_Z_down.SetDirectory(0)
    self.ES_lowpt_up.SetDirectory(0)
    self.ES_lowpt_down.SetDirectory(0)
    self.MUONS_up.SetDirectory(0)
    self.MUONS_down.SetDirectory(0)
    self.lep_pt.SetDirectory(0)
    self.lep_PS1_up.SetDirectory(0)
    self.lep_PS1_down.SetDirectory(0)
    self.lep_PS2_up.SetDirectory(0)
    self.lep_PS2_down.SetDirectory(0)
    self.lep_MAT1_up.SetDirectory(0)
    self.lep_MAT1_down.SetDirectory(0)
    self.lep_MAT2_up.SetDirectory(0)
    self.lep_MAT2_down.SetDirectory(0)
    self.lep_ES_Z_up.SetDirectory(0)
    self.lep_ES_Z_down.SetDirectory(0)
    self.lep_ES_lowpt_up.SetDirectory(0)
    self.lep_ES_lowpt_down.SetDirectory(0)
    self.lep_MUONS_up.SetDirectory(0)
    self.lep_MUONS_down.SetDirectory(0)

  def set_ownership(self, dir):
    self.nominal.SetDirectory(dir)
    self.PS1_up.SetDirectory(dir)
    self.PS1_down.SetDirectory(dir)
    self.PS2_up.SetDirectory(dir)
    self.PS2_down.SetDirectory(dir)
    self.MAT1_up.SetDirectory(dir)
    self.MAT1_down.SetDirectory(dir)
    self.MAT2_up.SetDirectory(dir)
    self.MAT2_down.SetDirectory(dir)
    self.ES_Z_up.SetDirectory(dir)
    self.ES_Z_down.SetDirectory(dir)
    self.ES_lowpt_up.SetDirectory(dir)
    self.ES_lowpt_down.SetDirectory(dir)
    self.MUONS_up.SetDirectory(dir)
    self.MUONS_down.SetDirectory(dir)
    self.lep_pt.SetDirectory(dir)
    self.lep_PS1_up.SetDirectory(dir)
    self.lep_PS1_down.SetDirectory(dir)
    self.lep_PS2_up.SetDirectory(dir)
    self.lep_PS2_down.SetDirectory(dir)
    self.lep_MAT1_up.SetDirectory(dir)
    self.lep_MAT1_down.SetDirectory(dir)
    self.lep_MAT2_up.SetDirectory(dir)
    self.lep_MAT2_down.SetDirectory(dir)
    self.lep_ES_Z_up.SetDirectory(dir)
    self.lep_ES_Z_down.SetDirectory(dir)
    self.lep_ES_lowpt_up.SetDirectory(dir)
    self.lep_ES_lowpt_down.SetDirectory(dir)
    self.lep_MUONS_up.SetDirectory(dir)
    self.lep_MUONS_down.SetDirectory(dir)

def fill_syst_histos(t, syst):
  poids = 1.

  # 4-momenta used for computing the new Higgs mass
  Z1_lepplus_4m_up = ROOT.TLorentzVector()
  Z1_lepminus_4m_up = ROOT.TLorentzVector()
  Z2_lepplus_4m_up = ROOT.TLorentzVector()
  Z2_lepminus_4m_up = ROOT.TLorentzVector()
  Z1_lepplus_4m_down = ROOT.TLorentzVector()
  Z1_lepminus_4m_down = ROOT.TLorentzVector()
  Z2_lepplus_4m_down = ROOT.TLorentzVector()
  Z2_lepminus_4m_down = ROOT.TLorentzVector()


  # nominal histogram
  syst.nominal.Fill(t.H_m/1000., poids)

  # PS1: use PSStatUp/Down for all and only the electrons between 0-1.5 (cl_eta)
  #
  Z1_lepplus_corr = False
  if (t.Z1_lepplus_m < 100):
    if (abs(t.Z1_lepplus_cl_eta) < 1.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS1_up.Fill(E_corrected_up / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS1_down.Fill(E_corrected_down / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))

      Z1_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)

      Z1_lepplus_corr = True
  if not Z1_lepplus_corr:
      Z1_lepplus_4m_up.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
  #
  Z1_lepminus_corr = False
  if (t.Z1_lepminus_m < 100):
    if (abs(t.Z1_lepminus_cl_eta) < 1.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS1_up.Fill(E_corrected_up / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS1_down.Fill(E_corrected_down / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))

      Z1_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)

      Z1_lepminus_corr = True
  if not Z1_lepminus_corr:
      Z1_lepminus_4m_up.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
  #
  Z2_lepplus_corr = False
  if (t.Z2_lepplus_m < 100):
    if (abs(t.Z2_lepplus_cl_eta) < 1.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS1_up.Fill(E_corrected_up / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS1_down.Fill(E_corrected_down / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))

      Z2_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)

      Z2_lepplus_corr = True
  if not Z2_lepplus_corr:
      Z2_lepplus_4m_up.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
  #
  Z2_lepminus_corr = False
  if (t.Z2_lepminus_m < 100):
    if (abs(t.Z2_lepminus_cl_eta) < 1.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS1_up.Fill(E_corrected_up / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS1_down.Fill(E_corrected_down / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))

      Z2_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

      Z2_lepminus_corr = True
  if not Z2_lepminus_corr:
      Z2_lepminus_4m_up.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

  syst.PS1_up.Fill((Z1_lepplus_4m_up + Z1_lepminus_4m_up + Z2_lepplus_4m_up + Z2_lepminus_4m_up).M()/1000., poids)
  syst.PS1_down.Fill((Z1_lepplus_4m_down + Z1_lepminus_4m_down + Z2_lepplus_4m_down + Z2_lepminus_4m_down).M()/1000., poids)


  # PS2: use PSStatUp/Down for all and only the electrons between 1.5-1.8 (cl_eta)
  #
  Z1_lepplus_corr = False
  if (t.Z1_lepplus_m < 100):
    if (abs(t.Z1_lepplus_cl_eta) > 1.5 and abs(t.Z1_lepplus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS2_up.Fill(E_corrected_up / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS2_down.Fill(E_corrected_down / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))

      Z1_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)

      Z1_lepplus_corr = True
  if not Z1_lepplus_corr:
      Z1_lepplus_4m_up.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
  #
  Z1_lepminus_corr = False
  if (t.Z1_lepminus_m < 100):
    if (abs(t.Z1_lepminus_cl_eta) > 1.5 and abs(t.Z1_lepminus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS2_up.Fill(E_corrected_up / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS2_down.Fill(E_corrected_down / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))

      Z1_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)

      Z1_lepminus_corr = True
  if not Z1_lepminus_corr:
      Z1_lepminus_4m_up.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
  #
  Z2_lepplus_corr = False
  if (t.Z2_lepplus_m < 100):
    if (abs(t.Z2_lepplus_cl_eta) > 1.5 and abs(t.Z2_lepplus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS2_up.Fill(E_corrected_up / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS2_down.Fill(E_corrected_down / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))

      Z2_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)

      Z2_lepplus_corr = True
  if not Z2_lepplus_corr:
      Z2_lepplus_4m_up.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
  #
  Z2_lepminus_corr = False
  if (t.Z2_lepminus_m < 100):
    if (abs(t.Z2_lepminus_cl_eta) > 1.5 and abs(t.Z2_lepminus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.PSStatDown)
      syst.lep_PS2_up.Fill(E_corrected_up / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_PS2_down.Fill(E_corrected_down / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))

      Z2_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

      Z2_lepminus_corr = True
  if not Z2_lepminus_corr:
      Z2_lepminus_4m_up.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

  syst.PS2_up.Fill((Z1_lepplus_4m_up + Z1_lepminus_4m_up + Z2_lepplus_4m_up + Z2_lepminus_4m_up).M()/1000., poids)
  syst.PS2_down.Fill((Z1_lepplus_4m_down + Z1_lepminus_4m_down + Z2_lepplus_4m_down + Z2_lepminus_4m_down).M()/1000., poids)
 
  # MAT1: use R12StatUp/Down for all and only the electrons between 0-1.8 (cl_eta)
  #
  Z1_lepplus_corr = False
  if (t.Z1_lepplus_m < 100):
    if (abs(t.Z1_lepplus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT1_up.Fill(E_corrected_up / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT1_down.Fill(E_corrected_down / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))

      Z1_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)

      Z1_lepplus_corr = True
  if not Z1_lepplus_corr:
      Z1_lepplus_4m_up.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
  #
  Z1_lepminus_corr = False
  if (t.Z1_lepminus_m < 100):
    if (abs(t.Z1_lepminus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT1_up.Fill(E_corrected_up / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT1_down.Fill(E_corrected_down / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))

      Z1_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)

      Z1_lepminus_corr = True
  if not Z1_lepminus_corr:
      Z1_lepminus_4m_up.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
  #
  Z2_lepplus_corr = False
  if (t.Z2_lepplus_m < 100):
    if (abs(t.Z2_lepplus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT1_up.Fill(E_corrected_up / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT1_down.Fill(E_corrected_down / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))

      Z2_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)

      Z2_lepplus_corr = True
  if not Z2_lepplus_corr:
      Z2_lepplus_4m_up.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
  #
  Z2_lepminus_corr = False
  if (t.Z2_lepminus_m < 100):
    if (abs(t.Z2_lepminus_cl_eta) < 1.8):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT1_up.Fill(E_corrected_up / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT1_down.Fill(E_corrected_down / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))

      Z2_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

      Z2_lepminus_corr = True
  if not Z2_lepminus_corr:
      Z2_lepminus_4m_up.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

  syst.MAT1_up.Fill((Z1_lepplus_4m_up + Z1_lepminus_4m_up + Z2_lepplus_4m_up + Z2_lepminus_4m_up).M()/1000., poids)
  syst.MAT1_down.Fill((Z1_lepplus_4m_down + Z1_lepminus_4m_down + Z2_lepplus_4m_down + Z2_lepminus_4m_down).M()/1000., poids)


  # MAT2: use R12StatUp/Down for all and only the electrons between 1.8-2.5 (cl_eta)
  #
  Z1_lepplus_corr = False
  if (t.Z1_lepplus_m < 100):
    if (abs(t.Z1_lepplus_cl_eta) > 1.8 and abs(t.Z1_lepplus_cl_eta) < 2.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT2_up.Fill(E_corrected_up / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT2_down.Fill(E_corrected_down / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))

      Z1_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)

      Z1_lepplus_corr = True
  if not Z1_lepplus_corr:
      Z1_lepplus_4m_up.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
  #
  Z1_lepminus_corr = False
  if (t.Z1_lepminus_m < 100):
    if (abs(t.Z1_lepminus_cl_eta) > 1.8 and abs(t.Z1_lepminus_cl_eta) < 2.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT2_up.Fill(E_corrected_up / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT2_down.Fill(E_corrected_down / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))

      Z1_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)

      Z1_lepminus_corr = True
  if not Z1_lepminus_corr:
      Z1_lepminus_4m_up.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
  #
  Z2_lepplus_corr = False
  if (t.Z2_lepplus_m < 100):
    if (abs(t.Z2_lepplus_cl_eta) > 1.8 and abs(t.Z2_lepplus_cl_eta) < 2.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT2_up.Fill(E_corrected_up / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT2_down.Fill(E_corrected_down / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))

      Z2_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)

      Z2_lepplus_corr = True
  if not Z2_lepplus_corr:
      Z2_lepplus_4m_up.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
  #
  Z2_lepminus_corr = False
  if (t.Z2_lepminus_m < 100):
    if (abs(t.Z2_lepminus_cl_eta) > 1.8 and abs(t.Z2_lepminus_cl_eta) < 2.5):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.R12StatDown)
      syst.lep_MAT2_up.Fill(E_corrected_up / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_MAT2_down.Fill(E_corrected_down / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))

      Z2_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

      Z2_lepminus_corr = True
  if not Z2_lepminus_corr:
      Z2_lepminus_4m_up.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

  syst.MAT2_up.Fill((Z1_lepplus_4m_up + Z1_lepminus_4m_up + Z2_lepplus_4m_up + Z2_lepminus_4m_up).M()/1000., poids)
  syst.MAT2_down.Fill((Z1_lepplus_4m_down + Z1_lepminus_4m_down + Z2_lepplus_4m_down + Z2_lepminus_4m_down).M()/1000., poids)



  # ES_Z: use ZeeMethodUp/Down for all the electrons
  #
  Z1_lepplus_corr = False
  if (t.Z1_lepplus_m < 100):
    E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodUp)
    E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodDown)
    syst.lep_pt.Fill(t.Z1_lepplus_cl_E_calibsmearnoscale/1000.) # filled here since it's the only one for all electrons :)
    syst.lep_ES_Z_up.Fill(E_corrected_up / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)
    syst.lep_ES_Z_down.Fill(E_corrected_down / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)

    pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))
    pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))

    Z1_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
    Z1_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)

    Z1_lepplus_corr = True
  if not Z1_lepplus_corr:
    Z1_lepplus_4m_up.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
    Z1_lepplus_4m_down.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
  #
  Z1_lepminus_corr = False
  if (t.Z1_lepminus_m < 100):
    E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodUp)
    E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodDown)
    syst.lep_pt.Fill(t.Z1_lepminus_cl_E_calibsmearnoscale/1000.)
    syst.lep_ES_Z_up.Fill(E_corrected_up / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)
    syst.lep_ES_Z_down.Fill(E_corrected_down / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)

    pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))
    pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))

    Z1_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
    Z1_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)

    Z1_lepminus_corr = True
  if not Z1_lepminus_corr:
    Z1_lepminus_4m_up.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
    Z1_lepminus_4m_down.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
  #
  Z2_lepplus_corr = False
  if (t.Z2_lepplus_m < 100):
    E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodUp)
    E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodDown)
    syst.lep_pt.Fill(t.Z2_lepplus_cl_E_calibsmearnoscale/1000.)
    syst.lep_ES_Z_up.Fill(E_corrected_up / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)
    syst.lep_ES_Z_down.Fill(E_corrected_down / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)

    pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))
    pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))

    Z2_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
    Z2_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)

    Z2_lepplus_corr = True
  if not Z2_lepplus_corr:
    Z2_lepplus_4m_up.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
    Z2_lepplus_4m_down.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
  #
  Z2_lepminus_corr = False
  if (t.Z2_lepminus_m < 100):
    E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodUp)
    E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.ZeeMethodDown)
    syst.lep_pt.Fill(t.Z2_lepminus_cl_E_calibsmearnoscale/1000.)
    syst.lep_ES_Z_up.Fill(E_corrected_up / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)
    syst.lep_ES_Z_down.Fill(E_corrected_down / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)

    pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))
    pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))

    Z2_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
    Z2_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

    Z2_lepminus_corr = True
  if not Z2_lepminus_corr:
    Z2_lepminus_4m_up.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
    Z2_lepminus_4m_down.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

  syst.ES_Z_up.Fill((Z1_lepplus_4m_up + Z1_lepminus_4m_up + Z2_lepplus_4m_up + Z2_lepminus_4m_up).M()/1000., poids)
  syst.ES_Z_down.Fill((Z1_lepplus_4m_down + Z1_lepminus_4m_down + Z2_lepplus_4m_down + Z2_lepminus_4m_down).M()/1000., poids)


  # ES_lowpt: use PSStatUp/Down for all and only the electrons between 1.5-1.8 (cl_eta)
  #
  Z1_lepplus_corr = False
  if (t.Z1_lepplus_m < 100):
    if (t.Z1_lepplus_cl_E_calibsmearnoscale < 20000):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepplus_cl_eta, t.Z1_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtDown)
      syst.lep_ES_lowpt_up.Fill(E_corrected_up / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_ES_lowpt_down.Fill(E_corrected_down / t.Z1_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepplus_eta)))

      Z1_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)

      Z1_lepplus_corr = True
  if not Z1_lepplus_corr:
      Z1_lepplus_4m_up.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
      Z1_lepplus_4m_down.SetPtEtaPhiM(t.Z1_lepplus_pt, t.Z1_lepplus_eta, t.Z1_lepplus_phi, t.Z1_lepplus_m)
  #
  Z1_lepminus_corr = False
  if (t.Z1_lepminus_m < 100):
    if (t.Z1_lepminus_cl_E_calibsmearnoscale < 20000):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z1_lepminus_cl_eta, t.Z1_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtDown)
      syst.lep_ES_lowpt_up.Fill(E_corrected_up / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_ES_lowpt_down.Fill(E_corrected_down / t.Z1_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z1_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z1_lepminus_eta)))

      Z1_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)

      Z1_lepminus_corr = True
  if not Z1_lepminus_corr:
      Z1_lepminus_4m_up.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
      Z1_lepminus_4m_down.SetPtEtaPhiM(t.Z1_lepminus_pt, t.Z1_lepminus_eta, t.Z1_lepminus_phi, t.Z1_lepminus_m)
  #
  Z2_lepplus_corr = False
  if (t.Z2_lepplus_m < 100):
    if (t.Z2_lepplus_cl_E_calibsmearnoscale < 20000):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepplus_cl_eta, t.Z2_lepplus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtDown)
      syst.lep_ES_lowpt_up.Fill(E_corrected_up / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)
      syst.lep_ES_lowpt_down.Fill(E_corrected_down / t.Z2_lepplus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepplus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepplus_eta)))

      Z2_lepplus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)

      Z2_lepplus_corr = True
  if not Z2_lepplus_corr:
      Z2_lepplus_4m_up.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
      Z2_lepplus_4m_down.SetPtEtaPhiM(t.Z2_lepplus_pt, t.Z2_lepplus_eta, t.Z2_lepplus_phi, t.Z2_lepplus_m)
  #
  Z2_lepminus_corr = False
  if (t.Z2_lepminus_m < 100):
    if (t.Z2_lepminus_cl_E_calibsmearnoscale < 20000):
      E_corrected_up = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtUp)
      E_corrected_down = energyRescaler.applyEnergyCorrection(t.Z2_lepminus_cl_eta, t.Z2_lepminus_cl_E_calibsmearnoscale, 0, ROOT.egRescaler.EnergyRescalerUpgrade.LowPtDown)
      syst.lep_ES_lowpt_up.Fill(E_corrected_up / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)
      syst.lep_ES_lowpt_down.Fill(E_corrected_down / t.Z2_lepminus_cl_E_calibsmearnoscale - 1.)

      pt_corrected_up = math.sqrt(math.pow(E_corrected_up, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))
      pt_corrected_down = math.sqrt(math.pow(E_corrected_down, 2) - math.pow(t.Z2_lepminus_m, 2)) * math.sin(2 * math.atan(math.exp(-t.Z2_lepminus_eta)))

      Z2_lepminus_4m_up.SetPtEtaPhiM(pt_corrected_up, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(pt_corrected_down, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

      Z2_lepminus_corr = True
  if not Z2_lepminus_corr:
      Z2_lepminus_4m_up.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)
      Z2_lepminus_4m_down.SetPtEtaPhiM(t.Z2_lepminus_pt, t.Z2_lepminus_eta, t.Z2_lepminus_phi, t.Z2_lepminus_m)

  syst.ES_lowpt_up.Fill((Z1_lepplus_4m_up + Z1_lepminus_4m_up + Z2_lepplus_4m_up + Z2_lepminus_4m_up).M()/1000., poids)
  syst.ES_lowpt_down.Fill((Z1_lepplus_4m_down + Z1_lepminus_4m_down + Z2_lepplus_4m_down + Z2_lepminus_4m_down).M()/1000., poids)
 




#### THE ACTUAL SEQUENCE

syst_run = {}

channels = ['4mu', '2mu2e', '2e2mu', '4e']

for chan in channels:
  syst_run[chan] = {}


crack_counter = 0

for file in input:
  print 'VERBOSE: opening %s' % file
  f = ROOT.TFile.Open(file)

  if (f):
    t = f.Get('candidates')
     
    if (t):
      for entry in xrange(t.GetEntries()):
        t.GetEntry(entry)

        if (t.run not in syst_run[channels[t.type]].keys()):
	  print 'initializing for %d' % t.run
          syst_run[channels[t.type]][t.run] = syst_histos('ggH130_%s' % channels[t.type])
          syst_run[channels[t.type]][t.run].set_ownership(outfile)

	# crack removal
        Z1_lepplus_e_crk = (t.Z1_lepplus_m < 100 and 1.37 < abs(t.Z1_lepplus_cl_eta) < 1.52)
        Z1_lepminus_e_crk = (t.Z1_lepminus_m < 100 and 1.37 < abs(t.Z1_lepminus_cl_eta) < 1.52)
        Z2_lepplus_e_crk = (t.Z2_lepplus_m < 100 and 1.37 < abs(t.Z2_lepplus_cl_eta) < 1.52)
        Z2_lepminus_e_crk = (t.Z2_lepminus_m < 100 and 1.37 < abs(t.Z2_lepminus_cl_eta) < 1.52)

	if (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
	  if (t.selected == 1):
	    crack_counter += 1
	    if (removeCrack):
	      continue
	else: # REMOVE ME: if we don't want to remove the crack, we keep the crack alone
	  if (not removeCrack):
	    continue

	
	# fill histos

	if (t.selected == 1):
          fill_syst_histos(t, syst_run[channels[t.type]][t.run])
    else:
      print 'WARNING: tree not found'
  else:
    print 'WARNING: file not found'

print ''
if (removeCrack):
  print 'WARNING: removed %d selected events since at least one electron was in the crack' % crack_counter
else:
  print 'INFO: found %d selected events with at least an electron in the crack' % crack_counter

outfile.Write()
outfile.Close()
