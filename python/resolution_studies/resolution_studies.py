# This script is used to perform resolution studies.
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", "--filename", dest="filename", help="write output on ROOT file starting by FILE", metavar="FILE", default="test")
parser.add_option("-c", "--constraint", action="store_true", default=False, help="apply mass constraint", dest="doConstraint")
parser.add_option("-s", "--nocategories", action="store_false", default=True, help="simply run on a whole category", dest="doCategories")
parser.add_option("-b", "--brem", action="store_true", default=False, help="use brem categories for H->4e", dest="doBrem")
parser.add_option("-t", "--notruth", action="store_false", default=True, help="do not fit with RooFit the m4l-mH distributions", dest="doTruthFit")
parser.add_option("-k", "--truth_categories", action="store_true", default=False, help="apply TruthFit in categories", dest="doTruthFitInCategories")
parser.add_option("-r", "--no_rescale", action="store_false", default=True, help="do not apply Higgs mass rescale to sample mass", dest="doRescaleToSampleMass")
parser.add_option("-v", "--version", dest="version", help="use configuration for MC tag TAG (default is mc11c)", metavar="TAG", default="mc11c")
parser.add_option("-m", "--merge-2e2mu", action="store_true", dest="merge_2e2mu", help="merge 2e2mu and 2mu2e", default=False)

(options, args) = parser.parse_args()

doApplyMassConstraint = options.doConstraint
doCategories = options.doCategories
doBrem = options.doBrem and doCategories
doTruthFit = options.doTruthFit
doTruthFitInCategories = options.doTruthFitInCategories and doTruthFit
doRescaleToSampleMass = options.doRescaleToSampleMass
tag = options.version
merge_2e2mu = options.merge_2e2mu

output_filename_base = options.filename

import sys
sys.path.append('./Higgs4lepAnalysis/python/resolution_studies')

from ROOT import *


from common_utils import *

if (merge_2e2mu):
  types = ['4mu', '2e2mu', '4e']

# options
doUseApproxFractionOutside2Sigma = True
doNormalize = True
doSeparateOutputs = False

if merge_2e2mu:
  output_filename_base = output_filename_base + '_merge2e2mu'
if (doApplyMassConstraint):
  output_filename = output_filename_base + '.constrained.root'
else:
  output_filename = output_filename_base + '.unconstrained.root'

# print options
print '**********************'
print '  output_filename        = ' + str(output_filename)
print '  doApplyMassConstraint  = ' + str(doApplyMassConstraint)
print '  doCategories           = ' + str(doCategories)
print '  doBrem                 = ' + str(doBrem)
print '  doRescaleToSampleMass  = ' + str(doRescaleToSampleMass)
print '  doNormalize            = ' + str(doNormalize)
print '  doSeparateOutputs      = ' + str(doSeparateOutputs)
print '  doTruthFit             = ' + str(doTruthFit)
print '  doTruthFitInCategories = ' + str(doTruthFitInCategories)
print '  tag                    = ' + str(tag)
print '  merge_2e2mu            = ' + str(merge_2e2mu)
print '**********************'

# storage elements
histo_m4l = allocate_space()
histo_m4l_rescaled = allocate_space()
histo_m4l_truth = allocate_space()
histo_m4l_rescaled_truth = allocate_space()
count_plain = allocate_space()
count_weighted = allocate_space()
fit_result = allocate_space()
fit_result_rescaled = allocate_space()
fit_result_truth = allocate_space()
fit_result_rescaled_truth = allocate_space()
population = allocate_space()
histo_mZ1 = allocate_space()
histo_mZ2 = allocate_space()
histo_m4l_minus_truth = allocate_space()
histo_mZ1_minus_truth = allocate_space()
histo_mZ2_minus_truth = allocate_space()

# summary plots for each channel+category
rescale_vs_mH = {}
sigma_vs_mH = {}
fwhm_vs_mH = {}
rescale_vs_mH_rescaled = {}
sigma_vs_mH_rescaled = {}
fwhm_vs_mH_rescaled = {}
fwhm_vs_mH_truth = {}
fwhm_vs_mH_rescaled_truth = {}

# open the output
output_file = TFile(output_filename, 'RECREATE')

# define categories
categories = {}
if (doCategories):
  categories["4mu"] = ["all", "bbbb", "bbb", "bb", "other"]
  if (not doBrem):
    categories["4e"] =  ["all", "bbbb", "onecrk", "bbb", "other"]
  else:
    categories["4e"] =  ["all",
                         "lblb_bbbb", "lblb_onecrk", "lblb_bbb", "lblb_other",
                         "hbhb_bbbb", "hbhb_onecrk", "hbhb_bbb", "hbhb_other",
                         "lbhb_bbbb", "lbhb_onecrk", "lbhb_bbb", "lbhb_other"]
  categories["2e2mu"] = ["all",
                         "onecrk_any",
                         "bb_bb", "bb_other",
          	       "other_other"]
  categories["2mu2e"] = ["all",
                         "any_onecrk",
                         "bb_bb",
                         "other_bb", "other_other"]
else:
  categories["4mu"] = ["all"]
  categories["4e"]  = ["all"]
  categories["2e2mu"]  = ["all"]
  categories["2mu2e"]  = ["all"]


# implement categories
def get_category(tree):
  # locations
  Z1_lepplus_mu_b = (abs(tree.Z1_lepplus_eta) < 1.05)
  Z1_lepplus_mu_e = (abs(tree.Z1_lepplus_eta) > 1.05)
  Z1_lepplus_e_b = (abs(tree.Z1_lepplus_eta) < 1.37)
  Z1_lepplus_e_crk = (1.37 < abs(tree.Z1_lepplus_eta) < 1.52)
  Z1_lepplus_e_e = (1.52 < abs(tree.Z1_lepplus_eta) < 2.47)
  Z1_lepminus_mu_b = (abs(tree.Z1_lepminus_eta) < 1.05)
  Z1_lepminus_mu_e = (abs(tree.Z1_lepminus_eta) > 1.05)
  Z1_lepminus_e_b = (abs(tree.Z1_lepminus_eta) < 1.37)
  Z1_lepminus_e_crk = (1.37 < abs(tree.Z1_lepminus_eta) < 1.52)
  Z1_lepminus_e_e = (1.52 < abs(tree.Z1_lepminus_eta) < 2.47)
  Z2_lepplus_mu_b = (abs(tree.Z2_lepplus_eta) < 1.05)
  Z2_lepplus_mu_e = (abs(tree.Z2_lepplus_eta) > 1.05)
  Z2_lepplus_e_b = (abs(tree.Z2_lepplus_eta) < 1.37)
  Z2_lepplus_e_crk = (1.37 < abs(tree.Z2_lepplus_eta) < 1.52)
  Z2_lepplus_e_e = (1.52 < abs(tree.Z2_lepplus_eta) < 2.47)
  Z2_lepminus_mu_b = (abs(tree.Z2_lepminus_eta) < 1.05)
  Z2_lepminus_mu_e = (abs(tree.Z2_lepminus_eta) > 1.05)
  Z2_lepminus_e_b = (abs(tree.Z2_lepminus_eta) < 1.37)
  Z2_lepminus_e_crk = (1.37 < abs(tree.Z2_lepminus_eta) < 1.52)
  Z2_lepminus_e_e = (1.52 < abs(tree.Z2_lepminus_eta) < 2.47)

  # brem
  Z1_lepplus_e_lb = (tree.Z1_lepplus_GSF_dp < 0.3)
  Z1_lepplus_e_hb = (tree.Z1_lepplus_GSF_dp > 0.3)
  Z1_lepminus_e_lb = (tree.Z1_lepminus_GSF_dp < 0.3)
  Z1_lepminus_e_hb = (tree.Z1_lepminus_GSF_dp > 0.3)
  Z2_lepplus_e_lb = (tree.Z2_lepplus_GSF_dp < 0.3)
  Z2_lepplus_e_hb = (tree.Z2_lepplus_GSF_dp > 0.3)
  Z2_lepminus_e_lb = (tree.Z2_lepminus_GSF_dp < 0.3)
  Z2_lepminus_e_hb = (tree.Z2_lepminus_GSF_dp > 0.3)

  if (tree.type == 0): # 4mu
    if (Z1_lepplus_mu_b and Z1_lepminus_mu_b and Z2_lepplus_mu_b and Z2_lepminus_mu_b):
      return 'bbbb'
    elif (Z1_lepplus_mu_b + Z1_lepminus_mu_b + Z2_lepplus_mu_b + Z2_lepminus_mu_b == 3):
      return 'bbb'
    elif (Z1_lepplus_mu_b + Z1_lepminus_mu_b + Z2_lepplus_mu_b + Z2_lepminus_mu_b == 2):
      return 'bb'
    else:
      return 'other'
  elif (tree.type == 2): # 2e2mu
    if (Z1_lepplus_e_crk or Z1_lepminus_e_crk):
      return 'onecrk_any'
    elif (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_mu_b and Z2_lepminus_mu_b):
      return 'bb_bb'
    elif (Z1_lepplus_e_b and Z1_lepminus_e_b):
      return 'bb_other'
    else:
      return 'other_other'
  elif (tree.type == 1): # 2mu2e
    if (Z2_lepplus_e_crk or Z2_lepminus_e_crk):
      return 'any_onecrk'
    elif (Z2_lepplus_e_b and Z2_lepminus_e_b and Z1_lepplus_mu_b and Z1_lepminus_mu_b):
      return 'bb_bb'
    elif (Z2_lepplus_e_b and Z2_lepminus_e_b):
      return 'other_bb'
    else:
      return 'other_other'
  elif (tree.type == 3): # 4e
    if (not doBrem):
      if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
        return 'bbbb'
      elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
        return 'onecrk'
      elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
        return 'bbb'
      else:
        return 'other'
    else:
      if (Z1_lepplus_e_lb and Z1_lepminus_e_lb and Z2_lepplus_e_lb and Z2_lepminus_e_lb):
        if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
          return 'lblb_bbbb'
        elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
          return 'lblb_onecrk'
        elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
          return 'lblb_bbb'
        else:
          return 'lblb_other'
      elif ((Z1_lepplus_e_hb or Z1_lepminus_e_hb) and (Z2_lepplus_e_hb or Z2_lepminus_e_hb)):
        if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
          return 'hbhb_bbbb'
        elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
          return 'hbhb_onecrk'
        elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
          return 'hbhb_bbb'
        else:
          return 'hbhb_other'
      else:
        if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
          return 'lbhb_bbbb'
        elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
          return 'lbhb_onecrk'
        elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
          return 'lbhb_bbb'
        else:
          return 'lbhb_other'
  else:
    return 'bug'


# initialize storage elements
for sample in samples:
  for type in types:
    for cat in categories[type]:
      count_plain[sample.mass][type][cat] = 0.
      count_weighted[sample.mass][type][cat] = 0.

      histo_m4l[sample.mass][type][cat] = TH1D('histo_m4l_%d_%s_%s' % (sample.mass, type, cat), '', 140, sample.min, sample.max)
      histo_m4l[sample.mass][type][cat].Sumw2()

      histo_m4l_rescaled[sample.mass][type][cat] = TH1D('histo_m4l_rescaled_%d_%s_%s' % (sample.mass, type, cat), '', 140, sample.min, sample.max)
      histo_m4l_rescaled[sample.mass][type][cat].Sumw2()

      if (cat != 'all'):
        if (sample.mass < 360):
          true_min = -30
          true_max = 15
        else:
          true_min = -30
          true_max = 30
      else:
        if (sample.mass < 130):
          true_min = -10
          true_max = 5
        elif (sample.mass < 180):
          true_min = -10
          true_max = 7
        elif (sample.mass < 360):
          true_min = -20
          true_max = 10
        elif (sample.mass < 460):
          true_min = -30
          true_max = 20
        else:
          true_min = -30
          true_max = 30
      true_nbins = int((true_max - true_min) / 0.5) # 500 MeV binning
      histo_m4l_truth[sample.mass][type][cat] = TH1D('histo_m4l_truth_%d_%s_%s' % (sample.mass, type, cat), '', true_nbins, true_min, true_max)
      histo_m4l_truth[sample.mass][type][cat].Sumw2()

      histo_m4l_rescaled_truth[sample.mass][type][cat] = TH1D('histo_m4l_rescaled_truth_%d_%s_%s' % (sample.mass, type, cat), '', true_nbins, true_min, true_max)
      histo_m4l_rescaled_truth[sample.mass][type][cat].Sumw2()

      histo_mZ1[sample.mass][type][cat] = TH1D('histo_mZ1_%d_%s_%s' % (sample.mass, type, cat), '', 100, 66, 116)
      histo_mZ1[sample.mass][type][cat].Sumw2()
      histo_mZ2[sample.mass][type][cat] = TH1D('histo_mZ2_%d_%s_%s' % (sample.mass, type, cat), '', 100, 66, 116)
      histo_mZ2[sample.mass][type][cat].Sumw2()
      histo_m4l_minus_truth[sample.mass][type][cat] = TH1D('histo_m4l_minus_truth_%d_%s_%s' % (sample.mass, type, cat), '', 200, -100, 100)
      histo_m4l_minus_truth[sample.mass][type][cat].Sumw2()
      histo_mZ1_minus_truth[sample.mass][type][cat] = TH1D('histo_mZ1_minus_truth_%d_%s_%s' % (sample.mass, type, cat), '', 120, -30, 30)
      histo_mZ1_minus_truth[sample.mass][type][cat].Sumw2()
      histo_mZ2_minus_truth[sample.mass][type][cat] = TH1D('histo_mZ2_minus_truth_%d_%s_%s' % (sample.mass, type, cat), '', 120, -30, 30)
      histo_mZ2_minus_truth[sample.mass][type][cat].Sumw2()

# initialize summary plots
for type in types:
  color_index = 40
  rescale_vs_mH[type] = {}
  sigma_vs_mH[type] = {}
  fwhm_vs_mH[type] = {}
  rescale_vs_mH_rescaled[type] = {}
  sigma_vs_mH_rescaled[type] = {}
  fwhm_vs_mH_rescaled[type] = {}
  fwhm_vs_mH_truth[type] = {}
  fwhm_vs_mH_rescaled_truth[type] = {}
  for cat in categories[type]:
    rescale_vs_mH[type][cat] = TGraphErrors()
    sigma_vs_mH[type][cat] = TGraphErrors()
    fwhm_vs_mH[type][cat] = TGraphErrors()
    rescale_vs_mH_rescaled[type][cat] = TGraphErrors()
    sigma_vs_mH_rescaled[type][cat] = TGraphErrors()
    fwhm_vs_mH_rescaled[type][cat] = TGraphErrors()
    fwhm_vs_mH_truth[type][cat] = TGraphErrors()
    fwhm_vs_mH_rescaled_truth[type][cat] = TGraphErrors()

    rescale_vs_mH[type][cat].SetName('rescale_vs_mH_' + type + '_' + cat)
    sigma_vs_mH[type][cat].SetName('sigma_vs_mH_' + type + '_' + cat)
    fwhm_vs_mH[type][cat].SetName('fwhm_vs_mH_' + type + '_' + cat)
    rescale_vs_mH_rescaled[type][cat].SetName('rescale_vs_mH_rescaled_' + type + '_' + cat)
    sigma_vs_mH_rescaled[type][cat].SetName('sigma_vs_mH_rescaled_' + type + '_' + cat)
    fwhm_vs_mH_rescaled[type][cat].SetName('fwhm_vs_mH_rescaled_' + type + '_' + cat)
    fwhm_vs_mH_truth[type][cat].SetName('fwhm_vs_mH_truth_' + type + '_' + cat)
    fwhm_vs_mH_rescaled_truth[type][cat].SetName('fwhm_vs_mH_rescaled_truth_' + type + '_' + cat)
    
    rescale_vs_mH[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')
    sigma_vs_mH[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')
    fwhm_vs_mH[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')
    rescale_vs_mH_rescaled[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')
    sigma_vs_mH_rescaled[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')
    fwhm_vs_mH_rescaled[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')
    fwhm_vs_mH_truth[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')
    fwhm_vs_mH_rescaled_truth[type][cat].GetXaxis().SetTitle('m_{H} [GeV]')

    rescale_vs_mH[type][cat].GetYaxis().SetTitle('m_{H} / m')
    sigma_vs_mH[type][cat].GetYaxis().SetTitle('#sigma [GeV]')
    fwhm_vs_mH[type][cat].GetYaxis().SetTitle('FWHM [GeV]')
    rescale_vs_mH_rescaled[type][cat].GetYaxis().SetTitle('m_{H} / m')
    sigma_vs_mH_rescaled[type][cat].GetYaxis().SetTitle('#sigma [GeV]')
    fwhm_vs_mH_rescaled[type][cat].GetYaxis().SetTitle('FWHM [GeV]')
    fwhm_vs_mH_truth[type][cat].GetYaxis().SetTitle('FWHM [GeV]')
    fwhm_vs_mH_rescaled_truth[type][cat].GetYaxis().SetTitle('FWHM [GeV]')

    rescale_vs_mH[type][cat].SetMarkerColor(color_index)
    sigma_vs_mH[type][cat].SetMarkerColor(color_index)
    fwhm_vs_mH[type][cat].SetMarkerColor(color_index)
    rescale_vs_mH_rescaled[type][cat].SetMarkerColor(color_index)
    sigma_vs_mH_rescaled[type][cat].SetMarkerColor(color_index)
    fwhm_vs_mH_rescaled[type][cat].SetMarkerColor(color_index)
    fwhm_vs_mH_truth[type][cat].SetMarkerColor(color_index)
    fwhm_vs_mH_rescaled_truth[type][cat].SetMarkerColor(color_index)
    color_index = color_index + 1



# fits
def CSC_fit(input):
  output_file.cd()
  print '' 
  print '' 
  print 'sample with m_H = %d' % input.sample.mass
  print 'type = %s' % input.type
  print 'cat = %s' % input.region

  # initialize fit results
  fit_histo = input.histo
  fit_m = observable('m', 2, 0)
  fit_sigma = observable('sigma', 2, 0)
  fit_bulk = observable_range(0, 0)
  fit_frac = observable('frac', 0, 0)
  fit_fwhm = observable('fwhm', fwhm(input.histo), 0)
  fit_range = observable_range(input.sample.min, input.sample.max)

  # perform the fit
  if (input.histo.GetEntries() > 0):
    n_itr = 0
    tolerance = 0.00001
    prev_m = observable('previous_m', 1, 0)
    prev_sigma = observable('previous_sigma', 1, 0)

    while (n_itr < 100 and not (abs(get_distance(fit_m, prev_m)) < tolerance and abs(get_distance(fit_sigma, prev_sigma)) < tolerance)):
      n_itr = n_itr + 1

      prev_m.copyFrom(fit_m)
      prev_sigma.copyFrom(fit_sigma)

      fitres = fit_histo.Fit('gaus', 'SR', '', fit_range.min, fit_range.max)

      fitfunction = fit_histo.GetFunction('gaus')

      fit_m.val     = fitfunction.GetParameter(1)
      fit_m.err     = fitfunction.GetParError(1)
      fit_sigma.val = fitfunction.GetParameter(2)
      fit_sigma.err = fitfunction.GetParError(2)

      if (input.type == '4mu'):
        fit_range.min = fit_m.val - 2.0 * fit_sigma.val
        fit_range.max = fit_m.val + 2.0 * fit_sigma.val
      else:
        fit_range.min = fit_m.val - 1.5 * fit_sigma.val
        fit_range.max = fit_m.val + 2.5 * fit_sigma.val
        
    print 'did %d iterations' % n_itr
    print fit_m
    print fit_sigma

    # compute the fraction
    fit_bulk   = observable_range(fit_m.val - 2.0 * fit_sigma.val, fit_m.val + 2.0 * fit_sigma.val)
    left_part  = fit_histo.Integral(0, fit_histo.FindBin(fit_bulk.min))
    right_part = fit_histo.Integral(fit_histo.FindBin(fit_bulk.max), fit_histo.GetNbinsX() + 1)
    fit_frac.val = (left_part + right_part) / fit_histo.Integral(0, fit_histo.GetNbinsX() + 1)
    fit_frac.err = 0
  else:
    print 'empty histogram, unable to fit it...'

  # print the canvas
  c1 = TCanvas(fit_histo.GetName() + '_fit', m4l_name[input.type], 600, 600)
  c1.cd()
  fit_histo.Draw()
  fit_histo.GetXaxis().SetTitle(m4l_name[input.type] + ' [GeV]')
  fit_histo.GetXaxis().SetNdivisions(505);
  fit_histo.GetYaxis().SetTitle('a.u. / 0.5 GeV')
  fit_histo.GetYaxis().SetRangeUser(0, fit_histo.GetMaximum()*1.2)
  fit_histo.GetYaxis().SetNdivisions(505);
  fitfunction = fit_histo.GetFunction('gaus')
  if (fitfunction):
    fitfunction.SetLineColor(kRed)
    fitfunction.Draw("same")

  c1.cd()
  ATLASLabel(.2,.85, 'Simulation Internal')
  region_label = TLatex()
  region_label.SetName(fit_histo.GetName() + '_region_label')
  region_label.SetNDC()
  region_label.SetTextSize(0.033)
  region_label.SetText(.2, .6, input.region)
  region_label.Draw()
  mass_label = TLatex()
  mass_label.SetName(fit_histo.GetName() + '_mass_label')
  mass_label.SetNDC()
  mass_label.SetTextSize(0.033)
  mass_label.SetText(0.2, 0.5, "m = (%.2f #pm %.2f) GeV" % (fit_m.val, fit_m.err))
  mass_label.Draw()
  sigma_label = TLatex()
  sigma_label.SetName(fit_histo.GetName() + '_sigma_label')
  sigma_label.SetNDC()
  sigma_label.SetTextSize(0.033)
  sigma_label.SetText(0.2, 0.5-.05, "#sigma = (%.2f #pm %.2f) GeV" % (fit_sigma.val, fit_sigma.err))
  sigma_label.Draw()

  plot_legend = TLegend()
  plot_legend.SetName(fit_histo.GetName() + '_plot_legend')
  plot_legend.SetX1NDC(0.2)
  plot_legend.SetY1NDC(0.7)
  plot_legend.SetX2NDC(0.65)
  plot_legend.SetY2NDC(0.78)
  plot_legend.SetTextSize(0.033)
  plot_legend.AddEntry(fit_histo, "m_{H} = " + str(sample.mass) + " GeV", "p")
  plot_legend.AddEntry(fitfunction, "gaussian fit", "l")
  plot_legend.SetBorderSize(0)
  plot_legend.SetFillColor(0)
  plot_legend.Draw()

  placeholder_frac = TLatex()
  placeholder_frac.SetName(fit_histo.GetName() + '_placeholder_frac')
  placeholder_frac.SetNDC()
  placeholder_frac.SetTextSize(0.033)
  placeholder_frac.SetText(0.2, 0.5-.1, "fraction outside #pm 2#sigma: %.2f" % fit_frac.val)
  placeholder_frac.Draw()

  c1.Modified()
  c1.Update()

  # save it here (ROOT bug?)
  c1.Write()

  if (doSeparateOutputs):
    c1.SaveAs(c1.GetName() + '.root')

  # return results
  result = fitResult(c1, fit_m, fit_sigma, fit_bulk, fit_frac, fit_fwhm)
  return result

# TRUTH fitter
def truth_fit(input):
   # define observables and models
   if (cat != 'all'):
     if (input.sample.mass < 360):
       min = -30
       max = 15
     else:
       min = -30
       max = 30
   else:
     if (input.sample.mass < 130):
       min = -10
       max = 5
     elif (input.sample.mass < 180):
       min = -10
       max = 7
     elif (input.sample.mass < 360):
       min = -20
       max = 10
     elif (input.sample.mass < 460):
       min = -30
       max = 20
     else:
       min = -30
       max = 30

   m4l = RooRealVar("m4l", m4l_name[input.type] + ' - m_{H}^{(truth)}', min, max, "GeV");
   data = RooDataHist("data", "data", RooArgList(m4l), input.histo);

   cb_m0 = RooRealVar("cb_m0", "cb_m0", -4, 4, "GeV");
   cb_sigma = RooRealVar("cb_sigma", "cb_sigma", 0, 30, "GeV");
   cb_alpha = RooRealVar("cb_alpha", "cb_alpha", 0, 3);
   cb_n = RooRealVar("cb_n", "cb_n", 0, 200);
   cb = RooCBShape("cb", "cb", m4l, cb_m0, cb_sigma, cb_alpha, cb_n);

  #intrinsic_mean = RooRealVar("intrinsic_mean", "intrinsic_mean", sample.mass, sample.mass, sample.mass, "GeV")
  #intrinsic_width = RooRealVar("intrinsic_width", "intrinsic_width", higgs_width[sample.mass], higgs_width[sample.mass], higgs_width[sample.mass], "GeV")
  #intrinsic = RooBreitWigner("intrinsic", "intrinsic", m4l, intrinsic_mean, intrinsic_width);

   # set the signal model
   signal = cb

   # perform the fit
   if (input.histo.GetEntries() > 0):
     signal.fitTo(data, RooFit.SumW2Error(True))

   # graphic output
   c2 = TCanvas(input.histo.GetName() + '_truthfit', m4l_name[input.type], 600, 600)
   c2.cd()
   frame = m4l.frame(RooFit.Name(input.histo.GetName() + '_truthfitframe'), RooFit.Title(''))
   data.plotOn(frame)
   if (input.histo.GetEntries() > 0):
     signal.plotOn(frame, RooFit.LineColor(kRed))
     signal.paramOn(frame)
   frame.Draw()


   # obtain fwhm and uncertainty
   fit_fwhm = observable('truth_fwhm', 0, 0)

   if (input.histo.GetEntries() > 0 and doTruthFit):
     justcb = cb

     gen_dataset = justcb.generate(RooArgSet(m4l), 1000000)
     gen_hist = TH1D('tmp_datahist', '', 1000, min, max)
     gen_dataset.fillHistogram(gen_hist,  RooArgList(m4l))

     fit_fwhm.val = fwhm(gen_hist)

     del gen_dataset
     del gen_hist

     # generate 100 toys with the same number of entries as the input dataset
     # to get a reasonable uncertainty on this number
     sum_fwhm = 0
     sum_fwhm_sq = 0

     n_toys = 100
     n_events = int(input.histo.GetEntries())
     for i in range(n_toys):
       gen_dataset = justcb.generate(RooArgSet(m4l), n_events)
       gen_hist = TH1D('tmp_datahist', '', 1000, min, max)
       gen_dataset.fillHistogram(gen_hist,  RooArgList(m4l))

       this_fwhm = fwhm(gen_hist)
       sum_fwhm = sum_fwhm + this_fwhm
       sum_fwhm_sq = sum_fwhm_sq + this_fwhm * this_fwhm

       del gen_dataset
       del gen_hist

     fwhm_err = sqrt(sum_fwhm_sq/n_toys - pow(sum_fwhm/n_toys, 2))
     fit_fwhm.err = fwhm_err

   c2.Modified()
   c2.Update()

  #result = [c2, fit_fwhm]
   result = [frame, fit_fwhm]

   if (doSeparateOutputs):
    c2.SaveAs(c2.GetName() + '.root')

   del m4l
   del data
   del cb_m0
   del cb_sigma
   del cb_alpha
   del cb_n
   del cb
  #del intrinsic_mean
  #del intrinsic_width
  #del intrinsic

   return result

# compute the rescaled higgs mass
def get_rescaled_higgs_mass(tree, sample, type, cat):
  ptRescale = sample.mass / fit_result[sample.mass][type][cat].m.val

  if (doApplyMassConstraint):
    return tree.H_m_constrained * ptRescale

  v4_1 = TLorentzVector()
  v4_2 = TLorentzVector()
  v4_3 = TLorentzVector()
  v4_4 = TLorentzVector()
  v4_1.SetPtEtaPhiM(tree.Z1_lepplus_pt, tree.Z1_lepplus_eta, tree.Z1_lepplus_phi, tree.Z1_lepplus_m)
  v4_2.SetPtEtaPhiM(tree.Z1_lepminus_pt, tree.Z1_lepminus_eta, tree.Z1_lepminus_phi, tree.Z1_lepminus_m)
  v4_3.SetPtEtaPhiM(tree.Z2_lepplus_pt, tree.Z2_lepplus_eta, tree.Z2_lepplus_phi, tree.Z2_lepplus_m)
  v4_4.SetPtEtaPhiM(tree.Z2_lepminus_pt, tree.Z2_lepminus_eta, tree.Z2_lepminus_phi, tree.Z2_lepminus_m)
  pippo=(v4_1 + v4_2 + v4_3 + v4_4).M()
  v4_1.SetXYZM(ptRescale*v4_1.Px(), ptRescale*v4_1.Py(), ptRescale*v4_1.Pz(), v4_1.M())
  v4_2.SetXYZM(ptRescale*v4_2.Px(), ptRescale*v4_2.Py(), ptRescale*v4_2.Pz(), v4_2.M())
  v4_3.SetXYZM(ptRescale*v4_3.Px(), ptRescale*v4_3.Py(), ptRescale*v4_3.Pz(), v4_3.M())
  v4_4.SetXYZM(ptRescale*v4_4.Px(), ptRescale*v4_4.Py(), ptRescale*v4_4.Pz(), v4_4.M())
  res = (v4_1 + v4_2 + v4_3 + v4_4).M()
  if (res/1000 != sample.mass*ptRescale):
    pass
   #print 'fixme: ptRescale=%f recomputedOld=%f oldMass=%f expectedNew=%f newMass=%f' % (ptRescale, pippo, tree.H_m, tree.H_m * ptRescale, res)
  return res


# fills histograms and counters
def fill_histo_counters(target_histo, target_histo_truth, doPtRescale):
  # loop over samples
  for sample in samples:
    print 'looking for samples with mH = ' + str(sample.mass) + ' GeV'
  
    # cross section weighting stuff
    gg_xsec_weight = 0.
  
    # loop over the merging of gg and vbf files
    for filename in mergelists(sample.list_gg, sample.list_vbf):
      print '  opening ' + filename
      file = TFile.Open(filename)
      tree = file.Get("candidates")
      
      for entry in range(tree.GetEntries()):
        tree.GetEntry(entry)
      
        if (tree.selected != 1):
     	  continue
  
	if merge_2e2mu:
          if (tree.type > 1):
            type = types[tree.type - 1]
          else:
            type = types[tree.type]
	else:
          type = types[tree.type]
      
        # compute overall weight
        weight  = 1.0
        weight *= tree.pu_weight
        if filename in sample.list_gg:
  	  gg_xsec_weight = tree.xsec_weight
        else:
  	  if gg_xsec_weight == 0.:
  	    print 'WARNING: gluon-gluon fusion sample not yet retrieved, ABORT!!!'
  	  print 'ERROR: reweight not yet implemented, PLEASE do it! (use an entries map)'
  	weight *= tree.xsec_weight / gg_xsec_weight 
	if tag == 'mc11c':
          weight *= tree.ggF_weight
	elif tag == 'mc12a':
          weight *= tree.vxz_weight
        weight *= tree.Z1_lepplus_weight
        weight *= tree.Z1_lepminus_weight
        weight *= tree.Z2_lepplus_weight
        weight *= tree.Z2_lepminus_weight
        weight *= tree.trigSF_weight

	# find the event category (needed for scale correction and optional category study)
        cat = get_category(tree)

	# compute the higgs mass
	if (not doPtRescale):
	  if (not doApplyMassConstraint):
            my_4l_mass = tree.H_m / 1000.
            my_Z1_mass = tree.Z1_m / 1000.
            my_Z2_mass = tree.Z2_m / 1000.
	  else:
            my_4l_mass = tree.H_m_constrained / 1000.
            my_Z1_mass = tree.Z1_m_constrained / 1000.
            my_Z2_mass = tree.Z2_m_constrained / 1000.
	else:
	  my_4l_mass = get_rescaled_higgs_mass(tree, sample, type, cat) / 1000
          my_Z1_mass = tree.Z1_m_constrained / 1000.
          my_Z2_mass = tree.Z2_m_constrained / 1000.
     
        # fill the standard plot
        count_plain[sample.mass][type]['all'] += 1
        count_weighted[sample.mass][type]['all'] += weight
        target_histo[sample.mass][type]['all'].Fill(my_4l_mass, weight)
        target_histo_truth[sample.mass][type]['all'].Fill(my_4l_mass - tree.H_m_truth/1000, weight)
	# fill (ONCE ONLY) the reco-truth plots
	if (not doPtRescale):
          histo_mZ1[sample.mass][type]['all'].Fill(my_Z1_mass, weight)
          histo_mZ2[sample.mass][type]['all'].Fill(my_Z2_mass, weight)
	  histo_m4l_minus_truth[sample.mass][type]['all'].Fill(my_4l_mass - tree.H_m_truth/1000, weight)
	  histo_mZ1_minus_truth[sample.mass][type]['all'].Fill(my_Z1_mass - tree.Z1_m_truth/1000, weight)
	  histo_mZ2_minus_truth[sample.mass][type]['all'].Fill(my_Z2_mass - tree.Z2_m_truth/1000, weight)
  
        if (doCategories):
   	  # fill the right category
  	  count_plain[sample.mass][type][cat] += 1
  	  count_weighted[sample.mass][type][cat] += weight
  	  target_histo[sample.mass][type][cat].Fill(my_4l_mass, weight)
  	  target_histo_truth[sample.mass][type][cat].Fill(my_4l_mass - tree.H_m_truth/1000, weight)
          # fill (ONCE ONLY) the reco-truth plots
          if (not doPtRescale):
            histo_mZ1[sample.mass][type][cat].Fill(my_Z1_mass, weight)
            histo_mZ2[sample.mass][type][cat].Fill(my_Z2_mass, weight)
            histo_m4l_minus_truth[sample.mass][type][cat].Fill(my_4l_mass - tree.H_m_truth/1000, weight)
            histo_mZ1_minus_truth[sample.mass][type][cat].Fill(my_Z1_mass - tree.Z1_m_truth/1000, weight)
            histo_mZ2_minus_truth[sample.mass][type][cat].Fill(my_Z2_mass - tree.Z2_m_truth/1000, weight)
      
      file.Close()
      






####### THE ACTUAL SEQUENCE #########

# retrieve the histograms
fill_histo_counters(histo_m4l, histo_m4l_truth, False)

# perform the first fit
output_file.cd()

for sample in samples:
  for type in types:
    for cat in categories[type]:
      # normalize to unity if requested
      if (doNormalize):
        # ordinary histo
        area = histo_m4l[sample.mass][type][cat].Integral(0, histo_m4l[sample.mass][type][cat].GetNbinsX()+1)
	if (area != 0):
          histo_m4l[sample.mass][type][cat].Scale(1/area)

        # truth histo
        area = histo_m4l_truth[sample.mass][type][cat].Integral(0, histo_m4l_truth[sample.mass][type][cat].GetNbinsX()+1)
	if (area != 0):
          histo_m4l_truth[sample.mass][type][cat].Scale(1/area)

      # do the preliminary fit
      myInput = fitInput(sample, histo_m4l[sample.mass][type][cat], type, cat)
      fit_result[sample.mass][type][cat] = CSC_fit(myInput)

      # do the corresponding truth fit
      if (doTruthFitInCategories or cat == 'all'):
        myInput2 = fitInput(sample, histo_m4l_truth[sample.mass][type][cat], type, cat)
        fit_result_truth[sample.mass][type][cat] = truth_fit(myInput2)

      # fill the rescaling map
      print 'fixme: {m_H = %d GeV} [%-5s] (%-5s) %f' % (sample.mass, type, cat, sample.mass/fit_result[sample.mass][type][cat].m.val)

# retrieve the rescaled histograms
if (doRescaleToSampleMass):
  # reset stuff
  for sample in samples:
    for type in types:
      for cat in categories[type]:
        count_plain[sample.mass][type][cat] = 0.
        count_weighted[sample.mass][type][cat] = 0.
	# shouldn't these be different counters, in principle?

  # compute everything again after the mass rescale
  fill_histo_counters(histo_m4l_rescaled, histo_m4l_rescaled_truth, True)

  # perform the second fit
  output_file.cd()

  for sample in samples:
    for type in types:
      for cat in categories[type]:
        # normalize to unity if requested
        if (doNormalize):
	  # ordinary histo
          area = histo_m4l_rescaled[sample.mass][type][cat].Integral(0, histo_m4l_rescaled[sample.mass][type][cat].GetNbinsX()+1)
          if (area != 0):
            histo_m4l_rescaled[sample.mass][type][cat].Scale(1/area)

	  # truth histo
          area = histo_m4l_rescaled_truth[sample.mass][type][cat].Integral(0, histo_m4l_rescaled_truth[sample.mass][type][cat].GetNbinsX()+1)
          if (area != 0):
            histo_m4l_rescaled_truth[sample.mass][type][cat].Scale(1/area)

        # do the second fit
        myInput = fitInput(sample, histo_m4l_rescaled[sample.mass][type][cat], type, cat)
        fit_result_rescaled[sample.mass][type][cat] = CSC_fit(myInput)

        # do the corresponding truth fit
        if (doTruthFitInCategories or cat == 'all'):
          myInput2 = fitInput(sample, histo_m4l_rescaled_truth[sample.mass][type][cat], type, cat)
          fit_result_rescaled_truth[sample.mass][type][cat] = truth_fit(myInput2)

        # fill the rescaling map
        print 'fixme: {m_H = %d GeV} [%-5s] (%-5s) %f' % (sample.mass, type, cat, sample.mass/fit_result[sample.mass][type][cat].m.val)
        print 'fixme: {m_H = %d GeV} [%-5s] (%-5s) %f' % (sample.mass, type, cat, sample.mass/fit_result_rescaled[sample.mass][type][cat].m.val)

# print event summary
for sample in samples:
  print ''
  print 'm_H = %d GeV' % sample.mass
  for type in types:
    print '[%s]' % type
    for cat in categories[type]:
      print '  %-6s: %f / %f' % (cat, count_plain[sample.mass][type][cat], count_weighted[sample.mass][type][cat])

      # fill event counter
      if (count_weighted[sample.mass][type]['all'] != 0):
        population[sample.mass][type][cat] = count_weighted[sample.mass][type][cat] / count_weighted[sample.mass][type]['all']
      else:
        population[sample.mass][type][cat] = 0


# print detailed LaTeX output
print ''
print ''
print ''
print '*** CONFNOTE OUTPUT ***'
print ''

import re

confnote_descr = {}
for type in types:
  confnote_descr[type] = {}

confnote_descr['4mu']['all'] = 'all events'
confnote_descr['4mu']['bbbb'] = 'all muons in the barrel'
confnote_descr['4mu']['bbb'] = 'three muons in the barrel'
confnote_descr['4mu']['bb'] = 'two muons in the barrel'
confnote_descr['4mu']['other'] = 'any other event'
if not merge_2e2mu:
  confnote_descr['2mu2e']['all'] = 'all events'
  confnote_descr['2mu2e']['any_onecrk'] = 'at least one electron in the crack region'
  confnote_descr['2mu2e']['bb_bb'] = 'all leptons in the barrel'
  confnote_descr['2mu2e']['other_bb'] = 'electrons in the barrel, at least a muon in the endcap'
  confnote_descr['2mu2e']['other_other'] = 'any other event'
confnote_descr['2e2mu']['all'] = 'all events'
confnote_descr['2e2mu']['onecrk_any'] = 'at least one electron in the crack region'
confnote_descr['2e2mu']['bb_bb'] = 'all leptons in the barrel'
confnote_descr['2e2mu']['bb_other'] = 'electrons in the barrel, at least a muon in the endcap'
confnote_descr['2e2mu']['other_other'] = 'any other event'
confnote_descr['4e']['all'] = 'all events'
confnote_descr['4e']['bbbb'] = 'all electrons in the barrel'
confnote_descr['4e']['onecrk'] = 'at least one electron in the crack region'
confnote_descr['4e']['bbb'] = 'three electrons in the barrel (none in the crack)'
confnote_descr['4e']['other'] = 'any other event'


# default fit
print '[USUAL FIT]'
print ''
n_sample = -1
for sample in samples:
  n_sample = n_sample + 1
  print ''
  print '% MASS = ' + str(sample.mass)
  print '\\hline'
  print 'channel & name & description & frequency & $m$ [GeV] & $\\sigma$ [GeV] & events outside $\\pm2\\sigma$  & FWHM [GeV]\\\\\\hline\\hline'
  for type in types:
    for cat in categories[type]:
      mytype     = type_latex[type]
      mypop      = population[sample.mass][type][cat]
      mymassval  = fit_result[sample.mass][type][cat].m.val
      mymasserr  = fit_result[sample.mass][type][cat].m.err
      mysigmaval = fit_result[sample.mass][type][cat].sigma.val
      mysigmaerr = fit_result[sample.mass][type][cat].sigma.err
      mytails    = fit_result[sample.mass][type][cat].frac.val
      myfwhm     = fit_result[sample.mass][type][cat].fwhm.val

      cat_withoutunderscore = re.sub('_', '\_', cat)
      print '$%s$ & \\texttt{%s} & %s & $%.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$ & $%.2f$ & $%.1f$ \\\\\\hline' % (mytype, cat_withoutunderscore, confnote_descr[type][cat], mypop, mymassval, mymasserr, mysigmaval, mysigmaerr, mytails, myfwhm)

      histo_m4l[sample.mass][type][cat].SetDirectory(output_file)
      histo_mZ1[sample.mass][type][cat].SetDirectory(output_file)
      histo_mZ2[sample.mass][type][cat].SetDirectory(output_file)
      histo_m4l_minus_truth[sample.mass][type][cat].SetDirectory(output_file)
      histo_mZ1_minus_truth[sample.mass][type][cat].SetDirectory(output_file)
      histo_mZ2_minus_truth[sample.mass][type][cat].SetDirectory(output_file)
     #fit_result[sample.mass][type][cat].histo.Write()

      rescale_vs_mH[type][cat].SetPoint(n_sample, sample.mass, sample.mass/mymassval)
      rescale_vs_mH[type][cat].SetPointError(n_sample, 0, abs(sample.mass/(mymassval*mymassval)) * mymasserr)
      sigma_vs_mH[type][cat].SetPoint(n_sample, sample.mass, mysigmaval)
      sigma_vs_mH[type][cat].SetPointError(n_sample, 0, mysigmaerr)
      fwhm_vs_mH[type][cat].SetPoint(n_sample, sample.mass, myfwhm)
      if (doTruthFitInCategories or cat == 'all'):
        fit_result_truth[sample.mass][type][cat][0].Write()

        fwhm_vs_mH_truth[type][cat].SetPoint(n_sample, sample.mass, fit_result_truth[sample.mass][type][cat][1].val)
        fwhm_vs_mH_truth[type][cat].SetPointError(n_sample, 0, fit_result_truth[sample.mass][type][cat][1].err)

# after rescale
if (doRescaleToSampleMass):
  print ''
  print ''
  print '[AFTER pT RESCALE]'
  print ''

  n_sample = -1
  for sample in samples:
    n_sample = n_sample + 1
    print ''
    print '% MASS = ' + str(sample.mass)
    print '\\hline'
    print 'channel & name & description & frequency & $m$ [GeV] & $\\sigma$ [GeV] & events outside $\\pm2\\sigma$  & FWHM [GeV]\\\\\\hline\\hline'
    for type in types:
      for cat in categories[type]:
        mytype     = type_latex[type]
        mypop      = population[sample.mass][type][cat]
        mymassval  = fit_result_rescaled[sample.mass][type][cat].m.val
        mymasserr  = fit_result_rescaled[sample.mass][type][cat].m.err
        mysigmaval = fit_result_rescaled[sample.mass][type][cat].sigma.val
        mysigmaerr = fit_result_rescaled[sample.mass][type][cat].sigma.err
        mytails    = fit_result_rescaled[sample.mass][type][cat].frac.val
        myfwhm     = fit_result_rescaled[sample.mass][type][cat].fwhm.val
  
        cat_withoutunderscore = re.sub('_', '\_', cat)
        print '$%s$ & \\texttt{%s} & %s & $%.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$ & $%.2f$ & $%.1f$ \\\\\\hline' % (mytype, cat_withoutunderscore, confnote_descr[type][cat], mypop, mymassval, mymasserr, mysigmaval, mysigmaerr, mytails, myfwhm)
  
        histo_m4l_rescaled[sample.mass][type][cat].SetDirectory(output_file)
       #fit_result_rescaled[sample.mass][type][cat].histo.Write()

        rescale_vs_mH_rescaled[type][cat].SetPoint(n_sample, sample.mass, sample.mass/mymassval)
        rescale_vs_mH_rescaled[type][cat].SetPointError(n_sample, 0, abs(sample.mass/mymassval) * mymasserr)
        sigma_vs_mH_rescaled[type][cat].SetPoint(n_sample, sample.mass, mysigmaval)
        sigma_vs_mH_rescaled[type][cat].SetPointError(n_sample, 0, mysigmaerr)
        fwhm_vs_mH_rescaled[type][cat].SetPoint(n_sample, sample.mass, myfwhm)
        if (doTruthFitInCategories or cat == 'all'):
          fit_result_rescaled_truth[sample.mass][type][cat][0].Write()
	 
          fwhm_vs_mH_rescaled_truth[type][cat].SetPoint(n_sample, sample.mass, fit_result_rescaled_truth[sample.mass][type][cat][1].val)
          fwhm_vs_mH_rescaled_truth[type][cat].SetPointError(n_sample, 0, fit_result_rescaled_truth[sample.mass][type][cat][1].err)
  

# save the summary histograms
for type in types:
  for cat in categories[type]:
    rescale_vs_mH[type][cat].Write()
    sigma_vs_mH[type][cat].Write()
    fwhm_vs_mH[type][cat].Write()
    rescale_vs_mH_rescaled[type][cat].Write()
    sigma_vs_mH_rescaled[type][cat].Write()
    fwhm_vs_mH_rescaled[type][cat].Write()
    fwhm_vs_mH_truth[type][cat].Write()
    fwhm_vs_mH_rescaled_truth[type][cat].Write()


output_file.Write()
output_file.Close()
