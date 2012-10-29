# efficiency of electron isolation cuts [0.15 for ptcone20, 0.20-0.25-0.30 for etcone20 and topoetcone20] vs nvx, pt, eta, phi
# for higgs, Z->mumu (divided in truth categories, i.e. electron-fake-conversion)
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *

ROOT.gROOT.LoadMacro("AtlasStyle.C") 
ROOT.gROOT.LoadMacro("AtlasLabels.C")
SetAtlasStyle()

breakdown_bkg = False

file_sig = TFile('output_electron_topoiso_ggH125.root')
file_bkg = TFile('output_electron_topoiso_Zmumu.root')
file_bkg_conversions = TFile('output_electron_topoiso_Zmumu_conversions.root')
file_bkg_fakes = TFile('output_electron_topoiso_Zmumu_fakes.root')

interesting_cuts = ['02', '025', '03']
interesting_x = ['nvx', 'pt', 'eta', 'phi']
interesting_x_label = ['n_{VX}', 'p_{T} [GeV]', '#eta', '#phi']
interesting_y = ['topoetcone20']#'etcone20', 'topoetcone20']

canvas = {}
leg = {}

for y in interesting_y:
  for x in interesting_x:
    eff_surname = 'eff_%s_vs_%s' % (y, x)
    canvas[eff_surname] = TCanvas('c_%s' % eff_surname, eff_surname)

    for i, cut in enumerate(interesting_cuts):
      eff_name = 'eff_%s_lt_%s_vs_%s' % (y, cut, x)

      h_sig = file_sig.Get(eff_name)
      h_bkg = file_bkg.Get(eff_name)
      h_bkg_conversions = file_bkg_conversions.Get(eff_name)
      h_bkg_fakes = file_bkg_fakes.Get(eff_name)

      h_sig.SetDirectory(0)
      h_bkg.SetDirectory(0)
      h_bkg_conversions.SetDirectory(0)
      h_bkg_fakes.SetDirectory(0)



      h_sig.SetLineColor(kRed + i)
      h_bkg.SetLineColor(kViolet + i)
      h_bkg_conversions.SetLineColor(kGreen-8 + i)
      h_bkg_fakes.SetLineColor(kBlue+1 + i)

      h_sig.SetLineWidth(2)
      h_bkg.SetLineWidth(2)
      h_bkg_conversions.SetLineWidth(1)
      h_bkg_fakes.SetLineWidth(1)

      h_sig.SetFillColor(0)
      h_bkg.SetFillColor(0)
      h_bkg_conversions.SetFillColor(0)
      h_bkg_fakes.SetFillColor(0)

      h_sig.SetMarkerSize(0)
      h_bkg.SetMarkerSize(0)
      h_bkg_conversions.SetMarkerSize(0)
      h_bkg_fakes.SetMarkerSize(0)

      h_sig.Draw('al' if i == 0 else 'lsame')
      h_bkg.Draw('lsame')
      if breakdown_bkg:
        h_bkg_conversions.Draw('lsame')
        h_bkg_fakes.Draw('lsame')

      if (i == 0):
        canvas[eff_surname].Update()
	h_sig.GetPaintedGraph().GetXaxis().SetTitle('electron %s' % interesting_x_label[interesting_x.index(x)])
	h_sig.GetPaintedGraph().GetYaxis().SetTitle('efficiency')

      h_sig.SetTitle('%s, %s/pt < 0.%s' % ('sig', y, cut[1:]))
      h_bkg.SetTitle('%s, %s/pt < 0.%s' % ('conv+had', y, cut[1:]))
      h_bkg_conversions.SetTitle('%s, %s/pt < 0.%s' % ('conv', y, cut[1:]))
      h_bkg_fakes.SetTitle('%s, %s/pt < 0.%s' % ('had', y, cut[1:]))

    leg[eff_surname] = canvas[eff_surname].BuildLegend()
    leg[eff_surname].SetBorderSize(0)
    leg[eff_surname].SetFillColor(0)

file_sig.Close()
file_bkg.Close()
file_bkg_conversions.Close()
file_bkg_fakes.Close()
