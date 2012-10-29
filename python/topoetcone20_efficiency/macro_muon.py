# produce final plots comparing certain sets of muon isolation variables
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *

ROOT.gROOT.LoadMacro("AtlasStyle.C") 
ROOT.gROOT.LoadMacro("AtlasLabels.C")
SetAtlasStyle()


class plotoptions:
  def __init__(self, logy, xtitle, ytitle, xmin = -9999, xmax = -9999):
    self.logy = logy
    self.xtitle = xtitle
    self.ytitle = ytitle
    self.xmin = xmin
    self.xmax = xmax

colors = [kBlack, kBlue+1, kRed+1, kGreen-8, kViolet, kOrange-3, kAzure+2]
markers = [20, 23, 25, 27, 24, 26]


compare_sets = [['etcone20', 'etcone20_final'], ['ptcone20_final']]

plots = ['cut_scan_%s_pt', 'histo_%s', 'histo_%s_vs_nvx_pfx']
plotopts = [plotoptions(0, 'muon isolation/p_{T} cut', 'efficiency', 0, 0.5), plotoptions(1, 'muon isolation [MeV]', 'entries'), plotoptions(0, 'number of vertices with #geq 2 tracks', 'average muon isolation [MeV]', 0, 30)]


filename = 'output_muon_topoiso.root'

f = TFile(filename)

canvases = {}

for plot in plots:
  for myset in compare_sets:
    canv_name = '%s_%d' % (plot, compare_sets.index(myset))
    canvases[canv_name] = TCanvas(canv_name, canv_name)
    canvases[canv_name].SetLogy(plotopts[plots.index(plot)].logy)

    for obs in myset:
      plotname = plot % obs
      the_plot = f.Get(plotname)

      if (the_plot):
	# draw option: prefix
	drawopt = ''
        if (myset.index(obs) == 0):
          if (the_plot.ClassName() != 'TGraph'):
	    drawopt = '' 
          elif (the_plot.ClassName() == 'TGraph'):
	    drawopt = 'a'

	# draw option: middle
        if (the_plot.ClassName() == 'TGraph'):
	  drawopt = drawopt + 'pl'
        else:
	  drawopt = drawopt + ''

	# draw option: end
        if (myset.index(obs) == 0):
	  drawopt = drawopt + ''
        else:
	  drawopt = drawopt + 'same'

        the_plot.SetMarkerStyle(markers[myset.index(obs)])
	if (the_plot.ClassName() == 'TGraph'):
          the_plot.SetMarkerColor(colors[myset.index(obs)])
	if (the_plot.ClassName() == 'TGraph' or the_plot.ClassName() == 'TH1F'):
          the_plot.SetLineColor(colors[myset.index(obs)])

	the_plot.Draw(drawopt)

	if (myset.index(obs) == 0):
	  the_plot.GetXaxis().SetTitle(plotopts[plots.index(plot)].xtitle)
	  the_plot.GetYaxis().SetTitle(plotopts[plots.index(plot)].ytitle)

	  if (plotopts[plots.index(plot)].xmin != -9999 or plotopts[plots.index(plot)].xmax != -9999):
	    the_plot.GetXaxis().SetRangeUser(plotopts[plots.index(plot)].xmin, plotopts[plots.index(plot)].xmax)

	print '%s drawn over %s with option %s' % (plotname, canv_name, drawopt)
      else:
        print 'unable to retrieve %s from %s' % (plotname, f.GetName())

    leg = canvases[canv_name].BuildLegend()
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
