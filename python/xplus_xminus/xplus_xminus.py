from ROOT import *
import glob
import re

gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro("AtlasStyle.C") 
ROOT.gROOT.LoadMacro("AtlasLabels.C")
SetAtlasStyle()

plotdir = './plots'

samples = glob.glob('*/*root*')
samples.sort()
names   = []

#user.vippolit.H4l000000.mc11_7TeV.116765.PowHegPythia_ggH135_ZZ4lep.merge.p869.xplus_xminus

regexp_content = 'user\.(\w+)\.(\w+)\.(\w+)\.(\w+)\.(\w+)\..*'

regexp = re.compile(regexp_content)


colors = [
           kRed+1,
	   kBlue+1,
	   kGreen-8,
	   kViolet,
	   kBlack,
	   kOrange-3,
	   kAzure+2,
         ]

histograms = [
               'xplus_reco', 'xminus_reco',
               'xplus_truth', 'xminus_truth',
             ]
histo = {}

for h in histograms:
  histo[h] = {}


outfile = TFile('outfile_x.root', 'RECREATE')


for sample in samples:
  regexp_result = regexp.search(sample)

  if (regexp_result):
    _user = regexp_result.group(1)
    _code = regexp_result.group(2)
    _type = regexp_result.group(3)
    _runnumber = regexp_result.group(4)
    _name = regexp_result.group(5)

    if (_type == 'mc11_7TeV'):
      print 'processing %s' % (_name)

      outfile.cd()
      histo['xplus_reco'][_name] = TH1D('xplus_reco_%s' % _name, 'xplus_reco_%s;reco x_{+};arbitrary units / 0.02' % _name, 200, -2, 2)
      histo['xminus_reco'][_name] = TH1D('xminus_reco_%s' % _name, 'xminus_reco_%s;reco x_{-};arbitrary units / 0.02' % _name, 200, -2, 2)
      histo['xplus_truth'][_name] = TH1D('xplus_truth_%s' % _name, 'xplus_truth_%s;truth (reco matched) x_{+};arbitrary units / 0.02' % _name, 200, -2, 2)
      histo['xminus_truth'][_name] = TH1D('xminus_truth_%s' % _name, 'xminus_truth_%s;truth (reco matched) x_{-};arbitrary units / 0.02' % _name, 200, -2, 2)

      f = TFile(sample)

      if (f):
        t = f.Get('candidates')
        
	if (t):
	  for entry in range(t.GetEntries()):
	    t.GetEntry(entry)

	    if (t.selected == 1):
	      xplus_reco = 1. - (t.Z1_m + t.Z2_m) / t.H_m
	      xminus_reco = 1. - (t.Z1_m - t.Z2_m) / t.H_m
	      if (t.H_m_truth > 0 and t.Z1_m_truth > 0 and t.Z2_m_truth > 0):
                xplus_truth = 1. - (t.Z1_m_truth + t.Z2_m_truth) / t.H_m_truth
	        xminus_truth = 1. - (t.Z1_m_truth - t.Z2_m_truth) / t.H_m_truth
	      else:
                xplus_truth = -9999.9
	        xminus_truth = -9999.9

	     #weight = t.xsec_weight / t.processed_entries * t.pu_weight * t.ggF_weight
	     #weight = t.xsec_weight / t.processed_entries * t.pu_weight
	     #weight = t.xsec_weight / t.processed_entries * t.ggF_weight
	      weight = 1.

	      histo['xplus_reco'][_name].Fill(xplus_reco, weight)
	      histo['xminus_reco'][_name].Fill(xminus_reco, weight)
	      histo['xplus_truth'][_name].Fill(xplus_truth, weight)
	      histo['xminus_truth'][_name].Fill(xminus_truth, weight)
	else:
          print 'troubles opening candidates\' tree in %s' % (sample)

        f.Close()
      else:
        print 'troubles opening %s' % (sample)

outfile.Write()

for h in histograms:
  c = TCanvas('c_%s' % h, 'c_%s' % h)
  for sample in histo[h].keys():
    i = histo[h].keys().index(sample)

    histo[h][sample].SetLineColor(colors[i])
    histo[h][sample].SetMarkerColor(colors[i])

    if (histo[h][sample].Integral() > 0):
      if (i == 0):
        histo[h][sample].DrawNormalized()
      else:
        histo[h][sample].DrawNormalized('same')

  c.BuildLegend()
  c.SetLogy()
  outfile.cd()
  c.Write()#'%s/plot_%s.png' % (plotdir, h))
    
outfile.Close()
