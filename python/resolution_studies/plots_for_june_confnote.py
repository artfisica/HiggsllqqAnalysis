# produces some histograms involving truth and reco, to validate mass constraint
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import ROOT

input = 'output_constraint_mc12.root'
output = 'plots_constraintstudies.root'

f_out = ROOT.TFile(output, 'RECREATE')

channels = ['4mu', '2mu2e', '2e2mu', '4e']
histo_names_1D = [
                  'm12_minus_truth_over_truth',
                  'm12_constrained_minus_truth_over_truth',
	          'm4l_minus_truth_over_truth',
	          'm4l_constrained_minus_truth_over_truth',
	          'lep12_pt_minus_truth_over_truth'
	          'lep12_pt_constrained_minus_truth_over_truth'
		 ]
histo_names_2D = [
                  'm12_constrained_vs_m12',
	          'm4l_constrained_vs_m4l',
	         ]
histos = {}

for chan in channels:
  histos[chan] = {}

  # create
  histos[chan]['m12_minus_truth_over_truth'] = ROOT.TH1F('h_%s_m12_minus_truth_over_truth' % chan, ';m_{12}(reco)/m_{12}(truth) - 1;a.u. / 0.001', 600, -0.3, 0.3)
  histos[chan]['m12_constrained_minus_truth_over_truth'] = ROOT.TH1F('h_%s_m12_constrained_minus_truth_over_truth' % chan, ';m_{12}(constrained)/m_{12}(truth) - 1;a.u. / 0.001', 600, -0.3, 0.3)
  histos[chan]['m4l_minus_truth_over_truth'] = ROOT.TH1F('h_%s_m4l_minus_truth_over_truth' % chan, ';m_{4l}(reco)/m_{4l}(truth) - 1;a.u. / 0.001', 600, -0.3, 0.3)
  histos[chan]['m4l_constrained_minus_truth_over_truth'] = ROOT.TH1F('h_%s_m4l_constrained_minus_truth_over_truth' % chan, ';m_{4l}(constrained)/m_{4l}(truth) - 1;a.u. / 0.001', 600, -0.3, 0.3)
  histos[chan]['lep12_pt_minus_truth_over_truth'] = ROOT.TH1F('h_%s_lep12_pt_minus_truth_over_truth' % chan, ';lepton p_{T}(constrained)/lepton p_{T}(truth) - 1;a.u. / 0.001', 600, -0.3, 0.3)
  histos[chan]['lep12_pt_constrained_minus_truth_over_truth'] = ROOT.TH1F('h_%s_lep12_pt_constrained_minus_truth_over_truth' % chan, ';lepton p_{T}(constrained)/lepton p_{T}(truth) - 1;a.u. / 0.001', 600, -0.3, 0.3)

  histos[chan]['m12_constrained_vs_m12'] = ROOT.TH2F('h_%s_m12_constrained_vs_m12' % chan, ';m_{12}(reco) [GeV];m_{12}(constrained) [GeV]', 70, 50, 120, 70, 50, 120)
  histos[chan]['m4l_constrained_vs_m4l'] = ROOT.TH2F('h_%s_m4l_constrained_vs_m4l' % chan, ';m_{4l}(reco) [GeV];m_{4l}(constrained) [GeV]', 80, 100, 180, 80, 100, 180)


  # set directory
  histos[chan]['m12_minus_truth_over_truth'].SetDirectory(f_out)
  histos[chan]['m12_constrained_minus_truth_over_truth'].SetDirectory(f_out)
  histos[chan]['m4l_minus_truth_over_truth'].SetDirectory(f_out)
  histos[chan]['m4l_constrained_minus_truth_over_truth'].SetDirectory(f_out)
  histos[chan]['lep12_pt_minus_truth_over_truth'].SetDirectory(f_out)
  histos[chan]['lep12_pt_constrained_minus_truth_over_truth'].SetDirectory(f_out)

  histos[chan]['m12_constrained_vs_m12'].SetDirectory(f_out)
  histos[chan]['m4l_constrained_vs_m4l'].SetDirectory(f_out)

f = ROOT.TFile.Open(input)

t = f.Get('candidates')

if (t):
  for entry in xrange(t.GetEntries()):
    t.GetEntry(entry)

    if (t.selected == 1):

      # no cross-section weight for such a plot (we use only ggH130 :) )
      poids = t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.top_weight * t.powhegbug_weight * t.vxz_weight

      if (t.Z1_m_truth > 0):
        histos[channels[t.type]]['m12_minus_truth_over_truth'].Fill(t.Z1_m/t.Z1_m_truth - 1., poids)
        histos[channels[t.type]]['m12_constrained_minus_truth_over_truth'].Fill(t.Z1_m_constrained/t.Z1_m_truth - 1., poids)
        histos[channels[t.type]]['m4l_minus_truth_over_truth'].Fill(t.H_m/t.H_m_truth - 1., poids)
        histos[channels[t.type]]['m4l_constrained_minus_truth_over_truth'].Fill(t.H_m_constrained/t.H_m_truth - 1., poids)
      if (t.Z1_lepplus_pt_truth > 0):     
        histos[channels[t.type]]['lep12_pt_minus_truth_over_truth'].Fill(t.Z1_lepplus_pt/t.Z1_lepplus_pt_truth - 1., poids)
        histos[channels[t.type]]['lep12_pt_constrained_minus_truth_over_truth'].Fill(t.Z1_lepplus_pt_constrained/t.Z1_lepplus_pt_truth - 1., poids)
      if (t.Z1_lepminus_pt_truth > 0):     
        histos[channels[t.type]]['lep12_pt_minus_truth_over_truth'].Fill(t.Z1_lepminus_pt/t.Z1_lepminus_pt_truth - 1., poids)
        histos[channels[t.type]]['lep12_pt_constrained_minus_truth_over_truth'].Fill(t.Z1_lepminus_pt_constrained/t.Z1_lepminus_pt_truth - 1., poids)

      histos[channels[t.type]]['m12_constrained_vs_m12'].Fill(t.Z1_m/1000., t.Z1_m_constrained/1000., poids)
      histos[channels[t.type]]['m4l_constrained_vs_m4l'].Fill(t.H_m/1000.,  t.H_m_constrained/1000.,  poids)

f.Close()

f_out.Write()
f_out.Close()
