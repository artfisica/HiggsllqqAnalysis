# makes the cutflow in a python way (we all love python <3)
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import os

the_dir='atlasscratchdisk/user/vippolit/H4l000009'

#list_of_files = os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "ggH120|Alpgen|ggH200" > __TMP_CUTFLOWMAKER__' % the_dir)
#list_of_files = os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "campaign" > __TMP_CUTFLOWMAKER__' % the_dir)
list_of_files = os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2012_nuove/*/*root* | grep -E "PowhegPythia8_AU2CT10_ZZ" > __TMP_CUTFLOWMAKER__')
#list_of_files = os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2011/*/*root* | grep -E "Pythiazz4l_3MultiLeptonFilterElecM" > __TMP_CUTFLOWMAKER__')

generated = {}
passing = {}
passing[0] = {}
passing[1] = {}
passing[2] = {}
passing[3] = {}

#poids = 'pu_weight * trigSF_weight * Z1_lepplus_weight * Z1_lepminus_weight * Z2_lepplus_weight * Z2_lepminus_weight * top_weight * powhegbug_weight * vxz_weight'
poids = '1'#pu_weight * trigSF_weight * Z1_lepplus_weight * Z1_lepminus_weight * Z2_lepplus_weight * Z2_lepminus_weight * top_weight * powhegbug_weight * vxz_weight'

for filename in open('__TMP_CUTFLOWMAKER__'):
  print 'opening %s' % filename
  f = TFile.Open(filename[:-1])

  e = f.Get('cutflow')
  h = f.Get('generatedEntriesHisto')
  e.GetEntry(0)

  if e.run not in generated.keys():
    generated[e.run] = 0.
    for chan in range(0, 4):
      passing[chan][e.run] = {}
      for i in range(12):
        passing[chan][e.run][i] = 0.

 #generated[e.run] = generated[e.run] + e.GetEntries('analysis==0')
  generated[e.run] += h.GetBinContent(5) # mc12!!!

  c = f.Get('candidates')

  for chan in range(0, 4):
    for i in range(11):
     #c.Draw('H_m>>htemp', '(type==%d&&last>=%d&&H_m<160000)*(%s)' % (chan, i, poids), 'goff')
      c.Draw('H_m>>htemp', '(type==%d&&last>=%d&&closesttozz&&H_m<160000)*(%s)' % (chan, i, poids), 'goff')
     #c.Draw('H_m>>htemp', '(type==%d&&last>=%d&&H_m<160000&&(Z1_m_truth>0&&Z2_m_truth>0))*(%s)' % (chan, i, poids), 'goff')
      tmp_histo = gROOT.FindObject('htemp')
      this_contribution = tmp_histo.Integral(0, tmp_histo.GetNbinsX()+1)
     #passing[chan][e.run][i] = passing[chan][e.run][i] + c.GetEntries('(type==%d&&last>=%d&&closesttozz)*(%s)' % (chan, i, poids))
      passing[chan][e.run][i] += this_contribution
    i = 11
   #c.Draw('H_m>>htemp', '(type==%d&&last>=%d&&H_m<160000&&abs(Z1_m-91000)<15000)*(%s)' % (chan, i-1, poids), 'goff')
    c.Draw('H_m>>htemp', '(type==%d&&last>=%d&&closesttozz&&H_m<160000&&abs(Z1_m-91000)<15000)*(%s)' % (chan, i-1, poids), 'goff')
   #c.Draw('H_m>>htemp', '(type==%d&&last>=%d&&H_m<160000&&abs(Z1_m-91000)<15000&&(Z1_m_truth>0&&Z2_m_truth>0))*(%s)' % (chan, i-1, poids), 'goff')
    tmp_histo = gROOT.FindObject('htemp')
    this_contribution = tmp_histo.Integral(0, tmp_histo.GetNbinsX()+1)
   #passing[chan][e.run][i] = passing[chan][e.run][i] + c.GetEntries('(type==%d&&last>=%d&&closesttozz)*(%s)' % (chan, i, poids))
    passing[chan][e.run][i] += this_contribution

  f.Close()

chan_names = ['4mu', '2mu2e', '2e2mu', '4e']
cut_names = ['mass_2e2mu', 'opposite_charge', 'kinematics', 'trigmatch', 'Z1_mass', 'Z2_mass', 'DeltaR', 'best', 'track_iso', 'calo_iso', 'd0_sig', 'mass_window']

for run in generated.keys():
  print 'generated[' + str(run) + '] = ' + str(str(generated[run]))
  for chan in range(0, 4):
    for i in range(12):
      if (i == 0): rel_eff = "NA"
      elif (passing[chan][run][i-1] == 0): rel_eff = "NaN"
      else: rel_eff = str(passing[chan][run][i]/passing[chan][run][i-1])

      print '%s [%s] %s = %s (%s rel eff) (%s abs eff)' % (str(run), chan_names[chan], cut_names[i], str(passing[chan][run][i]), rel_eff, str(passing[chan][run][i]/generated[run]))
