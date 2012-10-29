# makes the cutflow in a python way (we all love python <3), adding the deltar > 0.2 cut to evaluate impact
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import os

list_of_files = []
the_dir='atlaslocalgroupdisk/user/vippolit/H4l000009'

# not together, in the input, since mc_channel_number is filled badly for these validation samples -> will be fixed
#os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "H120.*topoetcone" > __TMP_CUTFLOWMAKER__' % the_dir)
#os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "H130.*topoetcone" > __TMP_CUTFLOWMAKER__' % the_dir)
#os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "Zee.*topoetcone" > __TMP_CUTFLOWMAKER__' % the_dir)
#os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "H130.*2012_campaign" > __TMP_CUTFLOWMAKER__' % the_dir)
os.system('ls output_test.ZZ_2e2mu.noFSRcut.root > __TMP_CUTFLOWMAKER__')
generated = {}
passing = {}
passing[0] = {}
passing[1] = {}
passing[2] = {}
passing[3] = {}

#for filename in list_of_files:
for filename in open('__TMP_CUTFLOWMAKER__'):
  print 'opening %s' % filename
 #f = TFile.Open(filename)
  f = TFile.Open(filename[:-1])

  e = f.Get('cutflow')
  e.GetEntry(0)

  if e.run not in generated.keys():
    generated[e.run] = 0.
    for chan in range(0, 4):
      passing[chan][e.run] = {}
      for i in range(12):
        passing[chan][e.run][i] = 0.

  generated[e.run] = generated[e.run] + e.GetEntries('analysis==0')

  c = f.Get('candidates')

  for chan in range(0, 4):
    for i in range(11):
      passing[chan][e.run][i] = passing[chan][e.run][i] + c.GetEntries('type==%d&&last>=%d&&closesttozz' % (chan, i))

  for entry in range(c.GetEntries()):
    c.GetEntry(entry)

    if not (c.selected == 1):
      continue

    tlv_Z1_lepplus = TLorentzVector(); tlv_Z1_lepplus.SetPtEtaPhiM(c.Z1_lepplus_pt, c.Z1_lepplus_eta, c.Z1_lepplus_phi, c.Z1_lepplus_m)
    tlv_Z1_lepminus = TLorentzVector(); tlv_Z1_lepminus.SetPtEtaPhiM(c.Z1_lepminus_pt, c.Z1_lepminus_eta, c.Z1_lepminus_phi, c.Z1_lepminus_m)
    tlv_Z2_lepplus = TLorentzVector(); tlv_Z2_lepplus.SetPtEtaPhiM(c.Z2_lepplus_pt, c.Z2_lepplus_eta, c.Z2_lepplus_phi, c.Z2_lepplus_m)
    tlv_Z2_lepminus = TLorentzVector(); tlv_Z2_lepminus.SetPtEtaPhiM(c.Z2_lepminus_pt, c.Z2_lepminus_eta, c.Z2_lepminus_phi, c.Z2_lepminus_m)
    tlv = [tlv_Z1_lepplus,tlv_Z1_lepminus,tlv_Z2_lepplus,tlv_Z2_lepminus]

    howmany = 1
    for i in range(0, 4):
      for j in range(i+1, 4):
        if (tlv[i].DeltaR(tlv[j]) < 0.2): howmany = 0
    passing[c.type][c.run][11] = passing[c.type][c.run][11] + howmany

  f.Close()

chan_names = ['4mu', '2mu2e', '2e2mu', '4e']
cut_names = ['mass_2e2mu', 'opposite_charge', 'kinematics', 'trigmatch', 'Z1_mass', 'Z2_mass', 'DeltaR', 'best', 'track_iso', 'calo_iso', 'd0_sig', 'mass_window', 'dr_less_02']

for run in generated.keys():
  print 'generated[' + str(run) + '] = ' + str(str(generated[run]))
  for chan in range(0, 4):
    for i in range(12):
      print '%s [%s] %s = %s' % (str(run), chan_names[chan], cut_names[i], str(passing[chan][run][i]))
