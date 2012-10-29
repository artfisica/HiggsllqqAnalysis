# check the effect of using the d0 cut on all the leptons
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import os

the_dir='atlasscratchdisk/user/vippolit/H4l000009'

os.system('ls /tmp/vippolit/user*/*root* > __TMP_CUTFLOWMAKER__')

selected = {}
passing_additional = {}

for i in range(0, 4):
  selected[i] = {}
  passing_additional[i] = {}

for filename in open('__TMP_CUTFLOWMAKER__'):
  print 'opening %s' % filename
  f = TFile.Open(filename.rstrip('\n'))

  t = f.Get('candidates')

  if (not t):
    print 'tree not found here.'
    continue

  t.GetEntry(0)
  if t.run not in selected[0].keys():
    for i in range(0, 4):
      selected[i][t.run] = 0.
      passing_additional[i][t.run] = 0.

  for entry in range(t.GetEntries()):
    t.GetEntry(entry)

    poids = 1.#= t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.ggF_weight * t.top_weight * t.powhegbug_weight * t.vxz_weight
    if poids != poids:
      print '!!!'
      print t.Z1_lepplus_weight
      print t.Z1_lepminus_weight
      print t.Z2_lepplus_weight
      print t.Z2_lepminus_weight
      print t.ggF_weight
      print t.top_weight
      print t.powhegbug_weight
      print t.vxz_weight
    #poids = t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.ggF_weight * t.top_weight * t.powhegbug_weight * t.vxz_weight
   #print poids

   #poids = t.vxz_weight * t.pu_weight

    if (t.selected == 1):
      selected[t.type][t.run] += poids
      
      Z1_lepplus_d0sign_cut = 3.5 if t.Z1_lepplus_m > 100 else 6.5
      Z1_lepminus_d0sign_cut = 3.5 if t.Z1_lepminus_m > 100 else 6.5
      Z2_lepplus_d0sign_cut = 3.5 if t.Z2_lepplus_m > 100 else 6.5
      Z2_lepminus_d0sign_cut = 3.5 if t.Z2_lepminus_m > 100 else 6.5

      if (abs(t.Z1_lepplus_d0/t.Z1_lepplus_d0_sig) > Z1_lepplus_d0sign_cut) : continue
      if (abs(t.Z1_lepminus_d0/t.Z1_lepminus_d0_sig) > Z1_lepminus_d0sign_cut) : continue
      if (abs(t.Z2_lepplus_d0/t.Z2_lepplus_d0_sig) > Z2_lepplus_d0sign_cut) : continue
      if (abs(t.Z2_lepminus_d0/t.Z2_lepminus_d0_sig) > Z2_lepminus_d0sign_cut) : continue

      passing_additional[t.type][t.run] += poids

  f.Close()

chan_names = ['4mu', '2mu2e', '2e2mu', '4e']

for run in selected[0].keys():
  for chan in range(0, 4):
    print '%s [%s] %lf (%lf/%lf)' % (str(run), chan_names[chan], passing_additional[chan][run]/selected[chan][run], passing_additional[chan][run], selected[chan][run])
