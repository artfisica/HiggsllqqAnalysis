# imitate the SM analysis cuts and see what happens
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import array
import PyCintex
PyCintex.Cintex.Enable()
from ROOT import *

yield_h4l = {}
yield_sm = {}

channels = ['4mu', '2mu2e', '2e2mu', '4e']
regions = ['all', '<160', '>160']

for region in regions:
  yield_h4l[region] = {}
  yield_sm[region] = {}
  for chan in channels:
    yield_h4l[region][chan] = 0.
    yield_sm[region][chan] = 0.

lumi = {}
# SM lumi
#lumi['4mu'] = 2.1
#lumi['2mu2e'] = 2.1
#lumi['2e2mu'] = 2.1
#lumi['4e'] = 2.1
lumi['4mu'] = 3.2
lumi['2mu2e'] = 3.2
lumi['2e2mu'] = 3.2
lumi['4e'] = 3.2

xsec = {}
xsec[126937] = 69.75
xsec[126938] = 145.37
xsec[126939] = 103.06
xsec[126940] = 69.75
xsec[126941] = 103.06
xsec[126942] = 8.15

xsec[116601] = 0.570000
xsec[116602] = 0.570000
xsec[116603] = 1.140000

gen = {}
gen[126937] = 148659.516158
gen[126938] = 148440.078446
gen[126939] = 148547.713514
gen[126940] = 148613.073800
gen[126941] = 148730.803571
gen[126942] = 148719.853615

gen[116601] = 44509.572256
gen[116602] = 44375.866890
gen[116603] = 44396.331425
  
import glob

input_list_qq = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_nuove/*PowhegPythia8_AU2CT10_ZZ*/*root*')
input_list_gg = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_nuove/*gg2ZZ*/*root*')

input_list = [y for x in [input_list_qq, input_list_gg] for y in x]


for filename in input_list:
  f = TFile.Open(filename)
  print 'opened %s' % filename

  if not f: continue

  t = f.Get("candidates")

  if not t: continue

  for entry in range(t.GetEntries()):
    t.GetEntry(entry)

    # apply h4l cuts
    if (t.selected != 1): continue

    # remove the unused sample
    if (t.run == 116600): continue

    poids = t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.top_weight * t.powhegbug_weight * t.vxz_weight
    poids = poids * xsec[t.run] / gen[t.run] * lumi[channels[t.type]]

    yield_h4l['all'][channels[t.type]] += poids
    if (t.H_m < 160000): yield_h4l['<160'][channels[t.type]] += poids
    if (t.H_m >= 160000): yield_h4l['>160'][channels[t.type]] += poids


    # apply also SM cuts
    willBeRemoved = False

    if (t.Z1_lepplus_m < 100) and not (t.Z1_lepplus_pt > 15000 and (t.Z1_lepplus_z0) < 2 and abs(t.Z1_lepplus_d0/t.Z1_lepplus_d0_sig) < 6): willBeRemoved = True
    if (t.Z1_lepminus_m < 100) and not (t.Z1_lepminus_pt > 15000 and (t.Z1_lepminus_z0) < 2 and abs(t.Z1_lepminus_d0/t.Z1_lepminus_d0_sig) < 6): willBeRemoved = True
    if (t.Z2_lepplus_m < 100) and not (t.Z2_lepplus_pt > 15000 and (t.Z2_lepplus_z0) < 2 and abs(t.Z2_lepplus_d0/t.Z2_lepplus_d0_sig) < 6): willBeRemoved = True
    if (t.Z2_lepminus_m < 100) and not (t.Z2_lepminus_pt > 15000 and (t.Z2_lepminus_z0) < 2 and abs(t.Z2_lepminus_d0/t.Z2_lepminus_d0_sig) < 6): willBeRemoved = True

    if (t.Z1_lepplus_m > 100) and not (t.Z1_lepplus_author < 8 and t.Z1_lepplus_pt > 15000 and abs(t.Z1_lepplus_z0) < 2 and t.Z1_lepplus_isSA != 1): willBeRemoved = True
    if (t.Z1_lepminus_m > 100) and not (t.Z1_lepminus_author < 8 and t.Z1_lepminus_pt > 15000 and abs(t.Z1_lepminus_z0) < 2 and t.Z1_lepminus_isSA != 1): willBeRemoved = True
    if (t.Z2_lepplus_m > 100) and not (t.Z2_lepplus_author < 8 and t.Z2_lepplus_pt > 15000 and abs(t.Z2_lepplus_z0) < 2 and t.Z2_lepplus_isSA != 1): willBeRemoved = True
    if (t.Z2_lepminus_m > 100) and not (t.Z2_lepminus_author < 8 and t.Z2_lepminus_pt > 15000 and abs(t.Z2_lepminus_z0) < 2 and t.Z2_lepminus_isSA != 1): willBeRemoved = True

    if (not(t.Z1_m > 66000 and t.Z1_m < 116000)): willBeRemoved = True
    if (not(t.Z2_m > 66000 and t.Z2_m < 116000)): willBeRemoved = True

    if (willBeRemoved): continue

    # event has been accepted by SM

    yield_sm['all'][channels[t.type]] += poids
    if (t.H_m < 160000): yield_sm['<160'][channels[t.type]] += poids
    if (t.H_m >= 160000): yield_sm['>160'][channels[t.type]] += poids

  f.Close()


print ''
print ''

print 'h4l rates'
for region in regions:
  for chan in channels:
    print '[%s] (%s) %lf' % (region, chan, yield_h4l[region][chan])

print ''
print ''

print 'sm rates'
for region in regions:
  for chan in channels:
    print '[%s] (%s) %lf' % (region, chan, yield_sm[region][chan])
