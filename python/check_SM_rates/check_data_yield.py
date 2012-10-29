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
lumi['4mu'] = 2.1
lumi['2mu2e'] = 2.1
lumi['2e2mu'] = 2.1
lumi['4e'] = 2.1

import glob

input_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012/*physics*/*root*')


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

    poids = 1. # it is data!

    yield_h4l['all'][channels[t.type]] += poids
    if (t.H_m < 160000): yield_h4l['<160'][channels[t.type]] += poids
    if (t.H_m >= 160000): yield_h4l['>160'][channels[t.type]] += poids


    # apply also SM cuts
    willBeRemoved = False

    if (t.Z1_lepplus_m < 100) and not (t.Z1_lepplus_pt > 15000 and (t.Z1_lepplus_z0) < 2 and abs(t.Z1_lepplus_d0/t.Z1_lepplus_d0_sig) < 6): willBeRemoved = True
    if (t.Z1_lepminus_m < 100) and not (t.Z1_lepminus_pt > 15000 and (t.Z1_lepminus_z0) < 2 and abs(t.Z1_lepminus_d0/t.Z1_lepminus_d0_sig) < 6): willBeRemoved = True
    if (t.Z2_lepplus_m < 100) and not (t.Z2_lepplus_pt > 15000 and (t.Z2_lepplus_z0) < 2 and abs(t.Z2_lepplus_d0/t.Z2_lepplus_d0_sig) < 6): willBeRemoved = True
    if (t.Z2_lepminus_m < 100) and not (t.Z2_lepminus_pt > 15000 and (t.Z2_lepminus_z0) < 2 and abs(t.Z2_lepminus_d0/t.Z2_lepminus_d0_sig) < 6): willBeRemoved = True

    if (t.Z1_lepplus_m > 100) and not (t.Z1_lepplus_author != 18 and t.Z1_lepplus_pt > 15000 and abs(t.Z1_lepplus_z0) < 2 and t.Z1_lepplus_isSA != 1): willBeRemoved = True
    if (t.Z1_lepminus_m > 100) and not (t.Z1_lepminus_author != 18 and t.Z1_lepminus_pt > 15000 and abs(t.Z1_lepminus_z0) < 2 and t.Z1_lepminus_isSA != 1): willBeRemoved = True
    if (t.Z2_lepplus_m > 100) and not (t.Z2_lepplus_author != 18 and t.Z2_lepplus_pt > 15000 and abs(t.Z2_lepplus_z0) < 2 and t.Z2_lepplus_isSA != 1): willBeRemoved = True
    if (t.Z2_lepminus_m > 100) and not (t.Z2_lepminus_author != 18 and t.Z2_lepminus_pt > 15000 and abs(t.Z2_lepminus_z0) < 2 and t.Z2_lepminus_isSA != 1): willBeRemoved = True

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
