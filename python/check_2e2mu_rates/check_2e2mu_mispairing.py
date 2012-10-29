# check the effect of artificially introducing (e,mu)+(e,mu) mispairing in 2e2mu events
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import array
import PyCintex
PyCintex.Cintex.Enable()
from ROOT import *


import os, re
import glob

import sys
sys.path.append('./datasets')

interesting_files = []

def passesZ2MassCut(H, Z_2):
  
   nbins = 7
   mass_4lep = [   120,   130, 150, 160, 165, 180, 190 ]
   cut_Z2    = [  17.5,  22.5,  30,  30,  35,  40,  50 ]

   my_cut = 0
   my_4lep_mass = H.M() / 1000.

   index = -1
   for j in range(0, nbins):
     if (my_4lep_mass > mass_4lep[j]): index = j

   if (index == -1): my_cut = 17.5 * 1000.
   elif (index == nbins - 1): my_cut = 50.0 * 1000.
   else: my_cut = 1000. * (cut_Z2[index] + (my_4lep_mass - mass_4lep[index]) * (cut_Z2[index + 1] - cut_Z2[index]) / (mass_4lep[index + 1] - mass_4lep[index]))

   return (Z_2.M() > my_cut and Z_2.M() < 115000)

def passesDeltaRCut(leptons):
  for i in range(0, 4):
    for j in range(i+1, 4):
      dr_cut = 0.2 if (leptons[i].M() != leptons[j].M()) else 0.1
      if (leptons[i].DeltaR(leptons[j]) < dr_cut): return False

  # ONLY 4mu and 4e
  if (leptons[0].M() == leptons[2].M()):
    if ((leptons[0]+leptons[3]).M() < 5000 or (leptons[1]+leptons[2]).M() < 5000): return False

  return True

def passesIsolationD0(momentum, ptcone20, etcone20, d0sign):
  calo_cut = 0.15 if (momentum.M() > 100 and abs(momentum.Eta()) > 2.5) else 0.30
  if (momentum.M() < 100): calo_cut = 0.20

  if (momentum.M() < 100):
    if (ptcone20/momentum.Pt() > 0.15 or etcone20/momentum.Pt() > calo_cut or abs(d0sign) > 6.5): return False
  else:
    if (ptcone20/momentum.Pt() > 0.15 or etcone20/momentum.Pt() > calo_cut or abs(d0sign) > 3.5): return False
  return True
    

#output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_nuove/*PowhegPythia8_AU2CT10_ZZ*2e2mu*/*root*')
#output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_nuove/*PowhegPythia8_AU2CT10_ZZ*4mu*/*root*')
output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2011/*109292*/*root*')

format = '%-6s | %6d | %10d | % 5d | % 4.2f | % 4.2f | % 4.2f | % 4.2f |'

channel_names = ['4mu', '2mu2e', '2e2mu', '4e']

candidates_in_event = {}
candidates_in_event_raw = {}
events_selected = {}
events_selected_raw = {}
for chan in channel_names:
  candidates_in_event[chan] = {}
  candidates_in_event_raw[chan] = {}
  events_selected[chan] = {}
  events_selected_raw[chan] = {}
  candidates_in_event[chan]['<160'] = {}
  candidates_in_event[chan]['>160'] = {}
  candidates_in_event_raw[chan]['<160'] = {}
  candidates_in_event_raw[chan]['>160'] = {}
  events_selected[chan]['<160'] = 0.
  events_selected[chan]['>160'] = 0.
  events_selected_raw[chan]['<160'] = 0.
  events_selected_raw[chan]['>160'] = 0.



for filename in output_list:
  file = TFile.Open(filename.rstrip('\n'))
  print file.GetName()

  if not file: continue

  tree = file.Get("candidates")

  if not tree: continue;

  for entry in range(tree.GetEntries()):
    tree.GetEntry(entry)

    type = "none"
    if (tree.type == 0):
      type = "4mu"
    elif (tree.type == 1):
      type = "2mu2e"
    elif (tree.type == 2):
      type = "2e2mu"
    elif (tree.type == 3):
      type = "4e"

    # determines the low/high mass region
    region = '<160' if tree.H_m < 160000 else '>160'

    # reset the candidate counter if this event has not been processed yet
    if tree.event not in candidates_in_event[type]['<160'].keys():
      candidates_in_event[type]['<160'][tree.event] = []
      candidates_in_event_raw[type]['<160'][tree.event] = []
    if tree.event not in candidates_in_event[type]['>160'].keys():
      candidates_in_event[type]['>160'][tree.event] = []
      candidates_in_event_raw[type]['>160'][tree.event] = []

    if tree.last >= 6: # events passing deltaR cut are the only ones which might be chosen against others (they undergo the "best" candidate selection)
      flag = 1 if tree.selected == 1 else 0
      candidates_in_event[type][region][tree.event].append([tree.Z1_m, tree.Z2_m, flag]) # Z1_m, Z2_m, passes_all_selection
      candidates_in_event_raw[type][region][tree.event].append([tree.Z1_m, tree.Z2_m, flag]) # Z1_m, Z2_m, passes_all_selection

    # mispair artificially: take events just after the Z1 cut (the idea is to flip e/mu and reapply the Z1,Z2 cuts to the new pairing
    if (tree.last >= 3):# and not (tree.selected == 1)):
      Z1_lepplus_4m = TLorentzVector()
      Z1_lepminus_4m = TLorentzVector()
      Z2_lepplus_4m = TLorentzVector()
      Z2_lepminus_4m = TLorentzVector()

      Z1_lepplus_4m.SetPtEtaPhiM(tree.Z1_lepplus_pt, tree.Z1_lepplus_eta, tree.Z1_lepplus_phi, tree.Z1_lepplus_m)
      Z1_lepminus_4m.SetPtEtaPhiM(tree.Z1_lepminus_pt, tree.Z1_lepminus_eta, tree.Z1_lepminus_phi, tree.Z1_lepminus_m)
      Z2_lepplus_4m.SetPtEtaPhiM(tree.Z2_lepplus_pt, tree.Z2_lepplus_eta, tree.Z2_lepplus_phi, tree.Z2_lepplus_m)
      Z2_lepminus_4m.SetPtEtaPhiM(tree.Z2_lepminus_pt, tree.Z2_lepminus_eta, tree.Z2_lepminus_phi, tree.Z2_lepminus_m)

      quadrimomenta = [Z1_lepplus_4m, Z1_lepminus_4m, Z2_lepplus_4m, Z2_lepminus_4m]

      if (not passesDeltaRCut(quadrimomenta)): continue # deltaR cut does not depend on lepton ordering, so we do it now

      # new Z1, Z2
      doInvertPairs = False # set to False to check consistency

      if (doInvertPairs):
        Z_a = Z1_lepplus_4m + Z2_lepminus_4m
        Z_b = Z2_lepplus_4m + Z1_lepminus_4m
      else:
        Z_a = Z1_lepplus_4m + Z1_lepminus_4m
        Z_b = Z2_lepplus_4m + Z2_lepminus_4m

      Z_1 = Z_a if (abs(Z_a.M() - 91187.6) < abs(Z_b.M() - 91187.6)) else Z_b
      Z_2 = Z_b if (abs(Z_a.M() - 91187.6) < abs(Z_b.M() - 91187.6)) else Z_a
      H = Z_1+Z_2

      if (abs(H.M() - tree.H_m) > 1):
        print 'ERROR: unexpected mispair mass %lf versus tree mass %lf in event %d, they should be equal!' % (H.M(), tree.H_m, tree.event)

      status = 0 # tells if the candidate, which passes all cuts up to deltaR of course, has been selected

      # Z1 mass cut
      if (Z_1.M() > 50000 and Z_1.M() < 106000):
        # Z2 mass cut
	if (passesZ2MassCut(H, Z_2)):
          if (passesIsolationD0(Z1_lepplus_4m, tree.Z1_lepplus_ptcone20_final, tree.Z1_lepplus_etcone20_final, tree.Z1_lepplus_d0/tree.Z1_lepplus_d0_sig)): 
            if (passesIsolationD0(Z1_lepminus_4m, tree.Z1_lepminus_ptcone20_final, tree.Z1_lepminus_etcone20_final, tree.Z1_lepminus_d0/tree.Z1_lepminus_d0_sig)): 
              if (passesIsolationD0(Z2_lepplus_4m, tree.Z2_lepplus_ptcone20_final, tree.Z2_lepplus_etcone20_final, tree.Z2_lepplus_d0/tree.Z2_lepplus_d0_sig)): 
                if (passesIsolationD0(Z2_lepminus_4m, tree.Z2_lepminus_ptcone20_final, tree.Z2_lepminus_etcone20_final, tree.Z2_lepminus_d0/tree.Z2_lepminus_d0_sig)): 
	          status = 1

      candidates_in_event[type][region][tree.event].append([Z_1.M(), Z_2.M(), status])

  file.Close()


print 'Computing...'


for chan in channel_names:
  for region in candidates_in_event[chan]:
    for event in candidates_in_event[chan][region]:
      closest = [-99999999, -99999999, 0] # closesttoZZ candidate: [Z1_m, Z2_m, selected]
      
      for candidate in candidates_in_event[chan][region][event]:
        if (abs(candidate[0] - 91187.6) < abs(closest[0] - 91187.6)): # case 1: this candidate has Z1 mass closer to PDG
	  closest = candidate
        elif (abs(candidate[0] - 91187.6) == abs(closest[0] - 91187.6) and candidate[1] > closest[1]): # case 2: same Z1 mass but higher Z2 mass
	  closest = candidate
        
      if closest[2] == 1:
     #if closest[2] == 1 and abs(closest[0] - 91187.6) < 15000:
        events_selected[chan][region] += 1.

for chan in channel_names:
  for region in candidates_in_event_raw[chan]:
    for event in candidates_in_event_raw[chan][region]:
      closest = [-99999999, -99999999, 0] # closesttoZZ candidate: [Z1_m, Z2_m, selected]
      
      for candidate in candidates_in_event_raw[chan][region][event]:
        if (abs(candidate[0] - 91187.6) < abs(closest[0] - 91187.6)): # case 1: this candidate has Z1 mass closer to PDG
	  closest = candidate
        elif (abs(candidate[0] - 91187.6) == abs(closest[0] - 91187.6) and candidate[1] > closest[1]): # case 2: same Z1 mass but higher Z2 mass
	  closest = candidate
        
      if closest[2] == 1:
     #if closest[2] == 1 and abs(closest[0] - 91187.6) < 15000:
        events_selected_raw[chan][region] += 1.


print ''
print ''
print ' ### RAW OUTPUT (BASELINE-SELECTED) ##### '

for chan in channel_names:
  print '[%s] (%lf < 160 - %lf > 160)' % (chan, events_selected_raw[chan]['<160'], events_selected_raw[chan]['>160'])

print ''
print ''
print ' ### AFTER CROSS-PAIRING (BASELINE-SELECTED (+) MISPAIRING) ##### '

for chan in channel_names:
  print '[%s] (%lf < 160 - %lf > 160)' % (chan, events_selected[chan]['<160'], events_selected[chan]['>160'])
