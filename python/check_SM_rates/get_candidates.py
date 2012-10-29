# create a table with the candidate list
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import array
import PyCintex
PyCintex.Cintex.Enable()
from ROOT import *

import numpy

output_file = TFile('data_candidates.root', 'RECREATE')

histo_min_dr = TH1F('histo_min_dr', 'minimum dr(l,l) per candidate', 100, 0, 5)

newtree_type = numpy.zeros(1, dtype=int)
newtree_mass = numpy.zeros(1, dtype=float)
newtree_mass_constrained = numpy.zeros(1, dtype=float)
newtree_Z1_m = numpy.zeros(1, dtype=float)
newtree_Z2_m = numpy.zeros(1, dtype=float)
newtree = TTree('data_candidates', 'minimal overlap-removed data candidates ntuple')
newtree.Branch('type', newtree_type, 'type/i')
newtree.Branch('mass', newtree_mass, 'mass/D')
newtree.Branch('mass_constrained', newtree_mass_constrained, 'mass_constrained/D')
newtree.Branch('Z1_m', newtree_Z1_m, 'Z1_m/D')
newtree.Branch('Z2_m', newtree_Z2_m, 'Z2_m/D')

import os, re
import glob

import sys
sys.path.append('./datasets')

interesting_files = []

class file_location:
  def __init__(self, run, file):
    self.run = run
    self.file = file

def find_processed_dataset(run):
  result = []
  look_for = 'data11_7TeV.%08d' % int(run)

  import sys
  sys.path.append('./Higgs4lepAnalysis/python')
  import samples_data11

  for sample in samples_data11.physics_Muons:
    if sample.startswith(look_for):
      result.append(sample)
  for sample in samples_data11.physics_Egamma:
    if sample.startswith(look_for):
      result.append(sample)
  for sample in samples_data11.debugrec_hltacc:
    if sample.startswith(look_for):
      result.append(sample)

  return result

def remove_duplicates(seq, idfun=None):  
  # order preserving 
  if idfun is None: 
    def idfun(x): return x 
  seen = {} 
  result = [] 
  for item in seq: 
    marker = idfun(item) 
    if marker in seen: continue 
    seen[marker] = 1 
    result.append(item) 
  return result
 


output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_v20/*phys*/*root*')
#output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2011/*phys*/*root*')
#output_list = glob.glob('output_test.root')

format = '%-6s | %6d | %10d | % 5d | % 4.2f | % 4.2f | % 4.2f | % 4.2f |'

list_4mu = []
list_2e2mu = []
list_4e = []

partial_runmax = 204158
counter_lowmass_partial = {}
channel_names = ['4mu', '2mu2e', '2e2mu', '4e']
for chan in channel_names:
  counter_lowmass_partial[chan] = {}
  counter_lowmass_partial[chan]['<160'] = 0.
  counter_lowmass_partial[chan]['>160'] = 0.

for filename in output_list:
#for filename in open('lista_file_dati_2011'):
  file = TFile.Open(filename.rstrip('\n'))

  if not file: continue

  tree = file.Get("candidates")

  if not tree: continue;

  for entry in range(tree.GetEntries()):
    tree.GetEntry(entry)

    if (tree.selected != 1):
      continue;

    type = "none"
    if (tree.type == 0):
      type = "4mu"
    elif (tree.type == 1):
      type = "2mu2e"
    elif (tree.type == 2):
      type = "2e2mu"
    elif (tree.type == 3):
      type = "4e"


    run = int(tree.run)

    event = int(tree.event)

    lbn = int(tree.lbn)

    H_m = float(tree.H_m) / 1000.

    H_m_constrained = float(tree.H_m_constrained) / 1000.

    Z1_m = float(tree.Z1_m) / 1000.

    Z2_m = float(tree.Z2_m) / 1000.

    line = format % (type, run, event, lbn, H_m, Z1_m, Z2_m, H_m_constrained)

    willBeRemoved = False
    if (tree.Z1_lepplus_m < 100 and tree.Z1_lepplus_etcone20_final/tree.Z1_lepplus_pt > 0.2): willBeRemoved = True
    if (tree.Z1_lepminus_m < 100 and tree.Z1_lepminus_etcone20_final/tree.Z1_lepminus_pt > 0.2): willBeRemoved = True
    if (tree.Z2_lepplus_m < 100 and tree.Z2_lepplus_etcone20_final/tree.Z2_lepplus_pt > 0.2): willBeRemoved = True
    if (tree.Z2_lepminus_m < 100 and tree.Z2_lepminus_etcone20_final/tree.Z2_lepminus_pt > 0.2): willBeRemoved = True
    if (willBeRemoved): line = '%s *** removed ***' % line

    # actually apply 2012 cuts (topoetcone20 is not in ntuple_2012 folder, just from ntuple_2012_nuove)
    if (willBeRemoved):
      line = '%s [no 2012 iso]' % line
     #continue

    if line not in list_4mu and line not in list_4e and line not in list_2e2mu:
      newtree_type[0] = tree.type
      newtree_mass[0] = tree.H_m
      newtree_mass_constrained[0] = tree.H_m_constrained
      newtree_Z1_m[0] = tree.Z1_m
      newtree_Z2_m[0] = tree.Z2_m
      newtree.Fill()

      # apply also SM cuts
      willBeRemoved = False
      alreadyPrinted = False
      
      if (tree.Z1_lepplus_m < 100) and not (tree.Z1_lepplus_pt > 15000 and (tree.Z1_lepplus_z0) < 2 and abs(tree.Z1_lepplus_d0/tree.Z1_lepplus_d0_sig) < 6): willBeRemoved = True
      if (tree.Z1_lepminus_m < 100) and not (tree.Z1_lepminus_pt > 15000 and (tree.Z1_lepminus_z0) < 2 and abs(tree.Z1_lepminus_d0/tree.Z1_lepminus_d0_sig) < 6): willBeRemoved = True
      if (tree.Z2_lepplus_m < 100) and not (tree.Z2_lepplus_pt > 15000 and (tree.Z2_lepplus_z0) < 2 and abs(tree.Z2_lepplus_d0/tree.Z2_lepplus_d0_sig) < 6): willBeRemoved = True
      if (tree.Z2_lepminus_m < 100) and not (tree.Z2_lepminus_pt > 15000 and (tree.Z2_lepminus_z0) < 2 and abs(tree.Z2_lepminus_d0/tree.Z2_lepminus_d0_sig) < 6): willBeRemoved = True

      if (willBeRemoved and not alreadyPrinted):
        line = '%s [no el pt/z0/d0]' % line
	alreadyPrinted = True
      
      if (tree.Z1_lepplus_m > 100) and not (tree.Z1_lepplus_author < 8 and tree.Z1_lepplus_pt > 15000 and abs(tree.Z1_lepplus_z0) < 2 and tree.Z1_lepplus_isSA != 1): willBeRemoved = True
      if (tree.Z1_lepminus_m > 100) and not (tree.Z1_lepminus_author < 8 and tree.Z1_lepminus_pt > 15000 and abs(tree.Z1_lepminus_z0) < 2 and tree.Z1_lepminus_isSA != 1): willBeRemoved = True
      if (tree.Z2_lepplus_m > 100) and not (tree.Z2_lepplus_author < 8 and tree.Z2_lepplus_pt > 15000 and abs(tree.Z2_lepplus_z0) < 2 and tree.Z2_lepplus_isSA != 1): willBeRemoved = True
      if (tree.Z2_lepminus_m > 100) and not (tree.Z2_lepminus_author < 8 and tree.Z2_lepminus_pt > 15000 and abs(tree.Z2_lepminus_z0) < 2 and tree.Z2_lepminus_isSA != 1): willBeRemoved = True

      if (willBeRemoved and not alreadyPrinted):
        line = '%s [no mu author/pt/z0]' % line
	alreadyPrinted = True
      
      if (not(tree.Z1_m > 66000 and tree.Z1_m < 116000)): willBeRemoved = True

      if (willBeRemoved and not alreadyPrinted):
        line = '%s [no Z1 mass]' % line
	alreadyPrinted = True

      if (not(tree.Z2_m > 66000 and tree.Z2_m < 116000)): willBeRemoved = True

      if (willBeRemoved and not alreadyPrinted):
        line = '%s [no Z2 mass]' % line
	alreadyPrinted = True
      

      # do not apply SM cuts
     #willBeRemoved = False
      if (willBeRemoved):
        print '%s (%.2lf %.2lf %.2lf %.2lf)' % (line, tree.Z1_lepplus_pt/1000., tree.Z1_lepminus_pt/1000., tree.Z2_lepplus_pt/1000., tree.Z2_lepminus_pt/1000.)
        continue

     #if (tree.H_m > 100000):
     #if (tree.run <= 203876 and run not in [201489, 203258, 203336, 203636, 203739]): # use 2.1/fb
      if (1): # use full stats
        if (tree.H_m < 160000):
          counter_lowmass_partial[type]['<160'] += 1
        else:
          counter_lowmass_partial[type]['>160'] += 1
      else:
        continue;

      # computes dr
      Z1_lepplus_4m = TLorentzVector()
      Z1_lepminus_4m = TLorentzVector()
      Z2_lepplus_4m = TLorentzVector()
      Z2_lepminus_4m = TLorentzVector()

      Z1_lepplus_4m.SetPtEtaPhiM(tree.Z1_lepplus_pt, tree.Z1_lepplus_eta, tree.Z1_lepplus_phi, tree.Z1_lepplus_m)
      Z1_lepminus_4m.SetPtEtaPhiM(tree.Z1_lepminus_pt, tree.Z1_lepminus_eta, tree.Z1_lepminus_phi, tree.Z1_lepminus_m)
      Z2_lepplus_4m.SetPtEtaPhiM(tree.Z2_lepplus_pt, tree.Z2_lepplus_eta, tree.Z2_lepplus_phi, tree.Z2_lepplus_m)
      Z2_lepminus_4m.SetPtEtaPhiM(tree.Z2_lepminus_pt, tree.Z2_lepminus_eta, tree.Z2_lepminus_phi, tree.Z2_lepminus_m)

      quadrimomenta = [Z1_lepplus_4m, Z1_lepminus_4m, Z2_lepplus_4m, Z2_lepminus_4m]

      min_dr = 9999.

      for i in range(len(quadrimomenta)):
        for j in range(len(quadrimomenta)):
	    if quadrimomenta[i].DeltaR(quadrimomenta[j]) < min_dr and i != j:
	      min_dr = quadrimomenta[i].DeltaR(quadrimomenta[j])

      histo_min_dr.Fill(min_dr)


    if (tree.type == 0):
      list_4mu.append(line)
    elif (tree.type == 3):
      list_4e.append(line)
    else:
      list_2e2mu.append(line)

    interesting_files.append(file_location(tree.run, tree.filename[0]))

  file.Close()

list_4mu.sort()
list_4e.sort()
list_2e2mu.sort()

output_file.Write()
output_file.Close()

list_4mu = remove_duplicates(list_4mu)
list_4e = remove_duplicates(list_4e)
list_2e2mu = remove_duplicates(list_2e2mu)

for item in list_4mu:
  print item
for item in list_4e:
  print item
for item in list_2e2mu:
  print item


print ''
print ''
print ' ######## '

for my_file in interesting_files:
  for its_dataset in find_processed_dataset(my_file.run):
    print 'dq2-get -f %s %s' % (my_file.file, its_dataset)


for chan in channel_names:
  print '[%s] (%lf < 160 - %lf > 160)' % (chan, counter_lowmass_partial[chan]['<160'], counter_lowmass_partial[chan]['>160'])
