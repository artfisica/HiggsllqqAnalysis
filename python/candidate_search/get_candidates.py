# create a table with the candidate list
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import array
import PyCintex
PyCintex.Cintex.Enable()
from ROOT import *

import numpy

#tag = 'data12'
tag = 'data11'
apply100GeVcut = False

applyRunLimit = False
max_runnumber = 205017

remove_weird_events = False
weird_events = [[190256, 56600537]]

output_file = TFile('data_candidates.root', 'RECREATE')
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
sys.path.append('./Higgs4lepAnalysis/python/resolution_studies')

def get_category(tree):
  doBrem = False

  # locations
  Z1_lepplus_mu_b = (abs(tree.Z1_lepplus_eta) < 1.05)
  Z1_lepplus_mu_e = (abs(tree.Z1_lepplus_eta) > 1.05)
  Z1_lepplus_e_b = (abs(tree.Z1_lepplus_eta) < 1.37)
  Z1_lepplus_e_crk = (1.37 < abs(tree.Z1_lepplus_eta) < 1.52)
  Z1_lepplus_e_e = (1.52 < abs(tree.Z1_lepplus_eta) < 2.47)
  Z1_lepminus_mu_b = (abs(tree.Z1_lepminus_eta) < 1.05)
  Z1_lepminus_mu_e = (abs(tree.Z1_lepminus_eta) > 1.05)
  Z1_lepminus_e_b = (abs(tree.Z1_lepminus_eta) < 1.37)
  Z1_lepminus_e_crk = (1.37 < abs(tree.Z1_lepminus_eta) < 1.52)
  Z1_lepminus_e_e = (1.52 < abs(tree.Z1_lepminus_eta) < 2.47)
  Z2_lepplus_mu_b = (abs(tree.Z2_lepplus_eta) < 1.05)
  Z2_lepplus_mu_e = (abs(tree.Z2_lepplus_eta) > 1.05)
  Z2_lepplus_e_b = (abs(tree.Z2_lepplus_eta) < 1.37)
  Z2_lepplus_e_crk = (1.37 < abs(tree.Z2_lepplus_eta) < 1.52)
  Z2_lepplus_e_e = (1.52 < abs(tree.Z2_lepplus_eta) < 2.47)
  Z2_lepminus_mu_b = (abs(tree.Z2_lepminus_eta) < 1.05)
  Z2_lepminus_mu_e = (abs(tree.Z2_lepminus_eta) > 1.05)
  Z2_lepminus_e_b = (abs(tree.Z2_lepminus_eta) < 1.37)
  Z2_lepminus_e_crk = (1.37 < abs(tree.Z2_lepminus_eta) < 1.52)
  Z2_lepminus_e_e = (1.52 < abs(tree.Z2_lepminus_eta) < 2.47)

  # brem
  Z1_lepplus_e_lb = (tree.Z1_lepplus_GSF_dp < 0.3)
  Z1_lepplus_e_hb = (tree.Z1_lepplus_GSF_dp > 0.3)
  Z1_lepminus_e_lb = (tree.Z1_lepminus_GSF_dp < 0.3)
  Z1_lepminus_e_hb = (tree.Z1_lepminus_GSF_dp > 0.3)
  Z2_lepplus_e_lb = (tree.Z2_lepplus_GSF_dp < 0.3)
  Z2_lepplus_e_hb = (tree.Z2_lepplus_GSF_dp > 0.3)
  Z2_lepminus_e_lb = (tree.Z2_lepminus_GSF_dp < 0.3)
  Z2_lepminus_e_hb = (tree.Z2_lepminus_GSF_dp > 0.3)

  if (tree.type == 0): # 4mu
    if (Z1_lepplus_mu_b and Z1_lepminus_mu_b and Z2_lepplus_mu_b and Z2_lepminus_mu_b):
      return 'bbbb'
    elif (Z1_lepplus_mu_b + Z1_lepminus_mu_b + Z2_lepplus_mu_b + Z2_lepminus_mu_b == 3):
      return 'bbb'
    elif (Z1_lepplus_mu_b + Z1_lepminus_mu_b + Z2_lepplus_mu_b + Z2_lepminus_mu_b == 2):
      return 'bb'
    else:
      return 'other'
  elif (tree.type == 2): # 2e2mu
    if (Z1_lepplus_e_crk or Z1_lepminus_e_crk):
      return 'onecrk_any'
    elif (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_mu_b and Z2_lepminus_mu_b):
      return 'bb_bb'
    elif (Z1_lepplus_e_b and Z1_lepminus_e_b):
      return 'bb_other'
    else:
      return 'other_other'
  elif (tree.type == 1): # 2mu2e
    if (Z2_lepplus_e_crk or Z2_lepminus_e_crk):
      return 'any_onecrk'
    elif (Z2_lepplus_e_b and Z2_lepminus_e_b and Z1_lepplus_mu_b and Z1_lepminus_mu_b):
      return 'bb_bb'
    elif (Z2_lepplus_e_b and Z2_lepminus_e_b):
      return 'other_bb'
    else:
      return 'other_other'
  elif (tree.type == 3): # 4e
    if (not doBrem):
      if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
        return 'bbbb'
      elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
        return 'onecrk'
      elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
        return 'bbb'
      else:
        return 'other'
    else:
      if (Z1_lepplus_e_lb and Z1_lepminus_e_lb and Z2_lepplus_e_lb and Z2_lepminus_e_lb):
        if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
          return 'lblb_bbbb'
        elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
          return 'lblb_onecrk'
        elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
          return 'lblb_bbb'
        else:
          return 'lblb_other'
      elif ((Z1_lepplus_e_hb or Z1_lepminus_e_hb) and (Z2_lepplus_e_hb or Z2_lepminus_e_hb)):
        if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
          return 'hbhb_bbbb'
        elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
          return 'hbhb_onecrk'
        elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
          return 'hbhb_bbb'
        else:
          return 'hbhb_other'
      else:
        if (Z1_lepplus_e_b and Z1_lepminus_e_b and Z2_lepplus_e_b and Z2_lepminus_e_b):
          return 'lbhb_bbbb'
        elif (Z1_lepplus_e_crk or Z1_lepminus_e_crk or Z2_lepplus_e_crk or Z2_lepminus_e_crk):
          return 'lbhb_onecrk'
        elif (Z1_lepplus_e_b + Z1_lepminus_e_b + Z2_lepplus_e_b + Z2_lepminus_e_b == 3):
          return 'lbhb_bbb'
        else:
          return 'lbhb_other'
  else:
    return 'bug'

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
 


#output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_nuove/*phys*/*root*')
if (tag == 'data12'):
  output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2012_v20/*phys*/*root*')
elif (tag == 'data11'):
  output_list = glob.glob('/afs/cern.ch/work/v/vippolit/ntuple_2011_v20/*phys*/*root*')

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
  file = TFile.Open(filename.rstrip('\n'))

  if not file: continue

  tree = file.Get("candidates")

  if not tree: continue;

  for entry in range(tree.GetEntries()):
    tree.GetEntry(entry)

    if (tree.selected != 1):
      continue;

    if (remove_weird_events and [tree.run, tree.event] in weird_events):
      continue

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
    if (tag == 'data12'):
      if (tree.Z1_lepplus_m < 100 and tree.Z1_lepplus_etcone20_final/tree.Z1_lepplus_pt > 0.2): willBeRemoved = True
      if (tree.Z1_lepminus_m < 100 and tree.Z1_lepminus_etcone20_final/tree.Z1_lepminus_pt > 0.2): willBeRemoved = True
      if (tree.Z2_lepplus_m < 100 and tree.Z2_lepplus_etcone20_final/tree.Z2_lepplus_pt > 0.2): willBeRemoved = True
      if (tree.Z2_lepminus_m < 100 and tree.Z2_lepminus_etcone20_final/tree.Z2_lepminus_pt > 0.2): willBeRemoved = True
     #### TRIGMATCH AND TRIGGER ON TOP  
     #if (tree.last < 11):
     #  line = '%s *** no trigger ***' % line
     #elif (tree.last < 11):
     #  line = '%s *** no trigmatch ***' % line
    if (willBeRemoved): line = '%s *** removed ***' % line

    if (tree.Z1_lepplus_author == 16 or tree.Z1_lepminus_author == 16 or tree.Z2_lepplus_author == 16 or tree.Z2_lepminus_author == 16): line = '%s *** CALO muons ***' % line
    if (tree.Z1_lepplus_isSA == 1 or tree.Z1_lepminus_isSA == 1 or tree.Z2_lepplus_isSA == 1 or tree.Z2_lepminus_isSA == 1): line = '%s *** SA muons ***' % line

    if (willBeRemoved):
      pass
      continue

    line = '%s %s' % (line, get_category(tree)) # moved it here, otherwise duplicate removal does not work for tree and counters


    if line not in list_4mu and line not in list_4e and line not in list_2e2mu:
      if ((applyRunLimit == False or tree.run <= max_runnumber) and (apply100GeVcut == False or tree.H_m > 100000)):
        newtree_type[0] = tree.type
        newtree_mass[0] = tree.H_m
        newtree_mass_constrained[0] = tree.H_m_constrained
        newtree_Z1_m[0] = tree.Z1_m
        newtree_Z2_m[0] = tree.Z2_m
        newtree.Fill()

        if (tree.H_m < 160000):
          counter_lowmass_partial[type]['<160'] += 1
        else:
          counter_lowmass_partial[type]['>160'] += 1
      else:
        continue;

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
