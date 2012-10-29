# create the same table as get_candidates.py but from the HSG2 twiki table text
# (it's just a parser to reduce data to the same format)
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import array
import os, re
import glob

from ROOT import TFile, TTree
import numpy

import sys
sys.path.append('./datasets')

correct_bug = True # set to true for summer2012

input_file = 'old' # summer2012
#input_file = 'old_old' # winter
input_file_constraint = 'eleni' # this is a fix, to add kostas' mass constraint

list_constrained = {}

for line in open(input_file_constraint):
  tokens = line.split('*')
  run = int(tokens[2])
  event = int(tokens[3])
  H_m_constrained = float(tokens[5])

  print 'read %d %d %lf' % (run, event, H_m_constrained)

  list_constrained['%d_%d' % (run, event)] = H_m_constrained/1000.

print ''
print ''


format = '%-6s | %6d | %10d | % 5d | % 4.2f | % 4.2f | % 4.2f'


output_file = TFile('data_candidates.root', 'RECREATE')
tree_type = numpy.zeros(1, dtype=int)
tree_mass = numpy.zeros(1, dtype=float)
tree_mass_constrained = numpy.zeros(1, dtype=float)
tree_Z1_m = numpy.zeros(1, dtype=float)
tree_Z2_m = numpy.zeros(1, dtype=float)
tree = TTree('data_candidates', 'minimal overlap-removed data candidates ntuple')
tree.Branch('type', tree_type, 'type/i')
tree.Branch('mass', tree_mass, 'mass/D')
tree.Branch('mass_constrained', tree_mass_constrained, 'mass_constrained/D')
tree.Branch('Z1_m', tree_Z1_m, 'Z1_m/D')
tree.Branch('Z2_m', tree_Z2_m, 'Z2_m/D')


for line in open(input_file):

  field = line.split('|') # 2, 3, 4, 5, 6, 7, 8 are interesting fields

  if (len(field) < 11):
    continue;

  type = -1

  if (field[2]  == ' 4&#956; '):
    channel = '4mu'
    type = 0
  elif (field[2]  == ' 2e2&#956; '):
    channel = '2e2mu'
    type = 2
  elif (field[2]  == ' 2&#956;2e '):
    channel = '2mu2e'
    type = 1
  elif (field[2]  == ' 4e '):
    channel = '4e'
    type = 3
  else:
    continue;

  run = int(field[3])
  event = int(field[4])
  lbn = int(field[5])
  H_m = float(field[6])
  if type == 3 and correct_bug: Z1_m = float(field[7]) # bugfix in the twiki :)
  else: Z1_m = float(field[8])
  if (correct_bug): Z2_m = float(field[9])
  else: Z2_m = float(field[8])

  tree_type[0] = type
  tree_mass[0] = H_m
  #tree_mass_constrained[0] = list_constrained['%d_%d' % (run, event)]
  tree_Z1_m[0] = Z1_m
  tree_Z2_m[0] = Z2_m
  tree.Fill()
    
  print format % (channel, run, event, lbn, H_m, Z1_m, Z2_m)

output_file.Write()
output_file.Close()
