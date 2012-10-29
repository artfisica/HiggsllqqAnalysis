import array
import os, re

import PyCintex
PyCintex.Cintex.Enable()
from ROOT import *

gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro("AtlasStyle.C") 
ROOT.gROOT.LoadMacro("AtlasLabels.C")
SetAtlasStyle()

import sys
sys.path.append('./Higgs4lepAnalysis/python/resolution_studies')

from higgs_widths import * 

from sample_list import *

# PDG values
Z_pdg_mass = 91187.6
mu_pdg_mass = 105.65836
el_pdg_mass = 0.510998910
Z_pdg_width = 2495.2


# candidate types
types = ['4mu', '2mu2e', '2e2mu', '4e']
# UNCOMMENT TO MERGE 2E2MU AND 2MU2E (not needed in resolution_studies, there's a flag now)
#types = ['4mu', '2e2mu', '4e']

# type graphics
m4l_name = {}
m4l_name['4mu'] = 'm_{#mu#mu#mu#mu}'
m4l_name['4e'] = 'm_{eeee}'
m4l_name['2mu2e'] = 'm_{#mu#muee}'
m4l_name['2e2mu'] = 'm_{ee#mu#mu}'
type_latex = {}
type_latex['4mu']   = '\mu\mu\mu\mu'
type_latex['4e']    = 'eeee'
type_latex['2mu2e'] = '\mu\mu ee'
type_latex['2e2mu'] = 'ee \mu\mu'

# Z boson types
types_2l = ['mu', 'e']

# type graphics
m2l_name = {}
m2l_name['mu'] = 'm_{#mu#mu}'
m2l_name['e'] = 'm_{ee}'


# merging lists
def mergelists(list1, list2):
  result = []
  for item in list1:
    result.append(item)
  for item in list2:
    result.append(item)
  return result

# fwhm
def fwhm(histo):
  bin1 = histo.FindFirstBinAbove(histo.GetMaximum() / 2)
  bin2 = histo.FindLastBinAbove(histo.GetMaximum() / 2)
  res  = histo.GetBinCenter(bin2) - histo.GetBinCenter(bin1)
  return res

# observable
class observable:
  def __init__(self, _name, _val, _err):
    self.name = _name
    self.val = _val
    self.err = _err
  def __str__(self):
    res = '%s = %f +/- %f' % (self.name, self.val, self.err)
    return res
  def copyFrom(self, other):
    self.name = other.name
    self.val = other.val
    self.err = other.err

# range
class observable_range:
  def __init__(self, _min, _max):
    if (_min < _max):
      self.min = _min
      self.max = _max
    else:
      self.min = _max
      self.max = _min
  def __str__(self):
    res = '[%f, %f]' % (self.min, self.max)
    return res

# fit input
class fitInput:
  def __init__(self, _sample, _histo, _type, _region):
    self.sample = _sample
    self.histo = _histo
    self.type = _type
    self.region = _region

# fit results
class fitResult:
  def __init__(self, _histo, _m, _sigma, _bulk, _frac, _fwhm):
    self.histo = _histo
    self.m = _m
    self.sigma = _sigma
    self.bulk = _bulk
    self.frac = _frac
    self.fwhm = _fwhm

# distance between observables
def get_distance(new, old):
  res = 0.
  if (old.val != 0.):
    res = (new.val - old.val) / old.val
  else:
    res = 99999
  return res

# allocate storage element
def allocate_space():
  me = {}
  for sample in samples:
    me[sample.mass] = {}
    for type in types:
      me[sample.mass][type] = {}
  return me

# allocate storage element (2-lepton version)
def allocate_space_2l():
  me = {}
  for sample in samples:
    me[sample.mass] = {}
    for type in types_2l:
      me[sample.mass][type] = {}
  return me

# luminosity
lumi = {}
lumi['4mu'] = 4.81
lumi['2mu2e'] = 4.81
lumi['2e2mu'] = 4.81
lumi['4e'] = 4.91
