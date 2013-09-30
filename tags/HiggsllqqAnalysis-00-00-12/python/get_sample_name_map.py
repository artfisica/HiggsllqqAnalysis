# utility functions to retrieve name of a sample
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import sys
from pyAMI.pyAMI import AMI
ami = AMI()

import samples_mc11c
import samples_mc12

def getChannelNameMap(list_of_lists):
  result = {}

  for dataset_list in list_of_lists:
    for dataset in dataset_list:
      tokens = dataset.split('.')
      channel = tokens[1]
      name = tokens[2]
      result[channel] = name

  return result

def printChannelEntriesMap(my_map, type):
  if (type == 'python'):
    print 'sample_name = {}'
    for key in my_map.keys():
      print 'sample_name[%s] = \'%s\';' % (key, my_map[key])
  elif (type == 'C++'):
    print 'run_number = ntuple->eventinfo.mc_channel_number();'
    for key in my_map.keys():
      print 'if (this_channel == %s) sample_name = "%s";' % (key, my_map[key])


lists_to_process = [
# MC
  samples_mc11c.signals_ggF,
  samples_mc11c.signals_VBF,
  samples_mc11c.signals_VH,
  samples_mc11c.alpgen_Zbb,
  samples_mc11c.alpgen_Zjet,
  samples_mc11c.ZZ,
# samples_mc11c.qcd,
  samples_mc11c.top,

  samples_mc12.signals_ggF,
  samples_mc12.signals_VBF,
  samples_mc12.signals_VH,
  samples_mc12.alpgen_Zbb,
  samples_mc12.alpgen_Zjet,
  samples_mc12.ZZ,
# samples_mc12.qcd,
  samples_mc12.top,
]







my_map = getChannelNameMap(lists_to_process)

languages = ['C++', 'python']

for language in languages:
  print 'python output'
  printChannelEntriesMap(my_map, language)