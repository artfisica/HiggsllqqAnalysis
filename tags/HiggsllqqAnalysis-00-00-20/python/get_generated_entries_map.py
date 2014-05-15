# utility functions to retrieve number of generated entries
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

import sys
from pyAMI.pyAMI import AMI
ami = AMI()

import samples_data11
import samples_mc11c
import samples_mc12

def getNumberOfGeneratedEntries(dataset):
  command = ['GetDatasetInfo', 'logicalDatasetName=%s/' % dataset.rstrip('/')]

  result = ami.execute(command)
  output = result.output().split()

  totalEvents=-1

  for i in xrange(len(output)):
    if ('totalEvents' == output[i]):
      totalEvents = output[i+2]
  
  return str(totalEvents)

def getChannelEntriesMap(list_of_lists):
  result = {}

  for dataset_list in list_of_lists:
    for dataset in dataset_list:
      tokens = dataset.split('.')
      channel = tokens[1]
      entries = getNumberOfGeneratedEntries(dataset)
      if (entries != -1):
        result[channel] = entries

  return result

def printChannelEntriesMap(my_map, type):
  if (type == 'python'):
    print 'gen_entries = {}'
    for key in my_map.keys():
      print 'gen_entries[%s] = %s.;' % (key, my_map[key])
  elif (type == 'C++'):
    print 'run_number = ntuple->eventinfo.mc_channel_number();'
    for key in my_map.keys():
      print 'if (this_channel == %s) gen_entries = %s.;' % (key, my_map[key])


lists_to_process = [
# MC
 #samples_mc11c.signals_ggF,
 #samples_mc11c.signals_VBF,
 #samples_mc11c.signals_VH,
 #samples_mc11c.alpgen_Zbb,
 #samples_mc11c.alpgen_Zjet,
 #samples_mc11c.ZZ,
##samples_mc11c.qcd,
 #samples_mc11c.top,
  samples_mc12.signals_ggF,
  samples_mc12.signals_VBF,
  samples_mc12.signals_VH,
  samples_mc12.alpgen_Zbb,
  samples_mc12.alpgen_Zjet,
  samples_mc12.ZZ,
# samples_mc12.qcd,
  samples_mc12.top,
]







my_map = getChannelEntriesMap(lists_to_process)

languages = ['C++', 'python']

for language in languages:
  print 'python output'
  printChannelEntriesMap(my_map, language)
