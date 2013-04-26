#!/bin/bash
import os

import samples_data11
import samples_mc11c
import samples_data12
import samples_mc12


lists_to_process_full = [
## DATA
# # 2012
  samples_data12.physics_Muons,
  samples_data12.physics_Egamma,
# # 2011
  samples_data11.physics_Muons,
  samples_data11.physics_Egamma,
#
# MC
# # 2012
  samples_mc12.signals_ggF,
  samples_mc12.signals_VBF,
  samples_mc12.signals_VH,
  samples_mc12.alpgen_Zbb,
  samples_mc12.pythia_Zjet,
  samples_mc12.alpgen_Zjet,
  samples_mc12.WZ,
  samples_mc12.ZZ,
  samples_mc12.qcd,
  samples_mc12.top,
# # 2011
  samples_mc11c.signals_ggF,
  samples_mc11c.signals_VBF,
  samples_mc11c.signals_VH,
  samples_mc11c.alpgen_Zbb,
  samples_mc11c.alpgen_Zjet,
 #samples_mc11c.WZ,
  samples_mc11c.ZZ,
  samples_mc11c.qcd,
  samples_mc11c.top,
]

lists_to_process_reduced = [
]

lists_to_process = lists_to_process_full



### analysis options
command_outputname='output_test.root'

command_name='./HiggsllqqAnalysis/bin/test'

command_opts = {}
command_opts['data12_8TeV'] = '--analysis rel_17_2 --useTopoIso --input input.txt --output %s'  % (command_outputname) # rel. 17.2 standard, with topoiso - use for data12 and mc12a processing
command_opts['mc12_8TeV'] = command_opts['data12_8TeV']
command_opts['data11_7TeV'] ='--analysis rel_17 --GSF --input input.txt --output %s'  % (command_outputname) # rel. 17 standard - use for data11 and mc11c processing
command_opts['mc11_7TeV'] = command_opts['data11_7TeV']



### grid options
username='arturos' # who is running the job
mycodeversion='Hllqq'

suffix = {}
suffix['data12_8TeV'] = '2012.01'
suffix['mc12_8TeV'] = suffix['data12_8TeV']
suffix['data11_7TeV'] = '2011.01'
suffix['mc11_7TeV'] = suffix['data11_7TeV']
filename=command_outputname
excludedSite=''
useItCloud=False
speedUp=True # use tarball
temp_dir='/tmp/arturos'

data_transfer='INFN-ROMA1_LOCALGROUPDISK' # put here INFN-ROMA1_LOCALGROUPDISK to transfer output ntuples to the localgroupdisk; leave it empty not to do it
doMerge=False # put to true to use output merge - WARNING: never tested so far!

avoidBigJobs=True # the stupid 10 Gb limit

period_name = {}
period_name['data12_8TeV'] = 'run_200804_215643' # used for data streams (data12)
period_name['data11_7TeV'] = 'run_177986_191933' # used for data streams (data11)

printOnly=True # do not execute the command, just print it

# prun usually skips files you need! check its output and put them here
files_to_preserve = [
'HiggsllqqAnalysis/packages/files/pileup/*root*',
'HiggsllqqAnalysis/packages/files/ggFHiggsPtWeight/*root*',
'HiggsZZ4lUtils/share/*',
'TrigMuonEfficiency/share/*root*',
'egammaAnalysisUtils/share/*root*',
'MuonEfficiencyCorrections/share/*txt*',
'MuonMomentumCorrections/share/*',
'egammaFourMomentumError/share/*root*',
]


#####################
#####################

import re


_suffix = suffix

sample_regexp = '(\w+)\.(\w+)\.(\w+)(\..*)\.NTUP_\w*HSG2\.([^p]+)(p\w+)' # 1: production (data/MC) 2: channel 3: name 4: 'merge' 5: AOD tag 6: D3PD tag
sample_re = re.compile(sample_regexp)

genericsample_regexp = '(\w+)\.(\d+)\.(\w+)\.(.*)\.(\w+)' # 1: production (data/MC/perf*) 2: channel 3: name 4: AOD tag 5: last tag
genericsample_re = re.compile(genericsample_regexp)

tarball_already_created = False

# utility class, to store sample info obtained from regexp
class d3pd_dataset:
  def __init__(self):
    self.run = ''
    self.name = ''
    self.tag = ''
    self.fulltag = ''

# utility function, which composes the actual prun command, given the input datasets(') list and sample properties
def compose_prun_command(d3pd, dataset_string):
  tagfound = False
  mytagname = 'NONE'
  for tagname in command_opts.keys():
    if tagname in d3pd.project:
      if not tagfound:
        prun_exec_string = '--exec \"echo %%IN > input.txt; python HiggsllqqAnalysis/python/regolari.py; cat input.txt; %s %s\"' % (command_name, command_opts[tagname])
	tagfound = True
	mytagname = tagname
      else:
        print 'ERROR: command_opts.keys() contains multiple tag names which match d3pd named "%s"...' % d3pd.project
	return
  if not tagfound:
    print 'ERROR: unable to decode tag for d3pd named "%s"...' % d3pd.project
    return

  this_command = 'prun ' + prun_exec_string + ' --outDS user.' + username + '.' + mycodeversion + '.' + d3pd.run + '.' + d3pd.name + '.' + d3pd.tag + '.' + _suffix[mytagname] + ' --outputs ' + filename + ' --useAthenaPackages --inDS ' + dataset_string
  
  if (avoidBigJobs):
    this_command += ' --nGBPerJob=10'
  
  ## add the excluded sites
  if (excludedSite != ''):
    this_command += ' --excludedSite ' + excludedSite
  
  # add the rootcore stuff
  this_command += ' --useRootCore'
  
  # preserve files
  if (len(files_to_preserve) > 0):
    this_command += ' --extFile='
    for interesting_file in files_to_preserve:
      this_command += interesting_file + ','
  
  # use italian cloud?
  if (useItCloud):
    this_command += ' --cloud=IT'

  # tarball stuff
  if (speedUp):
    if not tarball_already_created:
      this_command += ' --outTarBall tarball.tar --noSubmit'
    else:
      this_command += ' --inTarBall tarball.tar'

  # tmp directory
  if (temp_dir != ''):
    this_command += ' --tmpDir=' + temp_dir

  # DaTRI
  if (data_transfer != ''):
    this_command += ' --destSE=' + data_transfer

  # output merging
  if (doMerge):
    this_command += ' --mergeOutput'

  return this_command


# actual loop over datasets
for dataset_list in lists_to_process:
  if (len(dataset_list) > 0):
    match_re = sample_re.search(dataset_list[0])
    match_generic_re = genericsample_re.search(dataset_list[0])

    if (match_re):
      d3pd = d3pd_dataset()
      d3pd.project = match_re.group(1)
      d3pd.run = match_re.group(2)
      d3pd.name = match_re.group(3)
      d3pd.fulltag = match_re.group(5) + match_re.group(6)
      d3pd.tag = re.sub('_tid.*', '', match_re.group(6)) # protect against tid (datasets and not containers)

      isData = d3pd.project.startswith('data')
    elif (match_generic_re):
      d3pd = d3pd_dataset()
      d3pd.project = match_generic_re.group(1)
      d3pd.run = match_generic_re.group(2)
      d3pd.name = match_generic_re.group(3)
      d3pd.fulltag = match_generic_re.group(4) + match_generic_re.group(5)
      d3pd.tag = re.sub('_tid.*', '', match_generic_re.group(5)) # protect against tid (datasets and not containers)

      isData = False

    if (match_re or match_generic_re):
      if (isData):
      # DATA: one job per stream
        tagfound = False
	mytagname = ''
        for tagname in period_name.keys():
          if tagname in d3pd.project:
            if not tagfound:
              mytagname = tagname
	      tagfound = True
            else:
              print 'ERROR: command_opts.keys() contains multiple tag names which match d3pd named "%s"...' % d3pd.project
	      continue
        if not tagfound:
          print 'ERROR: unable to decode tag for d3pd named "%s"...' % d3pd.project
	  continue

	d3pd.run = period_name[mytagname]
        dataset_string = ''

        for this_dataset in dataset_list:
          dataset_string += this_dataset + ','
	dataset_string = dataset_string[:-1]

	if (speedUp and not tarball_already_created):
	  print '### TARBALL CREATION COMMAND'
	  this_command = compose_prun_command(d3pd, dataset_string)
	  print this_command
          if (not printOnly):
            os.system(this_command)
	  print '### DONE'
          tarball_already_created = True

	this_command = compose_prun_command(d3pd, dataset_string)

        print this_command
        if (not printOnly):
          os.system(this_command)
      else:
      # MC: one job per sample
        for this_dataset in dataset_list:
          dataset_string = this_dataset

          match_re = sample_re.search(this_dataset)
          match_generic_re = genericsample_re.search(this_dataset)
       
          if (match_re):
            d3pd = d3pd_dataset()
            d3pd.project = match_re.group(1)
            d3pd.run = match_re.group(2)
            d3pd.name = match_re.group(3)
            d3pd.fulltag = match_re.group(5) + match_re.group(6)
            d3pd.tag = re.sub('_tid.*', '', match_re.group(6)) # protect against tid (datasets and not containers)
          elif (match_generic_re):
            d3pd = d3pd_dataset()
            d3pd.project = match_generic_re.group(1)
            d3pd.run = match_generic_re.group(2)
            d3pd.name = match_generic_re.group(3)
            d3pd.fulltag = match_generic_re.group(4) + match_generic_re.group(5)
            d3pd.tag = re.sub('_tid.*', '', match_generic_re.group(5)) # protect against tid (datasets and not containers)
       
          if (match_re or match_generic_re):
            if (speedUp and not tarball_already_created):
              print '### TARBALL CREATION COMMAND'
              this_command = compose_prun_command(d3pd, dataset_string)
              print this_command
              if (not printOnly):
                os.system(this_command)
              print '### DONE'
              tarball_already_created = True

            this_command = compose_prun_command(d3pd, dataset_string)
            
            print this_command
            if (not printOnly):
              os.system(this_command)
          else:
	    print 'error: sample %s does not match regular expression, skipping it...' % (this_dataset)
    else:
      print 'error: sample %s does not match regular expression, skipping this list...' % (dataset_list[0])
