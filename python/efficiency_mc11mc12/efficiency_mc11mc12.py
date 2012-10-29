# compute efficiency vs mH for various samples
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import sys
import glob
import numpy

sys.path.append('Higgs4lepAnalysis/python/yield_estimator')
from chan_name_map import *
from chan_genentries_map import *


import re

regexp_ggF = re.compile('ggH(\d\d\d)')

# update maps with 8 TeV stuff (not needed at some point hopefully)
sample_name[160152] = 'PowhegPythia8_AU2CT10_ggH110_ZZ4lep'
sample_name[160153] = 'PowhegPythia8_AU2CT10_ggH115_ZZ4lep'
sample_name[160154] = 'PowhegPythia8_AU2CT10_ggH120_ZZ4lep'
sample_name[160155] = 'PowhegPythia8_AU2CT10_ggH125_ZZ4lep'
sample_name[160156] = 'PowhegPythia8_AU2CT10_ggH130_ZZ4lep'
sample_name[160157] = 'PowhegPythia8_AU2CT10_ggH135_ZZ4lep'
sample_name[160159] = 'PowhegPythia8_AU2CT10_ggH145_ZZ4lep'
sample_name[160160] = 'PowhegPythia8_AU2CT10_ggH150_ZZ4lep'
sample_name[160161] = 'PowhegPythia8_AU2CT10_ggH155_ZZ4lep'
sample_name[160162] = 'PowhegPythia8_AU2CT10_ggH160_ZZ4lep'
sample_name[160163] = 'PowhegPythia8_AU2CT10_ggH165_ZZ4lep'
sample_name[160164] = 'PowhegPythia8_AU2CT10_ggH170_ZZ4lep'
sample_name[160165] = 'PowhegPythia8_AU2CT10_ggH175_ZZ4lep'
sample_name[160166] = 'PowhegPythia8_AU2CT10_ggH180_ZZ4lep'
sample_name[160167] = 'PowhegPythia8_AU2CT10_ggH185_ZZ4lep'
sample_name[160168] = 'PowhegPythia8_AU2CT10_ggH190_ZZ4lep'
sample_name[160169] = 'PowhegPythia8_AU2CT10_ggH195_ZZ4lep'
sample_name[160170] = 'PowhegPythia8_AU2CT10_ggH200_ZZ4lep'
sample_name[160171] = 'PowhegPythia8_AU2CT10_ggH220_ZZ4lep'
sample_name[160173] = 'PowhegPythia8_AU2CT10_ggH260_ZZ4lep'
sample_name[160174] = 'PowhegPythia8_AU2CT10_ggH280_ZZ4lep'
sample_name[160175] = 'PowhegPythia8_AU2CT10_ggH300_ZZ4lep'
sample_name[160180] = 'PowhegPythia8_AU2CT10_ggH400_ZZ4lep'

gen_entries[160173] = 49900.;
gen_entries[160155] = 150000.;
gen_entries[160154] = 149998.;
gen_entries[160157] = 149798.;
gen_entries[160156] = 149998.;
gen_entries[160153] = 149998.;
gen_entries[160152] = 149900.;
gen_entries[160175] = 50000.;
gen_entries[160174] = 49900.;
gen_entries[160159] = 150000.;
gen_entries[160171] = 50000.;
gen_entries[160170] = 50000.;
gen_entries[160180] = 50000.;
gen_entries[160160] = 150000.;
gen_entries[160161] = 50000.;
gen_entries[160162] = 50000.;
gen_entries[160163] = 50000.;
gen_entries[160164] = 49900.;
gen_entries[160165] = 50000.;
gen_entries[160166] = 50000.;
gen_entries[160167] = 50000.;
gen_entries[160168] = 50000.;
gen_entries[160169] = 50000.;


filelist = glob.glob('/tmp/vippolit/*/*root')
#filelist = ['output_test.root', 'output_test.17.root']

higgs_mass = {}
generated = {}
selected  = {}

analyses = ['4mu', '2mu2e', '4e']
sample_types = ['mc11c', 'mc12'] # used just for the final plot (run number is unique!)

for analysis in analyses:
  selected[analysis] = {}

has_alerted = {}

for filename in filelist:
  print 'considering %s' % filename


  print 'opening %s' % filename
  f = TFile.Open(filename)

  h = f.Get('generatedEntriesHisto')
  e = f.Get('candidates')

  if (h and e):
    # find out which run is it
    e.GetEntry(0)

    # new run: initialize counters
    if (e.run not in generated.keys()):
      this_sample_name = sample_name[e.run]
      regexp_ggF_result = regexp_ggF.search(this_sample_name)
      if (not regexp_ggF_result): continue

      generated[e.run] = 0.
      higgs_mass[e.run] = int(regexp_ggF_result.group(1))

      for analysis in analyses:
        selected[analysis][e.run] = 0.
 
      has_alerted[e.run] = False

    # update counter of generated entries
    generated[e.run] += h.GetBinContent(4) # 4 is after mc_weight and ggF if proper
    if h.GetBinContent(4) != h.GetBinContent(3):
      print 'WARNING: run %d has ggF reweight applied' % e.run

    for entry in range(e.GetEntries()):
      e.GetEntry(entry)

      analysis_map = ['4mu', '2mu2e', '2mu2e', '4e']
      analysis = analysis_map[e.type]

      if (e.selected == 1):
        poids = e.ggF_weight

        selected[analysis][e.run] += poids

        if poids != 1. and not has_alerted[e.run]:
	  print 'WARNING: run %d has at least an event with non-zero ggF weight' % e.run
	  has_alerted[e.run] = True
  
  f.Close()



# save efficiency 
outfile = TFile('outfile_mc11mc12.root', 'RECREATE')

eff_vs_mH = {}

for analysis in analyses:
  for sample_type in sample_types:
    eff_vs_mH['%s_%s' % (analysis, sample_type)] = TEfficiency('eff_vs_mH_%s_%s' % (analysis, sample_type), '%s;Higgs mass [GeV];# %s sel. evts / # evts after 4-lep filter' % (sample_type, analysis), 1000, 0, 1000)

  for run in generated.keys():
    sample_type = 'mc12' if (run > 160152) else 'mc11c'
    print 'WARNING: %d considered as part of %s' % (run, sample_type)
 
    # rescale non-selected events for 2mu2e only, to be able to compare
    if (analysis != '2mu2e'):
      not_selected_events =       generated[run] - selected[analysis][run]
    else:
      not_selected_events = 2.0 * generated[run] - selected[analysis][run]

    for i in range(int(selected[analysis][run])):
      eff_vs_mH['%s_%s' % (analysis, sample_type)].Fill(1, higgs_mass[run])
    for j in range(int(not_selected_events)):
      eff_vs_mH['%s_%s' % (analysis, sample_type)].Fill(0, higgs_mass[run])

    print '(%s) run=%d mH=%d eff=%s' % (analysis, run, higgs_mass[run], selected[analysis][run]/generated[run])

filter_efficiency = {}
for sample_type in sample_types:
  filter_efficiency[sample_type] = TH1F("filter_efficiency_%s" % sample_type, "efficiency of 4-lepton filter;m_{H} [GeV]; filter efficiency", 1000, 0, 1000)

for run in generated.keys():
  sample_type = 'mc12' if (run > 160152) else 'mc11c'

  filter_efficiency[sample_type].Fill(higgs_mass[run], generated[run] / gen_entries[run])

outfile.Write()
outfile.Close()
