import os
from ROOT import *
import common as YieldTools
import sample_chan_map as ChanMap
import re

gROOT.SetBatch(True)

# configuration
tag = 'mc11c'
bin_generated_entries = 5

base_output_dir = '.'
file_allhistos = 'all_histos.root'

channels = ['4mu', '2mu2e', '2e2mu', '4e']
lumi = {}
lumi['4mu'] = 4.8
lumi['2mu2e'] = 4.8
lumi['2e2mu'] = 4.8
lumi['4e'] = 4.9

regions = ['all']
region_criteria = {}
region_criteria['all'] = '1' # inclusive

mode_cut = {}
mode_cut[''] = 'selected==1'
mode_cut['_noisoHI'] = 'closesttozz==1 && Z1_lepplus_ptcone20_final/Z1_lepplus_pt < 0.15 && Z1_lepminus_ptcone20_final/Z1_lepminus_pt < 0.15 && ((Z1_lepplus_isSA != 1 && Z1_lepplus_etcone20_final/Z1_lepplus_pt < 0.30) || (Z1_lepplus_isSA == 1 && Z1_lepplus_etcone20_final/Z1_lepplus_pt < 0.15)) && ((Z1_lepminus_isSA != 1 && Z1_lepminus_etcone20_final/Z1_lepminus_pt < 0.30) || (Z1_lepminus_isSA == 1 && Z1_lepminus_etcone20_final/Z1_lepminus_pt < 0.15)) && ((Z1_lepplus_m < 100 && abs(Z1_lepplus_d0/Z1_lepplus_d0_sig) < 6.5) || (Z1_lepplus_m > 100 && abs(Z1_lepplus_d0/Z1_lepplus_d0_sig) < 3.5)) && ((Z1_lepminus_m < 100 && abs(Z1_lepminus_d0/Z1_lepminus_d0_sig) < 6.5) || (Z1_lepminus_m > 100 && abs(Z1_lepminus_d0/Z1_lepminus_d0_sig) < 3.5))'
mode_cut['_noisoVAL'] = 'Z2_lepplus_ptcone20/Z2_lepplus_pt < 0.30 && Z2_lepminus_ptcone20/Z2_lepminus_pt < 0.30 && %s' % mode_cut['_noisoHI']
mode_cut['_noisoLO'] = 'Z2_lepplus_ptcone20/Z2_lepplus_pt < 0.15 && Z2_lepminus_ptcone20/Z2_lepminus_pt < 0.15 && %s' % mode_cut['_noisoHI']

vars = ['H_m', 'H_m_constrained', 'Z1_m', 'Z2_m']

nbins = {}
minx = {}
maxx = {}
nbins['H_m'] = 2000; minx['H_m'] = 0.; maxx['H_m'] = 1000.
nbins['H_m_constrained'] = 2000; minx['H_m_constrained'] = 0.; maxx['H_m_constrained'] = 1000.
nbins['Z1_m'] = 240; minx['Z1_m'] = 30.; maxx['Z1_m'] = 150.
nbins['Z2_m'] = 150; minx['Z2_m'] = 0.; maxx['Z2_m'] = 150.


meta_samples = {} # e.g. Higgs @ 125 GeV (composed by: ggF, VBF, WH, ZH)
physics_samples = {} # e.g. Zbb

if (tag == 'mc11c'):
  meta_samples['PowHegPythia_H125'] = ['PowHegPythia_ggH125_ZZ4lep', 'PowHegPythia_VBFH125_ZZ4lep', 'PythiaWH125_ZZ4lep', 'PythiaZH125_ZZ4lep']
  meta_samples['PowHegPythia_H130'] = ['PowHegPythia_ggH130_ZZ4lep', 'PowHegPythia_VBFH130_ZZ4lep', 'PythiaWH130_ZZ4lep', 'PythiaZH130_ZZ4lep']
  meta_samples['PowHegPythia_H150'] = ['PowHegPythia_ggH150_ZZ4lep', 'PowHegPythia_VBFH150_ZZ4lep', 'PythiaWH150_ZZ4lep', 'PythiaZH150_ZZ4lep']
  meta_samples['PowHegPythia_H190'] = ['PowHegPythia_ggH190_ZZ4lep', 'PowHegPythia_VBFH190_ZZ4lep', 'PythiaWH190_ZZ4lep', 'PythiaZH190_ZZ4lep']
  meta_samples['PowHegPythia_H200'] = ['PowHegPythia_ggH200_ZZ4lep', 'PowHegPythia_VBFH200_ZZ4lep', 'PythiaWH200_ZZ4lep', 'PythiaZH200_ZZ4lep']
  meta_samples['PowHegPythia_H400'] = ['PowHegPythia_ggH400_ZZ4lep', 'PowHegPythia_VBFH400_ZZ4lep']
  meta_samples['PowHegPythia_H600'] = ['PowHegPythia_ggH600_ZZ4lep', 'PowHegPythia_VBFH600_ZZ4lep']
  meta_samples['AlpgenHWfZeebb_4LepM'] = ['AlpgenHWfZeebbNp0_4LepM', 'AlpgenHWfZeebbNp1_4LepM', 'AlpgenHWfZeebbNp2_4LepM', 'AlpgenHWfZeebbNp3_4LepM']
  meta_samples['AlpgenHWfZmumubb_4LepM'] = ['AlpgenHWfZmumubbNp0_4LepM', 'AlpgenHWfZmumubbNp1_4LepM', 'AlpgenHWfZmumubbNp2_4LepM', 'AlpgenHWfZmumubbNp3_4LepM']
  meta_samples['AlpgenHWfZeebb_Veto4LepM_Pass3Lep'] = ['AlpgenHWfZeebbNp0_Veto4LepM_Pass3Lep', 'AlpgenHWfZeebbNp1_Veto4LepM_Pass3Lep', 'AlpgenHWfZeebbNp2_Veto4LepM_Pass3Lep', 'AlpgenHWfZeebbNp3_Veto4LepM_Pass3Lep']
  meta_samples['AlpgenHWfZmumubb_Veto4LepM_Pass3Lep'] = ['AlpgenHWfZmumubbNp0_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubbNp1_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubbNp2_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubbNp3_Veto4LepM_Pass3Lep']
  meta_samples['AlpgenJimmyLowMassDYeebb_nofilter'] = ['AlpgenJimmyLowMassDYeebbNp0_nofilter', 'AlpgenJimmyLowMassDYeebbNp1_nofilter', 'AlpgenJimmyLowMassDYeebbNp2_nofilter', 'AlpgenJimmyLowMassDYeebbNp3_nofilter']
  meta_samples['AlpgenJimmyLowMassDYmumubb_nofilter'] = ['AlpgenJimmyLowMassDYmumubbNp0_nofilter', 'AlpgenJimmyLowMassDYmumubbNp1_nofilter', 'AlpgenJimmyLowMassDYmumubbNp2_nofilter', 'AlpgenJimmyLowMassDYmumubbNp3_nofilter']
  meta_samples['AlpgenJimmyLowMassDYtautaubb_nofilter'] = ['AlpgenJimmyLowMassDYtautaubbNp0_nofilter', 'AlpgenJimmyLowMassDYtautaubbNp1_nofilter', 'AlpgenJimmyLowMassDYtautaubbNp2_nofilter', 'AlpgenJimmyLowMassDYtautaubbNp3_nofilter']
  meta_samples['AlpgenJimmyZee_pt20'] = ['AlpgenJimmyZeeNp0_pt20', 'AlpgenJimmyZeeNp1_pt20', 'AlpgenJimmyZeeNp2_pt20', 'AlpgenJimmyZeeNp3_pt20', 'AlpgenJimmyZeeNp4_pt20', 'AlpgenJimmyZeeNp5_pt20']
  meta_samples['AlpgenJimmyZmumu_pt20'] = ['AlpgenJimmyZmumuNp0_pt20', 'AlpgenJimmyZmumuNp1_pt20', 'AlpgenJimmyZmumuNp2_pt20', 'AlpgenJimmyZmumuNp3_pt20', 'AlpgenJimmyZmumuNp4_pt20', 'AlpgenJimmyZmumuNp5_pt20']
  meta_samples['AlpgenJimmyZtautau_pt20'] = ['AlpgenJimmyZtautauNp0_pt20', 'AlpgenJimmyZtautauNp1_pt20', 'AlpgenJimmyZtautauNp2_pt20', 'AlpgenJimmyZtautauNp3_pt20', 'AlpgenJimmyZtautauNp4_pt20', 'AlpgenJimmyZtautauNp5_pt20']
  meta_samples['AlpgenJimmyZee_Mll10to40_pt20'] = ['AlpgenJimmyZeeNp0_Mll10to40_pt20', 'AlpgenJimmyZeeNp1_Mll10to40_pt20', 'AlpgenJimmyZeeNp2_Mll10to40_pt20', 'AlpgenJimmyZeeNp3_Mll10to40_pt20', 'AlpgenJimmyZeeNp4_Mll10to40_pt20', 'AlpgenJimmyZeeNp5_Mll10to40_pt20']
  meta_samples['AlpgenJimmyZmumu_Mll10to40_pt20'] = ['AlpgenJimmyZmumuNp0_Mll10to40_pt20', 'AlpgenJimmyZmumuNp1_Mll10to40_pt20', 'AlpgenJimmyZmumuNp2_Mll10to40_pt20', 'AlpgenJimmyZmumuNp3_Mll10to40_pt20', 'AlpgenJimmyZmumuNp4_Mll10to40_pt20', 'AlpgenJimmyZmumuNp5_Mll10to40_pt20']
  meta_samples['Pythiazz4l_3MultiLeptonFilterElecMu'] = ['Pythiazz4l_3MultiLeptonFilterElecMu']
  meta_samples['PowHegZZ_trilep5GeV_Pythia'] = ['PowHegZZ_4e_trilep5GeV_Pythia', 'PowHegZZ_4mu_trilep5GeV_Pythia', 'PowHegZZ_2e2mu_trilep5GeV_Pythia', 'PowHegZZ_2mu2tau_trilep5GeV_Pythia', 'PowHegZZ_2e2tau_trilep5GeV_Pythia', 'PowHegZZ_4tau_trilep5GeV_Pythia']
  meta_samples['gg2ZZ_JIMMY'] = ['gg2ZZ_JIMMY_ZZ4e', 'gg2ZZ_JIMMY_ZZ4mu', 'gg2ZZ_JIMMY_ZZ2e2mu']
  meta_samples['T1_McAtNlo_Jimmy'] = ['T1_McAtNlo_Jimmy', 'T1_McAtNlo_Jimmy_4LepMass_Mll60GeV12GeV']#, 'T1_McAtNlo_Jimmy_4LepMass_Mll60GeV12GeV_tarrade']

  physics_samples['ZZ'] =  ['PowHegZZ_trilep5GeV_Pythia', 'gg2ZZ_JIMMY']
  physics_samples['Z'] =   ['AlpgenJimmyZee_pt20', 'AlpgenJimmyZmumu_pt20', 'AlpgenJimmyZtautau_pt20', 'AlpgenJimmyZee_Mll10to40_pt20', 'AlpgenJimmyZmumu_Mll10to40_pt20']
  physics_samples['Zbb'] = ['AlpgenHWfZeebb_4LepM', 'AlpgenHWfZmumubb_4LepM', 'AlpgenHWfZeebb_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubb_Veto4LepM_Pass3Lep', 'AlpgenJimmyLowMassDYeebb_nofilter', 'AlpgenJimmyLowMassDYmumubb_nofilter', 'AlpgenJimmyLowMassDYtautaubb_nofilter']
  physics_samples['tt'] =  ['T1_McAtNlo_Jimmy']
  physics_samples['125'] = ['PowHegPythia_H125']
  physics_samples['130'] = ['PowHegPythia_H130']
  physics_samples['150'] = ['PowHegPythia_H150']
  physics_samples['190'] = ['PowHegPythia_H190']
  physics_samples['200'] = ['PowHegPythia_H200']
  physics_samples['400'] = ['PowHegPythia_H400']
  physics_samples['600'] = ['PowHegPythia_H600']
elif (tag == 'mc12a'):
  pass

# obtains the list of channel numbers which fill a given histogram
physics_sample_runs = {}
for physics_sample in physics_samples:
  physics_sample_runs[physics_sample] = []

  for meta_sample in physics_samples[physics_sample]:
    for sample in meta_samples[meta_sample]:
      physics_sample_runs[physics_sample].append(ChanMap.sample_number[sample])




# weights to be applied
if (tag == 'mc11c'):
 #os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2011/*/*root* | grep -v "physics" > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data)
  os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2011/*/*root* | grep -v "physics" | grep ggH130 > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data)
  poids = 'pu_weight * trigSF_weight * Z1_lepplus_weight * Z1_lepminus_weight * Z2_lepplus_weight * Z2_lepminus_weight * ggF_weight * top_weight * powhegbug_weight'

elif (tag == 'mc12a'):
  os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2012/*/*root* | grep -v "physics" > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data)
  poids = 'pu_weight * trigSF_weight * Z1_lepplus_weight * Z1_lepminus_weight * Z2_lepplus_weight * Z2_lepminus_weight * top_weight * powhegbug_weight * vxz_weight'

poids_pythiaZZ = 'xsec_weight * %s' % poids


# histograms (one per channel)
output = TFile('%s/%s' % (base_output_dir, file_allhistos), 'RECREATE')
samples = {}

# lists of run numbers
list_of_reducible_backgrounds = [y for x in ['Z', 'Zbb', 'tt'] for y in physics_sample_runs[x]]

# sample_name map
sample_name = dict((ChanMap.sample_number[k], k) for k in ChanMap.sample_number)



### VERBOSE
print 'list of runs within each physics sample:'
print physics_sample_runs
print 'list of reducible backgrounds:'
print list_of_reducible_backgrounds
print 'sample_name map:'
print sample_name


# read input files and retrieve histograms
for filename in open('%s/__TMP_CUTFLOWMAKER__' % base_output_dir):
  print 'opening %s' % filename.rstrip('\n')
  f = TFile.Open(filename.rstrip('\n'))

  if f:
    h = f.Get('generatedEntriesHisto')
    c = f.Get('cutflow')
    t = f.Get('candidates')

    if (h and t and c):
      t_Draw = t.Draw # book method [for speed]

      c.GetEntry(0)
      run = c.run

      is_redux = (run in list_of_reducible_backgrounds)

      # how must an histo be created? default: selected == 1; reducible: closesttozz and other stuff
      modes = [''] if not is_redux else ['', '_noisoLO', '_noisoHI', '_noisoVAL']

      if (run not in samples.keys()): 
	samples[run] = YieldTools.sample(sample_name[run])

        for mode in modes:
  	  for chan in channels:
            for region in regions:
              for var in vars:
	        this_histoname = 'h_%d_%s_%s_%s%s' % (run, var, chan, region, mode)
	        samples[run].histo[this_histoname] = TH1F(this_histoname, 'histogram of %s; %s; entries' % (var, var), nbins[var], minx[var], maxx[var])
	        samples[run].histo[this_histoname].SetDirectory(output)
	        samples[run].final[this_histoname] = 0.

      samples[run].generated += h.GetBinContent(bin_generated_entries)

      output.cd()

      for mode in modes:
        for chanid, chan in enumerate(channels):
          for region in regions:
            selection = 'type==%d && %s && %s' % (chanid, mode_cut[mode], region_criteria[region])

	    for var in vars:
	      this_histoname = 'h_%d_%s_%s_%s%s' % (run, var, chan, region, mode)

	      draw_what = '%s/1000.>>+%s' % (var, this_histoname)
	      draw_when = '(%s)*(%s)' % (selection, poids)

              print draw_what
              print draw_when

              t_Draw(draw_what, draw_when, 'goff')

  f.Close()


# rescale to luminosity
gSystem.AddIncludePath('-IHiggsZZ4lUtils')
gSystem.Load('libHiggsZZ4lUtils')
ROOT.gSystem.Load('libCint')
gROOT.ProcessLine('.L HiggsZZ4lUtils/HiggsZZ4lUtils/CrossSection.h+')
gROOT.ProcessLine('.L HiggsZZ4lUtils/HiggsZZ4lUtils/HiggsCrossSection.h+')
gROOT.ProcessLine('.L HiggsZZ4lUtils/HiggsZZ4lUtils/BkgCrossSection.h+')

regexp_mass_ggF = re.compile('ggH(\d\d\d)')
regexp_mass_VBF = re.compile('VBFH(\d\d\d)')
regexp_mass_WH = re.compile('WH(\d\d\d)')
regexp_mass_ZH = re.compile('ZH(\d\d\d)')

my_COM_energy = CrossSections.SevenTeV if tag == 'mc11c' else CrossSections.EightTeV

for run in samples.keys():
  my_xsec = 0.

  bkg_xsec = CrossSections.GetBkgCrossSection7TeV(run, False)

  if (bkg_xsec != -1):
    my_xsec = bkg_xsec
  else:
    my_name = samples[run].name

    # look for ggF
    parse_result = regexp_mass_ggF.search(my_name)
    if (parse_result):
      my_xsec = CrossSections.higgs4lxsecGGF(float(parse_result.group(1)), my_COM_energy)
    # look for VBF
    parse_result = regexp_mass_VBF.search(my_name)
    if (parse_result):
      my_xsec = CrossSections.higgs4lxsecVBF(float(parse_result.group(1)), my_COM_energy)
    # look for WH
    parse_result = regexp_mass_WH.search(my_name)
    if (parse_result):
      my_xsec = CrossSections.higgs4lxsecWH(float(parse_result.group(1)), my_COM_energy)
    # look for ZH
    parse_result = regexp_mass_ZH.search(my_name)
    if (parse_result):
      my_xsec = CrossSections.higgs4lxsecZH(float(parse_result.group(1)), my_COM_energy)

  samples[run].xsec = my_xsec

  # rescale
  for histo in samples[run].histo:
    print histo
    samples[run].histo[histo].Scale(lumi['4mu'] * samples[run].xsec / samples[run].generated)



###

for run in samples.keys():
  print samples[run]

output.Write()
output.Close()
