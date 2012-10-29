# extract the final yield in the four final states, and produce the plots which input the workspace
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

##################
#########
######################## LA NORMALIZZAZIONE DI QUESTI CAZZO DI m12 E m34 E' SBAGLIATA SE SI TAGLIA A 100 GeV!!!!!!! PORCA TROIA!!!!!!
##
###### CERCA MISSING ALPGEN MC12 PER I CAZZI DI ALPGEN CHE NON CI STA ANCORA
#
#   AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp4Excl_Mll10to60
###################################### THIS IS MISSING!!! (and has been commented out)

from ROOT import *
import os
from optparse import OptionParser
import numpy
import re
from chan_genentries_map import *
#from chan_name_map import * # REMOVED!!
from systematics_map import *
from uncertainty_computer import *
from sample_chan_map import *
from common import *

gROOT.SetBatch(True)

gSystem.AddIncludePath('-IHiggsZZ4lUtils')
gSystem.Load('libHiggsZZ4lUtils')
gSystem.Load('libCint')
gROOT.ProcessLine('.L HiggsZZ4lUtils/HiggsZZ4lUtils/CrossSection.h+')
gROOT.ProcessLine('.L HiggsZZ4lUtils/HiggsZZ4lUtils/HiggsCrossSection.h+')
gROOT.ProcessLine('.L HiggsZZ4lUtils/HiggsZZ4lUtils/BkgCrossSection.h+')

### CONFIGURATION

### command line options

parser = OptionParser()
parser.add_option("-t", "--tag", dest="tag", help="specify the year tag (mc11c for 2011, mc12a for 2012)", default="mc12a")
parser.add_option("-d", "--directory", dest="directory", help="write output in the directory DIR", metavar="DIR", default=".")
parser.add_option("-c", "--constraint", action="store_true", default=False, help="apply mass constraint", dest="use_mass_constraint")
parser.add_option("-l", "--lower-cut", action="store_true", default=False, help="apply lower mass cut (100 GeV)", dest="do_apply_100GeV_cut")
parser.add_option("-k", "--do-categories", action="store_true", default=False, help="run in Z2 position categories", dest="do_categories")
parser.add_option("-w", "--check-window", action="store_true", default=False, help="replace low mass bin with 125pm6 GeV window for yield", dest="test_125_plusminus_6")

(options, args) = parser.parse_args()

tag = options.tag
base_output_dir = options.directory
do_apply_100GeV_cut = options.do_apply_100GeV_cut
use_mass_constraint = options.use_mass_constraint # triggers mass_to_use (histo and ZZ tree filling only) to be the value after constraint
test_125_plusminus_6 = options.test_125_plusminus_6

print ''
print ''
print 'CONFIGURATION:'
print '   - tag = %s' % str(tag)
print '   - base_output_dir = %s' % str(base_output_dir)
print '   - do_apply_100GeV_cut = %s' % str(do_apply_100GeV_cut)
print '   - use_mass_constraint = %s' % str(use_mass_constraint)
print '   - test_125_plusminus_6 = %s' % str(test_125_plusminus_6)
print ''
print ''
print ''
print ''


### channels
channels = ['4mu', '2mu2e', '2e2mu', '4e']

map_chanid_chan = {}
map_chanid_chan['4mu'] = 0
map_chanid_chan['2mu2e'] = 1
map_chanid_chan['2e2mu'] = 3
map_chanid_chan['4e'] = 2

mass_regions = ['all', '<160', '>160']

### physics


if (tag == 'mc11c'):
 #bin_for_generated = 2 # valerio
  bin_for_generated = 5
elif (tag == 'mc12a'):
  bin_for_generated = 5

tmp_outfile = TFile('%s/__yields_tmp_outfile.root' % base_output_dir, 'RECREATE')

if (tag == 'mc11c'):
 #os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2011/*/*root* | grep -v "physics" > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data)
  os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2011_v20/*/*root* | grep -v "physics" > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data)

  COM_energy = CrossSections.SevenTeV

  lumi = {}
  lumi['4mu'] = 4.8
  lumi['2mu2e'] = 4.8
  lumi['2e2mu'] = 4.8
  lumi['4e'] = 4.9

  lumi_uncertainty_relative = 0.039

  ES_shift_up = {}
  ES_shift_down = {}
  ES_shift_up['4mu'] = 0.
  ES_shift_down['4mu'] = 0.
  ES_shift_up['2mu2e'] = -0.25646/130.
  ES_shift_down['2mu2e'] = -0.168402/130.
  ES_shift_up['2e2mu'] = -0.578834/130.
  ES_shift_down['2e2mu'] = 0.376814/130.
  ES_shift_up['4e'] = -0.842431/130.
  ES_shift_down['4e'] = 0.576647/130.# 1/GeV [it's the histo average shift normalized to 130 GeV]

elif (tag == 'mc12a'):
 #os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2012_urtime/*/*root* | grep -v "physics" > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data)
 #os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2012_v20/*/*root* | grep -v "physics" > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data)
  os.system('ls /afs/cern.ch/work/v/vippolit/ntuple_2012_v20/*/*root* | grep -v "physics" | grep -v "ggH" | grep -v "VBFH" | grep -v WH | grep -v ZH | grep -v ZZ_ | grep -v gg2ZZ > %s/__TMP_CUTFLOWMAKER__' % base_output_dir) # exludes physics streams (data) and signals and ZZ

  COM_energy = CrossSections.EightTeV

  lumi = {}
  lumi['4mu'] = 5.831
  lumi['2mu2e'] = 5.831
  lumi['2e2mu'] = 5.831
  lumi['4e'] = 5.858

  lumi_uncertainty_relative = 0.039

  ES_shift_up = {}
  ES_shift_down = {}
  ES_shift_up['4mu'] = 0.
  ES_shift_down['4mu'] = 0.
  ES_shift_up['2mu2e'] = -0.182181/130.
  ES_shift_down['2mu2e'] = 0.267371/130.
  ES_shift_up['2e2mu'] = -0.402527/130.
  ES_shift_down['2e2mu'] = 0.610161/130.
  ES_shift_up['4e'] = -0.600013/130.
  ES_shift_down['4e'] = 0.896267/130.# 1/GeV [it's the histo average shift normalized to 130 GeV]



### SAMPLE ORGANIZATION

### tools to obtain mass from sample name
sample_mass = {}
regexp_mass_ggF = re.compile('ggH(\d\d\d)')
regexp_mass_VBF = re.compile('VBFH(\d\d\d)')
regexp_mass_WH = re.compile('WH(\d\d\d)')
regexp_mass_ZH = re.compile('ZH(\d\d\d)')


### tools to obtain interesting combinations of samples
meta_samples = {}

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

  #list_of_ZZ_samples = ['Pythiazz4l_3MultiLeptonFilterElecMu']
  list_of_ZZ_samples = ['PowHegZZ_trilep5GeV_Pythia', 'gg2ZZ_JIMMY']
  list_of_ZZ_samples_nogg2ZZ = ['PowHegZZ_trilep5GeV_Pythia']

elif (tag == 'mc12a'):
  meta_samples['PowhegPythia8_AU2CT10_H125'] = ['PowhegPythia8_AU2CT10_ggH125_ZZ4lep', 'PowhegPythia8_AU2CT10_VBFH125_ZZ4lep', 'Pythia8_AU2CTEQ6L1_WH125_ZZ4lep', 'Pythia8_AU2CTEQ6L1_ZH125_ZZ4lep']
  meta_samples['PowhegPythia8_AU2CT10_H130'] = ['PowhegPythia8_AU2CT10_ggH130_ZZ4lep', 'PowhegPythia8_AU2CT10_VBFH130_ZZ4lep', 'Pythia8_AU2CTEQ6L1_WH130_ZZ4lep', 'Pythia8_AU2CTEQ6L1_ZH130_ZZ4lep']
  meta_samples['PowhegPythia8_AU2CT10_H150'] = ['PowhegPythia8_AU2CT10_ggH150_ZZ4lep', 'PowhegPythia8_AU2CT10_VBFH150_ZZ4lep', 'Pythia8_AU2CTEQ6L1_WH150_ZZ4lep', 'Pythia8_AU2CTEQ6L1_ZH150_ZZ4lep']
  meta_samples['PowhegPythia8_AU2CT10_H190'] = ['PowhegPythia8_AU2CT10_ggH190_ZZ4lep', 'PowhegPythia8_AU2CT10_VBFH190_ZZ4lep', 'Pythia8_AU2CTEQ6L1_WH190_ZZ4lep', 'Pythia8_AU2CTEQ6L1_ZH190_ZZ4lep']
  meta_samples['PowhegPythia8_AU2CT10_H200'] = ['PowhegPythia8_AU2CT10_ggH200_ZZ4lep', 'PowhegPythia8_AU2CT10_VBFH200_ZZ4lep', 'Pythia8_AU2CTEQ6L1_WH200_ZZ4lep', 'Pythia8_AU2CTEQ6L1_ZH200_ZZ4lep']
  meta_samples['PowhegPythia8_AU2CT10_H400'] = ['PowhegPythia8_AU2CT10_ggH400_ZZ4lep', 'PowhegPythia8_AU2CT10_VBFH400_ZZ4lep']
  meta_samples['PowhegPythia8_AU2CT10_H600'] = ['PowhegPythia8_AU2CT10_ggH600_ZZ4lep', 'PowhegPythia8_AU2CT10_VBFH600_ZZ4lep']
  meta_samples['AlpgenHWfZeebb_4LepM'] = ['AlpgenHWfZeebbNp0_4LepM', 'AlpgenHWfZeebbNp1_4LepM', 'AlpgenHWfZeebbNp2_4LepM', 'AlpgenHWfZeebbNp3_4LepM'] # MISSING ALPGEN MC12
  meta_samples['AlpgenHWfZmumubb_4LepM'] = ['AlpgenHWfZmumubbNp0_4LepM', 'AlpgenHWfZmumubbNp1_4LepM', 'AlpgenHWfZmumubbNp2_4LepM', 'AlpgenHWfZmumubbNp3_4LepM'] # MISSING ALPGEN MC12
  meta_samples['AlpgenHWfZeebb_Veto4LepM_Pass3Lep'] = ['AlpgenHWfZeebbNp0_Veto4LepM_Pass3Lep', 'AlpgenHWfZeebbNp1_Veto4LepM_Pass3Lep', 'AlpgenHWfZeebbNp2_Veto4LepM_Pass3Lep', 'AlpgenHWfZeebbNp3_Veto4LepM_Pass3Lep'] # MISSING ALPGEN MC12
  meta_samples['AlpgenHWfZmumubb_Veto4LepM_Pass3Lep'] = ['AlpgenHWfZmumubbNp0_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubbNp1_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubbNp2_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubbNp3_Veto4LepM_Pass3Lep'] # MISSING ALPGEN MC12
  meta_samples['AlpgenJimmyLowMassDYeebb_nofilter'] = ['AlpgenJimmyLowMassDYeebbNp0_nofilter', 'AlpgenJimmyLowMassDYeebbNp1_nofilter', 'AlpgenJimmyLowMassDYeebbNp2_nofilter', 'AlpgenJimmyLowMassDYeebbNp3_nofilter'] # MISSING ALPGEN MC12
  meta_samples['AlpgenJimmyLowMassDYmumubb_nofilter'] = ['AlpgenJimmyLowMassDYmumubbNp0_nofilter', 'AlpgenJimmyLowMassDYmumubbNp1_nofilter', 'AlpgenJimmyLowMassDYmumubbNp2_nofilter', 'AlpgenJimmyLowMassDYmumubbNp3_nofilter'] # MISSING ALPGEN MC12
  meta_samples['AlpgenJimmyLowMassDYtautaubb_nofilter'] = ['AlpgenJimmyLowMassDYtautaubbNp0_nofilter', 'AlpgenJimmyLowMassDYtautaubbNp1_nofilter', 'AlpgenJimmyLowMassDYtautaubbNp2_nofilter', 'AlpgenJimmyLowMassDYtautaubbNp3_nofilter'] # MISSING ALPGEN MC12
  meta_samples['AlpgenJimmy_AUET2CTEQ6L1_Zee'] = ['AlpgenJimmy_AUET2CTEQ6L1_ZeeNp0', 'AlpgenJimmy_AUET2CTEQ6L1_ZeeNp1', 'AlpgenJimmy_AUET2CTEQ6L1_ZeeNp2', 'AlpgenJimmy_AUET2CTEQ6L1_ZeeNp3', 'AlpgenJimmy_AUET2CTEQ6L1_ZeeNp4', 'AlpgenJimmy_AUET2CTEQ6L1_ZeeNp5']
  meta_samples['AlpgenJimmy_AUET2CTEQ6L1_Zmumu'] = ['AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp0', 'AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp1', 'AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp2', 'AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp3', 'AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp4', 'AlpgenJimmy_AUET2CTEQ6L1_ZmumuNp5']
  meta_samples['AlpgenJimmy_AUET2CTEQ6L1_Ztautau'] = ['AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp0', 'AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp1', 'AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp2', 'AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp3', 'AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp4', 'AlpgenJimmy_AUET2CTEQ6L1_ZtautauNp5']
  meta_samples['AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeExcl_Mll10to60'] = ['AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp0Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp1Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp2Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp3Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp4Excl_Mll10to60']
  meta_samples['AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuExcl_Mll10to60'] = ['AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp0Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp1Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp2Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp3Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuNp4Excl_Mll10to60']
  meta_samples['AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauExcl_Mll10to60'] = ['AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp0Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp1Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp2Excl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp3Excl_Mll10to60']#,'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauNp4Excl_Mll10to60']
  meta_samples['McAtNloJimmy_CT10_ttbar'] = ['McAtNloJimmy_CT10_ttbar_LeptonFilter', 'McAtNloJimmy_CT10_ttbar_4LepMass_Mll40GeV12GeV']

  meta_samples['PowhegPythia8_AU2CT10_ZZ_mll4_2pt5'] = ['PowhegPythia8_AU2CT10_ZZ_4e_mll4_2pt5', 'PowhegPythia8_AU2CT10_ZZ_2e2mu_mll4_2pt5', 'PowhegPythia8_AU2CT10_ZZ_2e2tau_mll4_2pt5', 'PowhegPythia8_AU2CT10_ZZ_4mu_mll4_2pt5', 'PowhegPythia8_AU2CT10_ZZ_4tau_mll4_2pt5', 'PowhegPythia8_AU2CT10_ZZ_2mu2tau_mll4_2pt5']
  meta_samples['gg2ZZ_JIMMY_AUET2CT10'] = ['gg2ZZJimmy_AUET2CT10_ZZ4e', 'gg2ZZJimmy_AUET2CT10_ZZ4mu', 'gg2ZZJimmy_AUET2CT10_ZZ2e2mu']

  list_of_ZZ_samples = ['PowhegPythia8_AU2CT10_ZZ_mll4_2pt5', 'gg2ZZ_JIMMY_AUET2CT10']
  list_of_ZZ_samples_nogg2ZZ = ['PowhegPythia8_AU2CT10_ZZ_mll4_2pt5']


### lists of physics samples
physics_samples = {}

if (tag == 'mc11c'):
  physics_samples['ZZ'] =  list_of_ZZ_samples
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
  physics_samples['ZZ'] =  list_of_ZZ_samples
  physics_samples['Z'] =   ['AlpgenJimmy_AUET2CTEQ6L1_Zee', 'AlpgenJimmy_AUET2CTEQ6L1_Zmumu', 'AlpgenJimmy_AUET2CTEQ6L1_Ztautau', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeExcl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZmumuExcl_Mll10to60', 'AlpgenJimmy_Auto_AUET2CTEQ6L1_ZtautauExcl_Mll10to60']
  physics_samples['Zbb'] = ['AlpgenHWfZeebb_4LepM', 'AlpgenHWfZmumubb_4LepM', 'AlpgenHWfZeebb_Veto4LepM_Pass3Lep', 'AlpgenHWfZmumubb_Veto4LepM_Pass3Lep', 'AlpgenJimmyLowMassDYeebb_nofilter', 'AlpgenJimmyLowMassDYmumubb_nofilter', 'AlpgenJimmyLowMassDYtautaubb_nofilter']
  physics_samples['tt'] =  ['McAtNloJimmy_CT10_ttbar']
  physics_samples['125'] = ['PowhegPythia8_AU2CT10_H125']
  physics_samples['130'] = ['PowhegPythia8_AU2CT10_H130']
  physics_samples['150'] = ['PowhegPythia8_AU2CT10_H150']
  physics_samples['190'] = ['PowhegPythia8_AU2CT10_H190']
  physics_samples['200'] = ['PowhegPythia8_AU2CT10_H200']
  physics_samples['400'] = ['PowhegPythia8_AU2CT10_H400']
  physics_samples['600'] = ['PowhegPythia8_AU2CT10_H600']



### YIELD ESTIMATION

raw_yield = {}
final_yield = {} # after the possible mass cut at 100 GeV
final_yield_complete = {}  # BEFORE the possible mass cut at 100 GeV (and in full mass range of course): needed to normalize histograms since they are ALWAYS obtained without the 100 GeV cut
generated = {}
cross_section = {}
for mass_region in mass_regions:
  raw_yield[mass_region] = {}
  final_yield[mass_region] = {}

  for chan in channels:
    raw_yield[mass_region][chan] = {}
    final_yield[mass_region][chan] = {}
for chan in channels:
  final_yield_complete[chan] = {}

data_yield = {}
for mass_region in mass_regions:
  data_yield[mass_region] = {}
  for chan in channels:
    data_yield[mass_region][chan] = 0


### OUTPUT PREPARATION

### prepare output trees
tree_ZZ_weight = numpy.zeros(1, dtype=float)
tree_ZZ_mass = numpy.zeros(1, dtype=float)
tree_ZZ_gg = {}
tree_ZZ_qq = {}

for chan in channels:
  # qq
  tree_ZZ_gg[chan] = TTree('tree_%s' % chan, 'ZZ gg selected events [%s final state]' % chan)
  tree_ZZ_gg[chan].Branch('weight', tree_ZZ_weight, 'weight/D')
  tree_ZZ_gg[chan].Branch('mass', tree_ZZ_mass, 'mass/D')

  # gg
  tree_ZZ_qq[chan] = TTree('tree_%s' % chan, 'ZZ qq selected events [%s final state]' % chan)
  tree_ZZ_qq[chan].Branch('weight', tree_ZZ_weight, 'weight/D')
  tree_ZZ_qq[chan].Branch('mass', tree_ZZ_mass, 'mass/D')

### prepare output trees
tree_redux_weight = numpy.zeros(1, dtype=float)
tree_redux_mass = numpy.zeros(1, dtype=float)
tree_redux_m12 = numpy.zeros(1, dtype=float) # kostas
tree_redux_m34 = numpy.zeros(1, dtype=float) # kostas
tree_redux_type = numpy.zeros(1, dtype=int) # kostas
tree_redux_run = numpy.zeros(1, dtype=int) # kostas
tree_redux = TTree('tree_redux', 'events for reducible backgrounds [all final states together]')
tree_redux.Branch('weight', tree_redux_weight, 'weight/D')
tree_redux.Branch('mass', tree_redux_mass, 'mass/D')
tree_redux.Branch('m12', tree_redux_m12, 'm12/D') # kostas
tree_redux.Branch('m34', tree_redux_m34, 'm34/D') # kostas
tree_redux.Branch('type', tree_redux_type, 'type/I') # kostas
tree_redux.Branch('run', tree_redux_run, 'run/I') # kostas

tree_redux_syst_weight = numpy.zeros(1, dtype=float)
tree_redux_syst_mass = numpy.zeros(1, dtype=float)
tree_redux_syst = TTree('tree_redux_syst', 'events for reducible backgrounds [all final states together, d0 cut]')
tree_redux_syst.Branch('weight', tree_redux_syst_weight, 'weight/D')
tree_redux_syst.Branch('mass', tree_redux_syst_mass, 'mass/D')

tree_redux_syst2_weight = numpy.zeros(1, dtype=float)
tree_redux_syst2_mass = numpy.zeros(1, dtype=float)
tree_redux_syst2 = TTree('tree_redux_syst2', 'events for reducible backgrounds [all final states together, d0 cut, ptcone20/pt cut]')
tree_redux_syst2.Branch('weight', tree_redux_syst2_weight, 'weight/D')
tree_redux_syst2.Branch('mass', tree_redux_syst2_mass, 'mass/D')

### prepare output histos

# binning
bin_width = {}
bin_width['1100'] = 0.500000
bin_width['1105'] = 0.500000
bin_width['1110'] = 0.500000
bin_width['1115'] = 0.500000
bin_width['1120'] = 0.500000
bin_width['1125'] = 0.500000
bin_width['1130'] = 0.500000
bin_width['1135'] = 0.500000
bin_width['1140'] = 0.500000
bin_width['1145'] = 0.500000
bin_width['1150'] = 0.500000
bin_width['1155'] = 0.500000
bin_width['1160'] = 0.500000
bin_width['1165'] = 0.500000
bin_width['1170'] = 0.500000
bin_width['1175'] = 0.500000
bin_width['1180'] = 0.500000
bin_width['1185'] = 0.500000
bin_width['1190'] = 0.500000
bin_width['1195'] = 0.500000
bin_width['1200'] = 0.500000
bin_width['1205'] = 0.500000
bin_width['1210'] = 0.500000
bin_width['1215'] = 0.500000
bin_width['1220'] = 0.500000
bin_width['1225'] = 0.500000
bin_width['1230'] = 0.500000
bin_width['1235'] = 0.500000
bin_width['1240'] = 0.500000
bin_width['1245'] = 0.500000
bin_width['1250'] = 0.500000
bin_width['1255'] = 0.500000
bin_width['1260'] = 0.500000
bin_width['1265'] = 0.500000
bin_width['1270'] = 0.500000
bin_width['1275'] = 0.500000
bin_width['1280'] = 0.500000
bin_width['1285'] = 0.500000
bin_width['1290'] = 0.500000
bin_width['1295'] = 0.500000
bin_width['1300'] = 0.500000
bin_width['1305'] = 0.500000
bin_width['1310'] = 0.500000
bin_width['1315'] = 0.500000
bin_width['1320'] = 0.500000
bin_width['1325'] = 0.500000
bin_width['1330'] = 0.500000
bin_width['1335'] = 0.500000
bin_width['1340'] = 0.500000
bin_width['1345'] = 0.500000
bin_width['1350'] = 0.500000
bin_width['1355'] = 0.500000
bin_width['1360'] = 0.500000
bin_width['1365'] = 0.500000
bin_width['1370'] = 0.500000
bin_width['1375'] = 0.500000
bin_width['1380'] = 0.500000
bin_width['1385'] = 0.500000
bin_width['1390'] = 0.500000
bin_width['1395'] = 0.500000
bin_width['1400'] = 0.500000
bin_width['1405'] = 0.500000
bin_width['1410'] = 0.500000
bin_width['1415'] = 0.500000
bin_width['1420'] = 0.500000
bin_width['1425'] = 0.500000
bin_width['1430'] = 0.500000
bin_width['1435'] = 0.500000
bin_width['1440'] = 0.500000
bin_width['1445'] = 0.500000
bin_width['1450'] = 0.500000
bin_width['1455'] = 0.500000
bin_width['1460'] = 0.500000
bin_width['1465'] = 0.500000
bin_width['1470'] = 0.500000
bin_width['1475'] = 0.500000
bin_width['1480'] = 0.500000
bin_width['1485'] = 0.500000
bin_width['1490'] = 0.500000
bin_width['1495'] = 0.500000
bin_width['1500'] = 0.500000
bin_width['1505'] = 0.500000
bin_width['1510'] = 0.500000
bin_width['1515'] = 0.500000
bin_width['1520'] = 0.500000
bin_width['1525'] = 0.500000
bin_width['1530'] = 0.500000
bin_width['1535'] = 0.500000
bin_width['1540'] = 0.500000
bin_width['1545'] = 0.500000
bin_width['1550'] = 0.500000
bin_width['1555'] = 0.500000
bin_width['1560'] = 0.500000
bin_width['1565'] = 0.500000
bin_width['1570'] = 0.500000
bin_width['1575'] = 0.500000
bin_width['1580'] = 0.500000
bin_width['1585'] = 0.500000
bin_width['1590'] = 0.500000
bin_width['1595'] = 0.500000
bin_width['1600'] = 0.500000
bin_width['1610'] = 0.500000
bin_width['1620'] = 0.500000
bin_width['1630'] = 0.500000
bin_width['1640'] = 0.500000
bin_width['1650'] = 0.500000
bin_width['1660'] = 0.500000
bin_width['1670'] = 0.500000
bin_width['1680'] = 0.500000
bin_width['1690'] = 0.500000
bin_width['1700'] = 0.500000
bin_width['1710'] = 0.500000
bin_width['1720'] = 0.500000
bin_width['1730'] = 0.500000
bin_width['1740'] = 0.500000
bin_width['1750'] = 0.500000
bin_width['1760'] = 0.500000
bin_width['1770'] = 0.500000
bin_width['1780'] = 0.500000
bin_width['1790'] = 0.500000
bin_width['1800'] = 1.000000
bin_width['1810'] = 1.000000
bin_width['1820'] = 1.000000
bin_width['1830'] = 1.000000
bin_width['1840'] = 1.000000
bin_width['1850'] = 1.000000
bin_width['1860'] = 1.000000
bin_width['1870'] = 1.000000
bin_width['1880'] = 1.000000
bin_width['1890'] = 1.000000
bin_width['1900'] = 1.000000
bin_width['1910'] = 1.000000
bin_width['1920'] = 1.000000
bin_width['1930'] = 1.000000
bin_width['1940'] = 1.000000
bin_width['1950'] = 1.000000
bin_width['1960'] = 1.000000
bin_width['1970'] = 1.000000
bin_width['1980'] = 1.000000
bin_width['1990'] = 1.000000
bin_width['2000'] = 1.000000
bin_width['2020'] = 1.000000
bin_width['2040'] = 1.000000
bin_width['2060'] = 1.000000
bin_width['2080'] = 1.000000
bin_width['2100'] = 1.000000
bin_width['2120'] = 1.000000
bin_width['2140'] = 1.000000
bin_width['2160'] = 1.000000
bin_width['2180'] = 1.000000
bin_width['2200'] = 1.000000
bin_width['2220'] = 1.000000
bin_width['2240'] = 1.000000
bin_width['2260'] = 1.000000
bin_width['2280'] = 1.000000
bin_width['2300'] = 1.000000
bin_width['2320'] = 1.000000
bin_width['2340'] = 1.000000
bin_width['2360'] = 1.000000
bin_width['2380'] = 1.000000
bin_width['2400'] = 1.000000
bin_width['2420'] = 1.000000
bin_width['2440'] = 1.000000
bin_width['2460'] = 1.000000
bin_width['2480'] = 1.000000
bin_width['2500'] = 1.000000
bin_width['2520'] = 1.000000
bin_width['2540'] = 1.000000
bin_width['2560'] = 1.000000
bin_width['2580'] = 1.000000
bin_width['2600'] = 1.000000
bin_width['2620'] = 1.000000
bin_width['2640'] = 1.000000
bin_width['2660'] = 1.000000
bin_width['2680'] = 1.000000
bin_width['2700'] = 1.000000
bin_width['2720'] = 1.000000
bin_width['2740'] = 1.000000
bin_width['2760'] = 1.000000
bin_width['2780'] = 1.000000
bin_width['2800'] = 1.000000
bin_width['2820'] = 1.000000
bin_width['2840'] = 1.000000
bin_width['2860'] = 1.000000
bin_width['2880'] = 1.000000
bin_width['2900'] = 1.000000
bin_width['2920'] = 1.000000
bin_width['2940'] = 1.000000
bin_width['2950'] = 1.000000
bin_width['2960'] = 1.000000
bin_width['2980'] = 1.000000
bin_width['3000'] = 2.500000
bin_width['3050'] = 2.500000
bin_width['3100'] = 2.500000
bin_width['3150'] = 2.500000
bin_width['3200'] = 2.500000
bin_width['3250'] = 2.500000
bin_width['3300'] = 2.500000
bin_width['3350'] = 2.500000
bin_width['3400'] = 2.500000
bin_width['3450'] = 2.500000
bin_width['3500'] = 2.500000
bin_width['3550'] = 2.500000
bin_width['3600'] = 2.500000
bin_width['3700'] = 2.500000
bin_width['3800'] = 2.500000
bin_width['3900'] = 2.500000
bin_width['4000'] = 5.000000
bin_width['4200'] = 5.000000
bin_width['4400'] = 5.000000
bin_width['4600'] = 5.000000
bin_width['4800'] = 5.000000
bin_width['5000'] = 10.000000
bin_width['5200'] = 10.000000
bin_width['5400'] = 10.000000
bin_width['5600'] = 10.000000
bin_width['5800'] = 10.000000
bin_width['6000'] = 10.000000

# naming
histo_range = [0., 1000.]
widths = [0.5, 1., 2.5, 5., 10]
histo_names = ['histott', 'histoZbb', 'histoZ', 'histoZZ', 'histogg2ZZ', 'higgs', 'higgsVBF', 'higgsWH', 'higgsZH']

# histograms (for Higgs mass and additional plots like m12, m34)
mass_histo_simple = {}
mass_histo_additional = {}

for histo in histo_names:
  for i in range(0, 4):
    if ('higgs' in histo): # prepare Higgs signal histos with fine binning
      for mass_hypo in bin_width.keys():
        mass_histo_simple['%s_%s_%i' % (mass_hypo, histo, i)] = TH1F('%s_%s_%i' % (mass_hypo, histo, i), '', int((histo_range[1]-histo_range[0])/min(widths)), histo_range[0], histo_range[1])
        mass_histo_simple['%s_%s_ES_lo_%i' % (mass_hypo, histo, i)] = TH1F('%s_%s_ES_lo_%i' % (mass_hypo, histo, i), '', int((histo_range[1]-histo_range[0])/min(widths)), histo_range[0], histo_range[1])
        mass_histo_simple['%s_%s_ES_hi_%i' % (mass_hypo, histo, i)] = TH1F('%s_%s_ES_hi_%i' % (mass_hypo, histo, i), '', int((histo_range[1]-histo_range[0])/min(widths)), histo_range[0], histo_range[1])
        mass_histo_additional['%s_%sm12_%i' % (mass_hypo, histo, i)] = TH1F('%s_%sm12_%i' % (mass_hypo, histo, i), '', 240, 30, 150)
        mass_histo_additional['%s_%sm34_%i' % (mass_hypo, histo, i)] = TH1F('%s_%sm34_%i' % (mass_hypo, histo, i), '', 150,  0, 150)
    else: # prepare generic background histos with fine binning
        mass_histo_simple['%s_%s_%i' % ('finest', histo, i)] = TH1F('%s_%s_%i' % ('finest', histo, i), '', int((histo_range[1]-histo_range[0])/min(widths)), histo_range[0], histo_range[1])
        mass_histo_simple['%s_%s_hi_%i' % ('finest', histo, i)] = TH1F('%s_%s_hi_%i' % ('finest', histo, i), '', int((histo_range[1]-histo_range[0])/min(widths)), histo_range[0], histo_range[1])
        mass_histo_simple['%s_%s_lo_%i' % ('finest', histo, i)] = TH1F('%s_%s_lo_%i' % ('finest', histo, i), '', int((histo_range[1]-histo_range[0])/min(widths)), histo_range[0], histo_range[1])
        mass_histo_additional['%s_%sm12_%i' % ('finest', histo, i)] = TH1F('%s_%sm12_%i' % ('finest', histo, i), '', 240, 30, 150)
        mass_histo_additional['%s_%sm34_%i' % ('finest', histo, i)] = TH1F('%s_%sm34_%i' % ('finest', histo, i), '', 150,  0, 150)



### HISTOGRAM-SAMPLE ASSOCIATION

### prepare a map to tell which histogram each runnumber must fill (sample -> runnumber but runnumber -> [sample1, sample2], when switching from 7 to 8 TeV)
chan_properhisto_map = {}

# list of actual run numbers contained in each background category
samples_within = {}
samples_within['Zbb'] = [y for x in physics_samples['Zbb'] for y in meta_samples[x]]
samples_within['Z'] = [y for x in physics_samples['Z'] for y in meta_samples[x]]
samples_within['tt'] = [y for x in physics_samples['tt'] for y in meta_samples[x]]
if (tag == 'mc11c'):
  samples_within['gg2ZZ'] = [y for y in meta_samples['gg2ZZ_JIMMY']]
elif (tag == 'mc12a'):
  samples_within['gg2ZZ'] = [y for y in meta_samples['gg2ZZ_JIMMY_AUET2CT10']]
samples_within['ZZ'] = [y for x in list_of_ZZ_samples_nogg2ZZ for y in meta_samples[x]]

### list of interesting samples (signals + backgrounds)
list_of_interesting_samples = []
list_of_interesting_runs = []

# add backgrounds and resume their composition
print ''
print '### EXPECTED BACKGROUNDS'
for sample_type in samples_within:
  print '--- samples within %s' % sample_type
  for element in samples_within[sample_type]:
    print element
    list_of_interesting_samples.append(element)
    list_of_interesting_runs.append(sample_number[element])

  print ''
print ''

# make sure all higgs signals are saved
for samplename in sample_number:
  if samplename not in list_of_interesting_samples:
    useful = False

    if (regexp_mass_ggF.search(samplename)):   useful = True
    elif (regexp_mass_VBF.search(samplename)): useful = True
    elif (regexp_mass_WH.search(samplename)):  useful = True
    elif (regexp_mass_ZH.search(samplename)):  useful = True

    if (useful):
      list_of_interesting_samples.append(samplename)
      list_of_interesting_runs.append(sample_number[samplename])


### association between sample and histogram filled by events in this sample

sample_name = {}

# loop over interesting samples
for samplename in list_of_interesting_samples:
  if (samplename not in sample_number.keys()):
    print 'ERROR: skipping %s since no sample/run number association is found in the map!!!' % samplename
    continue
  run = sample_number[samplename]
  sample_name[run] = samplename

  print 'INFO: run %d sample %s' % (run, samplename)

  proper_histo = '' # _0, _1, _2, _3 are appended automatically

  if (samplename in samples_within['Zbb']):
    proper_histo = 'finest_histoZbb'
  elif (samplename in samples_within['Z']):
    proper_histo = 'finest_histoZ'
  elif (samplename in samples_within['ZZ']):
    proper_histo = 'finest_histoZZ'
  elif (samplename in samples_within['gg2ZZ']):
    proper_histo = 'finest_histogg2ZZ'
  elif (samplename in samples_within['tt']):
    proper_histo = 'finest_histott'
  else: # check if it's a signal
    if (regexp_mass_ggF.search(samplename)):
      proper_histo = '%s0_higgs' % (regexp_mass_ggF.search(samplename).group(1))
    elif (regexp_mass_VBF.search(samplename)):
      proper_histo = '%s0_higgsVBF' % (regexp_mass_VBF.search(samplename).group(1))
    elif (regexp_mass_WH.search(samplename)):
      proper_histo = '%s0_higgsWH' % (regexp_mass_WH.search(samplename).group(1))
    elif (regexp_mass_ZH.search(samplename)):
      proper_histo = '%s0_higgsZH' % (regexp_mass_ZH.search(samplename).group(1))
    else:
      print 'ERROR: sample %s (run %d) will not be saved into any histogram (you will have just the yield)' % (samplename, run)

    # remove samples we are not interested in
    if (proper_histo.split('_')[0] not in bin_width.keys()):
      print 'WARNING: sample %s (run %d) has been REMOVED from the histogram list since it is marked as not needed: would have filled histogram \'%s\'' % (samplename, run, proper_histo)
      proper_histo = ''

  # fill the map
  if (run in chan_properhisto_map):
    print 'WARNING: sample name %s (proper_histo=\'%s\') corresponds to run %d which is already in the list: will update the map only if proper_histo != \'\'' % (samplename, proper_histo, run)
    if (proper_histo != ''):
      chan_properhisto_map[run] = proper_histo
      print 'WARNING: ... -> updated'
  else:
    chan_properhisto_map[run] = proper_histo
  print 'INFO: events from run %d will fill histogram %s' % (run, chan_properhisto_map[run])

print ''
print '### CHAN_PROPERHISTO_MAP'
print chan_properhisto_map
print ''
print ''
print ''


### DEAL WITH DATA

### prepare data histograms
data_histo = {}
data_histo_additional = {}

for chan in channels:
  data_histo[chan] = TH1F('finest_histoDATA_%d' % map_chanid_chan[chan], '', int((histo_range[1]-histo_range[0])/min(widths)), histo_range[0], histo_range[1])

  data_histo_additional[chan] = {}
  data_histo_additional[chan]['m12'] = TH1F('finest_histoDATAm12_%d' % map_chanid_chan[chan], '', 240, 30, 150)
  data_histo_additional[chan]['m34'] = TH1F('finest_histoDATAm34_%d' % map_chanid_chan[chan], '', 150, 0, 150)

  data_histo[chan].SetDirectory(0)
  data_histo_additional[chan]['m12'].SetDirectory(0)
  data_histo_additional[chan]['m34'].SetDirectory(0)


### read the data file

data_file = TFile.Open('%s/data_candidates.root' % base_output_dir)

if (data_file):
  t = data_file.Get('data_candidates')
  if (t):
    for entry in range(t.GetEntries()):
      t.GetEntry(entry)

      mass_to_use = t.mass if not use_mass_constraint else t.mass_constrained

      chan = channels[t.type]

      # save data histo (it's cut just afterwards!)
      data_histo[chan].Fill(mass_to_use/1000.)
      data_histo_additional[chan]['m12'].Fill(t.Z1_m/1000.)
      data_histo_additional[chan]['m34'].Fill(t.Z2_m/1000.)

      if (do_apply_100GeV_cut and t.mass < 100000.): continue;

      # save data counter
      in_lowmass_bin = (test_125_plusminus_6 and abs(t.mass - 125000.) < 3000.) or (not test_125_plusminus_6 and t.mass < 160000.);

      if (in_lowmass_bin):
        data_yield['<160'][chan] += 1 
      else:
        data_yield['>160'][chan] += 1 

    print 'INFO: data retrieved from data_candidates.root (make sure it did not contain duplicates)'
  else:
    print 'ERROR: unable to locate data_candidates tree in data_candidates.root'

  data_file.Close()
else:
  print 'ERROR: unable to locate data_candidates.root, no data histo will be provided'



### RUN OVER MONTE-CARLO SIMULATION

### first step: build the generated entries and cross-section maps

print '### BUILDING GENERATED ENTRIES MAP...'

# loop over the files for a first time, retrieve the generated entries and cross-section and save them properly
for filename in open('%s/__TMP_CUTFLOWMAKER__' % base_output_dir):
  print 'VERBOSE: opening %s' % filename.rstrip('\n')
  f = TFile.Open(filename.rstrip('\n'))

  if f:
    h = f.Get('generatedEntriesHisto')
    e = f.Get('cutflow')

    if (h and e):
      # find out which run is it
      if (e.GetEntries() != 0.):
        e.GetEntry(0)
 
        # skip runs which are not interesting
        if (e.run not in list_of_interesting_runs):
    	  its_name = 'name not found' if e.run not in sample_name else sample_name[e.run]
          print 'WARNING: run number %d (%s) does not correspond to an interesting sample, so it is not used!' % (e.run, its_name)
          f.Close()
	  continue

      # initialize the counters the first time each run is found
      if (e.run not in generated.keys()):
	generated[e.run] = 0.
	for mass_region in mass_regions:
          for chan in channels:
  	    raw_yield[mass_region][chan][e.run] = 0.
  	    final_yield[mass_region][chan][e.run] = 0.
	for chan in channels:
 	  final_yield_complete[chan][e.run] = 0.
	
	# retrieve the cross-section
	cross_section[e.run] = -123456;

	# first hypo is background (easier macro to retrieve xsec)
       #bkg_xsec = CrossSections.bkgCrossSection(e.run, COM_energy, False)
        if (tag == 'mc11c'):
          bkg_xsec = CrossSections.GetBkgCrossSection7TeV(e.run, False)
        elif (tag == 'mc12a' and chan_properhisto_map[e.run] == 'finest_histoZbb'): # MISSING ALPGEN MC12
          bkg_xsec = 1.3 * CrossSections.GetBkgCrossSection7TeV(e.run, False)
	elif (tag == 'mc12a'):
          bkg_xsec = CrossSections.GetBkgCrossSection8TeV(e.run, True)

	if (bkg_xsec != -1):
	  # it's a background!
	  cross_section[e.run] = bkg_xsec
	else:
	  # test if it is a signal
	  higgs_xsector = HiggsCrossSection()

	  # look for ggF
	  parse_result = regexp_mass_ggF.search(sample_name[e.run])
	  if (parse_result):
	    sample_mass[e.run] = float(parse_result.group(1))
	    cross_section[e.run] = higgs_xsector.higgs4lxsecGGF(sample_mass[e.run], COM_energy)
	  # look for VBF
	  parse_result = regexp_mass_VBF.search(sample_name[e.run])
	  if (parse_result):
	    sample_mass[e.run] = float(parse_result.group(1))
	    cross_section[e.run] = higgs_xsector.higgs4lxsecVBF(sample_mass[e.run], COM_energy)
	  # look for WH
	  parse_result = regexp_mass_WH.search(sample_name[e.run])
	  if (parse_result):
	    sample_mass[e.run] = float(parse_result.group(1))
	    cross_section[e.run] = higgs_xsector.higgs4lxsecWH(sample_mass[e.run], COM_energy)
	  # look for ZH
	  parse_result = regexp_mass_ZH.search(sample_name[e.run])
	  if (parse_result):
	    sample_mass[e.run] = float(parse_result.group(1))
	    cross_section[e.run] = higgs_xsector.higgs4lxsecZH(sample_mass[e.run], COM_energy)

	if (cross_section[e.run] == -123456): print 'ERROR: unable to retrieve cross section for sample %s (%d)' % (sample_name[e.run], e.run)
	else: print 'INFO: retrieved xsec for sample %s [run %d] (%lf fb)' % (sample_name[e.run], e.run, cross_section[e.run])

       # update the generated entries map
      generated[e.run] += h.GetBinContent(bin_for_generated)
    f.Close()

### compare number of entries expected from AMI and actually put as denominator

# crosscheck: alert if processed less than expected
for run in generated.keys():
  if run not in gen_entries:
    print 'WARNING: [%d] (%s) %lf but AMI is not available !!!' % (run, sample_name[run], generated[run])
  elif (generated[run] != gen_entries[run]):
    print 'WARNING: [%d] (%s) %lf differs from AMI\'s %lf !!!' % (run, sample_name[run], generated[run], gen_entries[run])



### second step: get the event yield and fill the histograms

for filename in open('%s/__TMP_CUTFLOWMAKER__' % base_output_dir):
  print 'VERBOSE: opening %s' % filename.rstrip('\n')
  f = TFile.Open(filename.rstrip('\n'))

  if f:
    t = f.Get('candidates')

    if (t):
      # skip runs which are not interesting
      if (t.GetEntries() != 0.):
        t.GetEntry(0)

        if (t.run not in list_of_interesting_runs):
  	  its_name = 'name not found' if t.run not in sample_name else sample_name[t.run]
          print 'WARNING: run number %d (%s) does not correspond to an interesting sample, so it is not used!' % (t.run, its_name)
          f.Close()
  	  continue


      # initialize utility stuff for this run (will save us time)
      samplename = ''
      proper_histo = ''
      is_missing_xsec = False
      is_irreducible_bkg = False
      is_irreducible_gg_bkg = False
      is_reducible_bkg = False

      # loop over events
      for entry in xrange(t.GetEntries()):
        t.GetEntry(entry)

        # mass to be used (just for plotting, yields are always using unconstrained)
	mass_to_use = t.H_m if not use_mass_constraint else t.H_m_constrained

        # remove Z+jets double-counting
	if (t.hfor == 4): continue # should act on Z+jet only, not on Z+bb

	# optional: apply a 100 GeV mass cut (only to yields)
	must_be_skipped = (do_apply_100GeV_cut and t.H_m < 100000)

        # identify sample nature (sig, bkg, proper histo)
        if (samplename == ''):
	  samplename = sample_name[t.run]
	  proper_histo = chan_properhisto_map[t.run]

	  is_signal = (proper_histo[0:6] != 'finest')
	  is_irreducible_bkg = (proper_histo == 'finest_histoZZ' or proper_histo == 'finest_histogg2ZZ') # it's in chan_properhisto_map that one can select PowHeg+gg2Z against Phythia for irreducible bkg
	  is_irreducible_gg_bkg = (proper_histo == 'finest_histogg2ZZ')
	  is_reducible_bkg = (not is_irreducible_bkg and not is_signal)
	  is_missing_xsec = (t.run not in cross_section.keys())
	  is_Zbb = (proper_histo == 'finest_histoZbb') # MISSING ALPGEN MC12

        # initialize weights for cutflow and histos
	poids = 1.
	simplepoids = 1.


	chan = channels[t.type]

	# use prioritarily the xsec saved in the tree, for Pythia reweighted ZZ (otherwise the proper event weight)
	if (t.run == 109292 or t.run == 109291 or t.run not in cross_section.keys()):
	  if (tag == 'mc11c' or is_Zbb): # MISSING ALPGEN MC12
            poids = lumi[chan] * t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.ggF_weight * t.top_weight * t.powhegbug_weight * t.xsec_weight / generated[t.run]
	  elif (tag == 'mc12a'):
            poids = lumi[chan] * t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.top_weight * t.powhegbug_weight * t.vxz_weight * t.xsec_weight / generated[t.run]
	else:
	  if (tag == 'mc11c' or is_Zbb): # MISSING ALPGEN MC12
	    poids = lumi[chan] * t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.ggF_weight * t.top_weight * t.powhegbug_weight * cross_section[t.run] / generated[t.run]
	  elif (tag == 'mc12a'):
	    poids = lumi[chan] * t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.top_weight * t.powhegbug_weight * t.vxz_weight * cross_section[t.run] / generated[t.run]

	# extract yields, obtain signal and ZZ background histograms
	if (t.selected == 1):

          # fill ZZ trees
	  if (is_irreducible_bkg):
	    tree_ZZ_weight[0] = poids
	    tree_ZZ_mass[0] = mass_to_use/1000.

	    if (is_irreducible_gg_bkg):
	      tree_ZZ_gg[chan].Fill()
	    else:
	      tree_ZZ_qq[chan].Fill()

	  # fill histograms (only signal and ZZ)
	  if (is_signal or is_irreducible_bkg):
	    if (proper_histo != ''): # there are signals (e.g. 105 GeV) we are not interested in, they must be skipped...
	      if (t.type == 0): chan_id = 0     # 4mu   is 0
	      elif (t.type == 1): chan_id = 1   # 2mu2e is 1
	      elif (t.type == 2): chan_id = 3   # 2e2mu is 3
	      elif (t.type == 3): chan_id = 2   # 4e    is 2

              mass_histo_simple['%s_%i' % (proper_histo, chan_id)].Fill(mass_to_use/1000, poids)
	      if (is_signal): # only 110-150
	        if (sample_mass[t.run] >= 110. and sample_mass[t.run] <= 150):
                  mass_histo_simple['%s_ES_lo_%i' % (proper_histo, chan_id)].Fill(mass_to_use/1000 + ES_shift_down[channels[t.type]] * sample_mass[t.run], poids) # H_m + Delta/130 * sample_mass
                  mass_histo_simple['%s_ES_hi_%i' % (proper_histo, chan_id)].Fill(mass_to_use/1000 + ES_shift_up[channels[t.type]] * sample_mass[t.run], poids) # H_m + Delta/130 * sample_mass
              mass_histo_additional['%sm12_%i' % (proper_histo, chan_id)].Fill(t.Z1_m/1000, poids)
              mass_histo_additional['%sm34_%i' % (proper_histo, chan_id)].Fill(t.Z2_m/1000, poids)


	  # update the yield before the possible mass cut (all samples)
	  final_yield_complete[chan][t.run] += poids

	  # update the yields after the possible mass cut (all samples)
	  if (not must_be_skipped):

	    # update the counters
	    simplepoids = 1 if (t.pu_weight > 0) else -1
	    if t.pu_weight == 0: simplepoids = 0

	    raw_yield['all'][chan][t.run] += simplepoids
	    final_yield['all'][chan][t.run] += poids


            # low/high mass bins
            in_lowmass_bin = (test_125_plusminus_6 and abs(t.H_m-125000.)<3000) or (not test_125_plusminus_6 and t.H_m < 160000.)

	    if (in_lowmass_bin):
	      raw_yield['<160'][chan][t.run] += simplepoids
	      final_yield['<160'][chan][t.run] += poids
	    else:
	      raw_yield['>160'][chan][t.run] += simplepoids
	      final_yield['>160'][chan][t.run] += poids
	# (is_selected)
	

        # obtain the shape of reduciblebackgrounds
	if (is_reducible_bkg and proper_histo != ''):

	  # start from DeltaR+J/psiVeto+best cut and fill combining the four final states
	  if (t.closesttozz == 1):
	    # apply the analysis cuts on Z1
	    if (t.Z1_lepplus_m < 100): # electrons
	      ptcone_cut = 0.15
	      etcone_cut_lepplus = 0.30 if (tag == 'mc11c') else 0.20
	      etcone_cut_lepminus = 0.30 if (tag == 'mc11c') else 0.20
	      d0sign_cut = 6.5
	    else: # muons
	      ptcone_cut = 0.15
	      etcone_cut_lepplus = 0.30 if (t.Z1_lepplus_isSA != 1) else 0.15
	      etcone_cut_lepminus = 0.30 if (t.Z1_lepminus_isSA != 1) else 0.15
	      d0sign_cut = 3.5

  	    is_Z1_good = (t.Z1_lepplus_ptcone20_final/t.Z1_lepplus_pt < ptcone_cut and t.Z1_lepminus_ptcone20_final/t.Z1_lepminus_pt < ptcone_cut)
	    is_Z1_good = (is_Z1_good and t.Z1_lepplus_etcone20_final/t.Z1_lepplus_pt < etcone_cut_lepplus and t.Z1_lepminus_etcone20_final/t.Z1_lepminus_pt < etcone_cut_lepminus)
	    is_Z1_good = (is_Z1_good and abs(t.Z1_lepplus_d0/t.Z1_lepplus_d0_sig) < d0sign_cut and abs(t.Z2_lepminus_d0/t.Z2_lepminus_d0_sig) < d0sign_cut)

	    if (is_Z1_good):
	      histos_to_sum_together = ['finest_histott', 'finest_histoZbb', 'finest_histoZ']

	      # first systematics: no cuts on Z2 ("hi" = "high statistics" :) )
              for generic_histo in histos_to_sum_together:
                for chan_id in range(0,4):
                  mass_histo_simple['%s_hi_%i' % (generic_histo, chan_id)].Fill(mass_to_use/1000, poids)
	      
	      # baseline: doubled isolation cuts
	      #if (1): # kostas
	      if (t.Z2_lepplus_ptcone20_final/t.Z2_lepplus_pt < 0.30 and t.Z2_lepminus_ptcone20_final/t.Z2_lepplus_pt < 0.30):
		for genhist_idx, generic_histo in enumerate(histos_to_sum_together): # all the reducible bkgs are combined
		  for chan_id in range(0,4):
                    mass_histo_simple['%s_%i' % (generic_histo, chan_id)].Fill(mass_to_use/1000, poids)

		    # fill the tree only once per candidate, of course
		    if (genhist_idx == 0 and chan_id == 0):
                      tree_redux_weight[0] = t.top_weight #poids
                      tree_redux_mass[0] = mass_to_use
                      tree_redux_m12[0] = t.Z1_m # kostas
                      tree_redux_m34[0] = t.Z2_m # kostas
                      tree_redux_type[0] = t.type # kostas
                      tree_redux_run[0] = t.run # kostas
                      tree_redux.Fill()

		      if (abs(t.Z2_lepplus_d0/t.Z2_lepplus_d0_sig) < 6.5 and abs(t.Z2_lepminus_d0/t.Z2_lepminus_d0_sig) < 6.5):
                        tree_redux_syst_weight[0] = t.top_weight
                        tree_redux_syst_mass[0] = mass_to_use
                        tree_redux_syst.Fill()

	                if (t.Z2_lepplus_ptcone20_final/t.Z2_lepplus_pt < 0.15 and t.Z2_lepminus_ptcone20_final/t.Z2_lepminus_pt < 0.15):
                          tree_redux_syst2_weight[0] = t.top_weight
                          tree_redux_syst2_mass[0] = mass_to_use
                          tree_redux_syst2.Fill()

		    if (not do_apply_100GeV_cut or t.H_m > 100000):
                      mass_histo_additional['%sm12_%i' % (generic_histo, chan_id)].Fill(t.Z1_m/1000, poids)
                      mass_histo_additional['%sm34_%i' % (generic_histo, chan_id)].Fill(t.Z2_m/1000, poids)
	        
	      # second systematics: main cut ("lo" = "low statistics" :) )
	      if (t.Z2_lepplus_ptcone20_final/t.Z2_lepplus_pt < 0.15 and t.Z2_lepminus_ptcone20_final/t.Z2_lepminus_pt < 0.15):
		for generic_histo in histos_to_sum_together: # all the reducible bkgs are combined
		  for chan_id in range(0,4):
                    mass_histo_simple['%s_lo_%i' % (generic_histo, chan_id)].Fill(mass_to_use/1000, poids)
	    # (Z1 is good)
	  # (closesttozz==1)

    f.Close()



### OUTPUT TABLES

### print detailed table
print ''
print ''
print '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
print '                                       sample ||              4mu               |               2mu2e             |              2e2mu              |                4e               ||'
print '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'

# use "generated" just to crosscheck (it is a subset of list_of_interesting_runs)

for run in sorted(generated.iterkeys()):
  print '%-45s ||          %10lf            |           %10lf            |           %10lf            |           %10lf            ||' % (sample_name[run], final_yield['all']['4mu'][run], final_yield['all']['2mu2e'][run], final_yield['all']['2e2mu'][run], final_yield['all']['4e'][run])
  print '%-45s ||  %07.2lf / %8d (%4lf) |  %07.2lf / %8d (%4lf)  |  %07.2lf / %8d (%4lf)  |  %07.2lf / %8d (%4lf)  ||' % ('', raw_yield['all']['4mu'][run], generated[run], raw_yield['all']['4mu'][run]/generated[run], raw_yield['all']['2mu2e'][run], generated[run], raw_yield['all']['2mu2e'][run]/generated[run], raw_yield['all']['2e2mu'][run], generated[run], raw_yield['all']['2e2mu'][run]/generated[run], raw_yield['all']['4e'][run], generated[run], raw_yield['all']['4e'][run]/generated[run])
print '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
  

print ''
print ''
print ''
print ''

### print summary table
print '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'
print '                                       sample ||              4mu               |               2mu2e             |              2e2mu              |                4e               ||'
print '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'

# combine the final yields for each meta_sample
meta_final_yield = {} # yield
for mass_region in mass_regions:
    meta_final_yield[mass_region] = {}
    for chan in channels:
      meta_final_yield[mass_region][chan] = {}
      for metasample in sorted(meta_samples.iterkeys()):
        meta_final_yield[mass_region][chan][metasample] = 0.

# do the same but for the histogram normalization (so that each histogram has a corresponding normalization; needed for reducible bkgs...)
histo_final_integral = {}
for my_name in mass_histo_simple.keys(): # Higgs mass
  histo_final_integral[my_name] = 0.
for my_name in mass_histo_additional.keys(): # other histos
  histo_final_integral[my_name] = 0.

# loop over metasamples, then on each sample, and use them 
for metasample in sorted(meta_samples.iterkeys()):
  for sample in meta_samples[metasample]:
    if sample not in sample_number.keys():
      print 'WARNING: no corresponding run found for sample %s, skipping from contribution computation' % sample
      continue

    run = sample_number[sample]

    # protect against not found samples
    if run not in generated.keys(): # use generated just for consistency checks (could use list_of_interesting_runs)
      print 'WARNING: sample %s (run %d) not found among the read samples, skipping from contribution computation' % (sample, run)
      continue

    # loop over mass regions and channels, and use them
    for mass_region in mass_regions:
      for chan in channels:
        this_contribution = final_yield[mass_region][chan][run]

        # update the final yield of this meta_sample
        meta_final_yield[mass_region][chan][metasample] += this_contribution

	# update the normalization for histograms (used reducible background only: it is before the 100 GeV cut!)
        if (mass_region == 'all'):
	  # normalization is given by the yield before the cut
          this_contribution_complete = final_yield_complete[chan][run]

          proper_histo = chan_properhisto_map[run]

	  if (proper_histo != ''):
            my_name        =    '%s_%d'    % (proper_histo, map_chanid_chan[chan]) # actual histo name
            my_name_lo     = '%s_lo_%d'    % (proper_histo, map_chanid_chan[chan]) # low syst
            my_name_hi     = '%s_hi_%d'    % (proper_histo, map_chanid_chan[chan]) # high syst
            my_name_ES_lo  = '%s_ES_lo_%d' % (proper_histo, map_chanid_chan[chan]) # low ES syst
            my_name_ES_hi  = '%s_ES_hi_%d' % (proper_histo, map_chanid_chan[chan]) # high ES syst
            my_name_m12    = '%sm12_%d'    % (proper_histo, map_chanid_chan[chan]) # m12
            my_name_m34    = '%sm34_%d'    % (proper_histo, map_chanid_chan[chan]) # m34

            histo_final_integral[my_name] += this_contribution_complete
            histo_final_integral[my_name_m12] += this_contribution_complete
            histo_final_integral[my_name_m34] += this_contribution_complete

	    # if it is a reducible background
	    if (my_name_lo in histo_final_integral.keys()):
              histo_final_integral[my_name_lo] += this_contribution_complete
              histo_final_integral[my_name_hi] += this_contribution_complete

	    # if it is a signal sample
	    if (run in sample_mass.keys()):
              histo_final_integral[my_name_ES_lo] += this_contribution_complete
              histo_final_integral[my_name_ES_hi] += this_contribution_complete

  print '%-45s ||          %10lf            |           %10lf            |           %10lf            |           %10lf            ||' % (metasample, meta_final_yield['all']['4mu'][metasample], meta_final_yield['all']['2mu2e'][metasample], meta_final_yield['all']['2e2mu'][metasample], meta_final_yield['all']['4e'][metasample])
print '----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'



### UNCERTAINTY COMPUTATION

# compute uncertainties
# for signal formula is simple,
#   STDDEV_OVER_VALUE(N(physics_sample))^2 = lumi_uncertainty_relative^2 + xsec_uncertainty_relative^2 + systematics_on_muel_relative^2
# for backgrounds formula is
#   VARIANCE[N(physics_sample)] = (lumi_uncertainty_relative * final_integral(physics_sample))^2
#                               + sum_over_samples( (sel_efficiency_relative * final_yield[physics_sample])^2 )
#                               + sum_over_meta_samples( sum_over_samples( sum_over_samples( sel_eff_rel_i * sel_eff_rel_j * final_yield[physics_sample_i] * final_yield[physics_sample_j] ) ) )

for mass_region in mass_regions:
  continue
  for chan in channels:
    for physics_sample in physics_samples:
      result = 0

      lumi_part = 0
      eff_part = 0
      xsec_part = 0
      sel_part = 0 

   
      global_yield = 0
   
      for meta_sample in physics_samples[physics_sample]:
        for sample in meta_samples[meta_sample]:
          # update global yield
          this_yield = meta_final_yield[mass_region][chan][sample]
          global_yield += this_yield

  	  # update selection efficiency part (only for backgrounds)
	  sel_eff = get_bayes_eff_unc_relative(this_yield, generated[sample_number[sample]])
	  sel_part += pow(sel_eff * this_yield, 2)

	  # update [correlated] cross section part (only for backgrounds)
	  

      
      
      # lumi part
      lumi_part = pow(lumi_uncertainty_relative * global_yield, 2)





### LATEX TABLE

# latex table with signal and background expectations
def get_number_with_error(val):
  return '%.2f' % (val)
 #return '%.2f $\\pm %.2f' % (val, sqrt(val))

# ZZ
ltxZZ_low_4mu = sum([meta_final_yield['<160']['4mu'][x] for x in physics_samples['ZZ']])
ltxZZ_high_4mu = sum([meta_final_yield['>160']['4mu'][x] for x in physics_samples['ZZ']])
ltxZZ_low_mix = sum([meta_final_yield['<160']['2e2mu'][x] + meta_final_yield['<160']['2mu2e'][x] for x in physics_samples['ZZ']])
ltxZZ_high_mix = sum([meta_final_yield['>160']['2e2mu'][x] + meta_final_yield['>160']['2mu2e'][x] for x in physics_samples['ZZ']])
ltxZZ_low_4e = sum([meta_final_yield['<160']['4e'][x] for x in physics_samples['ZZ']])
ltxZZ_high_4e = sum([meta_final_yield['>160']['4e'][x] for x in physics_samples['ZZ']])
# Z
ltxZ_low_4mu = sum([meta_final_yield['<160']['4mu'][x] for x in physics_samples['Z']])
ltxZ_high_4mu = sum([meta_final_yield['>160']['4mu'][x] for x in physics_samples['Z']])
ltxZ_low_mix = sum([meta_final_yield['<160']['2e2mu'][x] + meta_final_yield['<160']['2mu2e'][x] for x in physics_samples['Z']])
ltxZ_high_mix = sum([meta_final_yield['>160']['2e2mu'][x] + meta_final_yield['>160']['2mu2e'][x] for x in physics_samples['Z']])
ltxZ_low_4e = sum([meta_final_yield['<160']['4e'][x] for x in physics_samples['Z']])
ltxZ_high_4e = sum([meta_final_yield['>160']['4e'][x] for x in physics_samples['Z']])
# Zbb
ltxZbb_low_4mu = sum([meta_final_yield['<160']['4mu'][x] for x in physics_samples['Zbb']])
ltxZbb_high_4mu = sum([meta_final_yield['>160']['4mu'][x] for x in physics_samples['Zbb']])
ltxZbb_low_mix = sum([meta_final_yield['<160']['2e2mu'][x] + meta_final_yield['<160']['2mu2e'][x] for x in physics_samples['Zbb']])
ltxZbb_high_mix = sum([meta_final_yield['>160']['2e2mu'][x] + meta_final_yield['>160']['2mu2e'][x] for x in physics_samples['Zbb']])
ltxZbb_low_4e = sum([meta_final_yield['<160']['4e'][x] for x in physics_samples['Zbb']])
ltxZbb_high_4e = sum([meta_final_yield['>160']['4e'][x] for x in physics_samples['Zbb']])
# tt
ltxtt_low_4mu = sum([meta_final_yield['<160']['4mu'][x] for x in physics_samples['tt']])
ltxtt_high_4mu = sum([meta_final_yield['>160']['4mu'][x] for x in physics_samples['tt']])
ltxtt_low_mix = sum([meta_final_yield['<160']['2e2mu'][x] + meta_final_yield['<160']['2mu2e'][x] for x in physics_samples['tt']])
ltxtt_high_mix = sum([meta_final_yield['>160']['2e2mu'][x] + meta_final_yield['>160']['2mu2e'][x] for x in physics_samples['tt']])
ltxtt_low_4e = sum([meta_final_yield['<160']['4e'][x] for x in physics_samples['tt']])
ltxtt_high_4e = sum([meta_final_yield['>160']['4e'][x] for x in physics_samples['tt']])
# total reducible
ltxred_low_4mu = ltxZ_low_4mu + ltxZbb_low_4mu + ltxtt_low_4mu
ltxred_high_4mu = ltxZ_high_4mu + ltxZbb_high_4mu + ltxtt_high_4mu
ltxred_low_mix = ltxZ_low_mix + ltxZbb_low_mix + ltxtt_low_mix
ltxred_high_mix = ltxZ_high_mix + ltxZbb_high_mix + ltxtt_high_mix
ltxred_low_4e = ltxZ_low_4e + ltxZbb_low_4e + ltxtt_low_4e
ltxred_high_4e = ltxZ_high_4e + ltxZbb_high_4e + ltxtt_high_4e
# total 
ltxtot_low_4mu = ltxred_low_4mu + ltxZZ_low_4mu
ltxtot_high_4mu = ltxred_high_4mu + ltxZZ_high_4mu
ltxtot_low_mix = ltxred_low_mix + ltxZZ_low_mix
ltxtot_high_mix = ltxred_high_mix + ltxZZ_high_mix
ltxtot_low_4e = ltxred_low_4e + ltxZZ_low_4e
ltxtot_high_4e = ltxred_high_4e + ltxZZ_high_4e
# data
ltxDATA_low_4mu = data_yield['<160']['4mu']
ltxDATA_high_4mu = data_yield['>160']['4mu']
ltxDATA_low_mix = data_yield['<160']['2mu2e'] + data_yield['<160']['2e2mu']
ltxDATA_high_mix = data_yield['>160']['2mu2e'] + data_yield['>160']['2e2mu']
ltxDATA_low_4e = data_yield['<160']['4e']
ltxDATA_high_4e = data_yield['>160']['4e']
# 125
ltx125_all_4mu = sum([meta_final_yield['all']['4mu'][x] for x in physics_samples['125']])
ltx125_all_mix = sum([meta_final_yield['all']['2e2mu'][x] + meta_final_yield['all']['2mu2e'][x] for x in physics_samples['125']])
ltx125_all_4e = sum([meta_final_yield['all']['4e'][x] for x in physics_samples['125']])
# 130
ltx130_all_4mu = sum([meta_final_yield['all']['4mu'][x] for x in physics_samples['130']])
ltx130_all_mix = sum([meta_final_yield['all']['2e2mu'][x] + meta_final_yield['all']['2mu2e'][x] for x in physics_samples['130']])
ltx130_all_4e = sum([meta_final_yield['all']['4e'][x] for x in physics_samples['130']])
# 150
ltx150_all_4mu = sum([meta_final_yield['all']['4mu'][x] for x in physics_samples['150']])
ltx150_all_mix = sum([meta_final_yield['all']['2e2mu'][x] + meta_final_yield['all']['2mu2e'][x] for x in physics_samples['150']])
ltx150_all_4e = sum([meta_final_yield['all']['4e'][x] for x in physics_samples['150']])
# 190
ltx190_all_4mu = sum([meta_final_yield['all']['4mu'][x] for x in physics_samples['190']])
ltx190_all_mix = sum([meta_final_yield['all']['2e2mu'][x] + meta_final_yield['all']['2mu2e'][x] for x in physics_samples['190']])
ltx190_all_4e = sum([meta_final_yield['all']['4e'][x] for x in physics_samples['190']])
# 200
ltx200_all_4mu = sum([meta_final_yield['all']['4mu'][x] for x in physics_samples['200']])
ltx200_all_mix = sum([meta_final_yield['all']['2e2mu'][x] + meta_final_yield['all']['2mu2e'][x] for x in physics_samples['200']])
ltx200_all_4e = sum([meta_final_yield['all']['4e'][x] for x in physics_samples['200']])
# 400
ltx400_all_4mu = sum([meta_final_yield['all']['4mu'][x] for x in physics_samples['400']])
ltx400_all_mix = sum([meta_final_yield['all']['2e2mu'][x] + meta_final_yield['all']['2mu2e'][x] for x in physics_samples['400']])
ltx400_all_4e = sum([meta_final_yield['all']['4e'][x] for x in physics_samples['400']])
# 600
ltx600_all_4mu = sum([meta_final_yield['all']['4mu'][x] for x in physics_samples['600']])
ltx600_all_mix = sum([meta_final_yield['all']['2e2mu'][x] + meta_final_yield['all']['2mu2e'][x] for x in physics_samples['600']])
ltx600_all_4e = sum([meta_final_yield['all']['4e'][x] for x in physics_samples['600']])

latex_text = '''
      \begin{table*}[!htbp]
      \begin{center}
      \caption{The expected number of signal and background events, with their systematic uncertainty, separated into ``Low mass'' ($m_{4\ell}<160\:\gev$) and ``High mass'' ($m_{4\ell}\geq 160\:\gev$) regions. The observed numbers of events are also presented.\label{tab:EventYields}}
      \begin{tabular}{ccccccc}
      \hline\hline
       & \multicolumn{2}{c}{\mmmm}& \multicolumn{2}{c}{\eemm} & \multicolumn{2}{c}{\eeee} \\
       & Low mass & High mass& Low mass & High mass& Low mass & High mass\\
      \hline
      Int. Luminosity\rule[-1mm]{0mm}{4.7mm} & \multicolumn{2}{c}{ \lumiFourMuon{}}& \multicolumn{2}{c}{\lumiTwoMuonTwoElectron{}} & \multicolumn{2}{c}{\lumiFourElectron{}} \\
      \hline
      $\tabscript{ZZ}{(*)}{}$             & %s  & %s  & %s & %s & %s & %s \\
      $Z$                                 & %s  & %s  & %s & %s & %s & %s \\ 
      $Zbb$                               & %s  & %s  & %s & %s & %s & %s \\
      $tt$                                & %s  & %s  & %s & %s & %s & %s \\\hline
      $Z$, $Zb\bar{b}$, and $t\bar{t}$    & %s  & %s  & %s & %s & %s & %s \\\hline
      Total Background                    & %s  & %s  & %s & %s & %s & %s \\\hline 
      Data & %d & %d & %d & %d & %d & %d\\
      \hline 
       $m_{H}=125\:\gev$ &\multicolumn{2}{c}{%s} &\multicolumn{2}{c}{%s}&\multicolumn{2}{c}{%s}\\
       $m_{H}=130\:\gev$ &\multicolumn{2}{c}{%s} &\multicolumn{2}{c}{%s}&\multicolumn{2}{c}{%s}\\
       $m_{H}=150\:\gev$ &\multicolumn{2}{c}{%s} &\multicolumn{2}{c}{%s}&\multicolumn{2}{c}{%s}\\
       $m_{H}=190\:\gev$ &\multicolumn{2}{c}{%s} &\multicolumn{2}{c}{%s}&\multicolumn{2}{c}{%s}\\
       $m_{H}=200\:\gev$ &\multicolumn{2}{c}{%s} &\multicolumn{2}{c}{%s}&\multicolumn{2}{c}{%s} \\
       $m_{H}=400\:\gev$ &\multicolumn{2}{c}{%s} &\multicolumn{2}{c}{%s}&\multicolumn{2}{c}{%s} \\
       $m_{H}=600\:\gev$ &\multicolumn{2}{c}{%s} &\multicolumn{2}{c}{%s}&\multicolumn{2}{c}{%s} \\
      \hline\hline
      \hline\hline
      \end{tabular}
      \end{center}
      \end{table*}
'''.replace('\\', '\\\\').replace('\b', '\\b').replace('\t', '\\t').replace('\r', '\\r') % ( get_number_with_error(ltxZZ_low_4mu), get_number_with_error(ltxZZ_high_4mu), get_number_with_error(ltxZZ_low_mix), get_number_with_error(ltxZZ_high_mix), get_number_with_error(ltxZZ_low_4e), get_number_with_error(ltxZZ_high_4e), get_number_with_error(ltxZ_low_4mu), get_number_with_error(ltxZ_high_4mu), get_number_with_error(ltxZ_low_mix), get_number_with_error(ltxZ_high_mix), get_number_with_error(ltxZ_low_4e), get_number_with_error(ltxZ_high_4e), get_number_with_error(ltxZbb_low_4mu), get_number_with_error(ltxZbb_high_4mu), get_number_with_error(ltxZbb_low_mix), get_number_with_error(ltxZbb_high_mix), get_number_with_error(ltxZbb_low_4e), get_number_with_error(ltxZbb_high_4e), get_number_with_error(ltxtt_low_4mu), get_number_with_error(ltxtt_high_4mu), get_number_with_error(ltxtt_low_mix), get_number_with_error(ltxtt_high_mix), get_number_with_error(ltxtt_low_4e), get_number_with_error(ltxtt_high_4e), get_number_with_error(ltxred_low_4mu), get_number_with_error(ltxred_high_4mu), get_number_with_error(ltxred_low_mix), get_number_with_error(ltxred_high_mix), get_number_with_error(ltxred_low_4e), get_number_with_error(ltxred_high_4e), get_number_with_error(ltxtot_low_4mu), get_number_with_error(ltxtot_high_4mu), get_number_with_error(ltxtot_low_mix), get_number_with_error(ltxtot_high_mix), get_number_with_error(ltxtot_low_4e), get_number_with_error(ltxtot_high_4e), ltxDATA_low_4mu, ltxDATA_high_4mu, ltxDATA_low_mix, ltxDATA_high_mix, ltxDATA_low_4e, ltxDATA_high_4e, get_number_with_error(ltx125_all_4mu), get_number_with_error(ltx125_all_mix), get_number_with_error(ltx125_all_4e), get_number_with_error(ltx130_all_4mu), get_number_with_error(ltx130_all_mix), get_number_with_error(ltx130_all_4e), get_number_with_error(ltx150_all_4mu), get_number_with_error(ltx150_all_mix), get_number_with_error(ltx150_all_4e), get_number_with_error(ltx190_all_4mu), get_number_with_error(ltx190_all_mix), get_number_with_error(ltx190_all_4e), get_number_with_error(ltx200_all_4mu), get_number_with_error(ltx200_all_mix), get_number_with_error(ltx200_all_4e), get_number_with_error(ltx400_all_4mu), get_number_with_error(ltx400_all_mix), get_number_with_error(ltx400_all_4e), get_number_with_error(ltx600_all_4mu), get_number_with_error(ltx600_all_mix), get_number_with_error(ltx600_all_4e))

print ''
print ' %%% LaTeX output %%%'
print latex_text
print ' %%% done %%%'
print ''



### FINALIZATION OF OUTPUT FILES

print 'INFO: saving ZZ and redux trees...'

### save the ZZ trees
for chan in channels:
  outfile_tree_gg = TFile('%s/ZZ_gg_tree_%s.root' % (base_output_dir, chan), 'RECREATE')
  tmp_tree = tree_ZZ_gg[chan].CloneTree()
  tmp_tree.SetDirectory(outfile_tree_gg)
  tmp_tree.SetName('tree')
  outfile_tree_gg.Write()
  outfile_tree_gg.Close()

  outfile_tree_qq = TFile('%s/ZZ_qq_tree_%s.root' % (base_output_dir, chan), 'RECREATE')
  tmp_tree = tree_ZZ_qq[chan].CloneTree()
  tmp_tree.SetDirectory(outfile_tree_qq)
  tmp_tree.SetName('tree')
  outfile_tree_qq.Write()
  outfile_tree_qq.Close()

### save the reducible background trees
outfile_tree_redux = TFile('%s/reduxbkg_tree.root' % (base_output_dir), 'RECREATE')
tmp_tree = tree_redux.CloneTree()
tmp_tree.SetDirectory(outfile_tree_redux)
tmp_tree.SetName('tree')
outfile_tree_redux.Write()
outfile_tree_redux.Close()

outfile_tree_redux_syst = TFile('%s/reduxsyst_bkg_tree.root' % (base_output_dir), 'RECREATE')
tmp_tree = tree_redux_syst.CloneTree()
tmp_tree.SetDirectory(outfile_tree_redux_syst)
tmp_tree.SetName('tree')
outfile_tree_redux_syst.Write()
outfile_tree_redux_syst.Close()

outfile_tree_redux_syst2 = TFile('%s/reduxsyst2_bkg_tree.root' % (base_output_dir), 'RECREATE')
tmp_tree = tree_redux_syst2.CloneTree()
tmp_tree.SetDirectory(outfile_tree_redux_syst2)
tmp_tree.SetName('tree')
outfile_tree_redux_syst2.Write()
outfile_tree_redux_syst2.Close()

 
### finalize the output histos

print 'INFO: finalizing output histos...'

outfile = {}
h = {}

# loop on each mass file
for mass_hypo in bin_width.keys():
  outfile[mass_hypo] = TFile('%s/workspace_histos/ahistos_%s.root' % (base_output_dir, mass_hypo), 'RECREATE')

  # MC: higgs mass histograms
  for histo in mass_histo_simple.keys():
    if (histo.split('_')[0] == 'finest'): # background sample
      # leave histoZZ and histogg2ZZ intact, rebin wisely and smooth the others
       
      is_irredux = False
      if (histo[0:14] == 'finest_histoZZ' or histo[0:17] == 'finest_histogg2ZZ'):
	is_irredux = True

      if (not is_irredux):
        h[histo] = rebin_and_unbin(mass_histo_simple[histo], 20) # kostas
        h[histo].SetName(histo[7:])
        h[histo].Smooth()
      else:
        h[histo] = mass_histo_simple[histo].Clone(histo[7:])

      # rebin anyways to match sample width
      h[histo].Rebin(int(bin_width[mass_hypo]/h[histo].GetBinWidth(1)))

      # rescale to expectation in the full mass range [no 100 GeV cut]
      # (really needed only for reducible, which is filled in a wider range - histos are always filled before the mass cut!)
      if (h[histo].Integral(0, h[histo].GetNbinsX()+1) != 0):
        h[histo].Scale(histo_final_integral[histo]/h[histo].Integral(0, h[histo].GetNbinsX()+1))

      # remove the low-m bins if appropriate (after the normalizations of course) ### REMOVED since we do it afterwards
     #if (do_apply_100GeV_cut):
     #  remove_bins_below(h[histo], 100.) # this is actually not used since ZZ is taken from smoothing

      # set it to be saved in this sample's file (e.g. ahistos_1300.root)
      h[histo].SetDirectory(outfile[mass_hypo])

    elif (histo.split('_')[0] == mass_hypo): # signal sample: save it in the proper file (note we don't loop on masses to reduce system load due to many open files at the same time)
      h[histo] = mass_histo_simple[histo].Clone(histo[5:])

      # first of all, save a copy of the histogram with finest binning, since it is used by the interpolation
      tmp_histo_signal = mass_histo_simple[histo].Clone()
      tmp_histo_signal.SetDirectory(tmp_outfile)
      tmp_histo_signal.Write()

      # then, rebin to match sample width
      # (this plot is actually not used since the final histograms in data_driven are obtained from 0.5 gev histograms saved above and interpolated)
      h[histo].Rebin(int(bin_width[mass_hypo]/h[histo].GetBinWidth(1)))

      # remove the low-m bins if appropriate (does not really affect plots as it's not really used, one uses the interpolation) ### REMOVED since we do it afterwards
     #if (do_apply_100GeV_cut):
     #  remove_bins_below(h[histo], 100.)

      # set it to be saved in this sample's file (e.g. ahistos_1300.root)
      h[histo].SetDirectory(outfile[mass_hypo])

  # MC: additional histograms
  for histo in mass_histo_additional.keys():
    if (histo.split('_')[0] == 'finest'): # background sample
      # leave histograms intact, just rename
      h[histo] = mass_histo_additional[histo].Clone(histo[7:])

      # rescale to expectation (still, effective mainly for reducible backgrounds)
      if (h[histo].Integral(0, h[histo].GetNbinsX()+1) != 0):
        h[histo].Scale(histo_final_integral[histo]/h[histo].Integral(0, h[histo].GetNbinsX()+1))

      # set it to be saved in this sample's file
      h[histo].SetDirectory(outfile[mass_hypo])

    elif (histo.split('_')[0] == mass_hypo): # signal sample
      h[histo] = mass_histo_additional[histo].Clone(histo[5:])

      # set it to be saved in this sample's file
      h[histo].SetDirectory(outfile[mass_hypo])

  # data (higgs mass): just rebin
  for datahist in data_histo.keys():
    h[datahist] = data_histo[datahist].Clone(data_histo[datahist].GetName()[7:])
    h[datahist].Rebin(int(bin_width[mass_hypo]/h[datahist].GetBinWidth(1)))
    h[datahist].SetDirectory(outfile[mass_hypo])

    # remove the low-m bins if appropriate  ### REMOVED since we do it afterwards
   #if (do_apply_100GeV_cut):
   #  remove_bins_below(h[datahist], 100.)

  # data (additional histograms): just rebin
  for datahist in data_histo_additional.keys():
    h['%s_m12' % datahist] = data_histo_additional[datahist]['m12'].Clone(data_histo_additional[datahist]['m12'].GetName()[7:])
    h['%s_m12' % datahist].SetDirectory(outfile[mass_hypo])
    h['%s_m34' % datahist] = data_histo_additional[datahist]['m34'].Clone(data_histo_additional[datahist]['m34'].GetName()[7:])
    h['%s_m34' % datahist].SetDirectory(outfile[mass_hypo])

  outfile[mass_hypo].Write()
  outfile[mass_hypo].Close()


### save the signal mass histograms
tmp_outfile.Write()
tmp_outfile.Close()
