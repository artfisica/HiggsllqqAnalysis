# checks the change in yield over signal/background samples due to a different cut on electron etcone20
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import glob

gROOT.SetBatch(True)

class sample:
  def __init__(self, name, filelist):
    self.name = name
    self.filelist = filelist

#Etcone20
sample_H130_Etcone20 = sample('H130_Etcone20', glob.glob('user.vippolit.004302*root*'))
sample_H120_Etcone20 = sample('H120_Etcone20', glob.glob('user.vippolit.004304*root*'))
sample_Zee_Etcone20  = sample('Zee_Etcone20', glob.glob('user.vippolit.004306*root*'))
#topoEtcone20
sample_H130_topoEtcone20 = sample('H130_topoEtcone20', glob.glob('user.vippolit.004314*root*'))
sample_H120_topoEtcone20 = sample('H120_topoEtcone20', glob.glob('user.vippolit.004316*root*'))
sample_Zee_topoEtcone20  = sample('Zee_topoEtcone20', glob.glob('user.vippolit.004318*root*'))
# topoEtcone20, official production mc12
sample_ggH125_topoEtcone20 = sample('ggH125_topoEtcone20', glob.glob('output.mc12_ggh125.root'))

new_etcone20_pt_cut = 0.2

samples = [
#          sample_H130_Etcone20,
#          sample_H120_Etcone20,
#   sample_Zee_Etcone20,
#          sample_H130_topoEtcone20,
#          sample_H120_topoEtcone20,
#   sample_Zee_topoEtcone20,
           sample_ggH125_topoEtcone20,
	  ]

channels = ['4mu', '2mu2e', '2e2mu', '4e']

num = {}
num_new = {}
den = {}

for chan in channels:
  num[chan] = {}
  num_new[chan] = {}
  den[chan] = {}

  for sample in samples:
    num[chan][sample.name] = 0.
    num_new[chan][sample.name] = 0.
    den[chan][sample.name] = 0.

for sample in samples:
  for filename in sample.filelist:
    print filename
    f = TFile.Open(filename)

    t_cf = f.Get('cutflow')
    for chan in channels:
      den[chan][sample.name] += t_cf.GetEntries('analysis==0')

    t_cand = f.Get('candidates')
    for entry in range(t_cand.GetEntries()):
      t_cand.GetEntry(entry)

      poids = 1.

      if (t_cand.selected) == 1:
        num[channels[t_cand.type]][sample.name] = num[channels[t_cand.type]][sample.name] + poids

        Z1_lepplus_OK = (t_cand.Z1_lepplus_etcone20_final/t_cand.Z1_lepplus_pt < new_etcone20_pt_cut)  or t_cand.Z1_lepplus_m > 100
        Z1_lepminus_OK = (t_cand.Z1_lepminus_etcone20_final/t_cand.Z1_lepminus_pt < new_etcone20_pt_cut)  or t_cand.Z1_lepminus_m > 100
        Z2_lepplus_OK = (t_cand.Z2_lepplus_etcone20_final/t_cand.Z2_lepplus_pt < new_etcone20_pt_cut)  or t_cand.Z2_lepplus_m > 100
        Z2_lepminus_OK = (t_cand.Z2_lepminus_etcone20_final/t_cand.Z2_lepminus_pt < new_etcone20_pt_cut)  or t_cand.Z2_lepminus_m > 100

        if (Z1_lepplus_OK and Z1_lepminus_OK and Z2_lepplus_OK and Z2_lepminus_OK):
          num_new[channels[t_cand.type]][sample.name] = num_new[channels[t_cand.type]][sample.name] + poids
        
    f.Close()

print ''
print ''

for chan in channels:
  for sample in samples:
    print '[%s] %s: eff = %lf *** (gen/old/new) = (%.0f/%.0f/%.0f)' % (chan, sample.name, num_new[chan][sample.name]/num[chan][sample.name], den[chan][sample.name], num[chan][sample.name], num_new[chan][sample.name])
