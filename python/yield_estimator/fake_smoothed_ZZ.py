# places non-smoothed ZZ histograms in the files where smoothed are expected to be
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *

file_input = 'no_constraint_no100gevcut/MC_only/workspace_histos/ahistos_1300.root'

new_name = 'histoZZSmooth_%d'
old_name_qq = 'histoZZ_%d'
old_name_gg = 'histogg2ZZ_%d'

filename_smoothed_qqZZ = 'histos_qqZZSmooth_SR_2012.root'
filename_smoothed_ggZZ = 'histos_ggZZSmooth_SR_2012.root'

f_smoothed_qqZZ = TFile(filename_smoothed_qqZZ, 'RECREATE')
f_smoothed_ggZZ = TFile(filename_smoothed_ggZZ, 'RECREATE')

f_input = TFile(file_input)

qq = {}
qq_new = {}
gg = {}
gg_new = {}

for i in range(0, 4):
  qq[i] = f_input.Get(old_name_qq % i)
  f_smoothed_qqZZ.cd()
  qq_new[i] = qq[i].Clone(new_name %i)
  qq_new[i].SetDirectory(f_smoothed_qqZZ)
  qq_new[i].Write()

  gg[i] = f_input.Get(old_name_gg % i)
  f_smoothed_ggZZ.cd()
  gg_new[i] = gg[i].Clone(new_name %i)
  gg_new[i].SetDirectory(f_smoothed_ggZZ)
  gg_new[i].Write()
f_input.Close()

f_smoothed_qqZZ.Write()
f_smoothed_ggZZ.Write()

f_smoothed_qqZZ.Close()
f_smoothed_ggZZ.Close()
