# make efficiency and performance plots comparing electron ptcone20_final, etcone20, etcone20_final and topoetcone20_final and same for 0.4 dr cone
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import glob
import numpy

output_list = glob.glob('user.vippolit.004324*root*')
input_file_tag = 'Zee'
#inlist = open('files_iso_%s.txt' % input_file_tag)
#output_list = inlist.readlines()
#inlist.close()

outfile = TFile('output_electron_topoiso_%s.root' % input_file_tag, 'RECREATE')

rel_cut_values = numpy.linspace(0., 0.5, 10)
abs_cut_values = numpy.linspace(0, 4000, 10)
rel_cut_values_full = numpy.linspace(0., 5., 100)
abs_cut_values_full = numpy.linspace(-3000, 50000, 200)

observables = [
               'etcone20', 'etcone20_final', 'topoetcone20_final', 'ptcone20_final',
               'etcone20_pt', 'etcone20_final_pt', 'topoetcone20_final_pt', 'ptcone20_final_pt',
               'etcone40', 'etcone40_final', 'topoetcone40_final', 'ptcone40_final',
               'etcone40_pt', 'etcone40_final_pt', 'topoetcone40_final_pt', 'ptcone40_final_pt',
	      ]

# raw distribution for each observable
histo_obs = {}
# observable vs pileup distributions
histo_obs_vs_bcid = {}
histo_obs_vs_nvx = {}
histo_obs_vs_averageIntPerXing = {}
histo_obs_vs_actualIntPerXing = {}
# <observable> vs pileup distributions
profile_obs_vs_bcid = {}
profile_obs_vs_nvx = {}
profile_obs_vs_averageIntPerXing = {}
profile_obs_vs_actualIntPerXing = {}

# efficiency computation for each observable
den = {}
num = {}

for obs in observables:
  den[obs] = {}
  num[obs] = {}

  if (obs.endswith('_pt')):
    proper_cut_values = rel_cut_values
    proper_cut_values_for_range = rel_cut_values_full
  else:
    proper_cut_values = abs_cut_values
    proper_cut_values_for_range = abs_cut_values_full

  histo_obs[obs] = TH1F('histo_%s' % obs, 'histo_%s' % obs, 100, min(proper_cut_values_for_range), max(proper_cut_values))
  histo_obs_vs_bcid[obs] = TH2F('histo_%s_vs_bcid' % obs, 'histo_%s_vs_bcid' % obs, 3000, 0, 3000, 100, min(proper_cut_values_for_range), max(proper_cut_values))
  histo_obs_vs_nvx[obs] = TH2F('histo_%s_vs_nvx' % obs, 'histo_%s_vs_nvx' % obs, 50, 0, 50, 100, min(proper_cut_values_for_range), max(proper_cut_values))
  histo_obs_vs_averageIntPerXing[obs] = TH2F('histo_%s_vs_averageIntPerXing' % obs, 'histo_%s_vs_averageIntPerXing' % obs, 50, 0, 50, 100, min(proper_cut_values_for_range), max(proper_cut_values))
  histo_obs_vs_actualIntPerXing[obs] = TH2F('histo_%s_vs_actualIntPerXing' % obs, 'histo_%s_vs_actualIntPerXing' % obs, 50, 0, 50, 100, min(proper_cut_values_for_range), max(proper_cut_values))

  for cut_val in proper_cut_values:
    den[obs][cut_val] = 0.
    num[obs][cut_val] = 0.

for filename in output_list:
  f = TFile.Open(filename.rstrip('\n'))
  print f.GetName()
  t = f.Get('leptontree')

  if t:
    for entry in range(t.GetEntries()):
      t.GetEntry(entry)

      poids = 1#5. * t.pu_weight * t.lep_weight * t.ggF_weight * t.xsec_weight / t.processed_entries

      if (t.lep_m < 100): # electron
        histo_obs['ptcone20_final'].Fill(t.lep_ptcone20_final)
        histo_obs['etcone20'].Fill(t.lep_etcone20)
        histo_obs['etcone20_final'].Fill(t.lep_etcone20_final)
        histo_obs['topoetcone20_final'].Fill(t.lep_topoetcone20_final)
        histo_obs['ptcone20_final_pt'].Fill(t.lep_ptcone20_final/t.lep_pt)
        histo_obs['etcone20_pt'].Fill(t.lep_etcone20/t.lep_pt)
        histo_obs['etcone20_final_pt'].Fill(t.lep_etcone20_final/t.lep_pt)
        histo_obs['topoetcone20_final_pt'].Fill(t.lep_topoetcone20_final/t.lep_pt)
        histo_obs['ptcone40_final'].Fill(t.lep_ptcone40_final)
        histo_obs['etcone40'].Fill(t.lep_etcone40)
        histo_obs['etcone40_final'].Fill(t.lep_etcone40_final)
        histo_obs['topoetcone40_final'].Fill(t.lep_topoetcone40_final)
        histo_obs['ptcone40_final_pt'].Fill(t.lep_ptcone40_final/t.lep_pt)
        histo_obs['etcone40_pt'].Fill(t.lep_etcone40/t.lep_pt)
        histo_obs['etcone40_final_pt'].Fill(t.lep_etcone40_final/t.lep_pt)
        histo_obs['topoetcone40_final_pt'].Fill(t.lep_topoetcone40_final/t.lep_pt)

        histo_obs_vs_bcid['ptcone20_final'].Fill(t.bcid, t.lep_ptcone20_final)
        histo_obs_vs_bcid['etcone20'].Fill(t.bcid, t.lep_etcone20)
        histo_obs_vs_bcid['etcone20_final'].Fill(t.bcid, t.lep_etcone20_final)
        histo_obs_vs_bcid['topoetcone20_final'].Fill(t.bcid, t.lep_topoetcone20_final)
        histo_obs_vs_bcid['ptcone20_final_pt'].Fill(t.bcid, t.lep_ptcone20_final/t.lep_pt)
        histo_obs_vs_bcid['etcone20_pt'].Fill(t.bcid, t.lep_etcone20/t.lep_pt)
        histo_obs_vs_bcid['etcone20_final_pt'].Fill(t.bcid, t.lep_etcone20_final/t.lep_pt)
        histo_obs_vs_bcid['topoetcone20_final_pt'].Fill(t.bcid, t.lep_topoetcone20_final/t.lep_pt)
        histo_obs_vs_nvx['ptcone20_final'].Fill(t.nvx, t.lep_ptcone20_final)
        histo_obs_vs_nvx['etcone20'].Fill(t.nvx, t.lep_etcone20)
        histo_obs_vs_nvx['etcone20_final'].Fill(t.nvx, t.lep_etcone20_final)
        histo_obs_vs_nvx['topoetcone20_final'].Fill(t.nvx, t.lep_topoetcone20_final)
        histo_obs_vs_nvx['ptcone20_final_pt'].Fill(t.nvx, t.lep_ptcone20_final/t.lep_pt)
        histo_obs_vs_nvx['etcone20_pt'].Fill(t.nvx, t.lep_etcone20/t.lep_pt)
        histo_obs_vs_nvx['etcone20_final_pt'].Fill(t.nvx, t.lep_etcone20_final/t.lep_pt)
        histo_obs_vs_nvx['topoetcone20_final_pt'].Fill(t.nvx, t.lep_topoetcone20_final/t.lep_pt)
        histo_obs_vs_averageIntPerXing['ptcone20_final'].Fill(t.averageIntPerXing, t.lep_ptcone20_final)
        histo_obs_vs_averageIntPerXing['etcone20'].Fill(t.averageIntPerXing, t.lep_etcone20)
        histo_obs_vs_averageIntPerXing['etcone20_final'].Fill(t.averageIntPerXing, t.lep_etcone20_final)
        histo_obs_vs_averageIntPerXing['topoetcone20_final'].Fill(t.averageIntPerXing, t.lep_topoetcone20_final)
        histo_obs_vs_averageIntPerXing['ptcone20_final_pt'].Fill(t.averageIntPerXing, t.lep_ptcone20_final/t.lep_pt)
        histo_obs_vs_averageIntPerXing['etcone20_pt'].Fill(t.averageIntPerXing, t.lep_etcone20/t.lep_pt)
        histo_obs_vs_averageIntPerXing['etcone20_final_pt'].Fill(t.averageIntPerXing, t.lep_etcone20_final/t.lep_pt)
        histo_obs_vs_averageIntPerXing['topoetcone20_final_pt'].Fill(t.averageIntPerXing, t.lep_topoetcone20_final/t.lep_pt)
        histo_obs_vs_actualIntPerXing['ptcone20_final'].Fill(t.actualIntPerXing, t.lep_ptcone20_final)
        histo_obs_vs_actualIntPerXing['etcone20'].Fill(t.actualIntPerXing, t.lep_etcone20)
        histo_obs_vs_actualIntPerXing['etcone20_final'].Fill(t.actualIntPerXing, t.lep_etcone20_final)
        histo_obs_vs_actualIntPerXing['topoetcone20_final'].Fill(t.actualIntPerXing, t.lep_topoetcone20_final)
        histo_obs_vs_actualIntPerXing['ptcone20_final_pt'].Fill(t.actualIntPerXing, t.lep_ptcone20_final/t.lep_pt)
        histo_obs_vs_actualIntPerXing['etcone20_pt'].Fill(t.actualIntPerXing, t.lep_etcone20/t.lep_pt)
        histo_obs_vs_actualIntPerXing['etcone20_final_pt'].Fill(t.actualIntPerXing, t.lep_etcone20_final/t.lep_pt)
        histo_obs_vs_actualIntPerXing['topoetcone20_final_pt'].Fill(t.actualIntPerXing, t.lep_topoetcone20_final/t.lep_pt)
        histo_obs_vs_bcid['ptcone40_final'].Fill(t.bcid, t.lep_ptcone40_final)
        histo_obs_vs_bcid['etcone40'].Fill(t.bcid, t.lep_etcone40)
        histo_obs_vs_bcid['etcone40_final'].Fill(t.bcid, t.lep_etcone40_final)
        histo_obs_vs_bcid['topoetcone40_final'].Fill(t.bcid, t.lep_topoetcone40_final)
        histo_obs_vs_bcid['ptcone40_final_pt'].Fill(t.bcid, t.lep_ptcone40_final/t.lep_pt)
        histo_obs_vs_bcid['etcone40_pt'].Fill(t.bcid, t.lep_etcone40/t.lep_pt)
        histo_obs_vs_bcid['etcone40_final_pt'].Fill(t.bcid, t.lep_etcone40_final/t.lep_pt)
        histo_obs_vs_bcid['topoetcone40_final_pt'].Fill(t.bcid, t.lep_topoetcone40_final/t.lep_pt)
        histo_obs_vs_nvx['ptcone40_final'].Fill(t.nvx, t.lep_ptcone40_final)
        histo_obs_vs_nvx['etcone40'].Fill(t.nvx, t.lep_etcone40)
        histo_obs_vs_nvx['etcone40_final'].Fill(t.nvx, t.lep_etcone40_final)
        histo_obs_vs_nvx['topoetcone40_final'].Fill(t.nvx, t.lep_topoetcone40_final)
        histo_obs_vs_nvx['ptcone40_final_pt'].Fill(t.nvx, t.lep_ptcone40_final/t.lep_pt)
        histo_obs_vs_nvx['etcone40_pt'].Fill(t.nvx, t.lep_etcone40/t.lep_pt)
        histo_obs_vs_nvx['etcone40_final_pt'].Fill(t.nvx, t.lep_etcone40_final/t.lep_pt)
        histo_obs_vs_nvx['topoetcone40_final_pt'].Fill(t.nvx, t.lep_topoetcone40_final/t.lep_pt)
        histo_obs_vs_averageIntPerXing['ptcone40_final'].Fill(t.averageIntPerXing, t.lep_ptcone40_final)
        histo_obs_vs_averageIntPerXing['etcone40'].Fill(t.averageIntPerXing, t.lep_etcone40)
        histo_obs_vs_averageIntPerXing['etcone40_final'].Fill(t.averageIntPerXing, t.lep_etcone40_final)
        histo_obs_vs_averageIntPerXing['topoetcone40_final'].Fill(t.averageIntPerXing, t.lep_topoetcone40_final)
        histo_obs_vs_averageIntPerXing['ptcone40_final_pt'].Fill(t.averageIntPerXing, t.lep_ptcone40_final/t.lep_pt)
        histo_obs_vs_averageIntPerXing['etcone40_pt'].Fill(t.averageIntPerXing, t.lep_etcone40/t.lep_pt)
        histo_obs_vs_averageIntPerXing['etcone40_final_pt'].Fill(t.averageIntPerXing, t.lep_etcone40_final/t.lep_pt)
        histo_obs_vs_averageIntPerXing['topoetcone40_final_pt'].Fill(t.averageIntPerXing, t.lep_topoetcone40_final/t.lep_pt)
        histo_obs_vs_actualIntPerXing['ptcone40_final'].Fill(t.actualIntPerXing, t.lep_ptcone40_final)
        histo_obs_vs_actualIntPerXing['etcone40'].Fill(t.actualIntPerXing, t.lep_etcone40)
        histo_obs_vs_actualIntPerXing['etcone40_final'].Fill(t.actualIntPerXing, t.lep_etcone40_final)
        histo_obs_vs_actualIntPerXing['topoetcone40_final'].Fill(t.actualIntPerXing, t.lep_topoetcone40_final)
        histo_obs_vs_actualIntPerXing['ptcone40_final_pt'].Fill(t.actualIntPerXing, t.lep_ptcone40_final/t.lep_pt)
        histo_obs_vs_actualIntPerXing['etcone40_pt'].Fill(t.actualIntPerXing, t.lep_etcone40/t.lep_pt)
        histo_obs_vs_actualIntPerXing['etcone40_final_pt'].Fill(t.actualIntPerXing, t.lep_etcone40_final/t.lep_pt)
        histo_obs_vs_actualIntPerXing['topoetcone40_final_pt'].Fill(t.actualIntPerXing, t.lep_topoetcone40_final/t.lep_pt)

        for obs in observables:
          if (obs.endswith('_pt')):
            proper_list = rel_cut_values
          else:
            proper_list = abs_cut_values

	  if (obs == 'ptcone20_final'): x = t.lep_ptcone20_final
	  elif (obs == 'etcone20'): x = t.lep_etcone20
	  elif (obs == 'etcone20_final'): x = t.lep_etcone20_final
	  elif (obs == 'topoetcone20_final'): x = t.lep_topoetcone20_final
	  elif (obs == 'ptcone20_final_pt'): x = t.lep_ptcone20_final/t.lep_pt
	  elif (obs == 'etcone20_pt'): x = t.lep_etcone20/t.lep_pt
	  elif (obs == 'etcone20_final_pt'): x = t.lep_etcone20_final/t.lep_pt
	  elif (obs == 'topoetcone20_final_pt'): x = t.lep_topoetcone20_final/t.lep_pt
	  elif (obs == 'ptcone40_final'): x = t.lep_ptcone40_final
	  elif (obs == 'etcone40'): x = t.lep_etcone40
	  elif (obs == 'etcone40_final'): x = t.lep_etcone40_final
	  elif (obs == 'topoetcone40_final'): x = t.lep_topoetcone40_final
	  elif (obs == 'ptcone40_final_pt'): x = t.lep_ptcone40_final/t.lep_pt
	  elif (obs == 'etcone40_pt'): x = t.lep_etcone40/t.lep_pt
	  elif (obs == 'etcone40_final_pt'): x = t.lep_etcone40_final/t.lep_pt
	  elif (obs == 'topoetcone40_final_pt'): x = t.lep_topoetcone40_final/t.lep_pt
        
          for val in proper_list:
            den[obs][val] = den[obs][val] + poids
	    if (x < val): 
              num[obs][val] = num[obs][val] + poids

outfile.cd()

for obs in observables:
  profile_obs_vs_bcid[obs] = histo_obs_vs_bcid[obs].ProfileX()
  profile_obs_vs_nvx[obs] = histo_obs_vs_nvx[obs].ProfileX()
  profile_obs_vs_averageIntPerXing[obs] = histo_obs_vs_averageIntPerXing[obs].ProfileX()
  profile_obs_vs_actualIntPerXing[obs] = histo_obs_vs_actualIntPerXing[obs].ProfileX()

# print out

cut_scan = {}
for obs in observables:
  cut_scan[obs] = TGraph()
  cut_scan[obs].SetName('cut_scan_%s' % obs)
  cut_scan[obs].SetTitle('cut_scan_%s' % obs)
  cut_scan[obs].GetXaxis().SetTitle('cut on lepton %s' % obs)
  cut_scan[obs].GetYaxis().SetTitle('efficiency')

for obs in observables:
  i = 0

  if (obs.endswith('_pt')):
    proper_list = rel_cut_values
  else:
    proper_list = abs_cut_values

  for val in proper_list:
    cut_scan[obs].SetPoint(i, val, num[obs][val] / den[obs][val])
    i = i + 1

for my_scan in cut_scan.keys():
  cut_scan[my_scan].Write()

outfile.Write()
outfile.Close()
