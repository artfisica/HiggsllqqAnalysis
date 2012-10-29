# return a d0 cut efficiency plot for each flavor
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import glob
import numpy

#output_list = glob.glob('output_test.root')
output_list = [
                'root://grid-cert-03.roma1.infn.it//dpm/roma1.infn.it/home/atlas/atlaslocalgroupdisk/user/vippolit/H4l000006/user.vippolit.H4l000006.105300.PythiaH130zz4l.valid_17_2.01.120501162718/user.vippolit.004302._00001.output_test.root',
                'root://grid-cert-03.roma1.infn.it//dpm/roma1.infn.it/home/atlas/atlaslocalgroupdisk/user/vippolit/H4l000006/user.vippolit.H4l000006.105300.PythiaH130zz4l.valid_17_2.01.120501162718/user.vippolit.004302._00002.output_test.root',
              ]

d0_values = {}
d0_values['el'] = numpy.linspace(0., 10., 100)
d0_values['mu'] = [3.5]

flavors = d0_values.keys()

den = {}
num = {}

# init counters
dict_list = [den, num]
for dict in dict_list:
  for flavor in flavors:
    dict[flavor] = {}
    for d0_value in d0_values[flavor]:
      dict[flavor][d0_value] = 0.

for filename in output_list:
  f = TFile.Open(filename)
  t = f.Get('candidates')

  if t:
    for entry in range(t.GetEntries()):
      t.GetEntry(entry)

      if (t.last >= 9 and t.type == 3):
        poids = 1#5. * t.pu_weight * t.trigSF_weight * t.Z1_lepplus_weight * t.Z1_lepminus_weight * t.Z2_lepplus_weight * t.Z2_lepminus_weight * t.ggF_weight * t.xsec_weight / t.processed_entries

        lepton_d0 = {}
	lepton_d0['mu'] = []
	lepton_d0['el'] = []

	if (t.Z1_lepplus_m < 100):
	  lepton_d0['el'].append(abs(t.Z1_lepplus_d0/t.Z1_lepplus_d0_sig))
	else:
	  lepton_d0['mu'].append(abs(t.Z1_lepplus_d0/t.Z1_lepplus_d0_sig))

	if (t.Z1_lepminus_m < 100):
	  lepton_d0['el'].append(abs(t.Z1_lepminus_d0/t.Z1_lepminus_d0_sig))
	else:
	  lepton_d0['mu'].append(abs(t.Z1_lepminus_d0/t.Z1_lepminus_d0_sig))

	if (t.Z2_lepplus_m < 100):
	  lepton_d0['el'].append(abs(t.Z2_lepplus_d0/t.Z2_lepplus_d0_sig))
	else:
	  lepton_d0['mu'].append(abs(t.Z2_lepplus_d0/t.Z2_lepplus_d0_sig))

	if (t.Z2_lepminus_m < 100):
	  lepton_d0['el'].append(abs(t.Z2_lepminus_d0/t.Z2_lepminus_d0_sig))
	else:
	  lepton_d0['mu'].append(abs(t.Z2_lepminus_d0/t.Z2_lepminus_d0_sig))

	
	# apply cuts
	pass_cut = {}
	for flavor in d0_values.keys():
	  pass_cut[flavor] = {}
	  for d0_value in d0_values[flavor]:
	    pass_cut[flavor][d0_value] = True

	for flavor in d0_values.keys():
	  for d0 in lepton_d0[flavor]:
	    for d0_value in d0_values[flavor]:
	      if (d0 > d0_value):
	        pass_cut[flavor][d0_value] = False

	for flavor in d0_values.keys():
	  for d0_value in d0_values[flavor]:
	    den[flavor][d0_value] += poids
	    if (pass_cut[flavor][d0_value]):
	      num[flavor][d0_value] += poids
  f.Close()


# print out

outfile = TFile('output_optimization_d0.root', 'RECREATE')

cut_scan = {}
for flavor in flavors:
  cut_scan[flavor] = TGraph()
  cut_scan[flavor].SetName('cut_scan_d0_%s' % flavor)
  cut_scan[flavor].SetTitle('cut_scan_d0_%s' % flavor)
  cut_scan[flavor].GetXaxis().SetTitle('cut on lepton |d_{0}/#sigma(d_{0})|')
  cut_scan[flavor].GetYaxis().SetTitle('efficiency')

for flavor in d0_values.keys():
  i = 0
  for d0_value in d0_values[flavor]:
    if (den[flavor][d0_value] != 0):
      cut_scan[flavor].SetPoint(i, d0_value, num[flavor][d0_value] / den[flavor][d0_value])
      i = i + 1

for my_scan in cut_scan.keys():
  cut_scan[my_scan].Write()

outfile.Write()
outfile.Close()
