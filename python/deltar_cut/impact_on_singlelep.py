# perform lepton studies (uses leptontree generated from IsolationStudies) to check the efficiency of a deltar cut on
# leptons in Z + X -> (l+l-) + e + Y, with the electron coming from FSR
# author: Valerio Ippolito <valerio.ippolito@cern.ch>

from ROOT import *
import os

list_of_files = []
the_dir='atlaslocalgroupdisk/user/vippolit/Iso000004'
os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "Iso000004.*PythiaZmumu_no_filter" > __TMP_CUTFLOWMAKER__' % the_dir)
#os.system('./Higgs4lepAnalysis/python/cutflow_maker/list_dpm_content.sh %s | grep -E "PythiaZee_no_filter.*deltar" > __TMP_CUTFLOWMAKER__' % the_dir)
#os.system('ls -d /storage/data1/ippolitv/dr_studies/*/* > __TMP_CUTFLOWMAKER__')

do_Zmumu_plus_e = True # if False, does Zee plus e



if (do_Zmumu_plus_e):
  outfile = TFile('output_dr_studies_Zmumu.root', 'RECREATE')
else:
  outfile = TFile('output_dr_studies_Zee.root', 'RECREATE')

histo = {}
histo['deltar_fsr'] = TH1F('deltar_fsr', 'deltar between FSR electron and closest muon from Z;FSR electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['deltar_had'] = TH1F('deltar_had', 'deltar between HAD electron and closest muon from Z;HAD electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['deltar_uniso'] = TH1F('deltar_uniso', 'deltar between UNISO electron and closest muon from Z;UNISO electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['deltar_all'] = TH1F('deltar_all', 'deltar between generic electron and closest muon from Z;electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['deltar_afteriso_fsr'] = TH1F('deltar_afteriso_fsr', 'deltar between FSR isolated electron and closest muon from Z;FSR isolated electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['deltar_afteriso_had'] = TH1F('deltar_afteriso_had', 'deltar between HAD isolated electron and closest muon from Z;HAD isolated electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['deltar_afteriso_uniso'] = TH1F('deltar_afteriso_uniso', 'deltar between UNISO isolated electron and closest muon from Z;UNISO isolated electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['deltar_afteriso_all'] = TH1F('deltar_afteriso_all', 'deltar between generic isolated electron and closest muon from Z;isolated electron #DeltaR with closest Z#rightarrow#mu#mu muon;entries', 100, 0, 5)
histo['Et_fsr'] = TH1F('Et_fsr', 'Et of FSR electrons;FSR electron E_{T} [GeV];entries', 200, 0, 100)
histo['Et_had'] = TH1F('Et_had', 'Et of HAD electrons;HAD electron E_{T} [GeV];entries', 200, 0, 100)
histo['Et_uniso'] = TH1F('Et_uniso', 'Et of UNISO electrons;UNISO electron E_{T} [GeV];entries', 200, 0, 100)
histo['Et_all'] = TH1F('Et_all', 'Et of electrons;electron E_{T} [GeV];entries', 200, 0, 100)
histo['Et_afteriso_fsr'] = TH1F('Et_afteriso_fsr', 'Et of isolated FSR electrons;isolated FSR electron E_{T} [GeV];entries', 200, 0, 100)
histo['Et_afteriso_had'] = TH1F('Et_afteriso_had', 'Et of isolated HAD electrons;isolated HAD electron E_{T} [GeV];entries', 200, 0, 100)
histo['Et_afteriso_uniso'] = TH1F('Et_afteriso_uniso', 'Et of isolated UNISO electrons;isolated UNISO electron E_{T} [GeV];entries', 200, 0, 100)
histo['Et_afteriso_all'] = TH1F('Et_afteriso_all', 'Et of isolated electrons;isolated electron E_{T} [GeV];entries', 200, 0, 100)

eff = {}
eff['deltar_fsr'] = TEfficiency('eff_deltarcut_fsr_Et', 'efficiency of deltar cut among FSR electron and closest muon, as a function of its Et;electron E_{T} [GeV];efficiency of #DeltaR > 0.2 cut', 100, 0, 100)
eff['deltar_had'] = TEfficiency('eff_deltarcut_had_Et', 'efficiency of deltar cut among HAD electron and closest muon, as a function of its Et;electron E_{T} [GeV];efficiency of #DeltaR > 0.2 cut', 100, 0, 100)
eff['deltar_uniso'] = TEfficiency('eff_deltarcut_uniso_Et', 'efficiency of deltar cut among UNISO electron and closest muon, as a function of its Et;electron E_{T} [GeV];efficiency of #DeltaR > 0.2 cut', 100, 0, 100)
eff['deltar_all'] = TEfficiency('eff_deltarcut_all_Et', 'efficiency of deltar cut among electron and closest muon, as a function of its Et;electron E_{T} [GeV];efficiency of #DeltaR > 0.2 cut', 100, 0, 100)

counter = {}
counter['fsr'] = {}
counter['had'] = {}
counter['uniso'] = {}
counter['all'] = {}
counter['all']['all'] = 0.
counter['all']['dr_more_02'] = 0.
counter['all']['dr_more_01'] = 0.
counter['fsr']['all'] = 0.
counter['fsr']['dr_more_02'] = 0.
counter['fsr']['dr_more_01'] = 0.
counter['had']['all'] = 0.
counter['had']['dr_more_02'] = 0.
counter['had']['dr_more_01'] = 0.
counter['uniso']['all'] = 0.
counter['uniso']['dr_more_02'] = 0.
counter['uniso']['dr_more_01'] = 0.
counter['afteriso_fsr'] = {}
counter['afteriso_had'] = {}
counter['afteriso_uniso'] = {}
counter['afteriso_all'] = {}
counter['afteriso_all']['all'] = 0.
counter['afteriso_all']['dr_more_02'] = 0.
counter['afteriso_all']['dr_more_01'] = 0.
counter['afteriso_fsr']['all'] = 0.
counter['afteriso_fsr']['dr_more_02'] = 0.
counter['afteriso_fsr']['dr_more_01'] = 0.
counter['afteriso_had']['all'] = 0.
counter['afteriso_had']['dr_more_02'] = 0.
counter['afteriso_had']['dr_more_01'] = 0.
counter['afteriso_uniso']['all'] = 0.
counter['afteriso_uniso']['dr_more_02'] = 0.
counter['afteriso_uniso']['dr_more_01'] = 0.

#for filename in list_of_files:
for filename in open('__TMP_CUTFLOWMAKER__'):
  print 'opening %s' % filename
  f = TFile.Open(filename.rstrip('\n'))

  if (f):
    t = f.Get('leptontree')

    if (t):
      for entry in range(t.GetEntries()):
        t.GetEntry(entry)

	# flavor 1 is muon, 0 is electron
	good_for_Zmumu_plus_e = (t.lep_Z_flavor == 1 and t.lep_m < 100)
	good_for_Zee_plus_e = (t.lep_Z_flavor == 0 and t.lep_m < 100)

	if ((do_Zmumu_plus_e and good_for_Zmumu_plus_e) or (not do_Zmumu_plus_e and good_for_Zee_plus_e)):
	  is_fsr = (t.lep_type == 4 and t.lep_origin == 5 and t.lep_originbkg == 40)
	  is_had = (t.lep_type == 17)
	  is_uniso = (t.lep_type == 3)

	  # rel. 17.0
	  is_isolated = (t.lep_ptcone20_final/t.lep_pt < 0.15 and t.lep_etcone20_final/t.lep_pt < 0.3)

	  el = TLorentzVector()
	  el.SetPtEtaPhiM(t.lep_pt, t.lep_eta, t.lep_phi, t.lep_m)
	  mu_plus = TLorentzVector()
	  mu_plus.SetPtEtaPhiM(t.lep_Z_lepplus_pt, t.lep_Z_lepplus_eta, t.lep_Z_lepplus_phi, t.lep_Z_lepplus_m)
	  mu_minus = TLorentzVector()
	  mu_minus.SetPtEtaPhiM(t.lep_Z_lepminus_pt, t.lep_Z_lepminus_eta, t.lep_Z_lepminus_phi, t.lep_Z_lepminus_m)

	  its_dr = min(mu_plus.DeltaR(el), mu_minus.DeltaR(el))

  	  histo['Et_all'].Fill(el.Et()/1000.)
	  if (is_isolated):
  	    histo['Et_afteriso_all'].Fill(el.Et()/1000.)
	    histo['deltar_afteriso_all'].Fill(its_dr)
	  histo['deltar_all'].Fill(its_dr)
	  eff['deltar_all'].Fill(its_dr > 0.2, el.Et()/1000.)
	  if (is_fsr):
  	    histo['Et_fsr'].Fill(el.Et()/1000.)
	    if (is_isolated):
  	      histo['Et_afteriso_fsr'].Fill(el.Et()/1000.)
  	      histo['deltar_afteriso_fsr'].Fill(its_dr)
  	    histo['deltar_fsr'].Fill(its_dr)
  	    eff['deltar_fsr'].Fill(its_dr > 0.2, el.Et()/1000.)
	  if (is_had):
  	    histo['Et_had'].Fill(el.Et()/1000.)
	    if (is_isolated):
  	      histo['Et_afteriso_had'].Fill(el.Et()/1000.)
  	      histo['deltar_afteriso_had'].Fill(its_dr)
  	    histo['deltar_had'].Fill(its_dr)
  	    eff['deltar_had'].Fill(its_dr > 0.2, el.Et()/1000.)
	  if (is_uniso):
  	    histo['Et_uniso'].Fill(el.Et()/1000.)
	    if (is_isolated):
  	      histo['Et_afteriso_uniso'].Fill(el.Et()/1000.)
  	      histo['deltar_afteriso_uniso'].Fill(its_dr)
  	    histo['deltar_uniso'].Fill(its_dr)
  	    eff['deltar_uniso'].Fill(its_dr > 0.2, el.Et()/1000.)

	  counter['all']['all'] += 1.
	  if (its_dr > 0.2): counter['all']['dr_more_02'] += 1.
	  if (its_dr > 0.1): counter['all']['dr_more_01'] += 1.
	  if (is_fsr):
	    counter['fsr']['all'] += 1.
	    if (its_dr > 0.2): counter['fsr']['dr_more_02'] += 1.
	    if (its_dr > 0.1): counter['fsr']['dr_more_01'] += 1.
	  if (is_had):
	    counter['had']['all'] += 1.
	    if (its_dr > 0.2): counter['had']['dr_more_02'] += 1.
	    if (its_dr > 0.1): counter['had']['dr_more_01'] += 1.
	  if (is_uniso):
	    counter['uniso']['all'] += 1.
	    if (its_dr > 0.2): counter['uniso']['dr_more_02'] += 1.
	    if (its_dr > 0.1): counter['uniso']['dr_more_01'] += 1.

	  if (is_isolated):
	    counter['afteriso_all']['all'] += 1.
	    if (its_dr > 0.2): counter['afteriso_all']['dr_more_02'] += 1.
	    if (its_dr > 0.1): counter['afteriso_all']['dr_more_01'] += 1.
	    if (is_fsr):
	      counter['afteriso_fsr']['all'] += 1.
	      if (its_dr > 0.2): counter['afteriso_fsr']['dr_more_02'] += 1.
	      if (its_dr > 0.1): counter['afteriso_fsr']['dr_more_01'] += 1.
	    if (is_had):
	      counter['afteriso_had']['all'] += 1.
	      if (its_dr > 0.2): counter['afteriso_had']['dr_more_02'] += 1.
	      if (its_dr > 0.1): counter['afteriso_had']['dr_more_01'] += 1.
	    if (is_uniso):
	      counter['afteriso_uniso']['all'] += 1.
	      if (its_dr > 0.2): counter['afteriso_uniso']['dr_more_02'] += 1.
	      if (its_dr > 0.1): counter['afteriso_uniso']['dr_more_01'] += 1.


    else:
      print 'unable to find leptontree in %s' % filename
  else:
    print 'unable to open %s' % filename

outfile.Write()
outfile.Close()


for motherkey in counter.keys():
  print '[%s]' % motherkey
  for key in counter[motherkey].keys():
    print '%s: %lf events (%lf)' % (key, counter[motherkey][key], counter[motherkey][key]/counter[motherkey]['all'])
