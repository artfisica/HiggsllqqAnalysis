#define IsolationEfficiency_cxx
// The class definition in IsolationEfficiency.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("IsolationEfficiency.C")
// Root > T->Process("IsolationEfficiency.C","some options")
// Root > T->Process("IsolationEfficiency.C+")
//

#include "IsolationEfficiency.h"
#include <TH2.h>
#include <TStyle.h>

std::vector<Double_t> cut_vector(Double_t min, Double_t max, UInt_t n)
{
   Double_t step = (max - min) / (Double_t)n;

   std::vector<Double_t> result;

   for (UInt_t i = 0; i < n; i++) {
      result.push_back(min + (Double_t) i * step);
   }

   return result;
}

void IsolationEfficiency::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   Info("Begin", "Configuration:");
   Info("Begin", "    - m_isElectron = %s:", (m_isElectron) ? "kTRUE" : "kFALSE");
   Info("Begin", "    - m_useElectronTruth = %s:", (m_useElectronTruth) ? "kTRUE" : "kFALSE");
   Info("Begin", "    - m_useConversions = %s:", (m_useConversions) ? "kTRUE" : "kFALSE");
   Info("Begin", "    - m_useFakes = %s:", (m_useFakes) ? "kTRUE" : "kFALSE");
   Info("Begin", "    - m_useSignalElectrons = %s:", (m_useSignalElectrons) ? "kTRUE" : "kFALSE");

   Info("Begin", "Creating output file...");

   outfile = new TFile(TString::Format("output_%s_topoiso_%s.root", (m_isElectron) ? "electron" : "muon", GetTag()), "RECREATE");

   Info("Begin", "Done.");
}

void IsolationEfficiency::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   Info("SlaveBegin", "Initialization...");

   /*
   prooffile = new TProofOutputFile(TString::Format("output_electron_topoiso_%s.root", GetTag()), TProofOutputFile::kMerge);
   outfile = prooffile->OpenFile("RECREATE");
   */
 //outfile = new TFile("output_electron_topoiso_%s.root", GetTag(), "RECREATE");
 //outfile->cd();

   observables.clear();

   observables.push_back("etcone20");
   observables.push_back("etcone20_final");
   observables.push_back("topoetcone20_final");
   observables.push_back("ptcone20_final");

   observables.push_back("etcone20_pt");
   observables.push_back("etcone20_final_pt");
   observables.push_back("topoetcone20_final_pt");
   observables.push_back("ptcone20_final_pt");

   observables.push_back("etcone40");
   observables.push_back("etcone40_final");
   observables.push_back("topoetcone40_final");
   observables.push_back("ptcone40_final");

   observables.push_back("etcone40_pt");
   observables.push_back("etcone40_final_pt");
   observables.push_back("topoetcone40_final_pt");
   observables.push_back("ptcone40_final_pt");

 //rel_cut_values = cut_vector(0., 5, 100);
   rel_cut_values = cut_vector(0., 0.5, 10);
   abs_cut_values = cut_vector(0., 4000, 10);
   rel_minmax = std::pair<Double_t, Double_t>(0, 5);
   abs_minmax = std::pair<Double_t, Double_t>(-3000, 50000);

   // cut < 0.15
   eff_vs_nvx["ptcone20_lt_015"] = new TEfficiency("eff_ptcone20_lt_015_vs_nvx", "efficiency of ptcone20/pt < 0.15 cut vs n_{VX};n_{VX};efficiency of ptcone20/p_{T} < 0.15", 50, 0, 50);
   eff_vs_pt["ptcone20_lt_015"] = new TEfficiency("eff_ptcone20_lt_015_vs_pt", "efficiency of ptcone20/pt < 0.15 cut vs p_{T};p_{T};efficiency of ptcone20/p_{T} < 0.15", 100, 0, 100);
   eff_vs_eta["ptcone20_lt_015"] = new TEfficiency("eff_ptcone20_lt_015_vs_eta", "efficiency of ptcone20/pt < 0.15 cut vs #eta;#eta;efficiency of ptcone20/p_{T} < 0.15", 60, -3, 3);
   eff_vs_phi["ptcone20_lt_015"] = new TEfficiency("eff_ptcone20_lt_015_vs_phi", "efficiency of ptcone20/pt < 0.15 cut vs #phi;#phi;efficiency of ptcone20/p_{T} < 0.15", 64, -3.2, 3.2);
   eff_vs_etaphi["ptcone20_lt_015"] = new TEfficiency("eff_ptcone20_lt_015_vs_etaphi", "efficiency of ptcone20/pt < 0.15 cut vs #eta,#phi;#phi;#eta", 64, -3.2, 3.2, 60, -3, 3);

   // cut < 0.2
   eff_vs_nvx["topoetcone20_lt_02"] = new TEfficiency("eff_topoetcone20_lt_02_vs_nvx", "efficiency of topoEtcone20/pt < 0.2 cut vs n_{VX};n_{VX};efficiency of topoEtcone20/p_{T} < 0.2", 50, 0, 50);
   eff_vs_nvx["etcone20_lt_02"] = new TEfficiency("eff_etcone20_lt_02_vs_nvx", "efficiency of etcone20/pt < 0.2 cut vs n_{VX};n_{VX};efficiency of etcone20/p_{T} < 0.2", 50, 0, 50);
   eff_vs_pt["topoetcone20_lt_02"] = new TEfficiency("eff_topoetcone20_lt_02_vs_pt", "efficiency of topoEtcone20/pt < 0.2 cut vs p_{T};p_{T};efficiency of topoEtcone20/p_{T} < 0.2", 100, 0, 100);
   eff_vs_pt["etcone20_lt_02"] = new TEfficiency("eff_etcone20_lt_02_vs_pt", "efficiency of etcone20/pt < 0.2 cut vs p_{T};p_{T};efficiency of etcone20/p_{T} < 0.2", 100, 0, 100);
   eff_vs_eta["topoetcone20_lt_02"] = new TEfficiency("eff_topoetcone20_lt_02_vs_eta", "efficiency of topoEtcone20/pt < 0.2 cut vs #eta;#eta;efficiency of topoEtcone20/p_{T} < 0.2", 60, -3, 3);
   eff_vs_eta["etcone20_lt_02"] = new TEfficiency("eff_etcone20_lt_02_vs_eta", "efficiency of etcone20/pt < 0.2 cut vs #phi;#eta;efficiency of etcone20/p_{T} < 0.2", 60, -3, 3);
   eff_vs_phi["topoetcone20_lt_02"] = new TEfficiency("eff_topoetcone20_lt_02_vs_phi", "efficiency of topoEtcone20/pt < 0.2 cut vs #phi;#phi;efficiency of topoEtcone20/p_{T} < 0.2", 64, -3.2, 3.2);
   eff_vs_phi["etcone20_lt_02"] = new TEfficiency("eff_etcone20_lt_02_vs_phi", "efficiency of etcone20/pt < 0.2 cut vs #phi;#phi;efficiency of etcone20/p_{T} < 0.2", 64, -3.2, 3.2);
   eff_vs_etaphi["etcone20_lt_02"] = new TEfficiency("eff_etcone20_lt_02_vs_etaphi", "efficiency of etcone20/pt < 0.2 cut vs #eta,#phi;#phi;#eta", 64, -3.2, 3.2, 60, -3, 3);
   eff_vs_etaphi["topoetcone20_lt_02"] = new TEfficiency("eff_topoetcone20_lt_02_vs_etaphi", "efficiency of topoetcone20/pt < 0.2 cut vs #eta,#phi;#phi;#eta", 64, -3.2, 3.2, 60, -3, 3);

   // cut < 0.25
   eff_vs_nvx["topoetcone20_lt_025"] = new TEfficiency("eff_topoetcone20_lt_025_vs_nvx", "efficiency of topoEtcone20/pt < 0.25 cut vs n_{VX};n_{VX};efficiency of topoEtcone20/p_{T} < 0.25", 50, 0, 50);
   eff_vs_nvx["etcone20_lt_025"] = new TEfficiency("eff_etcone20_lt_025_vs_nvx", "efficiency of etcone20/pt < 0.25 cut vs n_{VX};n_{VX};efficiency of etcone20/p_{T} < 0.25", 50, 0, 50);
   eff_vs_pt["topoetcone20_lt_025"] = new TEfficiency("eff_topoetcone20_lt_025_vs_pt", "efficiency of topoEtcone20/pt < 0.25 cut vs p_{T};p_{T};efficiency of topoEtcone20/p_{T} < 0.25", 100, 0, 100);
   eff_vs_pt["etcone20_lt_025"] = new TEfficiency("eff_etcone20_lt_025_vs_pt", "efficiency of etcone20/pt < 0.25 cut vs p_{T};p_{T};efficiency of etcone20/p_{T} < 0.25", 100, 0, 100);
   eff_vs_eta["topoetcone20_lt_025"] = new TEfficiency("eff_topoetcone20_lt_025_vs_eta", "efficiency of topoEtcone20/pt < 0.25 cut vs #eta;#eta;efficiency of topoEtcone20/p_{T} < 0.25", 60, -3, 3);
   eff_vs_eta["etcone20_lt_025"] = new TEfficiency("eff_etcone20_lt_025_vs_eta", "efficiency of etcone20/pt < 0.25 cut vs #phi;#eta;efficiency of etcone20/p_{T} < 0.25", 60, -3, 3);
   eff_vs_phi["topoetcone20_lt_025"] = new TEfficiency("eff_topoetcone20_lt_025_vs_phi", "efficiency of topoEtcone20/pt < 0.25 cut vs #phi;#phi;efficiency of topoEtcone20/p_{T} < 0.25", 64, -3.2, 3.2);
   eff_vs_phi["etcone20_lt_025"] = new TEfficiency("eff_etcone20_lt_025_vs_phi", "efficiency of etcone20/pt < 0.25 cut vs #phi;#phi;efficiency of etcone20/p_{T} < 0.25", 64, -3.2, 3.2);
   eff_vs_etaphi["etcone20_lt_025"] = new TEfficiency("eff_etcone20_lt_025_vs_etaphi", "efficiency of etcone20/pt < 0.25 cut vs #eta,#phi;#phi;#eta", 64, -3.2, 3.2, 60, -3, 3);
   eff_vs_etaphi["topoetcone20_lt_025"] = new TEfficiency("eff_topoetcone20_lt_025_vs_etaphi", "efficiency of topoetcone20/pt < 0.25 cut vs #eta,#phi;#phi;#eta", 64, -3.2, 3.2, 60, -3, 3);

   // cut < 0.3
   eff_vs_nvx["topoetcone20_lt_03"] = new TEfficiency("eff_topoetcone20_lt_03_vs_nvx", "efficiency of topoEtcone20/pt < 0.3 cut vs n_{VX};n_{VX};efficiency of topoEtcone20/p_{T} < 0.3", 50, 0, 50);
   eff_vs_nvx["etcone20_lt_03"] = new TEfficiency("eff_etcone20_lt_03_vs_nvx", "efficiency of etcone20/pt < 0.3 cut vs n_{VX};n_{VX};efficiency of etcone20/p_{T} < 0.3", 50, 0, 50);
   eff_vs_pt["topoetcone20_lt_03"] = new TEfficiency("eff_topoetcone20_lt_03_vs_pt", "efficiency of topoEtcone20/pt < 0.3 cut vs p_{T};p_{T};efficiency of topoEtcone20/p_{T} < 0.3", 100, 0, 100);
   eff_vs_pt["etcone20_lt_03"] = new TEfficiency("eff_etcone20_lt_03_vs_pt", "efficiency of etcone20/pt < 0.3 cut vs p_{T};p_{T};efficiency of etcone20/p_{T} < 0.3", 100, 0, 100);
   eff_vs_eta["topoetcone20_lt_03"] = new TEfficiency("eff_topoetcone20_lt_03_vs_eta", "efficiency of topoEtcone20/pt < 0.3 cut vs #eta;#eta;efficiency of topoEtcone20/p_{T} < 0.3", 60, -3, 3);
   eff_vs_eta["etcone20_lt_03"] = new TEfficiency("eff_etcone20_lt_03_vs_eta", "efficiency of etcone20/pt < 0.3 cut vs #phi;#eta;efficiency of etcone20/p_{T} < 0.3", 60, -3, 3);
   eff_vs_phi["topoetcone20_lt_03"] = new TEfficiency("eff_topoetcone20_lt_03_vs_phi", "efficiency of topoEtcone20/pt < 0.3 cut vs #phi;#phi;efficiency of topoEtcone20/p_{T} < 0.3", 64, -3.2, 3.2);
   eff_vs_phi["etcone20_lt_03"] = new TEfficiency("eff_etcone20_lt_03_vs_phi", "efficiency of etcone20/pt < 0.3 cut vs #phi;#phi;efficiency of etcone20/p_{T} < 0.3", 64, -3.2, 3.2);
   eff_vs_etaphi["etcone20_lt_03"] = new TEfficiency("eff_etcone20_lt_03_vs_etaphi", "efficiency of etcone20/pt < 0.3 cut vs #eta,#phi;#phi;#eta", 64, -3.2, 3.2, 60, -3, 3);
   eff_vs_etaphi["topoetcone20_lt_03"] = new TEfficiency("eff_topoetcone20_lt_03_vs_etaphi", "efficiency of topoetcone20/pt < 0.3 cut vs #eta,#phi;#phi;#eta", 64, -3.2, 3.2, 60, -3, 3);

   for (std::vector<TString>::iterator obs = observables.begin(); obs != observables.end(); ++obs) {
      std::vector<Double_t> *proper_cut_values;
      std::pair<Double_t, Double_t> *proper_minmax;

      if (obs->EndsWith("_pt")) {
         proper_cut_values = &rel_cut_values;
         proper_minmax = &rel_minmax;
      } else {
         proper_cut_values = &abs_cut_values;
         proper_minmax = &abs_minmax;
      }

      histo_obs[*obs] = new TH1F(TString::Format("histo_%s", obs->Data()), TString::Format("histo_%s", obs->Data()), 100, proper_minmax->first, proper_minmax->second);
      histo_obs_vs_bcid[*obs] = new TH2F(TString::Format("histo_%s_vs_bcid", obs->Data()), TString::Format("histo_%s_vs_bcid", obs->Data()), 3000, 0, 3000, 100, proper_minmax->first, proper_minmax->second);
      histo_obs_vs_nvx[*obs] = new TH2F(TString::Format("histo_%s_vs_nvx", obs->Data()), TString::Format("histo_%s_vs_nvx", obs->Data()), 50, 0, 50, 100, proper_minmax->first, proper_minmax->second);
      histo_obs_vs_averageIntPerXing[*obs] = new TH2F(TString::Format("histo_%s_vs_averageIntPerXing", obs->Data()), TString::Format("histo_%s_vs_averageIntPerXing", obs->Data()), 3000, 0, 3000, 100, proper_minmax->first, proper_minmax->second);
      histo_obs_vs_actualIntPerXing[*obs] = new TH2F(TString::Format("histo_%s_vs_actualIntPerXing", obs->Data()), TString::Format("histo_%s_vs_actualIntPerXing", obs->Data()), 3000, 0, 3000, 100, proper_minmax->first, proper_minmax->second);

      /*
      fOutput->Add(histo_obs[*obs]);
      fOutput->Add(histo_obs_vs_bcid[*obs]);
      fOutput->Add(histo_obs_vs_nvx[*obs]);
      fOutput->Add(histo_obs_vs_averageIntPerXing[*obs]);
      fOutput->Add(histo_obs_vs_actualIntPerXing[*obs]);
      */

      for (std::vector<Double_t>::iterator cutval = proper_cut_values->begin(); cutval != proper_cut_values->end(); ++cutval) {
         den[*obs][*cutval] = 0.;
         num[*obs][*cutval] = 0.;
      }
   }

   Info("SlaveBegin", "Done.");
}

Bool_t IsolationEfficiency::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either IsolationEfficiency::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   if (entry % 1000 == 0)
     Info("Process", "Processing entry %lld", entry);

   b_lep_m->GetEntry(entry);

   Double_t poids = 1.;

   if ((lep_m < 100 && m_isElectron) || (lep_m >= 100 && !m_isElectron)) {
      // LEPTONS
      b_lep_ptcone20_final->GetEntry(entry);
      b_lep_etcone20->GetEntry(entry);
      b_lep_etcone20_final->GetEntry(entry);
      b_lep_topoetcone20_final->GetEntry(entry);
      b_lep_ptcone40_final->GetEntry(entry);
      b_lep_etcone40->GetEntry(entry);
      b_lep_etcone40_final->GetEntry(entry);
      b_lep_topoetcone40_final->GetEntry(entry);
      b_lep_pt->GetEntry(entry);
      b_lep_eta->GetEntry(entry);
      b_lep_phi->GetEntry(entry);
      b_bcid->GetEntry(entry);
      b_nvx->GetEntry(entry);
      b_averageIntPerXing->GetEntry(entry);
      b_actualIntPerXing->GetEntry(entry);

      if (m_useElectronTruth && lep_m < 100) {
        b_lep_type->GetEntry(entry);
        b_lep_origin->GetEntry(entry);
        b_lep_originbkg->GetEntry(entry);

	Bool_t toBeKept(kFALSE);

	toBeKept = (toBeKept || (m_useSignalElectrons && lep_type == 2));
	toBeKept = (toBeKept || (m_useConversions && lep_type == 4 && lep_origin == 5));
	toBeKept = (toBeKept || (m_useFakes && lep_type == 17));

	if (!toBeKept) return kTRUE;
      }

    //b_lep_d0->GetEntry(entry); // REMOVE ME
    //if (TMath::Abs(lep_d0) > 0.01) return kTRUE; // REMOVE ME


      histo_obs["ptcone20_final"]->Fill(lep_ptcone20_final);
      histo_obs["etcone20"]->Fill(lep_etcone20);
      histo_obs["etcone20_final"]->Fill(lep_etcone20_final);
      histo_obs["topoetcone20_final"]->Fill(lep_topoetcone20_final);
      histo_obs["ptcone20_final_pt"]->Fill(lep_ptcone20_final / lep_pt);
      histo_obs["etcone20_pt"]->Fill(lep_etcone20 / lep_pt);
      histo_obs["etcone20_final_pt"]->Fill(lep_etcone20_final / lep_pt);
      histo_obs["topoetcone20_final_pt"]->Fill(lep_topoetcone20_final / lep_pt);
      histo_obs["ptcone40_final"]->Fill(lep_ptcone40_final);
      histo_obs["etcone40"]->Fill(lep_etcone40);
      histo_obs["etcone40_final"]->Fill(lep_etcone40_final);
      histo_obs["topoetcone40_final"]->Fill(lep_topoetcone40_final);
      histo_obs["ptcone40_final_pt"]->Fill(lep_ptcone40_final / lep_pt);
      histo_obs["etcone40_pt"]->Fill(lep_etcone40 / lep_pt);
      histo_obs["etcone40_final_pt"]->Fill(lep_etcone40_final / lep_pt);
      histo_obs["topoetcone40_final_pt"]->Fill(lep_topoetcone40_final / lep_pt);

      histo_obs_vs_bcid["ptcone20_final"]->Fill(bcid, lep_ptcone20_final);
      histo_obs_vs_bcid["etcone20"]->Fill(bcid, lep_etcone20);
      histo_obs_vs_bcid["etcone20_final"]->Fill(bcid, lep_etcone20_final);
      histo_obs_vs_bcid["topoetcone20_final"]->Fill(bcid, lep_topoetcone20_final);
      histo_obs_vs_bcid["ptcone20_final_pt"]->Fill(bcid, lep_ptcone20_final / lep_pt);
      histo_obs_vs_bcid["etcone20_pt"]->Fill(bcid, lep_etcone20 / lep_pt);
      histo_obs_vs_bcid["etcone20_final_pt"]->Fill(bcid, lep_etcone20_final / lep_pt);
      histo_obs_vs_bcid["topoetcone20_final_pt"]->Fill(bcid, lep_topoetcone20_final / lep_pt);
      histo_obs_vs_nvx["ptcone20_final"]->Fill(nvx, lep_ptcone20_final);
      histo_obs_vs_nvx["etcone20"]->Fill(nvx, lep_etcone20);
      histo_obs_vs_nvx["etcone20_final"]->Fill(nvx, lep_etcone20_final);
      histo_obs_vs_nvx["topoetcone20_final"]->Fill(nvx, lep_topoetcone20_final);
      histo_obs_vs_nvx["ptcone20_final_pt"]->Fill(nvx, lep_ptcone20_final / lep_pt);
      histo_obs_vs_nvx["etcone20_pt"]->Fill(nvx, lep_etcone20 / lep_pt);
      histo_obs_vs_nvx["etcone20_final_pt"]->Fill(nvx, lep_etcone20_final / lep_pt);
      histo_obs_vs_nvx["topoetcone20_final_pt"]->Fill(nvx, lep_topoetcone20_final / lep_pt);
      histo_obs_vs_averageIntPerXing["ptcone20_final"]->Fill(averageIntPerXing, lep_ptcone20_final);
      histo_obs_vs_averageIntPerXing["etcone20"]->Fill(averageIntPerXing, lep_etcone20);
      histo_obs_vs_averageIntPerXing["etcone20_final"]->Fill(averageIntPerXing, lep_etcone20_final);
      histo_obs_vs_averageIntPerXing["topoetcone20_final"]->Fill(averageIntPerXing, lep_topoetcone20_final);
      histo_obs_vs_averageIntPerXing["ptcone20_final_pt"]->Fill(averageIntPerXing, lep_ptcone20_final / lep_pt);
      histo_obs_vs_averageIntPerXing["etcone20_pt"]->Fill(averageIntPerXing, lep_etcone20 / lep_pt);
      histo_obs_vs_averageIntPerXing["etcone20_final_pt"]->Fill(averageIntPerXing, lep_etcone20_final / lep_pt);
      histo_obs_vs_averageIntPerXing["topoetcone20_final_pt"]->Fill(averageIntPerXing, lep_topoetcone20_final / lep_pt);
      histo_obs_vs_actualIntPerXing["ptcone20_final"]->Fill(actualIntPerXing, lep_ptcone20_final);
      histo_obs_vs_actualIntPerXing["etcone20"]->Fill(actualIntPerXing, lep_etcone20);
      histo_obs_vs_actualIntPerXing["etcone20_final"]->Fill(actualIntPerXing, lep_etcone20_final);
      histo_obs_vs_actualIntPerXing["topoetcone20_final"]->Fill(actualIntPerXing, lep_topoetcone20_final);
      histo_obs_vs_actualIntPerXing["ptcone20_final_pt"]->Fill(actualIntPerXing, lep_ptcone20_final / lep_pt);
      histo_obs_vs_actualIntPerXing["etcone20_pt"]->Fill(actualIntPerXing, lep_etcone20 / lep_pt);
      histo_obs_vs_actualIntPerXing["etcone20_final_pt"]->Fill(actualIntPerXing, lep_etcone20_final / lep_pt);
      histo_obs_vs_actualIntPerXing["topoetcone20_final_pt"]->Fill(actualIntPerXing, lep_topoetcone20_final / lep_pt);
      histo_obs_vs_bcid["ptcone40_final"]->Fill(bcid, lep_ptcone40_final);
      histo_obs_vs_bcid["etcone40"]->Fill(bcid, lep_etcone40);
      histo_obs_vs_bcid["etcone40_final"]->Fill(bcid, lep_etcone40_final);
      histo_obs_vs_bcid["topoetcone40_final"]->Fill(bcid, lep_topoetcone40_final);
      histo_obs_vs_bcid["ptcone40_final_pt"]->Fill(bcid, lep_ptcone40_final / lep_pt);
      histo_obs_vs_bcid["etcone40_pt"]->Fill(bcid, lep_etcone40 / lep_pt);
      histo_obs_vs_bcid["etcone40_final_pt"]->Fill(bcid, lep_etcone40_final / lep_pt);
      histo_obs_vs_bcid["topoetcone40_final_pt"]->Fill(bcid, lep_topoetcone40_final / lep_pt);
      histo_obs_vs_nvx["ptcone40_final"]->Fill(nvx, lep_ptcone40_final);
      histo_obs_vs_nvx["etcone40"]->Fill(nvx, lep_etcone40);
      histo_obs_vs_nvx["etcone40_final"]->Fill(nvx, lep_etcone40_final);
      histo_obs_vs_nvx["topoetcone40_final"]->Fill(nvx, lep_topoetcone40_final);
      histo_obs_vs_nvx["ptcone40_final_pt"]->Fill(nvx, lep_ptcone40_final / lep_pt);
      histo_obs_vs_nvx["etcone40_pt"]->Fill(nvx, lep_etcone40 / lep_pt);
      histo_obs_vs_nvx["etcone40_final_pt"]->Fill(nvx, lep_etcone40_final / lep_pt);
      histo_obs_vs_nvx["topoetcone40_final_pt"]->Fill(nvx, lep_topoetcone40_final / lep_pt);
      histo_obs_vs_averageIntPerXing["ptcone40_final"]->Fill(averageIntPerXing, lep_ptcone40_final);
      histo_obs_vs_averageIntPerXing["etcone40"]->Fill(averageIntPerXing, lep_etcone40);
      histo_obs_vs_averageIntPerXing["etcone40_final"]->Fill(averageIntPerXing, lep_etcone40_final);
      histo_obs_vs_averageIntPerXing["topoetcone40_final"]->Fill(averageIntPerXing, lep_topoetcone40_final);
      histo_obs_vs_averageIntPerXing["ptcone40_final_pt"]->Fill(averageIntPerXing, lep_ptcone40_final / lep_pt);
      histo_obs_vs_averageIntPerXing["etcone40_pt"]->Fill(averageIntPerXing, lep_etcone40 / lep_pt);
      histo_obs_vs_averageIntPerXing["etcone40_final_pt"]->Fill(averageIntPerXing, lep_etcone40_final / lep_pt);
      histo_obs_vs_averageIntPerXing["topoetcone40_final_pt"]->Fill(averageIntPerXing, lep_topoetcone40_final / lep_pt);
      histo_obs_vs_actualIntPerXing["ptcone40_final"]->Fill(actualIntPerXing, lep_ptcone40_final);
      histo_obs_vs_actualIntPerXing["etcone40"]->Fill(actualIntPerXing, lep_etcone40);
      histo_obs_vs_actualIntPerXing["etcone40_final"]->Fill(actualIntPerXing, lep_etcone40_final);
      histo_obs_vs_actualIntPerXing["topoetcone40_final"]->Fill(actualIntPerXing, lep_topoetcone40_final);
      histo_obs_vs_actualIntPerXing["ptcone40_final_pt"]->Fill(actualIntPerXing, lep_ptcone40_final / lep_pt);
      histo_obs_vs_actualIntPerXing["etcone40_pt"]->Fill(actualIntPerXing, lep_etcone40 / lep_pt);
      histo_obs_vs_actualIntPerXing["etcone40_final_pt"]->Fill(actualIntPerXing, lep_etcone40_final / lep_pt);
      histo_obs_vs_actualIntPerXing["topoetcone40_final_pt"]->Fill(actualIntPerXing, lep_topoetcone40_final / lep_pt);


      // cut < 0.15
      eff_vs_nvx["ptcone20_lt_015"]->Fill(lep_ptcone20_final/lep_pt < 0.15, nvx);
      eff_vs_pt["ptcone20_lt_015"]->Fill(lep_ptcone20_final/lep_pt < 0.15, lep_pt/1000.);
      eff_vs_eta["ptcone20_lt_015"]->Fill(lep_ptcone20_final/lep_pt < 0.15, lep_eta);
      eff_vs_phi["ptcone20_lt_015"]->Fill(lep_ptcone20_final/lep_pt < 0.15, lep_phi);
      eff_vs_etaphi["ptcone20_lt_015"]->Fill(lep_ptcone20_final/lep_pt < 0.15, lep_phi, lep_eta);

      // cut < 0.2
      eff_vs_nvx["topoetcone20_lt_02"]->Fill(lep_topoetcone20_final/lep_pt < 0.2, nvx);
      eff_vs_pt["topoetcone20_lt_02"]->Fill(lep_topoetcone20_final/lep_pt < 0.2, lep_pt/1000.);
      eff_vs_eta["topoetcone20_lt_02"]->Fill(lep_topoetcone20_final/lep_pt < 0.2, lep_eta);
      eff_vs_phi["topoetcone20_lt_02"]->Fill(lep_topoetcone20_final/lep_pt < 0.2, lep_phi);
      eff_vs_etaphi["topoetcone20_lt_02"]->Fill(lep_topoetcone20_final/lep_pt < 0.2, lep_phi, lep_eta);

      if (lep_m < 100) { // electrons
        eff_vs_nvx["etcone20_lt_02"]->Fill(lep_etcone20_final/lep_pt < 0.2, nvx);
        eff_vs_pt["etcone20_lt_02"]->Fill(lep_etcone20_final/lep_pt < 0.2, lep_pt/1000.);
        eff_vs_eta["etcone20_lt_02"]->Fill(lep_etcone20_final/lep_pt < 0.2, lep_eta);
        eff_vs_phi["etcone20_lt_02"]->Fill(lep_etcone20_final/lep_pt < 0.2, lep_phi);
        eff_vs_etaphi["etcone20_lt_02"]->Fill(lep_etcone20_final/lep_pt < 0.2, lep_phi, lep_eta);
      } else { // muons
        eff_vs_nvx["etcone20_lt_02"]->Fill(lep_etcone20/lep_pt < 0.2, nvx);
        eff_vs_pt["etcone20_lt_02"]->Fill(lep_etcone20/lep_pt < 0.2, lep_pt/1000.);
        eff_vs_eta["etcone20_lt_02"]->Fill(lep_etcone20/lep_pt < 0.2, lep_eta);
        eff_vs_phi["etcone20_lt_02"]->Fill(lep_etcone20/lep_pt < 0.2, lep_phi);
        eff_vs_etaphi["etcone20_lt_02"]->Fill(lep_etcone20/lep_pt < 0.2, lep_phi, lep_eta);
      }

      // cut < 0.25
      eff_vs_nvx["topoetcone20_lt_025"]->Fill(lep_topoetcone20_final/lep_pt < 0.25, nvx);
      eff_vs_pt["topoetcone20_lt_025"]->Fill(lep_topoetcone20_final/lep_pt < 0.25, lep_pt/1000.);
      eff_vs_eta["topoetcone20_lt_025"]->Fill(lep_topoetcone20_final/lep_pt < 0.25, lep_eta);
      eff_vs_phi["topoetcone20_lt_025"]->Fill(lep_topoetcone20_final/lep_pt < 0.25, lep_phi);
      eff_vs_etaphi["topoetcone20_lt_025"]->Fill(lep_topoetcone20_final/lep_pt < 0.25, lep_phi, lep_eta);

      if (lep_m < 100) { // electrons
        eff_vs_nvx["etcone20_lt_025"]->Fill(lep_etcone20_final/lep_pt < 0.25, nvx);
        eff_vs_pt["etcone20_lt_025"]->Fill(lep_etcone20_final/lep_pt < 0.25, lep_pt/1000.);
        eff_vs_eta["etcone20_lt_025"]->Fill(lep_etcone20_final/lep_pt < 0.25, lep_eta);
        eff_vs_phi["etcone20_lt_025"]->Fill(lep_etcone20_final/lep_pt < 0.25, lep_phi);
        eff_vs_etaphi["etcone20_lt_025"]->Fill(lep_etcone20_final/lep_pt < 0.25, lep_phi, lep_eta);
      } else { // muons
        eff_vs_nvx["etcone20_lt_025"]->Fill(lep_etcone20/lep_pt < 0.25, nvx);
        eff_vs_pt["etcone20_lt_025"]->Fill(lep_etcone20/lep_pt < 0.25, lep_pt/1000.);
        eff_vs_eta["etcone20_lt_025"]->Fill(lep_etcone20/lep_pt < 0.25, lep_eta);
        eff_vs_phi["etcone20_lt_025"]->Fill(lep_etcone20/lep_pt < 0.25, lep_phi);
        eff_vs_etaphi["etcone20_lt_025"]->Fill(lep_etcone20/lep_pt < 0.25, lep_phi, lep_eta);
      }

      // cut < 0.3
      eff_vs_nvx["topoetcone20_lt_03"]->Fill(lep_topoetcone20_final/lep_pt < 0.3, nvx);
      eff_vs_pt["topoetcone20_lt_03"]->Fill(lep_topoetcone20_final/lep_pt < 0.3, lep_pt/1000.);
      eff_vs_eta["topoetcone20_lt_03"]->Fill(lep_topoetcone20_final/lep_pt < 0.3, lep_eta);
      eff_vs_phi["topoetcone20_lt_03"]->Fill(lep_topoetcone20_final/lep_pt < 0.3, lep_phi);
      eff_vs_etaphi["topoetcone20_lt_03"]->Fill(lep_topoetcone20_final/lep_pt < 0.3, lep_phi, lep_eta);

      if (lep_m < 100) { // electrons
        eff_vs_nvx["etcone20_lt_03"]->Fill(lep_etcone20_final/lep_pt < 0.3, nvx);
        eff_vs_pt["etcone20_lt_03"]->Fill(lep_etcone20_final/lep_pt < 0.3, lep_pt/1000.);
        eff_vs_eta["etcone20_lt_03"]->Fill(lep_etcone20_final/lep_pt < 0.3, lep_eta);
        eff_vs_phi["etcone20_lt_03"]->Fill(lep_etcone20_final/lep_pt < 0.3, lep_phi);
        eff_vs_etaphi["etcone20_lt_03"]->Fill(lep_etcone20_final/lep_pt < 0.3, lep_phi, lep_eta);
      } else { // muons
        eff_vs_nvx["etcone20_lt_03"]->Fill(lep_etcone20/lep_pt < 0.3, nvx);
        eff_vs_pt["etcone20_lt_03"]->Fill(lep_etcone20/lep_pt < 0.3, lep_pt/1000.);
        eff_vs_eta["etcone20_lt_03"]->Fill(lep_etcone20/lep_pt < 0.3, lep_eta);
        eff_vs_phi["etcone20_lt_03"]->Fill(lep_etcone20/lep_pt < 0.3, lep_phi);
        eff_vs_etaphi["etcone20_lt_03"]->Fill(lep_etcone20/lep_pt < 0.3, lep_phi, lep_eta);
      }


      for (std::vector<TString>::iterator obs = observables.begin(); obs != observables.end(); ++obs) {
         std::vector<Double_t> *proper_cut_values;

         if (obs->EndsWith("_pt")) {
            proper_cut_values = &rel_cut_values;
         } else {
            proper_cut_values = &abs_cut_values;
         }

         Double_t x(-99999);

         if (*obs == "ptcone20_final") x = lep_ptcone20_final;
         else if (*obs == "etcone20") x = lep_etcone20;
         else if (*obs == "etcone20_final") x = lep_etcone20_final;
         else if (*obs == "topoetcone20_final") x = lep_topoetcone20_final;
         else if (*obs == "ptcone20_final_pt") x = lep_ptcone20_final / lep_pt;
         else if (*obs == "etcone20_pt") x = lep_etcone20 / lep_pt;
         else if (*obs == "etcone20_final_pt") x = lep_etcone20_final / lep_pt;
         else if (*obs == "topoetcone20_final_pt") x = lep_topoetcone20_final / lep_pt;
         else if (*obs == "ptcone40_final") x = lep_ptcone40_final;
         else if (*obs == "etcone40") x = lep_etcone40;
         else if (*obs == "etcone40_final") x = lep_etcone40_final;
         else if (*obs == "topoetcone40_final") x = lep_topoetcone40_final;
         else if (*obs == "ptcone40_final_pt") x = lep_ptcone40_final / lep_pt;
         else if (*obs == "etcone40_pt") x = lep_etcone40 / lep_pt;
         else if (*obs == "etcone40_final_pt") x = lep_etcone40_final / lep_pt;
         else if (*obs == "topoetcone40_final_pt") x = lep_topoetcone40_final / lep_pt;

         for (std::vector<Double_t>::iterator val = proper_cut_values->begin(); val != proper_cut_values->end(); ++val) {
            den[*obs][*val] = den[*obs][*val] + poids;
            if (x < *val) {
               num[*obs][*val] = num[*obs][*val] + poids;
            }
         } // loop on cut values
      } // loop on observables
   } // lepton


   return kTRUE;
}

void IsolationEfficiency::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

   Info("SlaveTerminate", "Post-run coding...");

   outfile->cd();

   for (std::vector<TString>::iterator obs = observables.begin(); obs != observables.end(); ++obs) {
      profile_obs_vs_bcid[*obs] = histo_obs_vs_bcid[*obs]->ProfileX();
      profile_obs_vs_nvx[*obs] = histo_obs_vs_nvx[*obs]->ProfileX();
      profile_obs_vs_averageIntPerXing[*obs] = histo_obs_vs_averageIntPerXing[*obs]->ProfileX();
      profile_obs_vs_actualIntPerXing[*obs] = histo_obs_vs_actualIntPerXing[*obs]->ProfileX();
      /*
      fOutput->Add(profile_obs_vs_bcid[*obs]);
      fOutput->Add(profile_obs_vs_nvx[*obs]);
      fOutput->Add(profile_obs_vs_averageIntPerXing[*obs]);
      fOutput->Add(profile_obs_vs_actualIntPerXing[*obs]);
      */
   }

   for (std::vector<TString>::iterator obs = observables.begin(); obs != observables.end(); ++obs) {
      cut_scan[*obs] = new TGraph();
      cut_scan[*obs]->SetName(TString::Format("cut_scan_%s", obs->Data()));
      cut_scan[*obs]->SetTitle(TString::Format("cut_scan_%s", obs->Data()));
      cut_scan[*obs]->GetXaxis()->SetTitle(TString::Format("cut on lepton %s", obs->Data()));
      cut_scan[*obs]->GetYaxis()->SetTitle(TString::Format("efficiency"));
      /*
      fOutput->Add(cut_scan[*obs]);
      */
   }

   for (std::vector<TString>::iterator obs = observables.begin(); obs != observables.end(); ++obs) {
      Int_t i(0);

      std::vector<Double_t> *proper_cut_values;

      if (obs->EndsWith("_pt")) {
         proper_cut_values = &rel_cut_values;
      } else {
         proper_cut_values = &abs_cut_values;
      }

      for (std::vector<Double_t>::iterator val = proper_cut_values->begin(); val != proper_cut_values->end(); ++val) {
         cut_scan[*obs]->SetPoint(i, *val, num[*obs][*val] / den[*obs][*val]);
         i++;
      }
   }

   for (std::vector<TString>::iterator obs = observables.begin(); obs != observables.end(); ++obs) {
      cut_scan[*obs]->Write();
   }

   /*
   fOutput->Add(prooffile);
   */

   Info("SlaveTerminate", "Done.");
}

void IsolationEfficiency::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   Info("Terminate", "Saving to disk...");

   outfile->Write();
   outfile->Close();

   Info("Terminate", "Cleaning...");

  
   for (std::vector<TString>::iterator obs = observables.begin(); obs != observables.end(); ++obs) {
     den[*obs].clear();
     num[*obs].clear();
     if (cut_scan[*obs]) delete cut_scan[*obs];
   }

   den.clear();
   num.clear();
   cut_scan.clear();
   observables.clear();
   rel_cut_values.clear();
   abs_cut_values.clear();

   Info("Terminate", "Done.");
}
