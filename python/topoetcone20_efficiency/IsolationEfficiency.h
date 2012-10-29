//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May  2 13:02:05 2012 by ROOT version 5.32/00-rc2
// from TTree leptontree/lepton ntuple
// found on file: user.vippolit.004324._00001.output_isol.root
//////////////////////////////////////////////////////////

#ifndef IsolationEfficiency_h
#define IsolationEfficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TProofOutputFile.h>
#include <TSelector.h>

#include <TProfile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <map>

#include <TEfficiency.h>

//#include "HiggsZZ4lUtils/ElectronMCClassification.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class IsolationEfficiency : public TSelector {
private:
   TProofOutputFile *prooffile;
   TFile *outfile;
   std::vector<Double_t> rel_cut_values;
   std::vector<Double_t> abs_cut_values;
   TString m_tag;
   Bool_t m_isElectron;
   Bool_t m_useElectronTruth;
   Bool_t m_useConversions;
   Bool_t m_useFakes;
   Bool_t m_useSignalElectrons;
   std::pair<Double_t, Double_t> rel_minmax;
   std::pair<Double_t, Double_t> abs_minmax;
   std::vector<TString> observables;
   std::map<TString, TH1F*> histo_obs;
   std::map<TString, TH2F*> histo_obs_vs_bcid;
   std::map<TString, TH2F*> histo_obs_vs_nvx;
   std::map<TString, TH2F*> histo_obs_vs_averageIntPerXing;
   std::map<TString, TH2F*> histo_obs_vs_actualIntPerXing;
   std::map<TString, TProfile*> profile_obs_vs_bcid;
   std::map<TString, TProfile*> profile_obs_vs_nvx;
   std::map<TString, TProfile*> profile_obs_vs_averageIntPerXing;
   std::map<TString, TProfile*> profile_obs_vs_actualIntPerXing;
   std::map<TString, std::map<Double_t, Double_t> > den;
   std::map<TString, std::map<Double_t, Double_t> > num;
   std::map<TString, TGraph *> cut_scan;
   std::map<TString, TEfficiency *>eff_vs_nvx;
   std::map<TString, TEfficiency *>eff_vs_pt;
   std::map<TString, TEfficiency *>eff_vs_eta;
   std::map<TString, TEfficiency *>eff_vs_phi;
   std::map<TString, TEfficiency *>eff_vs_etaphi;

public:
   void SetTag(TString val) {
      m_tag = val;
   }
   void SetFlavor(TString val) {
      if (val == "ELECTRON") {
        m_isElectron = kTRUE;
      } else if (val == "MUON") {
        m_isElectron = kFALSE;
      } else {
        Error("", "wrong lepton flavor %s, assuming ELECTRON", val.Data());
	m_isElectron = kFALSE;
      }
   }
   void SetTruth(Bool_t val) {
     m_useElectronTruth = kTRUE;
   }
   void SetClassification(Bool_t signal, Bool_t conversions, Bool_t fakes) {
     m_useConversions = conversions;
     m_useFakes = fakes;
     m_useSignalElectrons = signal;
   }
   const char *GetTag() {
      return m_tag.Data();
   }

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          event;
   Int_t           last;
   Int_t           nvx;
   Int_t           bcid;
   Float_t         actualIntPerXing;
   Float_t         averageIntPerXing;
   Float_t         xsec_weight;
   Float_t         pu_weight;
   Float_t         ggF_weight;
   Float_t         lep_weight;
   Int_t           lep_author;
   Int_t           lep_isSA;
   Float_t         lep_pt;
   Float_t         lep_eta;
   Float_t         lep_phi;
   Float_t         lep_m;
   Float_t         lep_d0;
   Float_t         lep_z0;
   Float_t         lep_etcone20;
   Float_t         lep_etcone20_final;
   Float_t         lep_topoetcone20;
   Float_t         lep_topoetcone20_final;
   Float_t         lep_ptcone20;
   Float_t         lep_ptcone20_final;
   Float_t         lep_etcone40;
   Float_t         lep_etcone40_final;
   Float_t         lep_topoetcone40;
   Float_t         lep_topoetcone40_final;
   Float_t         lep_ptcone40;
   Float_t         lep_ptcone40_final;
   Int_t           lep_type;
   Int_t           lep_origin;
   Int_t           lep_originbkg;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_last;   //!
   TBranch        *b_nvx;   //!
   TBranch        *b_bcid;   //!
   TBranch        *b_actualIntPerXing;   //!
   TBranch        *b_averageIntPerXing;   //!
   TBranch        *b_xsec_weight;   //!
   TBranch        *b_pu_weight;   //!
   TBranch        *b_ggF_weight;   //!
   TBranch        *b_lep_weight;   //!
   TBranch        *b_lep_author;   //!
   TBranch        *b_lep_isSA;   //!
   TBranch        *b_lep_pt;   //!
   TBranch        *b_lep_eta;   //!
   TBranch        *b_lep_phi;   //!
   TBranch        *b_lep_m;   //!
   TBranch        *b_lep_d0;   //!
   TBranch        *b_lep_z0;   //!
   TBranch        *b_lep_etcone20;   //!
   TBranch        *b_lep_etcone20_final;   //!
   TBranch        *b_lep_topoetcone20;   //!
   TBranch        *b_lep_topoetcone20_final;   //!
   TBranch        *b_lep_ptcone20;   //!
   TBranch        *b_lep_ptcone20_final;   //!
   TBranch        *b_lep_etcone40;   //!
   TBranch        *b_lep_etcone40_final;   //!
   TBranch        *b_lep_topoetcone40;   //!
   TBranch        *b_lep_topoetcone40_final;   //!
   TBranch        *b_lep_ptcone40;   //!
   TBranch        *b_lep_ptcone40_final;   //!
   TBranch        *b_lep_type;   //!
   TBranch        *b_lep_origin;   //!
   TBranch        *b_lep_originbkg;   //!

   IsolationEfficiency(TTree * /*tree*/ = 0) : fChain(0), m_useConversions(kFALSE), m_useFakes(kFALSE), m_useSignalElectrons(kFALSE) { }
   virtual ~IsolationEfficiency() { }
   virtual Int_t   Version() const {
      return 2;
   }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) {
      return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
   }
   virtual void    SetOption(const char *option) {
      fOption = option;
   }
   virtual void    SetObject(TObject *obj) {
      fObject = obj;
   }
   virtual void    SetInputList(TList *input) {
      fInput = input;
   }
   virtual TList  *GetOutputList() const {
      return fOutput;
   }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(IsolationEfficiency, 0);
};

#endif

#ifdef IsolationEfficiency_cxx
void IsolationEfficiency::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("last", &last, &b_last);
   fChain->SetBranchAddress("nvx", &nvx, &b_nvx);
   fChain->SetBranchAddress("bcid", &bcid, &b_bcid);
   fChain->SetBranchAddress("actualIntPerXing", &actualIntPerXing, &b_actualIntPerXing);
   fChain->SetBranchAddress("averageIntPerXing", &averageIntPerXing, &b_averageIntPerXing);
   fChain->SetBranchAddress("xsec_weight", &xsec_weight, &b_xsec_weight);
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fChain->SetBranchAddress("ggF_weight", &ggF_weight, &b_ggF_weight);
   fChain->SetBranchAddress("lep_weight", &lep_weight, &b_lep_weight);
   fChain->SetBranchAddress("lep_author", &lep_author, &b_lep_author);
   fChain->SetBranchAddress("lep_isSA", &lep_isSA, &b_lep_isSA);
   fChain->SetBranchAddress("lep_pt", &lep_pt, &b_lep_pt);
   fChain->SetBranchAddress("lep_eta", &lep_eta, &b_lep_eta);
   fChain->SetBranchAddress("lep_phi", &lep_phi, &b_lep_phi);
   fChain->SetBranchAddress("lep_m", &lep_m, &b_lep_m);
   fChain->SetBranchAddress("lep_d0", &lep_d0, &b_lep_d0);
   fChain->SetBranchAddress("lep_z0", &lep_z0, &b_lep_z0);
   fChain->SetBranchAddress("lep_etcone20", &lep_etcone20, &b_lep_etcone20);
   fChain->SetBranchAddress("lep_etcone20_final", &lep_etcone20_final, &b_lep_etcone20_final);
   fChain->SetBranchAddress("lep_topoetcone20", &lep_topoetcone20, &b_lep_topoetcone20);
   fChain->SetBranchAddress("lep_topoetcone20_final", &lep_topoetcone20_final, &b_lep_topoetcone20_final);
   fChain->SetBranchAddress("lep_ptcone20", &lep_ptcone20, &b_lep_ptcone20);
   fChain->SetBranchAddress("lep_ptcone20_final", &lep_ptcone20_final, &b_lep_ptcone20_final);
   fChain->SetBranchAddress("lep_etcone40", &lep_etcone40, &b_lep_etcone40);
   fChain->SetBranchAddress("lep_etcone40_final", &lep_etcone40_final, &b_lep_etcone40_final);
   fChain->SetBranchAddress("lep_topoetcone40", &lep_topoetcone40, &b_lep_topoetcone40);
   fChain->SetBranchAddress("lep_topoetcone40_final", &lep_topoetcone40_final, &b_lep_topoetcone40_final);
   fChain->SetBranchAddress("lep_ptcone40", &lep_ptcone40, &b_lep_ptcone40);
   fChain->SetBranchAddress("lep_ptcone40_final", &lep_ptcone40_final, &b_lep_ptcone40_final);
   fChain->SetBranchAddress("lep_type", &lep_type, &b_lep_type);
   fChain->SetBranchAddress("lep_origin", &lep_origin, &b_lep_origin);
   fChain->SetBranchAddress("lep_originbkg", &lep_originbkg, &b_lep_originbkg);
}

Bool_t IsolationEfficiency::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef IsolationEfficiency_cxx
