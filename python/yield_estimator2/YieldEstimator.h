//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jun 25 13:46:51 2012 by ROOT version 5.32/00
// from TTree candidates/SFOS Higgs candidates
// found on file: /afs/cern.ch/work/v/vippolit/ntuple_2012_v20/user.vippolit.H4l000020.146832.AlpgenJimmy_Auto_AUET2CTEQ6L1_ZeeNp2Excl_Mll10to60.p1044.2012_campaign.01.120625051350/user.vippolit.009834._00001.output_test.root
//////////////////////////////////////////////////////////

#ifndef YieldEstimator_h
#define YieldEstimator_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

#define N_CHANNELS 4
#define BIN_FOR_GENERATED 5

#include "../../../HiggsZZ4lUtils/HiggsZZ4lUtils/HiggsCrossSection.h"
#include "../../../HiggsZZ4lUtils/HiggsZZ4lUtils/BkgCrossSection.h"

using namespace std;

namespace MC {
   typedef enum {
      none = 0,
      // signal
      ggF,
      VBF,
      WH,
      ZH,
      // reducible
      Z,
      Zbb,
      tt,
      // irreducible
      ZZ,
      gg2ZZ,
   } ProductionProcess;
};

class YieldSample {
public:
   TString name;
   Bool_t isSignal;
   Float_t mass;
   MC::ProductionProcess type;
   Float_t xsec;
   Float_t generated;
   std::vector<TH1F*> h_H_m;
   std::vector<TH1F*> h_H_m_constrained;
   std::vector<TH1F*> h_Z1_m;
   std::vector<TH1F*> h_Z2_m;

   YieldSample() : name(""), isSignal(kFALSE), mass(-1), type(MC::none), xsec(-1), generated(-1) { }
   YieldSample(TString _name, Bool_t _isSignal = kFALSE, Float_t _mass = -1., MC::ProductionProcess _type = MC::none) : name(_name), isSignal(_isSignal), mass(_mass), type(_type), xsec(-1), generated(-1) {
      std::cout << "INFO : creating sample " << name << std::endl;
   }
   ~YieldSample() {}

   void print() const {
      std::cout << "INFO :       name = " << name << std::endl
                << "         isSignal = " << isSignal << std::endl
                << "             mass = " << mass << std::endl
                << "             type = " << type << std::endl
                << "             xsec = " << xsec << std::endl
                << "        generated = " << generated << std::endl;
   }
};

class FinalHistograms {
public:
   Double_t bin_width;

   // signal
   std::vector<TH1F*> higgs;
   std::vector<TH1F*> higgsVBF;
   std::vector<TH1F*> higgsWH;
   std::vector<TH1F*> higgsZH;
   // reducible
   std::vector<TH1F*> histoZ;
   std::vector<TH1F*> histoZbb;
   std::vector<TH1F*> histott;
   // irreducible
   std::vector<TH1F*> histoZZ;
   std::vector<TH1F*> histogg2ZZ;

   Bool_t rebin() {
      Int_t rebin_factor = (Int_t)(bin_width / higgs[0]->GetBinWidth(1));
      if (((Int_t)(bin_width * 10)) % ((Int_t)(higgs[0]->GetBinWidth(1) * 10)) != 0)
         std::cout << "WARNING : bin width is " << higgs[0]->GetBinWidth(1) << " but rebin factor to " << bin_width << " is not correct!" << std::endl;

      for (UInt_t i = 0; i < N_CHANNELS; i++) {
         higgs[i]->Rebin(rebin_factor);
         higgsVBF[i]->Rebin(rebin_factor);
         higgsWH[i]->Rebin(rebin_factor);
         higgsZH[i]->Rebin(rebin_factor);
         histoZ[i]->Rebin(rebin_factor);
         histoZbb[i]->Rebin(rebin_factor);
         histott[i]->Rebin(rebin_factor);
         histoZZ[i]->Rebin(rebin_factor);
         histogg2ZZ[i]->Rebin(rebin_factor);
      }

      return kTRUE;
   }
   void setDirectory(TDirectory *dir) {
      for (UInt_t i = 0; i < N_CHANNELS; i++) {
         higgs[i]->SetDirectory(dir);
         higgsVBF[i]->SetDirectory(dir);
         higgsWH[i]->SetDirectory(dir);
         higgsZH[i]->SetDirectory(dir);
         histoZ[i]->SetDirectory(dir);
         histoZbb[i]->SetDirectory(dir);
         histott[i]->SetDirectory(dir);
         histoZZ[i]->SetDirectory(dir);
         histogg2ZZ[i]->SetDirectory(dir);
      }
   }
   FinalHistograms(Double_t _bin_width = 0.5) {
      bin_width = _bin_width; // bin width is applied only when rebin() is called

      for (UInt_t i = 0; i < N_CHANNELS; i++) {
         higgs.push_back(new TH1F(TString::Format("higgs_%d", i), "", 2000, 0, 1000));
         higgsVBF.push_back(new TH1F(TString::Format("higgsVBF_%d", i), "", 2000, 0, 1000));
         higgsWH.push_back(new TH1F(TString::Format("higgsWH_%d", i), "", 2000, 0, 1000));
         higgsZH.push_back(new TH1F(TString::Format("higgsZH_%d", i), "", 2000, 0, 1000));
         histoZ.push_back(new TH1F(TString::Format("histoZ_%d", i), "", 2000, 0, 1000));
         histoZbb.push_back(new TH1F(TString::Format("histoZbb_%d", i), "", 2000, 0, 1000));
         histott.push_back(new TH1F(TString::Format("histott_%d", i), "", 2000, 0, 1000));
         histoZZ.push_back(new TH1F(TString::Format("histoZZ_%d", i), "", 2000, 0, 1000));
         histogg2ZZ.push_back(new TH1F(TString::Format("histogg2ZZ_%d", i), "", 2000, 0, 1000));
      }
   }
   FinalHistograms(FinalHistograms *other, Double_t _bin_width = 0.5) {
      for (UInt_t i = 0; i < N_CHANNELS; i++) {
         bin_width = _bin_width;
         higgs.push_back(dynamic_cast<TH1F*>(other->higgs[i]->Clone()));
         higgsVBF.push_back(dynamic_cast<TH1F*>(other->higgsVBF[i]->Clone()));
         higgsWH.push_back(dynamic_cast<TH1F*>(other->higgsWH[i]->Clone()));
         higgsZH.push_back(dynamic_cast<TH1F*>(other->higgsZH[i]->Clone()));
         histoZ.push_back(dynamic_cast<TH1F*>(other->histoZ[i]->Clone()));
         histoZbb.push_back(dynamic_cast<TH1F*>(other->histoZbb[i]->Clone()));
         histott.push_back(dynamic_cast<TH1F*>(other->histott[i]->Clone()));
         histoZZ.push_back(dynamic_cast<TH1F*>(other->histoZZ[i]->Clone()));
         histogg2ZZ.push_back(dynamic_cast<TH1F*>(other->histogg2ZZ[i]->Clone()));
      }
   }
   ~FinalHistograms() {
      // does not maintain ownership
   }
};



// Fixed size dimensions of array or collections stored in the TTree if any.

class YieldEstimator : public TSelector {
public:
   HiggsCrossSection higgs_xsector;
   CrossSections::LHCEnergy m_COM_energy;

   Double_t lumi[N_CHANNELS];
   Double_t generated[N_CHANNELS];

   TFile *m_output_with_constraint;
   TFile *m_output_without_constraint;
   TFile *m_output_raw;
   std::map<UInt_t, YieldSample> m_sample;
   std::map<TString, FinalHistograms> m_histos;

   Bool_t fillMCMap(TString release);

   UInt_t m_current_run;

   Bool_t m_interesting;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           last;
   Int_t           closesttozz;
   Int_t           selected;
   Int_t           type;
   UInt_t          run;
   UInt_t          event;
   UInt_t          lbn;
   Int_t           n_vx;
   Int_t           hfor;
   Int_t           top_weight;
   Int_t           powhegbug_weight;
   Int_t           actualIntPerXing;
   Int_t           averageIntPerXing;
   Float_t         xsec_weight;
   Long64_t        processed_entries;
   Long64_t        generated_entries;
   Float_t         pu_weight;
   Float_t         vxz_weight;
   Float_t         OQ_weight;
   Float_t         ggF_weight;
   Float_t         trigSF_weight;
   Float_t         Z1_lepplus_weight;
   Float_t         Z1_lepminus_weight;
   Float_t         Z2_lepplus_weight;
   Float_t         Z2_lepminus_weight;
   Float_t         H_m_truth;
   Float_t         ZZ_m_truth;
   Float_t         Z1_m_truth;
   Float_t         Z2_m_truth;
   Float_t         H_m_constrained;
   Float_t         Z1_m_constrained;
   Float_t         Z2_m_constrained;
   Float_t         Z1_chiSq;
   Float_t         Z2_chiSq;
   Float_t         H_m;
   Float_t         H_pt;
   Float_t         H_eta;
   Float_t         H_phi;
   Float_t         Z1_m;
   Float_t         Z1_pt;
   Float_t         Z1_eta;
   Float_t         Z1_phi;
   Float_t         Z2_m;
   Float_t         Z2_pt;
   Float_t         Z2_eta;
   Float_t         Z2_phi;
   Float_t         Z1_lepplus_m;
   Float_t         Z1_lepplus_pt;
   Float_t         Z1_lepplus_eta;
   Float_t         Z1_lepplus_phi;
   Float_t         Z1_lepplus_etcone20;
   Float_t         Z1_lepplus_etcone20_corr;
   Float_t         Z1_lepplus_etcone20_final;
   Float_t         Z1_lepplus_ptcone20;
   Float_t         Z1_lepplus_ptcone20_final;
   Float_t         Z1_lepplus_d0;
   Float_t         Z1_lepplus_z0;
   Float_t         Z1_lepplus_d0_sig;
   Float_t         Z1_lepplus_GSF_dp;
   Float_t         Z1_lepplus_GSF_p;
   Float_t         Z1_lepplus_charge;
   Float_t         Z1_lepplus_ptCB_nosmearnoscale;
   Float_t         Z1_lepplus_ptME_nosmearnoscale;
   Float_t         Z1_lepplus_ptID_nosmearnoscale;
   Float_t         Z1_lepplus_cl_E_calibsmearnoscale;
   Float_t         Z1_lepplus_cl_eta;
   Float_t         Z1_lepplus_cl_phi;
   Int_t           Z1_lepplus_isSA;
   Int_t           Z1_lepplus_author;
   Int_t           Z1_lepplus_type;
   Int_t           Z1_lepplus_typebkg;
   Int_t           Z1_lepplus_origin;
   Int_t           Z1_lepplus_originbkg;
   Float_t         Z1_lepminus_m;
   Float_t         Z1_lepminus_pt;
   Float_t         Z1_lepminus_eta;
   Float_t         Z1_lepminus_phi;
   Float_t         Z1_lepminus_etcone20;
   Float_t         Z1_lepminus_etcone20_corr;
   Float_t         Z1_lepminus_etcone20_final;
   Float_t         Z1_lepminus_ptcone20;
   Float_t         Z1_lepminus_ptcone20_final;
   Float_t         Z1_lepminus_d0;
   Float_t         Z1_lepminus_z0;
   Float_t         Z1_lepminus_d0_sig;
   Float_t         Z1_lepminus_GSF_dp;
   Float_t         Z1_lepminus_GSF_p;
   Float_t         Z1_lepminus_charge;
   Float_t         Z1_lepminus_ptCB_nosmearnoscale;
   Float_t         Z1_lepminus_ptME_nosmearnoscale;
   Float_t         Z1_lepminus_ptID_nosmearnoscale;
   Float_t         Z1_lepminus_cl_E_calibsmearnoscale;
   Float_t         Z1_lepminus_cl_eta;
   Float_t         Z1_lepminus_cl_phi;
   Int_t           Z1_lepminus_isSA;
   Int_t           Z1_lepminus_author;
   Int_t           Z1_lepminus_type;
   Int_t           Z1_lepminus_typebkg;
   Int_t           Z1_lepminus_origin;
   Int_t           Z1_lepminus_originbkg;
   Float_t         Z2_lepplus_m;
   Float_t         Z2_lepplus_pt;
   Float_t         Z2_lepplus_eta;
   Float_t         Z2_lepplus_phi;
   Float_t         Z2_lepplus_etcone20;
   Float_t         Z2_lepplus_etcone20_corr;
   Float_t         Z2_lepplus_etcone20_final;
   Float_t         Z2_lepplus_ptcone20;
   Float_t         Z2_lepplus_ptcone20_final;
   Float_t         Z2_lepplus_d0;
   Float_t         Z2_lepplus_z0;
   Float_t         Z2_lepplus_d0_sig;
   Float_t         Z2_lepplus_GSF_dp;
   Float_t         Z2_lepplus_GSF_p;
   Float_t         Z2_lepplus_charge;
   Float_t         Z2_lepplus_ptCB_nosmearnoscale;
   Float_t         Z2_lepplus_ptME_nosmearnoscale;
   Float_t         Z2_lepplus_ptID_nosmearnoscale;
   Float_t         Z2_lepplus_cl_E_calibsmearnoscale;
   Float_t         Z2_lepplus_cl_eta;
   Float_t         Z2_lepplus_cl_phi;
   Int_t           Z2_lepplus_isSA;
   Int_t           Z2_lepplus_author;
   Int_t           Z2_lepplus_type;
   Int_t           Z2_lepplus_typebkg;
   Int_t           Z2_lepplus_origin;
   Int_t           Z2_lepplus_originbkg;
   Float_t         Z2_lepminus_m;
   Float_t         Z2_lepminus_pt;
   Float_t         Z2_lepminus_eta;
   Float_t         Z2_lepminus_phi;
   Float_t         Z2_lepminus_etcone20;
   Float_t         Z2_lepminus_etcone20_corr;
   Float_t         Z2_lepminus_etcone20_final;
   Float_t         Z2_lepminus_ptcone20;
   Float_t         Z2_lepminus_ptcone20_final;
   Float_t         Z2_lepminus_d0;
   Float_t         Z2_lepminus_z0;
   Float_t         Z2_lepminus_d0_sig;
   Float_t         Z2_lepminus_GSF_dp;
   Float_t         Z2_lepminus_GSF_p;
   Float_t         Z2_lepminus_charge;
   Float_t         Z2_lepminus_ptCB_nosmearnoscale;
   Float_t         Z2_lepminus_ptME_nosmearnoscale;
   Float_t         Z2_lepminus_ptID_nosmearnoscale;
   Float_t         Z2_lepminus_cl_E_calibsmearnoscale;
   Float_t         Z2_lepminus_cl_eta;
   Float_t         Z2_lepminus_cl_phi;
   Int_t           Z2_lepminus_isSA;
   Int_t           Z2_lepminus_author;
   Int_t           Z2_lepminus_type;
   Int_t           Z2_lepminus_typebkg;
   Int_t           Z2_lepminus_origin;
   Int_t           Z2_lepminus_originbkg;
   vector<string>  *filename;

   // List of branches
   TBranch        *b_last;   //!
   TBranch        *b_closesttozz;   //!
   TBranch        *b_selected;   //!
   TBranch        *b_type;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_n_vx;   //!
   TBranch        *b_hfor;   //!
   TBranch        *b_top_weight;   //!
   TBranch        *b_powhegbug_weight;   //!
   TBranch        *b_actualIntPerXing;   //!
   TBranch        *b_averageIntPerXing;   //!
   TBranch        *b_xsec_weight;   //!
   TBranch        *b_processed_entries;   //!
   TBranch        *b_generated_entries;   //!
   TBranch        *b_pu_weight;   //!
   TBranch        *b_vxz_weight;   //!
   TBranch        *b_OQ_weight;   //!
   TBranch        *b_ggF_weight;   //!
   TBranch        *b_trigSF_weight;   //!
   TBranch        *b_Z1_lepplus_weight;   //!
   TBranch        *b_Z1_lepminus_weight;   //!
   TBranch        *b_Z2_lepplus_weight;   //!
   TBranch        *b_Z2_lepminus_weight;   //!
   TBranch        *b_H_m_truth;   //!
   TBranch        *b_ZZ_m_truth;   //!
   TBranch        *b_Z1_m_truth;   //!
   TBranch        *b_Z2_m_truth;   //!
   TBranch        *b_H_m_constrained;   //!
   TBranch        *b_Z1_m_constrained;   //!
   TBranch        *b_Z2_m_constrained;   //!
   TBranch        *b_Z1_chiSq;   //!
   TBranch        *b_Z2_chiSq;   //!
   TBranch        *b_H_m;   //!
   TBranch        *b_H_pt;   //!
   TBranch        *b_H_eta;   //!
   TBranch        *b_H_phi;   //!
   TBranch        *b_Z1_m;   //!
   TBranch        *b_Z1_pt;   //!
   TBranch        *b_Z1_eta;   //!
   TBranch        *b_Z1_phi;   //!
   TBranch        *b_Z2_m;   //!
   TBranch        *b_Z2_pt;   //!
   TBranch        *b_Z2_eta;   //!
   TBranch        *b_Z2_phi;   //!
   TBranch        *b_Z1_lepplus_m;   //!
   TBranch        *b_Z1_lepplus_pt;   //!
   TBranch        *b_Z1_lepplus_eta;   //!
   TBranch        *b_Z1_lepplus_phi;   //!
   TBranch        *b_Z1_lepplus_etcone20;   //!
   TBranch        *b_Z1_lepplus_etcone20_corr;   //!
   TBranch        *b_Z1_lepplus_etcone20_final;   //!
   TBranch        *b_Z1_lepplus_ptcone20;   //!
   TBranch        *b_Z1_lepplus_ptcone20_final;   //!
   TBranch        *b_Z1_lepplus_d0;   //!
   TBranch        *b_Z1_lepplus_z0;   //!
   TBranch        *b_Z1_lepplus_d0_sig;   //!
   TBranch        *b_Z1_lepplus_GSF_dp;   //!
   TBranch        *b_Z1_lepplus_GSF_p;   //!
   TBranch        *b_Z1_lepplus_charge;   //!
   TBranch        *b_Z1_lepplus_ptCB_nosmearnoscale;   //!
   TBranch        *b_Z1_lepplus_ptME_nosmearnoscale;   //!
   TBranch        *b_Z1_lepplus_ptID_nosmearnoscale;   //!
   TBranch        *b_Z1_lepplus_cl_E_calibsmearnoscale;   //!
   TBranch        *b_Z1_lepplus_cl_eta;   //!
   TBranch        *b_Z1_lepplus_cl_phi;   //!
   TBranch        *b_Z1_lepplus_isSA;   //!
   TBranch        *b_Z1_lepplus_author;   //!
   TBranch        *b_Z1_lepplus_type;   //!
   TBranch        *b_Z1_lepplus_typebkg;   //!
   TBranch        *b_Z1_lepplus_origin;   //!
   TBranch        *b_Z1_lepplus_originbkg;   //!
   TBranch        *b_Z1_lepminus_m;   //!
   TBranch        *b_Z1_lepminus_pt;   //!
   TBranch        *b_Z1_lepminus_eta;   //!
   TBranch        *b_Z1_lepminus_phi;   //!
   TBranch        *b_Z1_lepminus_etcone20;   //!
   TBranch        *b_Z1_lepminus_etcone20_corr;   //!
   TBranch        *b_Z1_lepminus_etcone20_final;   //!
   TBranch        *b_Z1_lepminus_ptcone20;   //!
   TBranch        *b_Z1_lepminus_ptcone20_final;   //!
   TBranch        *b_Z1_lepminus_d0;   //!
   TBranch        *b_Z1_lepminus_z0;   //!
   TBranch        *b_Z1_lepminus_d0_sig;   //!
   TBranch        *b_Z1_lepminus_GSF_dp;   //!
   TBranch        *b_Z1_lepminus_GSF_p;   //!
   TBranch        *b_Z1_lepminus_charge;   //!
   TBranch        *b_Z1_lepminus_ptCB_nosmearnoscale;   //!
   TBranch        *b_Z1_lepminus_ptME_nosmearnoscale;   //!
   TBranch        *b_Z1_lepminus_ptID_nosmearnoscale;   //!
   TBranch        *b_Z1_lepminus_cl_E_calibsmearnoscale;   //!
   TBranch        *b_Z1_lepminus_cl_eta;   //!
   TBranch        *b_Z1_lepminus_cl_phi;   //!
   TBranch        *b_Z1_lepminus_isSA;   //!
   TBranch        *b_Z1_lepminus_author;   //!
   TBranch        *b_Z1_lepminus_type;   //!
   TBranch        *b_Z1_lepminus_typebkg;   //!
   TBranch        *b_Z1_lepminus_origin;   //!
   TBranch        *b_Z1_lepminus_originbkg;   //!
   TBranch        *b_Z2_lepplus_m;   //!
   TBranch        *b_Z2_lepplus_pt;   //!
   TBranch        *b_Z2_lepplus_eta;   //!
   TBranch        *b_Z2_lepplus_phi;   //!
   TBranch        *b_Z2_lepplus_etcone20;   //!
   TBranch        *b_Z2_lepplus_etcone20_corr;   //!
   TBranch        *b_Z2_lepplus_etcone20_final;   //!
   TBranch        *b_Z2_lepplus_ptcone20;   //!
   TBranch        *b_Z2_lepplus_ptcone20_final;   //!
   TBranch        *b_Z2_lepplus_d0;   //!
   TBranch        *b_Z2_lepplus_z0;   //!
   TBranch        *b_Z2_lepplus_d0_sig;   //!
   TBranch        *b_Z2_lepplus_GSF_dp;   //!
   TBranch        *b_Z2_lepplus_GSF_p;   //!
   TBranch        *b_Z2_lepplus_charge;   //!
   TBranch        *b_Z2_lepplus_ptCB_nosmearnoscale;   //!
   TBranch        *b_Z2_lepplus_ptME_nosmearnoscale;   //!
   TBranch        *b_Z2_lepplus_ptID_nosmearnoscale;   //!
   TBranch        *b_Z2_lepplus_cl_E_calibsmearnoscale;   //!
   TBranch        *b_Z2_lepplus_cl_eta;   //!
   TBranch        *b_Z2_lepplus_cl_phi;   //!
   TBranch        *b_Z2_lepplus_isSA;   //!
   TBranch        *b_Z2_lepplus_author;   //!
   TBranch        *b_Z2_lepplus_type;   //!
   TBranch        *b_Z2_lepplus_typebkg;   //!
   TBranch        *b_Z2_lepplus_origin;   //!
   TBranch        *b_Z2_lepplus_originbkg;   //!
   TBranch        *b_Z2_lepminus_m;   //!
   TBranch        *b_Z2_lepminus_pt;   //!
   TBranch        *b_Z2_lepminus_eta;   //!
   TBranch        *b_Z2_lepminus_phi;   //!
   TBranch        *b_Z2_lepminus_etcone20;   //!
   TBranch        *b_Z2_lepminus_etcone20_corr;   //!
   TBranch        *b_Z2_lepminus_etcone20_final;   //!
   TBranch        *b_Z2_lepminus_ptcone20;   //!
   TBranch        *b_Z2_lepminus_ptcone20_final;   //!
   TBranch        *b_Z2_lepminus_d0;   //!
   TBranch        *b_Z2_lepminus_z0;   //!
   TBranch        *b_Z2_lepminus_d0_sig;   //!
   TBranch        *b_Z2_lepminus_GSF_dp;   //!
   TBranch        *b_Z2_lepminus_GSF_p;   //!
   TBranch        *b_Z2_lepminus_charge;   //!
   TBranch        *b_Z2_lepminus_ptCB_nosmearnoscale;   //!
   TBranch        *b_Z2_lepminus_ptME_nosmearnoscale;   //!
   TBranch        *b_Z2_lepminus_ptID_nosmearnoscale;   //!
   TBranch        *b_Z2_lepminus_cl_E_calibsmearnoscale;   //!
   TBranch        *b_Z2_lepminus_cl_eta;   //!
   TBranch        *b_Z2_lepminus_cl_phi;   //!
   TBranch        *b_Z2_lepminus_isSA;   //!
   TBranch        *b_Z2_lepminus_author;   //!
   TBranch        *b_Z2_lepminus_type;   //!
   TBranch        *b_Z2_lepminus_typebkg;   //!
   TBranch        *b_Z2_lepminus_origin;   //!
   TBranch        *b_Z2_lepminus_originbkg;   //!
   TBranch        *b_filename;   //!

   YieldEstimator(TTree * /*tree*/ = 0) : m_current_run(0), fChain(0) { }
   virtual ~YieldEstimator() { }
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

   ClassDef(YieldEstimator, 0);
};

#endif

#ifdef YieldEstimator_cxx
void YieldEstimator::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   filename = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("last", &last, &b_last);
   fChain->SetBranchAddress("closesttozz", &closesttozz, &b_closesttozz);
   fChain->SetBranchAddress("selected", &selected, &b_selected);
   fChain->SetBranchAddress("type", &type, &b_type);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("n_vx", &n_vx, &b_n_vx);
   fChain->SetBranchAddress("hfor", &hfor, &b_hfor);
   fChain->SetBranchAddress("top_weight", &top_weight, &b_top_weight);
   fChain->SetBranchAddress("powhegbug_weight", &powhegbug_weight, &b_powhegbug_weight);
   fChain->SetBranchAddress("actualIntPerXing", &actualIntPerXing, &b_actualIntPerXing);
   fChain->SetBranchAddress("averageIntPerXing", &averageIntPerXing, &b_averageIntPerXing);
   fChain->SetBranchAddress("xsec_weight", &xsec_weight, &b_xsec_weight);
   fChain->SetBranchAddress("processed_entries", &processed_entries, &b_processed_entries);
   fChain->SetBranchAddress("generated_entries", &generated_entries, &b_generated_entries);
   fChain->SetBranchAddress("pu_weight", &pu_weight, &b_pu_weight);
   fChain->SetBranchAddress("vxz_weight", &vxz_weight, &b_vxz_weight);
   fChain->SetBranchAddress("OQ_weight", &OQ_weight, &b_OQ_weight);
   fChain->SetBranchAddress("ggF_weight", &ggF_weight, &b_ggF_weight);
   fChain->SetBranchAddress("trigSF_weight", &trigSF_weight, &b_trigSF_weight);
   fChain->SetBranchAddress("Z1_lepplus_weight", &Z1_lepplus_weight, &b_Z1_lepplus_weight);
   fChain->SetBranchAddress("Z1_lepminus_weight", &Z1_lepminus_weight, &b_Z1_lepminus_weight);
   fChain->SetBranchAddress("Z2_lepplus_weight", &Z2_lepplus_weight, &b_Z2_lepplus_weight);
   fChain->SetBranchAddress("Z2_lepminus_weight", &Z2_lepminus_weight, &b_Z2_lepminus_weight);
   fChain->SetBranchAddress("H_m_truth", &H_m_truth, &b_H_m_truth);
   fChain->SetBranchAddress("ZZ_m_truth", &ZZ_m_truth, &b_ZZ_m_truth);
   fChain->SetBranchAddress("Z1_m_truth", &Z1_m_truth, &b_Z1_m_truth);
   fChain->SetBranchAddress("Z2_m_truth", &Z2_m_truth, &b_Z2_m_truth);
   fChain->SetBranchAddress("H_m_constrained", &H_m_constrained, &b_H_m_constrained);
   fChain->SetBranchAddress("Z1_m_constrained", &Z1_m_constrained, &b_Z1_m_constrained);
   fChain->SetBranchAddress("Z2_m_constrained", &Z2_m_constrained, &b_Z2_m_constrained);
   fChain->SetBranchAddress("Z1_chiSq", &Z1_chiSq, &b_Z1_chiSq);
   fChain->SetBranchAddress("Z2_chiSq", &Z2_chiSq, &b_Z2_chiSq);
   fChain->SetBranchAddress("H_m", &H_m, &b_H_m);
   fChain->SetBranchAddress("H_pt", &H_pt, &b_H_pt);
   fChain->SetBranchAddress("H_eta", &H_eta, &b_H_eta);
   fChain->SetBranchAddress("H_phi", &H_phi, &b_H_phi);
   fChain->SetBranchAddress("Z1_m", &Z1_m, &b_Z1_m);
   fChain->SetBranchAddress("Z1_pt", &Z1_pt, &b_Z1_pt);
   fChain->SetBranchAddress("Z1_eta", &Z1_eta, &b_Z1_eta);
   fChain->SetBranchAddress("Z1_phi", &Z1_phi, &b_Z1_phi);
   fChain->SetBranchAddress("Z2_m", &Z2_m, &b_Z2_m);
   fChain->SetBranchAddress("Z2_pt", &Z2_pt, &b_Z2_pt);
   fChain->SetBranchAddress("Z2_eta", &Z2_eta, &b_Z2_eta);
   fChain->SetBranchAddress("Z2_phi", &Z2_phi, &b_Z2_phi);
   fChain->SetBranchAddress("Z1_lepplus_m", &Z1_lepplus_m, &b_Z1_lepplus_m);
   fChain->SetBranchAddress("Z1_lepplus_pt", &Z1_lepplus_pt, &b_Z1_lepplus_pt);
   fChain->SetBranchAddress("Z1_lepplus_eta", &Z1_lepplus_eta, &b_Z1_lepplus_eta);
   fChain->SetBranchAddress("Z1_lepplus_phi", &Z1_lepplus_phi, &b_Z1_lepplus_phi);
   fChain->SetBranchAddress("Z1_lepplus_etcone20", &Z1_lepplus_etcone20, &b_Z1_lepplus_etcone20);
   fChain->SetBranchAddress("Z1_lepplus_etcone20_corr", &Z1_lepplus_etcone20_corr, &b_Z1_lepplus_etcone20_corr);
   fChain->SetBranchAddress("Z1_lepplus_etcone20_final", &Z1_lepplus_etcone20_final, &b_Z1_lepplus_etcone20_final);
   fChain->SetBranchAddress("Z1_lepplus_ptcone20", &Z1_lepplus_ptcone20, &b_Z1_lepplus_ptcone20);
   fChain->SetBranchAddress("Z1_lepplus_ptcone20_final", &Z1_lepplus_ptcone20_final, &b_Z1_lepplus_ptcone20_final);
   fChain->SetBranchAddress("Z1_lepplus_d0", &Z1_lepplus_d0, &b_Z1_lepplus_d0);
   fChain->SetBranchAddress("Z1_lepplus_z0", &Z1_lepplus_z0, &b_Z1_lepplus_z0);
   fChain->SetBranchAddress("Z1_lepplus_d0_sig", &Z1_lepplus_d0_sig, &b_Z1_lepplus_d0_sig);
   fChain->SetBranchAddress("Z1_lepplus_GSF_dp", &Z1_lepplus_GSF_dp, &b_Z1_lepplus_GSF_dp);
   fChain->SetBranchAddress("Z1_lepplus_GSF_p", &Z1_lepplus_GSF_p, &b_Z1_lepplus_GSF_p);
   fChain->SetBranchAddress("Z1_lepplus_charge", &Z1_lepplus_charge, &b_Z1_lepplus_charge);
   fChain->SetBranchAddress("Z1_lepplus_ptCB_nosmearnoscale", &Z1_lepplus_ptCB_nosmearnoscale, &b_Z1_lepplus_ptCB_nosmearnoscale);
   fChain->SetBranchAddress("Z1_lepplus_ptME_nosmearnoscale", &Z1_lepplus_ptME_nosmearnoscale, &b_Z1_lepplus_ptME_nosmearnoscale);
   fChain->SetBranchAddress("Z1_lepplus_ptID_nosmearnoscale", &Z1_lepplus_ptID_nosmearnoscale, &b_Z1_lepplus_ptID_nosmearnoscale);
   fChain->SetBranchAddress("Z1_lepplus_cl_E_calibsmearnoscale", &Z1_lepplus_cl_E_calibsmearnoscale, &b_Z1_lepplus_cl_E_calibsmearnoscale);
   fChain->SetBranchAddress("Z1_lepplus_cl_eta", &Z1_lepplus_cl_eta, &b_Z1_lepplus_cl_eta);
   fChain->SetBranchAddress("Z1_lepplus_cl_phi", &Z1_lepplus_cl_phi, &b_Z1_lepplus_cl_phi);
   fChain->SetBranchAddress("Z1_lepplus_isSA", &Z1_lepplus_isSA, &b_Z1_lepplus_isSA);
   fChain->SetBranchAddress("Z1_lepplus_author", &Z1_lepplus_author, &b_Z1_lepplus_author);
   fChain->SetBranchAddress("Z1_lepplus_type", &Z1_lepplus_type, &b_Z1_lepplus_type);
   fChain->SetBranchAddress("Z1_lepplus_typebkg", &Z1_lepplus_typebkg, &b_Z1_lepplus_typebkg);
   fChain->SetBranchAddress("Z1_lepplus_origin", &Z1_lepplus_origin, &b_Z1_lepplus_origin);
   fChain->SetBranchAddress("Z1_lepplus_originbkg", &Z1_lepplus_originbkg, &b_Z1_lepplus_originbkg);
   fChain->SetBranchAddress("Z1_lepminus_m", &Z1_lepminus_m, &b_Z1_lepminus_m);
   fChain->SetBranchAddress("Z1_lepminus_pt", &Z1_lepminus_pt, &b_Z1_lepminus_pt);
   fChain->SetBranchAddress("Z1_lepminus_eta", &Z1_lepminus_eta, &b_Z1_lepminus_eta);
   fChain->SetBranchAddress("Z1_lepminus_phi", &Z1_lepminus_phi, &b_Z1_lepminus_phi);
   fChain->SetBranchAddress("Z1_lepminus_etcone20", &Z1_lepminus_etcone20, &b_Z1_lepminus_etcone20);
   fChain->SetBranchAddress("Z1_lepminus_etcone20_corr", &Z1_lepminus_etcone20_corr, &b_Z1_lepminus_etcone20_corr);
   fChain->SetBranchAddress("Z1_lepminus_etcone20_final", &Z1_lepminus_etcone20_final, &b_Z1_lepminus_etcone20_final);
   fChain->SetBranchAddress("Z1_lepminus_ptcone20", &Z1_lepminus_ptcone20, &b_Z1_lepminus_ptcone20);
   fChain->SetBranchAddress("Z1_lepminus_ptcone20_final", &Z1_lepminus_ptcone20_final, &b_Z1_lepminus_ptcone20_final);
   fChain->SetBranchAddress("Z1_lepminus_d0", &Z1_lepminus_d0, &b_Z1_lepminus_d0);
   fChain->SetBranchAddress("Z1_lepminus_z0", &Z1_lepminus_z0, &b_Z1_lepminus_z0);
   fChain->SetBranchAddress("Z1_lepminus_d0_sig", &Z1_lepminus_d0_sig, &b_Z1_lepminus_d0_sig);
   fChain->SetBranchAddress("Z1_lepminus_GSF_dp", &Z1_lepminus_GSF_dp, &b_Z1_lepminus_GSF_dp);
   fChain->SetBranchAddress("Z1_lepminus_GSF_p", &Z1_lepminus_GSF_p, &b_Z1_lepminus_GSF_p);
   fChain->SetBranchAddress("Z1_lepminus_charge", &Z1_lepminus_charge, &b_Z1_lepminus_charge);
   fChain->SetBranchAddress("Z1_lepminus_ptCB_nosmearnoscale", &Z1_lepminus_ptCB_nosmearnoscale, &b_Z1_lepminus_ptCB_nosmearnoscale);
   fChain->SetBranchAddress("Z1_lepminus_ptME_nosmearnoscale", &Z1_lepminus_ptME_nosmearnoscale, &b_Z1_lepminus_ptME_nosmearnoscale);
   fChain->SetBranchAddress("Z1_lepminus_ptID_nosmearnoscale", &Z1_lepminus_ptID_nosmearnoscale, &b_Z1_lepminus_ptID_nosmearnoscale);
   fChain->SetBranchAddress("Z1_lepminus_cl_E_calibsmearnoscale", &Z1_lepminus_cl_E_calibsmearnoscale, &b_Z1_lepminus_cl_E_calibsmearnoscale);
   fChain->SetBranchAddress("Z1_lepminus_cl_eta", &Z1_lepminus_cl_eta, &b_Z1_lepminus_cl_eta);
   fChain->SetBranchAddress("Z1_lepminus_cl_phi", &Z1_lepminus_cl_phi, &b_Z1_lepminus_cl_phi);
   fChain->SetBranchAddress("Z1_lepminus_isSA", &Z1_lepminus_isSA, &b_Z1_lepminus_isSA);
   fChain->SetBranchAddress("Z1_lepminus_author", &Z1_lepminus_author, &b_Z1_lepminus_author);
   fChain->SetBranchAddress("Z1_lepminus_type", &Z1_lepminus_type, &b_Z1_lepminus_type);
   fChain->SetBranchAddress("Z1_lepminus_typebkg", &Z1_lepminus_typebkg, &b_Z1_lepminus_typebkg);
   fChain->SetBranchAddress("Z1_lepminus_origin", &Z1_lepminus_origin, &b_Z1_lepminus_origin);
   fChain->SetBranchAddress("Z1_lepminus_originbkg", &Z1_lepminus_originbkg, &b_Z1_lepminus_originbkg);
   fChain->SetBranchAddress("Z2_lepplus_m", &Z2_lepplus_m, &b_Z2_lepplus_m);
   fChain->SetBranchAddress("Z2_lepplus_pt", &Z2_lepplus_pt, &b_Z2_lepplus_pt);
   fChain->SetBranchAddress("Z2_lepplus_eta", &Z2_lepplus_eta, &b_Z2_lepplus_eta);
   fChain->SetBranchAddress("Z2_lepplus_phi", &Z2_lepplus_phi, &b_Z2_lepplus_phi);
   fChain->SetBranchAddress("Z2_lepplus_etcone20", &Z2_lepplus_etcone20, &b_Z2_lepplus_etcone20);
   fChain->SetBranchAddress("Z2_lepplus_etcone20_corr", &Z2_lepplus_etcone20_corr, &b_Z2_lepplus_etcone20_corr);
   fChain->SetBranchAddress("Z2_lepplus_etcone20_final", &Z2_lepplus_etcone20_final, &b_Z2_lepplus_etcone20_final);
   fChain->SetBranchAddress("Z2_lepplus_ptcone20", &Z2_lepplus_ptcone20, &b_Z2_lepplus_ptcone20);
   fChain->SetBranchAddress("Z2_lepplus_ptcone20_final", &Z2_lepplus_ptcone20_final, &b_Z2_lepplus_ptcone20_final);
   fChain->SetBranchAddress("Z2_lepplus_d0", &Z2_lepplus_d0, &b_Z2_lepplus_d0);
   fChain->SetBranchAddress("Z2_lepplus_z0", &Z2_lepplus_z0, &b_Z2_lepplus_z0);
   fChain->SetBranchAddress("Z2_lepplus_d0_sig", &Z2_lepplus_d0_sig, &b_Z2_lepplus_d0_sig);
   fChain->SetBranchAddress("Z2_lepplus_GSF_dp", &Z2_lepplus_GSF_dp, &b_Z2_lepplus_GSF_dp);
   fChain->SetBranchAddress("Z2_lepplus_GSF_p", &Z2_lepplus_GSF_p, &b_Z2_lepplus_GSF_p);
   fChain->SetBranchAddress("Z2_lepplus_charge", &Z2_lepplus_charge, &b_Z2_lepplus_charge);
   fChain->SetBranchAddress("Z2_lepplus_ptCB_nosmearnoscale", &Z2_lepplus_ptCB_nosmearnoscale, &b_Z2_lepplus_ptCB_nosmearnoscale);
   fChain->SetBranchAddress("Z2_lepplus_ptME_nosmearnoscale", &Z2_lepplus_ptME_nosmearnoscale, &b_Z2_lepplus_ptME_nosmearnoscale);
   fChain->SetBranchAddress("Z2_lepplus_ptID_nosmearnoscale", &Z2_lepplus_ptID_nosmearnoscale, &b_Z2_lepplus_ptID_nosmearnoscale);
   fChain->SetBranchAddress("Z2_lepplus_cl_E_calibsmearnoscale", &Z2_lepplus_cl_E_calibsmearnoscale, &b_Z2_lepplus_cl_E_calibsmearnoscale);
   fChain->SetBranchAddress("Z2_lepplus_cl_eta", &Z2_lepplus_cl_eta, &b_Z2_lepplus_cl_eta);
   fChain->SetBranchAddress("Z2_lepplus_cl_phi", &Z2_lepplus_cl_phi, &b_Z2_lepplus_cl_phi);
   fChain->SetBranchAddress("Z2_lepplus_isSA", &Z2_lepplus_isSA, &b_Z2_lepplus_isSA);
   fChain->SetBranchAddress("Z2_lepplus_author", &Z2_lepplus_author, &b_Z2_lepplus_author);
   fChain->SetBranchAddress("Z2_lepplus_type", &Z2_lepplus_type, &b_Z2_lepplus_type);
   fChain->SetBranchAddress("Z2_lepplus_typebkg", &Z2_lepplus_typebkg, &b_Z2_lepplus_typebkg);
   fChain->SetBranchAddress("Z2_lepplus_origin", &Z2_lepplus_origin, &b_Z2_lepplus_origin);
   fChain->SetBranchAddress("Z2_lepplus_originbkg", &Z2_lepplus_originbkg, &b_Z2_lepplus_originbkg);
   fChain->SetBranchAddress("Z2_lepminus_m", &Z2_lepminus_m, &b_Z2_lepminus_m);
   fChain->SetBranchAddress("Z2_lepminus_pt", &Z2_lepminus_pt, &b_Z2_lepminus_pt);
   fChain->SetBranchAddress("Z2_lepminus_eta", &Z2_lepminus_eta, &b_Z2_lepminus_eta);
   fChain->SetBranchAddress("Z2_lepminus_phi", &Z2_lepminus_phi, &b_Z2_lepminus_phi);
   fChain->SetBranchAddress("Z2_lepminus_etcone20", &Z2_lepminus_etcone20, &b_Z2_lepminus_etcone20);
   fChain->SetBranchAddress("Z2_lepminus_etcone20_corr", &Z2_lepminus_etcone20_corr, &b_Z2_lepminus_etcone20_corr);
   fChain->SetBranchAddress("Z2_lepminus_etcone20_final", &Z2_lepminus_etcone20_final, &b_Z2_lepminus_etcone20_final);
   fChain->SetBranchAddress("Z2_lepminus_ptcone20", &Z2_lepminus_ptcone20, &b_Z2_lepminus_ptcone20);
   fChain->SetBranchAddress("Z2_lepminus_ptcone20_final", &Z2_lepminus_ptcone20_final, &b_Z2_lepminus_ptcone20_final);
   fChain->SetBranchAddress("Z2_lepminus_d0", &Z2_lepminus_d0, &b_Z2_lepminus_d0);
   fChain->SetBranchAddress("Z2_lepminus_z0", &Z2_lepminus_z0, &b_Z2_lepminus_z0);
   fChain->SetBranchAddress("Z2_lepminus_d0_sig", &Z2_lepminus_d0_sig, &b_Z2_lepminus_d0_sig);
   fChain->SetBranchAddress("Z2_lepminus_GSF_dp", &Z2_lepminus_GSF_dp, &b_Z2_lepminus_GSF_dp);
   fChain->SetBranchAddress("Z2_lepminus_GSF_p", &Z2_lepminus_GSF_p, &b_Z2_lepminus_GSF_p);
   fChain->SetBranchAddress("Z2_lepminus_charge", &Z2_lepminus_charge, &b_Z2_lepminus_charge);
   fChain->SetBranchAddress("Z2_lepminus_ptCB_nosmearnoscale", &Z2_lepminus_ptCB_nosmearnoscale, &b_Z2_lepminus_ptCB_nosmearnoscale);
   fChain->SetBranchAddress("Z2_lepminus_ptME_nosmearnoscale", &Z2_lepminus_ptME_nosmearnoscale, &b_Z2_lepminus_ptME_nosmearnoscale);
   fChain->SetBranchAddress("Z2_lepminus_ptID_nosmearnoscale", &Z2_lepminus_ptID_nosmearnoscale, &b_Z2_lepminus_ptID_nosmearnoscale);
   fChain->SetBranchAddress("Z2_lepminus_cl_E_calibsmearnoscale", &Z2_lepminus_cl_E_calibsmearnoscale, &b_Z2_lepminus_cl_E_calibsmearnoscale);
   fChain->SetBranchAddress("Z2_lepminus_cl_eta", &Z2_lepminus_cl_eta, &b_Z2_lepminus_cl_eta);
   fChain->SetBranchAddress("Z2_lepminus_cl_phi", &Z2_lepminus_cl_phi, &b_Z2_lepminus_cl_phi);
   fChain->SetBranchAddress("Z2_lepminus_isSA", &Z2_lepminus_isSA, &b_Z2_lepminus_isSA);
   fChain->SetBranchAddress("Z2_lepminus_author", &Z2_lepminus_author, &b_Z2_lepminus_author);
   fChain->SetBranchAddress("Z2_lepminus_type", &Z2_lepminus_type, &b_Z2_lepminus_type);
   fChain->SetBranchAddress("Z2_lepminus_typebkg", &Z2_lepminus_typebkg, &b_Z2_lepminus_typebkg);
   fChain->SetBranchAddress("Z2_lepminus_origin", &Z2_lepminus_origin, &b_Z2_lepminus_origin);
   fChain->SetBranchAddress("Z2_lepminus_originbkg", &Z2_lepminus_originbkg, &b_Z2_lepminus_originbkg);
   fChain->SetBranchAddress("filename", &filename, &b_filename);
}

Bool_t YieldEstimator::Notify()
{
   TFile *file = 0;
   TChain *chain = dynamic_cast<TChain*>(fChain);

   if (chain) {
      // We are running locally...
      file = chain->GetFile();
   } else {
      // We are running on PROOF:
      file = fChain->GetCurrentFile();
   }

   Info("Notify", "Opening %s", file->GetName());

   fChain->GetTree()->GetEntry(0);
   m_current_run = run;

   Info("Notify", "        (run = %d)", (Int_t)run);

   // check if it is an interesting run
   if (m_sample.find(run) == m_sample.end()) {
      Warning("Notify", "        ---> not interesting run, skipping!");
      m_interesting = kFALSE;

      return kTRUE;
   } else {
      m_interesting = kTRUE;
   }


   // add if new
   if (m_sample[run].generated < 0) {
      Info("Notify", "        ---> was a new run, adding to the list...");
      for (Int_t i = 0; i < N_CHANNELS; i++) {
         Int_t chan_id[4] = {0, 1, 3, 2};
         Int_t chanid = chan_id[i];
         m_sample[run].h_H_m.push_back(new TH1F(TString::Format("raw_%d_H_m_%d", (Int_t)run, chanid), "", 2000, 0, 1000));
         m_sample[run].h_H_m_constrained.push_back(new TH1F(TString::Format("raw_%d_H_m_constrained_%d", (Int_t)run, chanid), "", 2000, 0, 1000));
         m_sample[run].h_Z1_m.push_back(new TH1F(TString::Format("raw_%d_Z1_m_%d", (Int_t)run, chanid), "", 240, 30, 150));
         m_sample[run].h_Z2_m.push_back(new TH1F(TString::Format("raw_%d_Z2_m_%d", (Int_t)run, chanid), "", 150, 0, 150));

         m_sample[run].h_H_m[i]->SetDirectory(m_output_raw);
         m_sample[run].h_H_m_constrained[i]->SetDirectory(m_output_raw);
         m_sample[run].h_Z1_m[i]->SetDirectory(m_output_raw);
         m_sample[run].h_Z2_m[i]->SetDirectory(m_output_raw);
      }

      if (m_sample[run].isSignal == kFALSE) {
         m_sample[run].xsec = CrossSections::bkgCrossSection(run, m_COM_energy, kFALSE);
         if (m_sample[run].xsec <= 0.) {
            Warning("Notify", "        ---> no cross section retrieved (but it should be a background!!!)");
         }
      } else {
         if (m_sample[run].type == MC::ggF) {
            m_sample[run].xsec = higgs_xsector.higgs4lxsecGGF(m_sample[run].mass, m_COM_energy);
         } else if (m_sample[run].type == MC::VBF) {
            m_sample[run].xsec = higgs_xsector.higgs4lxsecVBF(m_sample[run].mass, m_COM_energy);
         } else if (m_sample[run].type == MC::WH) {
            m_sample[run].xsec = higgs_xsector.higgs4lxsecWH(m_sample[run].mass, m_COM_energy);
         } else if (m_sample[run].type == MC::ZH) {
            m_sample[run].xsec = higgs_xsector.higgs4lxsecZH(m_sample[run].mass, m_COM_energy);
         }

         if (m_sample[run].xsec <= 0.) {
            Warning("Notify", "        ---> no cross section retrieved (but it should be a signal with mass %lf GeV!!!", m_sample[run].mass);
         }
      }
   }

   // update generated entries
   TH1D *gen_histo = dynamic_cast<TH1D*>(file->Get("generatedEntriesHisto"));

   if (gen_histo) {
      m_sample[run].generated += gen_histo->GetBinContent(BIN_FOR_GENERATED);
   } else {
      Error("Notify", "        ---> unable to find the generated entries histogram !!!");
   }


   return kTRUE;
}

#endif // #ifdef YieldEstimator_cxx
