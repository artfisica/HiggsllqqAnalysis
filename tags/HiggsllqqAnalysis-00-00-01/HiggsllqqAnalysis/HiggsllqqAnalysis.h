#ifndef __HiggsAnalysis_h__
#define __HiggsAnalysis_h__

/** @class HiggsllqqAnalysis HiggsllqqAnalysis.h

    Code to perform SM H -> ZZ(*) -> qqll analysis.

    @start  date 15/08/2012 
    @update date 29/10/2012
*/

#include <TEfficiency.h>

#include "HiggsAnalysis/HiggsAnalysis.h"

// ATLAS tools
#include "GoodRunsLists/TGoodRunsListReader.h"
#include "TrigRootAnalysis/TrigDecisionToolD3PD.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "egammaAnalysisUtils/egammaSFclass.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "MuonEfficiencyCorrections/AnalysisMuonConfigurableScaleFactors.h"
#include "MuonMomentumCorrections/SmearingClass.h"
#include "TrigMuonEfficiency/LeptonTriggerSF.h"
#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include "HiggsZZ4lUtils/IsEMPlusPlusH4lDefs.h"
#include "HiggsZZ4lUtils/BkgCrossSection.h"
#include "HiggsZZ4lUtils/HiggsCrossSection.h"
#include "HiggsZZ4lUtils/MzzWeightFromMCFM.h"
#include "TrigMuonEfficiency/MuonTriggerMatching.h"
#include "TrigMuonEfficiency/ElectronTriggerMatching.h"
#include "egammaAnalysisUtils/MultiLeptonDefs.h"
#include "egammaAnalysisUtils/VertexPositionReweightingTool.h"

// Higgs[llqq-4lep]Analysis common tools
#include "HiggsllqqAnalysis/ggFReweighting.h"
#include "HiggsllqqAnalysis/CommonTools.h"
#include "HiggsllqqAnalysis/ChargedLepton.h"
#include "HiggsllqqAnalysis/Jet.h"
#include "HiggsllqqAnalysis/Dilepton.h"
#include "HiggsllqqAnalysis/CutFlowTool.h"
#include "HiggsllqqAnalysis/DataPeriodTool.h"


// HiggsqqllAnalysis needs
#include "egammaAnalysisUtils/IsEMPlusPlusDefs.h"
#include "JetTagAlgorithms/MV1.h"
#include "ApplyJetCalibration/ApplyJetCalibration.h"
#include "ApplyJetResolutionSmearing/ApplyJetSmearing.h"
#include "HiggsllqqAnalysis/JetKinematicFitter.h"
#include "HiggsllqqAnalysis/HforToolD3PD.h"

// TestSelection needs (redundancies to be removed)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TLorentzVector.h>
#include <vector>
#include <TStopwatch.h>
#include <TH1D.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
//#include "TestSelection/JetKinematicFitter.h"
//#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"


namespace DataPeriod {
  enum {
    y2011_A, // data11
    y2011_B,
    y2011_C,
    y2011_D,
    y2011_E,
    y2011_F,
    y2011_G,
    y2011_H,
    y2011_I,
    y2011_J,
    y2011_K,
    y2011_L,
    y2011_M,
    y2012_A, // data12
    y2012_B, // data12
    // DO NOT USE VALUES BELOW THIS LINE FOR COMPARISON (e.g. AI < JM makes no sense)
    y2011_BD, // MC
    y2011_EH,
    y2011_LM,
    y2012_AllYear,
  };
}

namespace HllqqCutFlow {
  enum {
    HFOR,
    GRL,
    larError,
    Trigger,
    Vertex,
    METcleaning,
    LArHole,
    NumberOfLeptons,
    OppositeSign,
    PtLeptons,
    TriggerConsistency,
    MET,
    TwoJets,
    NumTagJets,
    DileptonMass,
    DiJetMass,
  };
}

namespace HllqqMuonQuality {
  enum {
    family,
    quality,
    cosmic,
    eta,
    pt,
    MCP,
    z0,
    overlap,
  };
}

namespace HllqqElectronQuality {
  enum {
    family,
    algorithm,
    quality,
    eta,
    Et,
    objectquality,
    z0,
    overlap,
  };
}

namespace HllqqJetQuality {
  enum {
    jetCleaning,
    kinematics,
    Pileup,
    overlap,
  };
}

namespace HllqqGRL {
  enum {
    // offset in the vector of GRLs: this number, plus the output of getChannel(), gives the location of the GRL to be read for that channel and that year
    // the difference between two consecutive members of this enum it must be equal to the number of entries in the enum of the class HiggsllqqAnalysis (i.e. { MU2, MUE, E2 } -> 3)
    data11 = 0,
    data12 = 3,
  };
}

namespace HllqqSystematics {
  class ChargedLepton {
  public:
    Float_t ptCB_nosmearnoscale;
    Float_t ptME_nosmearnoscale;
    Float_t ptID_nosmearnoscale;
    Float_t cl_E_calibsmearnoscale;
    Float_t cl_eta;
    Float_t cl_phi;
    
    ChargedLepton() {
      ptCB_nosmearnoscale = -9999.9;
      ptME_nosmearnoscale = -9999.9;
      ptID_nosmearnoscale = -9999.9;
      cl_E_calibsmearnoscale = -9999.9;
      cl_eta = -9999.9;
      cl_phi = -9999.9;
    }
    ~ChargedLepton() {}
  };
};

class HiggsllqqAnalysis : public HiggsAnalysis {
 public:
  
  enum {
    MU2,
    MUE,
    E2,
  };
  
  // constructor, destructor
  HiggsllqqAnalysis(TTree * /*tree*/ = 0) {
  }
  ~HiggsllqqAnalysis();
  
  // options setters
  virtual void setAnalysisVersion(TString val) {
    m_analysis_version = val;
  }
  virtual void setSmearing(Bool_t val) {
    m_doSmearing = val;
  }
  virtual void setTopoIso(Bool_t val) {
    m_useTopoIso = val;
  }
  virtual void setMuonFamily(Int_t val) {
    m_muonFamily = val;
  }
  virtual void setElectronFamily(Int_t val) {
    m_electronFamily = val;
  }
  virtual void setOutputFile(TString val) {
    m_outputFileName = val;
  }
  
 protected:
  
  // pointers to ATLAS tools
  Root::TGoodRunsListReader *m_GRL;
  Root::TPileupReweighting *m_PileupReweighter;
  VertexPositionReweightingTool *m_VertexPositionReweighter;
  D3PD::TrigDecisionToolD3PD *m_TrigDecisionToolD3PD;
  eg2011::EnergyRescaler *m_ElectronEnergyRescaler;
  egammaSFclass *m_ElectronEffSF;
  Analysis::AnalysisMuonConfigurableScaleFactors *m_MuonEffSF;
  Analysis::AnalysisMuonConfigurableScaleFactors *m_MuonEffSFCalo;
  Analysis::AnalysisMuonConfigurableScaleFactors *m_MuonEffSFSA;
  MuonSmear::SmearingClass *m_MuonSmearer;
  LeptonTriggerSF *m_MuonTrigSF;
  ggFReweighting *m_ggFReweighter;
  TriggerNavigationVariables m_trigNavVar;
  MuonTriggerMatching *m_MuonTriggerMatchTool;
  ElectronTriggerMatching *m_ElectronTriggerMatchTool;
  TH2F *m_smearD0[3];
  TRandom3 m_smearD0_rand;
  TAxis *m_smearD0_x;
  TAxis *m_smearD0_y;
  
  // containers for leptons, Jets, dileptons, in the event
  std::vector<Analysis::ChargedLepton *> m_Muons;
  std::vector<Analysis::ChargedLepton *> m_Electrons;
  std::vector<Analysis::Jet *> m_Jets;
  std::vector<Analysis::ChargedLepton *> m_GoodMuons;
  std::vector<Analysis::ChargedLepton *> m_GoodElectrons;
  std::vector<Analysis::Jet *> m_GoodJets;
  std::vector<Analysis::Dilepton *> m_Dileptons;
  
  // utility maps
  std::map<UInt_t, Float_t> m_CrossSection;
  std::map<UInt_t, Int_t> m_SignalSampleMass;
  
  // channels to be processed and cutflows
  std::vector<Int_t> m_Channels;
  std::vector<Analysis::CutFlowTool> m_EventCutflow;
  std::vector<Analysis::CutFlowTool> m_EventCutflow_rw;
  std::vector<Analysis::CutFlowTool> m_ElectronCutflow;
  std::vector<Analysis::CutFlowTool> m_MuonCutflow;
  std::vector<Analysis::CutFlowTool> m_JetCutflow;
  
  
  // functions dealing with tools' initialization / config
  virtual Bool_t initialize_tools();
  virtual Bool_t change_input();
  virtual Bool_t execute_tools(Long64_t entry);
  virtual Bool_t finalize_tools();
  
  // functions implementing the actual event loop
  virtual Bool_t initialize_analysis();
  virtual Bool_t execute_analysis();
  virtual Bool_t finalize_analysis();
  
  // event selection
  virtual Bool_t hasPowHegZZBug(); // to be removed and/or change!!
  virtual Int_t getLastCutPassed();
  virtual Bool_t passesGRL();
  virtual Bool_t hasGoodVertex();
  virtual Int_t getNumberOfGoodVertices();
  virtual std::vector<TString> getListOfAlternativeTriggers(TString sequence);
  virtual Bool_t passesTrigger();
  virtual Bool_t passesSingleMuonTrigger();
  virtual Bool_t passesDiMuonTrigger();
  virtual Bool_t passesSingleElectronTrigger();
  virtual Bool_t passesDiElectronTrigger();
  virtual Bool_t passesElectronMuonTrigger();
  virtual TString getSingleMuonTriggerName();
  virtual TString getDiMuonTriggerName();
  virtual TString getSingleElectronTriggerName();
  virtual TString getDiElectronTriggerName();
  virtual TString getElectronMuonTriggerName();
  virtual void applyChanges(Analysis::ChargedLepton *lep);
  virtual void applyChanges(Analysis::Jet *jet);
  virtual void getMuons(D3PDReader::MuonD3PDObject *mu_branch, Int_t family);
  virtual void getElectrons(D3PDReader::ElectronD3PDObject *el_branch, Int_t family);
  virtual void getJets(D3PDReader::JetD3PDObject *jet_branch); // To include Jet Family??
  virtual void getGoodMuons();
  virtual void getGoodElectrons();
  virtual void getGoodJets();
  virtual void getGoodLeptons();
  virtual void getDileptons();
  virtual void getGoodObjects();
  
  // object selection
  virtual Bool_t isMCPMuon(Analysis::ChargedLepton *lep);
  virtual Bool_t isGood(Analysis::ChargedLepton *lep);
  Bool_t IsGoodElectron(Analysis::ChargedLepton *lep);
  Bool_t IsGoodMuon(Analysis::ChargedLepton *lep);
  Bool_t isGoodJet(Analysis::Jet *jet);
  
  // utility functions for the selection
  virtual Bool_t isMC() {
    return ntuple->mc.n.IsAvailable();
  }
  virtual TString analysis_version() {
    return m_analysis_version;
  }
  virtual Bool_t doSmearing() {
    return m_doSmearing;
  }
  virtual Bool_t useTopoIso() {
    return m_useTopoIso;
  }
  virtual Int_t getMuonFamily() {
    return m_muonFamily;
  }
  virtual Int_t getElectronFamily() {
    return m_electronFamily;
  }
  virtual Int_t getTriggerInfo(TString chain);
  virtual Int_t getPeriod();
  virtual void setChannel(Int_t chan) {
    m_thisChannel = chan;
  }
  virtual Int_t getChannel() {
    return m_thisChannel;
  }
  virtual Float_t getD0SmearSigma(Int_t index_of_lepton, Int_t nBL, Float_t pt, Float_t eta);
  
  // event weighting helpers
  virtual Float_t getEventWeight();
  virtual Float_t getPileupWeight();
  virtual Float_t getVertexZWeight();
  virtual Float_t getLeptonWeight(Analysis::ChargedLepton *lep);
  virtual Float_t getSFWeight();
  virtual Float_t getggFWeight();
  
  // tools for cross section computation
  void initCrossSections(); // to be change!!
  Float_t getCrossSectionWeight();
  Float_t getTruthHiggsMass();
  virtual Float_t getTruthHiggsPt(); // to be change!!
    
  // option resume
  virtual void printAllOptions() {
    Info("printAllOptions", "======== options =======");
    Info("printAllOptions", "analysis_version = %s", m_analysis_version.Data());
    Info("printAllOptions", "doSmearing       = %s", m_doSmearing ? "kTRUE" : "kFALSE");
    Info("printAllOptions", "useTopoIso       = %s", m_useTopoIso ? "kTRUE" : "kFALSE");
    Info("printAllOptions", "muonFamily       = %d", m_muonFamily);
    Info("printAllOptions", "electronFamily   = %d", m_electronFamily);
    Info("printAllOptions", "outputFileName   = %s", m_outputFileName.Data());
    Info("printAllOptions", "========================");
  }
  
  
  // Methods from the qqll 2011 analysis
  Bool_t ApplyChangesMuon(Analysis::ChargedLepton *lep);
  Bool_t ApplyChangesElectron(Analysis::ChargedLepton *lep);
  Bool_t ApplyChangesJet(Analysis::Jet *jet);
  
  Int_t GetLastCutPassed();
  Bool_t GetGoodObjects();
  void LoadGRL();
  Bool_t PassesGRL();
  Bool_t JetInHole();
  Bool_t NotMETclean();
  Float_t getCorrectMETValue();
  Bool_t hasGoodMET();

    void InitReducedNtuple();
    void ResetReducedNtupleMembers();
    void FillReducedNtuple(Int_t cut, UInt_t channel);
    
  Bool_t GetDoLowMass() { return m_dolowmass; }
  Bool_t GetSysStudy() { return m_sysstudy; }
  Bool_t IsConsistentPt();
  Bool_t IsConsistentWithTrigger();
  Float_t GetMV1value(Analysis::Jet *jet);
  pair <Int_t,Double_t> GetFlavour(Analysis::Jet *jet);
  Bool_t isHeavyJet(Int_t pdg);
  Bool_t isLightJet(Int_t pdg);
  Bool_t isGluonJet(Int_t pdg);
  Float_t GetggFWeight();
  void InitMasses();
  Float_t GetTruthHiggsPt();
  
  // Trigger SF 
  std::pair<double, double> getCandidateTriggerSF(Int_t option);
  
  //Method to count the number of that events into the GoodJets vector.
  Int_t GetNumOfTags();
  
  //2 Methods to evaluate the impact of the Jet kinematic cuts in the selection
  Bool_t GoodPtJets();
  Bool_t GoodEtaJets();
  
  //Method to evaluate the isolation mu-jet
  Bool_t MuonJetOR();
  
  //Method to find the best DiJets using JetKinematicFitter
  Bool_t JetKinematicFitterResult();
  
  //Method to calculate the DiJet invariant mass for the tagged Jets!
  Bool_t JetDimassTagged();
  
  Bool_t GetDoQCDSelection() { 
    if(!isMC()) return m_doqcdselection; 
    return kFALSE;
  }
  
  void SetDoLowMass(Bool_t val) { m_dolowmass = val; }
  void SetSysStudy(Bool_t val) { m_sysstudy = val; }
  void SetDoQCDSelection(Bool_t val) { m_doqcdselection = val; }
  void SortIndex(std::vector<Analysis::Jet *> &vec);
  void SortIndex(std::vector<Analysis::ChargedLepton *> &vec);
  
  //Printing the Event tables
  void PrintCutFlowEvents(Int_t cut, Int_t event, Int_t count);
  
  //Method to calculate the Cleaning of the jet
  Bool_t isBadLooser(Analysis::Jet *jet);
  
  // Method to fill vectors to calculate the HFOR value
  void FillHFORvariables();
  
  
 protected:
  Bool_t m_called_getGoodLeptons;
  Bool_t m_called_getGoodJets;
  Bool_t m_called_getGoodObjects;
  Long64_t m_entriesInChain;
  Long64_t m_processedEntries;
  
 private:
  TString m_analysis_version;
  Bool_t m_doSmearing;
  Bool_t m_useTopoIso;
  Int_t m_muonFamily;
  Int_t m_electronFamily;
  Int_t m_thisChannel;
  TH1D *m_generatedEntriesHisto;
  std::map<TString, TH1F*> m_truthHistos; // To be update or remove
  TEfficiency *m_selectionEfficiencyVsNvx[4]; // To be update or remove
  TTree *m_TreeCutflow;
  TString m_outputFileName;
  TFile *m_outputFile;
  
  
  //Beginning
 protected:
  
  //JES AND JER TOOL
  JetCalibrationTool *myJES;
  JetSmearingTool *myJER;
  
  /// ggF reweighting tool
  ggFReweighting *fggFReweighter;
  
  // Jet kinematic fittr
  JetKinematicFitter *m_jetkinematicfitter;
  
  //HFOR TOOL
  HforToolD3PD *hforTool;
  
  // FLS tree
  TTree *m_reduced_ntuple;
  
  Int_t m_lep_chargeproduct;
  std::vector<Float_t> *m_lep_m;
  std::vector<Float_t> *m_lep_pt;
  std::vector<Float_t> *m_lep_eta;
  std::vector<Float_t> *m_lep_phi;
  std::vector<Int_t> *m_lep_charge;
  std::vector<Float_t> *m_lep_d0;
  std::vector<Float_t> *m_lep_sigd0;
  std::vector<Float_t> *m_lep_trackiso;
  std::vector<Float_t> *m_lep_caloiso;
  std::vector<Int_t> *m_lep_quality;
  std::vector<Float_t> *m_jets_m;
  std::vector<Float_t> *m_jets_pt;
  std::vector<Float_t> *m_jets_eta;
  std::vector<Float_t> *m_jets_eta_det;
  std::vector<Float_t> *m_jets_phi;
  std::vector<Float_t> *m_jets_MV1;
  std::vector<Float_t> *m_jets_flavortruth;
  std::vector<Float_t> *m_jets_jvtxf;
  std::vector<Int_t>   *m_jets_nTrk;
  std::vector<Float_t> *m_jets_width;
  std::vector<Int_t>   *m_jets_flavorpdg;
  std::vector<Double_t>   *m_jets_Epdg;
  UInt_t m_run, m_event;
  Float_t m_weight, m_mu, m_truthH_pt, m_ggFweight;
  Int_t m_cut, m_channel, m_qcdevent, m_NPV, m_Entries, m_HFOR;
  Float_t m_met_met, m_met_sumet, m_met_phi;
  Float_t m_trig_SF,m_trig_SF2,m_trig_SFC;
  Int_t m_trig_flag;
  
  //Truth quark information
  std::vector<Float_t> *m_quark_m;
  std::vector<Float_t> *m_quark_pt;
  std::vector<Int_t> *m_quark_pdg;
  std::vector<Float_t> *m_quark_eta;
  std::vector<Float_t> *m_quark_phi;
  std::vector<Float_t> *m_quark_E;
  
  Int_t                mc_n;
  std::vector<float>   *mc_pt;
  std::vector<int>     *mc_pdgId;
  std::vector<float>   *mc_m;
  std::vector<float>   *mc_eta;
  std::vector<float>   *mc_phi;
  std::vector<int>     *mc_status;
  std::vector<int>     *mc_vx_barcode;
  std::vector<vector<int> > *mc_child_index;
  std::vector<vector<int> > *mc_parent_index;   
  
 private:
  Int_t fChannel;
  std::string fOutputFile;
  Bool_t met_type_RefFinal;
  Bool_t met_type_LocHadTopo;
  Bool_t m_dolowmass;
  Bool_t m_sysstudy;
  Bool_t m_doqcdselection;
  std::map<UInt_t, Int_t> fSampleMass;

  float Mean_jets;
  float good_events;
  
 public:
  ClassDef(HiggsllqqAnalysis, 0);
};

#endif
