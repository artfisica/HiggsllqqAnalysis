#ifndef __HiggsAnalysis_h__
#define __HiggsAnalysis_h__

/** @class HiggsllqqAnalysis HiggsllqqAnalysis.h
    
    Code to perform SM H -> ZZ(*) -> qqll analysis.

    @start  date 08/15/2012 
    @update date 09/09/2013
*/

#include <TEfficiency.h>

#include "HiggsAnalysis/HiggsAnalysis.h"

// ATLAS tools
#include "GoodRunsLists/TGoodRunsListReader.h"
#include "TrigRootAnalysis/TrigDecisionToolD3PD.h"
#include "PileupReweighting/TPileupReweighting.h"
#include "egammaAnalysisUtils/egammaSFclass.h"
#include "MuonEfficiencyCorrections/AnalysisMuonConfigurableScaleFactors.h"
#include "MuonMomentumCorrections/SmearingClass.h"
#include "TrigMuonEfficiency/LeptonTriggerSF.h"
#include "egammaAnalysisUtils/CaloIsoCorrection.h"
#include "HiggsZZ4lUtils/BkgCrossSection.h"
#include "HiggsZZ4lUtils/HiggsCrossSection.h"
#include "HiggsZZ4lUtils/MzzWeightFromMCFM.h"
#include "TrigMuonEfficiency/MuonTriggerMatching.h"
#include "TrigMuonEfficiency/ElectronTriggerMatching.h"
#include "egammaAnalysisUtils/VertexPositionReweightingTool.h"
#include "HiggsZZ4lUtils/McOverlapRemoval.h"
#include "egammaEvent/egammaPIDdefs.h"
#include "egammaAnalysisUtils/EnergyRescalerUpgrade.h"

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
#include "JetTagAlgorithms/MV1c.h"
#include "ApplyJetCalibration/ApplyJetCalibration.h"
#include "ApplyJetResolutionSmearing/ApplyJetSmearing.h"
#include "HiggsllqqAnalysis/JetKinematicFitter.h"
#include "HiggsllqqAnalysis/HforToolD3PD.h"


// TestSelection needs (redundancies to be removed)
#include <fstream>
#include <TStyle.h>
#include <TTree.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TLorentzVector.h>
#include <vector>
#include <TStopwatch.h>
#include <TH1.h>
#include <TH1D.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Config.h"
#include <TMath.h>
#include "Math/Interpolator.h"
#include "CalibrationDataInterface/CalibrationDataInterfaceROOT.h"

using namespace std;


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
    y2012_C, // data12
    y2012_D, // data12
    y2012_E, // data12
    // DO NOT USE VALUES BELOW THIS LINE FOR COMPARISON (e.g. AI < JM makes no sense)
    y2011_BD, // MC
    y2011_EH,
    y2011_LM,
    y2012_AllYear,
  };
}


namespace HllqqCutFlow {
  enum {
    Entries,
    HFOR,
    GRL,
    larError,
    Trigger,
    Vertex,
    METcleaning,
    LArHole,
    NumberOfLeptons,
    TriggerConsistency,
    OppositeSign,
    TwoJets,
    DileptonMass,
    MET,
    NumTagJets,
    DiJetMass,
  };
}


namespace HllqqCutFlow0tag {
  enum {
    NumTagJets0,
    PtLeadingJet0,
    DiJetMass0,
    //    JetPt0,
    //    Ptleptons0,
    //    Dphileptons0,
  };
}


namespace HllqqCutFlow1tag {
  enum {
    NumTagJets1,
    PtLeadingJet1,
    DiJetMass1,
    //    JetPt0,
    //    Ptleptons0,
    //    Dphileptons0,
  };
}


namespace HllqqCutFlow2tag {
  enum {
    NumTagJets2,
    PtLeadingJet2,
    DiJetMass2,
    //    JetPt0,
    //    Ptleptons0,
    //    Dphileptons0,
  };
}


namespace HllqqMuonQuality {
  enum {
    family,
    quality,
    eta,
    pt,
    MCP,
    cosmic,
    z0,
    d0Sig, 
    Isolation,
    overlap,
    medium,
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
    d0Sig, 
    Isolation,
    overlap,
    medium,
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
      ptCB_nosmearnoscale    = -9999.9;
      ptME_nosmearnoscale    = -9999.9;
      ptID_nosmearnoscale    = -9999.9;
      cl_E_calibsmearnoscale = -9999.9;
      cl_eta = -9999.9;
      cl_phi = -9999.9;
    }
    ~ChargedLepton() {}
  };
};


//TestSelection Candidate Structures
typedef struct {
  unsigned int runnumber;
  unsigned int eventnumber;
  int   istagged;
  int   channel;
  int   isqcdevent;
  int   low_event;
  int   n_jets;
  int   n_b_jets;
  //////////////////////
  float lep1_m;
  float lep1_pt;
  float lep1_eta;
  float lep1_phi;
  float lep1_charge;
  float lep1_caloiso;
  float lep1_trackiso;
  float lep1_d0;
  float lep1_sigd0;
  int   lep1_quality;
  float lep2_m;
  float lep2_pt;
  float lep2_eta;
  float lep2_phi;
  float lep2_charge;
  float lep2_caloiso;
  float lep2_trackiso;
  float lep2_d0;
  float lep2_sigd0;
  int   lep2_quality;
  //
  float lepZ_m;
  float lepZ_pt;
  float lepZ_eta;
  float lepZ_phi;
  //////////////////////
  int   Jet1_KF_index;
  int   Jet2_KF_index;
  int   Jet1_LJ_index;
  int   Jet2_LJ_index;
  int   Jet1_BP_index;
  int   Jet2_BP_index;
  //////////////////////
  float realJ1_KF_m;
  float realJ1_KF_pt;
  float realJ1_KF_eta;
  float realJ1_KF_eta_det;
  float realJ1_KF_phi;
  float realJ1_KF_flavortruth;
  int   realJ1_KF_pdg;
  float realJ1_KF_jvf;
  int   realJ1_KF_ntrk;
  float realJ1_KF_width;
  int   realJ1_KF_ntrk12;
  float realJ1_KF_width12;
  float realJ1_KF_MV1;
  float realJ1_KF_Fisher;
  float realJ1_KF_LL;
  float realJ1_KF_LLMIX;
  float realJ2_KF_m;
  float realJ2_KF_pt;
  float realJ2_KF_eta;
  float realJ2_KF_eta_det;
  float realJ2_KF_phi;
  float realJ2_KF_flavortruth;
  int   realJ2_KF_pdg;
  float realJ2_KF_jvf;
  int   realJ2_KF_ntrk;
  float realJ2_KF_width;
  int   realJ2_KF_ntrk12;
  float realJ2_KF_width12;
  float realJ2_KF_MV1;
  float realJ2_KF_Fisher;
  float realJ2_KF_LL;
  float realJ2_KF_LLMIX;
  //
  float ll_2_KF_jets;
  float realZ_KF_m;
  float realZ_KF_pt;
  float realZ_KF_eta;
  float realZ_KF_phi;
  float realH_KF_m;
  float realH_KF_pt;
  float realH_KF_eta;
  float realH_KF_phi;
  //
  float corrJ1_KF_m;
  float corrJ1_KF_pt;
  float corrJ1_KF_eta;
  float corrJ1_KF_eta_det;
  float corrJ1_KF_phi;
  float corrJ1_KF_flavortruth;
  float corrJ1_KF_Fisher;
  float corrJ2_KF_m;
  float corrJ2_KF_pt;
  float corrJ2_KF_eta;
  float corrJ2_KF_eta_det;
  float corrJ2_KF_phi;
  float corrJ2_KF_flavortruth;
  float corrJ2_KF_Fisher;
  //
  float ll_2_KF_jets_corr;
  float corrZ_KF_m;
  float corrZ_KF_pt;
  float corrZ_KF_eta;
  float corrZ_KF_phi;
  float corrH_KF_m;
  float corrH_KF_pt;
  float corrH_KF_eta;
  float corrH_KF_phi;
  //////////////////////
  float realJ1_LJ_m;
  float realJ1_LJ_pt;
  float realJ1_LJ_eta;
  float realJ1_LJ_eta_det;
  float realJ1_LJ_phi;
  float realJ1_LJ_flavortruth;
  int   realJ1_LJ_pdg;
  float realJ1_LJ_jvf;
  int   realJ1_LJ_ntrk;
  float realJ1_LJ_width;
  int   realJ1_LJ_ntrk12;
  float realJ1_LJ_width12;
  float realJ1_LJ_MV1;
  float realJ1_LJ_Fisher;
  float realJ1_LJ_LL;
  float realJ1_LJ_LLMIX;
  float realJ2_LJ_m;
  float realJ2_LJ_pt;
  float realJ2_LJ_eta;
  float realJ2_LJ_eta_det;
  float realJ2_LJ_phi;
  float realJ2_LJ_flavortruth;
  int   realJ2_LJ_pdg;
  float realJ2_LJ_jvf;
  int   realJ2_LJ_ntrk;
  float realJ2_LJ_width;
  int   realJ2_LJ_ntrk12;
  float realJ2_LJ_width12;
  float realJ2_LJ_MV1;
  float realJ2_LJ_Fisher;
  float realJ2_LJ_LL;
  float realJ2_LJ_LLMIX;
  //
  float realZ_LJ_m;
  float realZ_LJ_pt;
  float realZ_LJ_eta;
  float realZ_LJ_phi;
  float realH_LJ_m;
  float realH_LJ_pt;
  float realH_LJ_eta;
  float realH_LJ_phi; 
  //////////////////////
  float realJ1_BP_m;
  float realJ1_BP_pt;
  float realJ1_BP_eta;
  float realJ1_BP_eta_det;
  float realJ1_BP_phi;
  float realJ1_BP_flavortruth;
  int   realJ1_BP_pdg;
  float realJ1_BP_jvf;
  int   realJ1_BP_ntrk;
  float realJ1_BP_width;
  int   realJ1_BP_ntrk12;
  float realJ1_BP_width12;
  float realJ1_BP_MV1;
  float realJ1_BP_Fisher;
  float realJ1_BP_LL;
  float realJ1_BP_LLMIX;
  float realJ2_BP_m;
  float realJ2_BP_pt;
  float realJ2_BP_eta;
  float realJ2_BP_eta_det;
  float realJ2_BP_phi;
  float realJ2_BP_flavortruth;
  int   realJ2_BP_pdg;
  float realJ2_BP_jvf;
  int   realJ2_BP_ntrk;
  float realJ2_BP_width;
  int   realJ2_BP_ntrk12;
  float realJ2_BP_width12;
  float realJ2_BP_MV1;
  float realJ2_BP_Fisher;
  float realJ2_BP_LL;
  float realJ2_BP_LLMIX;
  //
  float realZ_BP_m;
  float realZ_BP_pt;
  float realZ_BP_eta;
  float realZ_BP_phi;
  float realH_BP_m;
  float realH_BP_pt;
  float realH_BP_eta;
  float realH_BP_phi; 
  ////////////////////////
  float dR_ll;
  float dPhi_ll;
  //
  float dPhi_KF_jj;
  float dPhi_LJ_jj;
  float dPhi_BP_jj;
  float dR_KF_jj;
  float dR_LJ_jj;
  float dR_BP_jj;
  //
  float dPhi_KF_ZZ;
  float dPhi_LJ_ZZ;
  float dPhi_BP_ZZ;
  float dR_KF_ZZ;
  float dR_LJ_ZZ;
  float dR_BP_ZZ;
  ///////////////////////
  float chisquare;
  float mu;
  float met;
  int   NPV;
  int   HFOR;
  float sumet;
  float btagSF;
  int   Entries;
  float truthH_pt;
  float weight;
  float SFWeight;
  float ggFweight;
  float EventWeight;
  float PileupWeight;
  float VertexZWeight;
  float DPhijjZWeight;
  float TriggerSFWeight;
  //////////////////////
  float xWin_44p_4var;
  float yWin_44p_4var;
  float zWin_44p_4var;
  float gWin_44p_4var;
  float xWin_44p_6var;
  float yWin_44p_6var;
  float zWin_44p_6var;
  float gWin_44p_6var;
  float xWin_44h_4var;
  float yWin_44h_4var;
  float zWin_44h_4var;
  float gWin_44h_4var;
  float xWin_44h_6var;
  float yWin_44h_6var;
  float zWin_44h_6var;
  float gWin_44h_6var;
  float xWin_64p_4var;
  float yWin_64p_4var;
  float zWin_64p_4var;
  float gWin_64p_4var;
  float xWin_64p_6var;
  float yWin_64p_6var;
  float zWin_64p_6var;
  float gWin_64p_6var;
  float xWin_64h_4var;
  float yWin_64h_4var;
  float zWin_64h_4var;
  float gWin_64h_4var;
  float xWin_64h_6var;
  float yWin_64h_6var;
  float zWin_64h_6var;
  float gWin_64h_6var;
  //////////////////////
  float trig_SF;
  float trig_SF2;
  float trig_SFC;
  int   trig_flag;
  //////////////////////
  int   total_jet_ntrk1;
  float total_jet_width1;
  int   total_jet_ntrk2;
  float total_jet_width2;
  int   total_jet_ntrk3;
  float total_jet_width3;
  //Flavour Composition Variables
  float AllJet_MV1_1;    // 1 good jet 
  float AllJet_MV1_2;    // 2 good jet
  float AllJet_MV1_3;    // 3 good jet
} analysis_output_struct;


typedef struct SOMVar{
  float Ntrk, width, JetMass;
} SOMVar;


class HiggsllqqAnalysis : public HiggsAnalysis {
 public:
  
  enum {
    MU2,
    MUE,
    E2,
  };
  
  
  // constructor, destructor
  HiggsllqqAnalysis(TTree * /*tree*/ = 0)
    {
    }
  ~HiggsllqqAnalysis();
  
  
  // options setters
  virtual void setAnalysisVersion(TString val)
  {
    m_analysis_version = val;
  }
  
  virtual void setSmearing(Bool_t val)
  {
    m_doSmearing = val;
  }
  
  virtual void setTopoIso(Bool_t val)
  {
    m_useTopoIso = val;
  }
  
  virtual void setMuonFamily(Int_t val)
  {
    m_muonFamily = val;
  }
  
  virtual void setElectronFamily(Int_t val)
  {
    m_electronFamily = val;
  }
  
  virtual void setJetFamily(Int_t val)
  {
    m_jetFamily = val;
  }

  virtual void setJetbTagger(Int_t val)
  {
    m_JetbTagger = val;
  }
  
  virtual void setOutputFile(TString val)
  {
    m_outputFileName = val;
  }
  
  // Method to setup the Systematic Jet Studies (JES-JER)
  virtual void SetSysStudy(Bool_t val) 
  {
    m_sysstudy = val;
  }
  
  
 protected:
  
  // pointers to ATLAS tools
  Root::TGoodRunsListReader          *m_GRL;
  Root::TPileupReweighting           *m_PileupReweighter;
  D3PD::TrigDecisionToolD3PD         *m_TrigDecisionToolD3PD;
  egRescaler::EnergyRescalerUpgrade  *m_ElectronEnergyRescaler;
  
  Analysis::AnalysisMuonConfigurableScaleFactors  *m_MuonEffSF;
  Analysis::AnalysisMuonConfigurableScaleFactors  *m_MuonEffSFCalo;
  Analysis::AnalysisMuonConfigurableScaleFactors  *m_MuonEffSFSA;
  
  MuonSmear::SmearingClass *m_MuonSmearer;
  
  egammaSFclass                  *m_ElectronEffSF;
  VertexPositionReweightingTool  *m_VertexPositionReweighter;
  LeptonTriggerSF                *m_MuonTrigSF;
  ggFReweighting                 *m_ggFReweighter;
  TriggerNavigationVariables      m_trigNavVar;
  MuonTriggerMatching            *m_MuonTriggerMatchTool;
  ElectronTriggerMatching        *m_ElectronTriggerMatchTool;
  TH2F      *m_smearD0[3];
  TRandom3   m_smearD0_rand;
  TAxis     *m_smearD0_x;
  TAxis     *m_smearD0_y;
  
  
  // containers for leptons, Jets, dileptons, in the event
  std::vector<Analysis::ChargedLepton *>  m_Muons;
  std::vector<Analysis::ChargedLepton *>  m_Electrons;
  std::vector<Analysis::Jet *>            m_Jets;
  std::vector<Analysis::ChargedLepton *>  m_GoodMuons;
  std::vector<Analysis::ChargedLepton *>  m_GoodElectrons;
  std::vector<Analysis::Jet *>            m_GoodJets;
  std::vector<Analysis::Dilepton *>       m_Dileptons;
  
  
  // utility maps
  std::map<UInt_t, Float_t>  m_CrossSection;
  std::map<UInt_t, Int_t>    m_SignalSampleMass;
  
  
  // channels to be processed and cutflows
  std::vector<Int_t> m_Channels;
  std::vector<Analysis::CutFlowTool> m_EventCutflow;
  std::vector<Analysis::CutFlowTool> m_EventCutflow_rw;
  std::vector<Analysis::CutFlowTool> m_EventCutflow0tag;
  std::vector<Analysis::CutFlowTool> m_EventCutflow0tag_rw;
  std::vector<Analysis::CutFlowTool> m_EventCutflow1tag;
  std::vector<Analysis::CutFlowTool> m_EventCutflow1tag_rw;
  std::vector<Analysis::CutFlowTool> m_EventCutflow2tag;
  std::vector<Analysis::CutFlowTool> m_EventCutflow2tag_rw;
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
  virtual Int_t   getLastCutPassed();
  virtual Int_t   getNumberOfGoodVertices();
  virtual Bool_t  SherpaPt0Veto();
  virtual Bool_t  passesGRL();
  virtual Bool_t  hasGoodVertex();
  virtual Bool_t  passesTrigger();
  virtual Bool_t  passesSingleMuonTrigger();
  virtual Bool_t  passesDiMuonTrigger();
  virtual Bool_t  passesSingleElectronTrigger();
  virtual Bool_t  passesDiElectronTrigger();
  virtual Bool_t  passesElectronMuonTrigger();
  virtual TString getSingleMuonTriggerName();
  virtual TString getDiMuonTriggerName();
  virtual TString getSingleElectronTriggerName();
  virtual TString getDiElectronTriggerName();
  virtual TString getElectronMuonTriggerName();
  virtual void    applyChanges( Analysis::ChargedLepton *lep );
  virtual void    applyChanges( Analysis::Jet           *jet );
  virtual void    getMuons( D3PDReader::MuonD3PDObject     *mu_branch, Int_t family );
  virtual void    getElectrons( D3PDReader::ElectronD3PDObject *el_branch, Int_t family );
  virtual void    getJets( D3PDReader::JetD3PDObject       *jet_branch );
  virtual void    getGoodMuons();
  virtual void    getGoodElectrons();
  virtual void    getGoodJets();
  virtual void    getGoodLeptons();
  virtual Float_t getDiLeptonMass();
  
  virtual std::vector<TString> getListOfAlternativeTriggers(TString sequence);
  
  
  // object selection
  virtual Bool_t isMCPMuon(Analysis::ChargedLepton *lep);
  virtual Bool_t isGood(Analysis::ChargedLepton *lep);
  Bool_t isGoodJet(Analysis::Jet *jet);
  
  //object quality
  Bool_t Pair_Quality();
  
  // utility functions for the selection
  virtual Bool_t isMC()
  {
    return ntuple->mc.n.IsAvailable();
  }
  
  virtual TString analysis_version()
  {
    return m_analysis_version;
  }
  
  virtual Bool_t doSmearing()
  {
    return m_doSmearing;
  }
  
  virtual Bool_t useTopoIso()
  {
    return m_useTopoIso;
  }
  
  virtual Int_t getMuonFamily()
  {
    return m_muonFamily;
  }
  
  virtual Int_t getElectronFamily()
  {
    return m_electronFamily;
  }
  
  virtual Int_t getJetFamily()
  {
    return m_jetFamily;
  }
 

  virtual Int_t getJetbTagger()
  {
    return m_JetbTagger;
  }
 
  virtual Int_t getTriggerInfo(TString chain);
  virtual Int_t getPeriod();
  
  virtual void setChannel(Int_t chan)
  {
    m_thisChannel = chan;
  }
  
  virtual Int_t getChannel()
  {
    return m_thisChannel;
  }
  
  
  // event weighting helpers
  virtual Float_t getEventWeight();
  virtual Float_t getPileupWeight();
  virtual Float_t getVertexZWeight();
  virtual Float_t getLeptonWeight(Analysis::ChargedLepton *lep);
  virtual Float_t getSFWeight();
  virtual Float_t getggFWeight();
  virtual Float_t getDPhijjZWeight();
  
  
  // Trigger SF 
  Float_t getCandidateTriggerSF(TString syst = "");
  
  // tools for cross section computation
  void            initCrossSections();  // to be change!! error
  Float_t         getCrossSectionWeight();
  Float_t         getTruthHiggsMass();
  virtual Float_t getTruthHiggsPt();    // to be change!! error
  
  void FillHllqqCutFlowXtag(int last_event,UInt_t chan);
  
  // option resume
  virtual void printAllOptions()
  {
    Info("printAllOptions", "======== options =======");
    Info("printAllOptions", "analysis_version = %s", m_analysis_version.Data());
    Info("printAllOptions", "doSmearing       = %s", m_doSmearing ? "kTRUE" : "kFALSE");
    Info("printAllOptions", "useTopoIso       = %s", m_useTopoIso ? "kTRUE" : "kFALSE");
    Info("printAllOptions", "muonFamily       = %d", m_muonFamily);
    Info("printAllOptions", "electronFamily   = %d", m_electronFamily);
    Info("printAllOptions", "jetFamily        = %d", m_jetFamily);
    Info("printAllOptions", "jet b Tagger     = %d", m_JetbTagger);
    Info("printAllOptions", "outputFileName   = %s", m_outputFileName.Data());
    Info("printAllOptions", "========================");
  }
  
  
  // Methods from the qqll analysis
  Bool_t  ApplyChangesMuon(Analysis::ChargedLepton *lep);
  Bool_t  ApplyChangesElectron(Analysis::ChargedLepton *lep);
  Bool_t  ApplyChangesJet(Analysis::Jet *jet);
  Bool_t  JetInHole();
  Bool_t  NotMETclean();
  Bool_t  hasGoodMET();
  Bool_t  GetDoLowMass() { return m_dolowmass; }
  Bool_t  GetSysStudy()  { return m_sysstudy;  }
  Bool_t  IsConsistentWithTrigger();
  
  void    LoadGRL();
  void    InitReducedNtuple();
  void    ResetReducedNtupleMembers();
  void    FillReducedNtuple(Int_t cut, UInt_t channel);
  
  Float_t getCorrectMETValue();
  Float_t GetMV1value(Analysis::Jet *jet);
  pair <Int_t,Double_t> GetFlavour(Analysis::Jet *jet);
  
  // Flavor Jet Methods
  Bool_t isHeavyJet(Int_t pdg);
  Bool_t isLightJet(Int_t pdg);
  Bool_t isGluonJet(Int_t pdg);
  
  // Higgs MC (weight) Methods
  void    InitMasses();
  Float_t GetTruthHiggsPt();
  Float_t GetggFWeight();
  
  
  // Tracks and Widths 2012 methods
  Float_t isInTheJet(Int_t Index, Int_t JetIndex, vector<float> *whatinjet_eta, vector<float> *whatinjet_phi);
  std::pair<Float_t, Float_t>InfoTracks(Int_t JetIndex);
  Bool_t  isGoodTrack(Int_t TrackIndex);
  
  
  // Method to count the number of that events into the GoodJets vector.
  Int_t GetNumOfTags();
  
  
  // 2 Methods to evaluate the impact of the Jet kinematic cuts in the selection
  Bool_t GoodPtJets();
  Bool_t GoodEtaJets();
  
  
  // Method to evaluate the isolation mu-jet
  Bool_t MuonJetOR();
  
  
  // Method to find the best DiJets using JetKinematicFitter
  Bool_t JetKinematicFitterResult();
  
  // Method to find the best DiJets using BestPair algo
  Bool_t JetBestPairResult();
  
  // Method to calculate the DiJet invariant mass for the tagged Jets!
  Bool_t JetDimassTagged();
  
  Bool_t JetDimassOneTagged();
  
  
  // Method to flag/control the QCD selection  
  Bool_t GetDoQCDSelection()
  { 
    if(!isMC()) return m_doqcdselection; 
    return kFALSE;
  }  
  
  void SetDoQCDSelection(Bool_t val) { m_doqcdselection = val; }
  
  
  // Method to flag/control the Low or High Mass selection    
  void SetDoLowMass(Bool_t val) { m_dolowmass = val; }
  
  
  // Methods to sort the Jets or Leptons by pt
  void SortIndex(std::vector<Analysis::Jet *> &vec);
  void SortIndex(std::vector<Analysis::ChargedLepton *> &vec);
  
  
  // Printing the Event tables
  void PrintCutFlowEvents(Int_t cut, Int_t event, Int_t count);
  
  // Method to calculate the Cleaning of the jet
  Bool_t isBadLooser(Analysis::Jet *jet);
  
  // Method to fill vectors to calculate the HFOR value
  void FillHFORvariables();
  
  
  // Methods to Set, Reset and Fill the TestSelection Struct
  void SetAnalysisOutputBranches(analysis_output_struct *str);
  void ResetAnalysisOutputBranches(analysis_output_struct *str);
  void FillAnalysisOutputTree(analysis_output_struct *str, Int_t cut, UInt_t channel);
  pair <double,double> GetJetSFsvalue(int jetindex);
  
  // MVA methods (June 2013)
  void    SetTmvaReaders(TMVA::Reader *reader[36],Float_t var1[36], Float_t var2[36]);
  Float_t getFisher_KF(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getLikelihood_KF(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getLikelihoodMIX_KF(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getFisher_BP(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getLikelihood_BP(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getLikelihoodMIX_BP(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getFisher_LJ(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getLikelihood_LJ(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  Float_t getLikelihoodMIX_LJ(TMVA::Reader *reader[4],Float_t var1[4], Float_t var2[4], Float_t pt_jet, Float_t ntrk_jet,Float_t width_jet);
  
  
  // SOM METHODS (Dec 2012)
  int GetSOMx(TString which_map, int jetidx1, int jetidx2) { return (GetSOMWinner(which_map, jetidx1, jetidx2)).first;}
  
  int GetSOMy(TString which_map, int jetidx1, int jetidx2) { return (GetSOMWinner(which_map, jetidx1, jetidx2)).second;}
  
  float GetSOMz(TString which_map, int jetidx1, int jetidx2) {return GetSOMWinnerDist(which_map, jetidx1, jetidx2);}
  
  int GetSOMg(TString which_map, int jetidx1, int jetidx2);
  
  pair <int, int> GetSOMWinner(TString which_map, int jetidx1, int jetidx2);
  
  float GetSOMWinnerDist(TString which_map, int jetidx1, int jetidx2);
  
  pair <SOMVar, SOMVar> GetSOMVariables(int jetidx1, int jetidx2);
  
  pair<int, float> GetWinnerInfo(std::vector< pair<SOMVar,SOMVar> > SOMVectors, pair <SOMVar,SOMVar> Event);
  
  std::vector< pair<SOMVar,SOMVar> > GetSOMVectors(TString which_map, int jetidx1, int jetidx2);
  
  TString GetCorrectMap (TString which_map, int jetidx1, int jetidx2);
  
  std::vector< pair <int,int> > GetCoordinates(TString which_map);
  
  Bool_t CheckMap(TString which_map, int jetidx1, int jetidx2);
  
  Float_t Rightcut(Int_t efficiency, Float_t pt_jet, Float_t eta_jet);
  
  
 protected:
  Bool_t m_called_getGoodLeptons;
  Bool_t m_called_getGoodJets;
  Bool_t m_called_getGoodObjects;
  Long64_t m_entriesInChain;
  Long64_t m_processedEntries;
  
 private: 
  TString m_analysis_version;
  Bool_t  m_doSmearing;
  Bool_t  m_useTopoIso;
  Int_t   m_muonFamily;
  Int_t   m_electronFamily;
  Int_t   m_jetFamily;
  Int_t   m_JetbTagger;
  Int_t   m_thisChannel;

  TH1D *m_generatedEntriesHisto;
  TH1D *h_cutflow;
  TH1D *h_cutflow_weight;
  TH1D *h_cutflow_E2;
  TH1D *h_cutflow_weight_E2;
  TH1D *h_cutflow_MU2;
  TH1D *h_cutflow_weight_MU2;
  TH1D *h_cutflow_MUE;
  TH1D *h_cutflow_weight_MUE;
  
  TMVA::Reader *reader[36];
  Float_t var1[36],var2[36];
  
  TTree  *m_TreeCutflow;
  TString m_outputFileName;
  TFile  *m_outputFile;
  
  Analysis::CalibrationDataInterfaceROOT *calib;
  Analysis::CalibrationDataVariables ajet;
  Analysis::Uncertainty uncertainty;
  
  // Definition of the TestSelection Struct
  analysis_output_struct m_outevent;
  
  // Beginning
 protected:
  
  // JES AND JER TOOL
  JetCalibrationTool *myJES;
  JetSmearingTool    *myJER;
  
  // ggF reweighting tool
  ggFReweighting *fggFReweighter;
  
  // Jet kinematic fitter
  JetKinematicFitter *m_jetkinematicfitter;
  
  // HFOR TOOL
  HforToolD3PD *hforTool;
  
  // llqq Analysis tree
  TTree *m_reduced_ntuple;
  
  // TestSelection Tree
  TTree *analysistree;
  
  Int_t m_lep_chargeproduct;
  std::vector<Float_t>   *m_lep_m;
  std::vector<Float_t>   *m_lep_pt;
  std::vector<Float_t>   *m_lep_eta;
  std::vector<Float_t>   *m_lep_phi;
  std::vector<Int_t>     *m_lep_charge;
  std::vector<Float_t>   *m_lep_d0;
  std::vector<Float_t>   *m_lep_sigd0;
  std::vector<Float_t>   *m_lep_trackiso;
  std::vector<Float_t>   *m_lep_caloiso;
  std::vector<Int_t>     *m_lep_quality;
  std::vector<Float_t>   *m_jets_m;
  std::vector<Float_t>   *m_jets_pt;
  std::vector<Float_t>   *m_jets_eta;
  std::vector<Float_t>   *m_jets_eta_det;
  std::vector<Float_t>   *m_jets_phi;
  std::vector<Float_t>   *m_jets_MV1;
  std::vector<Float_t>   *m_jets_flavortruth;
  std::vector<Float_t>   *m_jets_jvtxf;
  std::vector<Int_t>     *m_jets_nTrk;
  std::vector<Float_t>   *m_jets_width;
  std::vector<Int_t>     *m_jets_flavorpdg;
  std::vector<Double_t>  *m_jets_Epdg;
  //
  UInt_t  m_run;
  UInt_t  m_event;
  Int_t   m_cut;
  Int_t   m_channel;
  Int_t   m_qcdevent;
  Int_t   m_NPV;
  Int_t   m_Entries;
  Int_t   m_HFOR;
  Int_t   m_low_event;
  Int_t   m_trig_flag;
  //
  Float_t m_weight;
  Float_t m_SFWeight;
  Float_t m_ggFweight;
  Float_t m_EventWeight;
  Float_t m_PileupWeight;
  Float_t m_VertexZWeight;
  Float_t m_DPhijjZWeight;
  Float_t m_TriggerSFWeight;
  Float_t m_mu;
  Float_t m_truthH_pt;
  Float_t m_met_met;
  Float_t m_met_sumet;
  Float_t m_met_phi;
  Float_t m_trig_SF;
  Float_t m_trig_SF2;
  Float_t m_trig_SFC;
  
  // Truth quark information
  std::vector<Float_t> *m_quark_m;
  std::vector<Float_t> *m_quark_pt;
  std::vector<Int_t>   *m_quark_pdg;
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
