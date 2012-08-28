#define HiggsllqqAnalysis_cxx

#include "HiggsllqqAnalysis/HiggsllqqAnalysis.h"

#include <fstream>
#include <TStyle.h>

//#define NEWTRIGSF


/*  when introducing changes depending on
    - 2011 vs 2012
    - 7 vs 8 TeV
    - rel. 17 vs rel. 17.2
    - mc11c vs mc12a
    please always check if just assuming "analysis_version()" (rel_17 or rel_17_2) is still reasonable, or one needs to handle with care these three differences

    !!! WARNING: it is expected to happen *AS SOON AS 2011 DATA IS REPROCESSED* !!!

    throughout the code, RELEASE17.2 means that the stuff there must or might be changed if running with all or particular 17.2 samples
    these changes are, unlike those relying on analysis_version(), to be done manually editing this code


  Information: dolowmass for muons     = True   if GetDoLowMass() == True
  dolowmass for electrons = True   if GetDoLowMass() == True
  
  August 10th 2012
*/

Bool_t cut_leptons   = kFALSE, 
  NotTrigConsistency = kFALSE, 
  NotPtConsistency   = kTRUE, 
  
//Smearing Options:
  MuonSmearing       = kFALSE, 
  JetSmearing        = kFALSE, 
  ElectronSmearing   = kFALSE,
  
  TaggedChannel      = kFALSE, 
  DoMETdataClean     = kFALSE, 
  is4lepGood         = kFALSE, 
  METtype_RefFinal   = kTRUE, 
  DoLowMass          = kTRUE, 
  DollqqAnalysis     = kTRUE, 
  DoCaloMuons        = kFALSE;

int count_events(0),eventNow(-1),overElectron(0),overMuon(0),overJet(0); int badevent=0, prebadevent = 0, ptchange=0, ptelecChange=0;
int periodBD(0),periodEH(0),periodI(0),periodJK(0),periodLM(0);
Float_t Muon0(0),Muon1(0),Muon2(0),Muon3(0),Muon4(0),Muon5(0),Muon6(0),Muon7(0),Muon8(0);
Float_t Electron0(0),Electron1(0),Electron2(0),Electron3(0),Electron4(0),Electron5(0),Electron6(0), b_rescaling = 1.05/*Fixed*/; 
int HFOR_value = -999;


HiggsllqqAnalysis::~HiggsllqqAnalysis()
{
}


Bool_t HiggsllqqAnalysis::ApplyChangesMuon(Analysis::ChargedLepton *lep)
{
  if (lep->flavor() == Analysis::ChargedLepton::MUON) {
    D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
    
    if (isMC() && doSmearing() && MuonSmearing) {
      double eta  = lep->Get4Momentum()->Eta();
      double ptcb = lep->Get4Momentum()->Pt();
      double ptme = lep->Get4Momentum_SA()->Pt();
      double ptid = lep->Get4Momentum_ID()->Pt();
      double mu_charge = lep->charge();
      double chargeflip = 1;
      
      m_MuonSmearer->SetSeed(ntuple->eventinfo.EventNumber(), (Int_t)mu->GetIndex());     
      
      if (mu->isCombinedMuon()) {
	m_MuonSmearer->Event(ptme, ptid, ptcb, eta, mu_charge);
	lep->Get4Momentum()->SetPtEtaPhiM(m_MuonSmearer->pTCB(), lep->Get4Momentum()->Eta(), lep->Get4Momentum()->Phi(), lep->Get4Momentum()->M());
	lep->Get4Momentum_SA()->SetPtEtaPhiM(m_MuonSmearer->pTMS(), lep->Get4Momentum_SA()->Eta(), lep->Get4Momentum_SA()->Phi(), lep->Get4Momentum_SA()->M());
	lep->Get4Momentum_ID()->SetPtEtaPhiM(m_MuonSmearer->pTID(), lep->Get4Momentum_ID()->Eta(), lep->Get4Momentum_ID()->Phi(), lep->Get4Momentum_ID()->M());
	chargeflip = (m_MuonSmearer->ChargeFlipCB());
      } else if (mu->isStandAloneMuon()) {
	m_MuonSmearer->Event(ptcb, eta, "MS", mu_charge);
	lep->Get4Momentum()->SetPtEtaPhiM(m_MuonSmearer->pTMS(), lep->Get4Momentum()->Eta(), lep->Get4Momentum()->Phi(), lep->Get4Momentum()->M());
	lep->Get4Momentum_SA()->SetPtEtaPhiM(m_MuonSmearer->pTMS(), lep->Get4Momentum_SA()->Eta(), lep->Get4Momentum_SA()->Phi(), lep->Get4Momentum_SA()->M());
	chargeflip = m_MuonSmearer->ChargeFlipMS();
      } else { // segment tagged, calo
	m_MuonSmearer->Event(ptcb, eta, "ID", mu_charge);
	lep->Get4Momentum()->SetPtEtaPhiM(m_MuonSmearer->pTID(), lep->Get4Momentum()->Eta(), lep->Get4Momentum()->Phi(), lep->Get4Momentum()->M());
	lep->Get4Momentum_ID()->SetPtEtaPhiM(m_MuonSmearer->pTID(), lep->Get4Momentum_ID()->Eta(), lep->Get4Momentum_ID()->Phi(), lep->Get4Momentum_ID()->M());
	chargeflip = m_MuonSmearer->ChargeFlipID();
      }
      
      if(chargeflip<0)
	{
	  cout<<"Charge Changed: "<<chargeflip<<endl;
	  lep->set_charge(mu->charge()*chargeflip);
	}
      
      if(TMath::Abs((ptcb - mu->pt())/ptcb)>=0.10) ptchange++;
    }
    return kTRUE;
  }
  else{
    cout<<" This is not a Muon! "<<endl;
    return kFALSE;
  }
}

Bool_t HiggsllqqAnalysis::ApplyChangesElectron(Analysis::ChargedLepton *lep)
{
  Int_t SYST_FLAG = 0; //SYST_FLAG is 0 for nominal scale, 1 or 2 for 1-sigma variations.
  
  if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {  
    
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    // first of all one must rescale energy in the crack! (both data and MC but only for 2011 so far)
    Float_t tmp_calibration(1.);
    if (analysis_version() == "rel_17") // rel. 17
      tmp_calibration = m_ElectronEnergyRescaler->applyMCCalibrationMeV(el->cl_eta(), el->cl_E() / TMath::CosH(el->tracketa()), "ELECTRON");
    Float_t tmp_E = el->cl_E() * TMath::Abs(tmp_calibration);
    Float_t tmp_Et = tmp_E / TMath::CosH(el->tracketa());
    Float_t tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    
    // then, apply the other corrections
    if (isMC() && doSmearing() && ElectronSmearing) {
      m_ElectronEnergyRescaler->SetRandomSeed(ntuple->eventinfo.EventNumber() + 100 * (Int_t)el->GetIndex());
      
      // false here means the MC is mc11c (no constant term)
      Float_t smearcorr = m_ElectronEnergyRescaler->getSmearingCorrectionMeV(el->cl_eta(), tmp_E, SYST_FLAG, false, "2011");
      
      tmp_E = tmp_E * smearcorr;
      tmp_Et = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    if (!isMC() && ElectronSmearing) {
      tmp_E = m_ElectronEnergyRescaler->applyEnergyCorrectionMeV(el->cl_eta(), el->cl_phi(), tmp_E, tmp_Et, SYST_FLAG, "ELECTRON");
      tmp_Et = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    
    Float_t tmp_pt = (TMath::Sqrt(TMath::Power(tmp_E, 2) - TMath::Power(lep->Get4Momentum()->M(), 2))) * TMath::Sin(2 * TMath::ATan(TMath::Exp(-el->tracketa())));
    Float_t tmp_pt_cl = (TMath::Sqrt(TMath::Power(tmp_E, 2) - TMath::Power(lep->Get4Momentum()->M(), 2))) * TMath::Sin(2 * TMath::ATan(TMath::Exp(-el->cl_eta())));
    
    // correct lepton 4-momenta
    lep->Get4Momentum()->SetPtEtaPhiM(tmp_pt, el->tracketa(), el->trackphi(), lep->Get4Momentum()->M());
    lep->Get4Momentum_ID()->SetPtEtaPhiM(el->trackpt(), el->tracketa(), el->trackphi(), lep->Get4Momentum()->M()); // untouched by energy corrections!
    lep->Get4Momentum_SA()->SetPtEtaPhiM(tmp_pt_cl, el->cl_eta(), el->cl_phi(), lep->Get4Momentum()->M());
    return kTRUE;
  }
  else{
    cout<<" This is not a Electron! "<<endl;
    return kFALSE;
  }
}

Bool_t HiggsllqqAnalysis::ApplyChangesJet(Analysis::Jet *jet)
{
  Float_t tmp_E(-9999.9), tmp_pt(-9999.9), tmp_eta(-9999.9), tmp_phi(-9999.9), tmp_Et(-9999.9);
  Bool_t sysstudy = GetSysStudy(); 
  
  D3PDReader::JetD3PDObjectElement *Jet = jet->GetJet();
  
  if (isMC() && doSmearing() && JetSmearing) { //NOT Jet Smearing Option! 14th June 2012
    double Eraw    = Jet->emscale_E();
    double eta_det = Jet->emscale_eta();
    double eta     = Jet->EtaOrigin();
    double phi     = Jet->PhiOrigin();
    double m       = Jet->MOrigin();
    
    double mu = ntuple->eventinfo.averageIntPerXing();
    int NPV=0;
    for (Int_t i = 0; i < ntuple->vxp.n(); i++) {
      if (ntuple->vxp[i].trk_n() >= 3) NPV++;
    }
    
    // Calibrate the jet!
    // Pile-up, origin, EtaJES correction applied, i.e. to OFFSET_ORIGIN_ETAJES scale
    TLorentzVector jet4v = myJES->ApplyOffsetEtaJES(Eraw,eta_det,eta,phi,m,mu,NPV);  
    
    // The below is systematic evaluation, and ONLY for MC
    // Smear the jet to match the MC resolution+1 sigma!
    if(isMC() && sysstudy){
      myJER->SetSeed(ntuple->eventinfo.EventNumber());
      myJER->SmearJet_SystRel17(jet4v);
    }
    
    tmp_E=jet4v.E();
    tmp_pt=jet4v.Pt();
    tmp_eta=jet4v.Eta();
    tmp_phi=jet4v.Phi();
    tmp_Et=jet4v.Et();
    
    //Inside the "if"
    jet->set_rightE(tmp_E);
    jet->set_rightpt(tmp_pt);
    jet->set_righteta(tmp_eta);
    jet->set_rightphi(tmp_phi);
    jet->set_rightEt(tmp_Et);
  }  
  
  return kTRUE;  
}

void HiggsllqqAnalysis::LoadGRL()
{
  // set GRL according to channel
  
  if (getChannel() == HiggsllqqAnalysis::MU2) {
    m_GRL->SetXMLFile("/afs/cern.ch/user/a/arturos/LoSterzo_llqq/HiggsllqqAnalysis/ANALYSIS/HiggsqqllAnalysis/grl/data_GRL_Hllqq.xml");
  } else if (getChannel() == HiggsllqqAnalysis::MU2) {
    m_GRL->SetXMLFile("/afs/cern.ch/user/a/arturos/LoSterzo_llqq/HiggsllqqAnalysis/ANALYSIS/HiggsqqllAnalysis/grl/data_GRL_Hllqq.xml");
  } else if (getChannel() == HiggsllqqAnalysis::MU2) {
    m_GRL->SetXMLFile("/afs/cern.ch/user/a/arturos/LoSterzo_llqq/HiggsllqqAnalysis/ANALYSIS/HiggsqqllAnalysis/grl/data_GRL_Hllqq.xml");
  } else {
    Error("doAnalysis", "unsupported channel, aborting.");
    exit(1);
  }
}

Bool_t HiggsllqqAnalysis::PassesGRL()
{
  if (isMC()) {
    return kTRUE;
  } else {
    if (ntuple->eventinfo.RunNumber() > 199999/*191517*/) { // recent runs for which no GRL is available
      return kTRUE;
    } else {
      return m_GRL->GetGoodRunsList(getChannel() + HllqqGRL::data11).HasRunLumiBlock(ntuple->eventinfo.RunNumber(), ntuple->eventinfo.lbn());
    }
  } 
}

Int_t HiggsllqqAnalysis::GetLastCutPassed()
{
  //reeplaced by getLastCutPastsed()
  return 0;
}



Bool_t HiggsllqqAnalysis::change_input()
{
  Info("change_input", "Processing a new file...");
  
  if (!fChain) {
    Error("change_input", "empty fChain pointer!");
    
    return kTRUE;
  }
  
  // initialize the trigger decision tool
  
  if (!m_TrigDecisionToolD3PD->SetEventTree(fEventTree)) {
    Error("change_input", "Problems with setting the event tree to the TDT");
  }
  if (!m_TrigDecisionToolD3PD->GetConfigSvc(kFALSE).SetConfigTree(fConfigTree)) {
    Error("change_input", "Problems with setting the config tree to the TDT");
  }
  
  // set TTree cache wisely
  if (m_processedEntries > 500) {
    ntuple->mu_staco.GetStatistics().AddToTreeCacheByEntries(fEventTree, 200);
    ntuple->mu_muid.GetStatistics().AddToTreeCacheByEntries(fEventTree, 200);
    ntuple->vxp.GetStatistics().AddToTreeCacheByEntries(fEventTree, 200);
    ntuple->mc.GetStatistics().AddToTreeCacheByEntries(fEventTree, 200);
  }
  
  // store number of entries to be processed
  m_entriesInChain = fChain->GetEntries();
  
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::initialize_tools()
{
  //Set the Low or High Mass Analysis
  SetDoLowMass(DoLowMass);
  
  // initiate the calibration tool
  TString jetAlgo="AntiKt4TopoEM";
  TString JES_config_file="ApplyJetCalibration/data/CalibrationConfigs/Rel17_JES.config";
  myJES = new JetCalibrationTool(jetAlgo,JES_config_file);
  myJER = new JetSmearingTool(jetAlgo);
  
  // initialize the kinematic fitter
  Info("doAnalysis", "Initializing JetKinematicFitter");
  int maxjet_KF = 5;
  m_jetkinematicfitter = new JetKinematicFitter(maxjet_KF,91187.6,2495.2); // mZ | gammaZ
  m_jetkinematicfitter->SetIsMC(isMC());
  
  // HFOR D3PD TOOL
  hforTool = new HforToolD3PD();
  
  printAllOptions();
  
  // trigger decision tool
  m_TrigDecisionToolD3PD = new D3PD::TrigDecisionToolD3PD();
  
  // goodrunslists
  m_GRL = new Root::TGoodRunsListReader();
  // put 2011 first
  m_GRL->AddXMLFile("HiggsllqqAnalysis/grl/data11_GRL_MU4.xml");////*//// Change with the real GRL for llqq analysis!
  m_GRL->AddXMLFile("HiggsllqqAnalysis/grl/data11_GRL_MU2E2.xml");
  m_GRL->AddXMLFile("HiggsllqqAnalysis/grl/data11_GRL_E4.xml");
  // then add 2012
  m_GRL->AddXMLFile("HiggsllqqAnalysis/grl/data12_GRL_MU4.xml");
  m_GRL->AddXMLFile("HiggsllqqAnalysis/grl/data12_GRL_MU2E2.xml");
  m_GRL->AddXMLFile("HiggsllqqAnalysis/grl/data12_GRL_E4.xml");
  
  m_GRL->Interpret();
  
  // pileup
  Info("initialize_tools", "Initializing the pile-up reweighting tool...");
  
  m_PileupReweighter = new Root::TPileupReweighting("m_PileupReweighter");
  m_PileupReweighter->SetUnrepresentedDataAction(2);
  
  if (analysis_version() == "rel_17") { // mc11c, 2011
    m_PileupReweighter->AddConfigFile("./HiggsllqqAnalysis/packages/files/pileup/MC11c.prw.root");
    m_PileupReweighter->AddLumiCalcFile("./HiggsllqqAnalysis/packages/files/pileup/ilumicalc_period_AllYear_Higgs_4l_2e2mu.root");
    m_PileupReweighter->SetDefaultChannel(109292);
  } else if (analysis_version() == "rel_17_2") { // mc12a, 2012
    m_PileupReweighter->AddConfigFile("./HiggsllqqAnalysis/packages/files/pileup/MC12a.prw.root");
    m_PileupReweighter->AddLumiCalcFile("./HiggsllqqAnalysis/packages/files/pileup/ilumicalc_2012_period_AllYear_Higgs_4l_2e2mu.root");
    m_PileupReweighter->SetDefaultChannel(160156);
  }
  m_PileupReweighter->Initialize();
  
  // vertex reweighting (MC12 only)
  Info("initialize_tools", "Initializing the vertex reweighting tool...");
  if (analysis_version() == "rel_17_2") {
    m_VertexPositionReweighter = new VertexPositionReweightingTool(VertexPositionReweightingTool::MC12a);
  } else {
    m_VertexPositionReweighter = 0;
  }
  
  // ElectronEnergyRescaler
  Info("initialize_tools", "Initializing the EnergyRescaler tool...");
  m_ElectronEnergyRescaler = new eg2011::EnergyRescaler();
  if (analysis_version() == "rel_17")
    m_ElectronEnergyRescaler->useDefaultCalibConstants("2011");
  else if (analysis_version() == "rel_17_2")
    m_ElectronEnergyRescaler->useDefaultCalibConstants("2012");
  
  
  // ElectronEffSF
  Info("initialize_tools", "Initializing the egammaSFclass tool...");
  m_ElectronEffSF = new egammaSFclass();
  
  // MuonEffSF
  Info("initialize_tools", "Initializing the Analysis::AnalysisMuonConfigurableScaleFactors tool...");
  
  Analysis::AnalysisMuonConfigurableScaleFactors::Configuration muonEffConfig;
  Analysis::AnalysisMuonConfigurableScaleFactors::Configuration muonEffConfigSA;
  
  if (analysis_version() == "rel_17") {
    muonEffConfig   = Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverPeriods;
    muonEffConfigSA = Analysis::AnalysisMuonConfigurableScaleFactors::Default;
  } else if (analysis_version() == "rel_17_2") {
    muonEffConfig   = Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverRuns;
    muonEffConfigSA = Analysis::AnalysisMuonConfigurableScaleFactors::AverageOverRuns;
  }
  
  std::vector<Double_t> int_lum = m_PileupReweighter->getIntegratedLumiVector();
  
  std::string mu_alg;
  std::string mu_alg_calo;
  std::string mu_alg_sa;
  
  if (analysis_version() == "rel_17") {
    if (getMuonFamily() == Muon::STACO) {
      mu_alg = "STACO_CB_plus_ST_2011_SF.txt";
      mu_alg_sa = "STACOHighEta.txt";
    } else if (getMuonFamily() == Muon::MUID) {
      mu_alg = "Muid_CB_plus_ST_2011_SF.txt";
      mu_alg_sa = "MuidHighEta.txt";
    }
    mu_alg_calo = "CaloTag_2011_SF.txt";
  } else if (analysis_version() == "rel_17_2") {
    if (getMuonFamily() == Muon::STACO) {
      mu_alg =  "STACO_CB_plus_ST_2012_SF.txt";
      mu_alg_sa = "STACO_CB_plus_ST_2012_SFms.txt";
    } else if (getMuonFamily() == Muon::MUID) {
      mu_alg = "Muid_CB_plus_ST_2012_SF.txt";
      mu_alg_sa = "Muid_CB_plus_ST_2012_SFms.txt";
    }
    mu_alg_calo =  "CaloTag_2012_SF.txt";
  }
  
  std::string unit("MeV");
  std::string muon_sf_directory("MuonEfficiencyCorrections/share/");
  
  m_MuonEffSF = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_directory, mu_alg, unit, muonEffConfig);
  m_MuonEffSFCalo = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_directory, mu_alg_calo, unit, muonEffConfig);
  m_MuonEffSFSA = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_directory, mu_alg_sa, unit, muonEffConfigSA);
  
  if (analysis_version() == "rel_17") {
    m_MuonEffSF->addPeriod("B", int_lum[0]);
    m_MuonEffSF->addPeriod("D", int_lum[1]);
    m_MuonEffSF->addPeriod("E", int_lum[2]);
    m_MuonEffSF->addPeriod("F", int_lum[3]);
    m_MuonEffSF->addPeriod("G", int_lum[4]);
    m_MuonEffSF->addPeriod("H", int_lum[5]);
    m_MuonEffSF->addPeriod("I", int_lum[6]);
    m_MuonEffSF->addPeriod("J", int_lum[7]);
    m_MuonEffSF->addPeriod("K", int_lum[8]);
    m_MuonEffSF->addPeriod("L", int_lum[9]);
    m_MuonEffSF->addPeriod("M", int_lum[10]);
    
    m_MuonEffSFCalo->addPeriod("B", int_lum[0]);
    m_MuonEffSFCalo->addPeriod("D", int_lum[1]);
    m_MuonEffSFCalo->addPeriod("E", int_lum[2]);
    m_MuonEffSFCalo->addPeriod("F", int_lum[3]);
    m_MuonEffSFCalo->addPeriod("G", int_lum[4]);
    m_MuonEffSFCalo->addPeriod("H", int_lum[5]);
    m_MuonEffSFCalo->addPeriod("I", int_lum[6]);
    m_MuonEffSFCalo->addPeriod("J", int_lum[7]);
    m_MuonEffSFCalo->addPeriod("K", int_lum[8]);
    m_MuonEffSFCalo->addPeriod("L", int_lum[9]);
    m_MuonEffSFCalo->addPeriod("M", int_lum[10]);
  }
  
  m_MuonEffSF->Initialise();
  m_MuonEffSFCalo->Initialise();
  m_MuonEffSFSA->Initialise();
  
  // MuonSmearer
  Info("initialize_tools", "Initializing the MuonSmear::SmearingClass tool...");
  if (getMuonFamily() == Muon::STACO) {
    mu_alg = "staco";
  } else if (getMuonFamily() == Muon::MUID) {
    mu_alg = "muid";
  }
  std::string muon_sm_directory("MuonMomentumCorrections/share/");
  if (analysis_version() == "rel_17")
    m_MuonSmearer = new MuonSmear::SmearingClass("Data11", mu_alg, "pT", "Rel17", muon_sm_directory);
  else if (analysis_version() == "rel_17_2")
    m_MuonSmearer = new MuonSmear::SmearingClass("Data12", mu_alg, "pT", "Rel17.2_preliminary", muon_sm_directory);
  m_MuonSmearer->UseScale(1);
  m_MuonSmearer->UseImprovedCombine();
  
  // MuonTrigSF
  Info("initialize_tools", "Initializing the LeptonTriggerSF tool...");
  if (analysis_version() == "rel_17") { // 2011
    m_MuonTrigSF = new LeptonTriggerSF("TrigMuonEfficiency/share/");
  } else if (analysis_version() == "rel_17_2") { // 2012
    m_MuonTrigSF = new LeptonTriggerSF(2012, "TrigMuonEfficiency/share", "muon_trigger_sf_2012.root");
  }
  
  // Muon and Electron trigger match
  m_MuonTriggerMatchTool = new MuonTriggerMatching(&m_trigNavVar);
  m_ElectronTriggerMatchTool = new ElectronTriggerMatching(&m_trigNavVar);
  
  // delay ggF tool initialization
  m_ggFReweighter = 0;
  
  // d0 smearing tool
  TFile *d0_smearing_file = new TFile("HiggsllqqAnalysis/packages/files/d0_reweight/impact_parameter_smearing.root");
  m_smearD0[0] = (TH2F*)d0_smearing_file->Get("smearD0_0");
  m_smearD0[1] = (TH2F*)d0_smearing_file->Get("smearD0_1");
  m_smearD0[2] = (TH2F*)d0_smearing_file->Get("smearD0_2");
  m_smearD0[0]->SetDirectory(0);
  m_smearD0[1]->SetDirectory(0);
  m_smearD0[2]->SetDirectory(0);
  m_smearD0_x = m_smearD0[0]->GetXaxis();
  m_smearD0_y = m_smearD0[0]->GetYaxis();
  d0_smearing_file->Close();
  
  return kTRUE;
}

Float_t HiggsllqqAnalysis::getD0SmearSigma(Int_t index_of_lepton, Int_t nBL, Float_t pt, Float_t eta)
{
  // from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/IPSmearing
  m_smearD0_rand.SetSeed(ntuple->eventinfo.EventNumber() + 100 * index_of_lepton);
  
  if (nBL >= 2) nBL = 2;
  
  float sinTheta = 1. / TMath::CosH(eta);
  float p = pt * TMath::CosH(eta);
  float p_quant = 1. / TMath::Sqrt(pt * pt * sinTheta) / 1000.; // 1./sqrt(p*p*sinTheta*sinTheta*sinTheta)
  
  int Xbin = m_smearD0_x->FindFixBin(eta);
  int Ybin = m_smearD0_y->FindFixBin(p_quant);
  float sigma = m_smearD0[nBL]->GetBinContent(Xbin, Ybin);
  
  return m_smearD0_rand.Gaus(0, sigma);
}

Bool_t HiggsllqqAnalysis::execute_tools(Long64_t entry)
{
  // read a new entry for the trigger decision tool
  m_TrigDecisionToolD3PD->GetEntry(entry, 0);
  
  // initialize the ggF reweighter with the first event info
  if (m_ggFReweighter == 0) {
    if (isMC()) {
      UInt_t chan(ntuple->eventinfo.mc_channel_number());
      if (m_SignalSampleMass.find(chan) != m_SignalSampleMass.end()) {
	Int_t samplemass = m_SignalSampleMass.find(chan)->second;
	Warning("execute_tools", "ggFReweighter being initialized with mass %d according to mc_channel_number = %u taken from the first event (!!!)", samplemass, chan);
	m_ggFReweighter = new ggFReweighting("PowHeg", samplemass, "Mean", "./HiggsllqqAnalysis/packages/files/ggFHiggsPtWeight/", "mc11");
      }
    }
  }
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::finalize_tools()
{
  delete m_GRL;
  delete m_PileupReweighter;
  if (m_VertexPositionReweighter) delete m_VertexPositionReweighter;
  delete m_TrigDecisionToolD3PD;
  delete m_ElectronEnergyRescaler;
  delete m_ElectronEffSF;
  delete m_MuonEffSF;
  delete m_MuonEffSFCalo;
  delete m_MuonEffSFSA;
  delete m_MuonSmearer;
  delete m_MuonTrigSF;
  delete m_MuonTriggerMatchTool;
  delete m_ElectronTriggerMatchTool;
  delete m_smearD0[0];
  delete m_smearD0[1];
  delete m_smearD0[2];
  
  delete m_ResolutionModel;
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::initialize_analysis()
{
  // open the output file
  m_outputFile = new TFile(m_outputFileName, "RECREATE");
  
  // prepare output histo
  m_generatedEntriesHisto = new TH1D("generatedEntriesHisto", "number of processed entries", 10, 0, 10);
  m_generatedEntriesHisto->Fill("raw", 0);
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix", 0);
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix_EvWeight", 0);
  m_generatedEntriesHisto->Fill("with_ggF_as_well", 0);
  m_generatedEntriesHisto->Fill("with_pu_and_vertex", 0);
  //m_truthHistos["truth_mu_pt"] = new TH1F("truth_mu_pt", "truth muon pT;truth muon p_{T} [GeV];entries", 100, 0, 100); // kostas
  //m_truthHistos["truth_mu_eta"] = new TH1F("truth_mu_eta", "truth muon #eta;truth muon #eta;entries", 100, -5, 5); // kostas
  //m_truthHistos["truth_H_pt"] = new TH1F("truth_H_pt", "truth Higgs pT;truth Higgs p_{T} [GeV];entries", 100, 0, 100); // kostas
  //m_truthHistos["truth_H_eta"] = new TH1F("truth_H_eta", "truth Higgs #eta;truth Higgs #eta;entries", 100, -5, 5); // kostas
  m_selectionEfficiencyVsNvx[0] = new TEfficiency("selectionEfficiencyVsNvx_MU2", "selectionEfficiencyVsNvx", 40, 0, 40); // kostas2
  m_selectionEfficiencyVsNvx[1] = new TEfficiency("selectionEfficiencyVsNvx_MUE", "selectionEfficiencyVsNvx", 40, 0, 40); // kostas2
  m_selectionEfficiencyVsNvx[2] = new TEfficiency("selectionEfficiencyVsNvx_E2", "selectionEfficiencyVsNvx", 40, 0, 40); // kostas2
  
  // prepare output trees
  m_TreeCutflow = new TTree("cutflow", "minimal event cutflow ntuple");
  m_TreeCutflow->Branch("analysis", &(m_cutflowStruct.analysis), "analysis/I");
  m_TreeCutflow->Branch("run", &(m_cutflowStruct.run), "run/i");
  m_TreeCutflow->Branch("event", &(m_cutflowStruct.event), "event/i");
  m_TreeCutflow->Branch("last", &(m_cutflowStruct.last), "last/I");
  
  
  m_TreeCandidates = new TTree("candidates", "SFOS Higgs candidates");
  setTreeCandidatesBranches(m_TreeCandidates, &m_candStruct);
  
  // reset the convenience event counter
  m_processedEntries = 0;
  
  // initialize cross-section and higgs-mass maps
  initCrossSections();
  
  // initialize the cutflows
  m_Channels.push_back(HiggsllqqAnalysis::E2);
  m_Channels.push_back(HiggsllqqAnalysis::MU2);
  m_Channels.push_back(HiggsllqqAnalysis::MUE);
  
  std::map<Int_t, TString> chan_name;
  chan_name.clear();
  chan_name[HiggsllqqAnalysis::E2]  = "E2";
  chan_name[HiggsllqqAnalysis::MU2] = "MU2";
  chan_name[HiggsllqqAnalysis::MUE] = "MUE";
  
  m_EventCutflow.clear();
  m_EventCutflow_rw.clear();
  m_CandidateCutflow.clear();
  m_ElectronCutflow.clear();
  m_MuonCutflow.clear();
  m_JetCutflow.clear();
  
  for (UInt_t i = 0; i < m_Channels.size(); i++) {
    
    m_EventCutflow.push_back(Analysis::CutFlowTool("Plain_" + chan_name[i]));
    m_EventCutflow[i].addCut("HFOR");
    m_EventCutflow[i].addCut("GRL");
    m_EventCutflow[i].addCut("larError");
    m_EventCutflow[i].addCut("Trigger");
    m_EventCutflow[i].addCut("Vertex");
    m_EventCutflow[i].addCut("METcleaning");
    m_EventCutflow[i].addCut("LArHole");
    m_EventCutflow[i].addCut("NumberOfLeptons");
    m_EventCutflow[i].addCut("OppositeSign");
    m_EventCutflow[i].addCut("PtLeptons");
    m_EventCutflow[i].addCut("TriggerConsistency");
    m_EventCutflow[i].addCut("MET");
    m_EventCutflow[i].addCut("TwoJets");
    m_EventCutflow[i].addCut("NumTagJets");
    m_EventCutflow[i].addCut("DileptonMass");
    m_EventCutflow[i].addCut("DiJetMass");
    m_EventCutflow[i].addCut("SelectedLeptons");
    m_EventCutflow[i].addCut("GoodCandidate");
    
    m_EventCutflow_rw.push_back(Analysis::CutFlowTool("Reweighted_" + chan_name[i]));
    m_EventCutflow_rw[i].addCut("HFOR");
    m_EventCutflow_rw[i].addCut("GRL");
    m_EventCutflow_rw[i].addCut("larError");
    m_EventCutflow_rw[i].addCut("Trigger");
    m_EventCutflow_rw[i].addCut("Vertex");
    m_EventCutflow_rw[i].addCut("METcleaning");
    m_EventCutflow_rw[i].addCut("LArHole");
    m_EventCutflow_rw[i].addCut("NumberOfLeptons");
    m_EventCutflow_rw[i].addCut("OppositeSign");
    m_EventCutflow_rw[i].addCut("PtLeptons");
    m_EventCutflow_rw[i].addCut("TriggerConsistency");
    m_EventCutflow_rw[i].addCut("MET");
    m_EventCutflow_rw[i].addCut("TwoJets");
    m_EventCutflow_rw[i].addCut("NumTagJets");
    m_EventCutflow_rw[i].addCut("DileptonMass");
    m_EventCutflow_rw[i].addCut("DiJetMass");
    m_EventCutflow_rw[i].addCut("SelectedLeptons");
    m_EventCutflow_rw[i].addCut("GoodCandidate");
    
    m_CandidateCutflow.push_back(Analysis::CutFlowTool("Candidates_" + chan_name[i]));
    m_CandidateCutflow[i].addCut("mass_2e2mu");
    m_CandidateCutflow[i].addCut("opposite_charge");
    m_CandidateCutflow[i].addCut("kinematics");
    m_CandidateCutflow[i].addCut("trigmatch");
    m_CandidateCutflow[i].addCut("Z1_mass");
    m_CandidateCutflow[i].addCut("Z2_mass");
    m_CandidateCutflow[i].addCut("DeltaR");
    m_CandidateCutflow[i].addCut("best");
    m_CandidateCutflow[i].addCut("track_iso");
    m_CandidateCutflow[i].addCut("calo_iso");
    m_CandidateCutflow[i].addCut("d0_sig");
    m_CandidateCutflow[i].addCut("mass_window");
    
    m_ElectronCutflow.push_back(Analysis::CutFlowTool("Electrons_" + chan_name[i]));
    m_ElectronCutflow[i].addCut("family");
    m_ElectronCutflow[i].addCut("algorithm");
    m_ElectronCutflow[i].addCut("quality");
    m_ElectronCutflow[i].addCut("eta");
    m_ElectronCutflow[i].addCut("Et");
    m_ElectronCutflow[i].addCut("objectquality");
    m_ElectronCutflow[i].addCut("z0");
    m_ElectronCutflow[i].addCut("overlap");
    
    m_MuonCutflow.push_back(Analysis::CutFlowTool("Muons_" + chan_name[i]));
    m_MuonCutflow[i].addCut("family");
    m_MuonCutflow[i].addCut("quality");
    m_MuonCutflow[i].addCut("cosmic");
    m_MuonCutflow[i].addCut("eta");
    m_MuonCutflow[i].addCut("pt");
    m_MuonCutflow[i].addCut("MCP");
    m_MuonCutflow[i].addCut("z0");
    m_MuonCutflow[i].addCut("overlap");
    
    m_JetCutflow.push_back(Analysis::CutFlowTool("Jets_" + chan_name[i]));
    m_JetCutflow[i].addCut("jetCleaning");
    m_JetCutflow[i].addCut("kinematics");
    m_JetCutflow[i].addCut("Pileup");
    m_JetCutflow[i].addCut("overlap");
  }
  
  // create the Z mass constraint resolution model
  m_ResolutionModel = new ResolutionModel(Z_pdg_mass, Z_pdg_width);
  m_ResolutionModel->LoadTemplatesFromFile("HiggsllqqAnalysis/packages/files/ResolutionModel/z_truth_mass_templates.root");
  
  return kTRUE;
}


Bool_t HiggsllqqAnalysis::hasPowHegZZBug()
{
  // check if an event in PowHeg ZZ sample has the bug
  Bool_t result(kFALSE);
  
  if (isMC()) {
    if (ntuple->eventinfo.mc_channel_number() == 126859 || ntuple->eventinfo.mc_channel_number() == 126860) {
      for (Int_t i = 0; i < ntuple->mc.n(); i++) {
	if (TMath::Abs(ntuple->mc[i].pdgId()) == 23 && ntuple->mc[i].m() < 0)
	  result = kTRUE;
      } // loop over truth particles
    } // PowHeg ZZ samples with bug
  }
  
  return result;
}

Bool_t HiggsllqqAnalysis::isWithinT1llqqPhaseSpace()
{
  // check if an event in T1 falls in the phase space of the 4-lepton filtered sample
  
  if (isMC()) {
    if (ntuple->eventinfo.mc_channel_number() == 105200) {
      std::vector<TLorentzVector> leptons;
      
      for (Int_t i = 0; i < ntuple->mc.n(); i++) {
	D3PDReader::TruthParticleD3PDObjectElement *p = &(ntuple->mc[i]);
	if (p->pt() > 5000 && p->status() == 1 && (TMath::Abs(p->pdgId()) == 13 || TMath::Abs(p->pdgId()) == 11)) {
	  leptons.push_back(CommonTools::getVector(p));
	}
      }
      
      std::vector<TLorentzVector> dileptons;
      
      for (UInt_t i = 0; i < leptons.size(); i++) {
	for (UInt_t j = i + 1; j < leptons.size(); j++) {
	  dileptons.push_back(leptons[i] + leptons[j]);
	}
      }
      
      for (UInt_t i = 0; i < dileptons.size(); i++) {
	for (UInt_t j = i + 1; j < dileptons.size(); j++) {
	  Double_t lowmass_cut = -1.;
	  if (analysis_version() == "rel_17") lowmass_cut = 60000.;
	  else if (analysis_version() == "rel_17_2") lowmass_cut = 40000.;
	  
	  if ((dileptons[i].M() > lowmass_cut && dileptons[j].M() > 12000) || (dileptons[j].M() > lowmass_cut && dileptons[i].M() > 12000)) {
	    return kTRUE;
	  }
	}
      }
    }
  }
  
  return kFALSE;
}

Int_t HiggsllqqAnalysis::getLastCutPassed()
{
  Int_t last(-1);
  
  // Including hfor veto 
  HFOR_value = hforTool->getDecision(ntuple->eventinfo.mc_channel_number(),mc_n, mc_pt, mc_eta, mc_phi, mc_m, mc_pdgId, mc_status, mc_vx_barcode, mc_parent_index, mc_child_index, HforToolD3PD::BBONLY);
  
  if (HFOR_value!=4) last = HllqqCutFlow::HFOR;
  else return last;   
  if (passesGRL()) last = HllqqCutFlow::GRL;
  else return last;
  if (ntuple->eventinfo.larError()  <= 1 /*!= 2*/) last = HllqqCutFlow::larError;
  else return last;
  if (passesTrigger()) last = HllqqCutFlow::Trigger;
  else return last;
  if (hasGoodVertex()) last = HllqqCutFlow::Vertex;
  else return last;
  
  getGoodLeptons();

  if(NotMETclean() /*&& !isMC()*/ && DoMETdataClean) return last;
  else last = HllqqCutFlow::METcleaning;
  
  if((JetInHole() && isMC() && getPeriod()==DataPeriod::y2011_EH) || (JetInHole() && !isMC())) return last;
  else last = HllqqCutFlow::LArHole;
  
  getGoodLeptons();
  
  //Number of Leptons Cut
  if ((getChannel() == HiggsllqqAnalysis::MU2 && m_GoodMuons.size() == 2 && m_GoodElectrons.size() == 0) ||
      (getChannel() == HiggsllqqAnalysis::E2 && m_GoodElectrons.size() == 2 && m_GoodMuons.size() == 0)  ||
      (getChannel() == HiggsllqqAnalysis::MUE && m_GoodMuons.size() == 1 && m_GoodElectrons.size() == 1)) last = HllqqCutFlow::NumberOfLeptons;  
  else return last;
  
  //OS selection on the 2 analysis channels
  int chargeprod = 0;
  if(getChannel() == HiggsllqqAnalysis::MU2 && !GetDoQCDSelection())
    chargeprod = (m_GoodMuons.at(1))->charge()*(m_GoodMuons.at(0))->charge(); 
  else if(getChannel() == HiggsllqqAnalysis::E2 && !GetDoQCDSelection())
    chargeprod = (m_GoodElectrons.at(1))->charge()*(m_GoodElectrons.at(0))->charge(); 
  if (chargeprod==-1 || getChannel() == HiggsllqqAnalysis::MUE || GetDoQCDSelection()) last = HllqqCutFlow::OppositeSign;
  else return last;
  
  //PtLeptons cut
  if(IsConsistentPt() || NotPtConsistency) last = HllqqCutFlow::PtLeptons;
  else return last;
  
  //Trigger Consistency Cut
  if(IsConsistentWithTrigger() || NotTrigConsistency) last = HllqqCutFlow::TriggerConsistency;
  else return last;
  
  //Missing Et cut
  if(hasGoodMET()) last = HllqqCutFlow::MET;
  else return last;
  
  //Minimum number of Jets cut
  if(m_GoodJets.size()>=2) last = HllqqCutFlow::TwoJets;
  else return last;
  
  //Minimun number of tag jets
  if((GetNumOfTags()==2 && TaggedChannel) || (GetNumOfTags()<2 && !TaggedChannel)) last = HllqqCutFlow::NumTagJets;
  else return last;
  
  //Dielpton Mass windows
  getDileptons();  
  std::vector<Analysis::Dilepton *>::iterator Z_itr_i= m_Dileptons.begin();
  Double_t DilepMass = (*Z_itr_i)->Get4Momentum()->M();
  
  if((GetDoLowMass() && DilepMass>20000. && DilepMass<76000.) ||
     (!GetDoLowMass() && DilepMass>83000. && DilepMass<99000.)) last = HllqqCutFlow::DileptonMass;
  else return last;
  
  //Invariant mass of the dijet
  if((JetKinematicFitterResult() && !TaggedChannel) || (JetDimassTagged() && TaggedChannel))last = HllqqCutFlow::DiJetMass;
  else return last;
  
  
  if ((getChannel() == HiggsllqqAnalysis::MU2 && m_GoodMuons.size() >= 4)    ||
      (getChannel() == HiggsllqqAnalysis::E2 && m_GoodElectrons.size() >= 4) ||
      (getChannel() == HiggsllqqAnalysis::MUE && m_GoodMuons.size() >= 2 && m_GoodElectrons.size() >= 2)) last = HllqqCutFlow::SelectedLeptons;
  else return last;
  
  getGoodObjects();
  
  if (m_GoodQuadrileptons.size() >= 1) last = HllqqCutFlow::GoodCandidate;
  else return last;
  
  return last;
}

Bool_t HiggsllqqAnalysis::passesGRL()
{
  if (isMC()) {
    return kTRUE;
  } else {
    if (ntuple->eventinfo.RunNumber() >= 205113) return kTRUE; // kostas (new data)
    
    return (m_GRL->GetGoodRunsList(getChannel() + HllqqGRL::data11).HasRunLumiBlock(ntuple->eventinfo.RunNumber(), ntuple->eventinfo.lbn()) ||
	    m_GRL->GetGoodRunsList(getChannel() + HllqqGRL::data12).HasRunLumiBlock(ntuple->eventinfo.RunNumber(), ntuple->eventinfo.lbn()));
  }
}

Bool_t HiggsllqqAnalysis::hasGoodVertex()
{
  for (Int_t i = 0; i < ntuple->vxp.n(); i++) {
    if (ntuple->vxp[i].trk_n() >= 3) return kTRUE;
  }
  return kFALSE;
}

Int_t HiggsllqqAnalysis::getNumberOfGoodVertices()
{
  Int_t result(0);
  
  for (Int_t i = 0; i < ntuple->vxp.n(); i++) {
    if (ntuple->vxp[i].trk_n() >= 3) result++;
  }
  
  return result;
}

std::vector<TString> HiggsllqqAnalysis::getListOfAlternativeTriggers(TString sequence)
{
  std::vector<TString> result;
  
  TObjArray *chains = sequence.Tokenize(";"); // expect e.g. EF_mu18;EF_mu18_MG (; means OR)
  
  if (chains->GetEntriesFast()) {
    TIter iChain(chains);
    TObjString *os = 0;
    
    while ((os = (TObjString*)iChain())) {
      result.push_back(os->GetString().Data());
    } // loop over chains
  } // non-zero input
  
  delete chains;
  
  return result;
}

Bool_t HiggsllqqAnalysis::passesTrigger()
{
  if (getChannel() == HiggsllqqAnalysis::MU2) {
    return (passesSingleMuonTrigger() || passesDiMuonTrigger());
  } else if (getChannel() == HiggsllqqAnalysis::MUE) {
    return (passesSingleMuonTrigger() || passesDiMuonTrigger() || passesSingleElectronTrigger() || passesDiElectronTrigger());// || passesElectronMuonTrigger());
  } else if (getChannel() == HiggsllqqAnalysis::E2) {
    return (passesSingleElectronTrigger() || passesDiElectronTrigger());
  } else {
    Abort("Unexpected channel passed to HiggsllqqAnalysis::passesTrigger()");
    return kFALSE;
  }
}

Bool_t HiggsllqqAnalysis::passesSingleMuonTrigger()
{
  return (getTriggerInfo(getSingleMuonTriggerName()) == 1);
}

Bool_t HiggsllqqAnalysis::passesDiMuonTrigger()
{
  return (getTriggerInfo(getDiMuonTriggerName()) == 1);
}

Bool_t HiggsllqqAnalysis::passesSingleElectronTrigger()
{
  Int_t period = getPeriod();
  
  if (period == DataPeriod::y2011_L || period == DataPeriod::y2011_M) {
    return (getTriggerInfo(getSingleElectronTriggerName()) == 1 || getTriggerInfo("EF_e45_medium1") == 1);
  }
  else {
    return (getTriggerInfo(getSingleElectronTriggerName()) == 1);
  }
}

Bool_t HiggsllqqAnalysis::passesDiElectronTrigger()
{
  return (getTriggerInfo(getDiElectronTriggerName()) == 1);
}

Bool_t HiggsllqqAnalysis::passesElectronMuonTrigger()
{
  return (getTriggerInfo(getElectronMuonTriggerName()) == 1);
}

TString HiggsllqqAnalysis::getSingleMuonTriggerName()
{
  TString chain_name("");
  Int_t period = getPeriod();
  
  if (!isMC()) {
    // DATA
    if (period <= DataPeriod::y2011_I) {
      // A-I
      chain_name = "EF_mu18_MG";
    } else if (period <= DataPeriod::y2011_M) {
      // L-M
      chain_name = "EF_mu18_MG_medium";
    } else if (period <= DataPeriod::y2012_B) {
      // 2012: A-B
      chain_name = "EF_mu24i_tight;EF_mu36_tight";
    }
  } else {
    // MC
    if (period == DataPeriod::y2011_BD) {
      // BD
      chain_name = "EF_mu18_MG";
    } else if (period == DataPeriod::y2011_EH) {
      // EH
      chain_name = "EF_mu18_MG";
    } else if (period == DataPeriod::y2011_I) {
      // I
      chain_name = "EF_mu18_MG";
    } else if (period == DataPeriod::y2011_J) {
      // J
      chain_name = "EF_mu18_MG_medium";
    } else if (period == DataPeriod::y2011_K) {
      // K
      chain_name = "EF_mu18_MG_medium";
    } else if (period == DataPeriod::y2011_LM) {
      // LM
      chain_name = "EF_mu18_MG_medium";
    } else if (period == DataPeriod::y2012_AllYear) {
      // 2012: AllYear
      
      chain_name = "EF_mu24i_tight;EF_mu36_tight";
    }
  }
  
  return chain_name;
}

TString HiggsllqqAnalysis::getDiMuonTriggerName()
{
  TString chain_name("");
  Int_t period = getPeriod();
  
  if (!isMC()) {
    // DATA
    if (period <= DataPeriod::y2011_M) {
      // A-M
      chain_name = "EF_2mu10_loose";
    } else {
      // 2012: A-B
      chain_name = "EF_2mu13;EF_mu18_tight_mu8_EFFS";
    }
  } else {
    // MC
    if (period == DataPeriod::y2011_BD) {
      chain_name = "EF_2mu10_loose";
    } else if (period == DataPeriod::y2011_EH) {
      chain_name = "EF_2mu10_loose";
    } else if (period == DataPeriod::y2011_I) {
      chain_name = "EF_2mu10_loose";
    } else if (period == DataPeriod::y2011_J) {
      chain_name = "EF_2mu10_loose";
    } else if (period == DataPeriod::y2011_K) {
      chain_name = "EF_2mu10_loose";
    } else if (period == DataPeriod::y2011_LM) {
      chain_name = "EF_2mu10_loose";
    } else if (period == DataPeriod::y2012_AllYear) {
      // 2012: AllYear
      chain_name = "EF_2mu13;EF_mu18_tight_mu8_EFFS";
    }
  }
  
  return chain_name;
}

TString HiggsllqqAnalysis::getSingleElectronTriggerName()
{
  TString chain_name("");
  Int_t period = getPeriod();
  
  if (!isMC()) {
    // DATA
    // A-J
    if (period <= DataPeriod::y2011_J) {
      chain_name = "EF_e20_medium";
    }
    // K
    else if (period == DataPeriod::y2011_K) {
      chain_name = "EF_e22_medium";
    }
    // L-M
    else if (period <= DataPeriod::y2011_M) {
      chain_name = "EF_e22vh_medium1";
    } else if (period == DataPeriod::y2012_A) {
      chain_name = "EF_e22vh_medium1;EF_e45_medium1";
    } else if (period == DataPeriod::y2012_B) {
      chain_name = "EF_e24vhi_medium1;EF_e60_medium1";
    }
  } else {
    // MC
    if (period == DataPeriod::y2011_BD) {
      chain_name = "EF_e20_medium";
    } else if (period == DataPeriod::y2011_EH) {
      chain_name = "EF_e20_medium";
    } else if (period == DataPeriod::y2011_I) {
      chain_name = "EF_e20_medium";
    } else if (period == DataPeriod::y2011_J) {
      chain_name = "EF_e20_medium";
    } else if (period == DataPeriod::y2011_K) {
      chain_name = "EF_e22_medium";
    } else if (period == DataPeriod::y2011_LM) {
      chain_name = "EF_e22vh_medium1";
    } else if (period == DataPeriod::y2012_AllYear) {
      // 2012: AllYear
      chain_name = "EF_e24vhi_medium1;EF_e60_medium1";
    }
  }
  
  return chain_name;
}

TString HiggsllqqAnalysis::getDiElectronTriggerName()
{
  TString chain_name("");
  Int_t period = getPeriod();
  
  if (!isMC()) {
    // DATA
    // A-J
    if (period <= DataPeriod::y2011_J) {
      chain_name = "EF_2e12_medium";
    }
    // K
    else if (period == DataPeriod::y2011_K) {
      chain_name = "EF_2e12T_medium";
    }
    // L-M
    else if (period <= DataPeriod::y2011_M) {
      chain_name = "EF_2e12Tvh_medium";
    }
    // 2012: A-B
    else if (period <= DataPeriod::y2012_B) {
      chain_name = "EF_2e12Tvh_loose1";
    }
  } else {
    // MC
    if (period == DataPeriod::y2011_BD) {
      chain_name = "EF_2e12_medium";
    } else if (period == DataPeriod::y2011_EH) {
      chain_name = "EF_2e12_medium";
    } else if (period == DataPeriod::y2011_I) {
      chain_name = "EF_2e12T_medium"; //"EF_2e12_medium";
    } else if (period == DataPeriod::y2011_J) {
      chain_name = "EF_2e12T_medium"; //"EF_2e12_medium";
    } else if (period == DataPeriod::y2011_K) {
      chain_name = "EF_2e12T_medium";
    } else if (period == DataPeriod::y2011_LM) {
      chain_name = "EF_2e12Tvh_medium";
    } else if (period == DataPeriod::y2012_AllYear) {
      // 2012: AllYear
      chain_name = "EF_2e12Tvh_loose1";
    }
  }
  
  return chain_name;
}

TString HiggsllqqAnalysis::getElectronMuonTriggerName()
{
  TString chain_name("");
  Int_t period = getPeriod();
  
  if (!isMC()) {
    // DATA
    // 2011
    if (period <= DataPeriod::y2011_M) {
      chain_name = "DO_NOT_USE_E_MU_TRIGGER";
    }
    // 2012: A-B
    else if (period <= DataPeriod::y2012_B) {
      chain_name = "EF_e12Tvh_medium1_mu8;EF_e24vhi_loose1_mu8";
    }
  } else {
    // MC
    if (period == DataPeriod::y2011_BD) {
      chain_name = "DO_NOT_USE_E_MU_TRIGGER";
    } else if (period == DataPeriod::y2011_EH) {
      chain_name = "DO_NOT_USE_E_MU_TRIGGER";
    } else if (period == DataPeriod::y2011_I) {
      chain_name = "DO_NOT_USE_E_MU_TRIGGER";
    } else if (period == DataPeriod::y2011_J) {
      chain_name = "DO_NOT_USE_E_MU_TRIGGER";
    } else if (period == DataPeriod::y2011_K) {
      chain_name = "DO_NOT_USE_E_MU_TRIGGER";
    } else if (period == DataPeriod::y2011_LM) {
      chain_name = "DO_NOT_USE_E_MU_TRIGGER";
    } else if (period == DataPeriod::y2012_AllYear) {
      // 2012: AllYear
      chain_name = "EF_e12Tvh_medium1_mu8;EF_e24vhi_loose1_mu8";
    }
  }
  
  return chain_name;
}

void HiggsllqqAnalysis::applyChanges(Analysis::ChargedLepton *lep)
{
  if (lep->flavor() == Analysis::ChargedLepton::MUON) {
    D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
    
    if (isMC() && doSmearing() && MuonSmearing) {
      double eta  = lep->Get4Momentum()->Eta();
      double ptcb = lep->Get4Momentum()->Pt();
      double ptme = lep->Get4Momentum_SA()->Pt();
      double ptid = lep->Get4Momentum_ID()->Pt();
      
      m_MuonSmearer->SetSeed(ntuple->eventinfo.EventNumber(), (Int_t)mu->GetIndex());
      
      if (mu->isCombinedMuon()) {
	m_MuonSmearer->Event(ptme, ptid, ptcb, eta);
	lep->Get4Momentum()->SetPtEtaPhiM(m_MuonSmearer->pTCB(), lep->Get4Momentum()->Eta(), lep->Get4Momentum()->Phi(), lep->Get4Momentum()->M());
	lep->Get4Momentum_SA()->SetPtEtaPhiM(m_MuonSmearer->pTMS(), lep->Get4Momentum_SA()->Eta(), lep->Get4Momentum_SA()->Phi(), lep->Get4Momentum_SA()->M());
	lep->Get4Momentum_ID()->SetPtEtaPhiM(m_MuonSmearer->pTID(), lep->Get4Momentum_ID()->Eta(), lep->Get4Momentum_ID()->Phi(), lep->Get4Momentum_ID()->M());
      } else if (mu->isStandAloneMuon()) {
	m_MuonSmearer->Event(ptcb, eta, "MS");
	lep->Get4Momentum()->SetPtEtaPhiM(m_MuonSmearer->pTMS(), lep->Get4Momentum()->Eta(), lep->Get4Momentum()->Phi(), lep->Get4Momentum()->M());
	lep->Get4Momentum_SA()->SetPtEtaPhiM(m_MuonSmearer->pTMS(), lep->Get4Momentum_SA()->Eta(), lep->Get4Momentum_SA()->Phi(), lep->Get4Momentum_SA()->M());
      } else { // segment tagged, calo
	m_MuonSmearer->Event(ptcb, eta, "ID");
	lep->Get4Momentum()->SetPtEtaPhiM(m_MuonSmearer->pTID(), lep->Get4Momentum()->Eta(), lep->Get4Momentum()->Phi(), lep->Get4Momentum()->M());
	lep->Get4Momentum_ID()->SetPtEtaPhiM(m_MuonSmearer->pTID(), lep->Get4Momentum_ID()->Eta(), lep->Get4Momentum_ID()->Phi(), lep->Get4Momentum_ID()->M());
      }
    }
  } else if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    Int_t SYST_FLAG = 0; //SYST_FLAG is 0 for nominal scale, 1 or 2 for 1-sigma variations.
    
    // first of all one must rescale energy in the crack! (both data and MC but only for 2011 so far)
    Float_t tmp_calibration(1.);
    if (analysis_version() == "rel_17") // rel. 17
      tmp_calibration = m_ElectronEnergyRescaler->applyMCCalibrationMeV(el->cl_eta(), el->cl_E() / TMath::CosH(el->tracketa()), "ELECTRON");
    Float_t tmp_E = el->cl_E() * TMath::Abs(tmp_calibration);
    Float_t tmp_Et = tmp_E / TMath::CosH(el->tracketa());
    Float_t tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    
    // then, apply the other corrections
    if (isMC() && doSmearing() && ElectronSmearing) {
      m_ElectronEnergyRescaler->SetRandomSeed(ntuple->eventinfo.EventNumber() + 100 * (Int_t)el->GetIndex());
      
      // false here means the MC is mc11c (no constant term)
      Float_t smearcorr = m_ElectronEnergyRescaler->getSmearingCorrectionMeV(el->cl_eta(), tmp_E, SYST_FLAG, false, "2011");
      
      tmp_E = tmp_E * smearcorr;
      tmp_Et = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    if (!isMC()) {
      tmp_E = m_ElectronEnergyRescaler->applyEnergyCorrectionMeV(el->cl_eta(), el->cl_phi(), tmp_E, tmp_Et, SYST_FLAG, "ELECTRON");
      tmp_Et = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    
    Float_t tmp_pt = (TMath::Sqrt(TMath::Power(tmp_E, 2) - TMath::Power(lep->Get4Momentum()->M(), 2))) * TMath::Sin(2 * TMath::ATan(TMath::Exp(-el->tracketa())));
    Float_t tmp_pt_cl = (TMath::Sqrt(TMath::Power(tmp_E, 2) - TMath::Power(lep->Get4Momentum()->M(), 2))) * TMath::Sin(2 * TMath::ATan(TMath::Exp(-el->cl_eta())));
    
    // correct lepton 4-momenta
    lep->Get4Momentum()->SetPtEtaPhiM(tmp_pt, el->tracketa(), el->trackphi(), lep->Get4Momentum()->M());
    lep->Get4Momentum_ID()->SetPtEtaPhiM(el->trackpt(), el->tracketa(), el->trackphi(), lep->Get4Momentum()->M()); // untouched by energy corrections!
    lep->Get4Momentum_SA()->SetPtEtaPhiM(tmp_pt_cl, el->cl_eta(), el->cl_phi(), lep->Get4Momentum()->M());
  } else {
    Abort("Unknown lepton flavor");
  }
}

void HiggsllqqAnalysis::applyChanges(Analysis::Jet *jet)
{
  Float_t tmp_E(-9999.9), tmp_pt(-9999.9), tmp_eta(-9999.9), tmp_phi(-9999.9), tmp_Et(-9999.9);
  Bool_t sysstudy = GetSysStudy(); 
  
  D3PDReader::JetD3PDObjectElement *this_jet = jet->GetJet();
  
  if (isMC() && doSmearing() && JetSmearing) {
    double eta     = this_jet->EtaOrigin();
    double phi     = this_jet->PhiOrigin();
    double Eraw    = this_jet->emscale_E();
    double eta_det = this_jet->emscale_eta();
    double m       = this_jet->MOrigin();
    
    double mu = ntuple->eventinfo.averageIntPerXing();
    int NPV=0;
    for (Int_t i = 0; i < ntuple->vxp.n(); i++) {
      if (ntuple->vxp[i].trk_n() >= 3) NPV++;
    }
    
    // Calibrate the jet!
    // Pile-up, origin, EtaJES correction applied, i.e. to OFFSET_ORIGIN_ETAJES scale
    TLorentzVector jet4v = myJES->ApplyOffsetEtaJES(Eraw,eta_det,eta,phi,m,mu,NPV);  
    
    // The below is systematic evaluation, and ONLY for MC
    // Smear the jet to match the MC resolution+1 sigma!
    if(isMC() && sysstudy){
      myJER->SetSeed(ntuple->eventinfo.EventNumber());
      myJER->SmearJet_SystRel17(jet4v);
    }
    
    tmp_E=jet4v.E();
    tmp_pt=jet4v.Pt();
    tmp_eta=jet4v.Eta();
    tmp_phi=jet4v.Phi();
    tmp_Et=jet4v.Et();
    
    //Inside the "if"
    jet->set_rightE(tmp_E);
    jet->set_rightpt(tmp_pt);
    jet->set_righteta(tmp_eta);
    jet->set_rightphi(tmp_phi);
    jet->set_rightEt(tmp_Et);
  }  
}

HllqqSystematics::ChargedLepton HiggsllqqAnalysis::getLeptonSystematics(Analysis::ChargedLepton *lep)
{
  HllqqSystematics::ChargedLepton result = HllqqSystematics::ChargedLepton();
  
  if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    // apply crack calibration if appropriate
    Float_t tmp_calibration(1.);
    if (analysis_version() == "rel_17") // rel. 17
      tmp_calibration = m_ElectronEnergyRescaler->applyMCCalibrationMeV(el->cl_eta(), el->cl_E() / TMath::CosH(el->tracketa()), "ELECTRON");
    Float_t tmp_E = el->cl_E() * TMath::Abs(tmp_calibration);
    
    // then, apply smearing
    if (isMC() && doSmearing()) {
      m_ElectronEnergyRescaler->SetRandomSeed(ntuple->eventinfo.EventNumber() + 100 * (Int_t)el->GetIndex());
      
      // false here means the MC is mc11c (no constant term)
      // SYST_FLAG is set to zero
      Float_t smearcorr = m_ElectronEnergyRescaler->getSmearingCorrectionMeV(el->cl_eta(), tmp_E, 0, false);
      
      tmp_E = tmp_E * smearcorr;
    }
    
    result.ptCB_nosmearnoscale = -9999.9;
    result.ptME_nosmearnoscale = -9999.9;
    result.ptID_nosmearnoscale = -9999.9;
    result.cl_E_calibsmearnoscale = tmp_E;
    result.cl_eta = el->cl_eta();
    result.cl_phi = el->cl_phi();
  } else {
    D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
    
    result.ptCB_nosmearnoscale = mu->pt();
    result.ptME_nosmearnoscale = TMath::Abs(1. / mu->me_qoverp() * TMath::Sin(mu->me_theta()));
    result.ptID_nosmearnoscale = TMath::Abs(1. / mu->id_qoverp() * TMath::Sin(mu->id_theta()));
    result.cl_E_calibsmearnoscale = -9999.9;
    result.cl_eta = -9999.9;
    result.cl_phi = -9999.9;
  }
  
  return result;
}

void HiggsllqqAnalysis::getMuons(D3PDReader::MuonD3PDObject *mu_branch, Int_t family)
{
  // fills m_Muons with those muons passing the one-lepton muon selection
  // (no overlap removal at this stage)
  
  for (Int_t i = 0; i < mu_branch->n(); i++) {
    Analysis::ChargedLepton *lep = new Analysis::ChargedLepton(&((*mu_branch)[i]), family);
    //applyChanges(lep);
    ApplyChangesMuon(lep);
    
    if(is4lepGood && isGood(lep)) m_Muons.push_back(lep);
    else if(!is4lepGood && IsGoodMuon(lep)) m_Muons.push_back(lep);
    else {
      m_MuonCutflow[getChannel()].addCutCounter(lep->lastcut(), 1);
      delete lep;
    }
  }
}

void HiggsllqqAnalysis::getElectrons(D3PDReader::ElectronD3PDObject *el_branch, Int_t family)
{
  // fills m_Electrons with those electrons passing the one-lepton electron selection
  // (no overlap removal at this stage)
  
  for (Int_t i = 0; i < el_branch->n(); i++) {
    Analysis::ChargedLepton *lep = new Analysis::ChargedLepton(&((*el_branch)[i]), family);
    //applyChanges(lep);
    ApplyChangesElectron(lep);
    
    if (is4lepGood && isGood(lep)) m_Electrons.push_back(lep);
    else if (!is4lepGood && IsGoodElectron(lep)) m_Electrons.push_back(lep);
    else {
      m_ElectronCutflow[getChannel()].addCutCounter(lep->lastcut(), 1);
      delete lep;
    }
  }
}

void HiggsllqqAnalysis::getJets(D3PDReader::JetD3PDObject *jet_branch)
{
  // fills m_Jets with those jets passing the one-jet selection
  // (no overlap removal at this stage)
  
  for (Int_t i = 0; i < jet_branch->n(); i++) {
    Analysis::Jet *jet = new Analysis::Jet(&((*jet_branch)[i]));
    //applyChanges(jet);
    ApplyChangesJet(jet);
    if (isGoodJet(jet)) m_Jets.push_back(jet);
    else {
      m_JetCutflow[getChannel()].addCutCounter(jet->lastcut(), 1);
      
      delete jet;
    }
  }  
}

void HiggsllqqAnalysis::getGoodMuons()
{
  // fills m_GoodMuons with those muons passing the full muon selection
  // (quality and overlap removal)
  
  // utility vector, tells if i-th muons must be skipped since overlapping with a jet
  std::vector<Bool_t> skip_muon;
  
  // remove overlap wrt jets
  skip_muon.clear();
  
  std::vector<Analysis::ChargedLepton *>::iterator muon_itr;
  for (muon_itr = m_Muons.begin(); muon_itr != m_Muons.end(); ++muon_itr) {
    skip_muon.push_back(kFALSE);
  }
  
  Int_t i(-1);
  
  for (std::vector<Analysis::ChargedLepton *>::iterator mu_itr_i = m_Muons.begin(); mu_itr_i != m_Muons.end(); ++mu_itr_i) {
    Analysis::ChargedLepton *mu_i = *mu_itr_i;
    
    Bool_t keepMe(kTRUE);
    
    // calo-staco
    if (mu_i->family() == Muon::CALO && DoCaloMuons) {
      for (std::vector<Analysis::ChargedLepton *>::iterator mu_itr_j = m_Muons.begin(); mu_itr_j != m_Muons.end(); ++mu_itr_j) {
	if (mu_itr_j == mu_itr_i) continue;
	
	Analysis::ChargedLepton *mu_j = *mu_itr_j;
	
	if (mu_j->family() != Muon::CALO && mu_i->GetMuon()->isStandAloneMuon() == 0) { // avoid pseudorapidity warning for SA muons [no ID track...]
	  Double_t dr = mu_i->Get4Momentum_ID()->DeltaR(*(mu_j->Get4Momentum_ID()));
	  
	  if (dr < 0.1) {
	    keepMe = kFALSE;
	  }
	}
      }
    }
    
    // SA-ST
    if (mu_i->family() != Muon::CALO && mu_i->GetMuon()->isStandAloneMuon() == 1) {
      for (std::vector<Analysis::ChargedLepton *>::iterator mu_itr_j = m_Muons.begin(); mu_itr_j != m_Muons.end(); ++mu_itr_j) {
	if (mu_itr_j == mu_itr_i) continue;
	
	Analysis::ChargedLepton *mu_j = *mu_itr_j;
	
	if (mu_j->family() != Muon::CALO && mu_i->GetMuon()->isSegmentTaggedMuon() == 1) { // avoid pseudorapidity warning for SA muons [no ID track...]
	  Double_t dr = mu_i->Get4Momentum_SA()->DeltaR(*(mu_j->Get4Momentum_ID()));
	  
	  if (dr < 0.2) {
	    keepMe = kFALSE;
	  }
	}
      }
    }
    
    i++;
    if (keepMe) {
      if (skip_muon[i]) continue;
      
      std::vector<Analysis::Jet*>::iterator jet_itr;
      for (jet_itr = m_GoodJets.begin(); jet_itr != m_GoodJets.end(); ++jet_itr) {
      	Analysis::Jet *jet = (*jet_itr);
	if (mu_i->Get4Momentum()->DeltaR(*(jet->Get4Momentum())) < 0.3) {
	  // found an jet overlapped to a jet
	  skip_muon[i] = kTRUE;
	} // overlapping muon/jet
      } // jet loop
    }
    else skip_muon[i] = kTRUE;  
  } // mu loop
  
  // fill the final tree with those jets which can be used in the analysis
  i = -1;
  for (muon_itr = m_Muons.begin(); muon_itr != m_Muons.end(); ++muon_itr) {
    i++;
    if (skip_muon[i]) continue;
    
    (*muon_itr)->set_lastcut(HllqqMuonQuality::overlap);
    m_GoodMuons.push_back((*muon_itr));
  } // jet loop
}

void HiggsllqqAnalysis::getGoodElectrons()
{
  // fills m_GoodElectrons with those electrons passing the full electron selection
  // (quality and overlap removal)
  
  // utility vector, tells if i-th electron must be skipped since overlapping with electron or muon
  std::vector<Bool_t> skip_electron;
  
  // remove overlap among electrons
  skip_electron.clear();
  
  std::vector<Analysis::ChargedLepton *>::iterator el_itr_i;
  std::vector<Analysis::ChargedLepton *>::iterator el_itr_j;
  for (el_itr_i = m_Electrons.begin(); el_itr_i != m_Electrons.end(); ++el_itr_i) {
    skip_electron.push_back(kFALSE);
  }
  
  Int_t i(-1);
  for (el_itr_i = m_Electrons.begin(); el_itr_i != m_Electrons.end(); ++el_itr_i) {
    i++;
    if (skip_electron[i]) continue;
    Analysis::ChargedLepton *el_i = (*el_itr_i);
    
    Int_t j = i;
    for (el_itr_j = el_itr_i + 1; el_itr_j != m_Electrons.end(); ++el_itr_j) {
      j++;
      if (skip_electron[j]) continue;
      Analysis::ChargedLepton *el_j = (*el_itr_j);
      
      // overlap is done at track level; we use directly D3PDReader::ElectronD3PDObjectElement variables
      // in order not to suffer from approximation errors in the building of TLorentzVectors
      if (analysis_version() == "rel_17") { // rel. 17
	if (el_j->GetElectron()->trackd0() == el_i->GetElectron()->trackd0() &&
	    el_j->GetElectron()->trackz0() == el_i->GetElectron()->trackz0() &&
	    el_j->GetElectron()->tracktheta() == el_i->GetElectron()->tracktheta() &&
	    el_j->GetElectron()->trackphi() == el_i->GetElectron()->trackphi() &&
	    el_j->GetElectron()->trackqoverp() == el_i->GetElectron()->trackqoverp()) {
	  // found two overlapped electrons, skip the lowest-Et one
	  Int_t to_be_skipped = (el_i->Get4Momentum()->Et() > el_j->Get4Momentum()->Et()) ? j : i;
	  skip_electron[to_be_skipped] = kTRUE;
	}
      } // overlapping electrons rel. 17
      else if (analysis_version() == "rel_17_2") { // rel. 17.2
	if (el_j->GetElectron()->Unrefittedtrack_d0() == el_i->GetElectron()->Unrefittedtrack_d0() &&
	    el_j->GetElectron()->Unrefittedtrack_z0() == el_i->GetElectron()->Unrefittedtrack_z0() &&
	    el_j->GetElectron()->Unrefittedtrack_theta() == el_i->GetElectron()->Unrefittedtrack_theta() &&
	    el_j->GetElectron()->Unrefittedtrack_phi() == el_i->GetElectron()->Unrefittedtrack_phi() &&
	    el_j->GetElectron()->Unrefittedtrack_qoverp() == el_i->GetElectron()->Unrefittedtrack_qoverp()) {
	  // found two overlapped electrons, skip the lowest-Et one
	  Int_t to_be_skipped = (el_i->Get4Momentum()->Et() > el_j->Get4Momentum()->Et()) ? j : i;
	  skip_electron[to_be_skipped] = kTRUE;
	} // tracks
	
	if (TMath::Abs(el_j->GetElectron()->cl_eta() - el_i->GetElectron()->cl_eta()) < 3 * 0.025 &&
	    TMath::Abs(TMath::ACos(TMath::Cos(el_j->GetElectron()->cl_phi() - el_i->GetElectron()->cl_phi()))) < 5 * 0.025) {
	  // found two overlapped electrons, skip the lowest-Et one
	  Int_t to_be_skipped = (el_i->Get4Momentum()->Et() > el_j->Get4Momentum()->Et()) ? j : i;
	  skip_electron[to_be_skipped] = kTRUE;
	} // clusters
      } // overlapping electrons rel. 17.2.2
    } // el_j loop
  } // el_i loop
  
  
  // now repeat with electron-muon overlap
  i = -1;
  for (el_itr_i = m_Electrons.begin(); el_itr_i != m_Electrons.end(); ++el_itr_i) {
    i++;
    if (skip_electron[i]) continue;
    Analysis::ChargedLepton *el = (*el_itr_i);
    
    std::vector<Analysis::ChargedLepton*>::iterator mu_itr;
    for (mu_itr = m_GoodMuons.begin(); mu_itr != m_GoodMuons.end(); ++mu_itr) {
      Analysis::ChargedLepton *mu = (*mu_itr);
      
      if (mu->GetMuon()->isStandAloneMuon() == 0) {
	if (el->Get4Momentum_ID()->DeltaR(*(mu->Get4Momentum_ID())) < 0.02) {
	  // found an electron overlapped to a muon
	  skip_electron[i] = kTRUE;
	} // overlapping electron/muon
      } // consider only non-SA muons
    } // good mu loop
  } // el loop
  
  
  // fill the final tree with those electrons which can be used in the analysis
  i = -1;
  for (el_itr_i = m_Electrons.begin(); el_itr_i != m_Electrons.end(); ++el_itr_i) {
    i++;
    if (skip_electron[i]) continue;
    
    (*el_itr_i)->set_lastcut(HllqqElectronQuality::overlap);
    m_GoodElectrons.push_back((*el_itr_i));
  } // el loop
}

void HiggsllqqAnalysis::getGoodJets()
{
  // fills m_GoodJets with those jets passing the full jet selection
  // (quality and jet-e overlap removal)
  
  // utility vector, tells if i-th jet must be skipped since overlapping with electron
  std::vector<Bool_t> skip_jet;
  
  // remove overlap wrt electrons
  skip_jet.clear();
  
  std::vector<Analysis::Jet *>::iterator jet_itr;
  for (jet_itr = m_Jets.begin(); jet_itr != m_Jets.end(); ++jet_itr) {
    skip_jet.push_back(kFALSE);
  }
  
  Int_t i(-1);
  for (jet_itr = m_Jets.begin(); jet_itr != m_Jets.end(); ++jet_itr) {
    i++;
    if (skip_jet[i]) continue;
    Analysis::Jet *jet = (*jet_itr);
    
    std::vector<Analysis::ChargedLepton*>::iterator el_itr;
    for (el_itr = m_Electrons.begin(); el_itr != m_Electrons.end(); ++el_itr) {
      Analysis::ChargedLepton *el = (*el_itr);
      if (jet->Get4Momentum()->DeltaR(*(el->Get4Momentum_ID())) < 0.4) {
	// found an jet overlapped to a electron
	skip_jet[i] = kTRUE;
      } // overlapping jet/electron
    } // el loop
  } // jet loop
  
  
  // fill the final tree with those jets which can be used in the analysis
  i = -1;
  for (jet_itr = m_Jets.begin(); jet_itr != m_Jets.end(); ++jet_itr) {
    i++;
    if (skip_jet[i]) continue;
    
    (*jet_itr)->set_lastcut(HllqqJetQuality::overlap);
    m_GoodJets.push_back((*jet_itr));
  } // jet loop
}

void HiggsllqqAnalysis::getDileptons()
{
  // fills m_Dileptons with those dileptons which can be used to build a Quadrilepton
  // removed charge requirement to enable filling quadruplets with same sign leptons for background studies
  
  std::vector<Analysis::ChargedLepton *>::iterator mu_itr_i;
  std::vector<Analysis::ChargedLepton *>::iterator mu_itr_j;
  for (mu_itr_i = m_GoodMuons.begin(); mu_itr_i != m_GoodMuons.end(); ++mu_itr_i) {
    for (mu_itr_j = mu_itr_i + 1; mu_itr_j != m_GoodMuons.end(); ++mu_itr_j) {
      if ((*mu_itr_i)->charge() * (*mu_itr_j)->charge() < 0.) {
	if ((*mu_itr_i)->charge() > 0) // dilepton expects (positive, negative)
	  m_Dileptons.push_back(new Analysis::Dilepton(*mu_itr_i, *mu_itr_j));
	else
	  m_Dileptons.push_back(new Analysis::Dilepton(*mu_itr_j, *mu_itr_i));
      } // opposite charge
      else
	m_Dileptons.push_back(new Analysis::Dilepton(*mu_itr_i, *mu_itr_j));
    } // mu_itr_j
  } // mu_itr_i
  
  std::vector<Analysis::ChargedLepton *>::iterator el_itr_i;
  std::vector<Analysis::ChargedLepton *>::iterator el_itr_j;
  for (el_itr_i = m_GoodElectrons.begin(); el_itr_i != m_GoodElectrons.end(); ++el_itr_i) {
    for (el_itr_j = el_itr_i + 1; el_itr_j != m_GoodElectrons.end(); ++el_itr_j) {
      if ((*el_itr_i)->charge() * (*el_itr_j)->charge() < 0.) {
	if ((*el_itr_i)->charge() > 0) // dilepton expects (positive, negative)
	  m_Dileptons.push_back(new Analysis::Dilepton(*el_itr_i, *el_itr_j));
	else
	  m_Dileptons.push_back(new Analysis::Dilepton(*el_itr_j, *el_itr_i));
      } // opposite charge
      else
	m_Dileptons.push_back(new Analysis::Dilepton(*el_itr_i, *el_itr_j));
    } // el_itr_j
  } // el_itr_i
}

void HiggsllqqAnalysis::getQuadrileptons()
{
  // fill m_Quadrileptons with those quadrileptions build out of m_Dileptons, ordering
  // the two Z candidates according to certain criteria (default analysis: Z1 is the closest to
  // PDG mass, Z2 the other)
  // isGood() is called for each new quadrilepton, to assign to Quadrilepton::lastcut() the ID
  // of the last quality cut passed
  
  std::vector<Analysis::Dilepton *>::iterator Z_itr_i;
  std::vector<Analysis::Dilepton *>::iterator Z_itr_j;
  for (Z_itr_i = m_Dileptons.begin(); Z_itr_i != m_Dileptons.end(); ++Z_itr_i) {
    for (Z_itr_j = Z_itr_i + 1; Z_itr_j != m_Dileptons.end(); ++Z_itr_j) {
      if ((*Z_itr_j)->OverlapsWith(*Z_itr_i) == kFALSE) { // Z1 and Z2 must not share a lepton
	Double_t mass_delta_i = TMath::Abs(Z_pdg_mass - (*Z_itr_i)->Get4Momentum()->M());
	Double_t mass_delta_j = TMath::Abs(Z_pdg_mass - (*Z_itr_j)->Get4Momentum()->M());
	
	Analysis::Quadrilepton *higgs(0);
	if ((*Z_itr_i)->IsNeutral() && (*Z_itr_j)->IsNeutral()) { // Standard Analysis, both Opposite Charged Pairs
	  if (mass_delta_i < mass_delta_j)
	    higgs = new Analysis::Quadrilepton(*Z_itr_i, *Z_itr_j);
	  else
	    higgs = new Analysis::Quadrilepton(*Z_itr_j, *Z_itr_i);
	} else { // Background studies, Z2 can have Same Sign Leptons
	  if ((*Z_itr_i)->IsNeutral())
	    higgs = new Analysis::Quadrilepton(*Z_itr_i, *Z_itr_j);
	  else if ((*Z_itr_j)->IsNeutral())
	    higgs = new Analysis::Quadrilepton(*Z_itr_j, *Z_itr_i);
	}
	if (higgs != 0) {
	  isGood(higgs);
	  m_Quadrileptons.push_back(higgs);
	}
      } // non-overlapping Z candidates
    } // Z_itr_j
  } // Z_itr_i
}

Analysis::Quadrilepton *HiggsllqqAnalysis::findClosestToZZ()
{
  // find the quadrilepton, among those passing all the cuts up to the DeltaR one included, with
  //  - Z1 mass closest to Z_pdg_mass
  //  - then, if needed, highest Z2 mass
  // and return a pointer to it (return 0 if none)
  // this function is used to ensure that up to 1 quadrilepton in the event can be selected
  
  Analysis::Quadrilepton *best(0);
  
  std::vector<Analysis::Quadrilepton *>::iterator H_itr;
  for (H_itr = m_Quadrileptons.begin(); H_itr != m_Quadrileptons.end(); ++H_itr) {
    if ((*H_itr)->lastcut() >= HllqqQuadrileptonQuality::DeltaR) {
      
      // first quadrilepton is always the best :)
      if (best == 0) {
	best = *H_itr;
      } else {
	Double_t mass_delta_this = TMath::Abs(Z_pdg_mass - (*H_itr)->GetZ1()->Get4Momentum()->M());
	Double_t mass_delta_best = TMath::Abs(Z_pdg_mass - best->GetZ1()->Get4Momentum()->M());
	
	if (mass_delta_this < mass_delta_best) {
	  best = (*H_itr);
	} else if (mass_delta_this == mass_delta_best) {
	  if ((*H_itr)->GetZ2()->Get4Momentum()->M() > best->GetZ2()->Get4Momentum()->M()) {
	    best = (*H_itr);
	  }
	}
      }
    } // quadrilepton passes at least the deltaR cut
  }
  
  return best;
}

Bool_t HiggsllqqAnalysis::isClosestToZZ(Analysis::Quadrilepton *higgs)
{
  return (higgs == findClosestToZZ());
}

Bool_t HiggsllqqAnalysis::isSelected(Analysis::Quadrilepton *higgs)
{
  // check if this quadrilepton is among the GoodQuadrileptons
  
  std::vector<Analysis::Quadrilepton *>::iterator sel_itr;
  for (sel_itr = m_GoodQuadrileptons.begin(); sel_itr != m_GoodQuadrileptons.end(); ++sel_itr) {
    if (higgs == (*sel_itr)) {
      return kTRUE;
    }
  }
  
  return kFALSE;
}

void HiggsllqqAnalysis::getGoodQuadrileptons()
{
  // fill m_GoodQuadrileptons with the quadrileptons passing the full selection (Higgs candidates)
  // official analysis will fill this vector up to 1 times per event
  
  // find the quadrilepton with Z1.M() closest to Z_pdg_mass (and highest Z2.M())
  Analysis::Quadrilepton *best = findClosestToZZ();
  
  // the selected quadruplet must of course pass the final cuts
  if (best != 0) {
    //if (best->lastcut() >= HllqqQuadrileptonQuality::mass_window) { // kostas
    if (best->lastcut() >= HllqqQuadrileptonQuality::d0_sig) {
      m_GoodQuadrileptons.push_back(best);
    }
  }
}

void HiggsllqqAnalysis::getGoodLeptons()
{
  // this function aims at collecting leptons to be used in the analysis
  
  // must be called once per event
  if (m_called_getGoodLeptons) return;
  else m_called_getGoodLeptons = kTRUE;
  
  /////////
  //  Preliminary step: collect leptons passing one-lepton cuts from the right branches
  //  (an one-lepton cut is a cut which is applied to a lepton regardless of the other leptons)
  /////////
  
  // add muons
  D3PDReader::MuonD3PDObject *mu_branch(0);
  if (getMuonFamily() == Muon::STACO) mu_branch = &(ntuple->mu_staco);
  else if (getMuonFamily() == Muon::MUID) mu_branch = &(ntuple->mu_muid);
  else Abort("Unknown muon family requested");
  
  getMuons(mu_branch, getMuonFamily());
  if (DoCaloMuons)
    getMuons(&(ntuple->mu_calo), Muon::CALO);
  
  
  // add electrons
  D3PDReader::ElectronD3PDObject *el_branch(0);
  if (getElectronFamily() == Electron::GSF) el_branch = &(ntuple->el_GSF);
  else if (getElectronFamily() == Electron::noGSF) el_branch = &(ntuple->el);
  else Abort("Unknown electron family requested");
  
  getElectrons(el_branch, getElectronFamily());
  
  
  // add jets     
  /*Include Jet Families???
    a)jet_antikt4truth 
    b)jet_akt4topoem
  */
  D3PDReader::JetD3PDObject *jet_branch(0);
  jet_branch = &(ntuple->jet_akt4topoem);
  getJets(jet_branch);
  
  
  ////////
  // Part 1: remove overlap among electrons and between muons and jets
  // To be implemented!!!
  ////////
  
  getGoodJets();
  
  /////////
  //  Part 2: remove overlap among muons
  /////////
  
  getGoodMuons();
  
  ////////
  // Part 3: remove overlap among electrons and between muons and electrons
  ////////
  
  getGoodElectrons();  
}

void HiggsllqqAnalysis::getGoodObjects()
{
  // the aim of this function is to build dileptons and quadrileptons for the analysis
  // using leptons selected by getGoodLeptons()
  
  // must be called once per event
  if (m_called_getGoodObjects) return;
  else m_called_getGoodObjects = kTRUE;
  
  ////////
  // Part 3: build dileptons
  ////////
  
  getDileptons();
  
  
  ////////
  // Part 5: build quadrileptons
  ////////
  
  getQuadrileptons();
  
  
  ////////
  // Part 6: find good quadrileptons (0 or 1, in the official analysis)
  ////////
  
  getGoodQuadrileptons();
  
  
  ////////
  // now we have filled
  //  - m_Muons, m_Electrons with pointers to muons and electrons passing selection cuts
  //  - m_GoodMuons, m_GoodElectrons with pointers to muons and electrons passing selection cuts and overlap removal
  //  - m_Dileptons with pointers to dileptons passing selection criteria (neutral charge in default analysis)
  //  - m_Quadrileptons with pointers to quadrileptons build from m_Dileptons
  //  - m_GoodQuadrileptons with pointers to those quadrileptons selected by the full analysis (up to 1 element in default analysis)
  // to be sure we are dealing correctly with object deletion, ALL and ONLY the objects in
  //      m_Muons      m_Electrons    m_Dileptons    m_Quadrileptons
  // must be deleted at the end of execute_analysis
  ////////
}

Bool_t HiggsllqqAnalysis::GetGoodObjects()
{
  /*
    See the HiggsqqllAnalysis Code for more information!!!    
  */
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::isMCPMuon(Analysis::ChargedLepton *lep)
{
  if (lep->flavor() != Analysis::ChargedLepton::MUON) {
    Error("isMCPMuon", "Calling MCP quality for a non-muon flavor (id=%d)...", lep->flavor());
    return kFALSE;
  }
  
  D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
  
  Double_t eta = TMath::Abs(mu->eta()); // shouldn't we use id_eta instead?
  
  if (mu->isStandAloneMuon()) {
    if (eta > 2.5) {
      Bool_t good_hits = (mu->nCSCEtaHits() + mu->nCSCPhiHits() > 0 && mu->nMDTEMHits() > 0 && mu->nMDTEOHits() > 0);
      if (!good_hits) return kFALSE;
    } else {
      return kFALSE;
    }
  } else {
    Bool_t blayer = (!mu->expectBLayerHit() || mu->nBLHits() > 0);
    if (!blayer) return kFALSE;
    if (analysis_version() == "rel_17") {
      Bool_t pixhits = (mu->nPixHits() + mu->nPixelDeadSensors() > 1);
      if (!pixhits) return kFALSE;
      Bool_t scthits = (mu->nSCTHits() + mu->nSCTDeadSensors() >= 6);
      if (!scthits) return kFALSE;
    } else if (analysis_version() == "rel_17_2") {
      Bool_t pixhits = (mu->nPixHits() + mu->nPixelDeadSensors() > 0);
      if (!pixhits) return kFALSE;
      Bool_t scthits = (mu->nSCTHits() + mu->nSCTDeadSensors() >= 5);
      if (!scthits) return kFALSE;
    }
    Bool_t holes = (mu->nPixHoles() + mu->nSCTHoles() < 3);
    if (!holes) return kFALSE;
    
    Int_t n = mu->nTRTHits() + mu->nTRTOutliers();
    
    if (lep->family() != Muon::CALO) {
      Bool_t case1 = ((n > 5) && ((Double_t) mu->nTRTOutliers()) < 0.9 * (Double_t)n);
      Bool_t case2 = (n > 5) ? ((Double_t) mu->nTRTOutliers() < 0.9 * (Double_t)n) : kTRUE;
      
      if (analysis_version() == "rel_17") {
	Bool_t trt_ext = (eta < 1.9) ? case1 : case2;
	if (!trt_ext) return kFALSE;
      } else if (analysis_version() == "rel_17_2") {
	Bool_t trt_ext = (0.1 < eta && eta < 1.9) ? case1 : case2;
	if (!trt_ext) return kFALSE;
      }
    } else {
      Bool_t good_trt = (eta < 0.1 && ((n < 6 || mu->nTRTOutliers() < 0.9 * (Double_t)n)));
      if (!good_trt) return kFALSE;
    }
  }
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::isGood(Analysis::ChargedLepton *lep)
{
  lep->set_lastcut(-1);
  
  //// Muons
  if (lep->flavor() == Analysis::ChargedLepton::MUON) {
    D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
    
    if (lep->family() == getMuonFamily() || (lep->family() == Muon::CALO && DoCaloMuons)) lep->set_lastcut(HllqqMuonQuality::family);
    else return kFALSE;
    if ((lep->family() == Muon::MUID && mu->tight() == 1) ||
	(lep->family() == Muon::STACO && (mu->author() == 6 || mu->author() == 7)) ||
	(lep->family() == Muon::CALO && DoCaloMuons && mu->author() == 16 && (mu->caloMuonIdTag() > 10 || mu->caloLRLikelihood() > 0.9))) lep->set_lastcut(HllqqMuonQuality::quality);
    else return kFALSE;
    if (TMath::Abs(mu->d0_exPV()) < 1. || mu->isStandAloneMuon()) lep->set_lastcut(HllqqMuonQuality::cosmic);
    else return kFALSE;
    if ((lep->family() != Muon::CALO && TMath::Abs(mu->eta()) < 2.5/*2.7*/) || (TMath::Abs(mu->eta()) < 0.1 && DoCaloMuons)) lep->set_lastcut(HllqqMuonQuality::eta);
    else return kFALSE;
    if ((lep->family() != Muon::CALO && lep->Get4Momentum()->Pt() > 7000/*6000*/) ||
	(lep->family() == Muon::CALO && DoCaloMuons && lep->Get4Momentum()->Pt() > 15000)) lep->set_lastcut(HllqqMuonQuality::pt);
    else return kFALSE;
    if (isMCPMuon(lep)) lep->set_lastcut(HllqqMuonQuality::MCP);
    else return kFALSE;
    if (TMath::Abs(mu->z0_exPV()) < 10. || mu->isStandAloneMuon()) lep->set_lastcut(HllqqMuonQuality::z0);
    else return kFALSE;
  }
  //// Electrons
  else if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    if (lep->family() == getElectronFamily()) lep->set_lastcut(HllqqElectronQuality::family);
    else return kFALSE;
    if (el->author() == 1 || el->author() == 3) lep->set_lastcut(HllqqElectronQuality::algorithm);
    else return kFALSE;
    if (analysis_version() == "rel_17") { // rel. 17
      if (el->tightPP() == 1) lep->set_lastcut(HllqqElectronQuality::quality);
      /*
      if (isLoosePlusPlusH4l(el->etas2(),
			     el->cl_E() / TMath::CosH(el->etas2()),
			     el->Ethad() / (el->cl_E() / TMath::CosH(el->etas2())),
			     el->Ethad1() / (el->cl_E() / TMath::CosH(el->etas2())),
			     el->reta(),
			     el->weta2(),
			     el->f1(),
			     el->wstot(),
			     (el->emaxs1() - el->Emax2()) / (el->emaxs1() + el->Emax2()),
			     el->deltaeta1(),
			     el->nSiHits(),
			     el->nSCTOutliers() + el->nPixelOutliers(),
			     el->nPixHits(),
			     el->nPixelOutliers(),
			     false,
			     false)) lep->set_lastcut(HllqqElectronQuality::quality);
      */
      else return kFALSE;
    } // rel. 17
    else if (analysis_version() == "rel_17_2") { // rel. 17.2
      if (el->tightPP() == 1) lep->set_lastcut(HllqqElectronQuality::quality);
      /*
      if (passMultiLepton(el->etas2(),
			  el->cl_E() / TMath::CosH(el->etas2()),
			  el->Ethad() / (el->cl_E() / TMath::CosH(el->etas2())),
			  el->Ethad1() / (el->cl_E() / TMath::CosH(el->etas2())),
			  el->reta(),
			  el->weta2(),
			  el->f1(),
			  el->f3(),
			  el->wstot(),
			  (el->emaxs1() - el->Emax2()) / (el->emaxs1() + el->Emax2()),
			  el->deltaeta1(),
			  el->nSiHits(),
			  el->nSCTDeadSensors() + el->nPixelDeadSensors(),
			  el->nPixHits(),
			  el->nPixelDeadSensors(),
			  el->deltaphiRescaled(),
			  CommonTools::getBremFitDp(el),
			  (el->nTRTHits() + el->nTRTOutliers() > 0) ? ((Double_t)el->nTRTHighTHits() + (Double_t)el->nTRTHighTOutliers()) / ((Double_t)el->nTRTHits() + (Double_t)el->nTRTOutliers()) : 0,
			  el->nTRTHits() + el->nTRTOutliers(),
			  el->nBLHits(),
			  el->expectBLayerHit(),
			  false)) lep->set_lastcut(HllqqElectronQuality::quality);
      */
      else return kFALSE;
    } // rel. 17.2
    if (TMath::Abs(el->cl_eta()) < 2.47) lep->set_lastcut(HllqqElectronQuality::eta);
    else return kFALSE;
    if (lep->Get4Momentum()->Et() > 7000) lep->set_lastcut(HllqqElectronQuality::Et);
    else return kFALSE;
    if ((el->OQ() & 1446) == 0) lep->set_lastcut(HllqqElectronQuality::objectquality);
    else return kFALSE;
    if (TMath::Abs(lep->z0()) < 10.) lep->set_lastcut(HllqqElectronQuality::z0);
    else return kFALSE;
    
  } else {
    Error("isGood", "Unknown lepton flavour %d", lep->flavor());
    return kFALSE;
  }
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::IsGoodMuon(Analysis::ChargedLepton *lep)
{
  /*
    family,
    quality,
    cosmic,
    eta,
    pt,
    MCP,
    z0,
    overlap,
  */

  lep->set_lastcut(-1);
  
  //// Muons
  if (lep->flavor() == Analysis::ChargedLepton::MUON) {
    D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
    Bool_t dolowmass = GetDoLowMass();
    if(cut_leptons) Muon0++;
    
    if (lep->family() == getMuonFamily() || (lep->family() == Muon::CALO && DoCaloMuons)) lep->set_lastcut(HllqqMuonQuality::family); 
    else return kFALSE;
    if(cut_leptons) Muon1++;
    
    if (((lep->family() != Muon::CALO && TMath::Abs(mu->eta()) < 2.5) || (TMath::Abs(mu->eta()) < 0.1 && DoCaloMuons)) && 
	((lep->family() != Muon::CALO && lep->Get4Momentum()->Pt() > 7000 && dolowmass)  || 
	 (lep->family() == Muon::CALO && lep->Get4Momentum()->Pt() > 15000 && dolowmass && DoCaloMuons) ||
	 (lep->family() != Muon::CALO && lep->Get4Momentum()->Pt() > 20000 && !dolowmass))) lep->set_lastcut(HllqqMuonQuality::quality/*Kinematics*/);
    else return kFALSE;
    if(cut_leptons) Muon2++;
    
    if ((lep->family() == Muon::MUID && mu->tight() == 1) ||
	(lep->family() == Muon::STACO && (mu->author() == 6 || mu->author() == 7)&& (mu->isCombinedMuon()==1 || mu->isSegmentTaggedMuon()==1)) ||
	(lep->family() == Muon::CALO && DoCaloMuons && mu->author() == 16 && (mu->caloMuonIdTag() > 10 || mu->caloLRLikelihood() > 0.9))) lep->set_lastcut(HllqqMuonQuality::cosmic/*author*/);
    else return kFALSE;
    if(cut_leptons) Muon3++;
    
    if (isMCPMuon(lep)) lep->set_lastcut(HllqqMuonQuality::eta/*MCP*/);
    else return kFALSE;
    if(cut_leptons) Muon4++;
    
    if ((TMath::Abs(mu->d0_exPV()) < 1. && TMath::Abs(mu->z0_exPV()) < 10.) || mu->isStandAloneMuon()) lep->set_lastcut(HllqqMuonQuality::pt/*cosmic*/);
    else return kFALSE;
    if(cut_leptons) Muon5++;
    
    if(!GetDoQCDSelection()) {
      
      //d0 Significance
      if (TMath::Abs(lep->d0() / lep->d0_sig()) < 3.5) lep->set_lastcut(HllqqMuonQuality::MCP/*d0Sig*/);
      else return kFALSE;
      if(cut_leptons) Muon6++;
      
      if (((mu->ptcone20()/mu->pt())<.1 && !dolowmass && mu->pt() > 20000.) || (dolowmass && (mu->ptcone20()/mu->pt())<.1)) lep->set_lastcut(HllqqMuonQuality::z0/*Isolation*/);
      else return kFALSE;
      if(cut_leptons) Muon7++;
    }
    return kTRUE;
  }
  else{
    cout<<" This is not a Muon! "<<endl;
    return kFALSE;
  }
}

Bool_t HiggsllqqAnalysis::IsGoodElectron(Analysis::ChargedLepton *lep)
{
  /*
    family,
    algorithm,
    quality,
    eta,
    Et,
    objectquality,
    z0,
    overlap,
  */
  
  if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    Bool_t dolowmass = GetDoLowMass();
    float lepton_et =0.;
    
    if(dolowmass)
      lepton_et = 7000.;
    else
      lepton_et = 20000.;
    
    if (lep->family() == getElectronFamily()) lep->set_lastcut(HllqqElectronQuality::family);
    else return kFALSE;
    if(cut_leptons) Electron0++;
    
    if (el->goodOQ()==1) lep->set_lastcut(HllqqElectronQuality::algorithm);//obj_quality
    else return kFALSE;
    if(cut_leptons) Electron1++;
    
    if((el->author() == 1 || el->author() == 3)) lep->set_lastcut(HllqqElectronQuality::quality);//family
    else return kFALSE;
    if(cut_leptons) Electron2++;
    
    if(!GetDoQCDSelection()) {
      if ((el->mediumPP() == 1 && !dolowmass) || (el->tightPP() == 1 && dolowmass)) lep->set_lastcut(HllqqElectronQuality::eta);//tight
      else return kFALSE;
    }
    else {
      if (el->loosePP() == 1) lep->set_lastcut(HllqqElectronQuality::eta);//tight
      else return kFALSE;
    }
    if(cut_leptons) Electron3++;
    
    //if (lep->Get4Momentum()->Et() > lepton_et) lep->set_lastcut(HllqqElectronQuality::);//Et
    //else return kFALSE;
    
    if (lep->Get4Momentum()->Et() > lepton_et && TMath::Abs(el->cl_eta()) < 2.47) lep->set_lastcut(HllqqElectronQuality::Et);//eta
    else return kFALSE;
    if(cut_leptons) Electron4++;
    
    if(!GetDoQCDSelection()) {
      
      if ((el->ptcone20()/el->pt()) < 0.1) lep->set_lastcut(HllqqElectronQuality::objectquality);//isolation
      else return kFALSE;
      if(cut_leptons) Electron5++;
      
      if((TMath::Abs(lep->d0()) / lep->d0_sig())<6.5) lep->set_lastcut(HllqqElectronQuality::z0);//d0Sig
      else return kFALSE;
      if(cut_leptons) Electron6++;
    }
    return kTRUE;
  }
  else{
    cout<<" This is not an Electron! "<<endl;
    return kFALSE;
  }
}

Bool_t HiggsllqqAnalysis::isGoodJet(Analysis::Jet *jet)
{
  // good jet definition (AntiKt4H1TopoJets)
  D3PDReader::JetD3PDObjectElement *Jet = jet->GetJet();
  
  Bool_t dolowmass = GetDoLowMass();
  
  if(/*!isBadLooser(jet)*/Jet->isBadLooseMinus() == 0) jet->set_lastcut(HllqqJetQuality::jetCleaning);
  else return kFALSE;
  
  if(dolowmass)
    if (jet->rightpt()/1000. > 20. && TMath::Abs(jet->righteta()) < 2.5 && jet->rightE()>0) jet->set_lastcut(HllqqJetQuality::kinematics);
    else return kFALSE;
  else
    if(jet->rightpt()/1000. > 25. && TMath::Abs(jet->righteta()) < 4.5 && jet->rightE()>0) jet->set_lastcut(HllqqJetQuality::kinematics);
    else return kFALSE;
  
  if(TMath::Abs(Jet->jvtxf()) > 0.75) jet->set_lastcut(HllqqJetQuality::Pileup);
  else return kFALSE;
  
  return kTRUE; 
}

Bool_t HiggsllqqAnalysis::isGood(Analysis::Quadrilepton *higgs)
{
  higgs->set_lastcut(-1);
  
  higgs->set_lastcut(HllqqQuadrileptonQuality::mass_2e2mu);
  
  if (isSFOS(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::opposite_charge);
  else return kFALSE;
  
  if (passesPtLeptonsCut(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::kinematics);
  else return kFALSE;
  
  if (isTriggerMatched(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::trigmatch);
  else return kFALSE;
  
  if (passesZ1MassCut(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::Z1_mass);
  else return kFALSE;
  
  if (passesZ2MassCut(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::Z2_mass);
  else return kFALSE;
  
  if (passesDeltaRCut(higgs) && passesJpsiVeto(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::DeltaR);
  else return kFALSE;
  
  if (passesTrackIsolationCut(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::track_iso);
  else return kFALSE;
  
  if (passesCaloIsolationCut(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::calo_iso);
  else return kFALSE;
  
  if (hasGoodImpactParameterSignificance(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::d0_sig);
  else return kFALSE;
  
  // higgs->set_lastcut(HllqqQuadrileptonQuality::mass_window); // kostas
  
  // if (passesTrigger()) higgs->set_lastcut(HllqqQuadrileptonQuality::trigger_on_top); // kostas
  // else return kFALSE;
  
  // if (isTriggerMatched(higgs)) higgs->set_lastcut(HllqqQuadrileptonQuality::trigmatch_on_top); // kostas
  // else return kFALSE;
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::isSFOS(Analysis::Quadrilepton * higgs)
{
  Bool_t charge = (higgs->IsNeutral());
  
  Int_t n_CALO_or_SA(0);
  
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  for (UInt_t i = 0; i < leptons.size(); i++) {
    if (leptons[i]->flavor() == Analysis::ChargedLepton::MUON) {
      if (leptons[i]->family() == Muon::CALO || leptons[i]->GetMuon()->isStandAloneMuon()) {
	n_CALO_or_SA++;
      }
    }
  }
  
  return (charge && n_CALO_or_SA <= 1);
}

Bool_t HiggsllqqAnalysis::passesPtLeptonsCut(Analysis::Quadrilepton * higgs)
{
  Int_t highpt_leptons(0);
  Int_t mediumpt_leptons(0);
  Int_t lowpt_leptons(0);
  
  Double_t highpt_threshold   = 20000;
  Double_t mediumpt_threshold = 15000;
  Double_t lowpt_threshold    = 10000;
  
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  for (UInt_t i = 0; i < leptons.size(); i++) {
    if (leptons[i]->Get4Momentum()->Pt() > highpt_threshold) highpt_leptons++;
    if (leptons[i]->Get4Momentum()->Pt() > mediumpt_threshold) mediumpt_leptons++;
    if (leptons[i]->Get4Momentum()->Pt() > lowpt_threshold) lowpt_leptons++;
  }
  
  if (highpt_leptons >= 1 && mediumpt_leptons >= 2 && lowpt_leptons >= 3) return kTRUE;
  return kFALSE;
}

Bool_t HiggsllqqAnalysis::isTriggerMatched(Analysis::Quadrilepton * higgs)
{
  // set the trigger navigation variables
  
  // we need non-const access to the variables
  D3PDReader::HSG2EventReader& ncevent = const_cast< D3PDReader::HSG2EventReader& >(*ntuple);
  
  // set the main navigation variables
  UInt_t smk = static_cast< UInt_t >(m_TrigDecisionToolD3PD->GetSMK());
  m_trigNavVar.set_trig_DB_SMK(smk);
  Int_t trig_Nav_n = ntuple->trig_Nav.n();
  m_trigNavVar.set_trig_Nav_n(trig_Nav_n);
  std::vector< short >* trig_Nav_chain_ChainId = ncevent.trig_Nav.chain_ChainId();
  m_trigNavVar.set_trig_Nav_chain_ChainId(trig_Nav_chain_ChainId);
  std::vector< std::vector< int > >* trig_Nav_chain_RoIType = ncevent.trig_Nav.chain_RoIType();
  m_trigNavVar.set_trig_Nav_chain_RoIType(trig_Nav_chain_RoIType);
  std::vector< std::vector< int > >* trig_Nav_chain_RoIIndex = ncevent.trig_Nav.chain_RoIIndex();
  m_trigNavVar.set_trig_Nav_chain_RoIIndex(trig_Nav_chain_RoIIndex);
  
  // set the electron variables
  std::vector< std::vector< int > >* trig_RoI_EF_e_egammaContainer_egamma_Electrons = ncevent.trig_RoI_EF_e.egammaContainer_egamma_Electrons();
  std::vector< std::vector< int > >* trig_RoI_EF_e_egammaContainer_egamma_ElectronsStatus = ncevent.trig_RoI_EF_e.egammaContainer_egamma_ElectronsStatus();
  Int_t trig_EF_el_n = ntuple->trig_EF_el.n();
  std::vector< float >* trig_EF_el_eta = ncevent.trig_EF_el.eta();
  std::vector< float >* trig_EF_el_phi = ncevent.trig_EF_el.phi();
  
  m_trigNavVar.set_trig_RoI_EF_e_egammaContainer_egamma_Electrons(trig_RoI_EF_e_egammaContainer_egamma_Electrons);
  m_trigNavVar.set_trig_RoI_EF_e_egammaContainer_egamma_ElectronsStatus(trig_RoI_EF_e_egammaContainer_egamma_ElectronsStatus);
  m_trigNavVar.set_trig_EF_el_n(trig_EF_el_n);
  m_trigNavVar.set_trig_EF_el_eta(trig_EF_el_eta);
  m_trigNavVar.set_trig_EF_el_phi(trig_EF_el_phi);
  
  // set the muon variables:
  std::vector< int >* trig_RoI_EF_mu_Muon_ROI = ncevent.trig_RoI_EF_mu.Muon_ROI();
  std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFInfoContainer = 0;
  if (analysis_version() == "rel_17_2")trig_RoI_EF_mu_TrigMuonEFInfoContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer();
  std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = 0;
  if (analysis_version() == "rel_17_2")trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainerStatus();
  std::vector< int >* trig_RoI_L2_mu_CombinedMuonFeature = ncevent.trig_RoI_L2_mu.CombinedMuonFeature();
  std::vector< int >* trig_RoI_L2_mu_CombinedMuonFeatureStatus = ncevent.trig_RoI_L2_mu.CombinedMuonFeatureStatus();
  std::vector< int >* trig_RoI_L2_mu_MuonFeature = ncevent.trig_RoI_L2_mu.MuonFeature();
  std::vector< int >* trig_RoI_L2_mu_Muon_ROI = ncevent.trig_RoI_L2_mu.Muon_ROI();
  std::vector< std::vector< float > >* trig_EF_trigmuonef_track_CB_pt = ncevent.trig_EF_trigmuonef.track_CB_pt();
  std::vector< std::vector< float > >* trig_EF_trigmuonef_track_CB_eta = ncevent.trig_EF_trigmuonef.track_CB_eta();
  std::vector< std::vector< float > >* trig_EF_trigmuonef_track_CB_phi = ncevent.trig_EF_trigmuonef.track_CB_phi();
  std::vector< std::vector< float > >* trig_EF_trigmuonef_track_SA_pt = ncevent.trig_EF_trigmuonef.track_SA_pt();
  std::vector< std::vector< float > >* trig_EF_trigmuonef_track_SA_eta = ncevent.trig_EF_trigmuonef.track_SA_eta();
  std::vector< std::vector< float > >* trig_EF_trigmuonef_track_SA_phi = ncevent.trig_EF_trigmuonef.track_SA_phi();
  std::vector< std::vector< float > >* trig_EF_trigmugirl_track_CB_pt = ncevent.trig_EF_trigmugirl.track_CB_pt();
  std::vector< std::vector< float > >* trig_EF_trigmugirl_track_CB_eta = ncevent.trig_EF_trigmugirl.track_CB_eta();
  std::vector< std::vector< float > >* trig_EF_trigmugirl_track_CB_phi = ncevent.trig_EF_trigmugirl.track_CB_phi();
  std::vector< float >* trig_L2_combmuonfeature_eta = ncevent.trig_L2_combmuonfeature.eta();
  std::vector< float >* trig_L2_combmuonfeature_phi = ncevent.trig_L2_combmuonfeature.phi();
  std::vector< float >* trig_L2_muonfeature_eta = ncevent.trig_L2_muonfeature.eta();
  std::vector< float >* trig_L2_muonfeature_phi = ncevent.trig_L2_muonfeature.phi();
  std::vector< float >* trig_L1_mu_eta = ncevent.trig_L1_mu.eta();
  std::vector< float >* trig_L1_mu_phi = ncevent.trig_L1_mu.phi();
  std::vector< std::string >* trig_L1_mu_thrName = ncevent.trig_L1_mu.thrName();
  
  std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFIsolationContainer = 0;
  if (ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainer.IsAvailable()) {
    trig_RoI_EF_mu_TrigMuonEFIsolationContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainer();
  } else {
    //Warning("isTriggerMatched", "unable to set TrigMuonEFIsolationContainer");
  }
  std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus = 0;
  if (ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainerStatus.IsAvailable()) {
    trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainerStatus();
  } else {
    //Warning("isTriggerMatched", "unable to set TrigMuonEFIsolationContainerStatus");
  }
  
  std::vector< int >* trig_EF_trigmuonef_EF_mu24i_tight = 0;
  if (ncevent.trig_EF_trigmuonef.EF_mu24i_tight.IsAvailable()) {
    trig_EF_trigmuonef_EF_mu24i_tight = ncevent.trig_EF_trigmuonef.EF_mu24i_tight();
  } else {
    //Warning("isTriggerMatched", "unable to set EF_mu24i_tight");
  }
  std::vector< int >* trig_EF_trigmuonef_EF_mu36_tight = 0;
  if (ncevent.trig_EF_trigmuonef.EF_mu36_tight.IsAvailable()) {
    trig_EF_trigmuonef_EF_mu36_tight = ncevent.trig_EF_trigmuonef.EF_mu36_tight();
  } else {
    //Warning("isTriggerMatched", "unable to set EF_mu36_tight");
  }
  
  m_trigNavVar.set_trig_RoI_EF_mu_Muon_ROI(trig_RoI_EF_mu_Muon_ROI);
  m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(trig_RoI_EF_mu_TrigMuonEFInfoContainer);
  m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus);
  m_trigNavVar.set_trig_RoI_L2_mu_CombinedMuonFeature(trig_RoI_L2_mu_CombinedMuonFeature);
  m_trigNavVar.set_trig_RoI_L2_mu_CombinedMuonFeatureStatus(trig_RoI_L2_mu_CombinedMuonFeatureStatus);
  m_trigNavVar.set_trig_RoI_L2_mu_MuonFeature(trig_RoI_L2_mu_MuonFeature);
  m_trigNavVar.set_trig_RoI_L2_mu_Muon_ROI(trig_RoI_L2_mu_Muon_ROI);
  m_trigNavVar.set_trig_EF_trigmuonef_track_CB_pt(trig_EF_trigmuonef_track_CB_pt);
  m_trigNavVar.set_trig_EF_trigmuonef_track_CB_eta(trig_EF_trigmuonef_track_CB_eta);
  m_trigNavVar.set_trig_EF_trigmuonef_track_CB_phi(trig_EF_trigmuonef_track_CB_phi);
  m_trigNavVar.set_trig_EF_trigmuonef_track_SA_pt(trig_EF_trigmuonef_track_SA_pt);
  m_trigNavVar.set_trig_EF_trigmuonef_track_SA_eta(trig_EF_trigmuonef_track_SA_eta);
  m_trigNavVar.set_trig_EF_trigmuonef_track_SA_phi(trig_EF_trigmuonef_track_SA_phi);
  m_trigNavVar.set_trig_EF_trigmugirl_track_CB_pt(trig_EF_trigmugirl_track_CB_pt);
  m_trigNavVar.set_trig_EF_trigmugirl_track_CB_eta(trig_EF_trigmugirl_track_CB_eta);
  m_trigNavVar.set_trig_EF_trigmugirl_track_CB_phi(trig_EF_trigmugirl_track_CB_phi);
  m_trigNavVar.set_trig_L2_combmuonfeature_eta(trig_L2_combmuonfeature_eta);
  m_trigNavVar.set_trig_L2_combmuonfeature_phi(trig_L2_combmuonfeature_phi);
  m_trigNavVar.set_trig_L2_muonfeature_eta(trig_L2_muonfeature_eta);
  m_trigNavVar.set_trig_L2_muonfeature_phi(trig_L2_muonfeature_phi);
  m_trigNavVar.set_trig_L1_mu_eta(trig_L1_mu_eta);
  m_trigNavVar.set_trig_L1_mu_phi(trig_L1_mu_phi);
  m_trigNavVar.set_trig_L1_mu_thrName(trig_L1_mu_thrName);
  m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFIsolationContainer(trig_RoI_EF_mu_TrigMuonEFIsolationContainer);
  m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus(trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus);
  m_trigNavVar.set_trig_EF_trigmuonef_EF_mu24i_tight(trig_EF_trigmuonef_EF_mu24i_tight);
  m_trigNavVar.set_trig_EF_trigmuonef_EF_mu36_tight(trig_EF_trigmuonef_EF_mu36_tight);
  
  /*
    if (!m_trigNavVar.isValid()) {
    Error("isTriggerMatched", "Wrong initialization of trigger navigation variables, this will affect trigger matching!!!");
    }
  */
  ///
  
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  // trigger matches
  Bool_t single_el(kFALSE);
  Bool_t single_mu(kFALSE);
  Bool_t di_el(kFALSE);
  Bool_t di_mu(kFALSE);
  Bool_t el_mu(kFALSE);
  
  // trigger chains
  std::vector<TString> chainlist_singlemu = getListOfAlternativeTriggers(getSingleMuonTriggerName());
  std::vector<TString> chainlist_dimu = getListOfAlternativeTriggers(getDiMuonTriggerName());
  std::vector<TString> chainlist_singleel = getListOfAlternativeTriggers(getSingleElectronTriggerName());
  std::vector<TString> chainlist_diel = getListOfAlternativeTriggers(getDiElectronTriggerName());
  std::vector<TString> chainlist_elmu = getListOfAlternativeTriggers(getElectronMuonTriggerName());
  
  // thresholds (manually applied to stay on plateau)
  std::map<TString, Float_t> singlemu_thr;
  std::map<TString, Float_t> singleel_thr;
  std::map<TString, std::pair<Float_t, Float_t> > dimu_thr;
  std::map<TString, std::pair<Float_t, Float_t> > diel_thr;
  std::map<TString, std::pair<Float_t, Float_t> > elmu_thr;
  
  /*
    singlemu_thr["EF_mu18_MG"] = 20000.;
    singlemu_thr["EF_mu18_MG_medium"] = 20000.;
    singlemu_thr["EF_mu24i_tight"] = 25000.;
    singlemu_thr["EF_mu36_tight"] = 37000.;
    
    singleel_thr["EF_e20_medium"] = 21000.;
    singleel_thr["EF_e22_medium"] = 23000.;
    singleel_thr["EF_e22vh_medium1"] = 23000.;
    singleel_thr["EF_e24vhi_medium1"] = 25000.;
    singleel_thr["EF_e60_medium1"] = 61000.;
    
    diel_thr["EF_2e12_medium"] = std::make_pair(13000., 13000.);
    diel_thr["EF_2e12T_medium"] = std::make_pair(13000., 13000.);
    diel_thr["EF_2e12Tvh_medium"] = std::make_pair(13000., 13000.);
    diel_thr["EF_2e12Tvh_loose1"] = std::make_pair(13000., 13000.);
    
    dimu_thr["EF_2mu10_loose"] = std::make_pair(12000., 12000.);
    dimu_thr["EF_2mu13"] = std::make_pair(14000., 14000.);
    dimu_thr["EF_mu18_tight_mu8_EFFS"] = std::make_pair(19000., 9000.);
    
    elmu_thr["EF_e12Tvh_medium1_mu8"] = std::make_pair(13000., 11000.);
    elmu_thr["EF_e24vhi_loose1_mu8"] = std::make_pair(25000., 9000.);
  */
  
  singlemu_thr["EF_mu18_MG"] = 20000.;
  singlemu_thr["EF_mu18_MG_medium"] = 20000.;
  singlemu_thr["EF_mu24i_tight"] = -1.;
  singlemu_thr["EF_mu36_tight"] = -1.;
  
  singleel_thr["EF_e20_medium"] = 21000.;
  singleel_thr["EF_e22_medium"] = 23000.;
  singleel_thr["EF_e22vh_medium1"] = 23000.;
  singleel_thr["EF_e24vhi_medium1"] = -1.;
  singleel_thr["EF_e60_medium1"] = -1.;
  
  diel_thr["EF_2e12_medium"] = std::make_pair(13000., 13000.);
  diel_thr["EF_2e12T_medium"] = std::make_pair(13000., 13000.);
  diel_thr["EF_2e12Tvh_medium"] = std::make_pair(13000., 13000.);
  diel_thr["EF_2e12Tvh_loose1"] = std::make_pair(-1., -1.);
  
  dimu_thr["EF_2mu10_loose"] = std::make_pair(12000., 12000.);
  dimu_thr["EF_2mu13"] = std::make_pair(-1., -1.);
  dimu_thr["EF_mu18_tight_mu8_EFFS"] = std::make_pair(-1., -1.);
  
  elmu_thr["EF_e12Tvh_medium1_mu8"] = std::make_pair(-1., -1.);
  elmu_thr["EF_e24vhi_loose1_mu8"] = std::make_pair(-1., -1.);
  
  
  for (UInt_t i = 0; i < leptons.size(); i++) {
    Int_t my_flavor = leptons[i]->flavor();
    
    if (my_flavor == Analysis::ChargedLepton::MUON && leptons[i]->family() != Muon::CALO) {
      // single muon trigger
      if (analysis_version() == "rel_17") { // bug in p956
	trig_RoI_EF_mu_TrigMuonEFInfoContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfo();
	trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfoStatus();
	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(trig_RoI_EF_mu_TrigMuonEFInfoContainer);
	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus);
      }
      
      for (std::vector<TString>::iterator singlemu = chainlist_singlemu.begin(); singlemu != chainlist_singlemu.end(); ++singlemu) {
	if (m_MuonTriggerMatchTool->match(leptons[i]->Get4Momentum()->Eta(), leptons[i]->Get4Momentum()->Phi(), singlemu->Data())) {
	  if (leptons[i]->Get4Momentum()->Pt() > singlemu_thr[*singlemu])
	    single_mu = kTRUE;
	} // trigger matched
      } // loop over chains in OR
      
      // dimuon trigger
      if (analysis_version() == "rel_17") { // bug in p956
	trig_RoI_EF_mu_TrigMuonEFInfoContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfo();
	trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfoStatus();
	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(trig_RoI_EF_mu_TrigMuonEFInfoContainer);
            m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus);
      }
      for (UInt_t j = i + 1; j < leptons.size(); j++) {
	if (leptons[j]->flavor() == my_flavor && leptons[j]->family() != Muon::CALO) {
	  std::pair<Bool_t, Bool_t> lep1(kFALSE, kFALSE), lep2(kFALSE, kFALSE);
	  
	  for (std::vector<TString>::iterator dimu = chainlist_dimu.begin(); dimu != chainlist_dimu.end(); ++dimu) {
	    m_MuonTriggerMatchTool->matchDimuon(*(leptons[i]->Get4Momentum()), *(leptons[j]->Get4Momentum()), dimu->Data(), lep1, lep2);
	    
	    Bool_t is_actually_matched(kFALSE);
	    Bool_t is_actually_above_threshold(kFALSE);
	    
	    if (analysis_version() == "rel_17") {
	      is_actually_matched = (lep1.first && lep2.first);
	      
	      Float_t max_pt = TMath::Max(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
	      Float_t min_pt = TMath::Min(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
	      is_actually_above_threshold = (max_pt > dimu_thr[*dimu].first && min_pt > dimu_thr[*dimu].second);
	    } else if (analysis_version() == "rel_17_2") {
	      is_actually_matched = ((lep1.first && lep2.second) || (lep1.second && lep2.first));
	      
	      Float_t pt_1 = leptons[i]->Get4Momentum()->Pt();
	      Float_t pt_2 = leptons[j]->Get4Momentum()->Pt();
	      is_actually_above_threshold = (lep1.first && lep2.second && pt_1 > dimu_thr[*dimu].first && pt_2 > dimu_thr[*dimu].second);
	      is_actually_above_threshold = is_actually_above_threshold || (lep2.first && lep1.second && pt_2 > dimu_thr[*dimu].first && pt_1 > dimu_thr[*dimu].second);
	    }
	    
	    if (is_actually_matched && is_actually_above_threshold) {
	      di_mu = kTRUE;
	    } // dimuon trigger matched
	  } // loop over chains in OR
	} // same flavor CB/ST/SA muons
      } // dimuon trigger
    } // CB/ST/SA muons
    else if (my_flavor == Analysis::ChargedLepton::ELECTRON) {
      // single electron trigger
      for (std::vector<TString>::iterator singleel = chainlist_singleel.begin(); singleel != chainlist_singleel.end(); ++singleel) {
	if (m_ElectronTriggerMatchTool->match(leptons[i]->Get4Momentum()->Eta(), leptons[i]->Get4Momentum()->Phi(), singleel->Data())) {
	  if (leptons[i]->Get4Momentum()->Pt() > singleel_thr[*singleel])
	    single_el = kTRUE;
	} // trigger matched
      } // loop over chains in OR
      
      // dielectron trigger
      for (UInt_t j = i + 1; j < leptons.size(); j++) {
	if (leptons[j]->flavor() == my_flavor) {
	  Bool_t lep1(kFALSE), lep2(kFALSE);
	  
	  for (std::vector<TString>::iterator diel = chainlist_diel.begin(); diel != chainlist_diel.end(); ++diel) {
	    m_ElectronTriggerMatchTool->matchDielectron(*(leptons[i]->Get4Momentum()), *(leptons[j]->Get4Momentum()), diel->Data(), lep1, lep2);
	    
	    if (lep1 && lep2) {
	      Float_t max_pt = TMath::Max(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
	      Float_t min_pt = TMath::Min(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
	      
	      if (max_pt > diel_thr[*diel].first && min_pt > diel_thr[*diel].second)
		di_el = kTRUE;
	    } // dielectron trigger matched
	  } // loop over chains in OR
	} // electron pair
      } // dielectron trigger
      
      // electron-muon trigger
      for (UInt_t j = 0; j < leptons.size(); j++) {
	if (leptons[j]->flavor() != my_flavor && leptons[j]->family() != Muon::CALO) {
	  Bool_t lep1(kFALSE), lep2(kFALSE);
	  
	  for (std::vector<TString>::iterator elmu = chainlist_elmu.begin(); elmu != chainlist_elmu.end(); ++elmu) {
	    // TO BE IMPLEMENTED
	    //m_ElectronTriggerMatchTool->matchDielectron(*(leptons[i]->Get4Momentum()), *(leptons[j]->Get4Momentum()), elmu->Data(), lep1, lep2);
	    
	    if (lep1 && lep2) {
	      if (leptons[i]->Get4Momentum()->Pt() > elmu_thr[*elmu].first && leptons[j]->Get4Momentum()->Pt() > elmu_thr[*elmu].second)
		el_mu = kTRUE;
	    } // el-mu trigger matched
	  } // loop over chains in OR
	} // electron-muon pair
      } // electron-muon trigger
    } // electrons
  }
  
  
  //if (passesSingleMuonTrigger() || passesDiMuonTrigger()) {
  //   return (single_mu || di_mu);
  //} else if (passesSingleMuonTrigger() || passesDiMuonTrigger() || passesSingleElectronTrigger() || passesDiElectronTrigger()) {
  return (single_mu || di_mu || single_el || di_el || el_mu); // kostas
  //} else if (passesSingleElectronTrigger() || passesDiElectronTrigger()) {
  //   return (single_el || di_el);
  //}
  
  // should have never been reached
  return kFALSE;
}

Bool_t HiggsllqqAnalysis::passesZ1MassCut(Analysis::Quadrilepton * higgs)
{
  if (50000 < higgs->GetZ1()->Get4Momentum()->M() && 106000 > higgs->GetZ1()->Get4Momentum()->M()) return kTRUE;
  return kFALSE;
}

Bool_t HiggsllqqAnalysis::passesZ2MassCut(Analysis::Quadrilepton * higgs)
{
  const UInt_t nbins(7);
  Double_t mass_llqq[nbins] = {   120,   130, 150, 160, 165, 180, 190 };
  Double_t cut_Z2[nbins]    = {  17.5,  22.5,  30,  30,  35,  40,  50 };
  
  Double_t my_cut(0);
  Double_t my_llqq_mass = higgs->Get4Momentum()->M() / 1000.; // input in GeV
  
  int index(-1);
  for (UInt_t j = 0; j < nbins; j++) {
    if (my_llqq_mass > mass_llqq[j]) index = j;
  }
  if (index == -1)     my_cut = 17.5 * 1000.;
  else if (index == nbins - 1) my_cut = 50.0 * 1000.;
  else my_cut = 1000. * (cut_Z2[index] + (my_llqq_mass - mass_llqq[index]) * (cut_Z2[index + 1] - cut_Z2[index]) / (mass_llqq[index + 1] - mass_llqq[index]));
  
  return (higgs->GetZ2()->Get4Momentum()->M() > my_cut && higgs->GetZ2()->Get4Momentum()->M() < 115000);
}

Bool_t HiggsllqqAnalysis::passesDeltaRCut(Analysis::Quadrilepton * higgs)
{
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  Double_t dr_cut_sameflavor(0.10);
  Double_t dr_cut_differentflavor(0.20);
  //Double_t dr_cut_differentflavor(0.10); // kostas
  
  std::vector<Analysis::ChargedLepton*>::iterator lep_itr_i;
  std::vector<Analysis::ChargedLepton*>::iterator lep_itr_j;
  
  for (lep_itr_i = leptons.begin(); lep_itr_i != leptons.end(); ++lep_itr_i) {
    for (lep_itr_j = lep_itr_i + 1; lep_itr_j != leptons.end(); ++lep_itr_j) {
      Double_t dr_cut = ((*lep_itr_i)->flavor() == (*lep_itr_j)->flavor()) ? dr_cut_sameflavor : dr_cut_differentflavor;
      
      if ((*lep_itr_i)->Get4Momentum()->DeltaR(*((*lep_itr_j)->Get4Momentum())) <= dr_cut) {
	return kFALSE;
      }
    } // loop on lep_j
  } // loop on lep_i
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::passesJpsiVeto(Analysis::Quadrilepton * higgs)
{
  Bool_t result(kTRUE);
  
  if (higgs->channel() == Analysis::Quadrilepton::MU2 || higgs->channel() == Analysis::Quadrilepton::E2) {
    TLorentzVector crossedpair_1 = *(higgs->GetZ1()->GetLepPlus()->Get4Momentum()) + *(higgs->GetZ2()->GetLepMinus()->Get4Momentum());
    TLorentzVector crossedpair_2 = *(higgs->GetZ2()->GetLepPlus()->Get4Momentum()) + *(higgs->GetZ1()->GetLepMinus()->Get4Momentum());
    
    if (crossedpair_1.M() < 5000 || crossedpair_2.M() < 5000) {
      result = kFALSE;
    }
  }
  
  return result;
}

std::vector<Float_t> HiggsllqqAnalysis::getFinalTrackIsoVector(Analysis::Quadrilepton * higgs)
{
  std::vector<Float_t> result;
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  Float_t dr_region(0.20);
  
  std::vector<Analysis::ChargedLepton*>::iterator lep_itr_i;
  std::vector<Analysis::ChargedLepton*>::iterator lep_itr_j;
  
  for (lep_itr_i = leptons.begin(); lep_itr_i != leptons.end(); ++lep_itr_i) {
    Double_t recalc_iso = (*lep_itr_i)->ptcone20();
    
    // subtract ID track pt from leptons within the cone
    for (lep_itr_j = leptons.begin(); lep_itr_j != leptons.end(); ++lep_itr_j) {
      if ((*lep_itr_i)->Get4Momentum()->DeltaR(*((*lep_itr_j)->Get4Momentum())) <= dr_region && lep_itr_i != lep_itr_j) {
	recalc_iso = recalc_iso - (*lep_itr_j)->Get4Momentum_ID()->Pt();
      }
    }
    
    result.push_back(recalc_iso);
  }
  
  return result;
}

Bool_t HiggsllqqAnalysis::passesTrackIsolationCut(Analysis::Quadrilepton * higgs)
{
  std::vector<Float_t> recalc_iso = getFinalTrackIsoVector(higgs);
  
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  Float_t isolation_cut(0.15);
  
  for (UInt_t i = 0; i < leptons.size(); i++) {
    if (recalc_iso[i] / leptons[i]->Get4Momentum()->Pt() > isolation_cut) return kFALSE;
  }
  
  return kTRUE;
}

Float_t HiggsllqqAnalysis::getCorrectedCaloIso(Analysis::ChargedLepton * lep)
{
  Int_t nPV(0);
  for (Int_t i = 0; i < ntuple->vxp.n(); i++) {
    if (ntuple->vxp[i].trk_n() >= 2)
      nPV++;
  }
  
  Float_t result(-1);
  
  if (lep->flavor() == Analysis::ChargedLepton::MUON) {
    result = lep->GetMuon()->etcone20();
  } else if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    if (analysis_version() == "rel_17") { // rel. 17
      result = CaloIsoCorrection::GetNPVCorrectedIsolation(nPV,
							   lep->GetElectron()->etas2(),
							   0.20,
							   isMC(),
							   lep->GetElectron()->Etcone20(),
							   CaloIsoCorrection::ELECTRON);
    } // rel. 17
    else if (analysis_version() == "rel_17_2") { // rel. 17.2
      if (useTopoIso()) {
	result = CaloIsoCorrection::GetPtEDCorrectedTopoIsolation(lep->GetElectron()->ED_median(),
								  lep->Get4Momentum_SA()->E(),
								  lep->GetElectron()->etas2(),
								  lep->GetElectron()->etap(),
								  lep->GetElectron()->cl_eta(),
								  0.20,
								  isMC(),
								  lep->GetElectron()->topoEtcone20(),
								  false,
								  CaloIsoCorrection::ELECTRON,
								  CaloIsoCorrection::REL17);
      } else {
	result = CaloIsoCorrection::GetPtNPVCorrectedIsolation(nPV,
							       lep->Get4Momentum_SA()->E(),
							       lep->GetElectron()->etas2(),
							       lep->GetElectron()->etap(),
							       lep->GetElectron()->cl_eta(),
							       0.20,
							       isMC(),
							       lep->GetElectron()->Etcone20(),
							       false,
							       CaloIsoCorrection::ELECTRON,
							       CaloIsoCorrection::REL17);
      }
    } // rel. 17.2
  }
  
  return result;
}

std::vector<Float_t> HiggsllqqAnalysis::getFinalCaloIsoVector(Analysis::Quadrilepton * higgs)
{
  std::vector<Float_t> result;
  result.clear();
  
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  Float_t dr_region(0.18);
  
  std::vector<Analysis::ChargedLepton*>::iterator lep_itr;
  
  // store pile-up corrected isolation
  Double_t corrected_calo_iso[4] = {0., 0., 0., 0.};
  Double_t calo_deposit[4] = {0., 0., 0., 0.};
  Int_t n(0);
  
  for (lep_itr = leptons.begin(); lep_itr != leptons.end(); ++lep_itr) {
    corrected_calo_iso[n] = getCorrectedCaloIso((*lep_itr));
    
    if ((*lep_itr)->flavor() == Analysis::ChargedLepton::MUON) {
      calo_deposit[n] = 0.;
    } else if ((*lep_itr)->flavor() == Analysis::ChargedLepton::ELECTRON) {
      calo_deposit[n] = (*lep_itr)->Get4Momentum()->Et();
    } else {
      Error("getFinalCaloIsoVector", "Unexpected lepton flavor %d", (*lep_itr)->flavor());
    }
    
    n++;
  }
  
  for (UInt_t i = 0; i < leptons.size(); i++) {
    Double_t recalc_iso = corrected_calo_iso[i];
    
    // subtract energy deposit from leptons within the cone
    for (UInt_t j = 0; j < leptons.size(); j++) {
      if (leptons[i]->Get4Momentum()->DeltaR(*(leptons[j]->Get4Momentum())) <= dr_region && i != j) {
	recalc_iso = recalc_iso - calo_deposit[j];
      }
    }
    
    result.push_back(recalc_iso);
  }
  
  return result;
}

Bool_t HiggsllqqAnalysis::passesCaloIsolationCut(Analysis::Quadrilepton * higgs)
{
  std::vector<Float_t> recalc_iso = getFinalCaloIsoVector(higgs);
  
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  Float_t isolation_cut_muon(0.30);
  Float_t isolation_cut_electron(0.30);
  if (analysis_version() == "rel_17_2") isolation_cut_electron = 0.20;
  
  Float_t isolation_cut_sa(0.15); // standalone muons
  
  for (UInt_t i = 0; i < leptons.size(); i++) {
    Float_t isolation_cut = (leptons[i]->flavor() == Analysis::ChargedLepton::MUON) ? isolation_cut_muon : isolation_cut_electron;
    
    if (leptons[i]->flavor() == Analysis::ChargedLepton::MUON) {
      if (leptons[i]->GetMuon()->isStandAloneMuon() == 1) {
	if (recalc_iso[i] / leptons[i]->Get4Momentum()->Pt() > isolation_cut_sa) return kFALSE;
      }
    }
    if (recalc_iso[i] / leptons[i]->Get4Momentum()->Pt() > isolation_cut) return kFALSE;
  }
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::hasGoodImpactParameterSignificance(Analysis::Quadrilepton * higgs)
{
  if (higgs->Get4Momentum()->M() > 190000 && analysis_version() == "rel_17")
    return kTRUE;
  
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  std::vector<Analysis::ChargedLepton*>::iterator lep_itr;
  
  for (lep_itr = leptons.begin(); lep_itr != leptons.end(); ++lep_itr) {
    Double_t d0 = (*lep_itr)->d0();
    
    if (isMC() && analysis_version() == "rel_17_2") {
      Bool_t apply_correction(kTRUE);
      
      if ((*lep_itr)->flavor() == Analysis::ChargedLepton::MUON) {
	if ((*lep_itr)->GetMuon()->isStandAloneMuon() == 1) {
	  apply_correction = kFALSE; // SA muons don't have nBL
	}
      }
      
      if (apply_correction) {
	Int_t lep_index = ((*lep_itr)->flavor() == Analysis::ChargedLepton::MUON) ? (*lep_itr)->GetMuon()->GetIndex() : (*lep_itr)->GetElectron()->GetIndex();
	Int_t nBL = ((*lep_itr)->flavor() == Analysis::ChargedLepton::MUON) ? (*lep_itr)->GetMuon()->nBLHits() : (*lep_itr)->GetElectron()->nBLHits();
	
	d0 += -0.002;
	
	d0 += getD0SmearSigma(lep_index, nBL, (*lep_itr)->Get4Momentum()->Pt(), (*lep_itr)->Get4Momentum()->Eta());; // https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsRecommendationsForSummer2012
      }
    }
    
    if ((*lep_itr)->flavor() == Analysis::ChargedLepton::MUON) {
      if (TMath::Abs(d0 / (*lep_itr)->d0_sig()) > 3.5) return kFALSE;
    } else if ((*lep_itr)->flavor() == Analysis::ChargedLepton::ELECTRON) {
      if (TMath::Abs(d0 / (*lep_itr)->d0_sig()) > 6.5) return kFALSE;
    }
  }
  
  return kTRUE;
}

Int_t HiggsllqqAnalysis::getBestCandidateCut()
{
  // return the ID of the last cut passed by the best candidate in this event (return -1 if not found)
  // this means that Quadrilepton::lastcut() is corrected for the output of isClosestToZZ()
  // in such a way that the ID returned can apply to the full event itself
  
  Int_t last(-9999);
  
  std::vector<Analysis::Quadrilepton *>::iterator H_itr;
  
  for (H_itr = m_Quadrileptons.begin(); H_itr != m_Quadrileptons.end(); ++H_itr) {
    // get the last cut passed by the quadrilepton
    Int_t its_cut = (*H_itr)->lastcut();
    
    // correct the last cut passed, taking into account the fact that the request that the quadrilepton is the best
    // in the event (closest to ZZ) must be done immediately after the DeltaR cut
    if (its_cut >= HllqqQuadrileptonQuality::best) {
      if (!isClosestToZZ((*H_itr))) {
	its_cut = HllqqQuadrileptonQuality::best - 1;
      }
    }
    
    // update, if it is the case, the result
    if (its_cut > last) {
      last = its_cut;
    }
  } // loop over saved quadrileptons
  
  return last;
}

Bool_t HiggsllqqAnalysis::execute_analysis()
{
  // internal utility counters/flags
  m_processedEntries++;
  m_called_getGoodLeptons = kFALSE;
  m_called_getGoodObjects = kFALSE;
  
  // bugfix for mc12a samples produced using FRONTIER
  if (analysis_version() == "rel_17_2" && isMC()) {
    if (getTriggerInfo("L1_RD0_FILLED") != 1) return kTRUE;
  }
  
  
  // update generated entries' histo
  m_generatedEntriesHisto->Fill("raw", 1);
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix", (!hasPowHegZZBug()) ? 1 : 0);
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix_EvWeight", (!hasPowHegZZBug() && isMC()) ? ntuple->eventinfo.mc_event_weight() : 0);
  m_generatedEntriesHisto->Fill("with_ggF_as_well", (!hasPowHegZZBug() && isMC()) ? ntuple->eventinfo.mc_event_weight() * getggFWeight() : 0);
  m_generatedEntriesHisto->Fill("with_pu_and_vertex", (!hasPowHegZZBug() && isMC()) ? getEventWeight() * getggFWeight() * getVertexZWeight() : 0);
  
  
  m_cutflowStruct.run = (isMC()) ? ntuple->eventinfo.mc_channel_number() : ntuple->eventinfo.RunNumber();
  m_cutflowStruct.event = ntuple->eventinfo.EventNumber();
  
  for (UInt_t i = 0; i < m_Channels.size(); i++) {
    UInt_t chan = m_Channels.at(i);
    setChannel(chan);
    
    Int_t last_event = getLastCutPassed(); //GetLastCutPassed();
    Int_t last_cand  = getBestCandidateCut();
    
    
    // update cutflows
    m_EventCutflow[chan].addCutCounter(last_event, 1);
    m_EventCutflow_rw[chan].addCutCounter(last_event, getEventWeight());//*GetSFWeight());
    if (last_event >= HllqqCutFlow::SelectedLeptons && last_cand != -9999) m_CandidateCutflow[chan].addCutCounter(last_cand, 1);
    
    std::vector<Analysis::ChargedLepton*>::iterator lep;
    for (lep = m_Muons.begin(); lep != m_Muons.end(); ++lep) {
      if ((*lep)->lastcut() != -1)
	m_MuonCutflow[chan].addCutCounter((*lep)->lastcut(), 1);
    }
    for (lep = m_Electrons.begin(); lep != m_Electrons.end(); ++lep) {
      if ((*lep)->lastcut() != -1)
	m_ElectronCutflow[chan].addCutCounter((*lep)->lastcut(), 1);
    }
    
    std::vector<Analysis::Jet*>::iterator jet;
    for (jet = m_Jets.begin(); jet != m_Jets.end(); ++jet) {
      if ((*jet)->lastcut() != -1)
	m_JetCutflow[chan].addCutCounter((*jet)->lastcut(), 1);
    }
    
    // fill the ntuple with the event cutflow
    m_cutflowStruct.analysis = chan;
    m_cutflowStruct.last = last_event;
    m_TreeCutflow->Fill();
    m_selectionEfficiencyVsNvx[chan]->Fill((last_event >= HllqqCutFlow::SelectedLeptons && m_GoodQuadrileptons.size() > 0), getNumberOfGoodVertices()); 
  }
  
  
  // fill the ntuple with the Higgs candidates in m_Quadrileptons
  for (std::vector<Analysis::Quadrilepton *>::iterator H_itr = m_Quadrileptons.begin(); H_itr != m_Quadrileptons.end(); ++H_itr) {
    fillCandidateStruct(&m_candStruct, (*H_itr));
    m_TreeCandidates->Fill();
  }
  
  // clear objects
  for (UInt_t k = 0; k < m_Quadrileptons.size(); k++)
    delete m_Quadrileptons.at(k);
  for (UInt_t k = 0; k < m_Dileptons.size(); k++)
    delete m_Dileptons.at(k);
  for (UInt_t k = 0; k < m_Muons.size(); k++)
    delete m_Muons.at(k);
  for (UInt_t k = 0; k < m_Electrons.size(); k++)
    delete m_Electrons.at(k);
  for (UInt_t k = 0; k < m_Jets.size(); k++)
    delete m_Jets.at(k);
  
  m_Muons.clear();
  m_Electrons.clear();
  m_Jets.clear();
  m_GoodMuons.clear();
  m_GoodElectrons.clear();
  m_GoodJets.clear();
  m_Dileptons.clear();
  m_Quadrileptons.clear();
  m_GoodQuadrileptons.clear();
  
  return kTRUE;
}

Bool_t HiggsllqqAnalysis::finalize_analysis()
{
  // print the number of processed entries
  Info("finalize_analysis", "Analysis terminated, processed %lld entries out of %lld available entries", m_processedEntries, m_entriesInChain);
  
  // print the cutflows
  for (UInt_t i = 0; i < m_Channels.size(); i++) {
    m_EventCutflow[i].print();
    m_EventCutflow_rw[i].print();
    m_CandidateCutflow[i].print();
    m_ElectronCutflow[i].print();
    m_MuonCutflow[i].print();
    m_JetCutflow[i].print();
  }
  
  // clear the maps
  m_CrossSection.clear();
  m_SignalSampleMass.clear();
  
  // clear the cutflows
  m_Channels.clear();
  m_EventCutflow.clear();
  m_EventCutflow_rw.clear();
  m_CandidateCutflow.clear();
  m_ElectronCutflow.clear();
  m_MuonCutflow.clear();
  m_JetCutflow.clear();
  
  // write the file
  m_outputFile->Write();
  m_outputFile->Close();
  
  return kTRUE;
}

Int_t HiggsllqqAnalysis::getTriggerInfo(TString chain)
{
  Int_t result(-1); // to distinguish between YES, NO and UNDEFINED
  Bool_t result_bool(kFALSE);
  
  std::vector<TString> all_triggers = getListOfAlternativeTriggers(chain);
  
  for (std::vector<TString>::iterator myChain = all_triggers.begin(); myChain != all_triggers.end(); ++myChain) {
    if (m_TrigDecisionToolD3PD->GetConfigSvc().IsConfigured(myChain->Data())) {
      if (m_TrigDecisionToolD3PD->IsPassed(myChain->Data())) result_bool = kTRUE;
      else result = 0; // needed to track chain existence
    } // defined chain
  } // loop over chains
  
  if (result_bool) result = 1; // otherwise keep result as it is (typically 0, might be -1 if no chain was found)
  
  return result;
}

Int_t HiggsllqqAnalysis::getPeriod()
{
  // MC in this function is mc11c-specific
  
  if (!isMC()) {
    // DATA
    if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_A").second) {
      // A
      return DataPeriod::y2011_A;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_B").second) {
      // B
      return DataPeriod::y2011_B;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_C").second) {
      // C
      return DataPeriod::y2011_C;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_D").second) {
      // D
      return DataPeriod::y2011_D;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_E").second) {
      // E
      return DataPeriod::y2011_E;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_F").second) {
      // F
      return DataPeriod::y2011_F;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_G").second) {
      // G
      return DataPeriod::y2011_G;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_H").second) {
      // H
      return DataPeriod::y2011_H;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_I").second) {
      // I
      return DataPeriod::y2011_I;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_J").second) {
      // J
      return DataPeriod::y2011_J;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_K").second) {
      // K
      return DataPeriod::y2011_K;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_L").second) {
      // L
         return DataPeriod::y2011_L;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2011_M").second) {
      // M
      return DataPeriod::y2011_M;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2012_A").second) {
      // 2012: A
      return DataPeriod::y2012_A;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2012_B").second) {
      // B
      return DataPeriod::y2012_B;
    } else {
      // SAFETY : consider last runs as belonging to the latest period
      return DataPeriod::y2012_B;
    }
  } else {
    // MC
    if (analysis_version() == "rel_17_2") { // mc12a
      if (ntuple->eventinfo.RunNumber() == 195847) {
	return DataPeriod::y2012_AllYear;
      }
    } else if (analysis_version() == "rel_17") { // mc11c
      if (ntuple->eventinfo.RunNumber() == 180164) {
	// B-D
	return DataPeriod::y2011_BD;
      } else if (ntuple->eventinfo.RunNumber() == 183003) {
	// E-H
	return DataPeriod::y2011_EH;
	//} else if (ntuple->eventinfo.RunNumber() == 185649) { // mc11a
      } else if (ntuple->eventinfo.RunNumber() == 186169 || ntuple->eventinfo.RunNumber() ==  185649) {
	// I-J-K
	TRandom3 random3;
	random3.SetSeed(ntuple->eventinfo.mc_channel_number() * ntuple->eventinfo.EventNumber());
	Double_t rd = random3.Uniform();
	Float_t periodI = m_PileupReweighter->GetIntegratedLumi(185353, 186493); //337.54; //399.206; //Changes Stelio's 29th May 2012
	Float_t periodJ = m_PileupReweighter->GetIntegratedLumi(186516, 186755); //226.39; //232.931;
	Float_t periodK = m_PileupReweighter->GetIntegratedLumi(186873, 187815); //590.36; //660.211;
	Double_t fracI = (periodI) / (periodI + periodJ + periodK);
	Double_t fracJ = (periodJ) / (periodI + periodJ + periodK);
	
	if (rd < fracI) {
	  // I
	  return DataPeriod::y2011_I;
	} else if (rd < fracJ + fracI) {
	  // J
	  return DataPeriod::y2011_J;
	} else {
	  // K
	  return DataPeriod::y2011_K;
	}
      } else if (ntuple->eventinfo.RunNumber() == 186275 || ntuple->eventinfo.RunNumber() == 185761 || ntuple->eventinfo.RunNumber() == 189751) {   //185761,186275,189751
	// L-M
	return DataPeriod::y2011_LM;
      }
    }
  }
  
  
  Warning("GetPeriod", "Unable to retrieve data period number for RunNumber %ud", ntuple->eventinfo.RunNumber());
  return -1;
}

Float_t HiggsllqqAnalysis::getEventWeight()
{
  Float_t result = getPileupWeight();
  
  if (isMC()) {
    result = result * ntuple->eventinfo.mc_event_weight();
  }
  
  return result;
}

Float_t HiggsllqqAnalysis::getPileupWeight()
{
  Float_t result = 1.;
  
  if (isMC()) {
    result = m_PileupReweighter->GetCombinedWeight(ntuple->eventinfo.RunNumber(), ntuple->eventinfo.mc_channel_number(), ntuple->eventinfo.averageIntPerXing());
  }
  
  if (result < 0) {
    Warning("getPileupWeight", "Negative weight from pile-up tool!");
  }
  
  return result;
}

Float_t HiggsllqqAnalysis::getVertexZWeight()
{
  Float_t result = 1.;
  
  if (isMC() && analysis_version() == "rel_17_2") {
    Float_t mc_vertex_z(0.);
    
    Int_t i(0);
    
    for (i = 0; i < ntuple->mc.n(); i++) {
      if (ntuple->mc[i].vx_z() != 0.) {
	mc_vertex_z = ntuple->mc[i].vx_z();
	break;
      }
    }
    
    if (mc_vertex_z == 0. || i > 20) {
      Warning("getVertexZWeight", "We expect a non-zero mc_vx_z within the first 10 entries, but it is not the case!");
    }
    
    result = m_VertexPositionReweighter->GetWeight(mc_vertex_z);
  }
  
  return result;
}

Float_t HiggsllqqAnalysis::getLeptonWeight(Analysis::ChargedLepton * lep)
{
  Float_t result(1);
  
  if (isMC()) {
    if (lep->flavor() == Analysis::ChargedLepton::MUON) {
      Float_t overall(1);
      
      // charge dependence is available for CB+ST (2011, 2012), SA (2012), but not for CALO
      
      if (lep->family() != Muon::CALO) {
	if (lep->GetMuon()->isStandAloneMuon() != 1)
	  overall = m_MuonEffSF->scaleFactor((Int_t)lep->charge(), *(lep->Get4Momentum()));
	else {
	  if (analysis_version() == "rel_17")
	    overall = m_MuonEffSFSA->scaleFactor(*(lep->Get4Momentum()));
	  else if (analysis_version() == "rel_17_2") 
	    overall = m_MuonEffSFSA->scaleFactor((Int_t)lep->charge(), *(lep->Get4Momentum()));
	}
      } else {
	overall = m_MuonEffSFCalo->scaleFactor(*(lep->Get4Momentum()));
      }
      
      result = overall;
    } else if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
      Float_t reco(1.), id(1.);
      
      Int_t representative_run_number = m_PileupReweighter->GetRandomRunNumber(ntuple->eventinfo.RunNumber());
      
      if (analysis_version() == "rel_17") {
	reco = m_ElectronEffSF->scaleFactor(lep->GetElectron()->cl_eta(), lep->Get4Momentum()->Et(), 4, 0, 6, true, representative_run_number).first;
	id   = m_ElectronEffSF->scaleFactor(lep->GetElectron()->cl_eta(), lep->Get4Momentum()->Et(), 5, 0, 6, true, representative_run_number).first;
      } else if (analysis_version() == "rel_17_2") {
	reco = m_ElectronEffSF->scaleFactor(lep->GetElectron()->cl_eta(), lep->Get4Momentum()->Et(), 4, 0, 8, true, representative_run_number).first;
	id   = m_ElectronEffSF->scaleFactor(lep->GetElectron()->cl_eta(), lep->Get4Momentum()->Et(), 30, 0, 8, true, representative_run_number).first;
      }
      
      result = reco * id;
    }
  }
  
  return result;
}

void HiggsllqqAnalysis::initCrossSections()
{
  m_CrossSection.clear();
  m_SignalSampleMass.clear();
  
  // values are in fb
  m_CrossSection[107650] = 827375;
  m_CrossSection[107651] = 166625;
  m_CrossSection[107652] = 50375;
  m_CrossSection[107653] = 14000;
  m_CrossSection[107654] = 3375;
  m_CrossSection[107655] = 1000;
  m_CrossSection[107660] = 822125;
  m_CrossSection[107661] = 166000;
  m_CrossSection[107662] = 49500;
  m_CrossSection[107663] = 13875;
  m_CrossSection[107664] = 3500;
  m_CrossSection[107665] = 1000;
  m_CrossSection[107670] = 828125;
  m_CrossSection[107671] = 167375;
  m_CrossSection[107672] = 50375;
  m_CrossSection[107673] = 13750;
  m_CrossSection[107674] = 3500;
  m_CrossSection[107675] = 1000;
  m_CrossSection[109291] = 146.8;
  m_CrossSection[109292] = 91.54;
  m_CrossSection[116960] = 20.701 * 1.4;
  m_CrossSection[116961] = 18.8029 * 1.4;
  m_CrossSection[116962] = 10.505 * 1.4;
  m_CrossSection[116963] = 7.30463 * 1.4;
  m_CrossSection[116965] = 21.516 * 1.4;
  m_CrossSection[116966] = 19.6674 * 1.4;
  m_CrossSection[116967] = 10.516 * 1.4;
  m_CrossSection[116968] = 7.93834 * 1.4;
  m_CrossSection[116950] = 1.4 * 756.32;
  m_CrossSection[116951] = 1.4 * 432.25;
  m_CrossSection[116952] = 1.4 * 176;
  m_CrossSection[116953] = 1.4 * 96.75;
  m_CrossSection[116955] = 1.4 * 730.24;
  m_CrossSection[116956] = 1.4 * 432.25;
  m_CrossSection[116957] = 1.4 * 179.3;
  m_CrossSection[116958] = 1.4 * 92.3962;
  m_CrossSection[105200] = 91550.6;
  m_CrossSection[109345] = 12707.2;
  m_CrossSection[109346] = 515.2;
  m_CrossSection[116611] = 5.89;
  m_CrossSection[116612] = 13.6;
  m_CrossSection[116761] = 0.944;
  m_CrossSection[116762] = 1.69;
  m_CrossSection[116763] = 2.81;
  m_CrossSection[116764] = 4.26;
  m_CrossSection[116765] = 7.46;
  m_CrossSection[116766] = 8.66;
  m_CrossSection[116767] = 9.24;
  m_CrossSection[116768] = 8.95;
  m_CrossSection[116769] = 7.41;
  m_CrossSection[116770] = 3.86;
  m_CrossSection[116771] = 1.92;
  m_CrossSection[116772] = 1.87;
  m_CrossSection[116773] = 2.4;
  m_CrossSection[116774] = 4.12;
  m_CrossSection[116775] = 9.51;
  m_CrossSection[116776] = 12.5;
  m_CrossSection[116777] = 13.4;
  m_CrossSection[116779] = 13.1;
  m_CrossSection[116780] = 12.3;
  m_CrossSection[116781] = 10.6;
  m_CrossSection[116782] = 9.31;
  m_CrossSection[116783] = 8.28;
  m_CrossSection[116784] = 7.53;
  m_CrossSection[116785] = 7.03;
  m_CrossSection[116786] = 6.92;
  m_CrossSection[116787] = 7.08;
  m_CrossSection[116788] = 6.41;
  m_CrossSection[116789] = 5.53;
  m_CrossSection[116790] = 4.69;
  m_CrossSection[116791] = 3.91;
  m_CrossSection[116792] = 3.25;
  m_CrossSection[116793] = 2.70;
  m_CrossSection[116794] = 2.24;
  m_CrossSection[116795] = 1.86;
  m_CrossSection[116796] = 1.56;
  m_CrossSection[116797] = 1.29;
  m_CrossSection[116798] = 1.08;
  m_CrossSection[116799] = 0.902;
  m_CrossSection[116621] = 0.481;
  m_CrossSection[116622] = 1.65;
  m_CrossSection[125063] = 0.0665;
  m_CrossSection[125064] = 0.124;
  m_CrossSection[125065] = 0.214;
  m_CrossSection[125066] = 0.337;
  m_CrossSection[125067] = 0.627;
  m_CrossSection[125068] = 0.751;
  m_CrossSection[125069] = 0.823;
  m_CrossSection[125070] = 0.819;
  m_CrossSection[125071] = 0.694;
  m_CrossSection[125072] = 0.373;
  m_CrossSection[125073] = 0.196;
  m_CrossSection[125074] = 0.198;
  m_CrossSection[125075] = 0.26;
  m_CrossSection[125076] = 0.458;
  m_CrossSection[125077] = 1.09;
  m_CrossSection[125078] = 1.47;
  m_CrossSection[125079] = 1.6;
  m_CrossSection[125081] = 1.63;
  m_CrossSection[125082] = 1.56;
  m_CrossSection[125083] = 1.38;
  m_CrossSection[125084] = 1.21;
  m_CrossSection[125085] = 1.06;
  m_CrossSection[125086] = 0.936;
  m_CrossSection[125087] = 0.822;
  m_CrossSection[125088] = 0.72;
  m_CrossSection[125089] = 0.605;
  m_CrossSection[125090] = 0.512;
  m_CrossSection[125091] = 0.441;
  m_CrossSection[125092] = 0.387;
  m_CrossSection[125093] = 0.344;
  m_CrossSection[125094] = 0.308;
  m_CrossSection[125095] = 0.278;
  m_CrossSection[125096] = 0.251;
  m_CrossSection[125097] = 0.228;
  m_CrossSection[125098] = 0.208;
  m_CrossSection[125099] = 0.19;
  m_CrossSection[125100] = 0.174;
  m_CrossSection[125101] = 0.159;
  
  // sample masses (values are in GeV)
  m_SignalSampleMass[116611] = 130;
  m_SignalSampleMass[116612] = 200;
  m_SignalSampleMass[116761] = 110;
  m_SignalSampleMass[116762] = 115;
  m_SignalSampleMass[116763] = 120;
  m_SignalSampleMass[116764] = 125;
  m_SignalSampleMass[116765] = 135;
  m_SignalSampleMass[116766] = 140;
  m_SignalSampleMass[116767] = 145;
  m_SignalSampleMass[116768] = 150;
  m_SignalSampleMass[116769] = 155;
  m_SignalSampleMass[116770] = 160;
  m_SignalSampleMass[116771] = 165;
  m_SignalSampleMass[116772] = 170;
  m_SignalSampleMass[116773] = 175;
  m_SignalSampleMass[116774] = 180;
  m_SignalSampleMass[116775] = 185;
  m_SignalSampleMass[116776] = 190;
  m_SignalSampleMass[116777] = 195;
  m_SignalSampleMass[116779] = 210;
  m_SignalSampleMass[116780] = 220;
  m_SignalSampleMass[116781] = 240;
  m_SignalSampleMass[116782] = 260;
  m_SignalSampleMass[116783] = 280;
  m_SignalSampleMass[116784] = 300;
  m_SignalSampleMass[116785] = 320;
  m_SignalSampleMass[116786] = 340;
  m_SignalSampleMass[116787] = 360;
  m_SignalSampleMass[116788] = 380;
  m_SignalSampleMass[116789] = 400;
  m_SignalSampleMass[116790] = 420;
  m_SignalSampleMass[116791] = 440;
  m_SignalSampleMass[116792] = 460;
  m_SignalSampleMass[116793] = 480;
  m_SignalSampleMass[116794] = 500;
  m_SignalSampleMass[116795] = 520;
  m_SignalSampleMass[116796] = 540;
  m_SignalSampleMass[116797] = 560;
  m_SignalSampleMass[116798] = 580;
  m_SignalSampleMass[116799] = 600;
}

Float_t HiggsllqqAnalysis::getCrossSectionWeight()
{
  Float_t result(-9999.9);
  
  if (isMC()) {
    UInt_t chan = ntuple->eventinfo.mc_channel_number();
    // ZZ cross section with M_ZZ dependent K factor
    if (chan == 109292 || chan == 109291) {
      Float_t xsec = CrossSections::GetBkgCrossSection(chan, kFALSE);
      Float_t mcfm = GetMzzWeightFromMCFM(getTruthZZMass() / 1000.);
      result = xsec * mcfm;
    }
    // Other cross sections
    else {
      Float_t xsec = -1.;
      std::map<UInt_t, Float_t>::iterator it;
      it = m_CrossSection.find(chan);
      if (it != m_CrossSection.end()) {
	xsec = it->second;
      }
      result = xsec;
    }
  }
  
  return result;
}

Float_t HiggsllqqAnalysis::getTruthZZMass()
{
  // return the mass of the ZZ system, using the first 2 Z bosons, if any
  TLorentzVector firstZ, secondZ;
  Int_t howMany(0); // how many Zs have been retrieved?
  
  if (isMC()) {
    for (Int_t i = 0; i < ntuple->mc.n() && howMany < 2; i++) {
      D3PDReader::TruthParticleD3PDObjectElement *particle = &(ntuple->mc[i]);
      
      if (TMath::Abs(ntuple->mc[i].pdgId()) == 23) {
	if (howMany == 0)
	  firstZ = CommonTools::getVector(particle);
	else if (howMany == 1)
	  secondZ = CommonTools::getVector(particle);
	howMany++;
      }
    }
  }
  
  return (firstZ + secondZ).M();
}

void HiggsllqqAnalysis::setTreeCandidatesBranches(TTree * tree, Analysis::CandidateStruct * cand)
{
  tree->Branch("last", &(cand->last), "last/I");
  tree->Branch("closesttozz", &(cand->closesttozz), "closesttozz/I");
  tree->Branch("selected", &(cand->selected), "selected/I");
  tree->Branch("type", &(cand->type), "type/I");
  tree->Branch("run", &(cand->run), "run/i");
  tree->Branch("event", &(cand->event), "event/i");
  tree->Branch("lbn", &(cand->lbn), "lbn/i");
  tree->Branch("n_vx", &(cand->n_vx), "n_vx/I");
  tree->Branch("hfor", &(cand->hfor), "hfor/I");
  tree->Branch("top_weight", &(cand->top_weight), "top_weight/I");
  tree->Branch("powhegbug_weight", &(cand->powhegbug_weight), "powhegbug_weight/I");
  tree->Branch("actualIntPerXing", &(cand->actualIntPerXing), "actualIntPerXing/I");
  tree->Branch("averageIntPerXing", &(cand->averageIntPerXing), "averageIntPerXing/I");
  tree->Branch("xsec_weight", &(cand->xsec_weight), "xsec_weight/F");
  tree->Branch("processed_entries", &(cand->processed_entries), "processed_entries/L");
  tree->Branch("generated_entries", &(cand->generated_entries), "generated_entries/L");
  tree->Branch("pu_weight", &(cand->pu_weight), "pu_weight/F");
  tree->Branch("vxz_weight", &(cand->vxz_weight), "vxz_weight/F");
  tree->Branch("OQ_weight", &(cand->OQ_weight), "OQ_weight/F");
  tree->Branch("ggF_weight", &(cand->ggF_weight), "ggF_weight/F");
  tree->Branch("trigSF_weight", &(cand->trigSF_weight), "trigSF_weight/F");
  tree->Branch("Z1_lepplus_weight", &(cand->Z1_lepplus_weight), "Z1_lepplus_weight/F");
  tree->Branch("Z1_lepminus_weight", &(cand->Z1_lepminus_weight), "Z1_lepminus_weight/F");
  tree->Branch("Z2_lepplus_weight", &(cand->Z2_lepplus_weight), "Z2_lepplus_weight/F");
  tree->Branch("Z2_lepminus_weight", &(cand->Z2_lepminus_weight), "Z2_lepminus_weight/F");
  tree->Branch("H_m_truth", &(cand->H_m_truth), "H_m_truth/F");
  tree->Branch("ZZ_m_truth", &(cand->ZZ_m_truth), "ZZ_m_truth/F");
  tree->Branch("Z1_m_truth", &(cand->Z1_m_truth), "Z1_m_truth/F");
  tree->Branch("Z2_m_truth", &(cand->Z2_m_truth), "Z2_m_truth/F");
  tree->Branch("H_m_constrained", &(cand->H_m_constrained), "H_m_constrained/F");
  tree->Branch("Z1_m_constrained", &(cand->Z1_m_constrained), "Z1_m_constrained/F");
  tree->Branch("Z2_m_constrained", &(cand->Z2_m_constrained), "Z2_m_constrained/F");
  tree->Branch("Z1_chiSq", &(cand->Z1_chiSq), "Z1_chiSq/F");
  tree->Branch("Z2_chiSq", &(cand->Z2_chiSq), "Z2_chiSq/F");
  tree->Branch("H_m", &(cand->H_m), "H_m/F");
  tree->Branch("H_pt", &(cand->H_pt), "H_pt/F");
  tree->Branch("H_eta", &(cand->H_eta), "H_eta/F");
  tree->Branch("H_phi", &(cand->H_phi), "H_phi/F");
  tree->Branch("Z1_m", &(cand->Z1_m), "Z1_m/F");
  tree->Branch("Z1_pt", &(cand->Z1_pt), "Z1_pt/F");
  tree->Branch("Z1_eta", &(cand->Z1_eta), "Z1_eta/F");
  tree->Branch("Z1_phi", &(cand->Z1_phi), "Z1_phi/F");
  tree->Branch("Z2_m", &(cand->Z2_m), "Z2_m/F");
  tree->Branch("Z2_pt", &(cand->Z2_pt), "Z2_pt/F");
  tree->Branch("Z2_eta", &(cand->Z2_eta), "Z2_eta/F");
  tree->Branch("Z2_phi", &(cand->Z2_phi), "Z2_phi/F");
  tree->Branch("Z1_lepplus_m", &(cand->Z1_lepplus_m), "Z1_lepplus_m/F");
  tree->Branch("Z1_lepplus_pt", &(cand->Z1_lepplus_pt), "Z1_lepplus_pt/F");
  tree->Branch("Z1_lepplus_pt_truth", &(cand->Z1_lepplus_pt_truth), "Z1_lepplus_pt_truth/F");
  tree->Branch("Z1_lepplus_pt_constrained", &(cand->Z1_lepplus_pt_constrained), "Z1_lepplus_pt_constrained/F");
  tree->Branch("Z1_lepplus_eta", &(cand->Z1_lepplus_eta), "Z1_lepplus_eta/F");
  tree->Branch("Z1_lepplus_phi", &(cand->Z1_lepplus_phi), "Z1_lepplus_phi/F");
  tree->Branch("Z1_lepplus_etcone20", &(cand->Z1_lepplus_etcone20), "Z1_lepplus_etcone20/F");
  tree->Branch("Z1_lepplus_etcone20_corr", &(cand->Z1_lepplus_etcone20_corr), "Z1_lepplus_etcone20_corr/F");
  tree->Branch("Z1_lepplus_etcone20_final", &(cand->Z1_lepplus_etcone20_final), "Z1_lepplus_etcone20_final/F");
  tree->Branch("Z1_lepplus_ptcone20", &(cand->Z1_lepplus_ptcone20), "Z1_lepplus_ptcone20/F");
  tree->Branch("Z1_lepplus_ptcone20_final", &(cand->Z1_lepplus_ptcone20_final), "Z1_lepplus_ptcone20_final/F");
  tree->Branch("Z1_lepplus_d0", &(cand->Z1_lepplus_d0), "Z1_lepplus_d0/F");
  tree->Branch("Z1_lepplus_z0", &(cand->Z1_lepplus_z0), "Z1_lepplus_z0/F");
  tree->Branch("Z1_lepplus_d0_sig", &(cand->Z1_lepplus_d0_sig), "Z1_lepplus_d0_sig/F");
  tree->Branch("Z1_lepplus_GSF_dp", &(cand->Z1_lepplus_GSF_dp), "Z1_lepplus_GSF_dp/F");
  tree->Branch("Z1_lepplus_GSF_p", &(cand->Z1_lepplus_GSF_p), "Z1_lepplus_GSF_p/F");
  tree->Branch("Z1_lepplus_charge", &(cand->Z1_lepplus_charge), "Z1_lepplus_charge/F");
  tree->Branch("Z1_lepplus_ptCB_nosmearnoscale", &(cand->Z1_lepplus_ptCB_nosmearnoscale), "Z1_lepplus_ptCB_nosmearnoscale/F");
  tree->Branch("Z1_lepplus_ptME_nosmearnoscale", &(cand->Z1_lepplus_ptME_nosmearnoscale), "Z1_lepplus_ptME_nosmearnoscale/F");
  tree->Branch("Z1_lepplus_ptID_nosmearnoscale", &(cand->Z1_lepplus_ptID_nosmearnoscale), "Z1_lepplus_ptID_nosmearnoscale/F");
  tree->Branch("Z1_lepplus_cl_E_calibsmearnoscale", &(cand->Z1_lepplus_cl_E_calibsmearnoscale), "Z1_lepplus_cl_E_calibsmearnoscale/F");
  tree->Branch("Z1_lepplus_cl_eta", &(cand->Z1_lepplus_cl_eta), "Z1_lepplus_cl_eta/F");
  tree->Branch("Z1_lepplus_cl_phi", &(cand->Z1_lepplus_cl_phi), "Z1_lepplus_cl_phi/F");
  tree->Branch("Z1_lepplus_isSA", &(cand->Z1_lepplus_isSA), "Z1_lepplus_isSA/I");
  tree->Branch("Z1_lepplus_author", &(cand->Z1_lepplus_author), "Z1_lepplus_author/I");
  tree->Branch("Z1_lepplus_type", &(cand->Z1_lepplus_type), "Z1_lepplus_type/I");
  tree->Branch("Z1_lepplus_typebkg", &(cand->Z1_lepplus_typebkg), "Z1_lepplus_typebkg/I");
  tree->Branch("Z1_lepplus_origin", &(cand->Z1_lepplus_origin), "Z1_lepplus_origin/I");
  tree->Branch("Z1_lepplus_originbkg", &(cand->Z1_lepplus_originbkg), "Z1_lepplus_originbkg/I");
  tree->Branch("Z1_lepminus_m", &(cand->Z1_lepminus_m), "Z1_lepminus_m/F");
  tree->Branch("Z1_lepminus_pt", &(cand->Z1_lepminus_pt), "Z1_lepminus_pt/F");
  tree->Branch("Z1_lepminus_pt_truth", &(cand->Z1_lepminus_pt_truth), "Z1_lepminus_pt_truth/F");
  tree->Branch("Z1_lepminus_pt_constrained", &(cand->Z1_lepminus_pt_constrained), "Z1_lepminus_pt_constrained/F");
  tree->Branch("Z1_lepminus_eta", &(cand->Z1_lepminus_eta), "Z1_lepminus_eta/F");
  tree->Branch("Z1_lepminus_phi", &(cand->Z1_lepminus_phi), "Z1_lepminus_phi/F");
  tree->Branch("Z1_lepminus_etcone20", &(cand->Z1_lepminus_etcone20), "Z1_lepminus_etcone20/F");
  tree->Branch("Z1_lepminus_etcone20_corr", &(cand->Z1_lepminus_etcone20_corr), "Z1_lepminus_etcone20_corr/F");
  tree->Branch("Z1_lepminus_etcone20_final", &(cand->Z1_lepminus_etcone20_final), "Z1_lepminus_etcone20_final/F");
  tree->Branch("Z1_lepminus_ptcone20", &(cand->Z1_lepminus_ptcone20), "Z1_lepminus_ptcone20/F");
  tree->Branch("Z1_lepminus_ptcone20_final", &(cand->Z1_lepminus_ptcone20_final), "Z1_lepminus_ptcone20_final/F");
  tree->Branch("Z1_lepminus_d0", &(cand->Z1_lepminus_d0), "Z1_lepminus_d0/F");
  tree->Branch("Z1_lepminus_z0", &(cand->Z1_lepminus_z0), "Z1_lepminus_z0/F");
  tree->Branch("Z1_lepminus_d0_sig", &(cand->Z1_lepminus_d0_sig), "Z1_lepminus_d0_sig/F");
  tree->Branch("Z1_lepminus_GSF_dp", &(cand->Z1_lepminus_GSF_dp), "Z1_lepminus_GSF_dp/F");
  tree->Branch("Z1_lepminus_GSF_p", &(cand->Z1_lepminus_GSF_p), "Z1_lepminus_GSF_p/F");
  tree->Branch("Z1_lepminus_charge", &(cand->Z1_lepminus_charge), "Z1_lepminus_charge/F");
  tree->Branch("Z1_lepminus_ptCB_nosmearnoscale", &(cand->Z1_lepminus_ptCB_nosmearnoscale), "Z1_lepminus_ptCB_nosmearnoscale/F");
  tree->Branch("Z1_lepminus_ptME_nosmearnoscale", &(cand->Z1_lepminus_ptME_nosmearnoscale), "Z1_lepminus_ptME_nosmearnoscale/F");
  tree->Branch("Z1_lepminus_ptID_nosmearnoscale", &(cand->Z1_lepminus_ptID_nosmearnoscale), "Z1_lepminus_ptID_nosmearnoscale/F");
  tree->Branch("Z1_lepminus_cl_E_calibsmearnoscale", &(cand->Z1_lepminus_cl_E_calibsmearnoscale), "Z1_lepminus_cl_E_calibsmearnoscale/F");
  tree->Branch("Z1_lepminus_cl_eta", &(cand->Z1_lepminus_cl_eta), "Z1_lepminus_cl_eta/F");
  tree->Branch("Z1_lepminus_cl_phi", &(cand->Z1_lepminus_cl_phi), "Z1_lepminus_cl_phi/F");
  tree->Branch("Z1_lepminus_isSA", &(cand->Z1_lepminus_isSA), "Z1_lepminus_isSA/I");
  tree->Branch("Z1_lepminus_author", &(cand->Z1_lepminus_author), "Z1_lepminus_author/I");
  tree->Branch("Z1_lepminus_type", &(cand->Z1_lepminus_type), "Z1_lepminus_type/I");
  tree->Branch("Z1_lepminus_typebkg", &(cand->Z1_lepminus_typebkg), "Z1_lepminus_typebkg/I");
  tree->Branch("Z1_lepminus_origin", &(cand->Z1_lepminus_origin), "Z1_lepminus_origin/I");
  tree->Branch("Z1_lepminus_originbkg", &(cand->Z1_lepminus_originbkg), "Z1_lepminus_originbkg/I");
  tree->Branch("Z2_lepplus_m", &(cand->Z2_lepplus_m), "Z2_lepplus_m/F");
  tree->Branch("Z2_lepplus_pt", &(cand->Z2_lepplus_pt), "Z2_lepplus_pt/F");
  tree->Branch("Z2_lepplus_pt_truth", &(cand->Z2_lepplus_pt_truth), "Z2_lepplus_pt_truth/F");
  tree->Branch("Z2_lepplus_pt_constrained", &(cand->Z2_lepplus_pt_constrained), "Z2_lepplus_pt_constrained/F");
  tree->Branch("Z2_lepplus_eta", &(cand->Z2_lepplus_eta), "Z2_lepplus_eta/F");
  tree->Branch("Z2_lepplus_phi", &(cand->Z2_lepplus_phi), "Z2_lepplus_phi/F");
  tree->Branch("Z2_lepplus_etcone20", &(cand->Z2_lepplus_etcone20), "Z2_lepplus_etcone20/F");
  tree->Branch("Z2_lepplus_etcone20_corr", &(cand->Z2_lepplus_etcone20_corr), "Z2_lepplus_etcone20_corr/F");
  tree->Branch("Z2_lepplus_etcone20_final", &(cand->Z2_lepplus_etcone20_final), "Z2_lepplus_etcone20_final/F");
  tree->Branch("Z2_lepplus_ptcone20", &(cand->Z2_lepplus_ptcone20), "Z2_lepplus_ptcone20/F");
  tree->Branch("Z2_lepplus_ptcone20_final", &(cand->Z2_lepplus_ptcone20_final), "Z2_lepplus_ptcone20_final/F");
  tree->Branch("Z2_lepplus_d0", &(cand->Z2_lepplus_d0), "Z2_lepplus_d0/F");
  tree->Branch("Z2_lepplus_z0", &(cand->Z2_lepplus_z0), "Z2_lepplus_z0/F");
  tree->Branch("Z2_lepplus_d0_sig", &(cand->Z2_lepplus_d0_sig), "Z2_lepplus_d0_sig/F");
  tree->Branch("Z2_lepplus_GSF_dp", &(cand->Z2_lepplus_GSF_dp), "Z2_lepplus_GSF_dp/F");
  tree->Branch("Z2_lepplus_GSF_p", &(cand->Z2_lepplus_GSF_p), "Z2_lepplus_GSF_p/F");
  tree->Branch("Z2_lepplus_charge", &(cand->Z2_lepplus_charge), "Z2_lepplus_charge/F");
  tree->Branch("Z2_lepplus_ptCB_nosmearnoscale", &(cand->Z2_lepplus_ptCB_nosmearnoscale), "Z2_lepplus_ptCB_nosmearnoscale/F");
  tree->Branch("Z2_lepplus_ptME_nosmearnoscale", &(cand->Z2_lepplus_ptME_nosmearnoscale), "Z2_lepplus_ptME_nosmearnoscale/F");
  tree->Branch("Z2_lepplus_ptID_nosmearnoscale", &(cand->Z2_lepplus_ptID_nosmearnoscale), "Z2_lepplus_ptID_nosmearnoscale/F");
  tree->Branch("Z2_lepplus_cl_E_calibsmearnoscale", &(cand->Z2_lepplus_cl_E_calibsmearnoscale), "Z2_lepplus_cl_E_calibsmearnoscale/F");
  tree->Branch("Z2_lepplus_cl_eta", &(cand->Z2_lepplus_cl_eta), "Z2_lepplus_cl_eta/F");
  tree->Branch("Z2_lepplus_cl_phi", &(cand->Z2_lepplus_cl_phi), "Z2_lepplus_cl_phi/F");
  tree->Branch("Z2_lepplus_isSA", &(cand->Z2_lepplus_isSA), "Z2_lepplus_isSA/I");
  tree->Branch("Z2_lepplus_author", &(cand->Z2_lepplus_author), "Z2_lepplus_author/I");
  tree->Branch("Z2_lepplus_type", &(cand->Z2_lepplus_type), "Z2_lepplus_type/I");
  tree->Branch("Z2_lepplus_typebkg", &(cand->Z2_lepplus_typebkg), "Z2_lepplus_typebkg/I");
  tree->Branch("Z2_lepplus_origin", &(cand->Z2_lepplus_origin), "Z2_lepplus_origin/I");
  tree->Branch("Z2_lepplus_originbkg", &(cand->Z2_lepplus_originbkg), "Z2_lepplus_originbkg/I");
  tree->Branch("Z2_lepminus_m", &(cand->Z2_lepminus_m), "Z2_lepminus_m/F");
  tree->Branch("Z2_lepminus_pt", &(cand->Z2_lepminus_pt), "Z2_lepminus_pt/F");
  tree->Branch("Z2_lepminus_pt_truth", &(cand->Z2_lepminus_pt_truth), "Z2_lepminus_pt_truth/F");
  tree->Branch("Z2_lepminus_pt_constrained", &(cand->Z2_lepminus_pt_constrained), "Z2_lepminus_pt_constrained/F");
  tree->Branch("Z2_lepminus_eta", &(cand->Z2_lepminus_eta), "Z2_lepminus_eta/F");
  tree->Branch("Z2_lepminus_phi", &(cand->Z2_lepminus_phi), "Z2_lepminus_phi/F");
  tree->Branch("Z2_lepminus_etcone20", &(cand->Z2_lepminus_etcone20), "Z2_lepminus_etcone20/F");
  tree->Branch("Z2_lepminus_etcone20_corr", &(cand->Z2_lepminus_etcone20_corr), "Z2_lepminus_etcone20_corr/F");
  tree->Branch("Z2_lepminus_etcone20_final", &(cand->Z2_lepminus_etcone20_final), "Z2_lepminus_etcone20_final/F");
  tree->Branch("Z2_lepminus_ptcone20", &(cand->Z2_lepminus_ptcone20), "Z2_lepminus_ptcone20/F");
  tree->Branch("Z2_lepminus_ptcone20_final", &(cand->Z2_lepminus_ptcone20_final), "Z2_lepminus_ptcone20_final/F");
  tree->Branch("Z2_lepminus_d0", &(cand->Z2_lepminus_d0), "Z2_lepminus_d0/F");
  tree->Branch("Z2_lepminus_z0", &(cand->Z2_lepminus_z0), "Z2_lepminus_z0/F");
  tree->Branch("Z2_lepminus_d0_sig", &(cand->Z2_lepminus_d0_sig), "Z2_lepminus_d0_sig/F");
  tree->Branch("Z2_lepminus_GSF_dp", &(cand->Z2_lepminus_GSF_dp), "Z2_lepminus_GSF_dp/F");
  tree->Branch("Z2_lepminus_GSF_p", &(cand->Z2_lepminus_GSF_p), "Z2_lepminus_GSF_p/F");
  tree->Branch("Z2_lepminus_charge", &(cand->Z2_lepminus_charge), "Z2_lepminus_charge/F");
  tree->Branch("Z2_lepminus_ptCB_nosmearnoscale", &(cand->Z2_lepminus_ptCB_nosmearnoscale), "Z2_lepminus_ptCB_nosmearnoscale/F");
  tree->Branch("Z2_lepminus_ptME_nosmearnoscale", &(cand->Z2_lepminus_ptME_nosmearnoscale), "Z2_lepminus_ptME_nosmearnoscale/F");
  tree->Branch("Z2_lepminus_ptID_nosmearnoscale", &(cand->Z2_lepminus_ptID_nosmearnoscale), "Z2_lepminus_ptID_nosmearnoscale/F");
  tree->Branch("Z2_lepminus_cl_E_calibsmearnoscale", &(cand->Z2_lepminus_cl_E_calibsmearnoscale), "Z2_lepminus_cl_E_calibsmearnoscale/F");
  tree->Branch("Z2_lepminus_cl_eta", &(cand->Z2_lepminus_cl_eta), "Z2_lepminus_cl_eta/F");
  tree->Branch("Z2_lepminus_cl_phi", &(cand->Z2_lepminus_cl_phi), "Z2_lepminus_cl_phi/F");
  tree->Branch("Z2_lepminus_isSA", &(cand->Z2_lepminus_isSA), "Z2_lepminus_isSA/I");
  tree->Branch("Z2_lepminus_author", &(cand->Z2_lepminus_author), "Z2_lepminus_author/I");
  tree->Branch("Z2_lepminus_type", &(cand->Z2_lepminus_type), "Z2_lepminus_type/I");
  tree->Branch("Z2_lepminus_typebkg", &(cand->Z2_lepminus_typebkg), "Z2_lepminus_typebkg/I");
  tree->Branch("Z2_lepminus_origin", &(cand->Z2_lepminus_origin), "Z2_lepminus_origin/I");
  tree->Branch("Z2_lepminus_originbkg", &(cand->Z2_lepminus_originbkg), "Z2_lepminus_originbkg/I");
  tree->Branch("filename", &(cand->filename));
}

void HiggsllqqAnalysis::fillCandidateStruct(Analysis::CandidateStruct * output, Analysis::Quadrilepton * higgs)
{
  // fill output->last with the last cut passed by the quadrilepton
  // independently on wether it was selected or not
  // output->closesttozz is then filled with a flag which tells if the candidate would
  // have been selected since it's the closest to ZZ among those passing the DeltaR cut
  output->last = higgs->lastcut();
  output->closesttozz = (isClosestToZZ(higgs)) ? 1 : 0;
  
  // check if this quadrilepton is among the GoodQuadrileptons
  output->selected = (isSelected(higgs)) ? 1 : 0;
  
  // fill other informations
  output->type = higgs->channel();
  output->run = (isMC()) ? ntuple->eventinfo.mc_channel_number() : ntuple->eventinfo.RunNumber();
  output->event = ntuple->eventinfo.EventNumber();
  output->lbn = ntuple->eventinfo.lbn();
  output->n_vx = getNumberOfGoodVertices();
  output->actualIntPerXing = ntuple->eventinfo.actualIntPerXing();
  output->averageIntPerXing = ntuple->eventinfo.averageIntPerXing();
  output->hfor = (isMC()) ? ntuple->top.hfor_type() : 0;
  output->top_weight = (isWithinT1llqqPhaseSpace()) ? 0 : 1; // multiplicative weight
  output->powhegbug_weight = (hasPowHegZZBug()) ? 0 : 1; // multiplicative weight
  
  output->xsec_weight = getCrossSectionWeight();
  output->processed_entries = m_entriesInChain;
  output->generated_entries = -9999;//no MC config file available, feature removed //(isMC()) ? m_PileupReweighter->GetNumberOfEvents(ntuple->eventinfo.mc_channel_number()) : -9999;
  output->pu_weight = getEventWeight();
  output->vxz_weight = getVertexZWeight();
  output->ggF_weight = getggFWeight();
  output->trigSF_weight = getCandidateTriggerSF(higgs);
  output->Z1_lepplus_weight = getLeptonWeight(higgs->GetZ1()->GetLepPlus());
  output->Z1_lepminus_weight = getLeptonWeight(higgs->GetZ1()->GetLepMinus());
  output->Z2_lepplus_weight = getLeptonWeight(higgs->GetZ2()->GetLepPlus());
  output->Z2_lepminus_weight = getLeptonWeight(higgs->GetZ2()->GetLepMinus());
  
  output->H_m_truth = getTruthHiggsMass();
  output->ZZ_m_truth = getTruthZZMass();
  
  //std::pair<Float_t, Float_t> truthZZ = getTruthZMass(higgs);
  //output->Z1_m_truth = truthZZ.first;
  //output->Z2_m_truth = truthZZ.second;
  output->Z1_m_truth = -1;
  output->Z2_m_truth = -1;
  output->Z1_lepplus_pt_truth = -1;
  output->Z1_lepminus_pt_truth = -1;
  output->Z2_lepplus_pt_truth = -1;
  output->Z2_lepminus_pt_truth = -1;
  std::vector<TLorentzVector *> truthLeptons = getTruthZLeptons(higgs);
  if (truthLeptons[0] && truthLeptons[1]) output->Z1_m_truth = (*truthLeptons[0] + *truthLeptons[1]).M();
  if (truthLeptons[2] && truthLeptons[3]) output->Z2_m_truth = (*truthLeptons[2] + *truthLeptons[3]).M();
  if (truthLeptons[0]) output->Z1_lepplus_pt_truth = truthLeptons[0]->Pt();
  if (truthLeptons[1]) output->Z1_lepminus_pt_truth = truthLeptons[1]->Pt();
  if (truthLeptons[2]) output->Z2_lepplus_pt_truth = truthLeptons[2]->Pt();
  if (truthLeptons[3]) output->Z2_lepminus_pt_truth = truthLeptons[3]->Pt();
  for (Int_t i = 0; i < 4; i++) if (truthLeptons[i]) delete truthLeptons[i];
  
  Analysis::ConstraintFitResult constraint_fit = getMassConstraintFit(higgs);
  output->H_m_constrained = constraint_fit.H_4m_updated.M();
  output->Z1_m_constrained = constraint_fit.Z1_4m_updated.M();
  output->Z2_m_constrained = constraint_fit.Z2_4m_updated.M();
  output->Z1_chiSq = constraint_fit.Z1_chiSq;
  output->Z2_chiSq = constraint_fit.Z2_chiSq;
  output->Z1_lepplus_pt_constrained = constraint_fit.Z1_lepplus_4m_updated.Pt();
  output->Z1_lepminus_pt_constrained = constraint_fit.Z1_lepminus_4m_updated.Pt();
  output->Z2_lepplus_pt_constrained = constraint_fit.Z2_lepplus_4m_updated.Pt();
  output->Z2_lepminus_pt_constrained = constraint_fit.Z2_lepminus_4m_updated.Pt();
  
  output->H_m = higgs->Get4Momentum()->M();
  output->H_pt = higgs->Get4Momentum()->Pt();
  output->H_eta = higgs->Get4Momentum()->Eta();
  output->H_phi = higgs->Get4Momentum()->Phi();
  
  output->Z1_m = higgs->GetZ1()->Get4Momentum()->M();
  output->Z1_pt = higgs->GetZ1()->Get4Momentum()->Pt();
  output->Z1_eta = higgs->GetZ1()->Get4Momentum()->Eta();
  output->Z1_phi = higgs->GetZ1()->Get4Momentum()->Phi();
  
  output->Z2_m = higgs->GetZ2()->Get4Momentum()->M();
  output->Z2_pt = higgs->GetZ2()->Get4Momentum()->Pt();
  output->Z2_eta = higgs->GetZ2()->Get4Momentum()->Eta();
  output->Z2_phi = higgs->GetZ2()->Get4Momentum()->Phi();
  
  std::vector<Float_t> trackIsolationVector = getFinalTrackIsoVector(higgs);
  std::vector<Float_t> caloIsolationVector = getFinalCaloIsoVector(higgs);
  HllqqSystematics::ChargedLepton Z1_lepplus_syst = getLeptonSystematics(higgs->GetZ1()->GetLepPlus());
  HllqqSystematics::ChargedLepton Z1_lepminus_syst = getLeptonSystematics(higgs->GetZ1()->GetLepMinus());
  HllqqSystematics::ChargedLepton Z2_lepplus_syst = getLeptonSystematics(higgs->GetZ2()->GetLepPlus());
  HllqqSystematics::ChargedLepton Z2_lepminus_syst = getLeptonSystematics(higgs->GetZ2()->GetLepMinus());
  
  output->Z1_lepplus_m = higgs->GetZ1()->GetLepPlus()->Get4Momentum()->M();
  output->Z1_lepplus_pt = higgs->GetZ1()->GetLepPlus()->Get4Momentum()->Pt();
  output->Z1_lepplus_eta = higgs->GetZ1()->GetLepPlus()->Get4Momentum()->Eta();
  output->Z1_lepplus_phi = higgs->GetZ1()->GetLepPlus()->Get4Momentum()->Phi();
  output->Z1_lepplus_etcone20 = higgs->GetZ1()->GetLepPlus()->etcone20();
  output->Z1_lepplus_etcone20_corr = getCorrectedCaloIso(higgs->GetZ1()->GetLepPlus());
  output->Z1_lepplus_etcone20_final = caloIsolationVector[0];
  output->Z1_lepplus_ptcone20 = higgs->GetZ1()->GetLepPlus()->ptcone20();
  output->Z1_lepplus_ptcone20_final = trackIsolationVector[0];
  output->Z1_lepplus_d0 = higgs->GetZ1()->GetLepPlus()->d0();
  
  if (isMC() && analysis_version() == "rel_17_2") output->Z1_lepplus_d0 -= 0.002;
  output->Z1_lepplus_z0 = higgs->GetZ1()->GetLepPlus()->z0();
  output->Z1_lepplus_d0_sig = higgs->GetZ1()->GetLepPlus()->d0_sig();
  output->Z1_lepplus_GSF_dp = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? CommonTools::getBremFitDp(higgs->GetZ1()->GetLepPlus()->GetElectron()) : -9999;
  output->Z1_lepplus_GSF_p = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? 1. / higgs->GetZ1()->GetLepPlus()->GetElectron()->trackqoverp() : -9999;
  output->Z1_lepplus_charge = higgs->GetZ1()->GetLepPlus()->charge();
  output->Z1_lepplus_ptCB_nosmearnoscale = Z1_lepplus_syst.ptCB_nosmearnoscale;
  output->Z1_lepplus_ptME_nosmearnoscale = Z1_lepplus_syst.ptME_nosmearnoscale;
  output->Z1_lepplus_ptID_nosmearnoscale = Z1_lepplus_syst.ptID_nosmearnoscale;
  output->Z1_lepplus_cl_E_calibsmearnoscale = Z1_lepplus_syst.cl_E_calibsmearnoscale;
  output->Z1_lepplus_cl_eta = Z1_lepplus_syst.cl_eta;
  output->Z1_lepplus_cl_phi = Z1_lepplus_syst.cl_phi;
  output->Z1_lepplus_isSA = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? -9999 : higgs->GetZ1()->GetLepPlus()->GetMuon()->isStandAloneMuon();
  output->Z1_lepplus_author = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? higgs->GetZ1()->GetLepPlus()->GetElectron()->author() : higgs->GetZ1()->GetLepPlus()->GetMuon()->author();
  output->Z1_lepplus_type = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepPlus()->GetElectron()->type() : -9999;
  output->Z1_lepplus_typebkg = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepPlus()->GetElectron()->typebkg() : -9999;
  output->Z1_lepplus_origin = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepPlus()->GetElectron()->origin() : -9999;
  output->Z1_lepplus_originbkg = (higgs->GetZ1()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepPlus()->GetElectron()->originbkg() : -9999;
  output->Z1_lepminus_m = higgs->GetZ1()->GetLepMinus()->Get4Momentum()->M();
  output->Z1_lepminus_pt = higgs->GetZ1()->GetLepMinus()->Get4Momentum()->Pt();
  output->Z1_lepminus_eta = higgs->GetZ1()->GetLepMinus()->Get4Momentum()->Eta();
  output->Z1_lepminus_phi = higgs->GetZ1()->GetLepMinus()->Get4Momentum()->Phi();
  output->Z1_lepminus_etcone20 = higgs->GetZ1()->GetLepMinus()->etcone20();
  output->Z1_lepminus_etcone20_corr = getCorrectedCaloIso(higgs->GetZ1()->GetLepMinus());
  output->Z1_lepminus_etcone20_final = caloIsolationVector[1];
  output->Z1_lepminus_ptcone20 = higgs->GetZ1()->GetLepMinus()->ptcone20();
  output->Z1_lepminus_ptcone20_final = trackIsolationVector[1];
  output->Z1_lepminus_d0 = higgs->GetZ1()->GetLepMinus()->d0();
  
  if (isMC() && analysis_version() == "rel_17_2") output->Z1_lepminus_d0 -= 0.002;
  output->Z1_lepminus_z0 = higgs->GetZ1()->GetLepMinus()->z0();
  output->Z1_lepminus_d0_sig = higgs->GetZ1()->GetLepMinus()->d0_sig();
  output->Z1_lepminus_GSF_dp = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? CommonTools::getBremFitDp(higgs->GetZ1()->GetLepMinus()->GetElectron()) : -9999;
  output->Z1_lepminus_GSF_p = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? 1. / higgs->GetZ1()->GetLepMinus()->GetElectron()->trackqoverp() : -9999;
  output->Z1_lepminus_charge = higgs->GetZ1()->GetLepMinus()->charge();
  output->Z1_lepminus_ptCB_nosmearnoscale = Z1_lepminus_syst.ptCB_nosmearnoscale;
  output->Z1_lepminus_ptME_nosmearnoscale = Z1_lepminus_syst.ptME_nosmearnoscale;
  output->Z1_lepminus_ptID_nosmearnoscale = Z1_lepminus_syst.ptID_nosmearnoscale;
  output->Z1_lepminus_cl_E_calibsmearnoscale = Z1_lepminus_syst.cl_E_calibsmearnoscale;
  output->Z1_lepminus_cl_eta = Z1_lepminus_syst.cl_eta;
  output->Z1_lepminus_cl_phi = Z1_lepminus_syst.cl_phi;
  output->Z1_lepminus_isSA = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? -9999  : higgs->GetZ1()->GetLepMinus()->GetMuon()->isStandAloneMuon();
  output->Z1_lepminus_author = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? higgs->GetZ1()->GetLepMinus()->GetElectron()->author()  :  higgs->GetZ1()->GetLepMinus()->GetMuon()->author();
  output->Z1_lepminus_type = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepMinus()->GetElectron()->type() : -9999;
  output->Z1_lepminus_typebkg = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepMinus()->GetElectron()->typebkg() : -9999;
  output->Z1_lepminus_origin = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepMinus()->GetElectron()->origin() : -9999;
  output->Z1_lepminus_originbkg = (higgs->GetZ1()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ1()->GetLepMinus()->GetElectron()->originbkg() : -9999;
  output->Z2_lepplus_m = higgs->GetZ2()->GetLepPlus()->Get4Momentum()->M();
  output->Z2_lepplus_pt = higgs->GetZ2()->GetLepPlus()->Get4Momentum()->Pt();
  output->Z2_lepplus_eta = higgs->GetZ2()->GetLepPlus()->Get4Momentum()->Eta();
  output->Z2_lepplus_phi = higgs->GetZ2()->GetLepPlus()->Get4Momentum()->Phi();
  output->Z2_lepplus_etcone20 = higgs->GetZ2()->GetLepPlus()->etcone20();
  output->Z2_lepplus_etcone20_corr = getCorrectedCaloIso(higgs->GetZ2()->GetLepPlus());
  output->Z2_lepplus_etcone20_final = caloIsolationVector[2];
  output->Z2_lepplus_ptcone20 = higgs->GetZ2()->GetLepPlus()->ptcone20();
  output->Z2_lepplus_ptcone20_final = trackIsolationVector[2];
  output->Z2_lepplus_d0 = higgs->GetZ2()->GetLepPlus()->d0();
  
  if (isMC() && analysis_version() == "rel_17_2") output->Z2_lepplus_d0 -= 0.002;
  output->Z2_lepplus_z0 = higgs->GetZ2()->GetLepPlus()->z0();
  output->Z2_lepplus_d0_sig = higgs->GetZ2()->GetLepPlus()->d0_sig();
  output->Z2_lepplus_GSF_dp = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? CommonTools::getBremFitDp(higgs->GetZ2()->GetLepPlus()->GetElectron()) : -9999;
  output->Z2_lepplus_GSF_p = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? 1. / higgs->GetZ2()->GetLepPlus()->GetElectron()->trackqoverp() : -9999;
  output->Z2_lepplus_charge = higgs->GetZ2()->GetLepPlus()->charge();
  output->Z2_lepplus_ptCB_nosmearnoscale = Z2_lepplus_syst.ptCB_nosmearnoscale;
  output->Z2_lepplus_ptME_nosmearnoscale = Z2_lepplus_syst.ptME_nosmearnoscale;
  output->Z2_lepplus_ptID_nosmearnoscale = Z2_lepplus_syst.ptID_nosmearnoscale;
  output->Z2_lepplus_cl_E_calibsmearnoscale = Z2_lepplus_syst.cl_E_calibsmearnoscale;
  output->Z2_lepplus_cl_eta = Z2_lepplus_syst.cl_eta;
  output->Z2_lepplus_cl_phi = Z2_lepplus_syst.cl_phi;
  output->Z2_lepplus_isSA = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? -9999 : higgs->GetZ2()->GetLepPlus()->GetMuon()->isStandAloneMuon();
  output->Z2_lepplus_author = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? higgs->GetZ2()->GetLepPlus()->GetElectron()->author() : higgs->GetZ2()->GetLepPlus()->GetMuon()->author();
  output->Z2_lepplus_type = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepPlus()->GetElectron()->type() : -9999;
  output->Z2_lepplus_typebkg = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepPlus()->GetElectron()->typebkg() : -9999;
  output->Z2_lepplus_origin = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepPlus()->GetElectron()->origin() : -9999;
  output->Z2_lepplus_originbkg = (higgs->GetZ2()->GetLepPlus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepPlus()->GetElectron()->originbkg() : -9999;
  output->Z2_lepminus_m = higgs->GetZ2()->GetLepMinus()->Get4Momentum()->M();
  output->Z2_lepminus_pt = higgs->GetZ2()->GetLepMinus()->Get4Momentum()->Pt();
  output->Z2_lepminus_eta = higgs->GetZ2()->GetLepMinus()->Get4Momentum()->Eta();
  output->Z2_lepminus_phi = higgs->GetZ2()->GetLepMinus()->Get4Momentum()->Phi();
  output->Z2_lepminus_etcone20 = higgs->GetZ2()->GetLepMinus()->etcone20();
  output->Z2_lepminus_etcone20_corr = getCorrectedCaloIso(higgs->GetZ2()->GetLepMinus());
  output->Z2_lepminus_etcone20_final = caloIsolationVector[3];
  output->Z2_lepminus_ptcone20 = higgs->GetZ2()->GetLepMinus()->ptcone20();
  output->Z2_lepminus_ptcone20_final = trackIsolationVector[3];
  output->Z2_lepminus_d0 = higgs->GetZ2()->GetLepMinus()->d0();
   
  if (isMC() && analysis_version() == "rel_17_2") output->Z2_lepminus_d0 -= 0.002;
  output->Z2_lepminus_z0 = higgs->GetZ2()->GetLepMinus()->z0();
  output->Z2_lepminus_d0_sig = higgs->GetZ2()->GetLepMinus()->d0_sig();
  output->Z2_lepminus_GSF_dp = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? CommonTools::getBremFitDp(higgs->GetZ2()->GetLepMinus()->GetElectron()) : -9999;
  output->Z2_lepminus_GSF_p = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? 1. / higgs->GetZ2()->GetLepMinus()->GetElectron()->trackqoverp() : -9999;
  output->Z2_lepminus_charge = higgs->GetZ2()->GetLepMinus()->charge();
  output->Z2_lepminus_ptCB_nosmearnoscale = Z2_lepminus_syst.ptCB_nosmearnoscale;
  output->Z2_lepminus_ptME_nosmearnoscale = Z2_lepminus_syst.ptME_nosmearnoscale;
  output->Z2_lepminus_ptID_nosmearnoscale = Z2_lepminus_syst.ptID_nosmearnoscale;
  output->Z2_lepminus_cl_E_calibsmearnoscale = Z2_lepminus_syst.cl_E_calibsmearnoscale;
  output->Z2_lepminus_cl_eta = Z2_lepminus_syst.cl_eta;
  output->Z2_lepminus_cl_phi = Z2_lepminus_syst.cl_phi;
  output->Z2_lepminus_isSA = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? -9999 : higgs->GetZ2()->GetLepMinus()->GetMuon()->isStandAloneMuon();
  output->Z2_lepminus_author = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON) ? higgs->GetZ2()->GetLepMinus()->GetElectron()->author() : higgs->GetZ2()->GetLepMinus()->GetMuon()->author();
  output->Z2_lepminus_type = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepMinus()->GetElectron()->type() : -9999;
  output->Z2_lepminus_typebkg = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepMinus()->GetElectron()->typebkg() : -9999;
  output->Z2_lepminus_origin = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepMinus()->GetElectron()->origin() : -9999;
  output->Z2_lepminus_originbkg = (higgs->GetZ2()->GetLepMinus()->flavor() == Analysis::ChargedLepton::ELECTRON && isMC()) ? higgs->GetZ2()->GetLepMinus()->GetElectron()->originbkg() : -9999;
  output->filename.clear();
  output->filename.push_back(std::string(fEventTree->GetCurrentFile()->GetName()));
}


Float_t HiggsllqqAnalysis::getggFWeight()
{
  Float_t result = 1.;
  
  if (isMC()) {
    // check if this is a signal PowHeg ggF sample
    if (m_SignalSampleMass.find(ntuple->eventinfo.mc_channel_number()) != m_SignalSampleMass.end()) {
      result = m_ggFReweighter->getWeight(getTruthHiggsPt() / 1000);
    }
  }
  
  return result;
}

Float_t HiggsllqqAnalysis::getTruthHiggsPt()
{
  Float_t result(-1);
  if (isMC()) {
    for (Int_t i = 0; i < ntuple->mc.n(); i++) {
      if (TMath::Abs(ntuple->mc[i].pdgId()) == 25) {
	result = ntuple->mc[i].pt();
      }
    }
  }
  
  return result;
}

Float_t HiggsllqqAnalysis::getTruthHiggsMass()
{
  // return the mass of the first Higgs in the truth bank
  
  Float_t result(-1);
  
  if (isMC()) {
    for (Int_t i = 0; i < ntuple->mc.n(); i++) {
      if (TMath::Abs(ntuple->mc[i].pdgId()) == 25) {
	if (result == -1)
	  result = ntuple->mc[i].m();
      }
    }
  }
  
  return result;
}

std::vector<TLorentzVector *> HiggsllqqAnalysis::getTruthZLeptons(Analysis::Quadrilepton * higgs) {
  // return, for each reco lepton, a pointer to a TLorentzVector representing the closest truth lepton, if it comes from Z (otherwise, null pointer)
  std::vector<TLorentzVector *> result;
  for (Int_t i = 0; i < 4; i++) result.push_back(0);
  
  std::vector<Analysis::ChargedLepton *> reco_leptons = higgs->GetLeptonVector();
  std::vector<D3PDReader::TruthParticleD3PDObjectElement *> truth_leptons;
  
  Double_t dr_cut(0.1);
  
  if (isMC()) {
    for (UInt_t i = 0; i < reco_leptons.size(); i++) {
      truth_leptons.push_back(0);
      
      Double_t best_dr(99999.);
      
      for (Int_t j = 0; j < ntuple->mc.n(); j++) {
	Int_t pdgid = TMath::Abs(ntuple->mc[j].pdgId());
	
	// take leptons
	if (ntuple->mc[j].status() == 1 && (pdgid == 11 || pdgid == 13 || pdgid == 15)) {
	  // find the closest
	  Double_t this_dr = reco_leptons[i]->Get4Momentum()->DeltaR(CommonTools::getVector(&(ntuple->mc[j])));
	  
	  if (this_dr  < dr_cut && this_dr < best_dr) {
	    truth_leptons[i] = &(ntuple->mc[j]);
	    best_dr = this_dr;
	  } // best match
	} // stable charged lepton
      } // loop over truth
      
      if (truth_leptons[i]) {
	// match found
	if (truth_leptons[i]->parent_index().size() > 0) {
	  Int_t index_parent = truth_leptons[i]->parent_index().at(0);
	  
	  if (ntuple->mc[index_parent].pdgId() == 23) {
	    result[i] = new TLorentzVector(CommonTools::getVector(truth_leptons[i]));
	  } // parent is Z!
	} // non-empty parent list
      } // reco-truth match found
    } // loop over reco
  } // MC
  
  
  return result;
}

std::pair<Float_t, Float_t> HiggsllqqAnalysis::getTruthZMass(Analysis::Quadrilepton * higgs)
{
  // return the mass of the truth leptons associated to the reco ones (if they come from Z of course, but any Z, so accounts for mis-pairing
  std::pair<Float_t, Float_t> result = std::make_pair<Float_t, Float_t>(-9999.9, -9999.9);
  
  // first two Z bosons in the event
  D3PDReader::TruthParticleD3PDObjectElement *truth_Z_a = 0;
  D3PDReader::TruthParticleD3PDObjectElement *truth_Z_b = 0;
  
  Int_t index_Z_a(-1);
  Int_t index_Z_b(-1);
  
  if (isMC()) {
    for (Int_t i = 0; i < ntuple->mc.n(); i++) {
      if (TMath::Abs(ntuple->mc[i].pdgId()) == 23) {
	if (truth_Z_a == 0) {
	  truth_Z_a = &(ntuple->mc[i]);
	  index_Z_a = i;
	} else if (truth_Z_b == 0) {
	  truth_Z_b = &(ntuple->mc[i]);
	  index_Z_b = i;
	}
      }
    }
  }
  
  D3PDReader::TruthParticleD3PDObjectElement *truth_Z_a_lep1 = 0;
  D3PDReader::TruthParticleD3PDObjectElement *truth_Z_a_lep2 = 0;
  D3PDReader::TruthParticleD3PDObjectElement *truth_Z_b_lep1 = 0;
  D3PDReader::TruthParticleD3PDObjectElement *truth_Z_b_lep2 = 0;
  
  if (isMC()) {
    for (Int_t i = 0; i < ntuple->mc.n(); i++) {
      Int_t pdgid = TMath::Abs(ntuple->mc[i].pdgId());
      
      // be sure they are leptons
      if (pdgid == 11 || pdgid == 13 || pdgid == 15) {
	if (ntuple->mc[i].parent_index().size() > 0) {
	  if (ntuple->mc[i].parent_index().at(0) == index_Z_a) {
	    if (truth_Z_a_lep1 == 0)
	      truth_Z_a_lep1 = &(ntuple->mc[i]);
	    else if (truth_Z_a_lep2 == 0)
	      truth_Z_a_lep2 = &(ntuple->mc[i]);
	  } else if (ntuple->mc[i].parent_index().at(0) == index_Z_b) {
	    if (truth_Z_b_lep1 == 0)
	      truth_Z_b_lep1 = &(ntuple->mc[i]);
	    else if (truth_Z_b_lep2 == 0)
	      truth_Z_b_lep2 = &(ntuple->mc[i]);
	  }
	}
      }
    }
  }
  
  Float_t dr_cut(0.1);
  
  Bool_t Z1_lepplus_matched_Z_a(kFALSE);
  Bool_t Z1_lepminus_matched_Z_a(kFALSE);
  Bool_t Z1_lepplus_matched_Z_b(kFALSE);
  Bool_t Z1_lepminus_matched_Z_b(kFALSE);
  
  Bool_t Z2_lepplus_matched_Z_a(kFALSE);
  Bool_t Z2_lepminus_matched_Z_a(kFALSE);
  Bool_t Z2_lepplus_matched_Z_b(kFALSE);
  Bool_t Z2_lepminus_matched_Z_b(kFALSE);
  
  if (truth_Z_a_lep1 && truth_Z_a_lep2 && truth_Z_b_lep1 && truth_Z_b_lep2) {
    Z1_lepplus_matched_Z_a = (higgs->GetZ1()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep1))  < dr_cut || higgs->GetZ1()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep2)) < dr_cut);
    Z1_lepminus_matched_Z_a = (higgs->GetZ1()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep1))  < dr_cut || higgs->GetZ1()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep2)) < dr_cut);
    Z1_lepplus_matched_Z_b = (higgs->GetZ1()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep1))  < dr_cut || higgs->GetZ1()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep2)) < dr_cut);
    Z1_lepminus_matched_Z_b = (higgs->GetZ1()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep1))  < dr_cut || higgs->GetZ1()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep2)) < dr_cut);
    Z2_lepplus_matched_Z_a = (higgs->GetZ2()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep1))  < dr_cut || higgs->GetZ2()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep2)) < dr_cut);
    Z2_lepminus_matched_Z_a = (higgs->GetZ2()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep1))  < dr_cut || higgs->GetZ2()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_a_lep2)) < dr_cut);
    Z2_lepplus_matched_Z_b = (higgs->GetZ2()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep1))  < dr_cut || higgs->GetZ2()->GetLepPlus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep2)) < dr_cut);
    Z2_lepminus_matched_Z_b = (higgs->GetZ2()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep1))  < dr_cut || higgs->GetZ2()->GetLepMinus()->Get4Momentum()->DeltaR(CommonTools::getVector(truth_Z_b_lep2)) < dr_cut);
  }
  
  result.first  = -9999999.;
  result.second = -9999999.;
  
  if (Z1_lepplus_matched_Z_a && Z1_lepminus_matched_Z_a) {
    result.first = truth_Z_a->m();
  }
  if (Z2_lepplus_matched_Z_a && Z2_lepminus_matched_Z_a) {
    result.second = truth_Z_a->m();
  }
  if (Z1_lepplus_matched_Z_b && Z1_lepminus_matched_Z_b) {
    result.first = truth_Z_b->m();
  }
  if (Z2_lepplus_matched_Z_b && Z2_lepminus_matched_Z_b) {
    result.second = truth_Z_b->m();
  }
  
  return result;
}

Analysis::ConstraintFitResult HiggsllqqAnalysis::getMassConstraintFit(Analysis::Quadrilepton * higgs)
{
  Analysis::ConstraintFitResult result(higgs);
  
  // adapted from https://svnweb.cern.ch/trac/atlasusr/browser/knikolop/analysis/trunk/analysis.c
  double quadrupletMass = higgs->Get4Momentum()->M();
  
  // define track parameters and uncertainties over the parameters (covariance matrix)
  const int ntracks = 4;
  double tracks[ntracks][5] = {0.};
  double etracks[ntracks][5][5] = {0.};
  double particle_mass[ntracks];
  
  // fill the vector of Leptons properly
  std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
  // fill track parameters and covariance matrix
  for (int i = 0; i < ntracks; i++) {
    // for convenience, save phi, theta and P
    double phi   = (leptons[i]->flavor() == Analysis::ChargedLepton::MUON) ? leptons[i]->GetMuon()->phi_exPV() : leptons[i]->GetElectron()->trackphi();
    double theta = (leptons[i]->flavor() == Analysis::ChargedLepton::MUON) ? leptons[i]->GetMuon()->theta_exPV() : leptons[i]->GetElectron()->tracktheta();
    double P     = leptons[i]->Get4Momentum()->P();
    particle_mass[i] = leptons[i]->Get4Momentum()->M();
    
    // fill the track
    tracks[i][0] = (leptons[i]->flavor() == Analysis::ChargedLepton::MUON) ? leptons[i]->GetMuon()->d0_exPV() : leptons[i]->GetElectron()->trackd0();
    tracks[i][1] = (leptons[i]->flavor() == Analysis::ChargedLepton::MUON) ? leptons[i]->GetMuon()->z0_exPV() : leptons[i]->GetElectron()->trackz0();
    tracks[i][2] = leptons[i]->Get4Momentum()->Px();
    tracks[i][3] = leptons[i]->Get4Momentum()->Py();
    tracks[i][4] = leptons[i]->Get4Momentum()->Pz();
    
    // fill its covariance matrix
    HepMatrix covmatrix = CommonTools::getCovarianceMatrix(m_ElectronEnergyRescaler, leptons[i]);
    
    //going from d0,z0,phi,theta,P --> d0,z0,px,py,pz
    HepMatrix Jacobian(5, 5, 0);
    Jacobian(1, 1) = 1.;
    Jacobian(2, 2) = 1.;
    Jacobian(3, 3) = -P * TMath::Sin(theta) * TMath::Sin(phi);
    Jacobian(3, 4) =  P * TMath::Sin(theta) * TMath::Cos(phi);
    Jacobian(4, 3) =  P * TMath::Cos(theta) * TMath::Cos(phi);
    Jacobian(4, 4) =  P * TMath::Cos(theta) * TMath::Sin(phi);
    Jacobian(4, 5) = -P * TMath::Sin(theta); // !!!
    Jacobian(5, 3) =      TMath::Sin(theta) * TMath::Cos(phi);
    Jacobian(5, 4) =      TMath::Sin(theta) * TMath::Sin(phi);
    Jacobian(5, 5) =      TMath::Cos(theta); // !!!
    
    HepMatrix newcovariance(Jacobian.T() * covmatrix * Jacobian);
    
    // save the covariance matrix
    for (int j = 0; j < 5; j++)
      for (int k = 0; k < 5; k++)
	etracks[i][j][k] = newcovariance[j][k];
  }
  
  // fill Z1 fit parameters and uncertainties with track parameters and covariance matrix elements
  double parameters[2][4] = {0.};
  double sigma[6][6] = {0.};
  for (Int_t i = 0; i < 2; i++) {
    for (Int_t j = 0; j < 3; j++)
      parameters[i][j] = tracks[i][j + 2]; //(*FittedTracks[i])(j+1,1);
    parameters[i][3] = particle_mass[i];
  }
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
	sigma[3 * i + j][3 * i + k] = etracks[i][j + 2][k + 2];
  
  // do the same for Z2
  double parameters2[2][4] = {0.};
  double sigma2[6][6] = {0.};
  for (Int_t i = 0; i < 2; i++) {
    for (Int_t j = 0; j < 3; j++) {
      parameters2[i][j] =  tracks[i + 2][j + 2]; //(*FittedTracks[i])(j+1,1);
    }
    parameters2[i][3] = particle_mass[i];
  }
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 3; j++)
      for (int k = 0; k < 3; k++)
	sigma2[3 * i + j][3 * i + k] = etracks[i + 2][j + 2][k + 2];
  
  // perform the actual Z1 refit
  bool hasWidth = true;
  m_ResolutionModel->Set(higgs, m_WhichResolutionModel, ResolutionModel::UseZ1);
  ConstraintFit *MassFit = new ConstraintFit(Z_pdg_mass, hasWidth, Z_pdg_width);
  MassFit->SetResolutionModel(m_ResolutionModel);
  MassFit->MassFitInterface(parameters, sigma, 2);
  result.Z1_chiSq = MassFit->MassFitRun(parameters, sigma);
  
  // repeat for Z2
  if (quadrupletMass > 190.e3) {
    m_ResolutionModel->Set(higgs, m_WhichResolutionModel, ResolutionModel::UseZ2);
    MassFit->MassFitInterface(parameters2, sigma2, 2);
    result.Z2_chiSq = MassFit->MassFitRun(parameters2, sigma2);
    // if fit is not performed, Z2_chiSq is anyway reset to a default value when constructor of ConstraintFitResult is called
  }
  
  result.Z1_lepplus_4m_updated.SetXYZM(parameters[0][0],  parameters[0][1], parameters[0][2], parameters[0][3]);
  result.Z1_lepminus_4m_updated.SetXYZM(parameters[1][0],  parameters[1][1], parameters[1][2], parameters[1][3]);
  result.Z2_lepplus_4m_updated.SetXYZM(parameters2[0][0], parameters2[0][1], parameters2[0][2], parameters2[0][3]);
  result.Z2_lepminus_4m_updated.SetXYZM(parameters2[1][0], parameters2[1][1], parameters2[1][2], parameters2[1][3]);
  
  // std::cout << "constrained Z1 lepplus: (" << result.Z1_lepplus_4m_updated.Pt() << ", " << result.Z1_lepplus_4m_updated.Eta() << ", " << result.Z1_lepplus_4m_updated.Phi() << ", " << result.Z1_lepplus_4m_updated.M() << ")" << std::endl; // kostas
  // std::cout << "constrained Z1 lepminus: (" << result.Z1_lepminus_4m_updated.Pt() << ", " << result.Z1_lepminus_4m_updated.Eta() << ", " << result.Z1_lepminus_4m_updated.Phi() << ", " << result.Z1_lepminus_4m_updated.M() << ")" << std::endl; // kostas
  
  result.Z1_4m_updated = result.Z1_lepplus_4m_updated + result.Z1_lepminus_4m_updated;
  result.Z2_4m_updated = result.Z2_lepplus_4m_updated + result.Z2_lepminus_4m_updated;
  result.H_4m_updated  = result.Z1_4m_updated + result.Z2_4m_updated;
  
  delete MassFit;
  
  return result;
}

Float_t HiggsllqqAnalysis::getCandidateTriggerSF(Analysis::Quadrilepton * higgs)
{
  Float_t result(1);
  
  if (isMC()) {
    
    // find a fake run number representing the data period to which this MC event is somewhat associated
    Int_t representative_run_number = m_PileupReweighter->GetRandomRunNumber(ntuple->eventinfo.RunNumber());
    
    electron_quality el_quality;
    
    if (analysis_version() == "rel_17")
      el_quality = loosepp;
    else if (analysis_version() == "rel_17_2")
      el_quality = ML;
    
    
    // use MeV (set to true to use GeV)
    //m_MuonTrigSF->setThresholds(false, representative_run_number);
    
    std::vector<TLorentzVector> muons_4m;
    std::vector<TLorentzVector> electrons_4m;
    
    std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
    
    for (UInt_t i = 0; i < leptons.size(); i++) {
      if (leptons[i]->flavor() == Analysis::ChargedLepton::MUON)
	muons_4m.push_back(*(leptons[i]->Get4Momentum()));
      else
	electrons_4m.push_back(*(leptons[i]->Get4Momentum()));
    }
    
    result = m_MuonTrigSF->GetTriggerSF(representative_run_number, false, muons_4m, loose, electrons_4m, el_quality, 0).first;
  }
  
  return result;
}

Float_t HiggsllqqAnalysis::getSFWeight()
{
  // to be implemented
  
  return 1;
}

Bool_t HiggsllqqAnalysis::JetInHole()
{
  // this method has to be called AFTER having filled 
  // the GoodJets vector with the GetGoodObjects() method
  
  Float_t ptthr = 40000.;
  
  std::vector<Analysis::Jet *>::iterator jet_itr_i;
  for (jet_itr_i = m_GoodJets.begin(); jet_itr_i != m_GoodJets.end(); ++jet_itr_i) {
    ptthr = 40000.;
    
    D3PDReader::JetD3PDObjectElement *Jet = (*jet_itr_i)->GetJet();
    
    if(!isMC()) ptthr *= (1-Jet->BCH_CORR_JET())/(1-Jet->BCH_CORR_CELL());
    if((*jet_itr_i)->righteta()>-.1 && (*jet_itr_i)->righteta()<1.5 && (*jet_itr_i)->rightphi()>-.9 && (*jet_itr_i)->rightphi()<-.1)
      {
	if((*jet_itr_i)->rightpt()>ptthr) 
	  return kTRUE;
      }
  }
  return kFALSE;
}

Bool_t HiggsllqqAnalysis::JetKinematicFitterResult() {
  
  CandidatePair *m_bestdijet;
  std::vector<int > m_jetindex;
  
  m_bestdijet = new CandidatePair(-1,-1,-9999.,-9999.,-9999.);
  
  std::vector<TLorentzVector > jetvector;
  jetvector.clear();
  
  std::vector<Analysis::Jet *>::iterator jet_itr;
  
  for (jet_itr = m_GoodJets.begin(); jet_itr != m_GoodJets.end(); ++jet_itr) {
    
    TLorentzVector thisJet;
    thisJet.SetPtEtaPhiM((*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
    jetvector.push_back(thisJet);
  }
  
  m_jetkinematicfitter->clearParticles();
  m_bestdijet->Reset();
  
  m_jetkinematicfitter->addParticles(jetvector);
  *m_bestdijet = m_jetkinematicfitter->findBestPair();
  int idx1 = m_bestdijet->index1!=-1 ? m_bestdijet->index1 : 0;
  int idx2 = m_bestdijet->index2!=-1 ? m_bestdijet->index2 : 1;
  
  TLorentzVector j1;
  TLorentzVector j2;
  int ii =-1;
  
  for(jet_itr = m_GoodJets.begin(); jet_itr != m_GoodJets.end(); ++jet_itr) {
    ii++;
    if(ii==idx1)
      j1.SetPtEtaPhiM((*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
    if(ii==idx2)
      j2.SetPtEtaPhiM((*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
  }
  
  TLorentzVector hadZ = j1 + j2;
  
  if((m_bestdijet->index1!=-1)&&(m_bestdijet->index2!=-1)&&(hadZ.M()>60000.)&&(hadZ.M()<115000.))
    return kTRUE;
  else
    return kFALSE;
}

Bool_t HiggsllqqAnalysis::IsConsistentPt()
{
  if(getChannel() == HiggsllqqAnalysis::MU2 && m_GoodMuons.size()==2 && m_GoodElectrons.size()==0) {
    // a bit redundant, but it's better to be safe	    
    if(m_GoodMuons.at(0)->Get4Momentum()->Pt()>20000. || m_GoodMuons.at(1)->Get4Momentum()->Pt()>20000.||
       (m_GoodMuons.at(0)->Get4Momentum()->Pt()>12000. && m_GoodMuons.at(1)->Get4Momentum()->Pt()>12000.))
      return kTRUE;
  }
  else if(getChannel() == HiggsllqqAnalysis::E2 && m_GoodElectrons.size()==2 && m_GoodMuons.size()==0) {
    // a bit redundant, but it's better to be safe
    if(m_GoodElectrons.at(0)->Get4Momentum()->Pt()>20000. || m_GoodElectrons.at(1)->Get4Momentum()->Pt()>20000. ||
       (m_GoodElectrons.at(0)->Get4Momentum()->Pt()>14000. && m_GoodElectrons.at(1)->Get4Momentum()->Pt()>14000.))
      return kTRUE;
  }
  return kFALSE;
}


Bool_t HiggsllqqAnalysis::IsConsistentWithTrigger()
{
  if(getChannel() == HiggsllqqAnalysis::MU2 && m_GoodMuons.size()==2 && m_GoodElectrons.size()==0) {	    
    if((passesSingleMuonTrigger() && m_GoodMuons.at(0)->Get4Momentum()->Pt()>20000. && m_GoodMuons.at(1)->Get4Momentum()->Pt()>7000.) ||
       (passesDiMuonTrigger() && m_GoodMuons.at(0)->Get4Momentum()->Pt()>12000. && m_GoodMuons.at(1)->Get4Momentum()->Pt()>12000.))
      return kTRUE;
  }
  else if(getChannel() == HiggsllqqAnalysis::E2 && m_GoodElectrons.size()==2 && m_GoodMuons.size()==0) {
    if( (passesSingleElectronTrigger() && m_GoodElectrons.at(0)->Get4Momentum()->Pt()>20000. && m_GoodElectrons.at(1)->Get4Momentum()->Pt()>7000.) ||
	(passesDiElectronTrigger() && m_GoodElectrons.at(0)->Get4Momentum()->Pt()>14000. && m_GoodElectrons.at(1)->Get4Momentum()->Pt()>14000.) )
      return kTRUE;
  }
  return kFALSE;
}

Bool_t HiggsllqqAnalysis::NotMETclean()
{  
  Float_t ptthr  = 20000.;
  bool Bad_event = kFALSE;
  D3PDReader::JetD3PDObject *jet_branch(0);
  jet_branch = &(ntuple->jet_akt4topoem);
  D3PDReader::ElectronD3PDObject *el_branch(0);
  if (getElectronFamily() == Electron::GSF) el_branch = &(ntuple->el_GSF);
  else if (getElectronFamily() == Electron::noGSF) el_branch = &(ntuple->el);
  
  for (Int_t i = 0; i < jet_branch->n(); i++) {
    Analysis::Jet *jet = new Analysis::Jet(&((*jet_branch)[i]));
    D3PDReader::JetD3PDObjectElement *Jet = jet->GetJet();
    ApplyChangesJet(jet);
    if((Jet->isBadLooseMinus()!=0) && jet->rightpt()>ptthr)
      {
	Bad_event = kTRUE;
	prebadevent++;
	
	for (Int_t i = 0; i < el_branch->n(); i++) {
	  Analysis::ChargedLepton *lep = new Analysis::ChargedLepton(&((*el_branch)[i]), getElectronFamily());
	  //applyChanges(lep);
	  ApplyChangesElectron(lep);
	  if (isGood(lep)) {  
	    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
	    Float_t dR_GoodEl_BadJet = TMath::Sqrt(TMath::Power(jet->righteta() - el->cl_eta(), 2) + TMath::Power(TVector2::Phi_mpi_pi(jet->rightphi() - el->cl_phi()), 2));
	    if(dR_GoodEl_BadJet<0.4)
	      Bad_event = kFALSE;
	  }
	}
	if(Bad_event)
	  {
	    badevent++;
	    return kTRUE;
	  }
      }
  }
  return kFALSE;
}

Float_t HiggsllqqAnalysis::getCorrectMETValue()
{
  if (METtype_RefFinal) 
    return  ntuple->MET_RefFinal.et();
  else
    {
      cout<<"WARNING!!! Retuning MET_Final NOT Topo Implemented yet!!"<<endl;
      return  ntuple->MET_RefFinal.et();
    }
}

Bool_t HiggsllqqAnalysis::hasGoodMET()
{
  Bool_t pass_cut;
  if(GetDoLowMass()) pass_cut = (getCorrectMETValue() < 30000.);
  else pass_cut = (getCorrectMETValue() < 50000.);
  
  return pass_cut;
}

Float_t HiggsllqqAnalysis::GetMV1value(Analysis::Jet *jet)
{
  D3PDReader::JetD3PDObjectElement *Jet = (jet)->GetJet();
  
  Float_t w_IP3D            = Jet->flavor_weight_IP3D();
  Float_t w_SV1             = Jet->flavor_weight_SV1();
  Float_t w_JetFitterCOMBNN = Jet->flavor_weight_JetFitterCOMBNN();
  Float_t MV1               = mv1Eval(w_IP3D,w_SV1,w_JetFitterCOMBNN,jet->rightpt(),jet->righteta());
  
  return MV1;
}

Int_t HiggsllqqAnalysis::GetNumOfTags()
{
  Int_t tags = 0;
  
  std::vector<Analysis::Jet *>::iterator jet_itr_i;
  
  for (jet_itr_i = m_GoodJets.begin(); jet_itr_i != m_GoodJets.end(); ++jet_itr_i) {
    if(GetMV1value(*jet_itr_i) > 0.601713){
      tags++;
    }
  }
  return tags;
}

//Method to calculate the DiJet invariant mass for the tagged Jets!
Bool_t HiggsllqqAnalysis::JetDimassTagged() {
  
  TLorentzVector j1;
  TLorentzVector j2;
  bool first =kTRUE;
  
  std::vector<Analysis::Jet *>::iterator jet_itr;
  
  for (jet_itr = m_GoodJets.begin(); jet_itr != m_GoodJets.end(); ++jet_itr) {
    
    if((GetMV1value(*jet_itr) > 0.601713) && first)
      {
	first =kFALSE;
	j1.SetPtEtaPhiM(b_rescaling*(*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
      }
    if((GetMV1value(*jet_itr) > 0.601713) && !first)
      {
	j2.SetPtEtaPhiM(b_rescaling*(*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
      }
  }
  
  TLorentzVector hadZ = j1 + j2;
  
  if((hadZ.M() > 60000.) && (hadZ.M() <115000.))
    {
      return kTRUE;
    }
  else
    return kFALSE;
}
