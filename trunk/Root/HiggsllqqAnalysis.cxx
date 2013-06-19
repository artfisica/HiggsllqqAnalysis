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
    
    
    Information: 
    dolowmass for muons     = True   if GetDoLowMass() == True
    dolowmass for electrons = True   if GetDoLowMass() == True
    
    April 22th 2013

    Authors:
       Arturo Sanchez <arturos@cern.ch> or <arturos@ula.ve>
       Valerio Ippolito <valerio.ippolito@cern.ch>
*/

Bool_t cut_leptons   = kFALSE, 
  NotTrigConsistency = kFALSE, 
  NotPtConsistency   = kTRUE, 
  
// Smearing Options:
  MuonSmearing       = kTRUE, 
  JetSmearing        = kTRUE, 
  ElectronSmearing   = kTRUE,
  
// MET cleaning methods aplication  
  DoMETdataClean     = kTRUE, 
  
// Tag channels granularity
  TaggedChannel      = kFALSE,
  TagOneJetChannel   = kFALSE,
  
  is4lepGood         = kFALSE, 
  
// Type of Missing Et
  METtype_RefFinal   = kTRUE,
  
  DoLowMass          = kTRUE, 
  DoCaloMuons        = kTRUE,
  
// Do Jet Kinematic Fitter OR USE THE 2 LEADING JETS
  DoKinematicFitter  = kTRUE,
  FillGluon          = kTRUE,
  Look_b_SFs         = kTRUE;


//Global Jets Variables.
int Pair_jet1(-1), Pair_jet2(-1), Jone(800), Jtwo(900);

float corr_jet_pt1(-1.), corr_jet_pt2(-1.), ChiSq(-1.);

int Print_low_OR_high = 1; // 0 for LowSelection ; 1 for HighSelection

int     count_events(0),eventNow(-1),overElectron(0),overMuon(0),overJet(0); 
int     badevent(0), prebadevent(0), ptchange(0), ptelecChange(0);
int     periodBD(0),periodEH(0),periodI(0),periodJK(0),periodLM(0);
Float_t Muon0(0),Muon1(0),Muon2(0),Muon3(0),Muon4(0),Muon5(0),Muon6(0),Muon7(0),Muon8(0);
Float_t Electron0(0),Electron1(0),Electron2(0),Electron3(0),Electron4(0),Electron5(0),Electron6(0), b_rescaling = 1.00/*Fixed  in 1.05 (2011)*/; 


// Definition of the Leptonic (dilepton) invariant mass window:
Float_t Mll_low_min  = 20000.;
Float_t Mll_low_max  = 83000.;
Float_t Mll_high_min = 70000.;
Float_t Mll_high_max = 99000.;


// Definition of the Hadronic (dijet) invariant mass window:  70 (60) GeV < Mjj < 105 (115) GeV 
Float_t Mjj_low_min       = 20000.; // 60000.;
Float_t Mjj_low_max       = 200000.; // 115000.;
Float_t Mjj_high_min      = 20000.; // 70000.;
Float_t Mjj_high_max      = 200000.; // 105000.;


//Definition of the MET cut:
Float_t MET_low_cut  = 50000.; // 30000.;
Float_t MET_high_cut = 60000.; // 50000.;

int HFOR_value       = -999;


//MV1 operating point 70%
Float_t MV1_OP70 = 0.795; //0.601713; //  


//Actual Jet Cone Size in used
Float_t Cone_size= 0.4;


//Actual (2012) JVF CUT
Float_t JVF_CUT = 0.5; // 0.75

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
  // Int_t SYST_FLAG = 0; //SYST_FLAG is 0 for nominal scale, 1 or 2 for 1-sigma variations.
  
  if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {  
    
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    // first of all one must rescale energy in the crack! (both data and MC but only for 2011 so far)
    Float_t tmp_calibration(1.);
    
    if (analysis_version() == "rel_17") // rel. 17
      tmp_calibration = m_ElectronEnergyRescaler->applyMCCalibration(el->cl_eta(), el->cl_E() / TMath::CosH(el->tracketa()), egRescaler::EnergyRescalerUpgrade::Electron);
    
    Float_t tmp_E = el->cl_E() * TMath::Abs(tmp_calibration);
    Float_t tmp_Et = tmp_E / TMath::CosH(el->tracketa());
    Float_t tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    
    // then, apply the other corrections
    if (isMC() && doSmearing() && ElectronSmearing) {
      m_ElectronEnergyRescaler->SetRandomSeed(ntuple->eventinfo.EventNumber() + 100 * (Int_t)el->GetIndex());
      
      // false here means the MC is mc11c (no constant term)
      Float_t smearcorr = m_ElectronEnergyRescaler->getSmearingCorrection(el->cl_eta(), tmp_E, egRescaler::EnergyRescalerUpgrade::NOMINAL);
      
      tmp_E = tmp_E * smearcorr;
      tmp_Et = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    
    
    if (!isMC() && doSmearing() && ElectronSmearing) {
      tmp_E = m_ElectronEnergyRescaler->applyEnergyCorrection(el->cl_eta(), tmp_E, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::ZeeAllDown, 1.0, ntuple->eventinfo.RunNumber());
      tmp_Et = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    
    Float_t tmp_pt    = (TMath::Sqrt(TMath::Power(tmp_E, 2) - TMath::Power(lep->Get4Momentum()->M(), 2))) * TMath::Sin(2 * TMath::ATan(TMath::Exp(-el->tracketa())));
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
  
  if (isMC() && doSmearing() && JetSmearing) { //Jet Smearing Option! 14th June 2012
    double Eraw    = Jet->emscale_E();
    double eta_det = Jet->emscale_eta();
    double eta     = Jet->EtaOrigin();
    double phi     = Jet->PhiOrigin();
    double m       = Jet->MOrigin();
    
    double mu = (isMC() && ntuple->eventinfo.lbn()==1 && int(ntuple->eventinfo.averageIntPerXing()+0.5)==1) ? 0. : ntuple->eventinfo.averageIntPerXing();
    
    
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
      myJER->SmearJet_Syst(jet4v);// myJER->SmearJet_SystRel17(jet4v);
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
    m_GRL->SetXMLFile("./HiggsllqqAnalysis/grl/grl_2012.xml");
  } else if (getChannel() == HiggsllqqAnalysis::MUE) {
    m_GRL->SetXMLFile("./HiggsllqqAnalysis/grl/grl_2012.xml");
  } else if (getChannel() == HiggsllqqAnalysis::E2) {
    m_GRL->SetXMLFile("./HiggsllqqAnalysis/grl/grl_2012.xml");
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
  //if (!m_TrigDecisionToolD3PD->GetConfigSvc(kFALSE).SetConfigTree(fConfigTree)) {
  if (!m_TrigDecisionToolD3PD->SetConfigTree(fConfigTree)) {
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
  
  printAllOptions();
  
  // initiate the calibration tool
  TString jetAlgo="";
  if (getJetFamily() == 0)      jetAlgo="AntiKt4TopoEM";
  else if (getJetFamily() == 1) jetAlgo="AntiKt4LCTopo";
  
  
  
  TString JES_config_file;
  
  if (analysis_version() == "rel_17_2") {
    JES_config_file="ApplyJetCalibration/data/CalibrationConfigs/JES_Full2012dataset_Preliminary_Jan13.config"; //JES_Full2012dataset_Preliminary_AFII_Jan13.config"; 
  } 
  else if (analysis_version() == "rel_17") {
    JES_config_file="ApplyJetCalibration/data/CalibrationConfigs/Rel17_JES.config";
  }
  
  
  bool isData(0);
  if (isMC())
    isData = false;
  else if (!isMC())
    isData = true;
  
  myJES = new JetCalibrationTool(jetAlgo,JES_config_file, isData);
  myJER = new JetSmearingTool(jetAlgo);
  
  
  // initialize the kinematic fitter
  Info("doAnalysis", "Initializing JetKinematicFitter");
  int maxjet_KF = 7;
  m_jetkinematicfitter = new JetKinematicFitter(maxjet_KF,91187.6,2495.2); // mZ | gammaZ
  m_jetkinematicfitter->SetIsMC(isMC());
  
  
  // HFOR D3PD TOOL
  hforTool = new HforToolD3PD();
  
  
  // trigger decision tool
  m_TrigDecisionToolD3PD = new D3PD::TrigDecisionToolD3PD();
  
  
  // goodrunslists
  m_GRL = new Root::TGoodRunsListReader();
  // put 2011 first
  m_GRL->AddXMLFile("./HiggsllqqAnalysis/grl/grl_2011.xml");////*//// Change with the real GRL for llqq analysis!
  m_GRL->AddXMLFile("./HiggsllqqAnalysis/grl/grl_2011.xml");
  m_GRL->AddXMLFile("./HiggsllqqAnalysis/grl/grl_2011.xml");
  // then add 2012
  m_GRL->AddXMLFile("./HiggsllqqAnalysis/grl/grl_2012.xml");
  m_GRL->AddXMLFile("./HiggsllqqAnalysis/grl/grl_2012.xml");
  m_GRL->AddXMLFile("./HiggsllqqAnalysis/grl/grl_2012.xml");
  
  m_GRL->Interpret();
  
  
  // pileup
  Info("initialize_tools", "Initializing the pile-up reweighting tool...");
  
  m_PileupReweighter = new Root::TPileupReweighting("m_PileupReweighter");
  m_PileupReweighter->SetUnrepresentedDataAction(2);
  
  if (analysis_version() == "rel_17") { // mc11c, 2011
    m_PileupReweighter->AddConfigFile("./HiggsllqqAnalysis/packages/files/pileup/MC11c.prw.root");
    //m_PileupReweighter->AddLumiCalcFile("./HiggsllqqAnalysis/packages/files/pileup/ilumicalc_period_AllYear_Higgs_4l_2e2mu.root");
    m_PileupReweighter->AddLumiCalcFile("./Higgs4lepAnalysis/packages/files/pileup/ilumicalc_2011_AllYear_All_Good.root");
    m_PileupReweighter->SetDefaultChannel(109292);
  } else if (analysis_version() == "rel_17_2") { // mc12a, 2012
    m_PileupReweighter->AddConfigFile("./HiggsllqqAnalysis/packages/files/pileup/MC12a.prw.root");
    //m_PileupReweighter->AddLumiCalcFile("./HiggsllqqAnalysis/packages/files/pileup/ilumicalc_2012_AllYear_All_Good.root");
    m_PileupReweighter->AddLumiCalcFile("./Higgs4lepAnalysis/packages/files/pileup/ilumicalc_2012_period_AllYear_Higgs_4l_2e2mu.root");
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
  m_ElectronEnergyRescaler = new egRescaler::EnergyRescalerUpgrade();
  
  if (analysis_version() == "rel_17")
    m_ElectronEnergyRescaler->Init("egammaAnalysisUtils/share/EnergyRescalerData.root", "2011", "es2011a");
  else if (analysis_version() == "rel_17_2")
    m_ElectronEnergyRescaler->Init("egammaAnalysisUtils/share/EnergyRescalerData.root", "2012", "es2012");
  
  
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
  
  
  m_MuonEffSF     = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_directory, mu_alg, unit, muonEffConfig);
  m_MuonEffSFCalo = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_directory, mu_alg_calo, unit, muonEffConfig);
  m_MuonEffSFSA   = new Analysis::AnalysisMuonConfigurableScaleFactors(muon_sf_directory, mu_alg_sa, unit, muonEffConfigSA);
  
  
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
    m_MuonSmearer = new MuonSmear::SmearingClass("Data12", mu_alg, "pT", "Rel17.2", muon_sm_directory);
  m_MuonSmearer->UseScale(1);
  m_MuonSmearer->UseImprovedCombine();
  
  
//   // MuonTrigSF
//   Info("initialize_tools", "Initializing the LeptonTriggerSF tool...");
//   if (analysis_version() == "rel_17") { // 2011
//     m_MuonTrigSF = new LeptonTriggerSF(2011, "TrigMuonEfficiency/share/", "muon_trigger_sf_mc11c.root");
//   } else if (analysis_version() == "rel_17_2") { // 2012
//     m_MuonTrigSF = new LeptonTriggerSF(2012, "TrigMuonEfficiency/share", "muon_trigger_sf_2012.root");
//   }
  
  // MuonTrigSF
  Info("initialize_tools", "Initializing the LeptonTriggerSF tool...");
  if (analysis_version() == "rel_17") { // 2011
    m_MuonTrigSF = new LeptonTriggerSF(2011, "TrigMuonEfficiency/share/", "muon_trigger_sf_mc11c.root", "ElectronEfficiencyCorrection/data/", "rel17p0.v01");
  } else if (analysis_version() == "rel_17_2") { // 2012
    m_MuonTrigSF = new LeptonTriggerSF(2012, "TrigMuonEfficiency/share", "muon_trigger_sf_2012.root", "ElectronEfficiencyCorrection/data/", "rel17p2.v02");
  }
  
  
  
  // Muon and Electron trigger match
  m_MuonTriggerMatchTool     = new MuonTriggerMatching(&m_trigNavVar);
  m_ElectronTriggerMatchTool = new ElectronTriggerMatching(&m_trigNavVar);
  
  
  // delay ggF tool initialization
  m_ggFReweighter = 0;
  
  
  // delay jet calibration tool initialization
  //m_jetCalibrationTool = 0;    // NOT! please, see up!!
  
  
  if(Look_b_SFs) //initiate the bTag SFs tool
    {
      calib = new Analysis::CalibrationDataInterfaceROOT("MV1","BTagCalibration.env","./HiggsllqqAnalysis/util/btagSF/"); 
      ajet.jetAuthor = jetAlgo;
      uncertainty = Analysis::Total;
    } 
  
  return kTRUE;
}


Float_t HiggsllqqAnalysis::getD0SmearSigma(Int_t index_of_lepton, Int_t nBL, Float_t pt, Float_t eta)
{
  
  // from https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/IPSmearing
  m_smearD0_rand.SetSeed(ntuple->eventinfo.EventNumber() + 100 * index_of_lepton);
  
  if (nBL >= 2) nBL = 2;
  
  float sinTheta = 1. / TMath::CosH(eta);
  //float p = pt * TMath::CosH(eta);
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
      
      //Printing the Cutflow for signal mass value!
      if(ntuple->eventinfo.mc_channel_number() < 160420)
	Print_low_OR_high=0;
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
  //delete m_jetCalibrationTool;
  
  return kTRUE;
}


Bool_t HiggsllqqAnalysis::initialize_analysis()
{
  
  //Mean of Jets
  Mean_jets =0.0;
  good_events =0.0;
  
  // open the output file
  m_outputFile = new TFile(m_outputFileName, "RECREATE");
  
  
  // prepare output histo
  m_generatedEntriesHisto = new TH1D("generatedEntriesHisto", "number of processed entries", 10, 0, 10);
  m_generatedEntriesHisto->Fill("raw", 0);
  m_generatedEntriesHisto->Fill("Sherpa_Veto", 0);  
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix", 0);
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix_EvWeight", 0);
  m_generatedEntriesHisto->Fill("with_ggF_as_well", 0);
  m_generatedEntriesHisto->Fill("with_pu_and_vertex", 0);
  
  
  //Definition of the different cutflow histogram: Event, and per channel.
  h_cutflow = new TH1D("cutflow_E2","",20,0,20);           //("cutflow","",20,0,20);
  h_cutflow->Fill("Entries",0);
  h_cutflow->Fill("HFOR",0);
  h_cutflow->Fill("GRL or SherpaVeto",0);
  h_cutflow->Fill("larError",0);
  h_cutflow->Fill("Trigger",0);
  h_cutflow->Fill("Vertex",0);
  h_cutflow->Fill("METcleaning",0);
  h_cutflow->Fill("LArHole",0);
  h_cutflow->Fill("NumberOfLeptons",0);
  h_cutflow->Fill("OppositeSign",0);
  h_cutflow->Fill("PtLeptons",0);
  h_cutflow->Fill("TriggerConsistency",0);
  h_cutflow->Fill("MET",0);
  h_cutflow->Fill("TwoJets",0);
  h_cutflow->Fill("NumTagJets",0);
  h_cutflow->Fill("DileptonMass",0);
  h_cutflow->Fill("DiJetMass",0);
  
  h_cutflow_weight = new TH1D("cutflow_MU2","",20,0,20); //("cutflow_weight","",20,0,20);
  h_cutflow_weight->Fill("Entries",0);
  h_cutflow_weight->Fill("HFOR",0);
  h_cutflow_weight->Fill("GRL or SherpaVeto",0);
  h_cutflow_weight->Fill("larError",0);
  h_cutflow_weight->Fill("Trigger",0);
  h_cutflow_weight->Fill("Vertex",0);
  h_cutflow_weight->Fill("METcleaning",0);
  h_cutflow_weight->Fill("LArHole",0);
  h_cutflow_weight->Fill("NumberOfLeptons",0);
  h_cutflow_weight->Fill("OppositeSign",0);
  h_cutflow_weight->Fill("PtLeptons",0);
  h_cutflow_weight->Fill("TriggerConsistency",0);
  h_cutflow_weight->Fill("MET",0);
  h_cutflow_weight->Fill("TwoJets",0);
  h_cutflow_weight->Fill("NumTagJets",0);
  h_cutflow_weight->Fill("DileptonMass",0);
  h_cutflow_weight->Fill("DiJetMass",0);
  
  
  //Initialization of the qqll ntuple
  InitReducedNtuple();
  
  
  // Set the TestSelection Tree
  analysistree = new TTree("tree","complete analysis output tree");
  SetAnalysisOutputBranches(&m_outevent);
  
  
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
  m_ElectronCutflow.clear();
  m_MuonCutflow.clear();
  m_JetCutflow.clear();
  
  
  
  for (UInt_t i = 0; i < m_Channels.size(); i++) {
    
    m_EventCutflow.push_back(Analysis::CutFlowTool("Plain_" + chan_name[i]));
    m_EventCutflow[i].addCut("HFOR");
    m_EventCutflow[i].addCut("GRL or SherpaVeto");
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
    
    m_EventCutflow_rw.push_back(Analysis::CutFlowTool("Reweighted_" + chan_name[i]));
    m_EventCutflow_rw[i].addCut("HFOR");
    m_EventCutflow_rw[i].addCut("GRL or SherpaVeto");
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
    
    m_ElectronCutflow.push_back(Analysis::CutFlowTool("Electrons_" + chan_name[i]));
    m_ElectronCutflow[i].addCut("family");
    m_ElectronCutflow[i].addCut("algorithm");
    m_ElectronCutflow[i].addCut("quality");
    m_ElectronCutflow[i].addCut("eta");
    m_ElectronCutflow[i].addCut("Et");
    m_ElectronCutflow[i].addCut("objectquality");
    m_ElectronCutflow[i].addCut("z0");
    m_ElectronCutflow[i].addCut("d0Sig");
    m_ElectronCutflow[i].addCut("Isolation");
    m_ElectronCutflow[i].addCut("overlap");
    
    m_MuonCutflow.push_back(Analysis::CutFlowTool("Muons_" + chan_name[i]));
    m_MuonCutflow[i].addCut("family");
    m_MuonCutflow[i].addCut("quality");
    m_MuonCutflow[i].addCut("cosmic");
    m_MuonCutflow[i].addCut("eta");
    m_MuonCutflow[i].addCut("pt");
    m_MuonCutflow[i].addCut("MCP");
    m_MuonCutflow[i].addCut("z0");
    m_MuonCutflow[i].addCut("d0Sig"); 
    m_MuonCutflow[i].addCut("Isolation");
    m_MuonCutflow[i].addCut("overlap");
    
    m_JetCutflow.push_back(Analysis::CutFlowTool("Jets_" + chan_name[i]));
    m_JetCutflow[i].addCut("jetCleaning");
    m_JetCutflow[i].addCut("kinematics");
    m_JetCutflow[i].addCut("Pileup");
    m_JetCutflow[i].addCut("overlap");
  }    
  
  
  
  //   if (analysis_version() == "rel_17_2") {
  //     if (isMC())
  //       m_MuonResolutionAndMomentumScaleFactors = new Analysis::MuonResolutionAndMomentumScaleFactors("./MuonMomentumCorrections/share/final_scale_factors_MC12_smearing.txt");
  //     else
  //       m_MuonResolutionAndMomentumScaleFactors = new Analysis::MuonResolutionAndMomentumScaleFactors("./MuonMomentumCorrections/share/final_scale_factors_data2012.txt");
  //   } else if (analysis_version() == "rel_17") {
  //     if (isMC())
  //       m_MuonResolutionAndMomentumScaleFactors = new Analysis::MuonResolutionAndMomentumScaleFactors("./MuonMomentumCorrections/share/final_scale_factors_MC11_smeared.txt");
  //     else
  //       m_MuonResolutionAndMomentumScaleFactors = new Analysis::MuonResolutionAndMomentumScaleFactors("./MuonMomentumCorrections/share/final_scale_factors_data2011.txt");
  //   }
  
  
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

Bool_t HiggsllqqAnalysis::SherpaPt0Veto()
{
  
  // check if an event in Sherpa boosted samples Z+jets pt0 should be remove (PT2>70GeV)
  Bool_t result(kFALSE);
  
  
  if (isMC()) {
    
    if (ntuple->eventinfo.mc_channel_number() > 167748 && ntuple->eventinfo.mc_channel_number() < 167758) {
      
      int pdg(0),one(0);
      TLorentzVector lepton_1;
      TLorentzVector lepton_2;
      TLorentzVector Z_12;
      
      for (Int_t i = 0; i < ntuple->mc.n(); i++) {
	
	if (ntuple->mc[i].status() != 3) { continue; }
	
	pdg = TMath::Abs(ntuple->mc[i].pdgId());
	
	if(pdg < 11) { continue; }
	if(pdg > 16) { continue; }
	
	one++;
	
	if(one==1)
	  lepton_1.SetPtEtaPhiM(ntuple->mc[i].pt(),
				ntuple->mc[i].eta(),
				ntuple->mc[i].phi(),
				ntuple->mc[i].m());
	
	if(one==2)
	  lepton_2.SetPtEtaPhiM(ntuple->mc[i].pt(),
				ntuple->mc[i].eta(),
				ntuple->mc[i].phi(),
				ntuple->mc[i].m());
	
      } // loop over truth particles
      
      Z_12 = lepton_1 + lepton_2;
      
      if(Z_12.Pt()>70000. && one==2) {
	//cout<<"Sherpa Event rejected. PT2 = "<<Z_12.Pt()<<endl;
	result = kTRUE;
      }
    } // Sherpa Z+Jets samples
  }
  
  return result;
}


Int_t HiggsllqqAnalysis::getLastCutPassed()
{
  
  Int_t last(-1);
  
  
  //Sherpa Pt0 samples event veto (PT2>70GeV)
  //if(SherpaPt0Veto()) return last;
  
  
  // Including hfor veto 
  if (isMC())
    HFOR_value = hforTool->getDecision(ntuple->eventinfo.mc_channel_number(),
				       mc_n,
				       mc_pt,
				       mc_eta,
				       mc_phi,
				       mc_m,
				       mc_pdgId,
				       mc_status,
				       mc_vx_barcode,
				       mc_parent_index,
				       mc_child_index,
				       HforToolD3PD::BBONLY);
  else
    HFOR_value = -99999;
  
  
  if (HFOR_value!=4) last = HllqqCutFlow::HFOR;
  else return last;   
  
  
  
  
  if (isMC()) {
    if (!SherpaPt0Veto()) last = HllqqCutFlow::GRL;
    else return last;
  }
  else {
    if (passesGRL()) last = HllqqCutFlow::GRL;
    else return last;
  }
  

  
  if (ntuple->eventinfo.larError() != 2 &&
      (isMC() ||  ntuple->eventinfo.tileError() != 2) &&
      (isMC() || (ntuple->eventinfo.coreFlags() & 0x40000) == 0)
      ) last = HllqqCutFlow::larError;
  else return last;	  
  
  
  
  if (passesTrigger()) last = HllqqCutFlow::Trigger;
  else return last;
  
  
  if (hasGoodVertex()) last = HllqqCutFlow::Vertex;
  else return last;
  
  
  //METHOD THAT GET ALL THE LEPTONS: First call to get the Jets and study the Cleaning of the Event  
  getGoodLeptons();
  
  
  if(NotMETclean() /*&& !isMC()*/ && DoMETdataClean) return last;
  else last = HllqqCutFlow::METcleaning;
  
  
  if(analysis_version() == "rel_17" && ((JetInHole() && isMC() && getPeriod()==DataPeriod::y2011_EH) || (JetInHole() && !isMC()))) return last;
  else last = HllqqCutFlow::LArHole;
  
  
  //METHOD THAT GET ALL THE LEPTONS: Second Call the really get all the Good Object that will be past the different analysis requests!!
  getGoodLeptons();
  
  
  //Number of Leptons Cut  
  if ((getChannel() == HiggsllqqAnalysis::MU2 && m_GoodMuons.size()    == 2 && m_GoodElectrons.size() == 0 && (GetDoLowMass() || Pair_Quality())) || /* Staco CB+ST muons pT>25 GeV, |eta|<2.5 */
      (getChannel() == HiggsllqqAnalysis::E2 && m_GoodElectrons.size() == 2 && m_GoodMuons.size()     == 0 && (GetDoLowMass() || Pair_Quality())) || /* el->mediumPP()==1, ET>25 GeV and |eta|<2.47 */
      (getChannel() == HiggsllqqAnalysis::MUE && m_GoodMuons.size()    == 1 && m_GoodElectrons.size() == 1 && !GetDoQCDSelection())) last = HllqqCutFlow::NumberOfLeptons;  
  else return last;
  
  
  //OS selection on the 2 analysis channels
  int chargeprod = 0;
  
  if(getChannel()      == HiggsllqqAnalysis::MU2 && !GetDoQCDSelection())
    chargeprod = (m_GoodMuons.at(1))->charge()     * (m_GoodMuons.at(0))->charge(); 
  
  else if(getChannel() == HiggsllqqAnalysis::E2  && !GetDoQCDSelection() && GetDoLowMass())
    chargeprod = (m_GoodElectrons.at(1))->charge() * (m_GoodElectrons.at(0))->charge(); 
  
  if (chargeprod==-1 || (getChannel() != HiggsllqqAnalysis::MU2 && !GetDoLowMass()) || GetDoQCDSelection()) last = HllqqCutFlow::OppositeSign;
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
  
  
  //Minimun number of tagged or untagged jets
  if((GetNumOfTags()==2 && TaggedChannel && !TagOneJetChannel) || 
     (GetNumOfTags()==1 && !TaggedChannel && TagOneJetChannel) || 
     (GetNumOfTags()==0 && !TaggedChannel && !TagOneJetChannel)) last = HllqqCutFlow::NumTagJets;
  else return last;
  
  
  //Dielpton Mass windows
  getDileptons(); //Create a single object using the pair of leptons  
  std::vector<Analysis::Dilepton *>::iterator Z_itr_i= m_Dileptons.begin();
  Double_t DilepMass = (*Z_itr_i)->Get4Momentum()->M();
  
  if((GetDoLowMass()  && DilepMass>Mll_low_min  && DilepMass<Mll_low_max )  ||
     (!GetDoLowMass() && DilepMass>Mll_high_min && DilepMass<Mll_high_max)) last = HllqqCutFlow::DileptonMass;
  else return last;
  
  
  //Invariant mass of the dijet
  if((JetKinematicFitterResult() && !TaggedChannel && !TagOneJetChannel) ||
     (JetDimassTagged()          &&  TaggedChannel && !TagOneJetChannel) ||
     (JetDimassOneTagged()       && !TaggedChannel &&  TagOneJetChannel))
    last = HllqqCutFlow::DiJetMass;
  
  else return last;
  
  
  //GET THE FINAL OBJECTS
  getGoodObjects();
  
  return last;
}


Bool_t HiggsllqqAnalysis::passesGRL()
{
  if (isMC()) {
    return kTRUE;
  } else {
    if (ntuple->eventinfo.RunNumber() >= 215643) return kTRUE;
    
    return (m_GRL->GetGoodRunsList(getChannel() + HllqqGRL::data11).HasRunLumiBlock(ntuple->eventinfo.RunNumber(), ntuple->eventinfo.lbn()) ||
	    m_GRL->GetGoodRunsList(getChannel() + HllqqGRL::data12).HasRunLumiBlock(ntuple->eventinfo.RunNumber(), ntuple->eventinfo.lbn()));
  }
}


Bool_t HiggsllqqAnalysis::hasGoodVertex()
{
  //for (Int_t i = 0; i < ntuple->vxp.n(); i++) {
  if (ntuple->vxp[0].trk_n() >= 3) return kTRUE;
  //}
  else return kFALSE;
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
    } else if (period <= DataPeriod::y2012_E) {
      // 2012: A-E
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
      if(GetDoLowMass())
	chain_name = "EF_2mu13;EF_mu18_tight_mu8_EFFS";
      else
	chain_name = "EF_2mu13";
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
      if(GetDoLowMass())
	chain_name = "EF_2mu13;EF_mu18_tight_mu8_EFFS";
      else
	chain_name = "EF_2mu13";
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
    } else if (period == DataPeriod::y2012_E) {
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
    else if (period <= DataPeriod::y2012_E) {
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
      chain_name = "EF_e10_medium_mu6";
    }
    // 2012: A-E
    else if (period <= DataPeriod::y2012_E) {
      chain_name = "EF_e12Tvh_medium1_mu8;EF_e24vhi_loose1_mu8";
    }
  } else {
    // MC
    if (period == DataPeriod::y2011_BD) {
      chain_name = "EF_e10_medium_mu6";
    } else if (period == DataPeriod::y2011_EH) {
      chain_name = "EF_e10_medium_mu6";
    } else if (period == DataPeriod::y2011_I) {
      chain_name = "EF_e10_medium_mu6";
    } else if (period == DataPeriod::y2011_J) {
      chain_name = "EF_e10_medium_mu6";
    } else if (period == DataPeriod::y2011_K) {
      chain_name = "EF_e10_medium_mu6";
    } else if (period == DataPeriod::y2011_LM) {
      chain_name = "EF_e10_medium_mu6";
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
  } 
  
  else if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    // Int_t SYST_FLAG = 0; //SYST_FLAG is 0 for nominal scale, 1 or 2 for 1-sigma variations.
    
    // first of all one must rescale energy in the crack! (both data and MC but only for 2011 so far)
    Float_t tmp_calibration(1.);
    if (analysis_version() == "rel_17") // rel. 17
      tmp_calibration = m_ElectronEnergyRescaler->applyMCCalibration(el->cl_eta(), el->cl_E() / TMath::CosH(el->tracketa()), egRescaler::EnergyRescalerUpgrade::Electron);
    Float_t tmp_E = el->cl_E() * TMath::Abs(tmp_calibration);
    Float_t tmp_Et = tmp_E / TMath::CosH(el->tracketa());
    Float_t tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    
    // then, apply the other corrections
    if (isMC() && doSmearing() && ElectronSmearing) {
      m_ElectronEnergyRescaler->SetRandomSeed(ntuple->eventinfo.EventNumber() + 100 * (Int_t)el->GetIndex());
      
      // false here means the MC is mc11c (no constant term)
      Float_t smearcorr = m_ElectronEnergyRescaler->getSmearingCorrection(el->cl_eta(), tmp_E, egRescaler::EnergyRescalerUpgrade::NOMINAL);
      
      tmp_E     = tmp_E * smearcorr;
      tmp_Et    = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    
    if (!isMC()) {
      tmp_E     = m_ElectronEnergyRescaler->applyEnergyCorrection(el->cl_eta(), tmp_E, egRescaler::EnergyRescalerUpgrade::Electron, egRescaler::EnergyRescalerUpgrade::Nominal, 1.0, ntuple->eventinfo.RunNumber());
      tmp_Et    = tmp_E / TMath::CosH(el->tracketa());
      tmp_Et_cl = tmp_E / TMath::CosH(el->cl_eta());
    }
    
    Float_t tmp_pt    = (TMath::Sqrt(TMath::Power(tmp_E, 2) - TMath::Power(lep->Get4Momentum()->M(), 2))) * TMath::Sin(2 * TMath::ATan(TMath::Exp(-el->tracketa())));
    Float_t tmp_pt_cl = (TMath::Sqrt(TMath::Power(tmp_E, 2) - TMath::Power(lep->Get4Momentum()->M(), 2))) * TMath::Sin(2 * TMath::ATan(TMath::Exp(-el->cl_eta())));
    
    // correct lepton 4-momenta
    lep->Get4Momentum()->SetPtEtaPhiM(tmp_pt, el->tracketa(), el->trackphi(), lep->Get4Momentum()->M());
    lep->Get4Momentum_ID()->SetPtEtaPhiM(el->trackpt(), el->tracketa(), el->trackphi(), lep->Get4Momentum()->M()); // untouched by energy corrections!
    lep->Get4Momentum_SA()->SetPtEtaPhiM(tmp_pt_cl, el->cl_eta(), el->cl_phi(), lep->Get4Momentum()->M());
  } 
  
  //If is not an Electron neither a Muon!!
  else {
    Abort("Unknown lepton flavor");
  }
}


void HiggsllqqAnalysis::applyChanges(Analysis::Jet *jet)
{
  Float_t tmp_E(-9999.9), tmp_pt(-9999.9), tmp_eta(-9999.9), tmp_phi(-9999.9), tmp_Et(-9999.9);
  Bool_t sysstudy = GetSysStudy(); 
  
  // For the pile-up correction, we need mu and NPV(2+ tracks)    
  double mu = (isMC() && ntuple->eventinfo.lbn()==1 && int(ntuple->eventinfo.averageIntPerXing()+0.5)==1) ? 0. : ntuple->eventinfo.averageIntPerXing();
  int NPV=0;
  
  D3PDReader::JetD3PDObjectElement *this_jet = jet->GetJet();
  
  if (/*isMC() &&*/ doSmearing() && JetSmearing) {
    
    //inside the event loop
    
    if (analysis_version() == "rel_17_2") {      
      
      double Eraw    = this_jet->emscale_E();
      double eta     = this_jet->emscale_eta(); //EtaOrigin();
      double phi     = this_jet->emscale_phi();
      double m       = this_jet->emscale_m();   //MOrigin();
      double Ax      = this_jet->ActiveAreaPx();
      double Ay      = this_jet->ActiveAreaPy();
      double Az      = this_jet->ActiveAreaPz();
      double Ae      = this_jet->ActiveAreaE();
      double rho     = ntuple->Eventshape.rhoKt4EM();
      
      for (Int_t i = 0; i < ntuple->vxp.n(); i++) {
	if (ntuple->vxp[i].trk_n() >= 2) NPV++;
      }
      
      //arturo
      //if(ntuple->eventinfo.EventNumber()==5553)      
      //    cout<<"E = "<<Eraw<<" eta= "<<eta<<" phi= "<<phi<<" m= "<<m<<" Ax= "<<Ax<<" Ay= "<<Ay<<" Az= "<<Az<<" Ae= "<<Ae<<" rho= "<<rho<<" mu= "<<mu<<" NPV= "<<NPV<<endl;
      
      // Calibrate the jet!
      // Pile-up, origin, EtaJES correction applied, i.e. to OFFSET_ORIGIN_ETAJES scale
      TLorentzVector jet4v = myJES->ApplyJetAreaOffsetEtaJES(Eraw,eta,phi,m,Ax,Ay,Az,Ae,rho,mu,NPV);
      
      
      // The below is systematic evaluation, and ONLY for MC
      // Smear the jet to match the MC resolution+1 sigma!
      if(isMC() && sysstudy){
	myJER->SetSeed(ntuple->eventinfo.EventNumber());
	myJER->SmearJet_Syst(jet4v);
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
    
    else if (analysis_version() == "rel_17") {
      
      double Eraw    = this_jet->emscale_E();
      double eta     = this_jet->EtaOrigin();
      double eta_det = this_jet->emscale_eta();
      double phi     = this_jet->emscale_phi();
      double m       = this_jet->MOrigin();
      
      // Calibrate the jet!
      // Pile-up, origin, EtaJES correction applied, i.e. to OFFSET_ORIGIN_ETAJES scale
      TLorentzVector jet4v = myJES->ApplyOffsetEtaJES(Eraw,eta_det,eta,phi,m,mu,NPV);  
      
      // The below is systematic evaluation, and ONLY for MC
      // Smear the jet to match the MC resolution+1 sigma!
      if(isMC() && sysstudy){
	myJER->SetSeed(ntuple->eventinfo.EventNumber());
	myJER->SmearJet_Syst(jet4v);
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
}


void HiggsllqqAnalysis::getMuons(D3PDReader::MuonD3PDObject *mu_branch, Int_t family)
{
  // fills m_Muons with those muons passing the one-lepton muon selection
  // (no overlap removal at this stage)
  
  for (Int_t i = 0; i < mu_branch->n(); i++) {
    Analysis::ChargedLepton *lep = new Analysis::ChargedLepton(&((*mu_branch)[i]), family);
    applyChanges(lep);
    //ApplyChangesMuon(lep);
    
    if(isGood(lep)) m_Muons.push_back(lep);
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
    applyChanges(lep);
    //ApplyChangesElectron(lep);
    
    if (isGood(lep)) m_Electrons.push_back(lep);
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
    applyChanges(jet);
    //ApplyChangesJet(jet);
    
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
  
  // utility vector, tells if i-th muons must be skipped if overlapping with a jet
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
	
	if ((mu_i->Get4Momentum()->DeltaR(*(jet->Get4Momentum()))<0.3 &&  GetDoLowMass()) ||
	    (mu_i->Get4Momentum()->DeltaR(*(jet->Get4Momentum()))<0.0 && !GetDoLowMass())) {
	  // found an jet overlapped to a jet
	  skip_muon[i] = kTRUE;
	} // overlapping muon/jet
      } // jet loop
    }
    else skip_muon[i] = kTRUE;  
  } // Jet loop
  
  
  
  // fill the final tree with those Muons which can be used in the analysis
  i = -1;
  for (muon_itr = m_Muons.begin(); muon_itr != m_Muons.end(); ++muon_itr) {
    i++;
    if (skip_muon[i]) continue;
    
    (*muon_itr)->set_lastcut(HllqqMuonQuality::overlap);
    m_GoodMuons.push_back((*muon_itr));
  } // Muon loop
}


void HiggsllqqAnalysis::getGoodElectrons()
{
  // fills m_GoodElectrons with those electrons passing the full electron selection
  // (quality and overlap removal)
  
  // utility vector, tells if i-th electron must be skipped since overlapping with electron or muon
  std::vector<Bool_t> skip_electron;
  
  // remove overlap among electrons
  skip_electron.clear();
 
  // utility vector, tells if i-th electron must be skipped since overlapping with electron or muon
  // same for muons, since CALO muons overlapping with electrons must be removed
  std::vector<Bool_t> skip_muon;
  
  // remove overlap among electrons
  skip_muon.clear();
 

 
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
	if (el_j->GetElectron()->trackd0()     == el_i->GetElectron()->trackd0()    &&
	    el_j->GetElectron()->trackz0()     == el_i->GetElectron()->trackz0()    &&
	    el_j->GetElectron()->tracktheta()  == el_i->GetElectron()->tracktheta() &&
	    el_j->GetElectron()->trackphi()    == el_i->GetElectron()->trackphi()   &&
	    el_j->GetElectron()->trackqoverp() == el_i->GetElectron()->trackqoverp()) {
	  // found two overlapped electrons, skip the lowest-Et one
	  Int_t to_be_skipped = (el_i->Get4Momentum()->Et() > el_j->Get4Momentum()->Et()) ? j : i;
	  skip_electron[to_be_skipped] = kTRUE;
	}
      } // overlapping electrons rel. 17
      else if (analysis_version() == "rel_17_2") { // rel. 17.2
	if (el_j->GetElectron()->Unrefittedtrack_d0()     == el_i->GetElectron()->Unrefittedtrack_d0()    &&
	    el_j->GetElectron()->Unrefittedtrack_z0()     == el_i->GetElectron()->Unrefittedtrack_z0()    &&
	    el_j->GetElectron()->Unrefittedtrack_theta()  == el_i->GetElectron()->Unrefittedtrack_theta() &&
	    el_j->GetElectron()->Unrefittedtrack_phi()    == el_i->GetElectron()->Unrefittedtrack_phi()   &&
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
	if (el->Get4Momentum_ID()->DeltaR(*(mu->Get4Momentum_ID())) < 0.2) {
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
  
  
//   // remove CALO muons overlapping with electrons
//   std::vector<Analysis::ChargedLepton *> tmpGoodMuons;
  
//   for (UInt_t j = 0; j < skip_muon.size(); j++) {
//     if (skip_muon[j]) {
//       // do not add it to the final good muon list, and retrocede its last cut
//       m_GoodMuons[j]->set_lastcut(HllqqMuonQuality::overlap - 1);
//     } else {
//       tmpGoodMuons.push_back(m_GoodMuons[j]);
//     }
//   } // good mu loop
  
//   // replace good muon vector with overlap-removed one
//   m_GoodMuons = tmpGoodMuons;
  
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
  
  
  //Sort the jets by pt
  if(m_GoodJets.size()>0)
    SortIndex(m_GoodJets);
}


void HiggsllqqAnalysis::SortIndex(std::vector<Analysis::Jet *> &vec)
{
  bool change_made = true;
  while (change_made) {
    change_made = false;
    
    for (unsigned int i = 0; i < vec.size() - 1; i++) {
      Analysis::Jet *jeta = vec.at(i);
      Analysis::Jet *jetb = vec.at(i+1);
      
      
      if (jeta->rightpt() < jetb->rightpt()) {
	Analysis::Jet *jettmp = vec.at(i);
        vec.at(i) = vec.at(i + 1);
        vec.at(i + 1) = jettmp;
        change_made = true;
      }
    }
  }
}


void HiggsllqqAnalysis::getDileptons()
{
  
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
  if (getMuonFamily() == Muon::STACO)     mu_branch = &(ntuple->mu_staco);
  else if (getMuonFamily() == Muon::MUID) mu_branch = &(ntuple->mu_muid);
  else Abort("Unknown muon family requested");
  
  getMuons(mu_branch, getMuonFamily());
  if (DoCaloMuons)
    getMuons(&(ntuple->mu_calo), Muon::CALO);
  
  
  // add electrons
  D3PDReader::ElectronD3PDObject *el_branch(0);
  if (getElectronFamily() == Electron::GSF)        el_branch = &(ntuple->el_GSF);
  else if (getElectronFamily() == Electron::noGSF) el_branch = &(ntuple->el);
  else Abort("Unknown electron family requested");
  
  getElectrons(el_branch, getElectronFamily());
  
  
  // add Jets     
  /*Include Jet Families??? error
    a)jet_AntiKt4LCTopo
    b)jet_akt4topoem
    c) jet FAT!?
  */
  D3PDReader::JetD3PDObject *jet_branch(0);
  if (getJetFamily() == 0)      jet_branch = &(ntuple->jet_akt4topoem);
  else if (getJetFamily() == 1) jet_branch = &(ntuple->jet_AntiKt4LCTopo);
  else Abort("Unknown jet family requested");
  
  getJets(jet_branch);
  
  
  ////////
  // Part 1: remove overlap among electrons and between muons and jets
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
  // must be called once per event
  if (m_called_getGoodObjects) return;
  else m_called_getGoodObjects = kTRUE;
  
  ////////
  // Build dileptons
  ////////
  
  getDileptons();
  
  ////////
  // now we have filled
  //  - m_Muons, m_Electrons with pointers to muons and electrons passing selection cuts
  //  - m_GoodMuons, m_GoodElectrons with pointers to muons and electrons passing selection cuts and overlap removal
  //  - m_Dileptons with pointers to dileptons passing selection criteria (neutral charge in default analysis)
  // to be sure we are dealing correctly with object deletion, ALL and ONLY the objects in
  //      m_Muons      m_Electrons    m_Dileptons
  // must be deleted at the end of execute_analysis
  ////////
}


Bool_t HiggsllqqAnalysis::GetGoodObjects()
{
  /*
    See the HiggsqqllAnalysis Code for more information!!! Not implemented here, see Method above!    
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
      Bool_t scthits = (mu->nSCTHits() + mu->nSCTDeadSensors() > 5);
      if (!scthits) return kFALSE;
    } else if (analysis_version() == "rel_17_2") {
      Bool_t pixhits = (mu->nPixHits() + mu->nPixelDeadSensors() > 0);
      if (!pixhits) return kFALSE;
      Bool_t scthits = (mu->nSCTHits() + mu->nSCTDeadSensors() > 4);
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
  
  Bool_t dolowmass = GetDoLowMass();
  
  //// Muons
  if (lep->flavor() == Analysis::ChargedLepton::MUON) {
    D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
    
    if (lep->family() == getMuonFamily() || (lep->family() == Muon::CALO && DoCaloMuons)) lep->set_lastcut(HllqqMuonQuality::family);
    else return kFALSE;
    
    
    if ((lep->family() == Muon::MUID  && mu->tight()   == 1) ||
	(lep->family() == Muon::STACO && (mu->author() == 6  || mu->author() == 7)) ||
	(lep->family() == Muon::CALO  && DoCaloMuons && mu->author() == 16 && (mu->caloMuonIdTag() > 10 || mu->caloLRLikelihood() > 0.9))) lep->set_lastcut(HllqqMuonQuality::quality);
    else return kFALSE;
    
    
    if ((lep->family() != Muon::CALO                              && TMath::Abs(mu->eta()) < 2.7) || 
	(lep->family() == Muon::CALO && DoCaloMuons               && TMath::Abs(mu->eta()) < 0.1) || 
	(lep->family() != Muon::CALO && mu->isStandAloneMuon()==1 && TMath::Abs(mu->eta()) >2.5 && TMath::Abs(mu->eta()) < 2.7)) lep->set_lastcut(HllqqMuonQuality::eta);
    else return kFALSE;
    
    
    if(dolowmass)
      {
	if ((lep->family() != Muon::CALO                && lep->Get4Momentum()->Pt() > 7000.) ||
	    (lep->family() == Muon::CALO && DoCaloMuons && lep->Get4Momentum()->Pt() > 15000.)) lep->set_lastcut(HllqqMuonQuality::pt);
	else return kFALSE;
      }
    else if (!dolowmass)
      {
	if ((lep->family() != Muon::CALO                && lep->Get4Momentum()->Pt() > 10000.) ||
	    (lep->family() == Muon::CALO && DoCaloMuons && lep->Get4Momentum()->Pt() > 20000.))   lep->set_lastcut(HllqqMuonQuality::pt);
	else return kFALSE;
      }
    
    
    if (isMCPMuon(lep)) lep->set_lastcut(HllqqMuonQuality::MCP);
    else return kFALSE;
    
    
    if (TMath::Abs(mu->d0_exPV()) < 1. || mu->isStandAloneMuon()) lep->set_lastcut(HllqqMuonQuality::cosmic);
    else return kFALSE;
    
    
    if (TMath::Abs(mu->z0_exPV()) < 10. || mu->isStandAloneMuon()) lep->set_lastcut(HllqqMuonQuality::z0);
    else return kFALSE;
    
    
    
    if(!GetDoQCDSelection()) {
      
      if(dolowmass)
	{
	  //d0 Significance
	  if (TMath::Abs(lep->d0() / lep->d0_sig()) < 3.5) lep->set_lastcut(HllqqMuonQuality::d0Sig);
	  else return kFALSE;
	}
      else lep->set_lastcut(HllqqMuonQuality::d0Sig);
      
      if ((mu->ptcone20()/mu->pt())<.1) lep->set_lastcut(HllqqMuonQuality::Isolation);
      else return kFALSE;
      
    }
    
    return kTRUE;
  }
  
  
  
  //// Electrons
  else if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
    
    if (lep->family() == getElectronFamily()) lep->set_lastcut(HllqqElectronQuality::family);
    else return kFALSE;
    
    
    if (el->author() == 1 || el->author() == 3) lep->set_lastcut(HllqqElectronQuality::algorithm);
    else return kFALSE;
    
    
    if (analysis_version() == "rel_17") { // rel. 17
      
      if(dolowmass)
	{
	  if (el->tightPP() == 1 && !is4lepGood) lep->set_lastcut(HllqqElectronQuality::quality);
	  
	  else return kFALSE;
	}
      else if(!dolowmass)
	{
	  if (el->mediumPP() == 1) lep->set_lastcut(HllqqElectronQuality::quality);
	  
	  else return kFALSE;
	}
    } // rel. 17
    
    
    else if (analysis_version() == "rel_17_2") { // rel. 17.2
      
      if(dolowmass)
	{
	  if (/*el->tightPP()== 1*/
	      isTightPlusPlus(el->etas2(),
			      el->cl_E() / TMath::CosH(el->etas2()),
			      el->f3(),
			      el->Ethad() / (el->cl_E() / TMath::CosH(el->etas2())),
			      el->Ethad1() / (el->cl_E() / TMath::CosH(el->etas2())),
			      el->reta(),
			      el->weta2(),
			      el->f1(),
			      el->wstot(),
			      el->emaxs1() + el->Emax2() > 0 ?(el->emaxs1() - el->Emax2()) / (el->emaxs1() + el->Emax2()):0.,
			      el->deltaeta1(),
			      el->trackd0_physics(),
			      el->TRTHighTOutliersRatio(),
			      el->nTRTHits(),
			      el->nTRTOutliers(),
			      el->nSiHits(),
			      el->nSCTOutliers() + el->nPixelOutliers(),
			      el->nPixHits(),
			      el->nPixelOutliers(),
			      el->nBLHits(),
			      el->nBLayerOutliers(),
			      el->expectBLayerHit(),
			      el->cl_E() * TMath::Abs(el->trackqoverp()),
			      el->deltaphi2(), 
			      el->isEM() & (1 << 1),
			      egammaMenu::eg2012,
			      false,
			      false)) lep->set_lastcut(HllqqElectronQuality::quality);
	  
	  else return kFALSE;
	}
      else if(!dolowmass)
	{
	  if (/*el->mediumPP() == 1*/ 
	      isLoosePlusPlus(el->etas2(),el->cl_E() / TMath::CosH(el->etas2()),el->Ethad() / (el->cl_E() / TMath::CosH(el->etas2())),
			      el->Ethad1() / (el->cl_E() / TMath::CosH(el->etas2())),
			      el->reta(),
			      el->weta2(),
			      el->f1(),
			      el->wstot(),
			      el->emaxs1() + el->Emax2() > 0 ?(el->emaxs1() - el->Emax2()) / (el->emaxs1() + el->Emax2()):0.,
			      el->deltaeta1(),
			      el->nSiHits(),
			      el->nSCTOutliers() + el->nPixelOutliers(),
			      el->nPixHits(),
			      el->nPixelOutliers(),
			      egammaMenu::eg2012,
			      false,
			      false)) lep->set_lastcut(HllqqElectronQuality::quality);
	  
	  else return kFALSE;
	}
    } // rel. 17.2
    
    
    if (TMath::Abs(el->cl_eta()) < 2.47) lep->set_lastcut(HllqqElectronQuality::eta);
    else return kFALSE;
    
    
    if(dolowmass)
      {	
	if (lep->Get4Momentum()->Et() > 7000) lep->set_lastcut(HllqqElectronQuality::Et);
	else return kFALSE;
      } 
    else if(!dolowmass)
      {
	if (lep->Get4Momentum()->Et() > 10000) lep->set_lastcut(HllqqElectronQuality::Et);
	else return kFALSE;	
      }   
    
    
    if ((el->OQ() & 1446) == 0) lep->set_lastcut(HllqqElectronQuality::objectquality);
    else return kFALSE;
    
    if(dolowmass) 
      {
	if (TMath::Abs(lep->z0()) < 10.) lep->set_lastcut(HllqqElectronQuality::z0);
	else return kFALSE;
      }
    else
      {
	lep->set_lastcut(HllqqElectronQuality::z0);
	//else return kFALSE;
      }
    
    
    if(!GetDoQCDSelection()) {
      
      if(dolowmass) {
	if((TMath::Abs(lep->d0()) / lep->d0_sig())<6.5) lep->set_lastcut(HllqqElectronQuality::d0Sig);
	else return kFALSE;
      }
      else lep->set_lastcut(HllqqElectronQuality::d0Sig);    
      
      if ((el->ptcone20()/el->pt()) < 0.1) lep->set_lastcut(HllqqElectronQuality::Isolation);
      else return kFALSE;    
    }
    
    return kTRUE;  
  }
  
  
  else {
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
	((lep->family() != Muon::CALO && lep->Get4Momentum()->Pt() > 7000  && dolowmass)                || 
	 (lep->family() == Muon::CALO && lep->Get4Momentum()->Pt() > 15000 && dolowmass && DoCaloMuons) ||
	 (lep->family() != Muon::CALO && lep->Get4Momentum()->Pt() > 20000 && !dolowmass)               ||
	 (lep->family() == Muon::CALO && lep->Get4Momentum()->Pt() > 20000 && !dolowmass && DoCaloMuons))) lep->set_lastcut(HllqqMuonQuality::quality/*Kinematics*/);
    else return kFALSE;
    
    if(cut_leptons) Muon2++;
    
    
    if ((lep->family() == Muon::MUID && mu->tight() == 1) ||
	(lep->family() == Muon::STACO && (mu->author() == 6 || mu->author() == 7) && (mu->isCombinedMuon()==1 || mu->isSegmentTaggedMuon()==1)) ||
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
  
  D3PDReader::JetD3PDObjectElement *Jet = jet->GetJet();
  
  Bool_t dolowmass = GetDoLowMass();
  
  if(Jet->isBadLooseMinus() == 0) jet->set_lastcut(HllqqJetQuality::jetCleaning);
  else return kFALSE;
  
  
  if(dolowmass)
    {
      if (jet->rightpt()/1000. > 20. && TMath::Abs(jet->righteta()) < 4.5 && jet->rightE()>0) jet->set_lastcut(HllqqJetQuality::kinematics);
      else return kFALSE;
    }  
  else if (!dolowmass)
    {
      if (jet->rightpt()/1000. > 20. && TMath::Abs(jet->righteta()) < 4.5 && jet->rightE()>0) jet->set_lastcut(HllqqJetQuality::kinematics);
      else return kFALSE;
    }
  
  Float_t jvtxf_cut = (analysis_version() == "rel_17") ? 0.75 : 0.50;
  
  if(dolowmass)
    {
      if((TMath::Abs(Jet->jvtxf()) > jvtxf_cut)/* || (jet->rightpt()/1000.< 50. && TMath::Abs(jet->righteta())<2.4)*/) jet->set_lastcut(HllqqJetQuality::Pileup);
      else return kFALSE;
    }  
  else if (!dolowmass)
    {
      if(TMath::Abs(Jet->jvtxf()) > jvtxf_cut) jet->set_lastcut(HllqqJetQuality::Pileup);
      else return kFALSE;
    }  
  
  
  if(dolowmass)
    {
      Bool_t hot_tile(kTRUE);
      
      if (!isMC()) {
	if (ntuple->eventinfo.RunNumber() == 202660
	    || ntuple->eventinfo.RunNumber() == 202668
	    || ntuple->eventinfo.RunNumber() == 202712
	    || ntuple->eventinfo.RunNumber() == 202740
	    || ntuple->eventinfo.RunNumber() == 202965
	    || ntuple->eventinfo.RunNumber() == 202987
	    || ntuple->eventinfo.RunNumber() == 202991
	    || ntuple->eventinfo.RunNumber() == 203027) {
	  Float_t j_fmax = Jet->fracSamplingMax();
	  Float_t j_smax = Jet->SamplingMax();
	  Float_t j_eta = jet->righteta();
	  Float_t j_phi = jet->rightphi();
	  
	  Bool_t etaphi28(kFALSE);
	  if (j_eta > -0.2 && j_eta < -0.1 && j_phi > 2.65 && j_phi < 2.75) etaphi28 = kTRUE;
	  
	  if (j_fmax > 0.6 && j_smax < 13 && etaphi28) hot_tile = kFALSE;
	}
      }
      if (!hot_tile)
	return kFALSE;
    }
  
  return kTRUE; 
}


Bool_t HiggsllqqAnalysis::execute_analysis()
{
  
  // internal utility counters/flags
  m_processedEntries++;
  
  
  // bugfix for mc12a samples produced using FRONTIER
  if (analysis_version() == "rel_17_2" && isMC()) {
    if (getTriggerInfo("L1_RD0_FILLED") != 1) return kTRUE;
  }
  
  
  // update generated entries' histo
  m_generatedEntriesHisto->Fill("raw", 1);
  m_generatedEntriesHisto->Fill("Sherpa_Veto",(!SherpaPt0Veto()) ? 1 : 0);  
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix", (!hasPowHegZZBug()) ? 1 : 0);
  m_generatedEntriesHisto->Fill("PowHegZZ_bugfix_EvWeight", (!hasPowHegZZBug() && isMC()) ? ntuple->eventinfo.mc_event_weight() : 0);
  m_generatedEntriesHisto->Fill("with_ggF_as_well", (!hasPowHegZZBug() && isMC()) ? ntuple->eventinfo.mc_event_weight() * getggFWeight() : 0);
  m_generatedEntriesHisto->Fill("with_pu_and_vertex", (!hasPowHegZZBug() && isMC()) ? getEventWeight() * getggFWeight() * getVertexZWeight() : 0);
  
  
  
  for (UInt_t i = 0; i < m_Channels.size(); i++) {    
    
    UInt_t chan = m_Channels.at(i);
    
    setChannel(chan);
    
    for(int sel=0; sel<2; sel++) //Looping on the Low/High Mass selections 
      {
	if(sel==0)
	  {  //Set the Low or High Mass Analysis
	    SetDoLowMass(kTRUE);
	  }
	else
	  {  //Set the Low or High Mass Analysis
	    SetDoLowMass(kFALSE);
	  }
	
	
	m_called_getGoodLeptons = kFALSE;      
	m_called_getGoodObjects = kFALSE;
	ResetReducedNtupleMembers();
	// TestSelection Reset 
	ResetAnalysisOutputBranches(&m_outevent);
	
	
	Int_t last_event = getLastCutPassed();	
	
	
	if(sel==Print_low_OR_high) //Printing of the cutflow shows Low==0 OR High==1 !!!!
	  {
	    
	    // update cutflows
	    m_EventCutflow[chan].addCutCounter(last_event, 1);
	    m_EventCutflow_rw[chan].addCutCounter(last_event, 1.*getSFWeight()*getPileupWeight()*ntuple->eventinfo.mc_event_weight());
	    
	    
	    // Fill the cutflow histograms	    
	    if(chan==2)
	      h_cutflow->Fill(last_event+2,1);
	    if(chan==0)
	      h_cutflow_weight->Fill(last_event+2,1/*getEventWeight()*/); //error: We are using the two histo to separate the cutflow per channel (MU2-E2)
	    
	    
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
	    
	  } //End Printing of the cutflow shows Low==0 OR High==1 !!!!
	
	
	
	//Set the QCD flag FALSE
	SetDoQCDSelection(kFALSE);
	
	if(chan!=1)// error: repare the mixing MUE channel
	  {
	    //Filling of the equivalent qqll tree (2011)
	    FillReducedNtuple(last_event,chan);    
	    
	    
	    // TestSelection Filling candidate struct
	    FillAnalysisOutputTree(&m_outevent,last_event,chan);
	  }	
	
	
	last_event = -1;
	
	
	m_Muons.clear();
	m_Electrons.clear();
	m_Jets.clear();
	m_GoodMuons.clear();
	m_GoodElectrons.clear();
	m_GoodJets.clear();
	
	
	/////////////////////
	// QCD selection - only for data
	/////////////////////
	if(!isMC()) {
	  
	  m_called_getGoodLeptons = kFALSE;      
	  m_called_getGoodObjects = kFALSE;
	  ResetReducedNtupleMembers();
	  // TestSelection Reset 
	  ResetAnalysisOutputBranches(&m_outevent);
	  SetDoQCDSelection(kTRUE);
	  last_event = getLastCutPassed();
	  
	  if(chan!=1)// error: repare the mixing MUE channel
	    {
	      // Filling of the equivalent qqll tree (2011)
	      FillReducedNtuple(last_event,chan);    
	      
	      
	      // TestSelection Filling candidate struct
	      FillAnalysisOutputTree(&m_outevent,last_event,chan);
	    }	
	  
	  
	  last_event = -1;
	  
	  m_Muons.clear();
	  m_Electrons.clear();
	  m_Jets.clear();
	  m_GoodMuons.clear();
	  m_GoodElectrons.clear();
	  m_GoodJets.clear();
	  
	} //end QCD selection
	
      } //End of loop into the Low/high Selection
    
  } // end of loop into the different channels
  
  
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
  
  
  return kTRUE;
}


Bool_t HiggsllqqAnalysis::finalize_analysis()
{  
  
  // print the number of processed entries
  Info("finalize_analysis", "Analysis terminated, processed %lld entries out of %lld available entries", m_processedEntries, m_entriesInChain);
  
  //Filling the First bin of the cutflows
  h_cutflow->Fill("Entries",m_processedEntries);
  h_cutflow_weight->Fill("Entries",m_processedEntries);
  
  
  // print the cutflows
  for (UInt_t i = 0; i < m_Channels.size(); i++) {
    m_EventCutflow[i].print();
    m_EventCutflow_rw[i].print();
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
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2012_C").second) {
      // C
      return DataPeriod::y2012_C;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2012_D").second) {
      // D
      return DataPeriod::y2012_D;
    } else if (ntuple->eventinfo.RunNumber() <= DataPeriodTool::GetPeriodRange("y2012_E").second) {
      // E
      return DataPeriod::y2012_E;
    } else {
      // SAFETY : consider last runs as belonging to the latest period
      return DataPeriod::y2012_E;
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
    double mu_i = (isMC() && ntuple->eventinfo.lbn()==1 && int(ntuple->eventinfo.averageIntPerXing()+0.5)==1) ? 0. : ntuple->eventinfo.averageIntPerXing();
    result = m_PileupReweighter->GetCombinedWeight(ntuple->eventinfo.RunNumber(), ntuple->eventinfo.mc_channel_number(), mu_i);
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
      
      m_PileupReweighter->SetRandomSeed(314159 + ntuple->eventinfo.mc_channel_number() * 2718 + ntuple->eventinfo.EventNumber());
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


void HiggsllqqAnalysis::initCrossSections()  //To be updated . Error
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
  m_SignalSampleMass[160402] = 110;
  m_SignalSampleMass[160403] = 115;
  m_SignalSampleMass[160404] = 120;
  m_SignalSampleMass[160405] = 125;
  m_SignalSampleMass[160406] = 130;
  m_SignalSampleMass[160407] = 135;
  m_SignalSampleMass[160408] = 140;
  m_SignalSampleMass[160409] = 145;
  m_SignalSampleMass[160410] = 150;
  m_SignalSampleMass[160411] = 155;
  m_SignalSampleMass[160412] = 160;
  m_SignalSampleMass[160413] = 165;
  m_SignalSampleMass[160414] = 170;
  m_SignalSampleMass[160415] = 175;
  m_SignalSampleMass[160416] = 180;
  m_SignalSampleMass[160417] = 185;
  m_SignalSampleMass[160418] = 190;
  m_SignalSampleMass[160419] = 195;
  m_SignalSampleMass[160420] = 200;
  m_SignalSampleMass[160421] = 220;
  m_SignalSampleMass[160422] = 240;
  m_SignalSampleMass[160423] = 260;
  m_SignalSampleMass[160424] = 280;
  m_SignalSampleMass[160425] = 300;
  m_SignalSampleMass[160426] = 320;
  m_SignalSampleMass[160427] = 340;
  m_SignalSampleMass[160428] = 360;
  m_SignalSampleMass[160429] = 380;
  m_SignalSampleMass[160430] = 400;
  /*
    m_SignalSampleMass[160431] = 420;
    m_SignalSampleMass[160432] = 440;
    m_SignalSampleMass[160433] = 460;
    m_SignalSampleMass[160434] = 480;
    m_SignalSampleMass[160435] = 500;
    m_SignalSampleMass[160436] = 520;
    m_SignalSampleMass[160437] = 540;
    m_SignalSampleMass[160438] = 560;
    m_SignalSampleMass[160439] = 580;
    m_SignalSampleMass[160440] = 600;
  */
}


Float_t HiggsllqqAnalysis::getCrossSectionWeight() //To be updated . Error
{
  Float_t result(-9999.9);
  
  if (isMC()) {
    UInt_t chan = ntuple->eventinfo.mc_channel_number();
    // ZZ cross section with M_ZZ dependent K factor
    if (chan == 109292 || chan == 109291) {
      Float_t xsec = CrossSections::GetBkgCrossSection(chan, kFALSE);
      Float_t mcfm = GetMzzWeightFromMCFM(126000/*getTruthZZMass()*/ / 1000.); //To be corrected!! Error
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


Float_t HiggsllqqAnalysis::getSFWeight()
{
  Float_t result = 1.;
  if (isMC()) 
    {
      if (getChannel()==HiggsllqqAnalysis::E2 && m_GoodElectrons.size()==2) {
	result *= getLeptonWeight(m_GoodElectrons.at(0));
	result *= getLeptonWeight(m_GoodElectrons.at(1));
      }
      else if (getChannel()==HiggsllqqAnalysis::MU2 && m_GoodMuons.size()==2) {
	result *= getLeptonWeight(m_GoodMuons.at(0));
	result *= getLeptonWeight(m_GoodMuons.at(1));
      }
      else if (getChannel()==HiggsllqqAnalysis::MUE && m_GoodMuons.size()==1 && m_GoodElectrons.size()==1) {
	result *= getLeptonWeight(m_GoodMuons.at(0));
	result *= getLeptonWeight(m_GoodElectrons.at(0));
      }
    }
  
  return result;
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
  
  Bool_t dolowmass = GetDoLowMass();
  
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
  
  
  int idx1 = 0;
  int idx2 = 1;
  
  if (DoKinematicFitter)
    {
      //cout<<".... Take care: Doing KF..."<<endl;
      idx1 = m_bestdijet->index1!=-1 ? m_bestdijet->index1 : 0;
      idx2 = m_bestdijet->index2!=-1 ? m_bestdijet->index2 : 1;
    }
  else
    {
      //cout<<".... Take care: Using two leading Jets..."<<endl;
      idx1 = 0;
      idx2 = 1;
    }
  
  
  //Filling Global variables with the index of the selected jets into the event to fill the TestSelection tree
  Pair_jet1 = idx1;
  Pair_jet2 = idx2;
  
  if(m_bestdijet->index1!=-1)
    {
      corr_jet_pt1 = m_bestdijet->refittedPt1;
      corr_jet_pt2 = m_bestdijet->refittedPt2;
      ChiSq        = m_bestdijet->chiSq;;
    }
  
  
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
  
  if((((m_bestdijet->index1!=-1)    && (m_bestdijet->index2!=-1) && (hadZ.M()>Mjj_low_min)  && (hadZ.M()<Mjj_low_max))  &&  dolowmass) 
     || (((m_bestdijet->index1!=-1) && (m_bestdijet->index2!=-1) && (hadZ.M()>Mjj_high_min) && (hadZ.M()<Mjj_high_max)) && !dolowmass))
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
  if(!GetDoLowMass())
    ptthr  = 20000.;
  
  bool Bad_event = kFALSE;
  
  Float_t dR_GoodEl_BadJet = -1.;
  
  D3PDReader::JetD3PDObject *jet_branch(0);
  if (getJetFamily() == 0)      jet_branch = &(ntuple->jet_akt4topoem);
  else if (getJetFamily() == 1) jet_branch = &(ntuple->jet_AntiKt4LCTopo);
  
  D3PDReader::ElectronD3PDObject *el_branch(0);
  if (getElectronFamily()      == Electron::GSF)   el_branch = &(ntuple->el_GSF);
  else if (getElectronFamily() == Electron::noGSF) el_branch = &(ntuple->el);
  
  
  for (Int_t i = 0; i < jet_branch->n(); i++) {

    Analysis::Jet *jet = new Analysis::Jet(&((*jet_branch)[i]));
    D3PDReader::JetD3PDObjectElement *Jet = jet->GetJet();
    
    applyChanges(jet);
    
    
    if((Jet->isBadLooseMinus()!=0) && jet->rightpt()>ptthr /*&& jet->righteta() < 2.5*/)
      {
	Bad_event = kTRUE;
	prebadevent++;
	
	for (Int_t i = 0; i < el_branch->n(); i++) 
	  {
	    Analysis::ChargedLepton *lep = new Analysis::ChargedLepton(&((*el_branch)[i]), getElectronFamily());
	    applyChanges(lep);
	    //ApplyChangesElectron(lep);
	    
	    if (isGood(lep)) {  
	      
	      D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();
	      
	      dR_GoodEl_BadJet = TMath::Sqrt(TMath::Power(jet->righteta() - el->cl_eta(), 2) + TMath::Power(TVector2::Phi_mpi_pi(jet->rightphi() - el->cl_phi()), 2));
	      
	      
	      /*if (dR_GoodEl_BadJet>-1. && dR_GoodEl_BadJet>0.4)
	      //index    Evtnum    pt      eta   isBadLooseMinus jvf   dr(j,e)
	      cout<<"\t"<<ntuple->eventinfo.EventNumber()
	      <<"\t"<<jet->rightpt()
	      <<"\t"<<jet->righteta()
	      <<"\t"<<Jet->isBadLooseMinus()
	      <<"\t"<<Jet->jvtxf()
	      <<"\t"<<dR_GoodEl_BadJet
	      <<"\t"<<endl;
	      */
	      
	      if(dR_GoodEl_BadJet<0.4 && dR_GoodEl_BadJet>-1.)
		Bad_event = kFALSE;
	    }
	  }
	
	if(Bad_event)
	  {
	    badevent++;
	    Bad_event = kFALSE;
	    
	    
	    return kTRUE;
	  }
      }
    
    dR_GoodEl_BadJet=-1.;
    
    
    Bool_t hot_tile(kTRUE);
    
    if (!isMC()) {
      if (ntuple->eventinfo.RunNumber()    == 202660
	  || ntuple->eventinfo.RunNumber() == 202668
	  || ntuple->eventinfo.RunNumber() == 202712
	  || ntuple->eventinfo.RunNumber() == 202740
	  || ntuple->eventinfo.RunNumber() == 202965
	  || ntuple->eventinfo.RunNumber() == 202987
	  || ntuple->eventinfo.RunNumber() == 202991
	  || ntuple->eventinfo.RunNumber() == 203027) {
	
	Float_t j_fmax = Jet->fracSamplingMax();
	Float_t j_smax = Jet->SamplingMax();
	Float_t j_eta = jet->righteta();
	Float_t j_phi = jet->rightphi();
	
	Bool_t etaphi28(kFALSE);
	if (j_eta > -0.2 && j_eta < -0.1 && j_phi > 2.65 && j_phi < 2.75) etaphi28 = kTRUE;
	
	if (j_fmax > 0.6 && j_smax < 13 && etaphi28) hot_tile = kFALSE;
      }
    }
    
    if (!hot_tile)
      return kFALSE;
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
  if(GetDoLowMass()) pass_cut = (getCorrectMETValue() < MET_low_cut);
  else pass_cut = (getCorrectMETValue() < MET_high_cut);
  
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
    if(GetMV1value(*jet_itr_i) > MV1_OP70){
      tags++;
    }
  }
  return tags;
}


//Method to calculate the DiJet invariant mass for the tagged Jets!
Bool_t HiggsllqqAnalysis::JetDimassTagged() {
  
  Bool_t dolowmass = GetDoLowMass();
  
  TLorentzVector j1;
  TLorentzVector j2;
  bool first =kTRUE;
  
  std::vector<Analysis::Jet *>::iterator jet_itr;
  
  for (jet_itr = m_GoodJets.begin(); jet_itr != m_GoodJets.end(); ++jet_itr) {
    
    if((GetMV1value(*jet_itr) > MV1_OP70) && first)
      {
	first =kFALSE;
	j1.SetPtEtaPhiM(b_rescaling*(*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
      }
    if((GetMV1value(*jet_itr) > MV1_OP70) && !first)
      {
	j2.SetPtEtaPhiM(b_rescaling*(*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
      }
  }
  
  TLorentzVector hadZ = j1 + j2;
  
  if(((hadZ.M() > Mjj_low_min) && (hadZ.M() < Mjj_low_max) && dolowmass) || ((hadZ.M() > Mjj_high_min) && (hadZ.M() < Mjj_high_max) && !dolowmass))
    {
      return kTRUE;
    }
  else
    return kFALSE;
}



//Method to calculate the DiJet invariant mass for the 1 tagged Jet and leading jet.
Bool_t HiggsllqqAnalysis::JetDimassOneTagged() {
  
  Bool_t dolowmass = GetDoLowMass();
  
  TLorentzVector j1;
  TLorentzVector j2;
  int bjet_index =-1;
  
  std::vector<Analysis::Jet *>::iterator jet_itr;
  int ii =-1;
  
  for (jet_itr = m_GoodJets.begin(); jet_itr != m_GoodJets.end(); ++jet_itr) {
    ii++;
    if((GetMV1value(*jet_itr) > MV1_OP70))
      {
	bjet_index = ii;
	j1.SetPtEtaPhiM(b_rescaling*(*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
      }
  }
  
  if(ii!=0)
    //    j2.SetPtEtaPhiM(b_rescaling*(*jet_itr)->rightpt(),(*jet_itr)->righteta(),(*jet_itr)->rightphi(),(*jet_itr)->Get4Momentum()->M());
    j2.SetPtEtaPhiM(b_rescaling*m_GoodJets.at(0)->rightpt(),m_GoodJets.at(0)->righteta(),m_GoodJets.at(0)->rightphi(),m_GoodJets.at(0)->Get4Momentum()->M());
  else
    j2.SetPtEtaPhiM(b_rescaling*m_GoodJets.at(1)->rightpt(),m_GoodJets.at(1)->righteta(),m_GoodJets.at(1)->rightphi(),m_GoodJets.at(1)->Get4Momentum()->M());
  
  
  TLorentzVector hadZ = j1 + j2;
  
  if(((hadZ.M() > Mjj_low_min) && (hadZ.M() < Mjj_low_max) && dolowmass) || ((hadZ.M() > Mjj_high_min) && (hadZ.M() < Mjj_high_max) && !dolowmass))
    {
      return kTRUE;
    }
  else
    return kFALSE;
}




//2nd September 2012
void HiggsllqqAnalysis::InitReducedNtuple() 
{
  // FLS tree
  m_jets_m           = new std::vector<Float_t>;
  m_jets_pt          = new std::vector<Float_t>;
  m_jets_eta         = new std::vector<Float_t>;
  m_jets_eta_det     = new std::vector<Float_t>;
  m_jets_phi         = new std::vector<Float_t>;
  m_jets_MV1         = new std::vector<Float_t>;
  m_jets_flavortruth = new std::vector<Float_t>;
  m_jets_jvtxf       = new std::vector<Float_t>;
  m_jets_nTrk        = new std::vector<Int_t>;
  m_jets_width       = new std::vector<Float_t>;
  m_jets_flavorpdg   = new std::vector<Int_t>;
  m_jets_Epdg        = new std::vector<Double_t>;
  m_lep_m            = new std::vector<Float_t>;
  m_lep_pt           = new std::vector<Float_t>;
  m_lep_eta          = new std::vector<Float_t>;
  m_lep_phi          = new std::vector<Float_t>;
  m_lep_charge       = new std::vector<Int_t>;
  m_lep_d0           = new std::vector<Float_t>;
  m_lep_sigd0        = new std::vector<Float_t>;
  m_lep_trackiso     = new std::vector<Float_t>;
  m_lep_caloiso      = new std::vector<Float_t>;
  m_lep_quality      = new std::vector<Int_t>;
  m_quark_m          = new std::vector<Float_t>;
  m_quark_pt         = new std::vector<Float_t>;
  m_quark_pdg        = new std::vector<Int_t>;
  m_quark_eta        = new std::vector<Float_t>;
  m_quark_phi        = new std::vector<Float_t>;
  m_quark_E          = new std::vector<Float_t>;
  
  
  // Just Intialization for HFOR calculation. Not to save into the ReducedNtuple
  if (isMC()) 
    {
      mc_n               = 0;
      mc_pt              = new std::vector<Float_t>;
      mc_pdgId           = new std::vector<Int_t>;
      mc_m               = new std::vector<Float_t>;
      mc_eta             = new std::vector<Float_t>;
      mc_phi             = new std::vector<Float_t>;
      mc_status          = new std::vector<Int_t>;
      mc_vx_barcode      = new std::vector<Int_t>;
      mc_child_index     = new std::vector<vector<Int_t> >;
      mc_parent_index    = new std::vector<vector<Int_t> >;
    }
  
  // tree
  m_reduced_ntuple = new TTree("higgs","reduced ntuple");
  m_reduced_ntuple->Branch("RunNumber",&m_run);
  m_reduced_ntuple->Branch("EventNumber",&m_event);
  m_reduced_ntuple->Branch("channel",&m_channel);
  m_reduced_ntuple->Branch("isqcdevent",&m_qcdevent);
  m_reduced_ntuple->Branch("islowevent",&m_low_event);
  m_reduced_ntuple->Branch("last_cut",&m_cut);
  m_reduced_ntuple->Branch("jet_m",&m_jets_m);
  m_reduced_ntuple->Branch("jet_pt",&m_jets_pt);
  m_reduced_ntuple->Branch("jet_eta",&m_jets_eta);
  m_reduced_ntuple->Branch("jet_eta_det",&m_jets_eta_det);
  m_reduced_ntuple->Branch("jet_phi",&m_jets_phi);
  m_reduced_ntuple->Branch("jet_MV1",&m_jets_MV1);
  m_reduced_ntuple->Branch("jet_flavortruth",&m_jets_flavortruth);
  m_reduced_ntuple->Branch("jet_jvtxf",&m_jets_jvtxf);
  m_reduced_ntuple->Branch("jet_nTrk",&m_jets_nTrk);
  m_reduced_ntuple->Branch("jet_width",&m_jets_width);
  m_reduced_ntuple->Branch("jet_flavorpdg",&m_jets_flavorpdg);
  m_reduced_ntuple->Branch("jet_Epdg",&m_jets_Epdg);
  m_reduced_ntuple->Branch("lep_m",&m_lep_m);
  m_reduced_ntuple->Branch("lep_pt",&m_lep_pt);
  m_reduced_ntuple->Branch("lep_eta",&m_lep_eta);
  m_reduced_ntuple->Branch("lep_phi",&m_lep_phi);
  m_reduced_ntuple->Branch("lep_charge",&m_lep_charge);
  m_reduced_ntuple->Branch("lep_d0",&m_lep_d0);
  m_reduced_ntuple->Branch("lep_sigd0",&m_lep_sigd0);
  m_reduced_ntuple->Branch("lep_trackiso",&m_lep_trackiso);
  m_reduced_ntuple->Branch("lep_caloiso",&m_lep_caloiso);
  m_reduced_ntuple->Branch("lep_quality",&m_lep_quality);
  m_reduced_ntuple->Branch("lep_chargeproduct",&m_lep_chargeproduct);
  m_reduced_ntuple->Branch("met_met",&m_met_met);
  m_reduced_ntuple->Branch("met_phi",&m_met_phi);
  m_reduced_ntuple->Branch("met_sumet",&m_met_sumet);
  m_reduced_ntuple->Branch("weight",&m_weight);
  m_reduced_ntuple->Branch("mu",&m_mu);
  m_reduced_ntuple->Branch("NPV",&m_NPV);
  m_reduced_ntuple->Branch("truthH_pt",&m_truthH_pt);
  m_reduced_ntuple->Branch("ggFweight",&m_ggFweight);
  
  //Trigger related scale factors and flag to tell single/dilepton triggers
  m_reduced_ntuple->Branch("trig_SF",&m_trig_SF);
  m_reduced_ntuple->Branch("trig_SF2",&m_trig_SF2);
  m_reduced_ntuple->Branch("trig_SFC",&m_trig_SFC);
  m_reduced_ntuple->Branch("trig_flag",&m_trig_flag);
  //
  m_reduced_ntuple->Branch("Entries",&m_Entries);
  m_reduced_ntuple->Branch("HFOR",&m_HFOR);
  m_reduced_ntuple->Branch("quark_m",&m_quark_m);
  m_reduced_ntuple->Branch("quark_pt",&m_quark_pt);
  m_reduced_ntuple->Branch("quark_pdg",&m_quark_pdg);
  m_reduced_ntuple->Branch("quark_eta",&m_quark_eta);
  m_reduced_ntuple->Branch("quark_phi",&m_quark_phi);
  m_reduced_ntuple->Branch("quark_E",&m_quark_E);
  
  ResetReducedNtupleMembers();  
}

void HiggsllqqAnalysis::ResetReducedNtupleMembers() 
{
  // FLS tree
  m_jets_m->clear();
  m_jets_pt->clear();
  m_jets_eta->clear();
  m_jets_eta_det->clear();
  m_jets_phi->clear();
  m_jets_MV1->clear();
  m_jets_flavortruth->clear();
  m_jets_jvtxf->clear();
  m_jets_nTrk->clear();
  m_jets_width->clear();
  m_jets_flavorpdg->clear();
  m_jets_Epdg->clear();
  m_lep_m->clear();
  m_lep_pt->clear();
  m_lep_eta->clear();
  m_lep_phi->clear();
  m_lep_charge->clear();
  m_lep_d0->clear();
  m_lep_sigd0->clear();
  m_lep_trackiso->clear();
  m_lep_caloiso->clear();
  m_lep_quality->clear();
  m_quark_m->clear();
  m_quark_pt->clear();
  m_quark_pdg->clear();
  m_quark_eta->clear();
  m_quark_phi->clear();
  m_quark_E->clear();
  
  m_lep_chargeproduct = 0;
  m_run       = -1;
  m_event     = -1;
  m_weight    = 1.;
  m_mu        = 1.;
  m_NPV       = 0;
  m_truthH_pt = -1.;
  m_ggFweight = 1.;
  m_cut       = -1;
  m_channel   = -1;
  m_qcdevent  = -1;
  m_low_event  = -1;
  m_met_met   = -9999.;
  m_met_phi   = -9999.;
  m_met_sumet = -9999.;
  m_Entries   = -9999;
  m_HFOR      = -9999;
}


void HiggsllqqAnalysis::FillReducedNtuple(Int_t cut, UInt_t channel)
{
  Int_t minimum_cut = 0;
  
  if(GetDoQCDSelection()) minimum_cut = HllqqCutFlow::PtLeptons;//MET;//NumberOfLeptons;//DiJetMass;
  else minimum_cut = HllqqCutFlow::PtLeptons;
  
  std::pair<float,float> jetsf;
  jetsf.first = 1.;
  jetsf.second = 1.;
  
  
  if(cut >= minimum_cut)
    {
      if (!isMC()) m_run         = ntuple->eventinfo.RunNumber();     
      if (isMC())  m_run         = ntuple->eventinfo.mc_channel_number();
      m_event       = ntuple->eventinfo.EventNumber();
      m_cut         = cut;
      m_weight      =  1.;
      m_mu          =  1.;
      m_truthH_pt   = -1.;
      m_ggFweight   =  1.;
      m_NPV         =  0;
      m_channel     = channel;
      
      
      if(channel == HiggsllqqAnalysis::MU2)
	{
	  m_lep_charge->push_back((m_GoodMuons.at(0))->charge());
	  m_lep_charge->push_back((m_GoodMuons.at(1))->charge());
	  
	  for (std::vector<Analysis::ChargedLepton*>::iterator mu_itr = m_GoodMuons.begin(); mu_itr != m_GoodMuons.end(); ++mu_itr)
	    {
	      Analysis::ChargedLepton *mu_a = (*mu_itr);
	      D3PDReader::MuonD3PDObjectElement *b_mu = mu_a->GetMuon();
	      
	      m_lep_m->push_back(mu_pdg_mass);
	      m_lep_pt->push_back(b_mu->pt());
	      m_lep_eta->push_back(b_mu->eta());
	      m_lep_phi->push_back(b_mu->phi());
	      m_lep_charge->push_back(b_mu->charge());
	      m_lep_d0->push_back(mu_a->d0());
	      m_lep_sigd0->push_back(mu_a->d0_sig());
	      m_lep_trackiso->push_back(b_mu->ptcone20()/b_mu->pt());
	      m_lep_caloiso->push_back(b_mu->etcone20()/b_mu->pt());
	      
	      if (mu_a->family() == Muon::MUID)
		m_lep_quality->push_back(1);
	      else if (mu_a->family() == Muon::STACO)
		m_lep_quality->push_back(2);
	      else if (mu_a->family() == Muon::CALO)
		m_lep_quality->push_back(3);
	    }
	}
      
      
      else if(channel == HiggsllqqAnalysis::E2)
	{
	  m_lep_charge->push_back((m_GoodElectrons.at(0))->charge());
	  m_lep_charge->push_back((m_GoodElectrons.at(1))->charge());
	  
	  for (std::vector<Analysis::ChargedLepton*>::iterator el_itr = m_GoodElectrons.begin(); el_itr != m_GoodElectrons.end(); ++el_itr)
	    {
	      Analysis::ChargedLepton *el_a = (*el_itr);
	      D3PDReader::ElectronD3PDObjectElement *b_el = el_a->GetElectron();
	      
	      m_lep_m->push_back(el_pdg_mass);
	      m_lep_pt->push_back(el_a->Get4Momentum()->Et());
	      m_lep_eta->push_back(b_el->cl_eta());
	      m_lep_phi->push_back(b_el->cl_phi());
	      m_lep_charge->push_back(b_el->charge());
	      m_lep_d0->push_back(el_a->d0());
	      m_lep_sigd0->push_back(el_a->d0_sig());
	      m_lep_trackiso->push_back(b_el->ptcone20()/b_el->pt());
	      m_lep_caloiso->push_back(b_el->Etcone20()/b_el->pt());
	      
	      if(b_el->tightPP() == 1)
		m_lep_quality->push_back(3);
	      else if(b_el->mediumPP() == 1)
		m_lep_quality->push_back(2);
	      else if(b_el->loosePP() == 1)
		m_lep_quality->push_back(1);
	    }
	}
      
      
      else if(channel == HiggsllqqAnalysis::MUE)
	{ //To be implemented!! error
	  for (std::vector<Analysis::ChargedLepton*>::iterator el_itr = m_GoodElectrons.begin(); el_itr != m_GoodElectrons.end(); ++el_itr)
	    {
	    }
	  
	  for (std::vector<Analysis::ChargedLepton*>::iterator mu_itr = m_GoodMuons.begin(); mu_itr != m_GoodMuons.end(); ++mu_itr)
	    {
	    }
	}
      
      
      std::vector<Analysis::Jet *>::iterator jet_itr;
      for (jet_itr = m_GoodJets.begin(); jet_itr != m_GoodJets.end(); ++jet_itr)
	{	
	  D3PDReader::JetD3PDObjectElement *Jet = (*jet_itr)->GetJet();
	  
	  m_jets_m->push_back((*jet_itr)->Get4Momentum()->M());
	  m_jets_pt->push_back((*jet_itr)->rightpt());
	  m_jets_eta->push_back((*jet_itr)->righteta());
	  m_jets_eta_det->push_back(Jet->emscale_eta());
	  m_jets_phi->push_back((*jet_itr)->rightphi());
	  m_jets_MV1->push_back(GetMV1value(*jet_itr));
	  m_jets_jvtxf->push_back(Jet->jvtxf());
	  m_jets_nTrk->push_back(Jet->nTrk());
	  m_jets_width->push_back(Jet->WIDTH());
	  
	  if (isMC()) 
	    {
	      m_jets_flavortruth->push_back(Jet->flavor_truth_label());
	      m_jets_flavorpdg->push_back(GetFlavour(*jet_itr).first);
	      m_jets_Epdg->push_back(GetFlavour(*jet_itr).second);
	    }
	}
      
      if (isMC())
	{
	  //Truth quark information
	  for (Int_t i = 0; i < ntuple->mc.n(); i++)
	    {
	      D3PDReader::TruthParticleD3PDObjectElement *p = &(ntuple->mc[i]);
	      TLorentzVector Tp = CommonTools::getVector(p);
	      
	      if(isHeavyJet(p->pdgId()))
		{
		  m_quark_m->push_back(Tp.M());
		  m_quark_pt->push_back(p->pt());
		  m_quark_pdg->push_back(p->pdgId());
		  m_quark_eta->push_back(p->eta());
		  m_quark_phi->push_back(p->phi());
		  m_quark_E->push_back(Tp.E());
		}
	    }
	}
      
      m_lep_chargeproduct = m_lep_charge->at(0)*m_lep_charge->at(1);
      m_met_met           = getCorrectMETValue();
      m_met_phi           = ntuple->MET_RefFinal.phi();
      m_met_sumet         = ntuple->MET_RefFinal.sumet(); 
      m_qcdevent          = GetDoQCDSelection() ? 1 : 0;
      m_low_event         = GetDoLowMass() ? 1 : 0;
      m_NPV               = getNumberOfGoodVertices();
      
      
      // all jet-related weights are missing
      
      if (isMC())
	{
	  m_weight *= getSFWeight()*getPileupWeight()*ntuple->eventinfo.mc_event_weight();
	  m_mu      = (isMC() && ntuple->eventinfo.lbn()==1 && int(ntuple->eventinfo.averageIntPerXing()+0.5)==1) ? 0. : ntuple->eventinfo.averageIntPerXing();
	}
      
      
      // reset and fill trigger flag word
      m_trig_flag=0;
      if(channel == HiggsllqqAnalysis::MU2)
	{	
	  //Fill flag to distinguish single & double lepton trigger
	  if (passesSingleMuonTrigger()) m_trig_flag |= 1<<0;
	  if (passesDiMuonTrigger())     m_trig_flag |= 1<<1;
	  
	} 
      else if(channel == HiggsllqqAnalysis::E2)
	{ 
	  //Fill flag to distinguish single & double lepton trigger
	  if (passesSingleElectronTrigger()) m_trig_flag |= 1<<0;
	  if (passesDiElectronTrigger())     m_trig_flag |= 1<<1;
	} 
      
      /*
	std::pair<double, double> triggerSF  (1.,0.);
	if (isMC()) { 
	triggerSF = getCandidateTriggerSF(0);
	m_trig_SF=triggerSF.first;
	triggerSF = getCandidateTriggerSF(2);
	m_trig_SF2=triggerSF.first;
	triggerSF = getCandidateTriggerSF(3);
	m_trig_SFC=triggerSF.first;
	}
      */
      
      // FILL ggF weight variables
      m_ggFweight = getggFWeight();
      m_truthH_pt = getTruthHiggsPt();
      m_Entries   = fChain->GetEntries();    
      m_HFOR      = HFOR_value; 
      
      m_reduced_ntuple->Fill();
      ResetReducedNtupleMembers();
    }
}


//Method to obtain the flavor of the MC jet, using the true information.
pair <Int_t,Double_t> HiggsllqqAnalysis::GetFlavour(Analysis::Jet *jet)
{
  
  Float_t pt_leading = 0;
  Float_t dr         = 0.4;
  Int_t flavour      = -999;
  Double_t E_leading = -1.;
  pair <Int_t,Double_t> mcinfo;
  
  
  for (Int_t i = 0; i < ntuple->mc.n(); i++) 
    {
      D3PDReader::TruthParticleD3PDObjectElement *p = &(ntuple->mc[i]);
      D3PDReader::JetD3PDObjectElement *this_jet    = jet->GetJet();
      Int_t pdg = p->pdgId();
      
      if(isGluonJet(pdg) || isHeavyJet(pdg) || isLightJet(pdg))
	{
	  TLorentzVector Tp = CommonTools::getVector(p);
	  Float_t tmp_dr = TMath::Sqrt(TMath::Power(this_jet->EtaOrigin() - p->eta(), 2) + TMath::Power(TVector2::Phi_mpi_pi(this_jet->PhiOrigin() - p->phi()), 2));
	  
	  if(tmp_dr<dr) 
	    {
	      if(isHeavyJet(pdg))
		{
		  pt_leading=p->pt();
		  flavour=pdg;
		  E_leading=Tp.E();
		}
	      else if(p->pt()>pt_leading)
		{
		  pt_leading=p->pt();
		  flavour=pdg;
		  E_leading=Tp.E();
		}
	    }
	}
    }
  
  mcinfo.first=flavour;
  mcinfo.second=E_leading;
  
  return mcinfo;
}


Bool_t HiggsllqqAnalysis::isHeavyJet(Int_t pdg)
{
  Bool_t flavor=kFALSE;
  if(pdg==4 || pdg==5 || pdg==6 || pdg==-4 || pdg==-5 || pdg==-6){flavor=kTRUE;}
  else{flavor=kFALSE;}
  return flavor;
}


Bool_t HiggsllqqAnalysis::isLightJet(Int_t pdg)
{
  Bool_t flavor=kFALSE;
  if(pdg==1 || pdg==2 || pdg==3 || pdg==-1 || pdg==-2 || pdg==-3){flavor=kTRUE;}
  else{flavor=kFALSE;}
  return flavor;
}


Bool_t HiggsllqqAnalysis::isGluonJet(Int_t pdg)
{
  Bool_t flavor=kFALSE;
  if(pdg==9 || pdg==21){flavor=kTRUE;}
  else{flavor=kFALSE;}
  return flavor;
}


void HiggsllqqAnalysis::InitMasses() {
  fSampleMass.clear();
  
  // sample masses (in GeV/c2)
  
  // gg samples ZZ->llqq
  fSampleMass[160402] = 110;
  fSampleMass[160403] = 115;
  fSampleMass[160404] = 120;
  fSampleMass[160405] = 125;
  fSampleMass[160406] = 130;
  fSampleMass[160407] = 135;
  fSampleMass[160408] = 140;
  fSampleMass[160409] = 145;
  fSampleMass[160410] = 150;
  fSampleMass[160411] = 155;
  fSampleMass[160412] = 160;
  fSampleMass[160413] = 165;
  fSampleMass[160414] = 170;
  fSampleMass[160415] = 175;
  fSampleMass[160416] = 180;
  fSampleMass[160417] = 185;
  fSampleMass[160418] = 190;
  fSampleMass[160419] = 195;
  fSampleMass[160420] = 200;
  fSampleMass[160421] = 220;
  fSampleMass[160422] = 240;
  fSampleMass[160423] = 260;
  fSampleMass[160424] = 280;
  fSampleMass[160425] = 300;
  fSampleMass[160426] = 320;
  fSampleMass[160427] = 340;
  fSampleMass[160428] = 360;
  fSampleMass[160429] = 380;
  fSampleMass[160430] = 400;
  
  // VBF samples not included since the reweight is only for ggH samples
}


/*
  std::pair<double, double> HiggsllqqAnalysis::getCandidateTriggerSF(Int_t option)
  // returns event SF,SF_err using GoodMuon or GoodElectron depending on channel
  // option:  
  //        0 standard code 
  //        1 single code (presently as 0)
  //        2 dilepton code 
  //        3 combined code 
  {
  std::pair<double, double> result;
  result.first=1.;
  result.second=0.;
  
  if (GetIsMC()) {
    
    // find a fake run number representing the data period to which this MC event is somewhat associated
    Int_t representative_run_number(-1);
    Int_t period = GetPeriod();
    
    
    // divide in per2011B2_I and per2011J_M  
    // period with bath muon trigger L3-L4 in principle need to be introduced yet
    if (period == DataPeriod::BD) {
      representative_run_number = 186493;
    } else if (period == DataPeriod::EH) {
      representative_run_number = 186493;
    } else if (period == DataPeriod::I) {
      representative_run_number = 186493;
    } else if (period == DataPeriod::J) {
      representative_run_number = 191933;
    } else if (period == DataPeriod::K) {
      representative_run_number = 191933;
    } else if (period == DataPeriod::LM) {
      representative_run_number = 191933;
    }
    
    // use MeV (set to true to use GeV)
    //m_MuonTrigSF->setThresholds(false, representative_run_number);
    
    std::vector<TLorentzVector> muons;
    std::vector<TLorentzVector> electrons;
    
    if(Channel() == HiggsqqllAnalysis::MU2) {
      
      muons.push_back(CommonTools::getVector(GoodMuons.at(0)));
      muons.push_back(CommonTools::getVector(GoodMuons.at(1)));
      
    } else if(Channel() == HiggsqqllAnalysis::E2) {
      
      electrons.push_back(CommonTools::getVector(GoodElectrons.at(0)));
      electrons.push_back(CommonTools::getVector(GoodElectrons.at(1)));
      
    } 
    
#ifdef NEWTRIGSF
    result = fMuonTrigSF->GetTriggerSF(representative_run_number, false, muons, combined, electrons, tightpp, 0, option);
#else
    result = fMuonTrigSF->GetTriggerSF(representative_run_number, false, muons, combined, electrons, tightpp, 0);
#endif
    
  }
  
  return result; 
}

*/


 // TestSelection ntuple Methods
void HiggsllqqAnalysis::ResetAnalysisOutputBranches(analysis_output_struct *str) {
  str->runnumber             = -999;
  str->eventnumber           = -999;
  str->istagged              = -999;
  str->channel               = -999;
  str->isqcdevent            = -999;
  str->low_event             = -999;
  str->n_jets                = -999;
  str->n_b_jets              = -999;
  str->weight                = 1.00; //Careful, this is a weight, initialization = 1.00
  //
  str->lep1_m                = -999;
  str->lep1_pt               = -999;
  str->lep1_eta              = -999;
  str->lep1_phi              = -999;
  str->lep1_charge           = -999;
  str->lep1_caloiso          = -999;
  str->lep1_trackiso         = -999;
  str->lep1_quality          = -999;
  str->lep2_m                = -999;
  str->lep2_m                = -999;
  str->lep2_pt               = -999;
  str->lep2_eta              = -999;
  str->lep2_phi              = -999;
  str->lep2_charge           = -999;
  str->lep2_caloiso          = -999;
  str->lep2_trackiso         = -999;
  str->lep2_quality          = -999;
  str->lepZ_m                = -999;
  str->lepZ_pt               = -999;
  str->lepZ_eta              = -999;
  str->lepZ_phi              = -999;
  //
  str->Jet1_KF_index         = -999;
  str->Jet2_KF_index         = -999;
  str->Jet1_LJ_index         = -999;
  str->Jet2_LJ_index         = -999;
  str->Jet1_BP_index         = -999;
  str->Jet2_BP_index         = -999;
  //
  str->realJ1_KF_m           = -999;
  str->realJ1_KF_pt          = -999;
  str->realJ1_KF_eta         = -999;
  str->realJ1_KF_eta_det     = -999;
  str->realJ1_KF_phi         = -999;
  str->realJ1_KF_flavortruth = -999;
  str->realJ1_KF_MV1         = -999;
  str->realJ1_KF_jvf         = -999;
  str->realJ1_KF_ntrk        = -999;
  str->realJ1_KF_width       = -999;
  str->realJ1_KF_ntrk12      = -999;
  str->realJ1_KF_width12     = -999;
  str->realJ1_KF_pdg         = -999;
  str->realJ2_KF_m           = -999;
  str->realJ2_KF_pt          = -999;
  str->realJ2_KF_eta         = -999;
  str->realJ2_KF_eta_det     = -999;
  str->realJ2_KF_phi         = -999;
  str->realJ2_KF_flavortruth = -999;
  str->realJ2_KF_MV1         = -999;
  str->realJ2_KF_jvf         = -999;
  str->realJ2_KF_ntrk        = -999;
  str->realJ2_KF_width       = -999;
  str->realJ2_KF_ntrk12      = -999;
  str->realJ2_KF_width12     = -999;
  str->realJ2_KF_pdg         = -999;
  //
  str->realZ_KF_m            = -999;
  str->realZ_KF_pt           = -999;
  str->realZ_KF_eta          = -999;
  str->realZ_KF_phi          = -999;
  str->realH_KF_m            = -999;
  str->realH_KF_pt           = -999;
  str->realH_KF_eta          = -999;
  str->realH_KF_phi          = -999;
  //
  str->corrJ1_KF_m           = -999;
  str->corrJ1_KF_pt          = -999;
  str->corrJ1_KF_eta         = -999;
  str->corrJ1_KF_eta_det     = -999;
  str->corrJ1_KF_phi         = -999;
  str->corrJ1_KF_flavortruth = -999;
  str->corrJ2_KF_m           = -999;
  str->corrJ2_KF_pt          = -999;
  str->corrJ2_KF_eta         = -999;
  str->corrJ2_KF_eta_det     = -999;
  str->corrJ2_KF_phi         = -999;
  str->corrJ2_KF_flavortruth = -999;
  //
  str->corrZ_KF_m            = -999;
  str->corrZ_KF_pt           = -999;
  str->corrZ_KF_eta          = -999;
  str->corrZ_KF_phi          = -999;
  str->corrH_KF_m            = -999;
  str->corrH_KF_pt           = -999;
  str->corrH_KF_eta          = -999;
  str->corrH_KF_phi          = -999;
  //
  str->realJ1_LJ_m           = -999;
  str->realJ1_LJ_pt          = -999;
  str->realJ1_LJ_eta         = -999;
  str->realJ1_LJ_eta_det     = -999;
  str->realJ1_LJ_phi         = -999;
  str->realJ1_LJ_flavortruth = -999;
  str->realJ1_LJ_MV1         = -999;
  str->realJ1_LJ_jvf         = -999;
  str->realJ1_LJ_ntrk        = -999;
  str->realJ1_LJ_width       = -999;
  str->realJ1_LJ_ntrk12      = -999;
  str->realJ1_LJ_width12     = -999;
  str->realJ1_LJ_pdg         = -999;
  str->realJ2_LJ_m           = -999;
  str->realJ2_LJ_pt          = -999;
  str->realJ2_LJ_eta         = -999;
  str->realJ2_LJ_eta_det     = -999;
  str->realJ2_LJ_phi         = -999;
  str->realJ2_LJ_flavortruth = -999;
  str->realJ2_LJ_MV1         = -999;
  str->realJ2_LJ_jvf         = -999;
  str->realJ2_LJ_ntrk        = -999;
  str->realJ2_LJ_width       = -999;
  str->realJ2_LJ_ntrk12      = -999;
  str->realJ2_LJ_width12     = -999;
  str->realJ2_LJ_pdg         = -999;
  //
  str->realZ_LJ_m            = -999;
  str->realZ_LJ_pt           = -999;
  str->realZ_LJ_eta          = -999;
  str->realZ_LJ_phi          = -999;
  str->realH_LJ_m            = -999;
  str->realH_LJ_pt           = -999;
  str->realH_LJ_eta          = -999;
  str->realH_LJ_phi          = -999;
  ///////////////////
  str->realJ1_BP_m           = -999;
  str->realJ1_BP_pt          = -999;
  str->realJ1_BP_eta         = -999;
  str->realJ1_BP_eta_det     = -999;
  str->realJ1_BP_phi         = -999;
  str->realJ1_BP_flavortruth = -999;
  str->realJ1_BP_MV1         = -999;
  str->realJ1_BP_jvf         = -999;
  str->realJ1_BP_ntrk        = -999;
  str->realJ1_BP_width       = -999;
  str->realJ1_BP_ntrk12      = -999;
  str->realJ1_BP_width12     = -999;
  str->realJ1_BP_pdg         = -999;
  str->realJ2_BP_m           = -999;
  str->realJ2_BP_pt          = -999;
  str->realJ2_BP_eta         = -999;
  str->realJ2_BP_eta_det     = -999;
  str->realJ2_BP_phi         = -999;
  str->realJ2_BP_flavortruth = -999;
  str->realJ2_BP_MV1         = -999;
  str->realJ2_BP_jvf         = -999;
  str->realJ2_BP_ntrk        = -999;
  str->realJ2_BP_width       = -999;
  str->realJ2_BP_ntrk12      = -999;
  str->realJ2_BP_width12     = -999;
  str->realJ2_BP_pdg         = -999;
  //
  str->realZ_BP_m            = -999;
  str->realZ_BP_pt           = -999;
  str->realZ_BP_eta          = -999;
  str->realZ_BP_phi          = -999;
  str->realH_BP_m            = -999;
  str->realH_BP_pt           = -999;
  str->realH_BP_eta          = -999;
  str->realH_BP_phi          = -999;
  ///////////////////////////////////
  str->chisquare             = -999;
  str->met                   = -999;
  str->sumet                 = -999;
  str->NPV                   = -999;
  str->mu                    = -999;
  str->HFOR                  =   -1;
  str->Entries               =   -1;
  str->truthH_pt             =   -1;
  str->ggFweight             = 1.00; //Careful, these are weights, initialization = 1.00
  str->btagSF                = 1.00; //Careful, these are weights, initialization = 1.00
  str->trig_SF               = 1.00; //Careful, these are weights, initialization = 1.00
  str->trig_SF2              = 1.00; //Careful, these are weights, initialization = 1.00
  str->trig_SFC              = 1.00; //Careful, these are weights, initialization = 1.00
  str->trig_flag             =   -1;
  /////////////
  str->dPhi_ll               =   -1;
  str->dR_ll                 =   -1;
  //
  str->dPhi_KF_jj            =   -1;
  str->dPhi_LJ_jj            =   -1;
  str->dPhi_BP_jj            =   -1;
  str->dR_KF_jj              =   -1;
  str->dR_LJ_jj              =   -1;
  str->dR_BP_jj              =   -1;
  //
  str->dPhi_KF_ZZ            =   -1;
  str->dPhi_LJ_ZZ            =   -1;
  str->dPhi_BP_ZZ            =   -1;
  str->dR_KF_ZZ              =   -1;
  str->dR_LJ_ZZ              =   -1;
  str->dR_BP_ZZ              =   -1;
  /////////////
  str->total_jet_ntrk1       =   -1;
  str->total_jet_width1      =   -1;
  str->total_jet_ntrk2       =   -1;
  str->total_jet_width2      =   -1;
  str->total_jet_ntrk3       =   -1;
  str->total_jet_width3      =   -1;
  //
  if(FillGluon){
    str->realJ1_KF_Fisher    = -999;
    str->realJ2_KF_Fisher    = -999;
    str->ll_2_KF_jets        = -999;
    str->corrJ1_KF_Fisher    = -999;
    str->corrJ2_KF_Fisher    = -999;
    str->ll_2_KF_jets_corr   = -999;
    
    //Flavour Composition Variables    
    str->AllJet_MV1_1        = -999; // 1 good jet
    str->AllJet_MV1_2        = -999; // 2 good jet
    str->AllJet_MV1_3        = -999; // 3 good jet
    ///////////////////////
    str->xWin_44p_4var       = -1;
    str->yWin_44p_4var       = -1;
    str->zWin_44p_4var       = -1;
    str->gWin_44p_4var       = -1;
    str->xWin_44p_6var       = -1;
    str->yWin_44p_6var       = -1;
    str->zWin_44p_6var       = -1;
    str->gWin_44p_6var       = -1;
    str->xWin_44h_4var       = -1;
    str->yWin_44h_4var       = -1;
    str->zWin_44h_4var       = -1;
    str->gWin_44h_4var       = -1;
    str->xWin_44h_6var       = -1;
    str->yWin_44h_6var       = -1;
    str->zWin_44h_6var       = -1;
    str->gWin_44h_6var       = -1;
    str->xWin_64p_4var       = -1;
    str->yWin_64p_4var       = -1;
    str->zWin_64p_4var       = -1;
    str->gWin_64p_4var       = -1;
    str->xWin_64p_6var       = -1;
    str->yWin_64p_6var       = -1;
    str->zWin_64p_6var       = -1;
    str->gWin_64p_6var       = -1;
    str->xWin_64h_4var       = -1;
    str->yWin_64h_4var       = -1;
    str->zWin_64h_4var       = -1;
    str->gWin_64h_4var       = -1;
    str->xWin_64h_6var       = -1;
    str->yWin_64h_6var       = -1;
    str->zWin_64h_6var       = -1;
    str->gWin_64h_6var       = -1;
  }
}


void HiggsllqqAnalysis::SetAnalysisOutputBranches(analysis_output_struct *str) 
{
  analysistree->Branch("RunNumber",             &(str->runnumber));
  analysistree->Branch("EventNumber",           &(str->eventnumber));
  analysistree->Branch("isTagged",              &(str->istagged));
  analysistree->Branch("channel",               &(str->channel));
  analysistree->Branch("isqcdevent",            &(str->isqcdevent));
  analysistree->Branch("low_event",             &(str->low_event));
  analysistree->Branch("HFOR",                  &(str->HFOR));
  analysistree->Branch("Entries",               &(str->Entries));
  //
  analysistree->Branch("n_jets",                &(str->n_jets));
  analysistree->Branch("n_b_jets",              &(str->n_b_jets));
  analysistree->Branch("weight",                &(str->weight));
  //  
  analysistree->Branch("lep1_m",                &(str->lep1_m));
  analysistree->Branch("lep1_pt",               &(str->lep1_pt));
  analysistree->Branch("lep1_eta",              &(str->lep1_eta));
  analysistree->Branch("lep1_phi",              &(str->lep1_phi));
  analysistree->Branch("lep1_charge",           &(str->lep1_charge));
  analysistree->Branch("lep1_caloiso",          &(str->lep1_caloiso));
  analysistree->Branch("lep1_trackiso",         &(str->lep1_trackiso));
  analysistree->Branch("lep1_d0",               &(str->lep1_d0));
  analysistree->Branch("lep1_sigd0",            &(str->lep1_sigd0));
  analysistree->Branch("lep1_quality",          &(str->lep1_quality));
  analysistree->Branch("lep2_m",                &(str->lep2_m));
  analysistree->Branch("lep2_pt",               &(str->lep2_pt));
  analysistree->Branch("lep2_eta",              &(str->lep2_eta));
  analysistree->Branch("lep2_phi",              &(str->lep2_phi));
  analysistree->Branch("lep2_charge",           &(str->lep2_charge));
  analysistree->Branch("lep2_caloiso",          &(str->lep2_caloiso));
  analysistree->Branch("lep2_trackiso",         &(str->lep2_trackiso));
  analysistree->Branch("lep2_d0",               &(str->lep2_d0));
  analysistree->Branch("lep2_sigd0",            &(str->lep2_sigd0));
  analysistree->Branch("lep2_quality",          &(str->lep2_quality));
  //
  analysistree->Branch("lepZ_m",                &(str->lepZ_m));
  analysistree->Branch("lepZ_pt",               &(str->lepZ_pt));
  analysistree->Branch("lepZ_eta",              &(str->lepZ_eta));
  analysistree->Branch("lepZ_phi",              &(str->lepZ_phi));
  ///////////////
  analysistree->Branch("realJ1_KF_m",           &(str->realJ1_KF_m));
  analysistree->Branch("realJ1_KF_pt",          &(str->realJ1_KF_pt));
  analysistree->Branch("realJ1_KF_eta",         &(str->realJ1_KF_eta));
  analysistree->Branch("realJ1_KF_eta_det",     &(str->realJ1_KF_eta_det));
  analysistree->Branch("realJ1_KF_phi",         &(str->realJ1_KF_phi));
  analysistree->Branch("realJ1_KF_flavortruth", &(str->realJ1_KF_flavortruth));
  analysistree->Branch("realJ1_KF_MV1",         &(str->realJ1_KF_MV1));
  analysistree->Branch("realJ1_KF_jvf",         &(str->realJ1_KF_jvf));
  analysistree->Branch("realJ1_KF_ntrk",        &(str->realJ1_KF_ntrk));
  analysistree->Branch("realJ1_KF_width",       &(str->realJ1_KF_width));
  analysistree->Branch("realJ1_KF_ntrk12",      &(str->realJ1_KF_ntrk12));
  analysistree->Branch("realJ1_KF_width12",     &(str->realJ1_KF_width12));
  analysistree->Branch("realJ1_KF_pdg",         &(str->realJ1_KF_pdg));
  analysistree->Branch("realJ2_KF_m",           &(str->realJ2_KF_m));
  analysistree->Branch("realJ2_KF_pt",          &(str->realJ2_KF_pt));
  analysistree->Branch("realJ2_KF_eta",         &(str->realJ2_KF_eta));
  analysistree->Branch("realJ2_KF_eta_det",     &(str->realJ2_KF_eta_det));
  analysistree->Branch("realJ2_KF_phi",         &(str->realJ2_KF_phi));
  analysistree->Branch("realJ2_KF_flavortruth", &(str->realJ2_KF_flavortruth));
  analysistree->Branch("realJ2_KF_MV1",         &(str->realJ2_KF_MV1));
  analysistree->Branch("realJ2_KF_jvf",         &(str->realJ2_KF_jvf));
  analysistree->Branch("realJ2_KF_ntrk",        &(str->realJ2_KF_ntrk));
  analysistree->Branch("realJ2_KF_width",       &(str->realJ2_KF_width));
  analysistree->Branch("realJ2_KF_ntrk12",      &(str->realJ2_KF_ntrk12));
  analysistree->Branch("realJ2_KF_width12",     &(str->realJ2_KF_width12)); 
  analysistree->Branch("realJ2_KF_pdg",         &(str->realJ2_KF_pdg));
  //
  analysistree->Branch("realZ_KF_m",            &(str->realZ_KF_m));
  analysistree->Branch("realZ_KF_pt",           &(str->realZ_KF_pt));
  analysistree->Branch("realZ_KF_eta",          &(str->realZ_KF_eta));
  analysistree->Branch("realZ_KF_phi",          &(str->realZ_KF_phi));
  analysistree->Branch("realH_KF_m",            &(str->realH_KF_m));
  analysistree->Branch("realH_KF_pt",           &(str->realH_KF_pt));
  analysistree->Branch("realH_KF_eta",          &(str->realH_KF_eta));
  analysistree->Branch("realH_KF_phi",          &(str->realH_KF_phi));
  ///////////
  analysistree->Branch("corrJ1_KF_m",           &(str->corrJ1_KF_m));
  analysistree->Branch("corrJ1_KF_pt",          &(str->corrJ1_KF_pt));
  analysistree->Branch("corrJ1_KF_eta",         &(str->corrJ1_KF_eta));
  analysistree->Branch("corrJ1_KF_eta_det",     &(str->corrJ1_KF_eta_det));
  analysistree->Branch("corrJ1_KF_phi",         &(str->corrJ1_KF_phi));
  analysistree->Branch("corrJ1_KF_flavortruth", &(str->corrJ1_KF_flavortruth));
  analysistree->Branch("corrJ2_KF_m",           &(str->corrJ2_KF_m));
  analysistree->Branch("corrJ2_KF_pt",          &(str->corrJ2_KF_pt));
  analysistree->Branch("corrJ2_KF_eta",         &(str->corrJ2_KF_eta));
  analysistree->Branch("corrJ2_KF_eta_det",     &(str->corrJ2_KF_eta_det));
  analysistree->Branch("corrJ2_KF_phi",         &(str->corrJ2_KF_phi));
  analysistree->Branch("corrJ2_KF_flavortruth", &(str->corrJ2_KF_flavortruth));
  analysistree->Branch("corrZ_KF_m",            &(str->corrZ_KF_m));
  analysistree->Branch("corrZ_KF_pt",           &(str->corrZ_KF_pt));
  analysistree->Branch("corrZ_KF_eta",          &(str->corrZ_KF_eta));
  analysistree->Branch("corrZ_KF_phi",          &(str->corrZ_KF_phi));
  analysistree->Branch("corrH_KF_m",            &(str->corrH_KF_m));
  analysistree->Branch("corrH_KF_pt",           &(str->corrH_KF_pt));
  analysistree->Branch("corrH_KF_eta",          &(str->corrH_KF_eta));
  analysistree->Branch("corrH_KF_phi",          &(str->corrH_KF_phi));
  ////////////////////
  analysistree->Branch("realJ1_LJ_m",           &(str->realJ1_LJ_m));
  analysistree->Branch("realJ1_LJ_pt",          &(str->realJ1_LJ_pt));
  analysistree->Branch("realJ1_LJ_eta",         &(str->realJ1_LJ_eta));
  analysistree->Branch("realJ1_LJ_eta_det",     &(str->realJ1_LJ_eta_det));
  analysistree->Branch("realJ1_LJ_phi",         &(str->realJ1_LJ_phi));
  analysistree->Branch("realJ1_LJ_flavortruth", &(str->realJ1_LJ_flavortruth));
  analysistree->Branch("realJ1_LJ_MV1",         &(str->realJ1_LJ_MV1));
  analysistree->Branch("realJ1_LJ_jvf",         &(str->realJ1_LJ_jvf));
  analysistree->Branch("realJ1_LJ_ntrk",        &(str->realJ1_LJ_ntrk));
  analysistree->Branch("realJ1_LJ_width",       &(str->realJ1_LJ_width));
  analysistree->Branch("realJ1_LJ_ntrk12",      &(str->realJ1_LJ_ntrk12));
  analysistree->Branch("realJ1_LJ_width12",     &(str->realJ1_LJ_width12));
  analysistree->Branch("realJ1_LJ_pdg",         &(str->realJ1_LJ_pdg));
  analysistree->Branch("realJ2_LJ_m",           &(str->realJ2_LJ_m));
  analysistree->Branch("realJ2_LJ_pt",          &(str->realJ2_LJ_pt));
  analysistree->Branch("realJ2_LJ_eta",         &(str->realJ2_LJ_eta));
  analysistree->Branch("realJ2_LJ_eta_det",     &(str->realJ2_LJ_eta_det));
  analysistree->Branch("realJ2_LJ_phi",         &(str->realJ2_LJ_phi));
  analysistree->Branch("realJ2_LJ_flavortruth", &(str->realJ2_LJ_flavortruth));
  analysistree->Branch("realJ2_LJ_MV1",         &(str->realJ2_LJ_MV1));
  analysistree->Branch("realJ2_LJ_jvf",         &(str->realJ2_LJ_jvf));
  analysistree->Branch("realJ2_LJ_ntrk",        &(str->realJ2_LJ_ntrk));
  analysistree->Branch("realJ2_LJ_width",       &(str->realJ2_LJ_width));
  analysistree->Branch("realJ2_LJ_ntrk12",      &(str->realJ2_LJ_ntrk12));
  analysistree->Branch("realJ2_LJ_width12",     &(str->realJ2_LJ_width12)); 
  analysistree->Branch("realJ2_LJ_pdg",         &(str->realJ2_LJ_pdg));
  //
  analysistree->Branch("realZ_LJ_m",            &(str->realZ_LJ_m));
  analysistree->Branch("realZ_LJ_pt",           &(str->realZ_LJ_pt));
  analysistree->Branch("realZ_LJ_eta",          &(str->realZ_LJ_eta));
  analysistree->Branch("realZ_LJ_phi",          &(str->realZ_LJ_phi));
  analysistree->Branch("realH_LJ_m",            &(str->realH_LJ_m));
  analysistree->Branch("realH_LJ_pt",           &(str->realH_LJ_pt));
  analysistree->Branch("realH_LJ_eta",          &(str->realH_LJ_eta));
  analysistree->Branch("realH_LJ_phi",          &(str->realH_LJ_phi));
  ///////////
  analysistree->Branch("realJ1_BP_m",           &(str->realJ1_BP_m));
  analysistree->Branch("realJ1_BP_pt",          &(str->realJ1_BP_pt));
  analysistree->Branch("realJ1_BP_eta",         &(str->realJ1_BP_eta));
  analysistree->Branch("realJ1_BP_eta_det",     &(str->realJ1_BP_eta_det));
  analysistree->Branch("realJ1_BP_phi",         &(str->realJ1_BP_phi));
  analysistree->Branch("realJ1_BP_flavortruth", &(str->realJ1_BP_flavortruth));
  analysistree->Branch("realJ1_BP_MV1",         &(str->realJ1_BP_MV1));
  analysistree->Branch("realJ1_BP_jvf",         &(str->realJ1_BP_jvf));
  analysistree->Branch("realJ1_BP_ntrk",        &(str->realJ1_BP_ntrk));
  analysistree->Branch("realJ1_BP_width",       &(str->realJ1_BP_width));
  analysistree->Branch("realJ1_BP_ntrk12",      &(str->realJ1_BP_ntrk12));
  analysistree->Branch("realJ1_BP_width12",     &(str->realJ1_BP_width12));
  analysistree->Branch("realJ1_BP_pdg",         &(str->realJ1_BP_pdg));
  analysistree->Branch("realJ2_BP_m",           &(str->realJ2_BP_m));
  analysistree->Branch("realJ2_BP_pt",          &(str->realJ2_BP_pt));
  analysistree->Branch("realJ2_BP_eta",         &(str->realJ2_BP_eta));
  analysistree->Branch("realJ2_BP_eta_det",     &(str->realJ2_BP_eta_det));
  analysistree->Branch("realJ2_BP_phi",         &(str->realJ2_BP_phi));
  analysistree->Branch("realJ2_BP_flavortruth", &(str->realJ2_BP_flavortruth));
  analysistree->Branch("realJ2_BP_MV1",         &(str->realJ2_BP_MV1));
  analysistree->Branch("realJ2_BP_jvf",         &(str->realJ2_BP_jvf));
  analysistree->Branch("realJ2_BP_ntrk",        &(str->realJ2_BP_ntrk));
  analysistree->Branch("realJ2_BP_width",       &(str->realJ2_BP_width));
  analysistree->Branch("realJ2_BP_ntrk12",      &(str->realJ2_BP_ntrk12));
  analysistree->Branch("realJ2_BP_width12",     &(str->realJ2_BP_width12)); 
  analysistree->Branch("realJ2_BP_pdg",         &(str->realJ2_BP_pdg));
  //
  analysistree->Branch("realZ_BP_m",            &(str->realZ_BP_m));
  analysistree->Branch("realZ_BP_pt",           &(str->realZ_BP_pt));
  analysistree->Branch("realZ_BP_eta",          &(str->realZ_BP_eta));
  analysistree->Branch("realZ_BP_phi",          &(str->realZ_BP_phi));
  analysistree->Branch("realH_BP_m",            &(str->realH_BP_m));
  analysistree->Branch("realH_BP_pt",           &(str->realH_BP_pt));
  analysistree->Branch("realH_BP_eta",          &(str->realH_BP_eta));
  analysistree->Branch("realH_BP_phi",          &(str->realH_BP_phi));
  ////////////////
  analysistree->Branch("Jet1_KF_index",         &(str->Jet1_KF_index));
  analysistree->Branch("Jet2_KF_index",         &(str->Jet2_KF_index));
  analysistree->Branch("Jet1_LJ_index",         &(str->Jet1_LJ_index));
  analysistree->Branch("Jet2_LJ_index",         &(str->Jet2_LJ_index));
  analysistree->Branch("Jet1_BP_index",         &(str->Jet1_BP_index));
  analysistree->Branch("Jet2_BP_index",         &(str->Jet2_BP_index));
  //
  analysistree->Branch("chisquare",             &(str->chisquare));
  //////////////
  analysistree->Branch("met",                   &(str->met));
  analysistree->Branch("sumet",                 &(str->sumet));
  analysistree->Branch("btagSF",                &(str->btagSF));
  analysistree->Branch("NPV",                   &(str->NPV));
  analysistree->Branch("truthH_pt",             &(str->truthH_pt));
  analysistree->Branch("ggFweight",             &(str->ggFweight));
  analysistree->Branch("mu",                    &(str->mu));
  analysistree->Branch("trig_SF",               &(str->trig_SF));
  analysistree->Branch("trig_SF2",              &(str->trig_SF2));
  analysistree->Branch("trig_SFC",              &(str->trig_SFC));
  analysistree->Branch("trig_flag",             &(str->trig_flag));
  //
  analysistree->Branch("dPhi_ll",               &(str->dPhi_ll));
  analysistree->Branch("dR_ll",                 &(str->dR_ll));
  //
  analysistree->Branch("dPhi_KF_jj",            &(str->dPhi_KF_jj));
  analysistree->Branch("dPhi_LJ_jj",            &(str->dPhi_LJ_jj));
  analysistree->Branch("dPhi_BP_jj",            &(str->dPhi_BP_jj));
  analysistree->Branch("dR_KF_jj",              &(str->dR_KF_jj));
  analysistree->Branch("dR_LJ_jj",              &(str->dR_LJ_jj));
  analysistree->Branch("dR_BP_jj",              &(str->dR_BP_jj));
  //
  analysistree->Branch("dPhi_KF_ZZ",            &(str->dPhi_KF_ZZ));
  analysistree->Branch("dPhi_LJ_ZZ",            &(str->dPhi_LJ_ZZ));
  analysistree->Branch("dPhi_BP_ZZ",            &(str->dPhi_BP_ZZ));
  analysistree->Branch("dR_KF_ZZ",              &(str->dR_KF_ZZ));
  analysistree->Branch("dR_LJ_ZZ",              &(str->dR_LJ_ZZ));
  analysistree->Branch("dR_BP_ZZ",              &(str->dR_BP_ZZ));
  //  
  analysistree->Branch("total_jet_ntrk1",       &(str->total_jet_ntrk1));
  analysistree->Branch("total_jet_width1",      &(str->total_jet_width1));
  analysistree->Branch("total_jet_ntrk2",       &(str->total_jet_ntrk2));
  analysistree->Branch("total_jet_width2",      &(str->total_jet_width2));
  analysistree->Branch("total_jet_ntrk3",       &(str->total_jet_ntrk3));
  analysistree->Branch("total_jet_width3",      &(str->total_jet_width3));
  
  if(FillGluon)
    {
      analysistree->Branch("realJ1_KF_Fisher",  &(str->realJ1_KF_Fisher));
      analysistree->Branch("realJ2_KF_Fisher",  &(str->realJ2_KF_Fisher));
      analysistree->Branch("ll_2_KF_jets",      &(str->ll_2_KF_jets));
      analysistree->Branch("corrJ1_KF_Fisher",  &(str->corrJ1_KF_Fisher));
      analysistree->Branch("corrJ2_KF_Fisher",  &(str->corrJ2_KF_Fisher));
      analysistree->Branch("ll_2_KF_jets_corr", &(str->ll_2_KF_jets_corr));
      
      //Flavour Composition Variables
      analysistree->Branch("AllJet_MV1_1",      &(str->AllJet_MV1_1));    // 1 good jet
      analysistree->Branch("AllJet_MV1_2",      &(str->AllJet_MV1_2));    // 2 good jet
      analysistree->Branch("AllJet_MV1_3",      &(str->AllJet_MV1_3));    // 3 good jet
      
      analysistree->Branch("xWin_44p_4var",     &(str->xWin_44p_4var));
      analysistree->Branch("yWin_44p_4var",     &(str->yWin_44p_4var));
      analysistree->Branch("zWin_44p_4var",     &(str->zWin_44p_4var));
      analysistree->Branch("gWin_44p_4var",     &(str->gWin_44p_4var));
      analysistree->Branch("xWin_44p_6var",     &(str->xWin_44p_6var));
      analysistree->Branch("yWin_44p_6var",     &(str->yWin_44p_6var));
      analysistree->Branch("zWin_44p_6var",     &(str->zWin_44p_6var));
      analysistree->Branch("gWin_44p_6var",     &(str->gWin_44p_6var));
      analysistree->Branch("xWin_44h_4var",     &(str->xWin_44h_4var));
      analysistree->Branch("yWin_44h_4var",     &(str->yWin_44h_4var));
      analysistree->Branch("zWin_44h_4var",     &(str->zWin_44h_4var));
      analysistree->Branch("gWin_44h_4var",     &(str->gWin_44h_4var));
      analysistree->Branch("xWin_44h_6var",     &(str->xWin_44h_6var));
      analysistree->Branch("yWin_44h_6var",     &(str->yWin_44h_6var));
      analysistree->Branch("zWin_44h_6var",     &(str->zWin_44h_6var));
      analysistree->Branch("gWin_44h_6var",     &(str->gWin_44h_6var));
      analysistree->Branch("xWin_64p_4var",     &(str->xWin_64p_4var));
      analysistree->Branch("yWin_64p_4var",     &(str->yWin_64p_4var));
      analysistree->Branch("zWin_64p_4var",     &(str->zWin_64p_4var));
      analysistree->Branch("gWin_64p_4var",     &(str->gWin_64p_4var));
      analysistree->Branch("xWin_64p_6var",     &(str->xWin_64p_6var));
      analysistree->Branch("yWin_64p_6var",     &(str->yWin_64p_6var));
      analysistree->Branch("zWin_64p_6var",     &(str->zWin_64p_6var));
      analysistree->Branch("gWin_64p_6var",     &(str->gWin_64p_6var));
      analysistree->Branch("xWin_64h_4var",     &(str->xWin_64h_4var));
      analysistree->Branch("yWin_64h_4var",     &(str->yWin_64h_4var));
      analysistree->Branch("zWin_64h_4var",     &(str->zWin_64h_4var));
      analysistree->Branch("gWin_64h_4var",     &(str->gWin_64h_4var));
      analysistree->Branch("xWin_64h_6var",     &(str->xWin_64h_6var));
      analysistree->Branch("yWin_64h_6var",     &(str->yWin_64h_6var));
      analysistree->Branch("zWin_64h_6var",     &(str->zWin_64h_6var));
      analysistree->Branch("gWin_64h_6var",     &(str->gWin_64h_6var));
    }
}


void HiggsllqqAnalysis::FillAnalysisOutputTree(analysis_output_struct *str, Int_t cut, UInt_t channel)
{  
  Int_t minimum_cut(0), second_cut(0); 
  
  if(GetDoQCDSelection()) 
    {
      minimum_cut = HllqqCutFlow::PtLeptons;//DiJetMass; //DileptonMass; //NumTagJets; //MET; //NumberOfLeptons;
      second_cut  = HllqqCutFlow::DiJetMass;
    }
  else
    { 
      minimum_cut = HllqqCutFlow::PtLeptons;
      second_cut  = HllqqCutFlow::DiJetMass;    
    }
  
  
  if(cut >= minimum_cut)
    {
      str->n_jets   = m_GoodJets.size();
      
      std::pair<float,float> jetsf;
      jetsf.first  = -1.;
      jetsf.second = -1.;
      
      TLorentzVector lep1, lep2, lepZ;
      
      
      if (!isMC()) str->runnumber   = ntuple->eventinfo.RunNumber();
      if (isMC())  str->runnumber   = ntuple->eventinfo.mc_channel_number();
      
      str->eventnumber = ntuple->eventinfo.EventNumber();
      str->HFOR        = HFOR_value;
      str->Entries     = fChain->GetEntries();      
      str->channel     = channel;
      str->low_event   = GetDoLowMass()      ? 1 : 0;
      str->isqcdevent  = GetDoQCDSelection() ? 1 : 0;
      str->met         = getCorrectMETValue();
      str->sumet       = ntuple->MET_RefFinal.sumet();
      str->NPV         = getNumberOfGoodVertices();
      str->truthH_pt   = getTruthHiggsPt();
      str->ggFweight   = getggFWeight();
      str->mu          = (isMC() && ntuple->eventinfo.lbn()==1 && int(ntuple->eventinfo.averageIntPerXing()+0.5)==1) ? 0. : ntuple->eventinfo.averageIntPerXing();
      str->weight      = 1.;
      if(isMC()) 
	str->weight   = 1.*getSFWeight()*getPileupWeight()*ntuple->eventinfo.mc_event_weight();
      
      
      //  Reset and fill trigger flag word
      //  str->trig_SF   = trig_SF;
      //  str->trig_SF2  = trig_SF2;
      //  str->trig_SFC  = trig_SFC;
      
      
      if(channel == HiggsllqqAnalysis::MU2)
	{
	  D3PDReader::MuonD3PDObjectElement *mu_1 = m_GoodMuons.at(0)->GetMuon();
	  D3PDReader::MuonD3PDObjectElement *mu_2 = m_GoodMuons.at(1)->GetMuon();
	  
	  //Fill flag to distinguish single & double lepton trigger
	  if (passesSingleMuonTrigger() && passesDiMuonTrigger()) str->trig_flag = 3;
	  else if (passesDiMuonTrigger())                         str->trig_flag = 2;
	  else if (passesSingleMuonTrigger())                     str->trig_flag = 1;
	  
	  
	  if(m_GoodMuons.at(0)->Get4Momentum()->Pt() >= m_GoodMuons.at(1)->Get4Momentum()->Pt())
	    {
	      //Filling leading Muon
	      str->lep1_m         = mu_pdg_mass;
	      str->lep1_pt        = mu_1->pt();
	      str->lep1_eta       = mu_1->eta();
	      str->lep1_phi       = mu_1->phi();
	      str->lep1_charge    = mu_1->charge();
	      str->lep1_caloiso   = mu_1->etcone20()/mu_1->pt();
	      str->lep1_trackiso  = mu_1->ptcone20()/mu_1->pt();
	      str->lep1_d0        = m_GoodMuons.at(0)->d0();
	      str->lep1_sigd0     = m_GoodMuons.at(0)->d0_sig();
	      
	      if (m_GoodMuons.at(0)->family() == Muon::MUID)
		str->lep1_quality = 1;
	      else if (m_GoodMuons.at(0)->family() == Muon::STACO)
		str->lep1_quality = 2;
	      else if (m_GoodMuons.at(0)->family() == Muon::CALO)
		str->lep1_quality = 3;
	      
	      
	      //Filling the Second Muon
	      str->lep2_m         = mu_pdg_mass;
	      str->lep2_pt        = mu_2->pt();
	      str->lep2_eta       = mu_2->eta();
	      str->lep2_phi       = mu_2->phi();
	      str->lep2_charge    = mu_2->charge();
	      str->lep2_caloiso   = mu_2->etcone20()/mu_2->pt();
	      str->lep2_trackiso  = mu_2->ptcone20()/mu_2->pt();
	      str->lep2_d0        = m_GoodMuons.at(1)->d0();
	      str->lep2_sigd0     = m_GoodMuons.at(1)->d0_sig();
	      
	      if (m_GoodMuons.at(1)->family() == Muon::MUID)
		str->lep2_quality = 1;
	      else if (m_GoodMuons.at(1)->family() == Muon::STACO)
		str->lep2_quality = 2;
	      else if (m_GoodMuons.at(1)->family() == Muon::CALO)
		str->lep2_quality = 3;
	    }
	  else
	    {
	      //Filling leading Muon
	      str->lep1_m         = mu_pdg_mass;
	      str->lep1_pt        = mu_2->pt();
	      str->lep1_eta       = mu_2->eta();
	      str->lep1_phi       = mu_2->phi();
	      str->lep1_charge    = mu_2->charge();
	      str->lep1_caloiso   = mu_2->etcone20()/mu_2->pt();
	      str->lep1_trackiso  = mu_2->ptcone20()/mu_2->pt();
	      str->lep1_d0        = m_GoodMuons.at(1)->d0();
	      str->lep1_sigd0     = m_GoodMuons.at(1)->d0_sig();
	      
	      if (m_GoodMuons.at(1)->family() == Muon::MUID)
		str->lep1_quality = 1;
	      else if (m_GoodMuons.at(1)->family() == Muon::STACO)
		str->lep1_quality = 2;
	      else if (m_GoodMuons.at(1)->family() == Muon::CALO)
		str->lep1_quality = 3;
	      
	      
	      //Filling the Second Muon
	      str->lep2_m         = mu_pdg_mass;
	      str->lep2_pt        = mu_1->pt();
	      str->lep2_eta       = mu_1->eta();
	      str->lep2_phi       = mu_1->phi();
	      str->lep2_charge    = mu_1->charge();
	      str->lep2_caloiso   = mu_1->etcone20()/mu_1->pt();
	      str->lep2_trackiso  = mu_1->ptcone20()/mu_1->pt();
	      str->lep2_d0        = m_GoodMuons.at(0)->d0();
	      str->lep2_sigd0     = m_GoodMuons.at(0)->d0_sig();
	      
	      if (m_GoodMuons.at(0)->family() == Muon::MUID)
		str->lep2_quality = 1;
	      else if (m_GoodMuons.at(0)->family() == Muon::STACO)
		str->lep2_quality = 2;
	      else if (m_GoodMuons.at(0)->family() == Muon::CALO)
		str->lep2_quality = 3;	
	    }
	  
	  lep1.SetPtEtaPhiM(str->lep1_pt,str->lep1_eta,str->lep1_phi,str->lep1_m);
	  lep2.SetPtEtaPhiM(str->lep2_pt,str->lep2_eta,str->lep2_phi,str->lep2_m);
	  lepZ = lep1 + lep2;
	  
	  str->lepZ_m           = lepZ.M();
	  str->lepZ_pt          = lepZ.Pt();
	  str->lepZ_eta         = lepZ.Eta();
	  str->lepZ_phi         = lepZ.Phi();
	  str->dPhi_ll          = TMath::Abs(lep1.DeltaPhi(lep2));
	  str->dR_ll            = lep1.DeltaR(lep2);
	  
	} //End of MU2 Channel if
      
      
      if(channel == HiggsllqqAnalysis::E2)
	{	
	  D3PDReader::ElectronD3PDObjectElement *el_1 = m_GoodElectrons.at(0)->GetElectron();
	  D3PDReader::ElectronD3PDObjectElement *el_2 = m_GoodElectrons.at(1)->GetElectron();
	  
	  
	  //Fill flag to distinguish single & double lepton trigger
	  if (passesSingleElectronTrigger() && passesDiElectronTrigger()) str->trig_flag = 3;
	  else if (passesDiElectronTrigger())                             str->trig_flag = 2;
	  else if (passesSingleElectronTrigger())                         str->trig_flag = 1;
	  
	  if(m_GoodElectrons.at(0)->Get4Momentum()->Pt() >= m_GoodElectrons.at(1)->Get4Momentum()->Pt())
	    {    
	      //Filling leading Electron
	      str->lep1_m         = el_pdg_mass;
	      str->lep1_pt        = el_1->pt();
	      str->lep1_eta       = el_1->eta();
	      str->lep1_phi       = el_1->phi();
	      str->lep1_charge    = el_1->charge();
	      str->lep1_caloiso   = el_1->Etcone20()/el_1->pt();
	      str->lep1_trackiso  = el_1->ptcone20()/el_1->pt();
	      str->lep1_d0        = m_GoodElectrons.at(0)->d0();
	      str->lep1_sigd0     = m_GoodElectrons.at(0)->d0_sig();
	      
	      if(el_1->tightPP() == 1)
		str->lep1_quality = 3;	  
	      else if(el_1->mediumPP() == 1)
		str->lep1_quality = 2;
	      else if(el_1->loosePP() == 1)
		str->lep1_quality = 1;
	      
	      
	      //Filling the Second Electron
	      str->lep2_m         = el_pdg_mass;
	      str->lep2_pt        = el_2->pt();
	      str->lep2_eta       = el_2->eta();
	      str->lep2_phi       = el_2->phi();
	      str->lep2_charge    = el_2->charge();
	      str->lep2_caloiso   = el_2->Etcone20()/el_2->pt();
	      str->lep2_trackiso  = el_2->ptcone20()/el_2->pt();
	      str->lep2_d0        = m_GoodElectrons.at(1)->d0();
	      str->lep2_sigd0     = m_GoodElectrons.at(1)->d0_sig();
	      
	      if(el_2->tightPP() == 1)
		str->lep2_quality = 3;	  
	      else if(el_2->mediumPP() == 1)
		str->lep2_quality = 2;
	      else if(el_2->loosePP() == 1)
		str->lep2_quality = 1;
	      
	    }
	  else
	    {	    
	      //Filling leading Electron
	      str->lep1_m         = el_pdg_mass;
	      str->lep1_pt        = el_2->pt();
	      str->lep1_eta       = el_2->eta();
	      str->lep1_phi       = el_2->phi();
	      str->lep1_charge    = el_2->charge();
	      str->lep1_caloiso   = el_2->Etcone20()/el_2->pt();
	      str->lep1_trackiso  = el_2->ptcone20()/el_2->pt();
	      str->lep1_d0        = m_GoodElectrons.at(1)->d0();
	      str->lep1_sigd0     = m_GoodElectrons.at(1)->d0_sig();
	      
	      if(el_2->tightPP() == 1)
		str->lep1_quality = 3;	  
	      else if(el_2->mediumPP() == 1)
		str->lep1_quality = 2;
	      else if(el_2->loosePP() == 1)
		str->lep1_quality = 1;
	      
	      
	      //Filling the Second Electron
	      str->lep2_m         = el_pdg_mass;
	      str->lep2_pt        = el_1->pt();
	      str->lep2_eta       = el_1->eta();
	      str->lep2_phi       = el_1->phi();
	      str->lep2_charge    = el_1->charge();
	      str->lep2_caloiso   = el_1->Etcone20()/el_1->pt();
	      str->lep2_trackiso  = el_1->ptcone20()/el_1->pt();
	      str->lep2_d0        = m_GoodElectrons.at(0)->d0();
	      str->lep2_sigd0     = m_GoodElectrons.at(0)->d0_sig();
	      
	      if(el_1->tightPP() == 1)
		str->lep2_quality = 3;	  
	      else if(el_1->mediumPP() == 1)
		str->lep2_quality = 2;
	      else if(el_1->loosePP() == 1)
		str->lep2_quality = 1;
	    }
	  
	  lep1.SetPtEtaPhiM(str->lep1_pt,str->lep1_eta,str->lep1_phi,str->lep1_m);
	  lep2.SetPtEtaPhiM(str->lep2_pt,str->lep2_eta,str->lep2_phi,str->lep2_m);
	  lepZ = lep1 + lep2;
	  
	  str->lepZ_m           = lepZ.M();
	  str->lepZ_pt          = lepZ.Pt();
	  str->lepZ_eta         = lepZ.Eta();
	  str->lepZ_phi         = lepZ.Phi();
	  str->dPhi_ll          = TMath::Abs(lep1.DeltaPhi(lep2));
	  str->dR_ll            = lep1.DeltaR(lep2);
	  
	} //End of E2 Channel if
      
      

      
      ///////////////////    JET FILLING!    ////////////////////      
      
      float tmpbtagsf = 1.;
      int howmanytags = GetNumOfTags();
      str->n_b_jets   = howmanytags;
      
      if(howmanytags==2)
	str->istagged = 1;
      if(howmanytags<2)
	str->istagged = 0;
      if(howmanytags>2)
	str->istagged = -1;
      
      // Using the Global Jets index variables to fill the reduced TestSelection ntuple,Warning: be sure that you call at least DiJetMass Cut if is OUT of the above "if"!
      
      if (cut >= second_cut /*m_GoodJets.size()>=2*/) 
	{  
	  std::pair<int,int> SelectedJets;
	  SelectedJets.first  = Pair_jet1;
	  SelectedJets.second = Pair_jet2;
	  JetBestPairResult();
	  
	  
	  // Getting the jet pair (KF or LJ or BP)
	  D3PDReader::JetD3PDObjectElement *Jet_1_KF = m_GoodJets.at(SelectedJets.first)->GetJet();
	  D3PDReader::JetD3PDObjectElement *Jet_2_KF = m_GoodJets.at(SelectedJets.second)->GetJet();
	  D3PDReader::JetD3PDObjectElement *Jet_1_LJ = m_GoodJets.at(0)->GetJet();
	  D3PDReader::JetD3PDObjectElement *Jet_2_LJ = m_GoodJets.at(1)->GetJet();
	  D3PDReader::JetD3PDObjectElement *Jet_1_BP = m_GoodJets.at(Jone)->GetJet();
	  D3PDReader::JetD3PDObjectElement *Jet_2_BP = m_GoodJets.at(Jtwo)->GetJet();
	  
	  
	  if (Look_b_SFs) // Looping to find the btagging SFs into the GoodJets selected into the event
	    {
	      for (UInt_t w = 0; w < m_GoodJets.size(); w++)
		{	
		  pair<double,double> thisjetSF = GetJetSFsvalue(w);
		  tmpbtagsf *= thisjetSF.first;
		  
		  if(w==0)
		    {
		      D3PDReader::JetD3PDObjectElement *Jet0 = m_GoodJets.at(0)->GetJet();
		      str->total_jet_ntrk1  = Jet0->nTrk();
		      str->total_jet_width1 = Jet0->WIDTH();
		      str->AllJet_MV1_1     = GetMV1value(m_GoodJets.at(0));
		    }
		  if(w==1)
		    {
		      D3PDReader::JetD3PDObjectElement *Jet1 = m_GoodJets.at(1)->GetJet();
		      str->total_jet_ntrk2  = Jet1->nTrk();
		      str->total_jet_width2 = Jet1->WIDTH();
		      str->AllJet_MV1_2     = GetMV1value(m_GoodJets.at(1));
		    }
		  if(w==2)
		    {
		      D3PDReader::JetD3PDObjectElement *Jet2 = m_GoodJets.at(2)->GetJet();
		      str->total_jet_ntrk3  = Jet2->nTrk();
		      str->total_jet_width3 = Jet2->WIDTH();
		      str->AllJet_MV1_3     = GetMV1value(m_GoodJets.at(2));
		    }
		}
	    }
	  
	  
	  //Filling the new definitions of Tracks and Width
	  std::pair<Float_t, Float_t> InfoNtracksWidthJ1_KF = InfoTracks(SelectedJets.first);
	  std::pair<Float_t, Float_t> InfoNtracksWidthJ2_KF = InfoTracks(SelectedJets.second);
	  std::pair<Float_t, Float_t> InfoNtracksWidthJ1_LJ = InfoTracks(0);
	  std::pair<Float_t, Float_t> InfoNtracksWidthJ2_LJ = InfoTracks(1);
	  std::pair<Float_t, Float_t> InfoNtracksWidthJ1_BP = InfoTracks(Jone);
	  std::pair<Float_t, Float_t> InfoNtracksWidthJ2_BP = InfoTracks(Jtwo);
	  
	  str->Jet1_KF_index         = SelectedJets.first;
	  str->Jet2_KF_index         = SelectedJets.second;
	  str->Jet1_LJ_index         = 0;//obvious!
	  str->Jet2_LJ_index         = 1;//obvious!
	  str->Jet1_BP_index         = Jone;
	  str->Jet2_BP_index         = Jtwo;
	  
	  
	  str->realJ1_KF_m           = m_GoodJets.at(SelectedJets.first)->Get4Momentum()->M();
	  str->realJ1_KF_pt          = m_GoodJets.at(SelectedJets.first)->rightpt();
	  str->realJ1_KF_eta         = m_GoodJets.at(SelectedJets.first)->righteta();
	  str->realJ1_KF_phi         = m_GoodJets.at(SelectedJets.first)->rightphi();
	  str->realJ1_KF_eta_det     = Jet_1_KF->emscale_eta();
	  str->realJ1_KF_flavortruth = Jet_1_KF->flavor_truth_label();
	  str->realJ1_KF_jvf         = Jet_1_KF->jvtxf();
	  str->realJ1_KF_ntrk        = Jet_1_KF->nTrk();
	  str->realJ1_KF_width       = Jet_1_KF->WIDTH();
	  str->realJ1_KF_MV1         = GetMV1value(m_GoodJets.at(SelectedJets.first));
	  str->realJ1_KF_ntrk12      = InfoNtracksWidthJ1_KF.first;
	  str->realJ1_KF_width12     = InfoNtracksWidthJ1_KF.second;
	  
	  str->realJ2_KF_m           = m_GoodJets.at(SelectedJets.second)->Get4Momentum()->M();
	  str->realJ2_KF_pt          = m_GoodJets.at(SelectedJets.second)->rightpt();
	  str->realJ2_KF_eta         = m_GoodJets.at(SelectedJets.second)->righteta();
	  str->realJ2_KF_phi         = m_GoodJets.at(SelectedJets.second)->rightphi();
	  str->realJ2_KF_eta_det     = Jet_2_KF->emscale_eta();
	  str->realJ2_KF_flavortruth = Jet_2_KF->flavor_truth_label();
	  str->realJ2_KF_jvf         = Jet_2_KF->jvtxf();
	  str->realJ2_KF_ntrk        = Jet_2_KF->nTrk();
	  str->realJ2_KF_width       = Jet_2_KF->WIDTH();
	  str->realJ2_KF_MV1         = GetMV1value(m_GoodJets.at(SelectedJets.second));
	  str->realJ2_KF_ntrk12      = InfoNtracksWidthJ2_KF.first;
	  str->realJ2_KF_width12     = InfoNtracksWidthJ2_KF.second;
	  
	  
	  
	  TLorentzVector j1;
	  j1.SetPtEtaPhiM(str->realJ1_KF_pt,str->realJ1_KF_eta,str->realJ1_KF_phi,str->realJ1_KF_m);
	  TLorentzVector j2;
	  j2.SetPtEtaPhiM(str->realJ2_KF_pt,str->realJ2_KF_eta,str->realJ2_KF_phi,str->realJ2_KF_m);
	  TLorentzVector hadZ = j1   +   j2;
	  TLorentzVector H    = lepZ + hadZ;
	  
	  
	  str->realZ_KF_m     = hadZ.M();
	  str->realZ_KF_pt    = hadZ.Pt();
	  str->realZ_KF_eta   = hadZ.Eta();
	  str->realZ_KF_phi   = hadZ.Phi();
	  
	  str->realH_KF_m     = H.M();
	  str->realH_KF_pt    = H.Pt();
	  str->realH_KF_eta   = H.Eta();
	  str->realH_KF_phi   = H.Phi();
	  
	  
	  
	  //ANGULAR JET VARIABLES
	  str->dPhi_KF_jj = TMath::Abs(j1.DeltaPhi(j2));
	  str->dR_KF_jj   = j1.DeltaR(j2);
	  str->dPhi_KF_ZZ = TMath::Abs(lepZ.DeltaPhi(hadZ));
	  str->dR_KF_ZZ   = lepZ.DeltaR(hadZ);
	  
	  
	  
	  if(corr_jet_pt1 > 0. && corr_jet_pt2 > 0.)
	    {	    
	      str->corrJ1_KF_m           = m_GoodJets.at(SelectedJets.first)->Get4Momentum()->M();
	      str->corrJ1_KF_pt          = corr_jet_pt1;
	      str->corrJ1_KF_eta         = m_GoodJets.at(SelectedJets.first)->righteta();
	      str->corrJ1_KF_eta_det     = Jet_1_KF->emscale_eta();
	      str->corrJ1_KF_phi         = m_GoodJets.at(SelectedJets.first)->rightphi();
	      str->corrJ1_KF_flavortruth = Jet_1_KF->flavor_truth_label();
	      
	      str->corrJ2_KF_m           = m_GoodJets.at(SelectedJets.second)->Get4Momentum()->M();
	      str->corrJ2_KF_pt          = corr_jet_pt2;
	      str->corrJ2_KF_eta         = m_GoodJets.at(SelectedJets.second)->righteta();
	      str->corrJ2_KF_eta_det     = Jet_2_KF->emscale_eta();
	      str->corrJ2_KF_phi         = m_GoodJets.at(SelectedJets.second)->rightphi();
	      str->corrJ2_KF_flavortruth = Jet_2_KF->flavor_truth_label();
	      
	      
	      j1.SetPtEtaPhiM(str->corrJ1_KF_pt,str->corrJ1_KF_eta,str->corrJ1_KF_phi,str->corrJ1_KF_m);
	      j2.SetPtEtaPhiM(str->corrJ2_KF_pt,str->corrJ2_KF_eta,str->corrJ2_KF_phi,str->corrJ2_KF_m);
	      hadZ = j1   +   j2;
	      H    = lepZ + hadZ;
	      
	      str->corrZ_KF_m   = hadZ.M();
	      str->corrZ_KF_pt  = hadZ.Pt();
	      str->corrZ_KF_eta = hadZ.Eta();
	      str->corrZ_KF_phi = hadZ.Phi();
	      
	      str->corrH_KF_m   = H.M();
	      str->corrH_KF_pt  = H.Pt();
	      str->corrH_KF_eta = H.Eta();
	      str->corrH_KF_phi = H.Phi();
	      
	      str->chisquare = ChiSq;
	    }
	  
	  
	  Analysis::Jet *jet_p = new Analysis::Jet(Jet_1_KF);
	  str->realJ1_KF_pdg = GetFlavour(jet_p).first;
	  Analysis::Jet *jet_s = new Analysis::Jet(Jet_2_KF);
	  str->realJ2_KF_pdg = GetFlavour(jet_s).first;
	  


	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  str->realJ1_LJ_m           = m_GoodJets.at(0)->Get4Momentum()->M();
	  str->realJ1_LJ_pt          = m_GoodJets.at(0)->rightpt();
	  str->realJ1_LJ_eta         = m_GoodJets.at(0)->righteta();
	  str->realJ1_LJ_phi         = m_GoodJets.at(0)->rightphi();
	  str->realJ1_LJ_eta_det     = Jet_1_LJ->emscale_eta();
	  str->realJ1_LJ_flavortruth = Jet_1_LJ->flavor_truth_label();
	  str->realJ1_LJ_jvf         = Jet_1_LJ->jvtxf();
	  str->realJ1_LJ_ntrk        = Jet_1_LJ->nTrk();
	  str->realJ1_LJ_width       = Jet_1_LJ->WIDTH();
	  str->realJ1_LJ_MV1         = GetMV1value(m_GoodJets.at(0));
	  str->realJ1_LJ_ntrk12      = InfoNtracksWidthJ1_LJ.first;
	  str->realJ1_LJ_width12     = InfoNtracksWidthJ1_LJ.second;
	  
	  str->realJ2_LJ_m           = m_GoodJets.at(1)->Get4Momentum()->M();
	  str->realJ2_LJ_pt          = m_GoodJets.at(1)->rightpt();
	  str->realJ2_LJ_eta         = m_GoodJets.at(1)->righteta();
	  str->realJ2_LJ_phi         = m_GoodJets.at(1)->rightphi();
	  str->realJ2_LJ_eta_det     = Jet_2_LJ->emscale_eta();
	  str->realJ2_LJ_flavortruth = Jet_2_LJ->flavor_truth_label();
	  str->realJ2_LJ_jvf         = Jet_2_LJ->jvtxf();
	  str->realJ2_LJ_ntrk        = Jet_2_LJ->nTrk();
	  str->realJ2_LJ_width       = Jet_2_LJ->WIDTH();
	  str->realJ2_LJ_MV1         = GetMV1value(m_GoodJets.at(1));
	  str->realJ2_LJ_ntrk12      = InfoNtracksWidthJ2_LJ.first;
	  str->realJ2_LJ_width12     = InfoNtracksWidthJ2_LJ.second;
	  
	  	  
	  TLorentzVector j1_LJ;
	  j1_LJ.SetPtEtaPhiM(str->realJ1_LJ_pt,str->realJ1_LJ_eta,str->realJ1_LJ_phi,str->realJ1_LJ_m);
	  TLorentzVector j2_LJ;
	  j2_LJ.SetPtEtaPhiM(str->realJ2_LJ_pt,str->realJ2_LJ_eta,str->realJ2_LJ_phi,str->realJ2_LJ_m);
	  TLorentzVector hadZ_LJ = j1_LJ   +   j2_LJ;
	  TLorentzVector H_LJ    = lepZ + hadZ_LJ;
	  
	  str->realZ_LJ_m     = hadZ_LJ.M();
	  str->realZ_LJ_pt    = hadZ_LJ.Pt();
	  str->realZ_LJ_eta   = hadZ_LJ.Eta();
	  str->realZ_LJ_phi   = hadZ_LJ.Phi();
	  
	  str->realH_LJ_m     = H_LJ.M();
	  str->realH_LJ_pt    = H_LJ.Pt();
	  str->realH_LJ_eta   = H_LJ.Eta();
	  str->realH_LJ_phi   = H_LJ.Phi();	  
	  
	  //ANGULAR JET VARIABLES
	  str->dPhi_LJ_jj = TMath::Abs(j1_LJ.DeltaPhi(j2_LJ));
	  str->dR_LJ_jj   = j1_LJ.DeltaR(j2_LJ);
	  str->dPhi_LJ_ZZ = TMath::Abs(lepZ.DeltaPhi(hadZ_LJ));
	  str->dR_LJ_ZZ   = lepZ.DeltaR(hadZ_LJ);	  
	  
	  Analysis::Jet *jet_p_LJ = new Analysis::Jet(Jet_1_LJ);
	  str->realJ1_LJ_pdg = GetFlavour(jet_p_LJ).first;
	  Analysis::Jet *jet_s_LJ = new Analysis::Jet(Jet_2_LJ);
	  str->realJ2_LJ_pdg = GetFlavour(jet_s_LJ).first;
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  str->realJ1_BP_m           = m_GoodJets.at(Jone)->Get4Momentum()->M();
	  str->realJ1_BP_pt          = m_GoodJets.at(Jone)->rightpt();
	  str->realJ1_BP_eta         = m_GoodJets.at(Jone)->righteta();
	  str->realJ1_BP_phi         = m_GoodJets.at(Jone)->rightphi();
	  str->realJ1_BP_eta_det     = Jet_1_BP->emscale_eta();
	  str->realJ1_BP_flavortruth = Jet_1_BP->flavor_truth_label();
	  str->realJ1_BP_jvf         = Jet_1_BP->jvtxf();
	  str->realJ1_BP_ntrk        = Jet_1_BP->nTrk();
	  str->realJ1_BP_width       = Jet_1_BP->WIDTH();
	  str->realJ1_BP_MV1         = GetMV1value(m_GoodJets.at(Jone));
	  str->realJ1_BP_ntrk12      = InfoNtracksWidthJ1_BP.first;
	  str->realJ1_BP_width12     = InfoNtracksWidthJ1_BP.second;
	  
	  str->realJ2_BP_m           = m_GoodJets.at(Jtwo)->Get4Momentum()->M();
	  str->realJ2_BP_pt          = m_GoodJets.at(Jtwo)->rightpt();
	  str->realJ2_BP_eta         = m_GoodJets.at(Jtwo)->righteta();
	  str->realJ2_BP_phi         = m_GoodJets.at(Jtwo)->rightphi();
	  str->realJ2_BP_eta_det     = Jet_2_BP->emscale_eta();
	  str->realJ2_BP_flavortruth = Jet_2_BP->flavor_truth_label();
	  str->realJ2_BP_jvf         = Jet_2_BP->jvtxf();
	  str->realJ2_BP_ntrk        = Jet_2_BP->nTrk();
	  str->realJ2_BP_width       = Jet_2_BP->WIDTH();
	  str->realJ2_BP_MV1         = GetMV1value(m_GoodJets.at(Jtwo));
	  str->realJ2_BP_ntrk12      = InfoNtracksWidthJ2_BP.first;
	  str->realJ2_BP_width12     = InfoNtracksWidthJ2_BP.second;
	  
	  	  
	  TLorentzVector j1_BP;
	  j1_BP.SetPtEtaPhiM(str->realJ1_BP_pt,str->realJ1_BP_eta,str->realJ1_BP_phi,str->realJ1_BP_m);
	  TLorentzVector j2_BP;
	  j2_BP.SetPtEtaPhiM(str->realJ2_BP_pt,str->realJ2_BP_eta,str->realJ2_BP_phi,str->realJ2_BP_m);
	  TLorentzVector hadZ_BP = j1_BP   +   j2_BP;
	  TLorentzVector H_BP    = lepZ + hadZ_BP;
	  
	  str->realZ_BP_m     = hadZ_BP.M();
	  str->realZ_BP_pt    = hadZ_BP.Pt();
	  str->realZ_BP_eta   = hadZ_BP.Eta();
	  str->realZ_BP_phi   = hadZ_BP.Phi();
	  
	  str->realH_BP_m     = H_BP.M();
	  str->realH_BP_pt    = H_BP.Pt();
	  str->realH_BP_eta   = H_BP.Eta();
	  str->realH_BP_phi   = H_BP.Phi();	  
	  
	  //ANGULAR JET VARIABLES
	  str->dPhi_BP_jj = TMath::Abs(j1_BP.DeltaPhi(j2_BP));
	  str->dR_BP_jj   = j1_BP.DeltaR(j2_BP);
	  str->dPhi_BP_ZZ = TMath::Abs(lepZ.DeltaPhi(hadZ_BP));
	  str->dR_BP_ZZ   = lepZ.DeltaR(hadZ_BP);	  
	  
	  Analysis::Jet *jet_p_BP = new Analysis::Jet(Jet_1_BP);
	  str->realJ1_BP_pdg = GetFlavour(jet_p_BP).first;
	  Analysis::Jet *jet_s_BP = new Analysis::Jet(Jet_2_BP);
	  str->realJ2_BP_pdg = GetFlavour(jet_s_BP).first;
	  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	  //////////QUARK/GLUON//////////	  
	  int idx1 = Pair_jet1;
	  int idx2 = Pair_jet2;
	  
	  if(FillGluon)
	    {	  	    
	      if (CheckMap("44p_4var",idx1,idx2)) {
		str->xWin_44p_4var = GetSOMx("44p_4var",idx1,idx2);
		str->yWin_44p_4var = GetSOMy("44p_4var",idx1,idx2);
		str->zWin_44p_4var = GetSOMz("44p_4var",idx1,idx2);
		str->gWin_44p_4var = GetSOMg("44p_4var",idx1,idx2);
	      }
	      if (CheckMap("44p_6var",idx1,idx2)) {
		str->xWin_44p_6var = GetSOMx("44p_6var",idx1,idx2);
		str->yWin_44p_6var = GetSOMy("44p_6var",idx1,idx2);
		str->zWin_44p_6var = GetSOMz("44p_6var",idx1,idx2);
		str->gWin_44p_6var = GetSOMg("44p_6var",idx1,idx2);
	      }
	      if (CheckMap("44h_4var",idx1,idx2)) {
		str->xWin_44h_4var = GetSOMx("44h_4var",idx1,idx2);
		str->yWin_44h_4var = GetSOMy("44h_4var",idx1,idx2);
		str->zWin_44h_4var = GetSOMz("44h_4var",idx1,idx2);
		str->gWin_44h_4var = GetSOMg("44h_4var",idx1,idx2);
	      }
	      if (CheckMap("44h_6var",idx1,idx2)) {
		str->xWin_44h_6var = GetSOMx("44h_6var",idx1,idx2);
		str->yWin_44h_6var = GetSOMy("44h_6var",idx1,idx2);
		str->zWin_44h_6var = GetSOMz("44h_6var",idx1,idx2);
		str->gWin_44h_6var = GetSOMg("44h_6var",idx1,idx2);
	      }
	      if (CheckMap("64p_4var",idx1,idx2)) {
		str->xWin_64p_4var = GetSOMx("64p_4var",idx1,idx2);
		str->yWin_64p_4var = GetSOMy("64p_4var",idx1,idx2);
		str->zWin_64p_4var = GetSOMz("64p_4var",idx1,idx2);
		str->gWin_64p_4var = GetSOMg("64p_4var",idx1,idx2);
	      }
	      if (CheckMap("64p_6var",idx1,idx2)) {
		str->xWin_64p_6var = GetSOMx("64p_6var",idx1,idx2);
		str->yWin_64p_6var = GetSOMy("64p_6var",idx1,idx2);
		str->zWin_64p_6var = GetSOMz("64p_6var",idx1,idx2);
		str->gWin_64p_6var = GetSOMg("64p_6var",idx1,idx2);
	      }
	      if (CheckMap("64h_4var",idx1,idx2)) {
		str->xWin_64h_4var = GetSOMx("64h_4var",idx1,idx2);
		str->yWin_64h_4var = GetSOMy("64h_4var",idx1,idx2);
		str->zWin_64h_4var = GetSOMz("64h_4var",idx1,idx2);
		str->gWin_64h_4var = GetSOMg("64h_4var",idx1,idx2);
	      }
	      if (CheckMap("64h_6var",idx1,idx2)) {
		str->xWin_64h_6var = GetSOMx("64h_6var",idx1,idx2);
		str->yWin_64h_6var = GetSOMy("64h_6var",idx1,idx2);
		str->zWin_64h_6var = GetSOMz("64h_6var",idx1,idx2);
		str->gWin_64h_6var = GetSOMg("64h_6var",idx1,idx2);
	      }
	    } // End of the FillGluon variables   
	  

	  /*
	    str->realJ1_Fisher = getFisher(reader,var1,var2,jet_pt->at(m_jetindex.at(idx1)),jet_eta->at(m_jetindex.at(idx1)),jet_nTrk->at(m_jetindex.at(idx1)),jet_width->at(m_jetindex.at(idx1)));
	    str->realJ2_Fisher = getFisher(reader,var1,var2,jet_pt->at(m_jetindex.at(idx2)),jet_eta->at(m_jetindex.at(idx2)),jet_nTrk->at(m_jetindex.at(idx2)),jet_width->at(m_jetindex.at(idx2)));
	    str->ll_2_jets = getLikelihood(reader,var1,var2,var3,var4,str->realJ1_pt,str->realJ2_pt,jet_nTrk->at(m_jetindex.at(idx1)),jet_width->at(m_jetindex.at(idx1)),jet_nTrk->at(m_jetindex.at(idx2)),jet_width->at(m_jetindex.at(idx2)));
	    
	    str->corrJ1_Fisher = getFisher(reader,var1,var2,jet_pt->at(m_jetindex.at(idx1)),jet_eta->at(m_jetindex.at(idx1)),jet_nTrk->at(m_jetindex.at(idx1)),jet_width->at(m_jetindex.at(idx1)));
	    str->corrJ2_Fisher = getFisher(reader,var1,var2,jet_pt->at(m_jetindex.at(idx2)),jet_eta->at(m_jetindex.at(idx2)),jet_nTrk->at(m_jetindex.at(idx2)),jet_width->at(m_jetindex.at(idx2)));
	    str->ll_2_jets_corr = getLikelihood(reader,var1,var2,var3,var4,str->corrJ1_pt,str->corrJ2_pt,jet_nTrk->at(m_jetindex.at(idx1)),jet_width->at(m_jetindex.at(idx1)),jet_nTrk->at(m_jetindex.at(idx2)),jet_width->at(m_jetindex.at(idx2)));
	    
	  */	
	  
	} //End of minimum number of Jets
      
      
      // btagging SF filling
      str->btagSF = tmpbtagsf;
      
      // Fill the Tree    
      analysistree->Fill();
    } // End of the If cut minimum!!!!  
}
///////////////////////////////////////////////////////////////////////////



pair <double,double> HiggsllqqAnalysis::GetJetSFsvalue(int jetindex)
{
  
  //Getting the jet Object
  D3PDReader::JetD3PDObjectElement *Jet_ = m_GoodJets.at(jetindex)->GetJet();
  
  ajet.jetPt  = m_GoodJets.at(jetindex)->rightpt();
  ajet.jetEta = m_GoodJets.at(jetindex)->righteta();
  
  Analysis::CalibResult res;
  
  /*
    std::string OperatingPoint;
    if(MV1_OP70==0.795)
    OperatingPoint = "0_795";
    
    else if(MV1_OP70==0.601713)
    OperatingPoint = "0_601713";
  */
  
  if (GetMV1value(m_GoodJets.at(jetindex)) > MV1_OP70)
    { 
      // jet flavor truth values checked on lxr
      if       (Jet_->flavor_truth_label() == 5)                                      res = calib->getScaleFactor(ajet,  "B"  , "0_795"/*OperatingPoint*/, uncertainty);
      else if  (Jet_->flavor_truth_label() == 4 || Jet_->flavor_truth_label() == 15)  res = calib->getScaleFactor(ajet,  "C"  , "0_795"/*OperatingPoint*/, uncertainty);
      else                                                                            res = calib->getScaleFactor(ajet,"Light", "0_795"/*OperatingPoint*/, uncertainty);
    } 
  else
    {
      if       (Jet_->flavor_truth_label() == 5)                                      res = calib->getInefficiencyScaleFactor(ajet,  "B"  , "0_795"/*OperatingPoint*/, uncertainty);
      else if  (Jet_->flavor_truth_label() == 4 || Jet_->flavor_truth_label() == 15)  res = calib->getInefficiencyScaleFactor(ajet,  "C"  , "0_795"/*OperatingPoint*/, uncertainty);
      else                                                                            res = calib->getInefficiencyScaleFactor(ajet,"Light", "0_795"/*OperatingPoint*/, uncertainty);
    }
  
  pair <double,double> result;
  
  result.first=res.first;
  result.second=res.second;
  
  return result;
}




/////////////////////////////////////////////////////////////////
// NTRK E WIDTH CALCULATION//
std::pair<Float_t, Float_t> HiggsllqqAnalysis::InfoTracks(Int_t JetIndex)
{  
  Float_t Ntracks=0.,W=0., Width=0., N=-1., dr=Cone_size; 
  
  for(Int_t i=0; i<ntuple->trk.n();i++){ //loop on tracks in the event    
    
    if(isGoodTrack(i)){//passes only good tracks      
      
      if(isInTheJet(i,JetIndex,ntuple->trk.eta(),ntuple->trk.phi_wrtPV())<dr){//check if track is inside jet 	
	
	Ntracks++;//count tracks	
	W+=ntuple->trk.pt()->at(i)*isInTheJet(i,JetIndex,ntuple->trk.eta(),ntuple->trk.phi_wrtPV());	
	N+=ntuple->trk.pt()->at(i);      
      }    
    }  
  }  
  
  Width=W/N;//track Width  
  
  pair<Float_t,Float_t> info;  
  info.first = Ntracks;  
  info.second = Width;  
  
  return info;
}


//GOOD TRACK//
Bool_t HiggsllqqAnalysis::isGoodTrack(Int_t TrackIndex)
{
  
  //  D3PDReader::TrackD3PDObject *track_branch(0);
  //  track_branch = &(ntuple->trk);
  
  Bool_t quality = ((ntuple->trk.chi2()->at(TrackIndex)/ntuple->trk.ndof()->at(TrackIndex)) <= 3);  
  Bool_t d0      = (ntuple->trk.d0_wrtPV()->at(TrackIndex) <= 1.0);  
  Bool_t pixel   = (ntuple->trk.nPixHits()->at(TrackIndex)>1);  
  Bool_t sct     = (ntuple->trk.nSCTHits()->at(TrackIndex)>6);  
  Bool_t pt      = (ntuple->trk.pt()->at(TrackIndex)>1000);  
  
  return (quality && d0 && pixel && sct && pt);
}


//IN THE JET..THIS METHOD CALCULATE DR...YOU CAN MODIFY WITH TLORENTZ VECTOR THAT HAS DR CALCULATION INSIDE!
Float_t HiggsllqqAnalysis::isInTheJet(Int_t Index, Int_t JetIndex, vector<float> *whatinjet_eta, vector<float> *whatinjet_phi)
{  
  Float_t deltar=0,deltaphi=0,deltaeta=0, deltaphi_tmp=0;  
  
  deltaphi_tmp = whatinjet_phi->at(Index) - m_GoodJets.at(JetIndex)->rightphi();  
  
  if(TMath::Abs(deltaphi_tmp)<TMath::Pi()) {
    deltaphi=deltaphi_tmp;
  }  
  else if(deltaphi_tmp>TMath::Pi()) {
    deltaphi=TMath::Pi()-deltaphi_tmp;
  }  
  else if(deltaphi_tmp<-TMath::Pi()) {
    deltaphi=-TMath::Pi()-deltaphi_tmp;
  }
  
  deltaeta = whatinjet_eta->at(Index)-m_GoodJets.at(JetIndex)->righteta();
  deltar   = sqrt((deltaeta*deltaeta)+(deltaphi*deltaphi));  
  
  return deltar;
}






// SOM METHODS (dec 2012)


int HiggsllqqAnalysis::GetSOMg(TString which_map, int jetidx1, int jetidx2)
{
  
  int g=0;
  pair <int, int> WinnerNode = GetSOMWinner(which_map, jetidx1,jetidx2); 
  std::vector< pair<int,int> > LinearMap;
  LinearMap.clear();
  
  LinearMap = GetCoordinates(which_map); //// change function for others linearization 
  
  for(UInt_t i=0; i<LinearMap.size(); i++) if (LinearMap.at(i)==WinnerNode) g=i;
  
  return g;
}


pair <int, int> HiggsllqqAnalysis::GetSOMWinner(TString which_map, int jetidx1, int jetidx2){
  
  pair <SOMVar, SOMVar> EventSOMVar = GetSOMVariables(jetidx1, jetidx2);
  std::vector< pair<SOMVar,SOMVar> > SOMVectors;
  std::vector< pair<int,int> > Map;
  SOMVectors.clear();
  Map.clear();
  
  SOMVectors = GetSOMVectors(which_map, jetidx1, jetidx2);
  Map = GetCoordinates(which_map);
  
  return Map.at(GetWinnerInfo(SOMVectors, EventSOMVar).first);
  
}

float HiggsllqqAnalysis::GetSOMWinnerDist(TString which_map, int jetidx1, int jetidx2){
  
  pair <SOMVar, SOMVar> EventSOMVar = GetSOMVariables(jetidx1, jetidx2);
  std::vector< pair <SOMVar, SOMVar> > SOMVectors;
  SOMVectors.clear();
  SOMVectors = GetSOMVectors(which_map, jetidx1, jetidx2);
  
  return GetWinnerInfo(SOMVectors, EventSOMVar).second;
  
}

pair <SOMVar, SOMVar> HiggsllqqAnalysis::GetSOMVariables(int jetidx1, int jetidx2){
  
  pair <SOMVar, SOMVar> result;
  SOMVar jet1SOM, jet2SOM;
  
  std::pair<Float_t, Float_t> InfoNtracksWidthJet1 = InfoTracks(jetidx1);
  std::pair<Float_t, Float_t> InfoNtracksWidthJet2 = InfoTracks(jetidx2);
  
  
  jet1SOM.Ntrk    = InfoNtracksWidthJet1.first/50.;
  jet1SOM.width   = InfoNtracksWidthJet1.second;
  jet1SOM.JetMass = GetMV1value(m_GoodJets.at(jetidx1));  //jet1SOM.JetMass = m_GoodJets.at(jetidx1)->Get4Momentum()->M();
  jet2SOM.Ntrk    = InfoNtracksWidthJet2.first/50.;
  jet2SOM.width   = InfoNtracksWidthJet2.second;
  jet2SOM.JetMass = GetMV1value(m_GoodJets.at(jetidx2));  //jet2SOM.JetMass = m_GoodJets.at(jetidx2)->Get4Momentum()->M();
  
  
  result.first = jet1SOM;
  result.second = jet2SOM;
  
  return result;
}

pair<int, float> HiggsllqqAnalysis::GetWinnerInfo(std::vector< pair<SOMVar, SOMVar> > SOMVectors, pair <SOMVar,SOMVar> Event){
  
  std::vector<float> distances;
  distances.clear();
  
  for (UInt_t i=0; i<SOMVectors.size(); i++){

    float dist2 = 
      pow(SOMVectors.at(i).first.Ntrk   - Event.first.Ntrk,2)  + 
      pow(SOMVectors.at(i).first.width  - Event.first.width,2) +
      pow(SOMVectors.at(i).second.Ntrk  - Event.second.Ntrk,2) +
      pow(SOMVectors.at(i).second.width - Event.second.width,2);
    if (SOMVectors.at(i).first.JetMass != -1. && SOMVectors.at(i).second.JetMass != -1.) dist2 += pow(SOMVectors.at(i).first.JetMass - Event.first.JetMass,2) + pow(SOMVectors.at(i).second.JetMass - Event.second.JetMass,2); 
    distances.push_back(sqrt(dist2));
  }
  
  float tmp_dist = distances.at(0);
  int idx_min    = 0;
  
  for (UInt_t i=1; i<distances.size(); i++) if (distances.at(i) < tmp_dist ) {idx_min = i; tmp_dist=distances.at(i);}
  
  pair <int, float> result;
  result.first  = idx_min;
  result.second = distances.at(idx_min);
  return result;
}



std::vector< pair<SOMVar,SOMVar> > HiggsllqqAnalysis::GetSOMVectors(TString which_map, int jetidx1, int jetidx2){
  
  TString maps_path    = "./HiggsllqqAnalysis/SOMmaps/";
  TString map_filename = GetCorrectMap(which_map, jetidx1, jetidx2);
  TString ps_dim;
  pair<SOMVar, SOMVar> tmp;
  std::vector< pair<SOMVar, SOMVar> > SOMVectors; 
  SOMVectors.clear();
  
  ifstream map;
  map.open(maps_path+map_filename);
  
  if      (which_map=="44p_4var" || which_map=="44h_4var" || which_map=="64p_4var" || which_map=="64h_4var") ps_dim = "4var";
  else if (which_map=="44p_6var" || which_map=="44h_6var" || which_map=="64p_6var" || which_map=="64h_6var") ps_dim = "6var";
  
  float Ntrk1, Ntrk2, width1, width2, JetMass1, JetMass2; 
  
  while (1){
    if      (ps_dim == "4var") map >> Ntrk1 >> width1 >> Ntrk2 >> width2;
    else if (ps_dim == "6var") map >> Ntrk1 >> width1 >> JetMass1 >> Ntrk2 >> width2 >> JetMass2;
    if (!map.good()) break;
    tmp.first.Ntrk   = Ntrk1;
    tmp.second.Ntrk  = Ntrk2;
    tmp.first.width  = width1;
    tmp.second.width = width2;
    if      (ps_dim == "6var") {tmp.first.JetMass = JetMass1; tmp.first.JetMass = JetMass2;}
    else if (ps_dim == "4var") {tmp.first.JetMass = -1.; tmp.first.JetMass = -1;}
    SOMVectors.push_back(tmp);
  }
  map.close();
  
  return SOMVectors;
}

TString HiggsllqqAnalysis::GetCorrectMap (TString which_map, int jetidx1, int jetidx2){
  
  float jet1_pt = m_GoodJets.at(jetidx1)->rightpt();
  float jet2_pt = m_GoodJets.at(jetidx2)->rightpt();
  TString name  = which_map;
  
  name += "/QGmap_";
  
  if ((20000.<jet1_pt) && (jet1_pt<40000.) && (20000.<jet2_pt) && (jet2_pt<40000.)){ name += "20_40";}
  else if ((40000.<jet1_pt) && (jet1_pt<75000.)) { name += "40_75"; }
  else if ((75000.<jet1_pt) && (jet1_pt<120000.)) { name += "75_120"; }
  else if ((120000.<jet1_pt) && (jet1_pt<180000.)){ name += "120_180";}
  else if ((180000.<jet1_pt) && (jet1_pt<240000.)){ name += "180_240";}
  else if (240000.<jet1_pt) {name += "240";}
  
  return name+".cod";
}

std::vector< pair<int,int> > HiggsllqqAnalysis::GetCoordinates(TString which_map){
  
  pair<int, int> point;
  std::vector< pair<int,int> > result;
  int xDimension=0, yDimension=0;
  
  if (which_map == "44p_4var" || which_map == "44h_4var" || which_map == "44p_6var" || which_map == "44h_6var") {xDimension = 4; yDimension = 4;}
  else if (which_map == "64p_4var" || which_map == "64h_4var" || which_map == "64p_6var" || which_map == "64h_6var") {xDimension = 6; yDimension = 4;}
  
  for (int j=0; j<yDimension; j++){
    for (int i=0; i<xDimension; i++){
      point.first = i;
      point.second = j;
      result.push_back(point);
    }
  }
  
  return result;
}

Bool_t HiggsllqqAnalysis::CheckMap(TString which_map, int jetidx1, int jetidx2){
  
  TString maps_path = "./HiggsllqqAnalysis/SOMmaps/";
  Bool_t opened = false;
  ifstream file;
  TString map_name = GetCorrectMap(which_map, jetidx1,jetidx2);
  TString name = maps_path+map_name;
  file.open(name);
  if (file.good()) opened = true; 
  file.close();
  
  return opened;
}


Float_t HiggsllqqAnalysis::Rightcut(Int_t efficiency, Float_t pt_jet, Float_t eta_jet)
{
  Int_t likebin = -1;
  Float_t cut[15]={0};
  Float_t cut0[15]={0.969,0.961,0.956,0.978,0.964,0.957,0.986,0.985,0.978,0.987,0.982,0.983,0.992,0.991,0.979};
  Float_t cut10[15]={0.895,0.88,0.875,0.914,0.899,0.89,0.949,0.93,0.914,0.961,0.949,0.94,0.966,0.954,0.938};
  Float_t cut20[15]={0.849,0.833,0.825,0.871,0.855,0.845,0.912,0.887,0.868,0.933,0.916,0.901,0.937,0.919,0.9};
  Float_t cut30[15]={0.802,0.788,0.772,0.823,0.805,0.795,0.867,0.837,0.815,0.895,0.874,0.855,0.895,0.875,0.852};
  Float_t cut40[15]={0.755,0.743,0.715,0.769,0.747,0.732,0.809,0.776,0.751,0.845,0.82,0.8,0.836,0.815,0.788};
  Float_t cut50[15]={0.704,0.695,0.658,0.702,0.686,0.663,0.731,0.7,0.676,0.772,0.752,0.727,0.752,0.732,0.71};
  Float_t cut60[15]={0.649,0.645,0.599,0.625,0.61,0.58,0.632,0.613,0.588,0.669,0.657,0.633,0.635,0.625,0.613};
  Float_t cut70[15]={0.588,0.588,0.537,0.533,0.518,0.487,0.512,0.502,0.492,0.537,0.537,0.52,0.488,0.498,0.495};
  Float_t cut80[15]={0.519,0.519,0.468,0.423,0.41,0.387,0.385,0.376,0.375,0.372,0.387,0.38,0.332,0.348,0.358};
  Float_t cut90[15]={0.425,0.426,0.376,0.285,0.279,0.265,0.243,0.235,0.241,0.211,0.221,0.229,0.191,0.203,0.213};
  
  Float_t right=0;
  
  if(efficiency==90)
    {
      memcpy(cut,cut90,sizeof(cut)); 
    }
  else if(efficiency==80)
    {
      memcpy(cut,cut80,sizeof(cut));  
    }
  else if(efficiency==70)
    {
      memcpy(cut,cut70,sizeof(cut)); 
    }
  else if(efficiency==60)
    {
      memcpy(cut,cut60,sizeof(cut)); 
    }
  else if(efficiency==50)
    {
      memcpy(cut,cut50,sizeof(cut)); 
    }
  else if(efficiency==40)
    {
      memcpy(cut,cut40,sizeof(cut)); 
    }
  else if(efficiency==30)
    {
      memcpy(cut,cut30,sizeof(cut)); 
    }
  else if(efficiency==20)
    {
      memcpy(cut,cut20,sizeof(cut)); 
    }
  else if(efficiency==10)
    {
      memcpy(cut,cut10,sizeof(cut)); 
    }
  else if(efficiency==0)
    {
      memcpy(cut,cut0,sizeof(cut)); 
    }

  if((pt_jet/1000.>19. && pt_jet/1000.<40.)&&(eta_jet>-0.8&& eta_jet<0.8))
    {
      likebin = 0;
    }
  else if((pt_jet/1000.>19. && pt_jet/1000.<40.)&&((eta_jet>-1.6&& eta_jet<-0.8)||(eta_jet>0.8&&eta_jet<1.6)))
    {
      likebin = 1;
    }
  else if((pt_jet/1000.>19. && pt_jet/1000.<40.)&&((eta_jet>-2.5&& eta_jet<-1.6)||(eta_jet>1.6&&eta_jet<2.5)))
    {
      likebin = 2;
    }
  else if((pt_jet/1000.>19. && pt_jet/1000.<40.)&&(eta_jet<-2.5&& eta_jet>2.5))
    {
      likebin = 2;
    }
  else if((pt_jet/1000.>40. && pt_jet/1000.<60.)&&(eta_jet>-0.8&& eta_jet<0.8))
    {
      likebin = 3;
    }
  else if((pt_jet/1000.>40. && pt_jet/1000.<60.)&&((eta_jet>-1.6&& eta_jet<-0.8)||(eta_jet>0.8&&eta_jet<1.6)))
    {
      likebin = 4;
    }
  else if((pt_jet/1000.>40. && pt_jet/1000.<60.)&&((eta_jet>-2.5&& eta_jet<-1.6)||(eta_jet>1.6&&eta_jet<2.5)))
    {
      likebin = 5;
    }
  else if((pt_jet/1000.>40. && pt_jet/1000.<60.)&&(eta_jet<-2.5&& eta_jet>2.5))
    {
      likebin = 5;
    }
  else if((pt_jet/1000.>60. && pt_jet/1000.<90.)&&(eta_jet>-0.8&& eta_jet<0.8))
    {
      likebin = 6;
    }
  else if((pt_jet/1000.>60. && pt_jet/1000.<90.)&&((eta_jet>-1.6&& eta_jet<-0.8)||(eta_jet>0.8&&eta_jet<1.6)))
    {
      likebin = 7;
    }
  else if((pt_jet/1000.>60. && pt_jet/1000.<90.)&&((eta_jet>-2.5&& eta_jet<-1.6)||(eta_jet>1.6&&eta_jet<2.5)))
    {
      likebin = 8;
    }
  else if((pt_jet/1000.>60. && pt_jet/1000.<90.)&&(eta_jet<-2.5&& eta_jet>2.5))
    {
      likebin = 8;
    }
  else if((pt_jet/1000.>90. && pt_jet/1000.<120.)&&(eta_jet>-0.8&& eta_jet<0.8))
    {
      likebin = 9;
    }
  else if((pt_jet/1000.>90. && pt_jet/1000.<120.)&&((eta_jet>-1.6&& eta_jet<-0.8)||(eta_jet>0.8&&eta_jet<1.6)))
    {
      likebin = 10;
    }
  else if((pt_jet/1000.>90. && pt_jet/1000.<120.)&&((eta_jet>-2.5&& eta_jet<-1.6)||(eta_jet>1.6&&eta_jet<2.5)))
    {
      likebin = 11;
    }
  else if((pt_jet/1000.>90. && pt_jet/1000.<120.)&&(eta_jet<-2.5&& eta_jet>2.5))
    {
      likebin = 11;
    }
  else if((pt_jet/1000.>120. && pt_jet/1000.<180.)&&(eta_jet>-0.8&& eta_jet<0.8))
    {
      likebin = 12;
    }
  else if((pt_jet/1000.>120. && pt_jet/1000.<180.)&&((eta_jet>-1.6&& eta_jet<-0.8)||(eta_jet>0.8&&eta_jet<1.6)))
    {
      likebin = 13;
    }
  else if((pt_jet/1000.>120. && pt_jet/1000.<180.)&&((eta_jet>-2.5&& eta_jet<-1.6)||(eta_jet>1.6&&eta_jet<2.5)))
    {
      likebin = 14;
    }
  else 
    {
      likebin = 14;
    }  
  right=cut[likebin];
  return right;
}


//END QUARK/GLOUN METHODS



Bool_t HiggsllqqAnalysis::Pair_Quality(){
  
  Bool_t GoodQ = false;
  
  D3PDReader::MuonD3PDObjectElement     *mu_1;
  D3PDReader::MuonD3PDObjectElement     *mu_2;
  D3PDReader::ElectronD3PDObjectElement *el_1;
  D3PDReader::ElectronD3PDObjectElement *el_2;
  
  
  if (getChannel() == HiggsllqqAnalysis::MU2 && m_GoodMuons.size() == 2) {
    mu_1 = m_GoodMuons.at(0)->GetMuon();
    mu_2 = m_GoodMuons.at(1)->GetMuon();     
  }
  if (getChannel() == HiggsllqqAnalysis::E2 && m_GoodElectrons.size() == 2) {
    el_1 = m_GoodElectrons.at(0)->GetElectron();
    el_2 = m_GoodElectrons.at(1)->GetElectron();
  }  
  
  
  if (getChannel() == HiggsllqqAnalysis::MU2 && m_GoodMuons.size() == 2 &&
      ((m_GoodMuons.at(0)->family() != Muon::CALO && mu_1->isStandAloneMuon()==0 && (mu_1->isCombinedMuon()==1 || mu_1->isSegmentTaggedMuon()==1) && mu_1->pt() > 25000. && TMath::Abs(mu_1->eta()) < 2.5) 
       ||
       (m_GoodMuons.at(1)->family() != Muon::CALO && mu_2->isStandAloneMuon()==0 && (mu_2->isCombinedMuon()==1 || mu_2->isSegmentTaggedMuon()==1) && mu_2->pt() > 25000. && TMath::Abs(mu_2->eta()) < 2.5)))
    GoodQ = true;
  
  
  if (getChannel() == HiggsllqqAnalysis::E2 && m_GoodElectrons.size() == 2 &&       
      ((el_1->pt()>25000. && isMediumPlusPlus(el_1->etas2(),
					      el_1->cl_E() / TMath::CosH(el_1->etas2()),
					      el_1->f3(),
					      el_1->Ethad() / (el_1->cl_E() / TMath::CosH(el_1->etas2())),
					      el_1->Ethad1() / (el_1->cl_E() / TMath::CosH(el_1->etas2())),
					      el_1->reta(),
					      el_1->weta2(),
					      el_1->f1(),
					      el_1->wstot(),
					      el_1->emaxs1() + el_1->Emax2() > 0 ?(el_1->emaxs1() - el_1->Emax2()) / (el_1->emaxs1() + el_1->Emax2()):0.,
					      el_1->deltaeta1(),
					      el_1->trackd0_physics(),
					      el_1->TRTHighTOutliersRatio(),
					      el_1->nTRTHits(),
					      el_1->nTRTOutliers(),
					      el_1->nSiHits(),
					      el_1->nSCTOutliers() + el_1->nPixelOutliers(),
					      el_1->nPixHits(),
					      el_1->nPixelOutliers(),
					      el_1->nBLHits(),
					      el_1->nBLayerOutliers(),
					      el_1->expectBLayerHit(),
					      egammaMenu::eg2012,
					      false,
					      false)) 
       ||
       (el_2->pt()>25000. && isMediumPlusPlus(el_2->etas2(),
					      el_2->cl_E() / TMath::CosH(el_2->etas2()),
					      el_2->f3(),
					      el_2->Ethad() / (el_2->cl_E() / TMath::CosH(el_2->etas2())),
					      el_2->Ethad1() / (el_2->cl_E() / TMath::CosH(el_2->etas2())),
					      el_2->reta(),
					      el_2->weta2(),
					      el_2->f1(),
					      el_2->wstot(),
					      el_2->emaxs1() + el_2->Emax2() > 0 ?(el_2->emaxs1() - el_2->Emax2()) / (el_2->emaxs1() + el_2->Emax2()):0.,
					      el_2->deltaeta1(),
					      el_2->trackd0_physics(),
					      el_2->TRTHighTOutliersRatio(),
					      el_2->nTRTHits(),
					      el_2->nTRTOutliers(),
					      el_2->nSiHits(),
					      el_2->nSCTOutliers() + el_2->nPixelOutliers(),
					      el_2->nPixHits(),
					      el_2->nPixelOutliers(),
					      el_2->nBLHits(),
					      el_2->nBLayerOutliers(),
					      el_2->expectBLayerHit(),
					      egammaMenu::eg2012,
					      false,
					      false))))
    GoodQ = true;
  
  return GoodQ;
}



//Best Pair method (May 2013)
void HiggsllqqAnalysis::JetBestPairResult()
{
  TLorentzVector j1;
  TLorentzVector j2;
  Jone =800;
  Jtwo =900;
  int ii(-1),jj(-1);
  float had(0);
  
  std::vector<Analysis::Jet *>::iterator jet_itr_a;
  std::vector<Analysis::Jet *>::iterator jet_itr_b;
  
  for(jet_itr_a = m_GoodJets.begin(); jet_itr_a != m_GoodJets.end(); ++jet_itr_a)
    {
      ii++;
      j1.SetPtEtaPhiM((*jet_itr_a)->rightpt(),(*jet_itr_a)->righteta(),(*jet_itr_a)->rightphi(),(*jet_itr_a)->Get4Momentum()->M());
      
      for(jet_itr_b = m_GoodJets.begin(); jet_itr_b != m_GoodJets.end(); ++jet_itr_b)
	{	
	  jj++;
	  
	  j2.SetPtEtaPhiM((*jet_itr_b)->rightpt(),(*jet_itr_b)->righteta(),(*jet_itr_b)->rightphi(),(*jet_itr_b)->Get4Momentum()->M());
	  
	  TLorentzVector hadZ = j1 + j2;
	  
	  if((hadZ.M()>Mjj_low_min  && hadZ.M()<Mjj_low_max  && GetDoLowMass() && ii<jj) || (hadZ.M()>Mjj_high_min && hadZ.M()<Mjj_high_max && !GetDoLowMass() && ii<jj)
	     )
	    {
	      int tmpJone = ii;
	      int tmpJtwo = jj;
	      
	      if(((tmpJtwo-tmpJone)<(Jtwo-Jone)) && ((Jtwo-Jone)>0))
		{
		  Jone = tmpJone;
		  Jtwo = tmpJtwo;
		  had = hadZ.M();
		}
	    }
	}
      jj=-1;
    }
  //if(Jone!=0 && Jtwo!=1)
  //cout<<"Jet 1= "<<Jone<<". Jet 2 = "<<Jtwo<<". Mass = "<<had<<endl;
}


//TRIGGER MATCHING AND SF's


// Bool_t HiggsllqqAnalysis::isTriggerMatched()
// {
//   // set the trigger navigation variables
  
//   // we need non-const access to the variables
//   D3PDReader::HSG2EventReader& ncevent = const_cast< D3PDReader::HSG2EventReader& >(*ntuple);
  
//   // set the main navigation variables
//   UInt_t smk = static_cast< UInt_t >(m_TrigDecisionToolD3PD->GetSMK());
//   m_trigNavVar.set_trig_DB_SMK(smk);
//   Int_t trig_Nav_n = ntuple->trig_Nav.n();
//   m_trigNavVar.set_trig_Nav_n(trig_Nav_n);
//   std::vector< short >* trig_Nav_chain_ChainId = ncevent.trig_Nav.chain_ChainId();
//   m_trigNavVar.set_trig_Nav_chain_ChainId(trig_Nav_chain_ChainId);
//   std::vector< std::vector< int > >* trig_Nav_chain_RoIType = ncevent.trig_Nav.chain_RoIType();
//   m_trigNavVar.set_trig_Nav_chain_RoIType(trig_Nav_chain_RoIType);
//   std::vector< std::vector< int > >* trig_Nav_chain_RoIIndex = ncevent.trig_Nav.chain_RoIIndex();
//   m_trigNavVar.set_trig_Nav_chain_RoIIndex(trig_Nav_chain_RoIIndex);
  
//   // set the electron variables
//   std::vector< std::vector< int > >* trig_RoI_EF_e_egammaContainer_egamma_Electrons = ncevent.trig_RoI_EF_e.egammaContainer_egamma_Electrons();
//   std::vector< std::vector< int > >* trig_RoI_EF_e_egammaContainer_egamma_ElectronsStatus = ncevent.trig_RoI_EF_e.egammaContainer_egamma_ElectronsStatus();
//   Int_t trig_EF_el_n = ntuple->trig_EF_el.n();
//   std::vector< float >* trig_EF_el_eta = ncevent.trig_EF_el.eta();
//   std::vector< float >* trig_EF_el_phi = ncevent.trig_EF_el.phi();
//   std::vector< float >* trig_EF_el_Et = ncevent.trig_EF_el.Et();
  
//   m_trigNavVar.set_trig_RoI_EF_e_egammaContainer_egamma_Electrons(trig_RoI_EF_e_egammaContainer_egamma_Electrons);
//   m_trigNavVar.set_trig_RoI_EF_e_egammaContainer_egamma_ElectronsStatus(trig_RoI_EF_e_egammaContainer_egamma_ElectronsStatus);
//   m_trigNavVar.set_trig_EF_el_n(trig_EF_el_n);
//   m_trigNavVar.set_trig_EF_el_eta(trig_EF_el_eta);
//   m_trigNavVar.set_trig_EF_el_phi(trig_EF_el_phi);
//   m_trigNavVar.set_trig_EF_el_Et(trig_EF_el_Et);
  
//   // set the muon variables:
//   std::vector< int >* trig_RoI_EF_mu_Muon_ROI = ncevent.trig_RoI_EF_mu.Muon_ROI();
//   std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFInfoContainer = 0;
//   if (analysis_version() == "rel_17_2")trig_RoI_EF_mu_TrigMuonEFInfoContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer();
//   std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = 0;
//   if (analysis_version() == "rel_17_2")trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainerStatus();
//   std::vector< int >* trig_RoI_L2_mu_CombinedMuonFeature = ncevent.trig_RoI_L2_mu.CombinedMuonFeature();
//   std::vector< int >* trig_RoI_L2_mu_CombinedMuonFeatureStatus = ncevent.trig_RoI_L2_mu.CombinedMuonFeatureStatus();
//   std::vector< int >* trig_RoI_L2_mu_MuonFeature = ncevent.trig_RoI_L2_mu.MuonFeature();
//   std::vector< int >* trig_RoI_L2_mu_Muon_ROI = ncevent.trig_RoI_L2_mu.Muon_ROI();
//   std::vector< std::vector< float > >* trig_EF_trigmuonef_track_CB_pt = ncevent.trig_EF_trigmuonef.track_CB_pt();
//   std::vector< std::vector< float > >* trig_EF_trigmuonef_track_CB_eta = ncevent.trig_EF_trigmuonef.track_CB_eta();
//   std::vector< std::vector< float > >* trig_EF_trigmuonef_track_CB_phi = ncevent.trig_EF_trigmuonef.track_CB_phi();
//   std::vector< std::vector< float > >* trig_EF_trigmuonef_track_SA_pt = ncevent.trig_EF_trigmuonef.track_SA_pt();
//   std::vector< std::vector< float > >* trig_EF_trigmuonef_track_SA_eta = ncevent.trig_EF_trigmuonef.track_SA_eta();
//   std::vector< std::vector< float > >* trig_EF_trigmuonef_track_SA_phi = ncevent.trig_EF_trigmuonef.track_SA_phi();
//   std::vector< std::vector< float > >* trig_EF_trigmugirl_track_CB_pt = ncevent.trig_EF_trigmugirl.track_CB_pt();
//   std::vector< std::vector< float > >* trig_EF_trigmugirl_track_CB_eta = ncevent.trig_EF_trigmugirl.track_CB_eta();
//   std::vector< std::vector< float > >* trig_EF_trigmugirl_track_CB_phi = ncevent.trig_EF_trigmugirl.track_CB_phi();
//   std::vector< float >* trig_L2_combmuonfeature_eta = ncevent.trig_L2_combmuonfeature.eta();
//   std::vector< float >* trig_L2_combmuonfeature_phi = ncevent.trig_L2_combmuonfeature.phi();
//   std::vector< float >* trig_L2_muonfeature_eta = ncevent.trig_L2_muonfeature.eta();
//   std::vector< float >* trig_L2_muonfeature_phi = ncevent.trig_L2_muonfeature.phi();
//   std::vector< float >* trig_L1_mu_eta = ncevent.trig_L1_mu.eta();
//   std::vector< float >* trig_L1_mu_phi = ncevent.trig_L1_mu.phi();
//   std::vector< std::string >* trig_L1_mu_thrName = ncevent.trig_L1_mu.thrName();
  
//   std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFIsolationContainer = 0;
//   if (ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainer.IsAvailable()) {
//     trig_RoI_EF_mu_TrigMuonEFIsolationContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainer();
//   } else {
//     //Warning("isTriggerMatched", "unable to set TrigMuonEFIsolationContainer");
//   }
//   std::vector< std::vector< int > >* trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus = 0;
//   if (ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainerStatus.IsAvailable()) {
//     trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFIsolationContainerStatus();
//   } else {
//     //Warning("isTriggerMatched", "unable to set TrigMuonEFIsolationContainerStatus");
//   }
  
//   std::vector< int >* trig_EF_trigmuonef_EF_mu24i_tight = 0;
//   if (ncevent.trig_EF_trigmuonef.EF_mu24i_tight.IsAvailable()) {
//     trig_EF_trigmuonef_EF_mu24i_tight = ncevent.trig_EF_trigmuonef.EF_mu24i_tight();
//   } else {
//     //Warning("isTriggerMatched", "unable to set EF_mu24i_tight");
//   }
//   std::vector< int >* trig_EF_trigmuonef_EF_mu36_tight = 0;
//   if (ncevent.trig_EF_trigmuonef.EF_mu36_tight.IsAvailable()) {
//     trig_EF_trigmuonef_EF_mu36_tight = ncevent.trig_EF_trigmuonef.EF_mu36_tight();
//   } else {
//     //Warning("isTriggerMatched", "unable to set EF_mu36_tight");
//   }
//   std::vector<std::vector< int > >* trig_EF_trigmuonef_track_MuonType = ncevent.trig_EF_trigmuonef.track_MuonType();
  
//   m_trigNavVar.set_trig_RoI_EF_mu_Muon_ROI(trig_RoI_EF_mu_Muon_ROI);
//   m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(trig_RoI_EF_mu_TrigMuonEFInfoContainer);
//   m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus);
//   m_trigNavVar.set_trig_RoI_L2_mu_CombinedMuonFeature(trig_RoI_L2_mu_CombinedMuonFeature);
//   m_trigNavVar.set_trig_RoI_L2_mu_CombinedMuonFeatureStatus(trig_RoI_L2_mu_CombinedMuonFeatureStatus);
//   m_trigNavVar.set_trig_RoI_L2_mu_MuonFeature(trig_RoI_L2_mu_MuonFeature);
//   m_trigNavVar.set_trig_RoI_L2_mu_Muon_ROI(trig_RoI_L2_mu_Muon_ROI);
//   m_trigNavVar.set_trig_EF_trigmuonef_track_CB_pt(trig_EF_trigmuonef_track_CB_pt);
//   m_trigNavVar.set_trig_EF_trigmuonef_track_CB_eta(trig_EF_trigmuonef_track_CB_eta);
//   m_trigNavVar.set_trig_EF_trigmuonef_track_CB_phi(trig_EF_trigmuonef_track_CB_phi);
//   m_trigNavVar.set_trig_EF_trigmuonef_track_SA_pt(trig_EF_trigmuonef_track_SA_pt);
//   m_trigNavVar.set_trig_EF_trigmuonef_track_SA_eta(trig_EF_trigmuonef_track_SA_eta);
//   m_trigNavVar.set_trig_EF_trigmuonef_track_SA_phi(trig_EF_trigmuonef_track_SA_phi);
//   m_trigNavVar.set_trig_EF_trigmugirl_track_CB_pt(trig_EF_trigmugirl_track_CB_pt);
//   m_trigNavVar.set_trig_EF_trigmugirl_track_CB_eta(trig_EF_trigmugirl_track_CB_eta);
//   m_trigNavVar.set_trig_EF_trigmugirl_track_CB_phi(trig_EF_trigmugirl_track_CB_phi);
//   m_trigNavVar.set_trig_L2_combmuonfeature_eta(trig_L2_combmuonfeature_eta);
//   m_trigNavVar.set_trig_L2_combmuonfeature_phi(trig_L2_combmuonfeature_phi);
//   m_trigNavVar.set_trig_L2_muonfeature_eta(trig_L2_muonfeature_eta);
//   m_trigNavVar.set_trig_L2_muonfeature_phi(trig_L2_muonfeature_phi);
//   m_trigNavVar.set_trig_L1_mu_eta(trig_L1_mu_eta);
//   m_trigNavVar.set_trig_L1_mu_phi(trig_L1_mu_phi);
//   m_trigNavVar.set_trig_L1_mu_thrName(trig_L1_mu_thrName);
//   m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFIsolationContainer(trig_RoI_EF_mu_TrigMuonEFIsolationContainer);
//   m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus(trig_RoI_EF_mu_TrigMuonEFIsolationContainerStatus);
//   m_trigNavVar.set_trig_EF_trigmuonef_EF_mu24i_tight(trig_EF_trigmuonef_EF_mu24i_tight);
//   m_trigNavVar.set_trig_EF_trigmuonef_EF_mu36_tight(trig_EF_trigmuonef_EF_mu36_tight);
//   m_trigNavVar.set_trig_EF_trigmuonef_track_MuonType(trig_EF_trigmuonef_track_MuonType);
  
  
  
// // aqui, cambiar por good 1 y 2  //arturo
//   std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
  
//   // trigger matches
//   Bool_t single_el(kFALSE);
//   Bool_t single_mu(kFALSE);
//   Bool_t di_el(kFALSE);
//   Bool_t di_mu(kFALSE);
//   Bool_t el_mu(kFALSE);
  
//   // trigger chains
//   std::vector<TString> chainlist_singlemu = getListOfAlternativeTriggers(getSingleMuonTriggerName());
//   std::vector<TString> chainlist_dimu = getListOfAlternativeTriggers(getDiMuonTriggerName());
//   std::vector<TString> chainlist_singleel = getListOfAlternativeTriggers(getSingleElectronTriggerName());
//   std::vector<TString> chainlist_diel = getListOfAlternativeTriggers(getDiElectronTriggerName());
//   std::vector<TString> chainlist_elmu = getListOfAlternativeTriggers(getElectronMuonTriggerName());
  
//   // thresholds (manually applied to stay on plateau)
//   std::map<TString, Float_t> singlemu_thr;
//   std::map<TString, Float_t> singleel_thr;
//   std::map<TString, std::pair<Float_t, Float_t> > dimu_thr;
//   std::map<TString, std::pair<Float_t, Float_t> > diel_thr;
//   std::map<TString, std::pair<Float_t, Float_t> > elmu_thr;
  
  
//   singlemu_thr["EF_mu18_MG"] = -1.;
//   singlemu_thr["EF_mu18_MG_medium"] = -1.;
//   singlemu_thr["EF_mu24i_tight"] = -1.;
//   singlemu_thr["EF_mu36_tight"] = -1.;
  
//   singleel_thr["EF_e_medium"] = -1.;
//   singleel_thr["EF_e22_medium"] = -1.;
//   singleel_thr["EF_e22vh_medium1"] = -1.;
//   singleel_thr["EF_e24vhi_medium1"] = -1.;
//   singleel_thr["EF_e60_medium1"] = -1.;
  
//   diel_thr["EF_2e12_medium"] = std::make_pair(-1., -1.);
//   diel_thr["EF_2e12T_medium"] = std::make_pair(-1., -1.);
//   diel_thr["EF_2e12Tvh_medium"] = std::make_pair(-1., -1.);
//   diel_thr["EF_2e12Tvh_loose1"] = std::make_pair(-1., -1.);
  
//   dimu_thr["EF_2mu10_loose"] = std::make_pair(-1., -1.);
//   dimu_thr["EF_2mu13"] = std::make_pair(-1., -1.);
//   dimu_thr["EF_mu18_tight_mu8_EFFS"] = std::make_pair(-1., -1.);
  
//   elmu_thr["EF_e12Tvh_medium1_mu8"] = std::make_pair(-1., -1.);
//   elmu_thr["EF_e24vhi_loose1_mu8"] = std::make_pair(-1., -1.);
  
  
//   for (UInt_t i = 0; i < leptons.size(); i++) {
//     Int_t my_flavor = leptons[i]->flavor();
    
//     if (my_flavor == Analysis::ChargedLepton::MUON && leptons[i]->family() != Muon::CALO) {
//       // single muon trigger
//       if (analysis_version() == "rel_17") { // bug in p956
// 	trig_RoI_EF_mu_TrigMuonEFInfoContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfo();
// 	trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_eMuonEFInfoStatus();
// 	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(trig_RoI_EF_mu_TrigMuonEFInfoContainer);
// 	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus);
//       }
      
//       for (std::vector<TString>::iterator singlemu = chainlist_singlemu.begin(); singlemu != chainlist_singlemu.end(); ++singlemu) {
// 	if (m_MuonTriggerMatchTool->match(leptons[i]->Get4Momentum()->Eta(), leptons[i]->Get4Momentum()->Phi(), singlemu->Data())) {
// 	  if (leptons[i]->Get4Momentum()->Pt() > singlemu_thr[*singlemu])
// 	    single_mu = kTRUE;
// 	} // trigger matched
//       } // loop over chains in OR
      
//       // dimuon trigger
//       if (analysis_version() == "rel_17") { // bug in p956
// 	trig_RoI_EF_mu_TrigMuonEFInfoContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfo();
// 	trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfoStatus();
// 	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(trig_RoI_EF_mu_TrigMuonEFInfoContainer);
// 	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus);
//       }
//       for (UInt_t j = i + 1; j < leptons.size(); j++) {
// 	if (leptons[j]->flavor() == my_flavor && leptons[j]->family() != Muon::CALO) {
// 	  std::pair<Bool_t, Bool_t> lep1(kFALSE, kFALSE), lep2(kFALSE, kFALSE);
	  
// 	  for (std::vector<TString>::iterator dimu = chainlist_dimu.begin(); dimu != chainlist_dimu.end(); ++dimu) {
// 	    m_MuonTriggerMatchTool->matchDimuon(*(leptons[i]->Get4Momentum()), *(leptons[j]->Get4Momentum()), dimu->Data(), lep1, lep2);
	    
// 	    Bool_t is_actually_matched(kFALSE);
// 	    Bool_t is_actually_above_threshold(kFALSE);
	    
// 	    if (analysis_version() == "rel_17") {
// 	      is_actually_matched = (lep1.first && lep2.first);
	      
// 	      Float_t max_pt = TMath::Max(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
// 	      Float_t min_pt = TMath::Min(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
// 	      is_actually_above_threshold = (max_pt > dimu_thr[*dimu].first && min_pt > dimu_thr[*dimu].second);
// 	    } else if (analysis_version() == "rel_17_2") {
// 	      is_actually_matched = ((lep1.first && lep2.second) || (lep1.second && lep2.first));
	      
// 	      Float_t pt_1 = leptons[i]->Get4Momentum()->Pt();
// 	      Float_t pt_2 = leptons[j]->Get4Momentum()->Pt();
// 	      is_actually_above_threshold = (lep1.first && lep2.second && pt_1 > dimu_thr[*dimu].first && pt_2 > dimu_thr[*dimu].second);
// 	      is_actually_above_threshold = is_actually_above_threshold || (lep2.first && lep1.second && pt_2 > dimu_thr[*dimu].first && pt_1 > dimu_thr[*dimu].second);
// 	    }
	    
// 	    if (is_actually_matched && is_actually_above_threshold) {
// 	      di_mu = kTRUE;
// 	    } // dimuon trigger matched
// 	  } // loop over chains in OR
// 	} // same flavor CB/ST/SA muons
//       } // dimuon trigger
//     } // CB/ST/SA muons
//     else if (my_flavor == Analysis::ChargedLepton::ELECTRON) {
//       // single electron trigger
//       for (std::vector<TString>::iterator singleel = chainlist_singleel.begin(); singleel != chainlist_singleel.end(); ++singleel) {
// 	if (m_ElectronTriggerMatchTool->match(leptons[i]->Get4Momentum()->Eta(), leptons[i]->Get4Momentum()->Phi(), singleel->Data())) {
// 	  if (leptons[i]->Get4Momentum()->Pt() > singleel_thr[*singleel])
// 	    single_el = kTRUE;
// 	} // trigger matched
//       } // loop over chains in OR
      
//       // dielectron trigger
//       for (UInt_t j = i + 1; j < leptons.size(); j++) {
// 	if (leptons[j]->flavor() == my_flavor) {
// 	  Bool_t lep1(kFALSE), lep2(kFALSE);
	  
// 	  for (std::vector<TString>::iterator diel = chainlist_diel.begin(); diel != chainlist_diel.end(); ++diel) {
// 	    m_ElectronTriggerMatchTool->matchDielectron(*(leptons[i]->Get4Momentum()), *(leptons[j]->Get4Momentum()), diel->Data(), lep1, lep2);
	    
// 	    if (lep1 && lep2) {
// 	      Float_t max_pt = TMath::Max(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
// 	      Float_t min_pt = TMath::Min(leptons[i]->Get4Momentum()->Pt(), leptons[j]->Get4Momentum()->Pt());
	      
// 	      if (max_pt > diel_thr[*diel].first && min_pt > diel_thr[*diel].second)
// 		di_el = kTRUE;
// 	    } // dielectron trigger matched
// 	  } // loop over chains in OR
// 	} // electron pair
//       } // dielectron trigger
      
//       // electron-muon trigger
//       if (analysis_version() == "rel_17") { // bug in p956
// 	trig_RoI_EF_mu_TrigMuonEFInfoContainer = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfo();
// 	trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus = ncevent.trig_RoI_EF_mu.TrigMuonEFInfoContainer_MuonEFInfoStatus();
// 	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainer(trig_RoI_EF_mu_TrigMuonEFInfoContainer);
// 	m_trigNavVar.set_trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus(trig_RoI_EF_mu_TrigMuonEFInfoContainerStatus);
//       }
//       for (UInt_t j = 0; j < leptons.size(); j++) {
// 	if (leptons[j]->flavor() != my_flavor && leptons[j]->family() != Muon::CALO) {
// 	  Bool_t lep1lep2(kFALSE);
	  
// 	  for (std::vector<TString>::iterator elmu = chainlist_elmu.begin(); elmu != chainlist_elmu.end(); ++elmu) {
// 	    // TO BE IMPLEMENTED
// 	    lep1lep2 = m_ElectronTriggerMatchTool->matchElectronMuon(*(leptons[i]->Get4Momentum()), *(leptons[j]->Get4Momentum()), elmu->Data());
	    
// 	    if (lep1lep2) {
// 	      if (leptons[i]->Get4Momentum()->Pt() > elmu_thr[*elmu].first && leptons[j]->Get4Momentum()->Pt() > elmu_thr[*elmu].second)
// 		el_mu = kTRUE;
// 	    } // el-mu trigger matched
// 	  } // loop over chains in OR
// 	} // electron-muon pair
//       } // electron-muon trigger
//     } // electrons
//   }
  
  
//   //if (passesSingleMuonTrigger() || passesDiMuonTrigger()) {
//   //   return (single_mu || di_mu);
//   //} else if (passesSingleMuonTrigger() || passesDiMuonTrigger() || passesSingleElectronTrigger() || passesDiElectronTrigger()) {
//   return (single_mu || di_mu || single_el || di_el || el_mu); // kostas
//   //} else if (passesSingleElectronTrigger() || passesDiElectronTrigger()) {
//   //   return (single_el || di_el);
//   //}
  
//   // should have never been reached
//   return kFALSE;
// }





// Float_t HiggsllqqAnalysis::getCandidateTriggerSF(Analysis::Quadrilepton * higgs)
// {
//   Float_t result(1);
  
//   if (isMC()) {
//     if (passesSingleMuonTrigger() || passesSingleElectronTrigger()) {
      
//       // find a fake run number representing the data period to which this MC event is somewhat associated
//       m_PileupReweighter->SetRandomSeed(314159 + ntuple->eventinfo.mc_channel_number() * 2718 + ntuple->eventinfo.EventNumber());
//       Int_t representative_run_number = m_PileupReweighter->GetRandomRunNumber(ntuple->eventinfo.RunNumber());
      
//       electron_quality el_quality;
      
//       if (analysis_version() == "rel_17")
// 	el_quality = loosepp;
//       else if (analysis_version() == "rel_17_2")
// 	el_quality = ML;
      
      
//       // use MeV (set to true to use GeV)
//       //m_MuonTrigSF->setThresholds(false, representative_run_number);
      
//       std::vector<TLorentzVector> muons_4m;
//       std::vector<TLorentzVector> electrons_4m;
      
//       std::vector<Analysis::ChargedLepton *> leptons = higgs->GetLeptonVector();
      
//       for (UInt_t i = 0; i < leptons.size(); i++) {
// 	if (leptons[i]->flavor() == Analysis::ChargedLepton::MUON)
// 	  muons_4m.push_back(*(leptons[i]->Get4Momentum()));
// 	else
// 	  electrons_4m.push_back(*(leptons[i]->Get4Momentum()));
//       }
      
//       result = m_MuonTrigSF->GetTriggerSF(representative_run_number, false, muons_4m, loose, electrons_4m, el_quality, 0).first;
//     } // passes single lepton trigger for 2011
//   }
  
//   return result;
// }
// 
