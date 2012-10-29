{
//gSystem.AddIncludePath("-IHiggsZZ4lUtils");
//gSystem.Load("libHiggsZZ4lUtils");
//ROOT.gSystem.Load("libCint");
//gROOT.ProcessLine(".L HiggsZZ4lUtils/HiggsZZ4lUtils/ElectronMCClassification.h+");


  gROOT->ProcessLine(".L Higgs4lepAnalysis/python/topoetcone20_efficiency/IsolationEfficiency.cxx+");

//TFileCollection *Zmumu_coll = new TFileCollection("Zmumu", "Pythia Z->mumu", "input.topoiso.Zmumu.txt");
//TFileCollection *ggH125_coll = new TFileCollection("ggH125", "Pythia ggH125", "input.topoiso.ggH125.txt");

  TChain Zmumu("leptontree");
//Zmumu.AddFileInfoList(Zmumu_coll->GetList());
  Zmumu.Add("/tmp/vippolit/user.vippolit.Iso000007.147807.PowhegPythia8_AU2CT10_Zmumu.p1044.2012_deltar_inclusivelep.01.120610025643/*root*");
  Zmumu.Add("/tmp/vippolit/user.vippolit.Iso000007.147807.PowhegPythia8_AU2CT10_Zmumu.p1044.2012_deltar_inclusivelep.01.120610025702/*root*");

  TChain ggH125("leptontree");
//ggH125.AddFileInfoList(ggH125_coll->GetList());
  ggH125.Add("/tmp/vippolit/user.vippolit.Iso000007.160155.PowhegPythia8_AU2CT10_ggH125_ZZ4lep.p1044.2012_deltar_inclusivelep.01.120610031705/*root*");


  // set the efficiency tool
  IsolationEfficiency eff;
  eff.SetFlavor("ELECTRON");
  eff.SetTruth(kTRUE);

  // signal electrons on ggH125
  eff.SetClassification(kTRUE, kFALSE, kFALSE); // e-c-f
  eff.SetTag("ggH125");
  ggH125.Process(&eff);

  // background electrons on Z->mumu
  eff.SetClassification(kFALSE, kTRUE, kTRUE); // e-c-f
  eff.SetTag("Zmumu");
  Zmumu.Process(&eff);

  // background electrons on Z->mumu (conversions only)
  eff.SetClassification(kFALSE, kTRUE, kFALSE); // e-c-f
  eff.SetTag("Zmumu_conversions");
  Zmumu.Process(&eff);

  // background electrons on Z->mumu (fakes only)
  eff.SetClassification(kFALSE, kFALSE, kTRUE); // e-c-f
  eff.SetTag("Zmumu_fakes");
  Zmumu.Process(&eff);
}
