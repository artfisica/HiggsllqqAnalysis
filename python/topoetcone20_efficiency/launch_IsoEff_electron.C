{
/*
  TChain H130("leptontree");
  TChain H120("leptontree");
  TChain Zee("leptontree");

  H130.Add("user.vippolit.004320*root*");
  H120.Add("user.vippolit.004322*root*");
  Zee.Add("user.vippolit.004324*root*");

  gROOT->ProcessLine(".L IsolationEfficiency.cxx+");

  Zee.SetProof();

  IsolationEfficiency eff;
//eff.SetTag("H130");
//H130.Process(&eff);
//eff.SetTag("H120");
//H120.Process(&eff);
  eff.SetTag("Zee");
  Zee.Process(&eff);
  */

//TProof *p = TProof::Open(gSystem->GetFromPipe("pod-info -c"));
//TProof *p = TProof::Open("");
  TFileCollection *Zee_coll = new TFileCollection("Zee", "Pythia Z->ee", "input.topoiso.Zee.txt");
  TFileCollection *H120_coll = new TFileCollection("H120", "Pythia H120", "input.topoiso.H120.txt");
  TFileCollection *H130_coll = new TFileCollection("H130", "Pythia H130", "input.topoiso.H130.txt");

  TChain Zee("leptontree");
  Zee.AddFileInfoList(Zee_coll->GetList());
  TChain H120("leptontree");
  H120.AddFileInfoList(H120_coll->GetList());
  TChain H130("leptontree");
  H130.AddFileInfoList(H130_coll->GetList());

  gROOT->ProcessLine(".L Higgs4lepAnalysis/python/topoetcone20_efficiency/IsolationEfficiency.cxx+");

  IsolationEfficiency eff;
  eff.SetFlavor("MUON");
//eff.SetFlavor("ELECTRON");

  eff.SetTag("Zee");
  Zee.Process(&eff);

  eff.SetTag("H120");
  H120.Process(&eff);

  eff.SetTag("H130");
  H130.Process(&eff);

//p->Process(Zee, "IsolationEfficiency.cxx+");
}
