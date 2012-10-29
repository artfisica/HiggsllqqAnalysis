#include "th1fmorph.C"

int Interpolate125(){
  std::cout << "copying __yields_tmp_outfile.root as *.backup" << std::endl;
  gSystem->Exec("cp __yields_tmp_outfile.root __yields_tmp_outfile.backup");
  
  TFile out("__yields_tmp_outfile.root", "UPDATE");
  TH1F *h120 = (TH1F*) out.Get("1200_higgs_2") ;
  TH1F *h130 = (TH1F*) out.Get("1300_higgs_2") ;
  TH1F* h125_old = (TH1F*) out.Get("1250_higgs_2");

  Double_t norm = h125_old->Integral();
  h125_old->Delete();

  TH1F* h125 =th1fmorph("1250_higgs_2", 
			"higgs 125",
			h120,h130,
			120,130,125,
			norm, 1);
  h125->SetDirectory(&out);
  h125->Write("", TObject::kOverwrite);
    
//out.Write();
  out.Close();

  return 0;
}
