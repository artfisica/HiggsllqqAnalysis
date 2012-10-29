#include <TString.h>
#include <TH1F.h>
#include <THStack.h>
#include <TFile.h>
#include <TCanvas.h>
#include <map>
#include <iostream>

void stampa_canale(TFile *file, int canale, int factor)
{
  std::map<TString, Float_t> datadriven_yield[4];

  datadriven_yield[0]["histoZ"] = 0.;
  datadriven_yield[0]["histoZbb"] = 0.149;
  datadriven_yield[0]["histott"] = 0.017;
  datadriven_yield[3]["histoZ"] = 0.;
  datadriven_yield[3]["histoZbb"] = 0.119;
  datadriven_yield[3]["histott"] = 0.015;
  datadriven_yield[1]["histoZ"] = 3.4;
  datadriven_yield[1]["histoZbb"] = 0.;
  datadriven_yield[1]["histott"] = 0.;
  datadriven_yield[2]["histoZ"] = 1.3;
  datadriven_yield[2]["histoZbb"] = 0.;
  datadriven_yield[2]["histott"] = 0.;

  

  file->cd(); 

  TH1F *ZZ = (TH1F*)file->Get(TString::Format("histoZZ_%d", canale));
  TH1F *gg2ZZ = (TH1F*)file->Get(TString::Format("histogg2ZZ_%d", canale));
  TH1F *Z = (TH1F*)file->Get(TString::Format("histoZ_%d", canale));
  TH1F *Zbb = (TH1F*)file->Get(TString::Format("histoZbb_%d", canale));
  TH1F *tt = (TH1F*)file->Get(TString::Format("histott_%d", canale));
  TH1F *DATA = (TH1F*)file->Get(TString::Format("histoDATA_%d", canale));
  TH1F *higgs = (TH1F*)file->Get(TString::Format("higgs_%d", canale));
  TH1F *higgsVBF = (TH1F*)file->Get(TString::Format("higgsVBF_%d", canale));
  TH1F *higgsWH = (TH1F*)file->Get(TString::Format("higgsWH_%d", canale));
  TH1F *higgsZH = (TH1F*)file->Get(TString::Format("higgsZH_%d", canale));

  if (Z->Integral(0, Z->GetNbinsX() + 1) != 0.) Z->Scale(datadriven_yield[canale]["histoZ"] / Z->Integral(0, Z->GetNbinsX() + 1));
  if (Zbb->Integral(0, Zbb->GetNbinsX() + 1) != 0.) Zbb->Scale(datadriven_yield[canale]["histoZbb"] / Zbb->Integral(0, Zbb->GetNbinsX() + 1));
  if (tt->Integral(0, tt->GetNbinsX() + 1) != 0.) tt->Scale(datadriven_yield[canale]["histott"] / tt->Integral(0, tt->GetNbinsX() + 1));

  ZZ->SetFillColor(kRed);
  gg2ZZ->SetFillColor(kViolet);
  Z->SetFillColor(kBlue+1);
  Zbb->SetFillColor(kGreen-5);
  tt->SetFillColor(kOrange);
  higgs->SetFillColor(kCyan);
  higgsVBF->SetFillColor(kCyan);
  higgsWH->SetFillColor(kCyan);
  higgsZH->SetFillColor(kCyan);

  DATA->Rebin(factor);
  ZZ->Rebin(factor);
  gg2ZZ->Rebin(factor);
  Z->Rebin(factor);
  Zbb->Rebin(factor);
  tt->Rebin(factor);
  higgs->Rebin(factor);
  higgsVBF->Rebin(factor);
  higgsWH->Rebin(factor);
  higgsZH->Rebin(factor);


  THStack *fondi = new THStack();
  fondi->Add(ZZ);
  fondi->Add(gg2ZZ);
  fondi->Add(Z);
  fondi->Add(Zbb);
  fondi->Add(tt);

  fondi->Add(higgs);
  fondi->Add(higgsVBF);
  fondi->Add(higgsWH);
  fondi->Add(higgsZH);

  TString chan_name[4]; TString chan_label[4];
  chan_name[0] = "4mu"; chan_label[0] = "m_{4#mu} [GeV]";
  chan_name[1] = "2mu2e"; chan_label[1] = "m_{2#mu2e} [GeV]";
  chan_name[3] = "2e2mu"; chan_label[3] = "m_{2e2#mu} [GeV]";
  chan_name[2] = "4e"; chan_label[2] = "m_{4e} [GeV]";

  TCanvas *c = new TCanvas(TString::Format("c_%d", canale), chan_name[canale].Data());

  DATA->SetMarkerStyle(20);
  DATA->GetXaxis()->SetTitle(chan_label[canale].Data());
  DATA->GetYaxis()->SetTitle(TString::Format("Events / %.1lf GeV", 0.5 * factor));

  DATA->Draw("pe");
  fondi->Draw("same");
  DATA->Draw("pesame");


  std::cout << "channel (kostas' convention) number " << canale << std::endl;
  std::cout << "ZZ: integral       = " << ZZ->Integral(0, ZZ->GetNbinsX() + 1) << std::endl;
  std::cout << "gg2ZZ: integral    = " << gg2ZZ->Integral(0, gg2ZZ->GetNbinsX() + 1) << std::endl;
  std::cout << "Z: integral        = " << Z->Integral(0, Z->GetNbinsX() + 1) << std::endl;
  std::cout << "Zbb: integral      = " << Zbb->Integral(0, Zbb->GetNbinsX() + 1) << std::endl;
  std::cout << "tt: integral       = " << tt->Integral(0, tt->GetNbinsX() + 1) << std::endl;
  std::cout << "DATA: integral     = " << DATA->Integral(0, DATA->GetNbinsX() + 1) << std::endl;
  std::cout << "higgs: integral    = " << higgs->Integral(0, higgs->GetNbinsX() + 1) << std::endl;
  std::cout << "higgsVBF: integral = " << higgsVBF->Integral(0, higgsVBF->GetNbinsX() + 1) << std::endl;
  std::cout << "higgsWH: integral  = " << higgsWH->Integral(0, higgsWH->GetNbinsX() + 1) << std::endl;
  std::cout << "higgsZH: integral  = " << higgsZH->Integral(0, higgsZH->GetNbinsX() + 1) << std::endl;

  c->SaveAs(TString::Format("%s.eps", chan_name[canale].Data()));
  DATA->GetXaxis()->SetRangeUser(100, 250);
  c->SaveAs(TString::Format("%s_zoom.eps", chan_name[canale].Data()));
}

void plot_mass() {
  TFile *file = new TFile("workspace_histos/ahistos_1250.root");
  if (!file) {
    std::cout << "file not found" << std::endl;
    return;
  }
  for (int i = 0; i < 4; i++) 
  //stampa_canale(file, i, 10);
    stampa_canale(file, i, 50);
}
