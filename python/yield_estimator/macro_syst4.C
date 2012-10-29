{
ggH130_4e_syst_nominal->Draw();
ggH130_4e_syst_ES_Z_up->Draw("same");
ggH130_4e_syst_ES_Z_down->Draw("same");
ggH130_4e_syst_ES_Z_up->SetLineColor(kRed);
ggH130_4e_syst_ES_Z_down->SetLineColor(kBlue+1);
ggH130_4e_syst_nominal->GetXaxis()->SetTitle("m_{eeee} [GeV]");
ggH130_4e_syst_nominal->GetYaxis()->SetTitle("a.u. / 0.5 GeV");
ggH130_4e_syst_nominal->GetBinWidth(3);
ggH130_4e_syst_nominal->GetXaxis()->SetRangeUser(80, 150);
ggH130_4e_syst_nominal->SetMarkerSize(0);
ggH130_4e_syst_ES_Z_up->SetMarkerSize(0);
ggH130_4e_syst_ES_Z_down->SetMarkerSize(0);
ggH130_4e_syst_ES_Z_down->SetTitle("ES_Z_down");
ggH130_4e_syst_ES_Z_up->SetTitle("ES_Z_up");
ggH130_4e_syst_nominal->SetTitle("nominal");
leg = c1->BuildLegend();
leg->SetFillColor(0);
leg->SetBorderSize(0);
}
