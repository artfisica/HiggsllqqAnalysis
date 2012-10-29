{

// brings histos on scope

h_4mu_m12_minus_truth_over_truth = h_4mu_m12_minus_truth_over_truth;
h_4mu_m12_constrained_minus_truth_over_truth = h_4mu_m12_constrained_minus_truth_over_truth;
h_4mu_m4l_minus_truth_over_truth = h_4mu_m4l_minus_truth_over_truth;
h_4mu_m4l_constrained_minus_truth_over_truth = h_4mu_m4l_constrained_minus_truth_over_truth;
h_4mu_lep12_pt_minus_truth_over_truth = h_4mu_lep12_pt_minus_truth_over_truth;
h_4mu_lep12_pt_constrained_minus_truth_over_truth = h_4mu_lep12_pt_constrained_minus_truth_over_truth;
h_4mu_m12_constrained_vs_m12 = h_4mu_m12_constrained_vs_m12;
h_4mu_m4l_constrained_vs_m4l = h_4mu_m4l_constrained_vs_m4l;
h_2mu2e_m12_minus_truth_over_truth = h_2mu2e_m12_minus_truth_over_truth;
h_2mu2e_m12_constrained_minus_truth_over_truth = h_2mu2e_m12_constrained_minus_truth_over_truth;
h_2mu2e_m4l_minus_truth_over_truth = h_2mu2e_m4l_minus_truth_over_truth;
h_2mu2e_m4l_constrained_minus_truth_over_truth = h_2mu2e_m4l_constrained_minus_truth_over_truth;
h_2mu2e_lep12_pt_minus_truth_over_truth = h_2mu2e_lep12_pt_minus_truth_over_truth;
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth = h_2mu2e_lep12_pt_constrained_minus_truth_over_truth;
h_2mu2e_m12_constrained_vs_m12 = h_2mu2e_m12_constrained_vs_m12;
h_2mu2e_m4l_constrained_vs_m4l = h_2mu2e_m4l_constrained_vs_m4l;
h_2e2mu_m12_minus_truth_over_truth = h_2e2mu_m12_minus_truth_over_truth;
h_2e2mu_m12_constrained_minus_truth_over_truth = h_2e2mu_m12_constrained_minus_truth_over_truth;
h_2e2mu_m4l_minus_truth_over_truth = h_2e2mu_m4l_minus_truth_over_truth;
h_2e2mu_m4l_constrained_minus_truth_over_truth = h_2e2mu_m4l_constrained_minus_truth_over_truth;
h_2e2mu_lep12_pt_minus_truth_over_truth = h_2e2mu_lep12_pt_minus_truth_over_truth;
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth = h_2e2mu_lep12_pt_constrained_minus_truth_over_truth;
h_2e2mu_m12_constrained_vs_m12 = h_2e2mu_m12_constrained_vs_m12;
h_2e2mu_m4l_constrained_vs_m4l = h_2e2mu_m4l_constrained_vs_m4l;
h_4e_m12_minus_truth_over_truth = h_4e_m12_minus_truth_over_truth;
h_4e_m12_constrained_minus_truth_over_truth = h_4e_m12_constrained_minus_truth_over_truth;
h_4e_m4l_minus_truth_over_truth = h_4e_m4l_minus_truth_over_truth;
h_4e_m4l_constrained_minus_truth_over_truth = h_4e_m4l_constrained_minus_truth_over_truth;
h_4e_lep12_pt_minus_truth_over_truth = h_4e_lep12_pt_minus_truth_over_truth;
h_4e_lep12_pt_constrained_minus_truth_over_truth = h_4e_lep12_pt_constrained_minus_truth_over_truth;
h_4e_m12_constrained_vs_m12 = h_4e_m12_constrained_vs_m12;
h_4e_m4l_constrained_vs_m4l = h_4e_m4l_constrained_vs_m4l;


// 4mu

h_4mu_m12_minus_truth_over_truth->Rebin(5);
h_4mu_m12_constrained_minus_truth_over_truth->Rebin(5);
h_4mu_m4l_minus_truth_over_truth->Rebin(5);
h_4mu_m4l_constrained_minus_truth_over_truth->Rebin(5);
h_4mu_lep12_pt_minus_truth_over_truth->Rebin(5);
h_4mu_lep12_pt_constrained_minus_truth_over_truth->Rebin(5);
h_4mu_m12_minus_truth_over_truth->SetMarkerSize(0);
h_4mu_m12_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_4mu_m4l_minus_truth_over_truth->SetMarkerSize(0);
h_4mu_m4l_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_4mu_lep12_pt_minus_truth_over_truth->SetMarkerSize(0);
h_4mu_lep12_pt_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_4mu_m12_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4mu_m12_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4mu_m4l_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4mu_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4mu_lep12_pt_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4mu_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4mu_m12_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4mu_m12_minus_truth_over_truth->GetMaximum() * 1.2);
h_4mu_m12_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4mu_m12_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_4mu_m4l_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4mu_m4l_minus_truth_over_truth->GetMaximum() * 1.2);
h_4mu_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4mu_m4l_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_4mu_lep12_pt_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4mu_lep12_pt_minus_truth_over_truth->GetMaximum() * 1.2);
h_4mu_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4mu_lep12_pt_constrained_minus_truth_over_truth->GetMaximum() * 1.2);

TCanvas *c_4mu[5];
TLegend *leg_4mu[5];

c_4mu[0] = new TCanvas("c_4mu_m12_pull");
c_4mu[0]->cd();
h_4mu_m12_minus_truth_over_truth->SetTitle("no constraint");
h_4mu_m12_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_4mu_m12_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_4mu_m12_minus_truth_over_truth->Draw();
h_4mu_m12_constrained_minus_truth_over_truth->Draw("same");
leg_4mu[0] = c_4mu[0]->BuildLegend(0.62, 0.79, 1, 1);
leg_4mu[0]->SetFillColor(0);
leg_4mu[0]->SetBorderSize(0);

c_4mu[1] = new TCanvas("c_4mu_m4l_pull");
c_4mu[1]->cd();
h_4mu_m4l_minus_truth_over_truth->SetTitle("no constraint");
h_4mu_m4l_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_4mu_m4l_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_4mu_m4l_minus_truth_over_truth->Draw();
h_4mu_m4l_constrained_minus_truth_over_truth->Draw("same");
leg_4mu[1] = c_4mu[1]->BuildLegend(0.62, 0.79, 1, 1);
leg_4mu[1]->SetFillColor(0);
leg_4mu[1]->SetBorderSize(0);

c_4mu[2] = new TCanvas("c_4mu_lep12_pt_pull");
c_4mu[2]->cd();
h_4mu_lep12_pt_minus_truth_over_truth->SetTitle("no constraint");
h_4mu_lep12_pt_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_4mu_lep12_pt_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_4mu_lep12_pt_minus_truth_over_truth->Draw();
h_4mu_lep12_pt_constrained_minus_truth_over_truth->Draw("same");
leg_4mu[2] = c_4mu[2]->BuildLegend(0.62, 0.79, 1, 1);
leg_4mu[2]->SetFillColor(0);
leg_4mu[2]->SetBorderSize(0);

c_4mu[3] = new TCanvas("c_4mu_m12_constr_vs_normal");
c_4mu[3]->cd();
h_4mu_m12_constrained_vs_m12->Draw("colz");

c_4mu[4] = new TCanvas("c_4mu_m4l_constr_vs_normal");
c_4mu[4]->cd();
h_4mu_m4l_constrained_vs_m4l->Draw("colz");

for (int i = 0; i < 5; i++) c_4mu[i].SaveAs(TString::Format("%s.eps", c_4mu[i].GetName()));

// 2mu2e

h_2mu2e_m12_minus_truth_over_truth->Rebin(5);
h_2mu2e_m12_constrained_minus_truth_over_truth->Rebin(5);
h_2mu2e_m4l_minus_truth_over_truth->Rebin(5);
h_2mu2e_m4l_constrained_minus_truth_over_truth->Rebin(5);
h_2mu2e_lep12_pt_minus_truth_over_truth->Rebin(5);
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->Rebin(5);
h_2mu2e_m12_minus_truth_over_truth->SetMarkerSize(0);
h_2mu2e_m12_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_2mu2e_m4l_minus_truth_over_truth->SetMarkerSize(0);
h_2mu2e_m4l_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_2mu2e_lep12_pt_minus_truth_over_truth->SetMarkerSize(0);
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_2mu2e_m12_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2mu2e_m12_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2mu2e_m4l_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2mu2e_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2mu2e_lep12_pt_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2mu2e_m12_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2mu2e_m12_minus_truth_over_truth->GetMaximum() * 1.2);
h_2mu2e_m12_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2mu2e_m12_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_2mu2e_m4l_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2mu2e_m4l_minus_truth_over_truth->GetMaximum() * 1.2);
h_2mu2e_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2mu2e_m4l_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_2mu2e_lep12_pt_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2mu2e_lep12_pt_minus_truth_over_truth->GetMaximum() * 1.2);
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->GetMaximum() * 1.2);

TCanvas *c_2mu2e[5];
TLegend *leg_2mu2e[5];

c_2mu2e[0] = new TCanvas("c_2mu2e_m12_pull");
c_2mu2e[0]->cd();
h_2mu2e_m12_minus_truth_over_truth->SetTitle("no constraint");
h_2mu2e_m12_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_2mu2e_m12_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_2mu2e_m12_minus_truth_over_truth->Draw();
h_2mu2e_m12_constrained_minus_truth_over_truth->Draw("same");
leg_2mu2e[0] = c_2mu2e[0]->BuildLegend(0.62, 0.79, 1, 1);
leg_2mu2e[0]->SetFillColor(0);
leg_2mu2e[0]->SetBorderSize(0);

c_2mu2e[1] = new TCanvas("c_2mu2e_m4l_pull");
c_2mu2e[1]->cd();
h_2mu2e_m4l_minus_truth_over_truth->SetTitle("no constraint");
h_2mu2e_m4l_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_2mu2e_m4l_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_2mu2e_m4l_minus_truth_over_truth->Draw();
h_2mu2e_m4l_constrained_minus_truth_over_truth->Draw("same");
leg_2mu2e[1] = c_2mu2e[1]->BuildLegend(0.62, 0.79, 1, 1);
leg_2mu2e[1]->SetFillColor(0);
leg_2mu2e[1]->SetBorderSize(0);

c_2mu2e[2] = new TCanvas("c_2mu2e_lep12_pt_pull");
c_2mu2e[2]->cd();
h_2mu2e_lep12_pt_minus_truth_over_truth->SetTitle("no constraint");
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_2mu2e_lep12_pt_minus_truth_over_truth->Draw();
h_2mu2e_lep12_pt_constrained_minus_truth_over_truth->Draw("same");
leg_2mu2e[2] = c_2mu2e[2]->BuildLegend(0.62, 0.79, 1, 1);
leg_2mu2e[2]->SetFillColor(0);
leg_2mu2e[2]->SetBorderSize(0);

c_2mu2e[3] = new TCanvas("c_2mu2e_m12_constr_vs_normal");
c_2mu2e[3]->cd();
h_2mu2e_m12_constrained_vs_m12->Draw("colz");

c_2mu2e[4] = new TCanvas("c_2mu2e_m4l_constr_vs_normal");
c_2mu2e[4]->cd();
h_2mu2e_m4l_constrained_vs_m4l->Draw("colz");

for (int i = 0; i < 5; i++) c_2mu2e[i].SaveAs(TString::Format("%s.eps", c_2mu2e[i].GetName()));

// 2e2mu

h_2e2mu_m12_minus_truth_over_truth->Rebin(5);
h_2e2mu_m12_constrained_minus_truth_over_truth->Rebin(5);
h_2e2mu_m4l_minus_truth_over_truth->Rebin(5);
h_2e2mu_m4l_constrained_minus_truth_over_truth->Rebin(5);
h_2e2mu_lep12_pt_minus_truth_over_truth->Rebin(5);
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->Rebin(5);
h_2e2mu_m12_minus_truth_over_truth->SetMarkerSize(0);
h_2e2mu_m12_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_2e2mu_m4l_minus_truth_over_truth->SetMarkerSize(0);
h_2e2mu_m4l_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_2e2mu_lep12_pt_minus_truth_over_truth->SetMarkerSize(0);
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_2e2mu_m12_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2e2mu_m12_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2e2mu_m4l_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2e2mu_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2e2mu_lep12_pt_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_2e2mu_m12_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2e2mu_m12_minus_truth_over_truth->GetMaximum() * 1.2);
h_2e2mu_m12_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2e2mu_m12_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_2e2mu_m4l_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2e2mu_m4l_minus_truth_over_truth->GetMaximum() * 1.2);
h_2e2mu_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2e2mu_m4l_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_2e2mu_lep12_pt_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2e2mu_lep12_pt_minus_truth_over_truth->GetMaximum() * 1.2);
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->GetMaximum() * 1.2);

TCanvas *c_2e2mu[5];
TLegend *leg_2e2mu[5];

c_2e2mu[0] = new TCanvas("c_2e2mu_m12_pull");
c_2e2mu[0]->cd();
h_2e2mu_m12_minus_truth_over_truth->SetTitle("no constraint");
h_2e2mu_m12_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_2e2mu_m12_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_2e2mu_m12_minus_truth_over_truth->Draw();
h_2e2mu_m12_constrained_minus_truth_over_truth->Draw("same");
leg_2e2mu[0] = c_2e2mu[0]->BuildLegend(0.62, 0.79, 1, 1);
leg_2e2mu[0]->SetFillColor(0);
leg_2e2mu[0]->SetBorderSize(0);

c_2e2mu[1] = new TCanvas("c_2e2mu_m4l_pull");
c_2e2mu[1]->cd();
h_2e2mu_m4l_minus_truth_over_truth->SetTitle("no constraint");
h_2e2mu_m4l_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_2e2mu_m4l_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_2e2mu_m4l_minus_truth_over_truth->Draw();
h_2e2mu_m4l_constrained_minus_truth_over_truth->Draw("same");
leg_2e2mu[1] = c_2e2mu[1]->BuildLegend(0.62, 0.79, 1, 1);
leg_2e2mu[1]->SetFillColor(0);
leg_2e2mu[1]->SetBorderSize(0);

c_2e2mu[2] = new TCanvas("c_2e2mu_lep12_pt_pull");
c_2e2mu[2]->cd();
h_2e2mu_lep12_pt_minus_truth_over_truth->SetTitle("no constraint");
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_2e2mu_lep12_pt_minus_truth_over_truth->Draw();
h_2e2mu_lep12_pt_constrained_minus_truth_over_truth->Draw("same");
leg_2e2mu[2] = c_2e2mu[2]->BuildLegend(0.62, 0.79, 1, 1);
leg_2e2mu[2]->SetFillColor(0);
leg_2e2mu[2]->SetBorderSize(0);

c_2e2mu[3] = new TCanvas("c_2e2mu_m12_constr_vs_normal");
c_2e2mu[3]->cd();
h_2e2mu_m12_constrained_vs_m12->Draw("colz");

c_2e2mu[4] = new TCanvas("c_2e2mu_m4l_constr_vs_normal");
c_2e2mu[4]->cd();
h_2e2mu_m4l_constrained_vs_m4l->Draw("colz");

for (int i = 0; i < 5; i++) c_2e2mu[i].SaveAs(TString::Format("%s.eps", c_2e2mu[i].GetName()));

// 4e

h_4e_m12_minus_truth_over_truth->Rebin(5);
h_4e_m12_constrained_minus_truth_over_truth->Rebin(5);
h_4e_m4l_minus_truth_over_truth->Rebin(5);
h_4e_m4l_constrained_minus_truth_over_truth->Rebin(5);
h_4e_lep12_pt_minus_truth_over_truth->Rebin(5);
h_4e_lep12_pt_constrained_minus_truth_over_truth->Rebin(5);
h_4e_m12_minus_truth_over_truth->SetMarkerSize(0);
h_4e_m12_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_4e_m4l_minus_truth_over_truth->SetMarkerSize(0);
h_4e_m4l_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_4e_lep12_pt_minus_truth_over_truth->SetMarkerSize(0);
h_4e_lep12_pt_constrained_minus_truth_over_truth->SetMarkerSize(0);
h_4e_m12_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4e_m12_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4e_m4l_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4e_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4e_lep12_pt_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4e_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetTitle("a.u. / 0.01");
h_4e_m12_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4e_m12_minus_truth_over_truth->GetMaximum() * 1.2);
h_4e_m12_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4e_m12_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_4e_m4l_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4e_m4l_minus_truth_over_truth->GetMaximum() * 1.2);
h_4e_m4l_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4e_m4l_constrained_minus_truth_over_truth->GetMaximum() * 1.2);
h_4e_lep12_pt_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4e_lep12_pt_minus_truth_over_truth->GetMaximum() * 1.2);
h_4e_lep12_pt_constrained_minus_truth_over_truth->GetYaxis()->SetRangeUser(0, h_4e_lep12_pt_constrained_minus_truth_over_truth->GetMaximum() * 1.2);

TCanvas *c_4e[5];
TLegend *leg_4e[5];

c_4e[0] = new TCanvas("c_4e_m12_pull");
c_4e[0]->cd();
h_4e_m12_minus_truth_over_truth->SetTitle("no constraint");
h_4e_m12_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_4e_m12_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_4e_m12_minus_truth_over_truth->Draw();
h_4e_m12_constrained_minus_truth_over_truth->Draw("same");
leg_4e[0] = c_4e[0]->BuildLegend(0.62, 0.79, 1, 1);
leg_4e[0]->SetFillColor(0);
leg_4e[0]->SetBorderSize(0);

c_4e[1] = new TCanvas("c_4e_m4l_pull");
c_4e[1]->cd();
h_4e_m4l_minus_truth_over_truth->SetTitle("no constraint");
h_4e_m4l_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_4e_m4l_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_4e_m4l_minus_truth_over_truth->Draw();
h_4e_m4l_constrained_minus_truth_over_truth->Draw("same");
leg_4e[1] = c_4e[1]->BuildLegend(0.62, 0.79, 1, 1);
leg_4e[1]->SetFillColor(0);
leg_4e[1]->SetBorderSize(0);

c_4e[2] = new TCanvas("c_4e_lep12_pt_pull");
c_4e[2]->cd();
h_4e_lep12_pt_minus_truth_over_truth->SetTitle("no constraint");
h_4e_lep12_pt_constrained_minus_truth_over_truth->SetTitle("with constraint");
h_4e_lep12_pt_constrained_minus_truth_over_truth->SetLineColor(kRed);
h_4e_lep12_pt_minus_truth_over_truth->Draw();
h_4e_lep12_pt_constrained_minus_truth_over_truth->Draw("same");
leg_4e[2] = c_4e[2]->BuildLegend(0.62, 0.79, 1, 1);
leg_4e[2]->SetFillColor(0);
leg_4e[2]->SetBorderSize(0);

c_4e[3] = new TCanvas("c_4e_m12_constr_vs_normal");
c_4e[3]->cd();
h_4e_m12_constrained_vs_m12->Draw("colz");

c_4e[4] = new TCanvas("c_4e_m4l_constr_vs_normal");
c_4e[4]->cd();
h_4e_m4l_constrained_vs_m4l->Draw("colz");

for (int i = 0; i < 5; i++) c_4e[i].SaveAs(TString::Format("%s.eps", c_4e[i].GetName()));

}
