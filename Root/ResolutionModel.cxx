#include "HiggsllqqAnalysis/ResolutionModel.h"

ResolutionModel::ResolutionModel()
{
   m_templates_file = 0;
   m_truth_template_var = 0;
   m_truth_template_pdf = 0;
   m_higgs = 0;
   m_type = -1;
   m_targetZ = 0;
   SetConstraint(0, 0);
}

ResolutionModel::ResolutionModel(Double_t constraint_mass, Double_t constraint_width)
{
   m_templates_file = 0;
   m_truth_template_var = 0;
   m_truth_template_pdf = 0;
   SetConstraint(constraint_mass, constraint_width);
}

ResolutionModel::ResolutionModel(Double_t constraint_mass, Double_t constraint_width, Analysis::Quadrilepton *higgs, Int_t type, Int_t target)
{
   m_templates_file = 0;
   m_truth_template_var = 0;
   m_truth_template_pdf = 0;
   SetConstraint(constraint_mass, constraint_width);
   Set(higgs, type, target);
}

void ResolutionModel::SetConstraint(Double_t constraint_mass, Double_t constraint_width)
{
   m_conMass = constraint_mass;
   m_conWidth = constraint_width;
}

void ResolutionModel::Set(Analysis::Quadrilepton *higgs, Int_t type, Int_t target)
{
   m_higgs = higgs;
   m_type = type;

   if (target == ResolutionModel::UseZ1) {
      m_targetZ = m_higgs->GetZ1();
   } else if (target == ResolutionModel::UseZ2) {
      m_targetZ = m_higgs->GetZ2();
   } else {
      Error("Set", "Unrecognized target option %d, using Z1 instead", target);
      m_targetZ = m_higgs->GetZ1();
   }

   if (m_templates_file) {
      if (m_templates_workspace) {
         m_truth_template_var = m_templates_workspace->var("m2l_truth");
         m_truth_template_pdf = m_templates_workspace->pdf(TString::Format("Z_truth_mass_pdf_%d", GuessMass()));
      } else {
         Error("Set", "Invalid Z truth mass templates workspace");
      }
   } else {
      Error("Set", "Invalid Z truth mass templates file");
   }
}

void ResolutionModel::LoadTemplatesFromFile(TString filename)
{
   Info("LoadTemplatesFromFile", "Loading Z truth mass templates from %s", filename.Data());
   m_templates_file = new TFile(filename);

   if (m_templates_file) {
      m_templates_workspace = (RooWorkspace *)m_templates_file->Get("Z_truth_templates_workspace");

      if (!m_templates_workspace) {
         Error("LoadTemplatesFromFile", "Unable to load workspace \"Z_truth_templates_workspace\" from file %s", filename.Data());
      }
   } else {
      Error("LoadTemplatesFromFile", "Unable to load templates from file named %s", filename.Data());
   }
}

Int_t ResolutionModel::GuessMass()
{
   // must return generated Higgs sample mass closest to the reco value
   Double_t masses[14] = {110, 115, 120, 125, 130, 135, 140, 145, 150, 180, 260, 360, 460, 600};

   Double_t best_distance(99999999);
   Double_t closest(-1);

   for (Int_t i = 0; i < 14; i++) {
      Double_t this_distance = TMath::Abs(masses[i] - m_higgs->Get4Momentum()->M() / 1000);
      if (this_distance < best_distance) {
         best_distance = this_distance;
         closest = masses[i];
      }
   }

   return (Int_t)closest;
}

Double_t ResolutionModel::EvaluateLogTruthTemplate(Double_t mass)
{
   if (m_truth_template_var && m_truth_template_pdf) {
      m_truth_template_var->setVal(mass / 1000); // templates are given in GeV

      Double_t pdf_value = m_truth_template_pdf->getVal();

      if (pdf_value == 0) {
         return -9.99e09;
      } else {
         return TMath::Log(pdf_value);
      }
   } else {
      // default PDF return value is 1 if PDF has not been found (log(1) = 0)
      return 0;
   }
}

Double_t ResolutionModel::EvaluateLogCrystalBall(Double_t mass, Double_t mass_error)
{
   // adapted from RooCBShape::evaluate()

   Double_t result(0.);

   Double_t t = (m_targetZ->Get4Momentum()->M() - mass - m_cb_m0) / mass_error;

   if (m_cb_alpha < 0) t = -t;

   Double_t absAlpha = TMath::Abs(m_cb_alpha);

   if (t >= -absAlpha) {
      result = -0.5 * t * t;
   } else {
      Double_t a = TMath::Power(m_cb_n / absAlpha, m_cb_n) * TMath::Exp(-0.5 * absAlpha * absAlpha);
      Double_t b = m_cb_n / absAlpha - absAlpha;

      result = TMath::Log(a) - m_cb_n * TMath::Log(b - t);
   }

   return result;
}

Int_t ResolutionModel::GetCategory(Analysis::Dilepton *dilepton)
{
   Bool_t lepplus_mu_b = (TMath::Abs(dilepton->GetLepPlus()->Get4Momentum()->Eta()) < 1.05);
//Bool_t lepplus_mu_e = (TMath::Abs(dilepton->GetLepPlus()->Get4Momentum()->Eta()) > 1.05);
   Bool_t lepplus_e_b = (TMath::Abs(dilepton->GetLepPlus()->Get4Momentum()->Eta()) < 1.37);
   Bool_t lepplus_e_crk = (TMath::Abs(dilepton->GetLepPlus()->Get4Momentum()->Eta()) > 1.37 && TMath::Abs(dilepton->GetLepPlus()->Get4Momentum()->Eta()) < 1.52);
//Bool_t lepplus_e_e = (TMath::Abs(dilepton->GetLepPlus()->Get4Momentum()->Eta()) > 1.52 && TMath::Abs(dilepton->GetLepPlus()->Get4Momentum()->Eta()) < 2.47);
   Bool_t lepminus_mu_b = (TMath::Abs(dilepton->GetLepMinus()->Get4Momentum()->Eta()) < 1.05);
//Bool_t lepminus_mu_e = (TMath::Abs(dilepton->GetLepMinus()->Get4Momentum()->Eta()) > 1.05);
   Bool_t lepminus_e_b = (TMath::Abs(dilepton->GetLepMinus()->Get4Momentum()->Eta()) < 1.37);
   Bool_t lepminus_e_crk = (TMath::Abs(dilepton->GetLepMinus()->Get4Momentum()->Eta()) > 1.37 && TMath::Abs(dilepton->GetLepMinus()->Get4Momentum()->Eta()) < 1.52);
//Bool_t lepminus_e_e = (TMath::Abs(dilepton->GetLepMinus()->Get4Momentum()->Eta()) > 1.52 && TMath::Abs(dilepton->GetLepMinus()->Get4Momentum()->Eta()) < 2.47);

   if (dilepton->flavor() == Analysis::ChargedLepton::MUON) {
      if (lepplus_mu_b && lepminus_mu_b) {
         return DileptonCategory::bb;
      } else if (lepplus_mu_b || lepminus_mu_b) {
         return DileptonCategory::b;
      } else {
         return DileptonCategory::other;
      }
   } else if (dilepton->flavor() == Analysis::ChargedLepton::ELECTRON) {
      if (lepplus_e_b && lepminus_e_b) {
         return DileptonCategory::bb;
      } else if (lepplus_e_crk || lepminus_e_crk) {
         return DileptonCategory::onecrk;
      } else if (lepplus_e_b || lepminus_e_b) {
         return DileptonCategory::b;
      } else {
         return DileptonCategory::other;
      }
   }

   return -1;
}

void ResolutionModel::InitializeCrystalBall()
{
   Int_t samplemass = GuessMass();

   Int_t channel = m_targetZ->flavor();
   Int_t cat = GetCategory(m_targetZ);

   m_cb_m0 = -1;
   m_cb_alpha = -1;
   m_cb_n = -1;

   // paste here the output of z_studies.py
   if (channel == Analysis::ChargedLepton::MUON) {
      if (cat == DileptonCategory::all) {
         if (samplemass == 110) {
            m_cb_m0    = 83.338725;
            //m_cb_sigma = 1759.898992;
            //m_cb_fwhm  = 4140.000000;
            m_cb_alpha = 1.970435;
            m_cb_n     = 1.552867;
         } else if (samplemass == 115) {
            m_cb_m0    = -66.878152;
            //m_cb_sigma = 1788.849586;
            //m_cb_fwhm  = 4140.000000;
            m_cb_alpha = 1.626374;
            m_cb_n     = 2.326747;
         } else if (samplemass == 120) {
            m_cb_m0    = -112.893252;
            //m_cb_sigma = 1809.890344;
            //m_cb_fwhm  = 4200.000000;
            m_cb_alpha = 1.628421;
            m_cb_n     = 2.609117;
         } else if (samplemass == 125) {
            m_cb_m0    = -116.345674;
            //m_cb_sigma = 1788.539492;
            //m_cb_fwhm  = 4080.000000;
            m_cb_alpha = 1.284711;
            m_cb_n     = 4.260946;
         } else if (samplemass == 130) {
            m_cb_m0    = -164.440402;
            //m_cb_sigma = 1834.081843;
            //m_cb_fwhm  = 4140.000000;
            m_cb_alpha = 1.627942;
            m_cb_n     = 1.751872;
         } else if (samplemass == 135) {
            m_cb_m0    = -153.292015;
            //m_cb_sigma = 1856.787681;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.639563;
            m_cb_n     = 1.818906;
         } else if (samplemass == 140) {
            m_cb_m0    = -160.459055;
            //m_cb_sigma = 1853.217337;
            //m_cb_fwhm  = 4200.000000;
            m_cb_alpha = 1.752527;
            m_cb_n     = 1.618816;
         } else if (samplemass == 145) {
            m_cb_m0    = -148.471963;
            //m_cb_sigma = 1862.668060;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.560487;
            m_cb_n     = 2.155463;
         } else if (samplemass == 150) {
            m_cb_m0    = -190.692435;
            //m_cb_sigma = 1874.141639;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.969737;
            m_cb_n     = 0.891731;
         } else if (samplemass == 180) {
            m_cb_m0    = 71.837233;
            //m_cb_sigma = 1783.690781;
            //m_cb_fwhm  = 4080.000000;
            m_cb_alpha = 1.368032;
            m_cb_n     = 18.383316;
         } else if (samplemass == 260) {
            m_cb_m0    = -95.350945;
            //m_cb_sigma = 1699.071430;
            //m_cb_fwhm  = 3960.000000;
            m_cb_alpha = 1.711272;
            m_cb_n     = 4.182116;
         } else if (samplemass == 360) {
            m_cb_m0    = -1.151551;
            //m_cb_sigma = 1834.818698;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.636070;
            m_cb_n     = 6.552992;
         } else if (samplemass == 460) {
            m_cb_m0    = -37.588822;
            //m_cb_sigma = 2013.863161;
            //m_cb_fwhm  = 4620.000000;
            m_cb_alpha = 1.478040;
            m_cb_n     = 89.243479;
         } else if (samplemass == 600) {
            m_cb_m0    = -49.434391;
            //m_cb_sigma = 2136.633914;
            //m_cb_fwhm  = 4860.000000;
            m_cb_alpha = 1.455291;
            m_cb_n     = 39.839344;
         }
      } else if (cat == DileptonCategory::bb) {
         if (samplemass == 110) {
            m_cb_m0    = 16.942569;
            //m_cb_sigma = 1550.540745;
            //m_cb_fwhm  = 3600.000000;
            m_cb_alpha = 1.212301;
            m_cb_n     = 39.115022;
         } else if (samplemass == 115) {
            m_cb_m0    = -98.891293;
            //m_cb_sigma = 1472.500119;
            //m_cb_fwhm  = 3360.000000;
            m_cb_alpha = 1.359399;
            m_cb_n     = 2.845198;
         } else if (samplemass == 120) {
            m_cb_m0    = -76.910136;
            //m_cb_sigma = 1424.192741;
            //m_cb_fwhm  = 3240.000000;
            m_cb_alpha = 1.267381;
            m_cb_n     = 3.901223;
         } else if (samplemass == 125) {
            m_cb_m0    = -180.459320;
            //m_cb_sigma = 1592.350294;
            //m_cb_fwhm  = 3600.000000;
            m_cb_alpha = 1.264711;
            m_cb_n     = 4.073286;
         } else if (samplemass == 130) {
            m_cb_m0    = -226.442789;
            //m_cb_sigma = 1494.137706;
            //m_cb_fwhm  = 3420.000000;
            m_cb_alpha = 1.634537;
            m_cb_n     = 1.691438;
         } else if (samplemass == 135) {
            m_cb_m0    = -184.401061;
            //m_cb_sigma = 1550.857315;
            //m_cb_fwhm  = 3600.000000;
            m_cb_alpha = 1.625806;
            m_cb_n     = 1.662485;
         } else if (samplemass == 140) {
            m_cb_m0    = -220.002149;
            //m_cb_sigma = 1521.024122;
            //m_cb_fwhm  = 3480.000000;
            m_cb_alpha = 1.852330;
            m_cb_n     = 1.219912;
         } else if (samplemass == 145) {
            m_cb_m0    = -225.280241;
            //m_cb_sigma = 1461.081645;
            //m_cb_fwhm  = 3360.000000;
            m_cb_alpha = 1.923728;
            m_cb_n     = 0.846741;
         } else if (samplemass == 150) {
            m_cb_m0    = -234.821908;
            //m_cb_sigma = 1509.971471;
            //m_cb_fwhm  = 3420.000000;
            m_cb_alpha = 1.655878;
            m_cb_n     = 1.467191;
         } else if (samplemass == 180) {
            m_cb_m0    = -8.417070;
            //m_cb_sigma = 1440.636756;
            //m_cb_fwhm  = 3300.000000;
            m_cb_alpha = 1.372532;
            m_cb_n     = 4.802734;
         } else if (samplemass == 260) {
            m_cb_m0    = -117.287234;
            //m_cb_sigma = 1449.962941;
            //m_cb_fwhm  = 3300.000000;
            m_cb_alpha = 1.592180;
            m_cb_n     = 4.592896;
         } else if (samplemass == 360) {
            m_cb_m0    = -27.872905;
            //m_cb_sigma = 1664.662750;
            //m_cb_fwhm  = 3840.000000;
            m_cb_alpha = 1.658073;
            m_cb_n     = 5.282399;
         } else if (samplemass == 460) {
            m_cb_m0    = -57.501413;
            //m_cb_sigma = 1849.068350;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.868581;
            m_cb_n     = 3.141333;
         } else if (samplemass == 600) {
            m_cb_m0    = -141.218335;
            //m_cb_sigma = 1913.264612;
            //m_cb_fwhm  = 4440.000000;
            m_cb_alpha = 1.531097;
            m_cb_n     = 9.622578;
         }
      } else if (cat == DileptonCategory::b) {
         if (samplemass == 110) {
            m_cb_m0    = 74.738505;
            //m_cb_sigma = 1694.245020;
            //m_cb_fwhm  = 3900.000000;
            m_cb_alpha = 1.650854;
            m_cb_n     = 2.314361;
         } else if (samplemass == 115) {
            m_cb_m0    = -85.326428;
            //m_cb_sigma = 1790.303598;
            //m_cb_fwhm  = 4140.000000;
            m_cb_alpha = 1.609774;
            m_cb_n     = 2.664098;
         } else if (samplemass == 120) {
            m_cb_m0    = -86.514910;
            //m_cb_sigma = 1863.189753;
            //m_cb_fwhm  = 4320.000000;
            m_cb_alpha = 1.466595;
            m_cb_n     = 4.130477;
         } else if (samplemass == 125) {
            m_cb_m0    = -120.638151;
            //m_cb_sigma = 1783.575879;
            //m_cb_fwhm  = 4080.000000;
            m_cb_alpha = 1.423805;
            m_cb_n     = 2.603929;
         } else if (samplemass == 130) {
            m_cb_m0    = -106.183460;
            //m_cb_sigma = 1815.246983;
            //m_cb_fwhm  = 4200.000000;
            m_cb_alpha = 1.523763;
            m_cb_n     = 2.079119;
         } else if (samplemass == 135) {
            m_cb_m0    = -169.842015;
            //m_cb_sigma = 1876.370208;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.723570;
            m_cb_n     = 1.592297;
         } else if (samplemass == 140) {
            m_cb_m0    = -186.161801;
            //m_cb_sigma = 1861.698580;
            //m_cb_fwhm  = 4320.000000;
            m_cb_alpha = 1.653779;
            m_cb_n     = 1.880290;
         } else if (samplemass == 145) {
            m_cb_m0    = -120.807814;
            //m_cb_sigma = 1982.871930;
            //m_cb_fwhm  = 4620.000000;
            m_cb_alpha = 1.469981;
            m_cb_n     = 2.739196;
         } else if (samplemass == 150) {
            m_cb_m0    = -174.455947;
            //m_cb_sigma = 1924.953307;
            //m_cb_fwhm  = 4440.000000;
            m_cb_alpha = 1.889823;
            m_cb_n     = 1.148414;
         } else if (samplemass == 180) {
            m_cb_m0    = 81.599899;
            //m_cb_sigma = 1845.127333;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.794708;
            m_cb_n     = 2.798884;
         } else if (samplemass == 260) {
            m_cb_m0    = -79.534866;
            //m_cb_sigma = 1717.744102;
            //m_cb_fwhm  = 4020.000000;
            m_cb_alpha = 1.520879;
            m_cb_n     = 9.499437;
         } else if (samplemass == 360) {
            m_cb_m0    = -133.639276;
            //m_cb_sigma = 2048.086633;
            //m_cb_fwhm  = 4740.000000;
            m_cb_alpha = 5.000000;
            m_cb_n     = 52.057446;
         } else if (samplemass == 460) {
            m_cb_m0    = -77.652443;
            //m_cb_sigma = 2096.932191;
            //m_cb_fwhm  = 4800.000000;
            m_cb_alpha = 2.031856;
            m_cb_n     = 2.572873;
         } else if (samplemass == 600) {
            m_cb_m0    = -51.630492;
            //m_cb_sigma = 2180.346278;
            //m_cb_fwhm  = 5100.000000;
            m_cb_alpha = 1.417003;
            m_cb_n     = 106.220254;
         }
      } else if (cat == DileptonCategory::other) {
         if (samplemass == 110) {
            m_cb_m0    = 395.888158;
            //m_cb_sigma = 2002.693916;
            //m_cb_fwhm  = 4620.000000;
            m_cb_alpha = 1.443703;
            m_cb_n     = 46.926071;
         } else if (samplemass == 115) {
            m_cb_m0    = 174.619687;
            //m_cb_sigma = 2215.943087;
            //m_cb_fwhm  = 5100.000000;
            m_cb_alpha = 1.456674;
            m_cb_n     = 6.689293;
         } else if (samplemass == 120) {
            m_cb_m0    = -12.095312;
            //m_cb_sigma = 2205.700716;
            //m_cb_fwhm  = 5100.000000;
            m_cb_alpha = 1.575268;
            m_cb_n     = 4.214354;
         } else if (samplemass == 125) {
            m_cb_m0    = 86.528226;
            //m_cb_sigma = 2116.069400;
            //m_cb_fwhm  = 4920.000000;
            m_cb_alpha = 1.301408;
            m_cb_n     = 6.743359;
         } else if (samplemass == 130) {
            m_cb_m0    = -107.072696;
            //m_cb_sigma = 2424.909705;
            //m_cb_fwhm  = 5640.000000;
            m_cb_alpha = 1.563471;
            m_cb_n     = 3.288962;
         } else if (samplemass == 135) {
            m_cb_m0    = -56.425759;
            //m_cb_sigma = 2333.856463;
            //m_cb_fwhm  = 5400.000000;
            m_cb_alpha = 1.520145;
            m_cb_n     = 2.705118;
         } else if (samplemass == 140) {
            m_cb_m0    = 31.847129;
            //m_cb_sigma = 2412.487970;
            //m_cb_fwhm  = 5640.000000;
            m_cb_alpha = 1.558331;
            m_cb_n     = 3.848597;
         } else if (samplemass == 145) {
            m_cb_m0    = -189.808543;
            //m_cb_sigma = 2335.255961;
            //m_cb_fwhm  = 5400.000000;
            m_cb_alpha = 1.952059;
            m_cb_n     = 2.104190;
         } else if (samplemass == 150) {
            m_cb_m0    = -63.205663;
            //m_cb_sigma = 2304.332372;
            //m_cb_fwhm  = 5340.000000;
            m_cb_alpha = 1.945404;
            m_cb_n     = 1.049522;
         } else if (samplemass == 180) {
            m_cb_m0    = 210.797657;
            //m_cb_sigma = 2164.042306;
            //m_cb_fwhm  = 5040.000000;
            m_cb_alpha = 2.068704;
            m_cb_n     = 1.559099;
         } else if (samplemass == 260) {
            m_cb_m0    = 24.317970;
            //m_cb_sigma = 2134.393138;
            //m_cb_fwhm  = 4920.000000;
            m_cb_alpha = 1.440831;
            m_cb_n     = 80.224762;
         } else if (samplemass == 360) {
            m_cb_m0    = 73.424893;
            //m_cb_sigma = 2200.979515;
            //m_cb_fwhm  = 4980.000000;
            m_cb_alpha = 1.509409;
            m_cb_n     = 152.189416;
         } else if (samplemass == 460) {
            m_cb_m0    = -36.106789;
            //m_cb_sigma = 2355.814793;
            //m_cb_fwhm  = 5460.000000;
            m_cb_alpha = 2.293606;
            m_cb_n     = 1.221477;
         } else if (samplemass == 600) {
            m_cb_m0    = 242.528760;
            //m_cb_sigma = 2530.390575;
            //m_cb_fwhm  = 5820.000000;
            m_cb_alpha = 1.763969;
            m_cb_n     = 2.064810;
         }
      }
   } else if (channel == Analysis::ChargedLepton::ELECTRON) {
      if (cat == DileptonCategory::all) {
         if (samplemass == 110) {
            m_cb_m0    = -702.348413;
            //m_cb_sigma = 1867.603458;
            //m_cb_fwhm  = 4320.000000;
            m_cb_alpha = 1.072333;
            m_cb_n     = 5.528733;
         } else if (samplemass == 115) {
            m_cb_m0    = -713.327966;
            //m_cb_sigma = 1946.494183;
            //m_cb_fwhm  = 4620.000000;
            m_cb_alpha = 0.872525;
            m_cb_n     = 71.625782;
         } else if (samplemass == 120) {
            m_cb_m0    = -728.794519;
            //m_cb_sigma = 1943.517709;
            //m_cb_fwhm  = 4680.000000;
            m_cb_alpha = 0.785115;
            m_cb_n     = 131.392040;
         } else if (samplemass == 125) {
            m_cb_m0    = -779.454345;
            //m_cb_sigma = 1913.839867;
            //m_cb_fwhm  = 4560.000000;
            m_cb_alpha = 0.815711;
            m_cb_n     = 119.527862;
         } else if (samplemass == 130) {
            m_cb_m0    = -924.281673;
            //m_cb_sigma = 2095.974860;
            //m_cb_fwhm  = 4980.000000;
            m_cb_alpha = 0.823554;
            m_cb_n     = 120.093386;
         } else if (samplemass == 135) {
            m_cb_m0    = -780.632072;
            //m_cb_sigma = 1983.696793;
            //m_cb_fwhm  = 4860.000000;
            m_cb_alpha = 0.737977;
            m_cb_n     = 128.415592;
         } else if (samplemass == 140) {
            m_cb_m0    = -832.257376;
            //m_cb_sigma = 1978.865402;
            //m_cb_fwhm  = 4860.000000;
            m_cb_alpha = 0.737363;
            m_cb_n     = 63.189604;
         } else if (samplemass == 145) {
            m_cb_m0    = -776.884115;
            //m_cb_sigma = 2011.667221;
            //m_cb_fwhm  = 4800.000000;
            m_cb_alpha = 0.767526;
            m_cb_n     = 130.435580;
         } else if (samplemass == 150) {
            m_cb_m0    = -857.262406;
            //m_cb_sigma = 2051.083155;
            //m_cb_fwhm  = 4980.000000;
            m_cb_alpha = 0.773731;
            m_cb_n     = 133.309990;
         } else if (samplemass == 180) {
            m_cb_m0    = -677.995309;
            //m_cb_sigma = 1983.649941;
            //m_cb_fwhm  = 4620.000000;
            m_cb_alpha = 1.196490;
            m_cb_n     = 3.731928;
         } else if (samplemass == 260) {
            m_cb_m0    = -515.207911;
            //m_cb_sigma = 1571.390054;
            //m_cb_fwhm  = 3600.000000;
            m_cb_alpha = 1.133162;
            m_cb_n     = 7.663120;
         } else if (samplemass == 360) {
            m_cb_m0    = -399.020299;
            //m_cb_sigma = 1416.271348;
            //m_cb_fwhm  = 3240.000000;
            m_cb_alpha = 1.169063;
            m_cb_n     = 6.325116;
         } else if (samplemass == 460) {
            m_cb_m0    = -430.915415;
            //m_cb_sigma = 1441.240573;
            //m_cb_fwhm  = 3300.000000;
            m_cb_alpha = 1.598916;
            m_cb_n     = 2.883179;
         } else if (samplemass == 600) {
            m_cb_m0    = -342.868319;
            //m_cb_sigma = 1331.016807;
            //m_cb_fwhm  = 3060.000000;
            m_cb_alpha = 1.098937;
            m_cb_n     = 84.435220;
         }
      } else if (cat == DileptonCategory::bb) {
         if (samplemass == 110) {
            m_cb_m0    = -637.655929;
            //m_cb_sigma = 1619.482813;
            //m_cb_fwhm  = 3780.000000;
            m_cb_alpha = 1.003046;
            m_cb_n     = 12.497988;
         } else if (samplemass == 115) {
            m_cb_m0    = -509.906654;
            //m_cb_sigma = 1584.995748;
            //m_cb_fwhm  = 3720.000000;
            m_cb_alpha = 0.886149;
            m_cb_n     = 19.471449;
         } else if (samplemass == 120) {
            m_cb_m0    = -566.550432;
            //m_cb_sigma = 1610.697842;
            //m_cb_fwhm  = 3840.000000;
            m_cb_alpha = 0.819942;
            m_cb_n     = 19.503604;
         } else if (samplemass == 125) {
            m_cb_m0    = -638.211212;
            //m_cb_sigma = 1626.442622;
            //m_cb_fwhm  = 3840.000000;
            m_cb_alpha = 0.837080;
            m_cb_n     = 70.541974;
         } else if (samplemass == 130) {
            m_cb_m0    = -709.368752;
            //m_cb_sigma = 1702.112537;
            //m_cb_fwhm  = 4020.000000;
            m_cb_alpha = 0.960046;
            m_cb_n     = 6.430758;
         } else if (samplemass == 135) {
            m_cb_m0    = -685.189391;
            //m_cb_sigma = 1705.056796;
            //m_cb_fwhm  = 3900.000000;
            m_cb_alpha = 0.951198;
            m_cb_n     = 6.206968;
         } else if (samplemass == 140) {
            m_cb_m0    = -658.189590;
            //m_cb_sigma = 1647.162453;
            //m_cb_fwhm  = 3960.000000;
            m_cb_alpha = 0.803084;
            m_cb_n     = 17.540297;
         } else if (samplemass == 145) {
            m_cb_m0    = -700.045385;
            //m_cb_sigma = 1636.153685;
            //m_cb_fwhm  = 3780.000000;
            m_cb_alpha = 1.186473;
            m_cb_n     = 2.471956;
         } else if (samplemass == 150) {
            m_cb_m0    = -743.787063;
            //m_cb_sigma = 1700.347762;
            //m_cb_fwhm  = 4020.000000;
            m_cb_alpha = 1.067634;
            m_cb_n     = 3.759846;
         } else if (samplemass == 180) {
            m_cb_m0    = -590.489159;
            //m_cb_sigma = 1695.949622;
            //m_cb_fwhm  = 3960.000000;
            m_cb_alpha = 1.024067;
            m_cb_n     = 13.571414;
         } else if (samplemass == 260) {
            m_cb_m0    = -417.122919;
            //m_cb_sigma = 1364.407771;
            //m_cb_fwhm  = 3120.000000;
            m_cb_alpha = 0.994014;
            m_cb_n     = 50.267415;
         } else if (samplemass == 360) {
            m_cb_m0    = -315.190514;
            //m_cb_sigma = 1236.703116;
            //m_cb_fwhm  = 2880.000000;
            m_cb_alpha = 0.954489;
            m_cb_n     = 140.323582;
         } else if (samplemass == 460) {
            m_cb_m0    = -430.874117;
            //m_cb_sigma = 1275.602149;
            //m_cb_fwhm  = 2940.000000;
            m_cb_alpha = 1.503591;
            m_cb_n     = 4.996123;
         } else if (samplemass == 600) {
            m_cb_m0    = -297.405989;
            //m_cb_sigma = 1136.583425;
            //m_cb_fwhm  = 2580.000000;
            m_cb_alpha = 1.090059;
            m_cb_n     = 96.281093;
         }
      } else if (cat == DileptonCategory::onecrk) {
         if (samplemass == 110) {
            m_cb_m0    = 764.977078;
            //m_cb_sigma = 1747.720532;
            //m_cb_fwhm  = 5520.000000;
            m_cb_alpha = 0.392909;
            m_cb_n     = 57.046439;
         } else if (samplemass == 115) {
            m_cb_m0    = -1270.446240;
            //m_cb_sigma = 3455.304902;
            //m_cb_fwhm  = 7980.000000;
            m_cb_alpha = 1.205204;
            m_cb_n     = 116.903541;
         } else if (samplemass == 120) {
            m_cb_m0    = 229.746486;
            //m_cb_sigma = 2479.587459;
            //m_cb_fwhm  = 8760.000000;
            m_cb_alpha = 0.303257;
            m_cb_n     = 105.179535;
         } else if (samplemass == 125) {
            m_cb_m0    = -1603.166013;
            //m_cb_sigma = 2918.806061;
            //m_cb_fwhm  = 6780.000000;
            m_cb_alpha = 1.210782;
            m_cb_n     = 114.011645;
         } else if (samplemass == 130) {
            m_cb_m0    = -2185.722099;
            //m_cb_sigma = 3971.666793;
            //m_cb_fwhm  = 9240.000000;
            m_cb_alpha = 1.873804;
            m_cb_n     = 96.778767;
         } else if (samplemass == 135) {
            m_cb_m0    = -832.006340;
            //m_cb_sigma = 3106.851476;
            //m_cb_fwhm  = 8820.000000;
            m_cb_alpha = 0.468540;
            m_cb_n     = 97.260091;
         } else if (samplemass == 140) {
            m_cb_m0    = -1219.150757;
            //m_cb_sigma = 3275.266855;
            //m_cb_fwhm  = 8580.000000;
            m_cb_alpha = 0.562347;
            m_cb_n     = 97.147994;
         } else if (samplemass == 145) {
            m_cb_m0    = -2.996015;
            //m_cb_sigma = 3009.474276;
            //m_cb_fwhm  = 7560.000000;
            m_cb_alpha = 0.625909;
            m_cb_n     = 132.156800;
         } else if (samplemass == 150) {
            m_cb_m0    = -1233.639457;
            //m_cb_sigma = 3547.929216;
            //m_cb_fwhm  = 8880.000000;
            m_cb_alpha = 0.668655;
            m_cb_n     = 65.951988;
         } else if (samplemass == 180) {
            m_cb_m0    = -60.291521;
            //m_cb_sigma = 2640.620165;
            //m_cb_fwhm  = 6480.000000;
            m_cb_alpha = 0.694988;
            m_cb_n     = 122.191427;
         } else if (samplemass == 260) {
            m_cb_m0    = -575.991700;
            //m_cb_sigma = 2494.405731;
            //m_cb_fwhm  = 5820.000000;
            m_cb_alpha = 0.962721;
            m_cb_n     = 136.640028;
         } else if (samplemass == 360) {
            m_cb_m0    = -97.932005;
            //m_cb_sigma = 1855.722612;
            //m_cb_fwhm  = 4260.000000;
            m_cb_alpha = 1.083659;
            m_cb_n     = 2.567155;
         } else if (samplemass == 460) {
            m_cb_m0    = 157.236884;
            //m_cb_sigma = 1939.763156;
            //m_cb_fwhm  = 4500.000000;
            m_cb_alpha = 1.014653;
            m_cb_n     = 2.859143;
         } else if (samplemass == 600) {
            m_cb_m0    = -168.546387;
            //m_cb_sigma = 2145.970952;
            //m_cb_fwhm  = 4980.000000;
            m_cb_alpha = 1.752844;
            m_cb_n     = 2.016750;
         }
      } else if (cat == DileptonCategory::b) {
         if (samplemass == 110) {
            m_cb_m0    = -1021.695587;
            //m_cb_sigma = 2275.974345;
            //m_cb_fwhm  = 5280.000000;
            m_cb_alpha = 0.919558;
            m_cb_n     = 133.306015;
         } else if (samplemass == 115) {
            m_cb_m0    = -1001.931960;
            //m_cb_sigma = 2127.431845;
            //m_cb_fwhm  = 5100.000000;
            m_cb_alpha = 0.776887;
            m_cb_n     = 123.905461;
         } else if (samplemass == 120) {
            m_cb_m0    = -1220.365535;
            //m_cb_sigma = 2349.187741;
            //m_cb_fwhm  = 5520.000000;
            m_cb_alpha = 0.952396;
            m_cb_n     = 117.615135;
         } else if (samplemass == 125) {
            m_cb_m0    = -1366.571568;
            //m_cb_sigma = 2632.862905;
            //m_cb_fwhm  = 6180.000000;
            m_cb_alpha = 0.943420;
            m_cb_n     = 116.375901;
         } else if (samplemass == 130) {
            m_cb_m0    = -1449.382006;
            //m_cb_sigma = 2423.953825;
            //m_cb_fwhm  = 5760.000000;
            m_cb_alpha = 0.855607;
            m_cb_n     = 119.030924;
         } else if (samplemass == 135) {
            m_cb_m0    = -1147.571070;
            //m_cb_sigma = 2246.478282;
            //m_cb_fwhm  = 5580.000000;
            m_cb_alpha = 0.697392;
            m_cb_n     = 118.313700;
         } else if (samplemass == 140) {
            m_cb_m0    = -1167.294293;
            //m_cb_sigma = 2244.427936;
            //m_cb_fwhm  = 5580.000000;
            m_cb_alpha = 0.665529;
            m_cb_n     = 56.968437;
         } else if (samplemass == 145) {
            m_cb_m0    = -1168.811533;
            //m_cb_sigma = 2193.108282;
            //m_cb_fwhm  = 5880.000000;
            m_cb_alpha = 0.559249;
            m_cb_n     = 111.834467;
         } else if (samplemass == 150) {
            m_cb_m0    = -1147.923620;
            //m_cb_sigma = 2363.281705;
            //m_cb_fwhm  = 5580.000000;
            m_cb_alpha = 0.796061;
            m_cb_n     = 4.782191;
         } else if (samplemass == 180) {
            m_cb_m0    = -754.208769;
            //m_cb_sigma = 2152.296599;
            //m_cb_fwhm  = 5160.000000;
            m_cb_alpha = 0.808961;
            m_cb_n     = 136.076437;
         } else if (samplemass == 260) {
            m_cb_m0    = -642.113962;
            //m_cb_sigma = 1745.698012;
            //m_cb_fwhm  = 4080.000000;
            m_cb_alpha = 0.903868;
            m_cb_n     = 124.779661;
         } else if (samplemass == 360) {
            m_cb_m0    = -655.140659;
            //m_cb_sigma = 1760.744420;
            //m_cb_fwhm  = 4020.000000;
            m_cb_alpha = 1.254523;
            m_cb_n     = 4.384580;
         } else if (samplemass == 460) {
            m_cb_m0    = -586.062337;
            //m_cb_sigma = 1734.049505;
            //m_cb_fwhm  = 3960.000000;
            m_cb_alpha = 1.398785;
            m_cb_n     = 3.584200;
         } else if (samplemass == 600) {
            m_cb_m0    = -688.327866;
            //m_cb_sigma = 1712.249402;
            //m_cb_fwhm  = 3900.000000;
            m_cb_alpha = 1.088195;
            m_cb_n     = 21.960387;
         }
      } else if (cat == DileptonCategory::other) {
         if (samplemass == 110) {
            m_cb_m0    = -934.885201;
            //m_cb_sigma = 1874.034319;
            //m_cb_fwhm  = 4320.000000;
            m_cb_alpha = 1.071391;
            m_cb_n     = 103.252137;
         } else if (samplemass == 115) {
            m_cb_m0    = -1511.711728;
            //m_cb_sigma = 2537.691198;
            //m_cb_fwhm  = 5760.000000;
            m_cb_alpha = 1.133619;
            m_cb_n     = 117.738809;
         } else if (samplemass == 120) {
            m_cb_m0    = -727.078723;
            //m_cb_sigma = 1871.242076;
            //m_cb_fwhm  = 4800.000000;
            m_cb_alpha = 0.621026;
            m_cb_n     = 122.096246;
         } else if (samplemass == 125) {
            m_cb_m0    = -578.261289;
            //m_cb_sigma = 1624.621587;
            //m_cb_fwhm  = 4320.000000;
            m_cb_alpha = 0.546583;
            m_cb_n     = 16.553067;
         } else if (samplemass == 130) {
            m_cb_m0    = -1350.302744;
            //m_cb_sigma = 2829.393890;
            //m_cb_fwhm  = 6660.000000;
            m_cb_alpha = 0.912280;
            m_cb_n     = 116.670113;
         } else if (samplemass == 135) {
            m_cb_m0    = -1058.060971;
            //m_cb_sigma = 2087.654440;
            //m_cb_fwhm  = 5340.000000;
            m_cb_alpha = 0.596419;
            m_cb_n     = 99.047606;
         } else if (samplemass == 140) {
            m_cb_m0    = -1213.274174;
            //m_cb_sigma = 2248.109909;
            //m_cb_fwhm  = 5460.000000;
            m_cb_alpha = 0.784316;
            m_cb_n     = 136.657233;
         } else if (samplemass == 145) {
            m_cb_m0    = -759.484393;
            //m_cb_sigma = 2498.939792;
            //m_cb_fwhm  = 5820.000000;
            m_cb_alpha = 1.397706;
            m_cb_n     = 0.587953;
         } else if (samplemass == 150) {
            m_cb_m0    = -1266.168544;
            //m_cb_sigma = 2531.114442;
            //m_cb_fwhm  = 5940.000000;
            m_cb_alpha = 0.883790;
            m_cb_n     = 105.949468;
         } else if (samplemass == 180) {
            m_cb_m0    = -774.049653;
            //m_cb_sigma = 2448.083529;
            //m_cb_fwhm  = 5760.000000;
            m_cb_alpha = 1.147287;
            m_cb_n     = 2.373943;
         } else if (samplemass == 260) {
            m_cb_m0    = -919.269362;
            //m_cb_sigma = 1961.334549;
            //m_cb_fwhm  = 4560.000000;
            m_cb_alpha = 5.000000;
            m_cb_n     = 52.057446;
         } else if (samplemass == 360) {
            m_cb_m0    = -820.487975;
            //m_cb_sigma = 1790.719402;
            //m_cb_fwhm  = 4200.000000;
            m_cb_alpha = 1.038052;
            m_cb_n     = 35.127166;
         } else if (samplemass == 460) {
            m_cb_m0    = -795.398250;
            //m_cb_sigma = 2108.526843;
            //m_cb_fwhm  = 4860.000000;
            m_cb_alpha = 5.000000;
            m_cb_n     = 52.057446;
         } else if (samplemass == 600) {
            m_cb_m0    = -680.814154;
            //m_cb_sigma = 1974.865820;
            //m_cb_fwhm  = 4560.000000;
            m_cb_alpha = 2.205243;
            m_cb_n     = 0.626546;
         }
      }
   }

   if (m_cb_m0 == -1 && m_cb_alpha == -1 && m_cb_n == -1) {
      Error("InitializeCrystalBall", "samplemass = %d not recognized in channel %d category %d", samplemass, channel, cat);
   }
}

Double_t ResolutionModel::EvaluateLogPdf(Double_t mass, Double_t mass_error)
{
   Double_t result(-1);

   if (GetType() == ResolutionModel::BW_GAUS) {
      result = - TMath::Power(mass - m_targetZ->Get4Momentum()->M(), 2) / (2. * TMath::Power(mass_error, 2)) - TMath::Log(TMath::Power(mass * mass - m_conMass * m_conMass, 2) + TMath::Power(m_conMass * m_conWidth, 2));
   } else if (GetType() == ResolutionModel::BW_CRYBALL) {
      InitializeCrystalBall();

      result = EvaluateLogCrystalBall(mass, mass_error) - TMath::Log(TMath::Power(mass * mass - m_conMass * m_conMass, 2) + TMath::Power(m_conMass * m_conWidth, 2));
   } else if (GetType() == ResolutionModel::TEMPLATE_CRYBALL) {
      InitializeCrystalBall();

      result = EvaluateLogCrystalBall(mass, mass_error) + EvaluateLogTruthTemplate(mass);
   } else if (GetType() == ResolutionModel::TEMPLATE_GAUS) {
      result = - TMath::Power(mass - m_targetZ->Get4Momentum()->M(), 2) / (2. * TMath::Power(mass_error, 2)) + EvaluateLogTruthTemplate(mass);
   }

   return result;
}

ResolutionModel::~ResolutionModel()
{
   // does not handle pointer deletion! this must be done by the code using this class
   if (m_templates_file) {
      m_templates_file->Close();
   }
}

void ResolutionModel::CopyFrom(ResolutionModel *theSource)
{
}
