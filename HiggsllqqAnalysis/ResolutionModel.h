#ifndef __ResolutionModel__
#define __ResolutionModel__

#include <iostream>
#include "HiggsllqqAnalysis/Quadrilepton.h"
#include <TMath.h>
#include <TFile.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>

class ResolutionModel {
private:
   Analysis::Quadrilepton *m_higgs;
   Analysis::Dilepton *m_targetZ;
   Int_t m_type;

   Double_t m_conMass;
   Double_t m_conWidth;
   Double_t m_cb_m0;
   Double_t m_cb_alpha;
   Double_t m_cb_n;

   TFile *m_templates_file;
   RooWorkspace *m_templates_workspace;

   RooRealVar *m_truth_template_var;
   RooAbsPdf *m_truth_template_pdf;


public:
   enum {
      BW_GAUS,
      BW_CRYBALL,
      TEMPLATE_GAUS,
      TEMPLATE_CRYBALL,
   };

   enum {
      UseZ1,
      UseZ2,
   };

   void Set(Analysis::Quadrilepton *higgs, Int_t type, Int_t target = ResolutionModel::UseZ1);
   void SetConstraint(Double_t constraint_mass, Double_t constraint_width);
   void LoadTemplatesFromFile(TString filename);
   Int_t GuessMass();
   Int_t GetType() {
      return m_type;
   }
   Int_t GetCategory(Analysis::Dilepton *dilepton);
   Double_t EvaluateLogPdf(Double_t mass, Double_t mass_error);
   Double_t EvaluateLogCrystalBall(Double_t mass, Double_t mass_error);
   Double_t EvaluateLogTruthTemplate(Double_t mass);
   void CopyFrom(ResolutionModel *theSource);

   ResolutionModel();
   ResolutionModel(Double_t constraint_mass, Double_t constraint_width);
   ResolutionModel(Double_t constraint_mass, Double_t constraint_width, Analysis::Quadrilepton *higgs, Int_t type, Int_t target = ResolutionModel::UseZ1);
   ~ResolutionModel();

protected:
   void InitializeCrystalBall();

};

namespace DileptonCategory {
   enum {
      bb,
      onecrk,
      b,
      other,
      all,
   };
};

#endif
