#ifndef __ConstraintFitResult_h__
#define __ConstraintFitResult_h__

#include "HiggsllqqAnalysis/Quadrilepton.h"

namespace Analysis {
   class ConstraintFitResult {
   public:
      TLorentzVector H_4m_updated;
      TLorentzVector Z1_4m_updated;
      TLorentzVector Z1_lepplus_4m_updated;
      TLorentzVector Z1_lepminus_4m_updated;
      TLorentzVector Z2_4m_updated;
      TLorentzVector Z2_lepplus_4m_updated;
      TLorentzVector Z2_lepminus_4m_updated;

      Double_t Z1_chiSq;
      Double_t Z2_chiSq;

      ConstraintFitResult(Quadrilepton *higgs) {
         H_4m_updated = *(higgs->Get4Momentum());
         Z1_4m_updated = *(higgs->GetZ1()->Get4Momentum());
         Z1_lepplus_4m_updated = *(higgs->GetZ1()->GetLepPlus()->Get4Momentum());
         Z1_lepminus_4m_updated = *(higgs->GetZ1()->GetLepMinus()->Get4Momentum());
         Z2_4m_updated = *(higgs->GetZ2()->Get4Momentum());
         Z2_lepplus_4m_updated = *(higgs->GetZ2()->GetLepPlus()->Get4Momentum());
         Z2_lepminus_4m_updated = *(higgs->GetZ2()->GetLepMinus()->Get4Momentum());
         Z1_chiSq = -9999.9;
         Z2_chiSq = -9999.9;
      }
      ~ConstraintFitResult() {}
   };
};

#endif
