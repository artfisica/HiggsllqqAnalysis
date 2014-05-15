#ifndef __Dilepton_h__
#define __Dilepton_h__

#include "HiggsllqqAnalysis/ChargedLepton.h"

namespace Analysis {
   class Dilepton : public TObject {
   public:
      Dilepton();
      Dilepton(Analysis::ChargedLepton *lepplus, Analysis::ChargedLepton *lepminus);
      ~Dilepton();
      TLorentzVector *Get4Momentum() {
         return m_momentum_cb;
      }
      TLorentzVector *Get4Momentum_ID() {
         return m_momentum_id;
      }
      TLorentzVector *Get4Momentum_SA() {
         return m_momentum_sa;
      }
      void set_flavor(Int_t val) {
         m_flavor = val;
      }
      Int_t flavor() {
         return m_flavor;
      }
      Bool_t IsNeutral() {
         return m_is_neutral;
      }

      void Set(Analysis::ChargedLepton *lepplus, Analysis::ChargedLepton *lepminus);
      Analysis::ChargedLepton *GetLepPlus() {
         return m_lepplus;
      }
      Analysis::ChargedLepton *GetLepMinus() {
         return m_lepminus;
      }

      Bool_t OverlapsWith(Dilepton *dilep);

   private:
      void Reset();

      Int_t   m_flavor;
      Bool_t  m_is_neutral;
      TLorentzVector *m_momentum_cb;
      TLorentzVector *m_momentum_id;
      TLorentzVector *m_momentum_sa; // sa = MS (muon), CALO (electron)

      Analysis::ChargedLepton *m_lepplus;
      Analysis::ChargedLepton *m_lepminus;

   public:

      ClassDef(Dilepton, 1);
   };
};
#endif
