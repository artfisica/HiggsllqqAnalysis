#ifndef __Quadrilepton_h__
#define __Quadrilepton_h__

#include "HiggsllqqAnalysis/Dilepton.h"

namespace Analysis {
   class Quadrilepton : public TObject {
   public:
      enum {
         MU2,
         MUE,
         EMU,
         E2,
      };

      Quadrilepton();
      Quadrilepton(Analysis::Dilepton *Z1, Analysis::Dilepton *Z2);
      ~Quadrilepton();
      TLorentzVector *Get4Momentum() {
         return m_momentum_cb;
      }
      TLorentzVector *Get4Momentum_ID() {
         return m_momentum_id;
      }
      TLorentzVector *Get4Momentum_SA() {
         return m_momentum_sa;
      }
      Bool_t IsNeutral() {
         return m_is_neutral;
      }
      void set_channel(Int_t val) {
         m_channel = val;
      }
      Int_t channel() {
         return m_channel;
      }
      void set_lastcut(Int_t val) {
         m_lastcut = val;
      }
      Int_t lastcut() {
         return m_lastcut;
      }

      void Set(Analysis::Dilepton *Z1, Analysis::Dilepton *Z2);
      Analysis::Dilepton *GetZ1() {
         return m_Z1;
      }
      Analysis::Dilepton *GetZ2() {
         return m_Z2;
      }

      std::vector<Analysis::ChargedLepton *> GetLeptonVector() {
         return m_lepVector;
      }

   private:
      void Reset();

      Bool_t  m_is_neutral;
      Int_t   m_channel;
      Int_t   m_lastcut;
      TLorentzVector *m_momentum_cb;
      TLorentzVector *m_momentum_id;
      TLorentzVector *m_momentum_sa; // sa = MS (muon), CALO (electron)

      Analysis::Dilepton *m_Z1;
      Analysis::Dilepton *m_Z2;

      std::vector<Analysis::ChargedLepton *> m_lepVector;

   public:

      ClassDef(Quadrilepton, 1);
   };
};
#endif
