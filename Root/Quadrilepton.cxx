// author: Valerio Ippolito <valerio.ippolito@cern.ch>
// see D3PDReadingInterface/Quadrilepton.h for details

#define __Quadrilepton_C__

#include "HiggsllqqAnalysis/Quadrilepton.h"

ClassImp(Analysis::Quadrilepton);

Analysis::Quadrilepton::Quadrilepton()
{
   Reset();
}

Analysis::Quadrilepton::Quadrilepton(Dilepton *Z1, Dilepton *Z2)
{
   Reset();
   Set(Z1, Z2);
}

Analysis::Quadrilepton::~Quadrilepton()
{
   delete m_momentum_cb;
   delete m_momentum_id;
   delete m_momentum_sa;
}

void Analysis::Quadrilepton::Reset()
{
   set_channel(-9999);
   set_lastcut(-9999);
   m_momentum_cb = new TLorentzVector();
   m_momentum_id = new TLorentzVector();
   m_momentum_sa = new TLorentzVector();
   m_Z1 = 0;
   m_Z2 = 0;
   m_lepVector.clear();
   m_is_neutral = kFALSE;
}

void Analysis::Quadrilepton::Set(Dilepton *Z1, Dilepton *Z2)
{
   Int_t mychannel(-1);

   if (Z1->flavor() == ChargedLepton::MUON) {
      if (Z2->flavor() == ChargedLepton::MUON) {
         mychannel = Quadrilepton::MU2;
      } else if (Z2->flavor() == ChargedLepton::ELECTRON) {
         mychannel = Quadrilepton::MUE;
      }
   } else if (Z1->flavor() == ChargedLepton::ELECTRON) {
      if (Z2->flavor() == ChargedLepton::MUON) {
         mychannel = Quadrilepton::EMU;
      } else if (Z2->flavor() == ChargedLepton::ELECTRON) {
         mychannel = Quadrilepton::E2;
      }
   }

   set_channel(mychannel);
   m_Z1 = Z1;
   m_Z2 = Z2;

   m_is_neutral = (m_Z1->IsNeutral() && m_Z2->IsNeutral());

   *m_momentum_cb = *m_Z1->Get4Momentum()    + *m_Z2->Get4Momentum();
   *m_momentum_id = *m_Z1->Get4Momentum_ID() + *m_Z2->Get4Momentum_ID();
   *m_momentum_sa = *m_Z1->Get4Momentum_SA() + *m_Z2->Get4Momentum_SA();

   m_lepVector.clear();
   m_lepVector.push_back(m_Z1->GetLepPlus());
   m_lepVector.push_back(m_Z1->GetLepMinus());
   m_lepVector.push_back(m_Z2->GetLepPlus());
   m_lepVector.push_back(m_Z2->GetLepMinus());
}
