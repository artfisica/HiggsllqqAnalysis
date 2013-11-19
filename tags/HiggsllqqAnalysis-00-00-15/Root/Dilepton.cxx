// author: Valerio Ippolito <valerio.ippolito@cern.ch>
// see D3PDReadingInterface/Dilepton.h for details

#define __Dilepton_C__

#include "HiggsllqqAnalysis/Dilepton.h"

ClassImp(Analysis::Dilepton);

Analysis::Dilepton::Dilepton()
{
   Reset();
}

Analysis::Dilepton::Dilepton(ChargedLepton *lepplus, ChargedLepton *lepminus)
{
   Reset();
   Set(lepplus, lepminus);
}

Analysis::Dilepton::~Dilepton()
{
   delete m_momentum_cb;
   delete m_momentum_id;
   delete m_momentum_sa;
}

void Analysis::Dilepton::Reset()
{
   set_flavor(-9999);
   m_momentum_cb = new TLorentzVector();
   m_momentum_id = new TLorentzVector();
   m_momentum_sa = new TLorentzVector();
   m_lepplus = 0;
   m_lepminus = 0;
   m_is_neutral = kFALSE;
}

void Analysis::Dilepton::Set(ChargedLepton *lepplus, ChargedLepton *lepminus)
{
   // flavor protection
   if (lepplus->flavor() != lepminus->flavor()) {
      return;
   }

   m_is_neutral = (lepplus->charge() * lepminus->charge() < 0);

   set_flavor(lepplus->flavor());
   m_lepplus = lepplus;
   m_lepminus = lepminus;

   *m_momentum_cb = *m_lepplus->Get4Momentum()    + *m_lepminus->Get4Momentum();
   *m_momentum_id = *m_lepplus->Get4Momentum_ID() + *m_lepminus->Get4Momentum_ID();
   *m_momentum_sa = *m_lepplus->Get4Momentum_SA() + *m_lepminus->Get4Momentum_SA();
}

Bool_t Analysis::Dilepton::OverlapsWith(Dilepton *dilep)
{
   return (GetLepPlus()  == dilep->GetLepPlus()  ||
           GetLepMinus() == dilep->GetLepMinus() ||
           GetLepPlus()  == dilep->GetLepMinus() ||
           GetLepMinus() == dilep->GetLepPlus());
}
