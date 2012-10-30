// author: Arturo Sanchez <arturos@cern.ch>
// see Jet.h for details

#define __Jet_C__

#include "HiggsllqqAnalysis/CommonTools.h"
#include "HiggsllqqAnalysis/Jet.h"

ClassImp(Analysis::Jet);

Analysis::Jet::Jet()
{
   Reset();
}

Analysis::Jet::Jet(D3PDReader::JetD3PDObjectElement *jet)
{
   Reset();
   Set(jet);
}

Analysis::Jet::~Jet()
{
   delete m_momentum;
}

void Analysis::Jet::Reset()
{
  set_lastcut(-9999);
  set_rightpt(-9999.9);
  set_righteta(-9999.9);
  set_rightphi(-9999.9);
  set_rightE(-9999.9);
  set_rightEt(-9999.9);
  m_momentum = new TLorentzVector();
  m_jet = 0;
}

void Analysis::Jet::Set(D3PDReader::JetD3PDObjectElement *jet)
{
   set_rightpt(jet->pt());
   set_righteta(jet->eta());
   set_rightphi(jet->phi());
   set_rightE(jet->E());
   set_rightEt(jet->emscale_pt()); //Correguir y calcular Et

   Get4Momentum()->SetPtEtaPhiM(jet->pt(), jet->eta(), jet->phi(), jet->m());
   m_jet = jet;
}
