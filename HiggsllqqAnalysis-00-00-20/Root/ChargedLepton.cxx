// author: Valerio Ippolito <valerio.ippolito@cern.ch>
// see D3PDReadingInterface/ChargedLepton.h for details

#define __ChargedLepton_C__

#include "HiggsllqqAnalysis/CommonTools.h"
#include "HiggsllqqAnalysis/ChargedLepton.h"

ClassImp(Analysis::ChargedLepton);

Analysis::ChargedLepton::ChargedLepton()
{
   Reset();
}

Analysis::ChargedLepton::ChargedLepton(D3PDReader::MuonD3PDObjectElement *mu, Int_t family)
{
   Reset();
   Set(mu, family);
}

Analysis::ChargedLepton::ChargedLepton(D3PDReader::ElectronD3PDObjectElement *el, Int_t family)
{
   Reset();
   Set(el, family);
}

Analysis::ChargedLepton::~ChargedLepton()
{
   delete m_momentum_cb;
   delete m_momentum_id;
   delete m_momentum_sa;
}

void Analysis::ChargedLepton::Reset()
{
   set_flavor(-9999);
   set_lastcut(-9999);
   set_family(-9999);
   set_charge(-9999.9);
   set_ptcone20(-9999.9);
   set_etcone20(-9999.9);
   set_d0(-9999.9);
   set_d0_sig(-9999.9);
   set_z0(-9999.9);
   set_z0_sig(-9999.9);
   m_momentum_cb = new TLorentzVector();
   m_momentum_id = new TLorentzVector();
   m_momentum_sa = new TLorentzVector();
   m_muon = 0;
   m_electron = 0;
}

void Analysis::ChargedLepton::Set(D3PDReader::MuonD3PDObjectElement *mu, Int_t family)
{
   set_flavor(ChargedLepton::MUON);
   set_family(family);
   set_charge(mu->charge());
   set_ptcone20(mu->ptcone20());
   set_etcone20(mu->etcone20());
   set_d0(mu->trackIPEstimate_d0_unbiasedpvunbiased());
   set_d0_sig(mu->trackIPEstimate_sigd0_unbiasedpvunbiased());
   set_z0(mu->trackIPEstimate_z0_unbiasedpvunbiased());
   set_z0_sig(mu->trackIPEstimate_sigz0_unbiasedpvunbiased());

   Get4Momentum()->SetPtEtaPhiM(mu->pt(), mu->eta(), mu->phi(), mu_pdg_mass);
   Get4Momentum_ID()->SetPtEtaPhiM(TMath::Abs(1. / mu->id_qoverp()*TMath::Sin(mu->id_theta())), -TMath::Log(TMath::Tan(0.5 * mu->id_theta())), mu->id_phi(), mu_pdg_mass);
   Get4Momentum_SA()->SetPtEtaPhiM(TMath::Abs(1. / mu->me_qoverp()*TMath::Sin(mu->me_theta())), -TMath::Log(TMath::Tan(0.5 * mu->me_theta())), mu->me_phi(), mu_pdg_mass);

   m_muon = mu;
   m_electron = 0;
}

void Analysis::ChargedLepton::Set(D3PDReader::ElectronD3PDObjectElement *el, Int_t family)
{
   set_flavor(ChargedLepton::ELECTRON);
   set_family(family);
   set_charge(el->charge());
   set_ptcone20(el->ptcone20());
   set_etcone20(el->Etcone20());
   set_d0(el->trackd0pvunbiased());
   set_d0_sig(el->tracksigd0pvunbiased());
   set_z0(el->trackz0pvunbiased());
   set_z0_sig(el->tracksigz0pvunbiased());

   Get4Momentum()->SetPtEtaPhiM(el->cl_E() / TMath::CosH(el->tracketa()), el->tracketa(), el->trackphi(), el_pdg_mass);
   Get4Momentum_ID()->SetPtEtaPhiM(el->trackpt(), el->tracketa(), el->trackphi(), el_pdg_mass);
   Get4Momentum_SA()->SetPtEtaPhiM(el->cl_pt(), el->cl_eta(), el->cl_phi(), el_pdg_mass);

   m_muon = 0;
   m_electron = el;
}
