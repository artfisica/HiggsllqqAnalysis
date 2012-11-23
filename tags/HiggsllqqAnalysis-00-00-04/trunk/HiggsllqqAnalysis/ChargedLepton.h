#ifndef __ChargedLepton_h__
#define __ChargedLepton_h__

#include <iostream>
#include <TObject.h>

#include <TLorentzVector.h>

#include "HiggsAnalysis/MuonD3PDObject.h"
#include "HiggsAnalysis/ElectronD3PDObject.h"

namespace Analysis {
   class ChargedLepton : public TObject {
   public:
      enum {
         ELECTRON,
         MUON,
      };

      ChargedLepton();
      ChargedLepton(D3PDReader::MuonD3PDObjectElement *mu, Int_t family);
      ChargedLepton(D3PDReader::ElectronD3PDObjectElement *el, Int_t family);
      ~ChargedLepton();
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
      void set_lastcut(Int_t val) {
         m_lastcut = val;
      }
      void set_family(Int_t val) {
         m_family = val;
      }
      void set_charge(Float_t val) {
         m_charge = val;
      }
      void set_ptcone20(Float_t val) {
         m_ptcone20 = val;
      }
      void set_etcone20(Float_t val) {
         m_etcone20 = val;
      }
      void set_d0(Float_t val) {
         m_d0 = val;
      }
      void set_d0_sig(Float_t val) {
         m_d0_sig = val;
      }
      void set_z0(Float_t val) {
         m_z0 = val;
      }
      void set_z0_sig(Float_t val) {
         m_z0_sig = val;
      }
      Int_t flavor() {
         return m_flavor;
      }
      Int_t lastcut() {
         return m_lastcut;
      }
      Int_t family() {
         return m_family;
      }
      Float_t charge() {
         return m_charge;
      }
      Float_t ptcone20() {
         return m_ptcone20;
      }
      Float_t etcone20() {
         return m_etcone20;
      }
      Float_t d0() {
         return m_d0;
      }
      Float_t d0_sig() {
         return m_d0_sig;
      }
      Float_t z0() {
         return m_z0;
      }
      Float_t z0_sig() {
         return m_z0_sig;
      }

      void Set(D3PDReader::MuonD3PDObjectElement *mu, Int_t family);
      void Set(D3PDReader::ElectronD3PDObjectElement *el, Int_t family);
      D3PDReader::MuonD3PDObjectElement *GetMuon() {
         return m_muon;
      }
      D3PDReader::ElectronD3PDObjectElement *GetElectron() {
         return m_electron;
      }


   private:
      void Reset();

      Int_t   m_flavor;
      Int_t   m_lastcut;
      Int_t   m_family;

      Float_t m_charge;
      Float_t m_ptcone20;
      Float_t m_etcone20;
      Float_t m_d0;
      Float_t m_d0_sig;
      Float_t m_z0;
      Float_t m_z0_sig;
      TLorentzVector *m_momentum_cb;
      TLorentzVector *m_momentum_id;
      TLorentzVector *m_momentum_sa; // sa = MS (muon), CALO (electron)

      D3PDReader::MuonD3PDObjectElement *m_muon;
      D3PDReader::ElectronD3PDObjectElement *m_electron;

   public:

      ClassDef(ChargedLepton, 1);
   };
};
#endif
