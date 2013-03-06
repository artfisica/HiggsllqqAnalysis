#ifndef __Jet_h__
#define __Jet_h__

#include <iostream>
#include <TObject.h>
#include <TMath.h>

#include <TLorentzVector.h>

#include "HiggsAnalysis/JetD3PDObject.h"

namespace Analysis {
  class Jet : public TObject {

  public:    
    Jet();
    Jet(D3PDReader::JetD3PDObjectElement *jet);
    ~Jet();

    TLorentzVector *Get4Momentum() {
      return m_momentum;
    }
    void set_lastcut(Int_t val) {
      m_lastcut = val;
    }
    void set_rightpt(Float_t val) {
      m_rightpt = val;
    }
    void set_righteta(Float_t val) {
      m_righteta = val;
    }
    void set_rightphi(Float_t val) {
      m_rightphi = val;
    }
    void set_rightE(Float_t val) {
      m_rightE = val;
    }
    void set_rightEt(Float_t val) {
      m_rightEt = val;
    }
    Int_t lastcut() {
      return m_lastcut;
    }
    Float_t rightpt() {
      return m_rightpt;
    }
    Float_t righteta() {
      return m_righteta;
    }
    Float_t rightphi() {
      return m_rightphi;
    }
    Float_t rightE() {
      return m_rightE;
    }
    Float_t rightEt() {
      return m_rightEt;
    }
    
    
    void Set(D3PDReader::JetD3PDObjectElement *jet);
    D3PDReader::JetD3PDObjectElement *GetJet() {
      return m_jet;
    }
    
    
  private:
    void Reset();

    Int_t   m_family;    
    Int_t   m_lastcut;
    Float_t m_rightpt;
    Float_t m_righteta;
    Float_t m_rightphi;
    Float_t m_rightE;
    Float_t m_rightEt;
    
    TLorentzVector *m_momentum;
    
    D3PDReader::JetD3PDObjectElement *m_jet;
    
  public:
    
    ClassDef(Jet, 1);
  };
};
#endif
