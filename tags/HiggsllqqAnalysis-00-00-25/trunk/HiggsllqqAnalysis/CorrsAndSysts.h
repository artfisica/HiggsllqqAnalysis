#ifndef CorrsAndSysts_h
#define CorrsAndSysts_h

/*
 * CorrsAndSysts class
 *
 * supply functions to apply corrections
 * and systematics needed for the H->bb
 * analysis
 *
 *  G. Facini, N. Morange & D. Buescher
 *  Wed Dec 12 13:07:00 CET 2012
 */

//
// Take all corrections and systematics that affect shapes out of Heather's script
//
// You should apply:
//
// NLO EW Higgs corrections
//
// W, Z, and top backgrounds shape corrections.
//
// pT and Mbb shape systematics. Cut-based values are here, for first estimates/checks.
// Should be replaced with systematics better tailored for MVA.
//
// Most pT systematics are binned in pT, each bin should be varied independently.
// The other are continuous parametrizations.
//


#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>

#include "TH1F.h"
#include "TF1.h"
#include "TString.h"
#include "TFile.h"
#include "TObject.h"

namespace CAS {

  enum EventType {
    WHlvbb=0,
    qqZHllbb,
    qqZHvvbb,
    ggZHllbb, //**added by Lei: 2014Feb01
    ggZHvvbb, //**added by Lei: 2014Feb01
    Wb,
    Wc,
    Wcc,
    Wl,
    Zb,
    Zc,
    Zcc,
    Zl,
    ttbar,
    stop_Wt,
    stop_s,
    stop_t,
    WW,
    WZ,
    ZZ,
    multijet,
    NONAME
  };
  enum DetailEventType {
    WWincl = 0,
    WlnuZhad,
    WhadZll,
    WhadZnunu,
    ZnunuZhad,
    ZllZhad,
    WWHerwig,
    WZHerwig,
    ZZHerwig,
    NODETAILNAME
  };
  enum SysVar { Do=0, Up=1 };
  enum SysBin {
    None=-2,                            // default
    Any=-1,                             // means no binning -> used for continuous systematics
    Bin0=0, Bin1, Bin2, Bin3, Bin4,     // pT bins for binned systematics
    NOTDEFINED
  };

  enum Systematic // only systematics affecting the shape are relevant to this tool
  {
    Nominal,
    // the following are binned in pT (5 bins, independent systematics)
    PTBINNED,
    // the following are continuous (thus correlated in pT)
    CONTINUOUS,
    SysTtbarPtWCont,
    SysTtbarMBBCont,
    SysTtbarPtB1Cont,
    SysTtbarMetCont,
    SysLepVeto,
    SysWtChanAcerMC,
    SysWtChanPythiaHerwig,
    SysTChanPtB2,
    SysHerwigPt,
    SysVVJetScalePt,
    SysVVJetPDFAlphaPt,
    SysVVJetScalePtST1,
    SysVVJetScalePtST2,
    SysTheoryHPt,
    SysTheoryVPtQCD, //** added by Lei: 2014jan31
    SysTopPt,
    SysTopMbb,
    SysWbbMbbGS, //inesochoa
    SysZMbb,
    SysWMbb,
    SysVVMbb,
    SysWDPhi,
    SysZDPhi,
    SysZPtV,
    SysWPtV,
    SysTruthTagDR,
    SysMJEleDR,
    SysMJElePtV,
    LAST
  };
} // end namespace CAS


class CorrsAndSysts {

  public:

    ~CorrsAndSysts();

  private:

    CorrsAndSysts(){};

    void Initialize();

    int m_debug;

    // string to enums or enums to string
    std::map<std::string, CAS::EventType> m_typeNames;
    std::map<CAS::Systematic, std::string> m_systNames;
    std::map<CAS::SysVar, std::string> m_varNames;
    std::map<CAS::SysBin, std::string> m_binNames;
    std::map<std::string, CAS::Systematic> m_systFromNames;
    std::map<std::string, CAS::SysVar> m_varFromNames;
    std::map<std::string, CAS::SysBin> m_binFromNames;

    Float_t pTbins[6];

    bool m_draw;
    bool m_zero;
    bool m_one;
    bool m_two;
    bool m_seven;
    bool m_eight;

    bool m_HZZ; //Flag for H->ZZ corrections

    // use TH* to store weights ; seems easier to maintain if ever we need e.g 2-d corrections
    // corrections
    TH1F* m_h_WHlvbbNLOEWKCorrection;
    TH1F* m_h_ZHllbbNLOEWKCorrection;
    TH1F* m_h_ZHvvbbNLOEWKCorrection;

    TH1F* m_h_ggZllH110_2JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZllH110_3JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZllH125_2JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZllH125_3JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZllH140_2JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZllH140_3JetCorrection; //** added by Lei: 2014jan31

    TH1F* m_h_ggZvvH110_2JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZvvH110_3JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZvvH125_2JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZvvH125_3JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZvvH140_2JetCorrection; //** added by Lei: 2014jan31
    TH1F* m_h_ggZvvH140_3JetCorrection; //** added by Lei: 2014jan31

    TH1F* m_h_pTbins;

    TH1F* m_h_topPtCorrection;

    // binned systematics
    TH1F* m_h_SysHerwigPt;

    TH1F* m_h_SysWW2JetInclScalePt;
    TH1F* m_h_SysWW2JetScalePt;
    TH1F* m_h_SysWW3JetScalePt;
    TH1F* m_h_SysWW3JetFor2JetScalePt;
    TH1F* m_h_SysWlnuZhad2JetInclScalePt;
    TH1F* m_h_SysWlnuZhad2JetScalePt;
    TH1F* m_h_SysWlnuZhad3JetScalePt;
    TH1F* m_h_SysWlnuZhad3JetFor2JetScalePt;
    TH1F* m_h_SysWhadZnunu2JetInclScalePt;
    TH1F* m_h_SysWhadZnunu2JetScalePt;
    TH1F* m_h_SysWhadZnunu3JetScalePt;
    TH1F* m_h_SysWhadZnunu3JetFor2JetScalePt;
    TH1F* m_h_SysWhadZll2JetInclScalePt;
    TH1F* m_h_SysWhadZll2JetScalePt;
    TH1F* m_h_SysWhadZll3JetScalePt;
    TH1F* m_h_SysWhadZll3JetFor2JetScalePt;
    TH1F* m_h_SysZllZhad2JetInclScalePt;
    TH1F* m_h_SysZllZhad2JetScalePt;
    TH1F* m_h_SysZllZhad3JetScalePt;
    TH1F* m_h_SysZllZhad3JetFor2JetScalePt;
    TH1F* m_h_SysZnunuZhad2JetInclScalePt;
    TH1F* m_h_SysZnunuZhad2JetScalePt;
    TH1F* m_h_SysZnunuZhad3JetScalePt;
    TH1F* m_h_SysZnunuZhad3JetFor2JetScalePt;

    TH1F* m_h_SysWW2JetPDFAlphaPt;
    TH1F* m_h_SysWW3JetPDFAlphaPt;
    TH1F* m_h_SysWlnuZhad2JetPDFAlphaPt;
    TH1F* m_h_SysWlnuZhad3JetPDFAlphaPt;
    TH1F* m_h_SysWhadZnunu2JetPDFAlphaPt;
    TH1F* m_h_SysWhadZnunu3JetPDFAlphaPt;
    TH1F* m_h_SysWhadZll2JetPDFAlphaPt;
    TH1F* m_h_SysWhadZll3JetPDFAlphaPt;
    TH1F* m_h_SysZllZhad2JetPDFAlphaPt;
    TH1F* m_h_SysZllZhad3JetPDFAlphaPt;
    TH1F* m_h_SysZnunuZhad2JetPDFAlphaPt;
    TH1F* m_h_SysZnunuZhad3JetPDFAlphaPt;

    TH1F* m_h_SysLepVeto2JetTop;
    TH1F* m_h_SysLepVeto2JetStop;
    TH1F* m_h_SysLepVeto2JetWbb;
    TH1F* m_h_SysLepVeto2JetWW;
    TH1F* m_h_SysLepVeto2JetWZ;
    TH1F* m_h_SysLepVeto3JetTop;
    TH1F* m_h_SysLepVeto3JetStop;
    TH1F* m_h_SysLepVeto3JetWbb;
    TH1F* m_h_SysLepVeto3JetWW;
    TH1F* m_h_SysLepVeto3JetWZ;

    // continuous systematics

    //new ttbar systematics
    TF1* m_f_SysTtbarPtVlpt;
    TF1* m_f_SysTtbarPtVhpt;
    TF1* m_f_SysTtbarMBB2jetlpt;
    TF1* m_f_SysTtbarMBB2jethpt;
    TF1* m_f_SysTtbarMBB3jetlpt;
    TF1* m_f_SysTtbarMBB3jethpt;
    TF1* m_f_SysTtbarPtB1lpt;
    TF1* m_f_SysTtbarPtB1hpt;
    TF1* m_f_SysTtbarMetlpt;
    TF1* m_f_SysTtbarMethpt;

    TH1F* m_h_SysTheoryWHlvbbPt;
    TH1F* m_h_SysTheoryZHllbbPt;
    TH1F* m_h_SysTheoryZHvvbbPt;

    //** added by Lei: 2014jan31
    TF1* m_f_SysTheoryqqVH2JetPtQCD;
    TF1* m_f_SysTheoryggZH2JetPtQCD;
    TF1* m_f_SysTheoryqqVH3JetPtQCD;
    TF1* m_f_SysTheoryggZH3JetPtQCD;

    //TF1* m_f_WDPhiCorr;
    TF1* m_f_WDPhiLowPt2JCorr;
    TF1* m_f_WDPhiHiPt2JCorr;
    TF1* m_f_WDPhiLowPt3JCorr;
    TF1* m_f_WDPhiHiPt3JCorr;
    //TF1* m_f_ZDPhiCorr;
    TF1* m_f_ZDPhiLowPt0TagCorr;
    TF1* m_f_ZDPhiHiPt0TagCorr;
    TF1* m_f_ZDPhiHiPt0TagHZZCorrErr;
    TF1* m_f_ZPtVCorr;
    TF1* m_f_ZbPtVCorr;
    TF1* m_f_SysWbbPtW; //inesochoa
    TF1* m_f_SysZbbMbb;
    TF1* m_f_SysWbbMbb;
    TF1* m_f_SysWbbMbbGS; //inesochoa
    TF1* m_f_SysTopMbb;
    TF1* m_f_SysWMbb;
    TF1* m_f_SysZMbb;
    TF1* m_f_SysWWMbb;
    TF1* m_f_SysWZMbb;
    TF1* m_f_SysZZMbb;

    TF1* m_f_SysTChanPt;
    TF1* m_f_SysWtPt2jLowPt;
    TF1* m_f_SysWtPt2jHighPt;
    TF1* m_f_SysWtPt3jLowPt;
    TF1* m_f_SysWtPt3jHighPt;

    TF1* m_f_dRjj; // truth-tagging deltaR correction/syst

    // multijet specific
    TH1F* m_h_SysMJ_dRBB_2lltag;
    TH1F* m_h_SysMJ_dRBB_2mmtag;
    TH1F* m_h_SysMJ_dRBB_2tttag;
    TH1F* m_h_SysMJ_pTV_2lltag;
    TH1F* m_h_SysMJ_pTV_2mmtag;
    TH1F* m_h_SysMJ_pTV_2tttag;

    std::map<TString, TObject*> m_allHists;

  public:

    CorrsAndSysts(TString name, bool draw=true, bool HZZ=false); // e.g OneLepton_8TeV
    CorrsAndSysts(int channel, int year, bool draw=true, bool HZZ=false); // channel: 0->0lepton, 1->1lepton, 2->2leptons

    inline void SetDebug(int i) { m_debug = i; }

    CAS::EventType GetEventType(TString name);
    CAS::SysBin GetSysBin(float vpt);

    void WriteHistsToFile(TString fname);

    // all values (VpT, Mbb) in MeV !

    // Higgs pT reweighting (NLO EW corrections)
    float Get_HiggsNLOEWKCorrection(CAS::EventType type, float VpT);
    inline float Get_HiggsNLOEWKCorrection(TString evtType, float VpT)
    { return Get_HiggsNLOEWKCorrection(m_typeNames[evtType.Data()], VpT); }

    //** added by Lei: 2014jan31
    // Z pT reweighting (gg->ZH corrections)
    float Get_ggZHCorrection(CAS::EventType type, float VpT, int njet);
    inline float Get_ggZHCorrection(TString evtType, float VpT, int njet)
    { return Get_ggZHCorrection(m_typeNames[evtType.Data()], VpT, njet); }

    // Top pT reweighting
    float Get_ToppTCorrection(CAS::EventType type, float avgTopPt);
    inline float Get_ToppTCorrection(TString evtType, float avgTopPt)
    { return Get_ToppTCorrection(m_typeNames[evtType.Data()], avgTopPt); }

    // Herwig diboson VpT correction
    float Get_DibosonCorrection(CAS::EventType type, float VpT, int mc_channel_number);
    float Get_DibosonCorrection(TString evtType, float VpT, int mc_channel_number);

    // same without mc_channel_number
    float Get_DibosonCorrection(CAS::EventType type, float VpT, CAS::DetailEventType detailtype);
    inline float Get_DibosonCorrection(TString evtType, float VpT, CAS::DetailEventType detailtype)
    { return Get_DibosonCorrection(m_typeNames[evtType.Data()], VpT, detailtype); }

    // Bkg DeltaPhi correction
    float Get_BkgDeltaPhiCorrection(CAS::EventType type, float DeltaPhi, int njet, float ptv);
    inline float Get_BkgDeltaPhiCorrection(TString evtType, float DeltaPhi, int njet, float ptv)
    { return Get_BkgDeltaPhiCorrection(m_typeNames[evtType.Data()], DeltaPhi, njet, ptv); }

    float Get_BkgPtVCorrection(CAS::EventType type, float ptv, int njet) ;
    inline float Get_BkgPtVCorrection(TString evtType, float ptv, int njet)
    { return Get_BkgPtVCorrection(m_typeNames[evtType.Data()], ptv, njet); }

    // DeltaR truth-tagging correction
    float Get_DeltaRTruthTagCorrection(CAS::EventType type, float deltaR, int ntag) ;
    inline float Get_DeltaRTruthTagCorrection(TString evtType, float deltaR, int ntag)
    { return Get_DeltaRTruthTagCorrection(m_typeNames[evtType.Data()], deltaR, ntag); }


    // Systematics on distributions
    float Get_SystematicWeight(CAS::EventType type, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR,
                               float pTB1, float pTB2, float met, int njet, int ntag, int mc_channel_number,
                               CAS::Systematic sys=CAS::Nominal, CAS::SysVar var=CAS::Up, CAS::SysBin bin=CAS::None);
    float Get_SystematicWeight(CAS::EventType type, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR,
                               float pTB1, float pTB2, float met, int njet, int ntag, int mc_channel_number, TString sysName);

    float Get_SystematicWeight(TString evtType, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR, float pTB1,
                               float pTB2, float met, int njet, int ntag, int mc_channel_number, CAS::Systematic sys,
                               CAS::SysVar var, CAS::SysBin bin);
    float Get_SystematicWeight(TString evtType, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR, float pTB1,
                               float pTB2, float met, int njet, int ntag, int mc_channel_number, TString sysName);


    // same without mc_channel_number
    float Get_SystematicWeight(CAS::EventType type, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR, float pTB1,
                               float pTB2, float met, int njet, int ntag, CAS::DetailEventType detailtype,
                               CAS::Systematic sys=CAS::Nominal, CAS::SysVar var=CAS::Up, CAS::SysBin bin=CAS::None);
    float Get_SystematicWeight(CAS::EventType type, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR, float pTB1,
                               float pTB2, float met, int njet, int ntag, CAS::DetailEventType detailtype, TString sysName);

    inline float Get_SystematicWeight(TString evtType, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR, float pTB1,
                                      float pTB2, float met, int njet, int ntag, CAS::DetailEventType detailtype,
                                      CAS::Systematic sys, CAS::SysVar var, CAS::SysBin bin)
    { return Get_SystematicWeight(m_typeNames[evtType.Data()], VpT, Mbb, truthPt, DeltaPhi, deltaR, pTB1, pTB2, met, njet, ntag, detailtype, sys, var, bin); }

    inline float Get_SystematicWeight(TString evtType, float VpT, float Mbb, float truthPt, float DeltaPhi, float deltaR, float pTB1,
                                      float pTB2, float met, int njet, int ntag, CAS::DetailEventType detailtype, TString sysName)
    { return Get_SystematicWeight(m_typeNames[evtType.Data()], VpT, Mbb, truthPt, DeltaPhi, deltaR, pTB1, pTB2, met, njet, ntag, detailtype, sysName); }

    // multijet specific RW and Syst
    float Get_dRBB_pTV_MJCorrection(CAS::EventType type, bool isEl, float dRBB, float pTV, int nTag, int nJet) const;
    float Get_MultijetSystematic(CAS::EventType type, bool isEl, float dRBB, float pTV, int nTag, int nJet, CAS::Systematic syst, CAS::SysVar direction) const;

    // forge the normalized syste name from the enums
    inline std::string GetSystName(CAS::Systematic sys, CAS::SysBin bin, CAS::SysVar var)
    { return m_systNames[sys]+m_binNames[bin]+m_varNames[var]; }

    // inverse function
    void GetSystFromName(TString name, CAS::Systematic& sys, CAS::SysBin& bin, CAS::SysVar& var);

};// close CorrsAndSysts class

namespace Utils {

  // utility
  TH1F* BuildTH1F(std::vector<Double_t> contents, TString hname, float min, float max, std::map<TString, TObject*>& hists);
  void  FillTH1F(std::vector<Float_t> contents, TH1F* h, std::map<TString, TObject*>& hists);
  void  FillTH1F(int len, Float_t* contents, TH1F* h, std::map<TString, TObject*>& hists);
  inline float GetScale(float value, TH1F* h);
  void SaveHist(TObject* h, std::map<TString, TObject*>& hists);
  void ArraySubstractOne(float* array, unsigned int length);
  CAS::DetailEventType GetDibosonType(int mc_channel_number, CAS::EventType type);

  // map<K,V> => map<V,K>
  template <typename T, typename U>
    std::map<U, T> reverseMap(const std::map<T, U>& m_in);

  // Implementation
  //
  template <typename T, typename U>
    std::map<U, T> reverseMap(const std::map<T, U>& m_in) {
      typedef typename std::map<T, U>::const_iterator map_it;
      map_it it=m_in.begin();
      std::map<U, T> m_out;
      while(it!=m_in.end()) {
        m_out[it->second]=it->first;
        it++;
      }
      return m_out;
    }

} // close Utils class

#endif //CorrsAndSysts_HPP_

