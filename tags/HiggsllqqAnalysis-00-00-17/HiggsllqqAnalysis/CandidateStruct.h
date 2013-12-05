#ifndef __CandidateStruct_h__
#define __CandidateStruct_h__

namespace Analysis {
  typedef struct {
    Int_t analysis;
    UInt_t run;
    UInt_t event;
    Int_t last;
  } CutflowStruct;
  
  typedef struct {
    Int_t last;
    Int_t closesttozz;
    Int_t selected;
    Int_t type;
    Int_t hfor;
    Int_t top_weight;
    Int_t powhegbug_weight;
    UInt_t run;
    UInt_t event;
    UInt_t lbn;
    Int_t n_vx;
    Int_t actualIntPerXing;
    Int_t averageIntPerXing;
    Float_t xsec_weight;
    Long64_t processed_entries;
    Long64_t generated_entries;
    Float_t pu_weight;
    Float_t vxz_weight;
    Float_t OQ_weight;
    Float_t ggF_weight;
    Float_t trigSF_weight;
    Float_t Z1_lepplus_weight;
    Float_t Z1_lepminus_weight;
    Float_t Z2_lepplus_weight;
    Float_t Z2_lepminus_weight;
    Float_t H_m_truth;
    Float_t ZZ_m_truth;
    Float_t Z1_m_truth;
    Float_t Z2_m_truth;
    Float_t H_m_constrained;
    Float_t Z1_m_constrained;
    Float_t Z2_m_constrained;
    Float_t Z1_chiSq;
    Float_t Z2_chiSq;
    Float_t H_m;
    Float_t H_pt;
    Float_t H_eta;
    Float_t H_phi;
    Float_t Z1_m;
    Float_t Z1_pt;
    Float_t Z1_eta;
    Float_t Z1_phi;
    Float_t Z2_m;
    Float_t Z2_pt;
    Float_t Z2_eta;
    Float_t Z2_phi;
    Float_t Z1_lepplus_m;
    Float_t Z1_lepplus_pt;
    Float_t Z1_lepplus_pt_truth;
    Float_t Z1_lepplus_pt_constrained;
    Float_t Z1_lepplus_eta;
    Float_t Z1_lepplus_phi;
    Float_t Z1_lepplus_etcone20;
    Float_t Z1_lepplus_etcone20_corr;
    Float_t Z1_lepplus_etcone20_final;
    Float_t Z1_lepplus_ptcone20;
    Float_t Z1_lepplus_ptcone20_final;
    Float_t Z1_lepplus_d0;
    Float_t Z1_lepplus_z0;
    Float_t Z1_lepplus_d0_sig;
    Float_t Z1_lepplus_GSF_dp;
    Float_t Z1_lepplus_GSF_p;
    Float_t Z1_lepplus_charge;
    Float_t Z1_lepplus_ptCB_nosmearnoscale;
    Float_t Z1_lepplus_ptME_nosmearnoscale;
    Float_t Z1_lepplus_ptID_nosmearnoscale;
    Float_t Z1_lepplus_cl_E_calibsmearnoscale;
    Float_t Z1_lepplus_cl_eta;
    Float_t Z1_lepplus_cl_phi;
    Int_t Z1_lepplus_isSA;
    Int_t Z1_lepplus_author;
    Int_t Z1_lepplus_type;
    Int_t Z1_lepplus_typebkg;
    Int_t Z1_lepplus_origin;
    Int_t Z1_lepplus_originbkg;
    Float_t Z1_lepminus_m;
    Float_t Z1_lepminus_pt;
    Float_t Z1_lepminus_pt_truth;
    Float_t Z1_lepminus_pt_constrained;
    Float_t Z1_lepminus_eta;
    Float_t Z1_lepminus_phi;
    Float_t Z1_lepminus_etcone20;
    Float_t Z1_lepminus_etcone20_corr;
    Float_t Z1_lepminus_etcone20_final;
    Float_t Z1_lepminus_ptcone20;
    Float_t Z1_lepminus_ptcone20_final;
    Float_t Z1_lepminus_d0;
    Float_t Z1_lepminus_z0;
    Float_t Z1_lepminus_d0_sig;
    Float_t Z1_lepminus_GSF_dp;
    Float_t Z1_lepminus_GSF_p;
    Float_t Z1_lepminus_charge;
    Float_t Z1_lepminus_ptCB_nosmearnoscale;
    Float_t Z1_lepminus_ptME_nosmearnoscale;
    Float_t Z1_lepminus_ptID_nosmearnoscale;
    Float_t Z1_lepminus_cl_E_calibsmearnoscale;
    Float_t Z1_lepminus_cl_eta;
    Float_t Z1_lepminus_cl_phi;
    Int_t Z1_lepminus_isSA;
    Int_t Z1_lepminus_author;
    Int_t Z1_lepminus_type;
    Int_t Z1_lepminus_typebkg;
    Int_t Z1_lepminus_origin;
    Int_t Z1_lepminus_originbkg;
    Float_t Z2_lepplus_m;
    Float_t Z2_lepplus_pt;
    Float_t Z2_lepplus_pt_truth;
    Float_t Z2_lepplus_pt_constrained;
    Float_t Z2_lepplus_eta;
    Float_t Z2_lepplus_phi;
    Float_t Z2_lepplus_etcone20;
    Float_t Z2_lepplus_etcone20_corr;
    Float_t Z2_lepplus_etcone20_final;
    Float_t Z2_lepplus_ptcone20;
    Float_t Z2_lepplus_ptcone20_final;
    Float_t Z2_lepplus_d0;
    Float_t Z2_lepplus_z0;
    Float_t Z2_lepplus_d0_sig;
    Float_t Z2_lepplus_GSF_dp;
    Float_t Z2_lepplus_GSF_p;
    Float_t Z2_lepplus_charge;
    Float_t Z2_lepplus_ptCB_nosmearnoscale;
    Float_t Z2_lepplus_ptME_nosmearnoscale;
    Float_t Z2_lepplus_ptID_nosmearnoscale;
    Float_t Z2_lepplus_cl_E_calibsmearnoscale;
    Float_t Z2_lepplus_cl_eta;
    Float_t Z2_lepplus_cl_phi;
    Int_t Z2_lepplus_isSA;
    Int_t Z2_lepplus_author;
    Int_t Z2_lepplus_type;
    Int_t Z2_lepplus_typebkg;
    Int_t Z2_lepplus_origin;
    Int_t Z2_lepplus_originbkg;
    Float_t Z2_lepminus_m;
    Float_t Z2_lepminus_pt;
    Float_t Z2_lepminus_pt_truth;
    Float_t Z2_lepminus_pt_constrained;
    Float_t Z2_lepminus_eta;
    Float_t Z2_lepminus_phi;
    Float_t Z2_lepminus_etcone20;
    Float_t Z2_lepminus_etcone20_corr;
    Float_t Z2_lepminus_etcone20_final;
    Float_t Z2_lepminus_ptcone20;
    Float_t Z2_lepminus_ptcone20_final;
    Float_t Z2_lepminus_d0;
    Float_t Z2_lepminus_z0;
    Float_t Z2_lepminus_d0_sig;
    Float_t Z2_lepminus_GSF_dp;
    Float_t Z2_lepminus_GSF_p;
    Float_t Z2_lepminus_charge;
    Float_t Z2_lepminus_ptCB_nosmearnoscale;
    Float_t Z2_lepminus_ptME_nosmearnoscale;
    Float_t Z2_lepminus_ptID_nosmearnoscale;
    Float_t Z2_lepminus_cl_E_calibsmearnoscale;
    Float_t Z2_lepminus_cl_eta;
    Float_t Z2_lepminus_cl_phi;
    Int_t Z2_lepminus_isSA;
    Int_t Z2_lepminus_author;
    Int_t Z2_lepminus_type;
    Int_t Z2_lepminus_typebkg;
    Int_t Z2_lepminus_origin;
    Int_t Z2_lepminus_originbkg;
    std::vector<std::string> filename;
  } CandidateStruct;
  
  //TestSelection Structures
  typedef struct {
    unsigned int runnumber;
    unsigned int eventnumber;
    int istagged;
    int channel;
    int isqcdevent;
    int n_jets;
    int n_b_jets;
    float weight;
    float lep1_m;
    float lep1_pt;
    float lep1_eta;
    float lep1_phi;
    float lep1_charge;
    float lep1_caloiso;
    float lep1_trackiso;
    float lep1_d0;
    float lep1_sigd0;
    int   lep1_quality;
    float lep2_m;
    float lep2_pt;
    float lep2_eta;
    float lep2_phi;
    float lep2_charge;
    float lep2_caloiso;
    float lep2_trackiso;
    float lep2_d0;
    float lep2_sigd0;
    int   lep2_quality;
    float lepZ_m;
    float lepZ_pt;
    float lepZ_eta;
    float lepZ_phi;
    float realJ1_m;
    float realJ1_pt;
    float realJ1_eta;
    float realJ1_eta_det;
    float realJ1_phi;
    float realJ1_flavortruth;
    float realJ1_pdg;
    float realJ1_jvf;
    float realJ1_ntrk;
    float realJ1_width;
    float realJ1_MV1;
    float realJ1_Fisher;
    float realJ2_m;
    float realJ2_pt;
    float realJ2_eta;
    float realJ2_eta_det;
    float realJ2_phi;
    float realJ2_flavortruth;
    float realJ2_pdg;
    float realJ2_jvf;
    float realJ2_ntrk;
    float realJ2_width;
    float realJ2_MV1;
    float realJ2_Fisher;
    float ll_2_jets;
    float realZ_m;
    float realZ_pt;
    float realZ_eta;
    float realZ_phi;
    float realH_m;
    float realH_pt;
    float realH_eta;
    float realH_phi;
    float corrJ1_m;
    float corrJ1_pt;
    float corrJ1_eta;
    float corrJ1_eta_det;
    float corrJ1_phi;
    float corrJ1_flavortruth;
    float corrJ1_Fisher;
    float corrJ2_m;
    float corrJ2_pt;
    float corrJ2_eta;
    float corrJ2_eta_det;
    float corrJ2_phi;
    float corrJ2_flavortruth;
    float corrJ2_Fisher;
    float ll_2_jets_corr;
    float corrZ_m;
    float corrZ_pt;
    float corrZ_eta;
    float corrZ_phi;
    float corrH_m;
    float corrH_pt;
    float corrH_eta;
    float corrH_phi;
    float chisquare;
    float met;
    float sumet;
    float btagSF;
    int NPV;
    float truthH_pt;
    float ggFweight;
    float xWin_44b;
    float yWin_44b;
    float zWin_44b;
    float gWin_44b;
    float xWin_44g;
    float yWin_44g;
    float zWin_44g;
    float gWin_44g;
    float xWin_55b;
    float yWin_55b;
    float zWin_55b;
    float gWin_55b;
    float xWin_55g;
    float yWin_55g;
    float zWin_55g;
    float gWin_55g;
    float mu;
    float trig_SF;
    float trig_SF2;
    float trig_SFC;
    int trig_flag;
    int HFOR;
    int Entries;
    //Flavour Composition Variables
    float SecondJet_MV1_b; // b jets 
    float SecondJet_MV1_c; // c jets
    float SecondJet_MV1_l; // light jets
    float AllJet_MV1_b; // b jets 
    float AllJet_MV1_c; // c jets
    float AllJet_MV1_l; // light jets
  } analysis_output_struct;
  
  typedef struct SOMVar{
    float Ntrk, width;
  } SOMVar;
};

#endif