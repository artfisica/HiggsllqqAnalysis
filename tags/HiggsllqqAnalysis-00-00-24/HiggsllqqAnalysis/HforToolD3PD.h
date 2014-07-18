
/**
   @class    HforToolD3PD.h

   @brief    tool for removal of heavy flavor overlap in Alpgen samples

   @author   Dominic Hirschbuehl, Michel Saners (original code in athena)
   @author   Takashi Yamanaka (converted the code for D3PD)

   taken from  Generators/GenAnalysisTools/TruthUtils/tags/tags/TruthUtils-00-01-00

   30 March 2012

   For the detail of usage, see https://twiki.cern.ch/twiki/bin/view/Main/HforToolD3PD
*/

#ifndef HFORTOOLD3PD_H
#define HFORTOOLD3PD_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include "TROOT.h"
#include <TParticle.h>
#include <TLorentzVector.h>

class HforToolD3PD {
 
 private:

  enum momentumSet {PtEtaPhiM=0, PtEtaPhiE, PxPyPzM, PxPyPzE};

 public:

  enum verbosity {DEBUG=0, INFO, WARNING, ERROR, FATAL};
  enum verbosity m_verbosity;

  enum removalmode {DEFAULT=0, ALL, BBONLY};

  HforToolD3PD();
  virtual ~HforToolD3PD();

  virtual void initialize(const std::string schema="angularbased",
			  const double matchingcone=0.4);
  void setVerbosity(const enum verbosity ver);

  // function to retrieve the decision
  // for datasetNumber use mc_channel_number in rel.17, use RunNumber in rel.16
  int getDecision(int datasetNumber, const int mc_n, 
		  const std::vector<float> *mc_pt,
		  const std::vector<float> *mc_eta,
		  const std::vector<float> *mc_phi,
		  const std::vector<float> *mc_m,
		  const std::vector<int> *mc_pdgId,
		  const std::vector<int> *mc_status,
		  const std::vector<int> *mc_vx_barcode,
		  const std::vector<std::vector<int> > *mc_parent_index,
		  const std::vector<std::vector<int> > *mc_child_index,
		  const enum removalmode mode=DEFAULT, const bool update=true);
   
 private:

  // functions to retrieves the b/c quarks
  const std::vector<TParticle> & get_bQuarks_MPI();
  const std::vector<TParticle> & get_bQuarks_GS();
  const std::vector<TParticle> & get_bQuarks_ME();
  const std::vector<TParticle> & get_bQuarks_MEin();
  const std::vector<TParticle> & get_bQuarks_PDF();
  const std::vector<TParticle> & get_bQuarks_unknown();
  const std::vector<TParticle> & get_bQuarks_topdecay();
  const std::vector<TParticle> & get_cQuarks_MPI();
  const std::vector<TParticle> & get_cQuarks_GS();
  const std::vector<TParticle> & get_cQuarks_ME();
  const std::vector<TParticle> & get_cQuarks_MEin();
  const std::vector<TParticle> & get_cQuarks_PDF();
  const std::vector<TParticle> & get_cQuarks_unknown();
  const std::vector<TParticle> & get_cQuarks_topdecay();

  // List of four-vectors for b/c quarks from MPI / gluon splitting /
  // MatrixElement (processed by parton shower) / MatrixElement (not processed) /
  // unknonwn origins
  // The int key is the absolute pdgId
  std::map< int,std::vector<TParticle> > m_Quarks_MPI;
  std::map< int,std::vector<TParticle> > m_Quarks_GS;
  std::map< int,std::vector<TParticle> > m_Quarks_ME;
  std::map< int,std::vector<TParticle> > m_Quarks_MEin;
  std::map< int,std::vector<TParticle> > m_Quarks_topdecay;
  std::map< int,std::vector<TParticle> > m_Quarks_PDF;
  std::map< int,std::vector<TParticle> > m_Quarks_unknown;

  void angularBasedRemoval(const enum removalmode mode=DEFAULT);
  void keepAllRemoval();

  std::string m_schema;
  double      m_matchingcone;
  
  // variables to classfy the sample
  std::string m_sampleType;
  bool        m_isZinclusive;

  // variable to keep the decision
  std::string m_result;

  // name of the genertor used for showeing (PYTHIA or HERWIG)
  std::string m_showerGenerator;

  // function to loop over the generated event record, find the b/c quarls,
  // and se the corresponding lists of four-vectors
  void findHFQuarks(const int mc_n, 
		    const std::vector<float> *mc_p1,
		    const std::vector<float> *mc_p2,
		    const std::vector<float> *mc_p3,
		    const std::vector<float> *mc_p4,
		    const std::vector<int> *mc_pdgId,
		    const std::vector<int> *mc_status,
		    const std::vector<int> *mc_vx_barcode,
		    const std::vector<std::vector<int> > *mc_parent_index,
		    const std::vector<std::vector<int> > *mc_child_index,
		    const enum momentumSet mode);
  
  void findHFQuarksHerwig(const std::map< int,std::vector<TParticle> > &finalstate_q,
			  const int mc_n,
			  const std::vector<int> *mc_pdgId,
			  const std::vector<int> *mc_status,
			  const std::vector<std::vector<int> > *mc_parent_index);
  
  void findHFQuarksPythia(const std::map< int,std::vector<TParticle> > &finalstate_q,
			  const int mc_n,
			  const std::vector<int> *mc_pdgId,
			  const std::vector<int> *mc_status,
			  const std::vector<int> *mc_vx_barcode,
			  const std::vector<std::vector<int> > *mc_parent_index);

  bool checkSampleType(int RunNumber);

  int matchdR(std::vector<TParticle>* quarks);

  inline double deltaR(const TParticle &v1, const TParticle &v2) {
    double dphi = std::fabs(v1.Phi() - v2.Phi());
    dphi = (dphi<=M_PI)? dphi : 2*M_PI-dphi;
    double deta = std::fabs(v1.Eta() - v2.Eta());
    return std::sqrt(dphi*dphi + deta*deta);
  }


};


#endif
