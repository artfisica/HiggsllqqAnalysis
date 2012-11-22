/**
   @class    HforToolD3PD.h

   @brief    tool for removal of heavy flavor overlap in Alpgen samples

   @author   Dominic Hirschbuehl, Michel Sanders (original code in athena)
   @author   Takashi Yamanaka (converted the code for D3PD)

   taken from  Generators/GenAnalysisTools/TruthUtils/tags/tags/TruthUtils-00-01-00

   30 March 2012
*/


#include "HiggsllqqAnalysis/HforToolD3PD.h"
#include <iostream>

HforToolD3PD::HforToolD3PD()
{
  initialize("angularbased", 0.4);
  m_verbosity = FATAL;
}


void HforToolD3PD::initialize(const std::string schema, const double matchingcone) 
{
  if (m_verbosity<=INFO) std::cout << "in initialize()" << std::endl;

  m_schema = schema;
  m_matchingcone = matchingcone;
}

void HforToolD3PD::setVerbosity(const enum verbosity ver)
{
  m_verbosity = ver;
}

///////////////////////////////////////////////////////////////////////
// return the q-quarks labeled as MPI, GS, ME, incoming ME or unknown
///////////////////////////////////////////////////////////////////////
const std::vector<TParticle> & HforToolD3PD::get_bQuarks_MPI() {
  return m_Quarks_MPI[5];
}

const std::vector<TParticle> & HforToolD3PD::get_bQuarks_GS() {
  return m_Quarks_GS[5];
}

const std::vector<TParticle> & HforToolD3PD::get_bQuarks_ME() {
  return m_Quarks_ME[5];
}

const std::vector<TParticle> & HforToolD3PD::get_bQuarks_MEin() {
  return m_Quarks_MEin[5];
}

const std::vector<TParticle> & HforToolD3PD::get_bQuarks_unknown() {
  return m_Quarks_unknown[5];
}

const std::vector<TParticle> & HforToolD3PD::get_bQuarks_topdecay() {
  return m_Quarks_topdecay[5];
}

const std::vector<TParticle> & HforToolD3PD::get_bQuarks_PDF() {
  return m_Quarks_PDF[5];
}

const std::vector<TParticle> & HforToolD3PD::get_cQuarks_MPI() {
  return m_Quarks_MPI[4];
}

const std::vector<TParticle> & HforToolD3PD::get_cQuarks_GS() {
  return m_Quarks_GS[4];
}

const std::vector<TParticle> & HforToolD3PD::get_cQuarks_ME() {
  return m_Quarks_ME[4];
}

const std::vector<TParticle> & HforToolD3PD::get_cQuarks_MEin() {
  return m_Quarks_MEin[4];
}

const std::vector<TParticle> & HforToolD3PD::get_cQuarks_unknown() {
  return m_Quarks_unknown[4];
}

const std::vector<TParticle> & HforToolD3PD::get_cQuarks_topdecay() {
  return m_Quarks_topdecay[4];
}

const std::vector<TParticle> & HforToolD3PD::get_cQuarks_PDF() {
  return m_Quarks_PDF[4];
}



///////////////////////////////////////////////////////////////
int HforToolD3PD::getDecision(const int datasetNumber,
			      const int mc_n, 
			      const std::vector<float> *mc_pt,
			      const std::vector<float> *mc_eta,
			      const std::vector<float> *mc_phi,
			      const std::vector<float> *mc_m,
			      const std::vector<int> *mc_pdgId,
			      const std::vector<int> *mc_status,
			      const std::vector<int> *mc_vx_barcode,
			      const std::vector<std::vector<int> > *mc_parent_index,
			      const std::vector<std::vector<int> > *mc_child_index,
			      const enum removalmode mode,
			      const bool update)
///////////////////////////////////////////////////////////////
{
  if (m_verbosity<=INFO) std::cout << "in getDecision()" << std::endl;
  
  if (!checkSampleType(datasetNumber)) {
    if(m_verbosity<=WARNING) std::cout << "This dataset does not need HforTool" << std::endl;
    return -1;
  }
  
  if (update) {
    findHFQuarks(mc_n, mc_pt, mc_eta, mc_phi, mc_m, mc_pdgId, mc_status, mc_vx_barcode, mc_parent_index, mc_child_index, PtEtaPhiM);
    
    if (m_verbosity<=INFO) {
      std::cout << "GS(b) = " << m_Quarks_GS[5].size() 
		<< " ME(b) = " << m_Quarks_ME[5].size() 
		<< " MPI(b) = " << m_Quarks_MPI[5].size() 
		<< " TOP(b) = " << m_Quarks_topdecay[5].size() << std::endl;
      std::cout << "GS(c) = " << m_Quarks_GS[4].size() 
		<< " ME(c) = " << m_Quarks_ME[4].size()
		<< " MPI(c) = " << m_Quarks_MPI[4].size() 
		<< " TOP(c) = " << m_Quarks_topdecay[4].size() << std::endl;
    }
  }

  if(m_schema == "angularbased")
    angularBasedRemoval(mode);
  else 
    keepAllRemoval();

  if( m_verbosity<=INFO ) std::cout << "result = " << m_result << std::endl;

  if (m_result == "isBB")               return 0;
  else if (m_result == "isCC")          return 1;
  else if (m_result == "isC")           return 2;
  else if (m_result == "isLightFlavor") return 3;
  else if (m_result == "kill")          return 4;
  else                                  return -1;  
}


/////////////////////////////////////////////////
// Find the heavy flavor quarks in this event
void HforToolD3PD::findHFQuarks(const int mc_n, 
				const std::vector<float> *mc_p1,
				const std::vector<float> *mc_p2,
				const std::vector<float> *mc_p3,
				const std::vector<float> *mc_p4,
				const std::vector<int> *mc_pdgId,
				const std::vector<int> *mc_status,
				const std::vector<int> *mc_vx_barcode,
				const std::vector<std::vector<int> > *mc_parent_index,
				const std::vector<std::vector<int> > *mc_child_index,
				const enum momentumSet mode)
/////////////////////////////////////////////////
{
  if(m_verbosity<=DEBUG) std::cout << "in findHFQuarks()" << std::endl;

  // clean up the previous event
  m_Quarks_MPI.clear();
  m_Quarks_GS.clear();
  m_Quarks_ME.clear();
  m_Quarks_MEin.clear();
  m_Quarks_topdecay.clear();
  m_Quarks_PDF.clear();
  m_Quarks_unknown.clear();

  if(m_verbosity<=DEBUG) std::cout << "m_Quarks cleared" << std::endl;

  // vectors with the intial and final state b/c quarks: 
  // ie. initial or final in the parton shower:
  // ignore b/c quarks from b/c-hadron decays
  std::map< int,std::vector<TParticle> > finalstate_q;

  if(m_verbosity<=DEBUG) std::cout << "mc_n=" << mc_n << std::endl;

  for (int i=0; i<mc_n; i++) {
    int pdg = mc_pdgId->at(i);
    int apdg = std::abs(pdg);
    int status = mc_status->at(i);
    int parent[2] = {-1, -1};
    int child[2] = {-1, -1};
    if ( apdg==5 || apdg==4 ) { // b or c quark
      
      if(m_verbosity<=DEBUG) std::cout << "pdg = " << pdg << ": " << i << " : stat " << status << " ( " << mc_p1->at(i) << " " << mc_p2->at(i) << " " << mc_p3->at(i) << " " << mc_p4->at(i) << " )" << std::endl;

      // find the production vertex and parents
      bool hasbchadronparent(false);
      bool hasmpiparent(false);
      bool hastopparent(false);
      bool hasWparent(false);
      if (mc_parent_index->at(i).size()) {

	if(m_verbosity<=DEBUG) std::cout << "  parents size: " << mc_parent_index->at(i).size() << std::endl;

	// check that there is no b/c-hadron as a parent, or a top quark as a parent
	// also find mpi parent
	for(unsigned int j=0; j<mc_parent_index->at(i).size(); j++) {
	  int parentId = mc_parent_index->at(i).at(j);
	  if(j<2) parent[j] = parentId;
     
	  if(m_verbosity<=DEBUG) std::cout << "   incoming: " << parentId << std::endl;

	  if (parentId < mc_n) {
	    int pdgin = std::abs(mc_pdgId->at( parentId ));
	    
	    if( pdgin==6 ) {
	      hastopparent = true;
	      if(m_verbosity<=DEBUG) std::cout << " b/c parton with a top parent" << std::endl;
	    }
	    if( pdgin==24 ) {
	      hasWparent = true;
	      if(m_verbosity<=DEBUG) std::cout << " b/c parton with a W parent" << std::endl;
	    }
	    if( (pdgin%10000)/1000 == apdg || (pdgin%1000)/100 == apdg )
	      hasbchadronparent = true;
	    
	    // also reject the c-quarks 
	    if ( apdg == 4 && ( pdgin == 5 || (pdgin%10000)/1000 == 5 ||
				(pdgin%1000)/100 == 5 ) )
	      hasbchadronparent = true;
	    if ( pdgin == 0 && mc_status->at( mc_parent_index->at(i).at(j) ) == 120 )
	      hasmpiparent = true;
	  } else {
	    if(m_verbosity<=DEBUG) std::cout << "  parent particle is out of range" << std::endl;
	  }
	}
      } else {
	// b/c parton without production vertex
	if(m_verbosity<=DEBUG) std::cout << "  b/c parton without production vertex" << std::endl;
      }

      if( hasbchadronparent ) {
	if(m_verbosity<=DEBUG) std::cout << " b/c parton with a b/c quark/hadron parent" << std::endl;
      }

      // find the decay vertex and children
      bool hasbcquarkdaughter(false);
      if ( mc_child_index->at(i).size() ) {
	
	if(m_verbosity<=DEBUG) std::cout << "  children size: " << mc_child_index->at(i).size() << std::endl;

	// check whether there are only non b/c-quark daughters
	for(unsigned int j=0; j<mc_child_index->at(i).size(); j++) {
	  int childId = mc_child_index->at(i).at(j); 
	  if(j<2) child[j] = childId;
	  if(m_verbosity<=DEBUG) std::cout << "   outgoing: " << childId << std::endl;

	  if (childId<mc_n) {
	    int pdgout = std::abs(mc_pdgId->at( childId ));
	    if ( pdgout==apdg )
	      hasbcquarkdaughter = true;
	  } else {
	    if(m_verbosity<=DEBUG) std::cout << "  child particle is out of range" << std::endl;
	  }
	}
      } else {
	// b/c parton without decay vertex
	if(m_verbosity<=DEBUG) std::cout << "  b/c parton without decay vertex" << std::endl;
      }

      TLorentzVector part_hlv;
      if (mode==PtEtaPhiM) {
	part_hlv.SetPtEtaPhiM(mc_p1->at(i), mc_p2->at(i), mc_p3->at(i), mc_p4->at(i));
      } else if (mode==PtEtaPhiE) {
	part_hlv.SetPtEtaPhiE(mc_p1->at(i), mc_p2->at(i), mc_p3->at(i), mc_p4->at(i));
      } else if (mode==PxPyPzM) {
	part_hlv.SetXYZM(mc_p1->at(i), mc_p2->at(i), mc_p3->at(i), mc_p4->at(i));
      } else if (mode==PxPyPzE) {
	part_hlv.SetPxPyPzE(mc_p1->at(i), mc_p2->at(i), mc_p3->at(i), mc_p4->at(i));
      } else {
	std::cerr << "unknown momentum set." << std::endl;
      }

      TParticle bcpart(pdg, i, parent[0], parent[1], child[0], child[1], part_hlv.Px(), part_hlv.Py(), part_hlv.Pz(), part_hlv.E(), 0., 0., 0., 0.);
     
      if (!hasbchadronparent && !hasbcquarkdaughter) {
	if(m_verbosity<=DEBUG) std::cout << " final state b/c quark, status=" << status << std::endl;
	// if no b/c-hadron parent and no b/c-quark daughter. keep it!
	finalstate_q[apdg].push_back(bcpart);
      }

      // if no b/c-hadron parent. check it to see whether it comes from the ME
      // but ignore the ones with an MPI parent
      if (!hasbchadronparent && !hasmpiparent &&
	  (status == 123 || status == 124 ) ) {
	if(m_verbosity<=DEBUG) std::cout << "  b/c-quark from ME" << std::endl;
	m_Quarks_MEin[apdg].push_back(bcpart);
      }
   
    } // particle is b or c quark

  } // loop over all particles in the event

  if(m_showerGenerator == "HERWIG")
    findHFQuarksHerwig(finalstate_q, mc_n, mc_pdgId, mc_status, mc_parent_index);
  else if(m_showerGenerator == "PYTHIA")
    findHFQuarksPythia(finalstate_q, mc_n, mc_pdgId, mc_status, mc_vx_barcode, mc_parent_index);
  else{
    // should never get here 
    if(m_verbosity<=WARNING) std::cout << "Shower generator type unknown" << std::endl;
  }
}



void HforToolD3PD::findHFQuarksHerwig(const std::map< int,std::vector<TParticle> > &finalstate_q,
				      const int mc_n,
				      const std::vector<int> *mc_pdgId,
				      const std::vector<int> *mc_status,
				      const std::vector<std::vector<int> > *mc_parent_index)
{
  if(m_verbosity<=DEBUG) std::cout << "findHFQuarksHerwig" << std::endl;

  // loop over all the final state b/c-quarks and find out where they come from
  // first loop over quarks flavors that were stored (b,c)
  for( std::map< int, std::vector<TParticle> >::const_iterator ipdg = finalstate_q.begin(); 
       ipdg != finalstate_q.end(); ipdg++ ) {
    int apdg(ipdg->first);
    if(m_verbosity<=DEBUG) std::cout << "looking for ancestors of pdg " << apdg << std::endl;

    // second loop over the final state quarks
    for( std::vector<TParticle>::const_iterator ibcpart = (ipdg->second).begin() ;
	 ibcpart != (ipdg->second).end(); ibcpart++) {
      const TParticle bcpart(*ibcpart);
      if(m_verbosity<=DEBUG) std::cout << "final state b/c " << bcpart.GetStatusCode() << std::endl;
      bool isMPI(false);
      bool isGS(false);
      bool isME(false);
      bool isPDF(false);
      bool isTopDecay(false);
      bool isWDecay(false); // subset of top-decays, for hadronic top-decays
      bool iscquarkfromb(false);
      // if the final state quark is a PDF parton, ignore it
      // (in AOD, descendants of these partons may be filtered out)
      if ( mc_status->at(bcpart.GetStatusCode()) == 141 || mc_status->at(bcpart.GetStatusCode()) == 142 ) {
	if(m_verbosity<=DEBUG) std::cout << "PDF !!" << std::endl;
	isPDF = true;
      }
      if ( !isPDF && bcpart.GetMother(0) != -1 ) {
	std::vector< std::vector<int> > ancestor_index;
	for(int i=0; i<2 && bcpart.GetMother(i)!=-1 && bcpart.GetMother(i)<mc_n; i++) {
	  std::vector<int> temp_vector;
	  temp_vector.push_back( bcpart.GetMother(i) );
	  ancestor_index.push_back( temp_vector );
	}
	for (unsigned int i=0; i<ancestor_index.size(); i++) {
	  std::vector<int> temp_vector;
	  for (unsigned int j=0; j<ancestor_index.at(i).size(); j++) {
	    for (unsigned int k=0; k<mc_parent_index->at(ancestor_index.at(i).at(j)).size(); k++) {
	      temp_vector.push_back( mc_parent_index->at(ancestor_index.at(i).at(j)).at(k) );
	    }
	  }
	  if(temp_vector.size()){
	    ancestor_index.push_back( temp_vector );
	  }
	}
 
	for(unsigned int i=0; i<ancestor_index.size(); i++) {
	  for(unsigned int j=0; j<ancestor_index.at(i).size(); j++) {
	    int apdgin = std::abs(mc_pdgId->at( ancestor_index.at(i).at(j) ));
	    int statusin = mc_status->at( ancestor_index.at(i).at(j) );
	    if (apdgin != apdg) {
	      // if MPI as a non-b parent, label it
	      if ( apdgin == 0 && statusin == 120 ) {
		isMPI = true;
	      }
	      // gluon splitting or ME origin: in evgen files.
	      // proton <id::2212> seem to be save in all events: not so in  AOD files
	      // Thus look for non-HF origin with status 121 or 122
	      if ( statusin == 121 || statusin == 122 ) {
		isGS = true;
	      }
	      // c quarks froma b quark (in b-hadron decays)
	      if (apdg == 4 && (apdgin == 5 || (apdgin&10000)/1000 == 5 ||
				(apdgin%1000)/100 == 5 ) ) {
		iscquarkfromb = true;
	      }
	      // b quark from a b-hadron decay
	      // b directly from b-hadron already rejected
	      if (apdg == 5 && ( (apdgin&10000)/1000 == 5 || 
				 (apdgin%1000)/100 == 5) ) {
		if(m_verbosity<=DEBUG) std::cout << "  b quark from b-hadron" << std::endl;
	      }
	      // top quark decay
	      if ( apdgin == 6 ) {
		if(m_verbosity<=DEBUG) std::cout << "  TOP !!" << std::endl;
		isTopDecay = true;
	      }
	      // W boson decay
	      if ( apdgin == 24 ) {
		if(m_verbosity<=DEBUG) std::cout << "  W !!" << std::endl;
		isWDecay = true;
	      }
	    } else {
	      if (m_verbosity<=DEBUG) std::cout << " b/c parent" << std::endl;
	      // if the status of a b-quark is 123 or 124, then it is a ME b-quark
	      if ( statusin == 123 || statusin == 124 ) {
		if (m_verbosity<=DEBUG) std::cout << "  ME !!" << std::endl;
		isME = true;
	      }
	      // if status 141 or 142 then it came from the PDF. ignore them!!
	      if ( statusin == 141 || statusin == 142 ) {
		if (m_verbosity<=DEBUG) std::cout << "  PDF !!" << std::endl;
		isPDF = true;
	      }
	    } // b/c or non-b/c quark
	  }
	} // loop over all ancestors
      } // final state b/c-quark with a production vertex

      // TopDecay does not depend on status code so it comes first
      // MPI output is also status 123,124 so MPI comes before anything else
      // ME parents have status 121 or 122 so ME comes before GS
      if ( !iscquarkfromb && !isPDF ) {
	if ( isTopDecay ) {
	  // If a b or c appers in the shower of a b-parton from a top decay,
	  // it should be counted as coming from the top decay.
	  // However, the b or c should not come from a W boson in a top decay
	  if ( !isWDecay ) {
	    m_Quarks_topdecay[apdg].push_back( bcpart );
	  }
	} else if ( isMPI ) {
	  if(m_verbosity<=DEBUG) std::cout << " come from MPI" << std::endl;
	  m_Quarks_MPI[apdg].push_back( bcpart );
	} else if ( isME ) {
	  if (m_verbosity<=DEBUG) std::cout << " come from ME" << std::endl;
	  m_Quarks_ME[apdg].push_back( bcpart );
	} else if ( isGS ) {
	  // in AOD, incoming ME partons may look like GS partons 
	  // if their descendants are filltered out
	  if ( !(mc_status->at(bcpart.GetStatusCode()) == 123 || mc_status->at(bcpart.GetStatusCode()) == 124) ) {
	    if (m_verbosity<=DEBUG) std::cout << " come from GS" << std::endl;
	    m_Quarks_GS[apdg].push_back( bcpart );
	  } else {
	    if(m_verbosity<=DEBUG) std::cout << "ME b/c-quark identified as GS" << std::endl;
	  }
	} 
	else {
	  if(m_verbosity<=DEBUG) std::cout << " Unidentified b/c-quark" << std::endl;
	  m_Quarks_unknown[apdg].push_back( bcpart );
	}
      } // not a c-quark from a b decay or a PDF c-quark
      else if ( !iscquarkfromb && isPDF ) {
	m_Quarks_PDF[apdg].push_back( bcpart );
      }
    
    } // loop over final state b/c-quarks
  } // loop over quark flavors
  if(m_verbosity<=DEBUG) std::cout << "Loop over quaurk flavors " << std::endl;

} // end of HforToolD3PD::findHFQuarksHerwig()
      
	  
///////////////////////////////////////////////////////
// Find the specifics for HF quarks in Pythia shower
void HforToolD3PD::findHFQuarksPythia(const std::map< int,std::vector<TParticle> > &finalstate_q,
				      const int mc_n,
				      const std::vector<int> *mc_pdgId,
				      const std::vector<int> *mc_status,
				      const std::vector<int> *mc_vx_barcode,
				      const std::vector<std::vector<int> > *mc_parent_index)
{
  if(m_verbosity<=DEBUG) std::cout << "findHFQuarksPythia" << std::endl;

  // loop over all the final state b/c-quarks and find out where they come from
  // first loop over quarks flavors that were stored (b,c)
  for ( std::map< int, std::vector<TParticle> >::const_iterator ipdg = finalstate_q.begin();
	ipdg != finalstate_q.end(); ipdg++) {
    int apdg(ipdg->first);
    if(m_verbosity<=DEBUG) std::cout << "looking for ancestors of pdg " << apdg << std::endl;

    // assume the Alpgen input (initial state -> final-stae) is completely
    // included in the event record with stat=3 particles 
    // the partons that we need are *not* these ones, but if 
    // these stat=3 partons exist, then the correct partons are in the event too
    std::vector<TParticle> MEParton;
    std::vector<TParticle> PDFParton;

    // loop over the stat=3 final state quarks
    for (std::vector<TParticle>::const_iterator ibcpart = (ipdg->second).begin();
	 ibcpart != (ipdg->second).end(); ibcpart++) {
      TParticle bcpart(*ibcpart);
      if ( mc_status->at(bcpart.GetStatusCode()) == 3 ) {
	if ( m_verbosity<=DEBUG ) std::cout << "final state b/c (stat=3) " << bcpart.GetStatusCode() << ", m = " << bcpart.GetMass() << std::endl;
	// if this parton has no descendants, then it's a ME parton
	if ( bcpart.GetDaughter(0)==-1 ) {
	  if(m_verbosity<=DEBUG) std::cout << "  ME parton" << std::endl;
	  m_Quarks_MEin[apdg].push_back( bcpart );
	  MEParton.push_back( bcpart );
	  // if there is a direct stat=3 ancestor with the same flavor,
	  // then there is PDF parton too (eg. qc->q'Wc)
	  if ( bcpart.GetMother(0) != -1 ) {
	    for ( int i=0; i<2 && bcpart.GetMother(i)!=-1 && bcpart.GetMother(i)<mc_n; i++) {
	      if ( m_verbosity<=DEBUG ) std::cout << "    incoming: " << bcpart.GetMother(i) << std::endl;
	      int pdgin(mc_pdgId->at(bcpart.GetMother(i)));
	      TParticle pin;
	      if(abs(pdgin)==apdg && mc_status->at(bcpart.GetMother(i))==3) {
		if(m_verbosity<=DEBUG) std::cout << "  PDF parton" << std::endl;
		PDFParton.push_back(pin);
	      }
	    } // ME parton with a production vertex
	  } // ME final state parton
	} // final state parton with no decay vertex
	// else it's a PDF parton that gets annihilated
	else {
	  if(m_verbosity<=DEBUG) std::cout << "  PDF parton" << std::endl;
	  PDFParton.push_back( bcpart );
	}
	
      } // stat=3 final state parton
    } // first loop over final state quarks

    int nMEPartons(MEParton.size());
    int nPDFPartons(PDFParton.size());
    
    // loop over the other final state quarks
    for ( std::vector<TParticle>::const_iterator ibcpart = (ipdg->second).begin(); ibcpart != (ipdg->second).end(); ibcpart++ ) {
      TParticle bcpart(*ibcpart);
      if ( mc_status->at(bcpart.GetStatusCode()) != 3 ) {
	if(m_verbosity<=DEBUG) std::cout << "final state b/c " << bcpart.GetStatusCode() << ", m = " << bcpart.GetMass() << std::endl;
	bool isTopDecay(false);
	bool isWDecay(false); // subset of top-decays, for hadronic top-decays
	bool iscquarkfromb(false);
	
	// separate GS/ME from MPI by looking at the ancestors:
	// one and only one proton ancestor -> MPI
	// two ancestors of which one a proton -> GS/ME
	//	bool hasPAncestor(false);
	bool hasPAncestor(false);
	int nAncestors(0);
	// if the parton has an origin, look at the ancestors
	if ( bcpart.GetMother(0) != -1 ) {
	  // check whether there is a proton ancestor,
	  // and how many ancestors there are
	  std::vector< std::vector<int> > ancestor_index;
	  for( int i=0; i<2 && bcpart.GetMother(i) != -1 && bcpart.GetMother(i)<mc_n; i++ ) {
	    std::vector<int> temp_vector;
	    temp_vector.push_back( bcpart.GetMother(i) );
	    ancestor_index.push_back( temp_vector );
	  }
	    
	  for ( unsigned int i=0; i<ancestor_index.size(); i++ ) {
	    std::vector<int> temp_vector;
	    for ( unsigned int j=0; j<ancestor_index.at(i).size(); j++) {
	      for (unsigned int k=0; k<mc_parent_index->at(ancestor_index.at(i).at(j)).size(); k++ ) {
		temp_vector.push_back( mc_parent_index->at(ancestor_index.at(i).at(j)).at(k) );
	      }
	    }
	    if(temp_vector.size()) {
	      ancestor_index.push_back( temp_vector );
	    }
	  }

	  for(unsigned int i=0; i<ancestor_index.size(); i++) {
	    for(unsigned int j=0; j<ancestor_index.at(i).size(); j++) {
	      int apdgin = std::abs( mc_pdgId->at(ancestor_index.at(i).at(j) ) );
	      if ( apdgin != apdg ) {
		if (m_verbosity<=DEBUG) std::cout << "  non b/c ancestor " << std::endl;
		// proton parent
		if ( apdgin == 2212 ) {
		  hasPAncestor = true;
		}
		// count number of ancestors
		nAncestors += 1;

		// c quark from a b quark (in b-hadron decays)
		if ( apdg == 4 && ( apdgin == 5 || (apdgin%10000)/1000 == 5 || (apdgin%1000)/100 == 5 ) ) {
		  if (m_verbosity<=DEBUG) std::cout << "  c quark from b quark or b hadron" << std::endl;
		  iscquarkfromb = true;
		}
		// b quark from a b-hadron decay 
		// (b directly from b-hadron already rejected)
		if ( apdg == 5 && ( (apdgin%10000)/1000 == 5 || (apdgin%1000)/100 == 5 ) ) {
		  if ( m_verbosity<=DEBUG ) std::cout << "  b quark from b hadron" << std::endl;
		  iscquarkfromb = true;
		}
		// top quark decay
		if ( apdgin == 6 ) {
		  if ( m_verbosity<=DEBUG ) std::cout << "  TOP !!" << std::endl;
		  isTopDecay = true;
		}
		// W boson decay
		if ( apdgin == 24 ) {
		  if ( m_verbosity<=DEBUG ) std::cout << "  W !!" << std::endl;
		  isWDecay = true;
		}
	      } else {
		if ( m_verbosity<=DEBUG ) std::cout << "  b/c or MC/PDF parent " << std::endl;
	      } // b/c or non-b/c quark as parent
	    }
	  } // loop over all ancestors
	
	  if ( m_verbosity<=DEBUG ) std::cout << "nAncstors = " << nAncestors << ", hasPAnacestor = " << hasPAncestor << std::endl;

	  if ( isTopDecay || isWDecay || iscquarkfromb ) {
	    if ( m_verbosity<=DEBUG ) std::cout << "topDecay or WDecay or cquarkfromb" << std::endl;
	    // M.P. Sanders 19.Mar.2012
	    // do nothing for the time being ==> DO NOT RUN THIS ON TTBAR EVENTS!
	  }
	  else if ( (nAncestors == 1 && hasPAncestor) || (nAncestors == 0 && !hasPAncestor) ) {
	    if ( m_verbosity<=DEBUG ) std::cout << "MPI parton" << std::endl;
	    m_Quarks_MPI[apdg].push_back( bcpart );
	  }
	  else if ( (nAncestors == 2 && hasPAncestor) || (nAncestors == 1 && !hasPAncestor) ) {
	    if ( m_verbosity<=DEBUG ) std::cout << "nAncestors == 2 && hasPAncestor" << std::endl;
	    // this can come from GS, ME or PDF
	    // exact matching not possible, e.g. c cbar -> s W cbar
	    // cbar in the final state can be ME or PDF (from c in initial state)
	    // or g c -> W s c cbar, cbar in final state?
	    // HF pair from ME and HF pair from GS: don't know which one in which
	    if ( (nMEPartons + nPDFPartons) == 0 ) {
	      if ( m_verbosity<=DEBUG ) std::cout << "GS parton" << std::endl;
	      m_Quarks_GS[apdg].push_back( bcpart );
	    }
	    else {
	      if ( m_verbosity<=DEBUG) std::cout << "  ME/PDf parton" << std::endl;
	      // need to check production vertices and compare them with
	      // prod. vertices of ME/PDF partons
	      int pdg(bcpart.GetPdgCode());

	      // vertex -3 / -4 have as output showered ME parton and link
	      // to vertex -5 with stat=3 ME parton; identical pdgid
	      // -> prod vtx of showered ME parton ==
	      //    prod vtx of one of the partons of stat=3 ME parton
		
	      // vertex -3 / -4 have as output showered PDF parton and link
	      // to vertex -5 with stat=5 PDF anti-parton; opposite pdgid
	      // -> prod vtx of showerd PDF parton ==
	      //    prod vtx of stat=3 PDF parton, and opposite pdgid

	      // first check that showered ME/PDF parton has prod. vtx -3 or -4
	      int pvtx34 = mc_vx_barcode->at(bcpart.GetStatusCode());
	      bool bc34( pvtx34==-3 || pvtx34==-4 );
	      if ( !bc34 ) {
		for(unsigned int i=0; i<ancestor_index.size(); i++) {
		  for ( unsigned int j=0; j<ancestor_index.at(i).size() && !bc34; j++ ) {
		    int bcpv(-1);
		    if ( mc_parent_index->at( ancestor_index.at(i).at(j) ).size() ) 
		      bcpv = mc_vx_barcode->at( ancestor_index.at(i).at(j) );

		    if ( mc_pdgId->at( ancestor_index.at(i).at(j) ) == pdg && (bcpv==-3 || bcpv==-4) ) {
		      pvtx34 = mc_vx_barcode->at(ancestor_index.at(i).at(j));
		      bc34 = true;
		    }
		  }
		} // loop over ancestors
	      } // prod-vetex has not barcode -3 or -4
	      if ( bc34 ) {
		// number of times this parton gets identified as ME or PDF
		int nid(0);
		// PDF parton
		for ( std::vector<TParticle>::const_iterator ipdf = PDFParton.begin(); ipdf != PDFParton.end(); ipdf++ ) {
		  if ( mc_vx_barcode->at(ipdf->GetStatusCode()) == pvtx34 && ipdf->GetPdgCode()==-pdg ) {
		    if ( m_verbosity<=DEBUG ) std::cout << "  -> PDF parton" << std::endl;
		    nid += 1;
		    m_Quarks_PDF[apdg].push_back( bcpart );
		  }
		} // loop over stat=3 PDF partons
		// ME parton
		bool isME(false);
		for ( std::vector<TParticle>::const_iterator ime = MEParton.begin(); !isME && ime != MEParton.end(); ime++ ) {
		  
		  if ( ime->GetMother(0) != -1 ) {
		    for ( int i=0; i<2 && ime->GetMother(i)!=-1 && i<mc_n; i++ ) {
		      if ( mc_vx_barcode->at(ime->GetMother(i)) == pvtx34 && ime->GetPdgCode() == pdg ) {
			if ( m_verbosity<=DEBUG ) std::cout << "  -> ME parton" << std::endl;
			nid += 1;
			isME = true;
			m_Quarks_ME[apdg].push_back( bcpart );
		      }
		    } // loop over partons
		  } // ME (stat=3) parton has a prod. vertex
		} // loop over ME (stat=3) partons	      
		if ( m_verbosity<=DEBUG ) std::cout << "this parton was id-d " << nid << " times" << std::endl;
		if ( !nid ) {  
		  if ( m_verbosity<=DEBUG ) std::cout << "  -> GS parton" << std::endl;
		  m_Quarks_GS[apdg].push_back( bcpart );
		} // not identified as ME or PDF -> GS
	      } // prod-vertex with barcode -3 or -4
	      else {
		if ( m_verbosity<=WARNING ) std::cout << "parton with no -3 or -4 prod. vertex" << std::endl;
		m_Quarks_unknown[apdg].push_back( bcpart );
	      }
	    } // event with ME or PDF stat=3 partons

	  } // parton with 2 ancestors of which one is proton
	  else {
	    if ( m_verbosity<=WARNING ) std::cout << "parton with strange set of ancestors" << std::endl;
	    m_Quarks_unknown[apdg].push_back( bcpart );
	  }  
	} // final state b/c-quark with a production vertex
	
	// if the final state parton has to be prod. vertex, then we
	// don't know what it is
	else {
	  if ( m_verbosity<=DEBUG ) std::cout << "b/c-quark without prod. vtx -> unindentified" << std::endl;
	  // m_Quarks_unknown[apdg].push_back( bcpart );
	  // in D3PD, MPI quark can loose prod. vtx
	  m_Quarks_MPI[apdg].push_back( bcpart );
	}
	
      } // stat!=3 final state b/c-quarks
      
    } // loop over final state b/c-quarks
    
    // print out wrongly identified PDF/GS partons
    if ( m_verbosity<=WARNING ) {
      if ( m_Quarks_PDF[apdg].size() != PDFParton.size() ) {
	if ( m_verbosity<=WARNING ) std::cout << "Mismatch in number of id-d PDF partons, apdg = " << apdg << std::endl;
	if ( m_verbosity<=DEBUG ) std::cout << "PDF partons from ME:" << std::endl;
	for ( std::vector<TParticle>::const_iterator ipdf = PDFParton.begin(); ipdf != PDFParton.end(); ipdf++ ) {
	  if ( m_verbosity<=DEBUG ) std::cout << ipdf->Px() << ", " << ipdf->Py() << ", " << ipdf->Pz() << " pT = " << ipdf->Pt()/1000. << ", eta = " << ipdf->Eta() << std::endl;
	}
	if ( m_verbosity<=DEBUG ) std::cout << "PDF partons identified:" << std::endl;;
	for ( std::vector<TParticle>::const_iterator iq = m_Quarks_PDF[apdg].begin(); iq!=m_Quarks_PDF[apdg].end(); iq++ ) {
	  if ( m_verbosity<=DEBUG ) std::cout << iq->Px() << ", " << iq->Py() << ", " << iq->Pz() << " pT = " << iq->Pt()/1000. << ", eta = " << iq->Eta() << std::endl;
	}
	
	// Could move the pdf parton with smaller eta to gs. but only if
	// some gs partons have been identified already
	// This is not fail-proof...
	
      }
    } // check number of identified PDF partons
    
  } // loop over quark flavors

} // end of HforToolD3PD::findHFQuarksPythia();


  
///////////////////////////////////////////////////////
// No overlap removal. only migration of events
// to be used only Wbb samples with phase space cuts
void HforToolD3PD::keepAllRemoval()
///////////////////////////////////////////////////////
{
  if ((m_Quarks_GS[5].size() + m_Quarks_ME[5].size())>0) 
    m_result = "isBB";
  else if ((m_Quarks_GS[4].size() + m_Quarks_ME[4].size())>0)
    m_result = "isCC";
  else 
    m_result = "isLightFlavor";
}




/////////////////////////////////////////
// Do the angular based removal
void HforToolD3PD::angularBasedRemoval(const enum removalmode mode)
// mode = DEFAULT : default overlap removal
// mode = ALL : overlap removal between light/cc/c/bb
// mode = BBONLY : overlap removal between light/bb only
//////////////////////////////////////////
{
  
  // container to keep matched quarks
  // 0:GS  1:ME
  bool hasCC[2] = {false, false};
  bool hasBB[2] = {false, false};

  int match_bGS = matchdR(&m_Quarks_GS[5]);
  if (match_bGS > 0) hasBB[0] = true;

  int match_bME = matchdR(&m_Quarks_ME[5]);
  if (match_bME > 0) hasBB[1] = true;

  int match_cGS = matchdR(&m_Quarks_GS[4]);
  if (match_cGS > 0) hasCC[0] = true;

  int match_cME = matchdR(&m_Quarks_ME[4]);
  if (match_cME > 0) hasCC[1] = true;


  // light flavor samples
  if (m_sampleType == "isLightFlavor") {
    m_result = "isLightFlavor";

    if ( (mode==DEFAULT && !m_isZinclusive)
	 || mode==ALL) {
      // remove ME HF
      if ( (m_Quarks_ME[5].size()>0) || (m_Quarks_ME[4].size()>0) ) {
	m_result = "kill";
      } else if (m_Quarks_GS[5].size()>0) {
	// remove unmatched HF from GS
	if (hasBB[0])
	  m_result = "isBB";
	else 
	  m_result = "kill";
      } else if (m_Quarks_GS[4].size()>0) {
	// remove unmatched HF from GS
	if (hasCC[0])
	  m_result = "isCC";
	else 
	  m_result = "kill";
      }
    } else {

      // ======================================
      // special case for the Z inclusive  sample
      // we only have to remove overlap with Zbb samples

      // remove ME HF
      if (m_Quarks_ME[5].size()>0) {
	m_result = "kill";
      } else if (m_Quarks_GS[5].size()>0) {
	// remove unmatched HF from GS
	if (hasBB[0])
	  m_result = "isBB";
	else 
	  m_result = "kill";
      } else if ( (m_Quarks_GS[4].size() + m_Quarks_ME[4].size()) > 0) {
	// nothing to remove in case of c quarks
	m_result = "isCC";
      }
      // ==========================================
    }
  }

  // cc samples
  if (m_sampleType == "isCC") {
    m_result = "isCC";
    
    // remove matched ME HF
    if (hasCC[1])
      m_result = "kill";
    else if (hasBB[0])
      m_result = "isBB";
    else if (m_Quarks_GS[5].size()>0)
      m_result = "kill";

  }
  
  // c samples
  if (m_sampleType == "isC") {
    m_result = "isC";
    
    // remove matched ME HF
    if (hasCC[1])
      m_result = "kill";
    else if (hasBB[0])
      m_result = "isBB";
    else if (m_Quarks_GS[5].size()>0)
      m_result = "kill";

  }
  
  // bb samples - we only "promote" events,
  // therefore not c-quarks have to be considered
  if (m_sampleType == "isBB") {
    m_result = "isBB";

    // remove matched ME HF
    if (hasBB[1])
      m_result = "kill";
  }

}



///////////////////////////////////////////////
// Perform DeltaR matching between two quarks
int HforToolD3PD::matchdR(std::vector<TParticle>* quarks)
///////////////////////////////////////////////
{
  int match = 0;
  if (quarks->size() > 1) {
    for (unsigned int i=0; i<quarks->size(); i++) {
      for (unsigned int j=i+1; j<quarks->size(); j++) {
	double dR = deltaR(quarks->at(i), quarks->at(j));
	if (m_verbosity<=DEBUG) std::cout << "deltaR( " << i << " , " << j << " )= " << dR << std::endl;
	if (dR < m_matchingcone) {
	  ++match;
	}
      }
    }
  }
  return match;
}


////////////////////////////////////////////////
// Check which sample we are running over
bool HforToolD3PD::checkSampleType(int datasetNumber)
////////////////////////////////////////////////
{
  if(m_verbosity<=DEBUG) std::cout << "in checkSampleType" << std::endl;

  m_isZinclusive = false;

  // W inclusive samples
  //   HERWIG
  if ( (datasetNumber >= 107680 && datasetNumber <= 107685) // enu
       || (datasetNumber >= 107690 && datasetNumber <= 107695) // munu
       || (datasetNumber >= 107700 && datasetNumber <= 107705) // taunu
       || (datasetNumber >= 144018 && datasetNumber <= 144020) // Np5_excl
       || (datasetNumber >= 144022 && datasetNumber <= 144024) // Np6
       || (datasetNumber >= 144196 && datasetNumber <= 144207) // susyfilt 
       ){ 
    m_sampleType = "isLightFlavor";
    m_showerGenerator = "HERWIG";
    return true;
  }

  //   PYTHIA
  if ( (datasetNumber >= 117680 && datasetNumber <= 117685) // enu 
       || (datasetNumber >= 117690 && datasetNumber <= 117695) // munu
       || (datasetNumber >= 117700 && datasetNumber <= 117705) // taunu
       ){
    m_sampleType = "isLightFlavor";
    m_showerGenerator = "PYTHIA";
    return true;
  }
 
  // Z inclusive samples
  //   HERWIG
  if ( (datasetNumber >= 107650 && datasetNumber <= 107655) // ee
       || (datasetNumber >= 107660 && datasetNumber <= 107665) // mumu
       || (datasetNumber >= 107670 && datasetNumber <= 107675) // tautau
       || (datasetNumber >= 107710 && datasetNumber <= 107715) // nunu
       || (datasetNumber == 144017) // nunuNp5_exl
       || (datasetNumber == 144021) // nunuNp6
       || (datasetNumber >= 144192 && datasetNumber <= 144195) // nunu_susyfilt
       || (datasetNumber >= 116250 && datasetNumber <= 116275) // DY Mll10to40
       ){ 
    m_sampleType = "isLightFlavor";
    m_isZinclusive = true;
    m_showerGenerator = "HERWIG";
    return true;
  }

  //   PYTHIA
  if ( (datasetNumber >= 117650 && datasetNumber <= 117655) // ee
       || (datasetNumber >= 117660 && datasetNumber <= 117665) // mumu
       || (datasetNumber >= 117670 && datasetNumber <= 117675) // tautau
       ){
    m_sampleType = "isLightFlavor";
    m_isZinclusive = true;
    m_showerGenerator = "PYTHIA";
     return true;  
   } 

  // Wc samples
  //   HERWIG
  if (datasetNumber >= 117288 && datasetNumber <= 117297) {
    m_sampleType = "isC";
    m_showerGenerator = "HERWIG";
    return true;
  }
  
  //   PYTHIA
  if (datasetNumber >= 126601 && datasetNumber <= 126605) {
    m_sampleType = "isC";
    m_showerGenerator = "PYTHIA";
    return true;
  }

  // Wcc samples
  //   HERWIG
  if (datasetNumber >= 117284 && datasetNumber <= 117287) {
    m_sampleType = "isCC";
    m_showerGenerator = "HERWIG";
    return true;
  }
  
  //   PYTHIA
  if (datasetNumber >= 126606 && datasetNumber <= 126609) {
  m_sampleType = "isCC";
    m_showerGenerator = "PYTHIA";
    return true;
  }
  
  // Wbb samples
  //   HERWIG
  if ( (datasetNumber >= 106280 && datasetNumber <= 106283)
       || (datasetNumber >= 107280 && datasetNumber <= 107283) ) {
    m_sampleType = "isBB";
    m_showerGenerator = "HERWIG";
    return true;
  }
  
  //   PYTHIA
  if( datasetNumber >= 126530 && datasetNumber <= 126533) {
    m_sampleType = "isBB";
    m_showerGenerator = "PYTHIA";
    return true;
  }

  // Zcc samples
  //   HERWIG
  if ( (datasetNumber >= 126414 && datasetNumber <= 126421) // ee or mumu
       ) {
    m_sampleType = "isCC";
    m_showerGenerator = "HERWIG";
    return true;
  }

  // Zbb samples
  //   HERWIG
  if ( (datasetNumber >= 109300 && datasetNumber <= 109313) // ee, mumu, tautau
       || (datasetNumber >= 118962 && datasetNumber <= 118965) // nunu
       || (datasetNumber >= 128130 && datasetNumber <= 128143) // DY Mll10to30
       ) {
    m_sampleType = "isBB";
    m_showerGenerator = "HERWIG";
    return true;
  }

  //   PYTHIA
  if ( datasetNumber >= 126560 && datasetNumber <= 126563 ) {
    m_sampleType = "isBB";
    m_showerGenerator = "PYTHIA";
    return true;  
  }

  // ttbar inclusive samples
  //   HERWIG
  if ( (datasetNumber >= 105890 && datasetNumber <= 105897) 
       || (datasetNumber >= 117887 && datasetNumber <= 117899) ){
    m_sampleType = "isLightFlavor";
    m_showerGenerator = "HERWIG";
    return true;
  }

  // ttbb sample
  //   HERWIG
  if ( datasetNumber == 116108 ) {
    m_sampleType = "isBB";
    m_showerGenerator = "HERWIG";
    return true;
  }

  // ttcc sample
  //  HERWIG
  if ( datasetNumber == 116109 ) {
    m_sampleType = "isCC";
    m_showerGenerator = "HERWIG";
    return true;
  }
  
  return false;

}
  


HforToolD3PD::~HforToolD3PD() {
}
