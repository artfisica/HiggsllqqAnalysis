#ifndef __CommonTools_h__
#define __CommonTools_h__

#define Z_pdg_mass 91187.6
#define mu_pdg_mass 105.65836
#define el_pdg_mass 0.510998910
#define Z_pdg_width 2495.2

#include "egammaAnalysisUtils/EnergyRescaler.h"

#include "HiggsAnalysis/HSG2EventReader.h"
#include "HiggsllqqAnalysis/ChargedLepton.h"
#include <TLorentzVector.h>
//#include "CLHEP/Matrix/Matrix.h"
#include <iostream>

class Matrix;

namespace Muon {
   enum {
      STACO,
      MUID,
      CALO,
   };
};

namespace Electron {
   enum {
      GSF,
      noGSF,
   };
};

namespace CommonTools {
   TLorentzVector getVector(D3PDReader::TruthParticleD3PDObjectElement *p);
   Float_t getBremFitDp(D3PDReader::ElectronD3PDObjectElement *el);
   Float_t getEnergyUncertainty(eg2011::EnergyRescaler *rescaler, Analysis::ChargedLepton *lep);
   //CLHEP::HepMatrix getCovarianceMatrixMuon(Analysis::ChargedLepton *lep);
   //CLHEP::HepMatrix getCovarianceMatrixElectron(eg2011::EnergyRescaler *rescaler, Analysis::ChargedLepton *lep);
   //CLHEP::HepMatrix getCovarianceMatrix(eg2011::EnergyRescaler *rescaler, Analysis::ChargedLepton *lep);
};

#endif
