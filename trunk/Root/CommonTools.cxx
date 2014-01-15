#define __CommonTools_cxx__

#include "HiggsllqqAnalysis/CommonTools.h"

TLorentzVector CommonTools::getVector(D3PDReader::TruthParticleD3PDObjectElement *p)
{
   TLorentzVector result;

   result.SetPtEtaPhiM(p->pt(), p->eta(), p->phi(), p->m());

   return result;
}

Float_t CommonTools::getBremFitDp(D3PDReader::ElectronD3PDObjectElement *el)
{
   Float_t result(-1);

   if (el->author() == 1 || el->author() == 3) {
     Int_t index = -1;

     // loop over refitted tracks and take the index of the GSF one
     // see https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=152105 for details
     for (Int_t i = 0; i < el->refittedTrack_n(); i++) {
        if ((el->refittedTrack_author())[i] == 4)
           index = i;
     } // loop over refitted tracks
     if (index != -1) {
        result = 1. - (el->refittedTrack_qoverp())[index] / (el->refittedTrack_LMqoverp())[index];
     } // found GSF track
   } // good author

   return result;
}

Float_t CommonTools::getEnergyUncertainty(eg2011::EnergyRescaler *rescaler, Analysis::ChargedLepton *lep) {
   // new recommendation
   return lep->Get4Momentum()->E() * rescaler->resolution(lep->Get4Momentum()->E()/1000., lep->GetElectron()->cl_eta(), kTRUE);
 //return el->cl_E() * rescaler->resolution(el->cl_E()/1000., el->cl_eta(), kTRUE);
 /*

   // old way

   Float_t abseta = TMath::Abs(el->cl_eta());
   Float_t energy = el->cl_E();

   Float_t sampling(0.);
   Float_t constant(0.);

   if (abseta < 0.8) {
     sampling = 0.091;
   } else if (abseta < 1.37) {
     sampling = 0.036 + 0.130 * abseta;
   } else if (abseta < 1.52) {
     sampling = 0.27;
   } else if (abseta < 2.0) {
     sampling = 0.85 - 0.36 * abseta;
   } else if (abseta < 2.3) {
     sampling = 0.16;
   } else if (abseta < 2.5) {
     sampling = -1.05 + 0.52 * abseta;
   } else { // not optimized
     sampling = -1.05 + 0.52 * abseta;
   }

   // table 4, http://arxiv.org/pdf/1110.3174v1
   if (abseta < 1.37) {
     constant = 0.012;
   } else if (abseta < 2.47) {
     constant = 0.018;
   } else if (abseta < 3.2) {
     constant = 0.033;
   } else {
     constant = 0.049;
   }

   Float_t sigmaE_over_E = TMath::Sqrt(TMath::Power(constant, 2) + TMath::Power(sampling/TMath::Sqrt(energy/1000.), 2));

   return energy * sigmaE_over_E;
   */
}

HepMatrix CommonTools::getCovarianceMatrixMuon(Analysis::ChargedLepton *lep) {
   D3PDReader::MuonD3PDObjectElement *mu = lep->GetMuon();
   double P     = lep->Get4Momentum()->P();

   // get the covariance matrix
   HepMatrix covmatrix(5, 5, 0);
   covmatrix(1, 1) = mu->cov_d0_exPV();
   covmatrix(2, 2) = mu->cov_z0_exPV();
   covmatrix(3, 3) = mu->cov_phi_exPV();
   covmatrix(4, 4) = mu->cov_theta_exPV();
   covmatrix(5, 5) = mu->cov_qoverp_exPV();
   covmatrix(1, 2) = mu->cov_d0_z0_exPV();
   covmatrix(1, 3) = mu->cov_d0_phi_exPV();
   covmatrix(1, 4) = mu->cov_d0_theta_exPV();
   covmatrix(1, 5) = mu->cov_d0_qoverp_exPV();
   covmatrix(2, 3) = mu->cov_z0_phi_exPV();
   covmatrix(2, 4) = mu->cov_z0_theta_exPV();
   covmatrix(2, 5) = mu->cov_z0_qoverp_exPV();
   covmatrix(3, 4) = mu->cov_phi_theta_exPV();
   covmatrix(3, 5) = mu->cov_phi_qoverp_exPV();
   covmatrix(4, 5) = mu->cov_theta_qoverp_exPV();

   // symmetrize it
   for (int ii = 1; ii < 6; ii++)
      for (int jj = ii + 1; jj < 6; jj++)
         covmatrix(jj, ii) = covmatrix(ii, jj);

   // go from d0,z0,phi,theta,1/P --> d0,z0,phi,theta,P
   HepMatrix Jacobian0(5, 5, 1);
   Jacobian0(5, 5) = -P * P;

   return HepMatrix(Jacobian0.T() * covmatrix * Jacobian0);
}

HepMatrix CommonTools::getCovarianceMatrixElectron(eg2011::EnergyRescaler *rescaler, Analysis::ChargedLepton *lep) {
   D3PDReader::ElectronD3PDObjectElement *el = lep->GetElectron();

   // let's assume zero correlation among CALO and ID measurement
   // get the covariance matrix d0,z0,phi,theta,P
   HepMatrix covmatrix(5, 5, 0);
   covmatrix(1, 1) = el->trackcov_d0();
   covmatrix(2, 2) = el->trackcov_z0();
   covmatrix(3, 3) = el->trackcov_phi();
   covmatrix(4, 4) = el->trackcov_theta();
   covmatrix(5, 5) = TMath::Power(getEnergyUncertainty(rescaler, lep), 2);
   covmatrix(1, 2) = el->trackcov_d0_z0();
   covmatrix(1, 3) = el->trackcov_d0_phi();
   covmatrix(1, 4) = el->trackcov_d0_theta();
   covmatrix(1, 5) = 0;
   covmatrix(2, 3) = el->trackcov_z0_phi();
   covmatrix(2, 4) = el->trackcov_z0_theta();
   covmatrix(2, 5) = 0;
   covmatrix(3, 4) = el->trackcov_phi_theta();
   covmatrix(3, 5) = 0;
   covmatrix(4, 5) = 0;

   // symmetrize it
   for (int ii = 1; ii < 6; ii++)
      for (int jj = ii + 1; jj < 6; jj++)
         covmatrix(jj, ii) = covmatrix(ii, jj);

   return covmatrix;
}

HepMatrix CommonTools::getCovarianceMatrix(eg2011::EnergyRescaler *rescaler, Analysis::ChargedLepton *lep) {
  if (lep->flavor() == Analysis::ChargedLepton::MUON) {
    return getCovarianceMatrixMuon(lep);
  }
  else if (lep->flavor() == Analysis::ChargedLepton::ELECTRON) {
    return getCovarianceMatrixElectron(rescaler, lep);
  } else {
    std::cerr << "ERROR: calling CommonTools::getCovarianceMatrix with an unknown flavour lepton!" << std::endl;
    HepMatrix covmatrix(5, 5, 0);
    return covmatrix;
  }
}
