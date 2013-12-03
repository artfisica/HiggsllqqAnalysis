// author: Valerio Ippolito <valerio.ippolito@cern.ch>
// see HllqqAnalysis.h for details

#ifndef __KinematicFitter_C__
#define __KinematicFitter_C__

#include "HiggsllqqAnalysis/KinematicFitter.h"

class chiSquare {
public:
   float m_mass;
   float m_width;
   TLorentzVector m_particle1;
   TLorentzVector m_particle2;
   float m_width1;
   float m_width2;

   double operator()(const double *xx) const {
      const double pt1_fit = xx[0];
      const double pt2_fit = xx[1];

      TLorentzVector particle1_fit;
      particle1_fit.SetPtEtaPhiM(xx[0], m_particle1.Eta(), m_particle1.Phi(), m_particle1.M());
      TLorentzVector particle2_fit;
      particle2_fit.SetPtEtaPhiM(xx[1], m_particle2.Eta(), m_particle2.Phi(), m_particle2.M());

      const double m_pp = (particle1_fit + particle2_fit).M();
      const double mass_constraint = TMath::Power(((m_pp - m_mass) / m_width), 2);
      const double pt_constraint1  = TMath::Power(((pt1_fit - m_particle1.Pt()) / m_width1), 2);
      const double pt_constraint2  = TMath::Power(((pt2_fit - m_particle2.Pt()) / m_width2), 2);

      return mass_constraint + pt_constraint1 + pt_constraint2;

   }
};

KinematicFitter::KinematicFitter(int n, double mass, double width)
{
   setDim(n);
   m_mass = mass;
   m_width = width;
   m_particle.clear();

   m_min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
}

KinematicFitter::~KinematicFitter()
{
   m_particle.clear();
   delete m_min;
}

CandidatePair KinematicFitter::findBestPair()
{
   std::vector<CandidatePair> candidate;
   candidate.clear();

   int max_idx = TMath::Min(m_dim, (int)m_particle.size()); // allow to restrict to just the first m_dim particles

   for (int i = 0; i < max_idx; i++) {
      for (int j = i + 1; j < max_idx; j++) {
         CandidatePair thisCandidate = minimizeChiSquare(i, j);
         candidate.push_back(thisCandidate);
      }
   }

   int best_pair_index(-1);
   double best_pair_chisq(99999.);
   for (int i = 0; i < (int)candidate.size(); i++) {

      double this_chisq = (candidate.at(i)).chiSq;
      if (this_chisq < best_pair_chisq) {
         best_pair_chisq = this_chisq;
         best_pair_index = i;
      }
   }

   if (best_pair_index == -1) {
//      std::cout << "WARNING: no candidate found, returning dummy pair" << std::endl;
      CandidatePair dummy(-1, -1, 99999);
      return dummy;
   } else {
      TLorentzVector uno = m_particle.at(candidate.at(best_pair_index).index1);
      TLorentzVector due = m_particle.at(candidate.at(best_pair_index).index2);
      return candidate.at(best_pair_index);
   }
}

CandidatePair KinematicFitter::minimizeChiSquare(int i, int j)
{
   chiSquare theChiSquare;

   theChiSquare.m_particle1 = m_particle.at(i);
   theChiSquare.m_particle2 = m_particle.at(j);
   theChiSquare.m_mass = this->m_mass;
   theChiSquare.m_width = this->m_width;
   theChiSquare.m_width1 = getResolution(theChiSquare.m_particle1);
   theChiSquare.m_width2 = getResolution(theChiSquare.m_particle2);

   m_min->SetMaxFunctionCalls(1000000);
   m_min->SetMaxIterations(100000);
   m_min->SetTolerance(0.001);

   ROOT::Math::Functor f(theChiSquare, 2);
   double step[2] = {0.01, 0.01};
   double variable[2];
   variable[0] = (theChiSquare.m_particle1).Pt();
   variable[1] = (theChiSquare.m_particle2).Pt();

   m_min->SetFunction(f);

   // Set the free variables to be minimized!
   m_min->SetVariable(0, "pt1_fit", variable[0], step[0]);
   m_min->SetVariable(1, "pt2_fit", variable[1], step[1]);

   m_min->Minimize();

   const double *xs = m_min->X();

   TLorentzVector particle1_refit, particle2_refit;
   particle1_refit.SetPtEtaPhiM(xs[0], (theChiSquare.m_particle1).Eta(), (theChiSquare.m_particle1).Phi(), (theChiSquare.m_particle1).M());
   particle2_refit.SetPtEtaPhiM(xs[1], (theChiSquare.m_particle2).Eta(), (theChiSquare.m_particle2).Phi(), (theChiSquare.m_particle2).M());
   double refittedMass = (particle1_refit + particle2_refit).M();

   CandidatePair result(i, j, theChiSquare(xs), refittedMass, particle1_refit.Pt(), particle2_refit.Pt());

   return result;
}


#endif
