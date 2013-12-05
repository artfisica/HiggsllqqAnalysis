////////////////////////////////////////////////////////////
// KinematicFitter
//
// A kinematic fitter for generic pairs of particles
// coming from a decay.
//
// author: Valerio Ippolito <valerio.ippolito@cern.ch>
//
////////////////////////////////////////////////////////////


#ifndef __KinematicFitter_h__
#define __KinematicFitter_h__

#include "TROOT.h"
#include "Math/GSLMinimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>
#include <TLorentzVector.h>
#include "JetResolution/JERProvider.h"

class CandidatePair {
public:
   CandidatePair(int i, int j, double chisquare, double refitMass = -9999., double refitPt1 = -9999., double refitPt2 = -9999.) {
      index1 = i;
      index2 = j;
      chiSq = chisquare;
      refittedMass = refitMass;
      refittedPt1  = refitPt1;
      refittedPt2  = refitPt2;
   }
   ~CandidatePair() {};
   void Reset() {
     index1 = -1;
     index2 = -1;
     chiSq = -1.;
     refittedMass = -1.;
     refittedPt1 = -1.;
     refittedPt2 = -1.;
   }

   int index1;
   int index2;
   double chiSq;
   double refittedMass;
   double refittedPt1;
   double refittedPt2;
};

class KinematicFitter {
public:
   KinematicFitter(int n, double mass, double width);
   ~KinematicFitter();
   void setDim(int n) {
      m_dim = n;
   }
   int dim(int n) {
      return m_dim;
   }
   CandidatePair findBestPair();
   bool addParticles(std::vector<TLorentzVector> &particle_list) {
      m_particle = particle_list;
      return true;
   }
   bool addParticle(TLorentzVector &particle) {
      m_particle.push_back(particle);
      return true;
   }
   bool clearParticles() {
      m_particle.clear();
      return true;
   }
   bool isMC() {
     return m_isMC;
   }
   void SetIsMC(bool val) {
     m_isMC = val;
   }

protected:
   CandidatePair minimizeChiSquare(int i, int j);
   virtual double getResolution(TLorentzVector &particle) {
      return 1;
   }
   int m_dim;
   double m_mass;
   double m_width;
   std::vector<TLorentzVector> m_particle;

private:
   ROOT::Math::Minimizer *m_min;
   bool m_isMC;
};

#endif
