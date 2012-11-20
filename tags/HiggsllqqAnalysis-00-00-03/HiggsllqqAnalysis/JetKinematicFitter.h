////////////////////////////////////////////////////////////
// JetKinematicFitter
//
// Implementation of KinematicFitter for dijet resonances
//
// author: Valerio Ippolito <valerio.ippolito@cern.ch>
//
////////////////////////////////////////////////////////////


#ifndef __JetKinematicFitter_h__
#define __JetKinematicFitter_h__

#include "HiggsllqqAnalysis/KinematicFitter.h"

class JetKinematicFitter : public KinematicFitter {
public:
   JetKinematicFitter(int n, double mass, double width);
   ~JetKinematicFitter();

protected:
   virtual double getResolution(TLorentzVector &particle);

   JERProvider *m_jer;
};

#endif
