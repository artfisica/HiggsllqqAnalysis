// author: Valerio Ippolito <valerio.ippolito@cern.ch>
// see HqqllAnalysis.h for details

#ifndef __JetKinematicFitter_C__
#define __JetKinematicFitter_C__

#include "HiggsllqqAnalysis/JetKinematicFitter.h"

JetKinematicFitter::JetKinematicFitter(int n, double mass, double width) : KinematicFitter(n, mass, width)
{
  m_jer = new JERProvider("AntiKt4TopoJES", "Truth", "JetResolution/share/JERProviderPlots.root");
  //m_jer = new JERProvider("AntiKt4TopoJES", "Truth", "JERProviderPlots.root");
  m_jer->init();
}

JetKinematicFitter::~JetKinematicFitter()
{
   delete m_jer;
}

double JetKinematicFitter::getResolution(TLorentzVector &particle)
{
   //std::cout << "DEBUG: correcting " << particle.Pt() / 1000. << " in eta " << particle.Eta() << " with fractional " << m_jer->getSigma(particle.Pt() / 1000, particle.Eta()) << std::endl;

   if(isMC()) return particle.Pt() * m_jer->getRelResolutionMC(particle.Pt() / 1000., particle.Eta()); // pt for the tool is in GeV!!!
   return particle.Pt() * m_jer->getRelResolutionData(particle.Pt() / 1000., particle.Eta()); // pt for the tool is in GeV!!!
}
#endif
