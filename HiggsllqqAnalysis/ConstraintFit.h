/////////////////////////////////////////
//
// Z mass constraint fit
// 
// original version by:
// Konstantinos Nikolopoulos
// Univ. of Athens
//
// changes in resolution model by:
// Valerio Ippolito
// valerio.ippolito@cern.ch
////////////////////////////////////////
#ifndef CONSTRAINTFIT_H
#define CONSTRAINTFIT_H

#include "HiggsllqqAnalysis/ResolutionModel.h"
#include "CLHEP/Matrix/Matrix.h"

class Matrix;

class ConstraintFit
{
 public:
  ConstraintFit(void);                   // default constructor
  ConstraintFit(double, bool, double);   // constructor setting the constraint characteristics.

  ~ConstraintFit(void); // default destructor 
  void  MassFitInterface(double (*)[4],double [6][6],int);
  double MassFitRun(double (*)[4],double [6][6]);
  double MassFitRun(double (*)[4],double [6][6],double);
  //double  MassFitInterface(double (*)[4],double (*)[3][3],int);

  void SetConstraintMass(double mass) { m_conMass = mass;};
  void SetConstraintWidth(double width) {m_conWidth = width;};
  void SetConstraintHasWidth(bool haswidth) {m_conHasWidth = haswidth;};
  double LikelihoodMass(void);
  double LikelihoodMass2(void);
  double LikelihoodMass(double);
  double LikelihoodMass2(double);

  void SetResolutionModel(ResolutionModel *resModel);

 private:
  // data members
  bool     m_conHasWidth;
  double   m_conMass;
  double   m_conWidth;

  int      m_parameters;
  int      m_obj;
  double * m_objmass;
  HepMatrix * m_parametersInit;
  HepMatrix * m_covarianceInit;
  HepMatrix * m_parametersFinal;
  HepMatrix * m_covarianceFinal;
  HepMatrix * m_chi2;

  ResolutionModel *m_resModel; // ownership is not held

  double CalculateChi2(const HepMatrix* p0, HepMatrix* var,const double Mass);
  double MassFitCalculation(HepMatrix*, HepMatrix*,const double);
  double MassFit(const HepMatrix*, HepMatrix*,const double, HepMatrix*,HepMatrix*);
  void ConstraintCalculation(const HepMatrix*,const double,HepMatrix*,HepMatrix*);
};
#endif // CONSTRAINTFIT_H
