///////////////////////////////////////////////////////////
// A C++ implementation of Mass constraint fitting
// 23/09/2006 
// K. Nikolopoulos
// --- * --- * --- * --- * ---* --- * --- * --- * ---
//
//
//
//
#include "HiggsllqqAnalysis/ConstraintFit.h"
#include "TMath.h"
double ConstraintFit::LikelihoodMass2(double MassResol)
{

    //--------------
  double p[4][2] = {0.};
  p[0][0] = (*m_parametersInit)(1,1);
  p[1][0] = (*m_parametersInit)(2,1);
  p[2][0] = (*m_parametersInit)(3,1);
  p[3][0] = m_objmass[0]*m_objmass[0];
  for(int i=0;i<3;i++)
    p[3][0] += p[i][0]*p[i][0];
  p[3][0] = sqrt(p[3][0]);
  
  p[0][1] = (*m_parametersInit)(4,1);
  p[1][1] = (*m_parametersInit)(5,1);
  p[2][1] = (*m_parametersInit)(6,1);
  p[3][1] = m_objmass[1]*m_objmass[1];
  for(int i=0;i<3;i++)
    p[3][1] += p[i][1]*p[i][1];
  p[3][1] = sqrt(p[3][1]);

  double Etot = p[3][0] + p[3][1];
  double Xtot = p[0][0] + p[0][1];
  double Ytot = p[1][0] + p[1][1];
  double Ztot = p[2][0] + p[2][1];
  double InitMass =  Etot*Etot-Xtot*Xtot-Ytot*Ytot-Ztot*Ztot;
  InitMass = sqrt(InitMass);  
  /*
  double maxmass = InitMass;
  double max = -(maxmass-InitMass)*(maxmass-InitMass)/2./MassResol/MassResol-TMath::Log((maxmass*maxmass-m_conMass*m_conMass)*(maxmass*maxmass-m_conMass*m_conMass)+m_conMass*m_conMass*m_conWidth*m_conWidth);
  //  std::cout << "** "<<maxmass << " "<< max<<std::endl;

  for(int i=1;i<401;i++)
    {
      double ytest = InitMass + (m_conMass-InitMass)/400*i;
      double val = -(ytest-InitMass)*(ytest-InitMass)/2./MassResol/MassResol-TMath::Log((ytest*ytest-m_conMass*m_conMass)*(ytest*ytest-m_conMass*m_conMass)+m_conMass*m_conMass*m_conWidth*m_conWidth);
      //std::cout << i << " "<<ytest << " "<< val<<std::endl;
      if(val>max)
	{
	  max=val;
	  maxmass=ytest;
	}
    }
  */
  //return maxmass;
  double sig = MassResol;
  double xLeft =InitMass;
  double xRight=m_conMass;
  if(m_conMass<InitMass)
    {
      xLeft = m_conMass;
      xRight= InitMass;
    }
  //double dLinitL = (InitMass-xLeft)/sig/sig-4.*(xLeft*xLeft-m_conMass*m_conMass)*xLeft/((xLeft*xLeft-m_conMass*m_conMass)*(xLeft*xLeft-m_conMass*m_conMass)+m_conMass*m_conMass*m_conWidth*m_conWidth);
  //double dLinitR = (InitMass-xRight)/sig/sig-4.*(xRight*xLeft-m_conMass*m_conMass)*xRight/((xRight*xRight-m_conMass*m_conMass)*(xRight*xRight-m_conMass*m_conMass)+m_conMass*m_conMass*m_conWidth*m_conWidth);
  double dLinitL = (InitMass-xLeft)/sig/sig-2.*(xLeft-m_conMass)/((xLeft-m_conMass)*(xLeft-m_conMass)+m_conWidth*m_conWidth);
  double dLinitR = (InitMass-xRight)/sig/sig-2.*(xRight-m_conMass)/((xRight-m_conMass)*(xRight-m_conMass)+m_conWidth*m_conWidth);
  //std::cout << dLinitL << " "<< dLinitR<<std::endl;
  if(dLinitL*dLinitR<0.)
    {
      while(xRight-xLeft>1.)//1 MeV
	{
	  double xM = (xRight+xLeft)/2.;
	  //double dL = (InitMass-xM)/sig/sig-4.*(xM*xLeft-m_conMass*m_conMass)*xM/((xM*xM-m_conMass*m_conMass)*(xM*xM-m_conMass*m_conMass)+m_conMass*m_conMass*m_conWidth*m_conWidth);
	  double dL = (InitMass-xM)/sig/sig-2.*(xM-m_conMass)/((xM-m_conMass)*(xM-m_conMass)+m_conWidth*m_conWidth); 

	  if(dL*dLinitL<0.)
	    {
	      xRight=xM;
	      dLinitR=dL;
	    }
	  else
	    {
	      xLeft=xM;
	      dLinitL=dL;
	    }
	}
      return (xLeft+xRight)/2.; 
    }
  else
    {
      if(dLinitL>dLinitR)
	return xLeft;
      else
	return xRight;
    }  
}
double ConstraintFit::LikelihoodMass(double MassResol)
{

    //--------------
  double p[4][2] = {0.};
  p[0][0] = (*m_parametersInit)(1,1);
  p[1][0] = (*m_parametersInit)(2,1);
  p[2][0] = (*m_parametersInit)(3,1);
  p[3][0] = m_objmass[0]*m_objmass[0];
  for(int i=0;i<3;i++)
    p[3][0] += p[i][0]*p[i][0];
  p[3][0] = sqrt(p[3][0]);
  
  p[0][1] = (*m_parametersInit)(4,1);
  p[1][1] = (*m_parametersInit)(5,1);
  p[2][1] = (*m_parametersInit)(6,1);
  p[3][1] = m_objmass[1]*m_objmass[1];
  for(int i=0;i<3;i++)
    p[3][1] += p[i][1]*p[i][1];
  p[3][1] = sqrt(p[3][1]);

  double Etot = p[3][0] + p[3][1];
  double Xtot = p[0][0] + p[0][1];
  double Ytot = p[1][0] + p[1][1];
  double Ztot = p[2][0] + p[2][1];
  double InitMass =  Etot*Etot-Xtot*Xtot-Ytot*Ytot-Ztot*Ztot;
  InitMass = sqrt(InitMass);

  double sig = MassResol;
  double xLeft =InitMass;
  double xRight=m_conMass;
  if(m_conMass<InitMass)
    {
      xLeft = m_conMass;
      xRight= InitMass;
    }
  double dLinitL = (InitMass-xLeft)/sig/sig-2.*(xLeft-m_conMass)/((xLeft-m_conMass)*(xLeft-m_conMass)+m_conWidth*m_conWidth);
  double dLinitR = (InitMass-xRight)/sig/sig-2.*(xRight-m_conMass)/((xRight-m_conMass)*(xRight-m_conMass)+m_conWidth*m_conWidth);
  if(dLinitL*dLinitR<0.)
    {
      while(xRight-xLeft>1.)//1 MeV
	{
	  double xM = (xRight+xLeft)/2.;
	  double dL = (InitMass-xM)/sig/sig-2.*(xM-m_conMass)/((xM-m_conMass)*(xM-m_conMass)+m_conWidth*m_conWidth); 
	  if(dL*dLinitL<0.)
	    {
	      xRight=xM;
	      dLinitR=dL;
	    }
	  else
	    {
	      xLeft=xM;
	      dLinitL=dL;
	    }
	}
      return (xLeft+xRight)/2.; 
    }
  else
    {
      if(dLinitL>dLinitR)
	return xLeft;
      else
	return xRight;
    }  
}
double ConstraintFit::LikelihoodMass2(void)
{
  //--------------
  double p[4][2] = {0.};
  p[0][0] = (*m_parametersInit)(1,1);// (*m_parametersInit)(1,1) * cos((*m_parametersInit)(2,1)) * sin((*m_parametersInit)(3,1));
  p[1][0] = (*m_parametersInit)(2,1);// (*m_parametersInit)(1,1) * sin((*m_parametersInit)(2,1)) * sin((*m_parametersInit)(3,1));
  p[2][0] = (*m_parametersInit)(3,1);//(*m_parametersInit)(1,1) * cos((*m_parametersInit)(3,1));
  p[3][0] = m_objmass[0]*m_objmass[0];
  for(int i=0;i<3;i++)
    p[3][0] += p[i][0]*p[i][0];
  p[3][0] = sqrt(p[3][0]);
  
  p[0][1] = (*m_parametersInit)(4,1);//(*m_parametersInit)(4,1) * cos((*m_parametersInit)(5,1)) * sin((*m_parametersInit)(6,1));
  p[1][1] = (*m_parametersInit)(5,1);//(*m_parametersInit)(4,1) * sin((*m_parametersInit)(5,1)) * sin((*m_parametersInit)(6,1));
  p[2][1] = (*m_parametersInit)(6,1);//(*m_parametersInit)(4,1) * cos((*m_parametersInit)(6,1));
  p[3][1] = m_objmass[1]*m_objmass[1];
  for(int i=0;i<3;i++)
    p[3][1] += p[i][1]*p[i][1];
  p[3][1] = sqrt(p[3][1]);

  double Etot = p[3][0] + p[3][1];
  double Xtot = p[0][0] + p[0][1];
  double Ytot = p[1][0] + p[1][1];
  double Ztot = p[2][0] + p[2][1];
  double InitMass =  Etot*Etot-Xtot*Xtot-Ytot*Ytot-Ztot*Ztot;
  InitMass = sqrt(InitMass);

  //// --- 
  //what is the uncertainty on InitMass?
  HepMatrix JacobianMass(6,1,0);
  JacobianMass(1,1) = (2.*(Etot)*p[0][0]/p[3][0] - 2.* Xtot );
  JacobianMass(2,1) = (2.*(Etot)*p[1][0]/p[3][0] - 2.* Ytot );
  JacobianMass(3,1) = (2.*(Etot)*p[2][0]/p[3][0] - 2.* Ztot );
  JacobianMass(4,1) = (2.*(Etot)*p[0][1]/p[3][1] - 2.* Xtot );
  JacobianMass(5,1) = (2.*(Etot)*p[1][1]/p[3][1] - 2.* Ytot );
  JacobianMass(6,1) = (2.*(Etot)*p[2][1]/p[3][1] - 2.* Ztot );
  double sig = (JacobianMass.T()*(*m_covarianceInit)*JacobianMass)(1,1);
  sig=sqrt(sig);
  sig/=2*InitMass;
//std::cout << "MASSA: " << InitMass << " incertezza " << sig << std::endl; // kostas
  double MassResol=sig;
  double maxmass = InitMass;


  double max = m_resModel->EvaluateLogPdf(maxmass, MassResol);
  // original version below
  // double max = -(maxmass-InitMass)*(maxmass-InitMass)/2./MassResol/MassResol-TMath::Log((maxmass*maxmass-m_conMass*m_conMass)*(maxmass*maxmass-m_conMass*m_conMass)+m_conMass*m_conMass*m_conWidth*m_conWidth);
  //  std::cout << "** "<<maxmass << " "<< max<<std::endl;


  for(int i=1;i<401;i++)
    {
      double ytest = InitMass + (m_conMass-InitMass)/400*i;
      double val = m_resModel->EvaluateLogPdf(ytest, MassResol);

      // original version below
      //double val = -(ytest-InitMass)*(ytest-InitMass)/2./MassResol/MassResol-TMath::Log((ytest*ytest-m_conMass*m_conMass)*(ytest*ytest-m_conMass*m_conMass)+m_conMass*m_conMass*m_conWidth*m_conWidth);
      //std::cout << i << " "<<ytest << " "<< val<<std::endl;
      if(val>max)
	{
	  max=val;
	  maxmass=ytest;
	}
    }

//std::cout << "MASSA DOPO: " << maxmass << std::endl; // kostas

  return maxmass;

  //std::cout<< sig <<std::endl;
  //std::cout << " Sig "<< sig<< " init mass " <<InitMass<<std::endl;
  //double A = 2.*m_conMass+InitMass;
  //double B = 2.*sig*sig+m_conWidth*m_conWidth+m_conMass*m_conMass+2.*InitMass*m_conMass;
  //double C = -(m_conMass*m_conMass*InitMass+m_conWidth*m_conWidth*InitMass+2.*sig*sig*m_conMass);
  //double dx = (m_conMass-InitMass)/100.;
  double xLeft =InitMass;
  double xRight=m_conMass;
  if(m_conMass<InitMass)
    {
      xLeft = m_conMass;
      xRight= InitMass;
    }
  /*xLeft*xLeft*xLeft-A*xLeft*xLeft+B*xLeft+C;/*/
  double dLinitL = (InitMass-xLeft)/sig/sig-2.*(xLeft-m_conMass)/((xLeft-m_conMass)*(xLeft-m_conMass)+m_conWidth*m_conWidth);
  /*xLeft*xLeft*xLeft-A*xLeft*xLeft+B*xLeft+C;/*/
  double dLinitR = (InitMass-xRight)/sig/sig-2.*(xRight-m_conMass)/((xRight-m_conMass)*(xRight-m_conMass)+m_conWidth*m_conWidth);
  if(dLinitL*dLinitR<0.)
    {
      while(xRight-xLeft>1.)//1 MeV
	{
	  double xM = (xRight+xLeft)/2.;
	  /*xM*xM*xM-A*xM*xM+B*xM+C;/*/
	  double dL = (InitMass-xM)/sig/sig-2.*(xM-m_conMass)/((xM-m_conMass)*(xM-m_conMass)+m_conWidth*m_conWidth); 
	  if(dL*dLinitL<0.)
	    {
	      xRight=xM;
	      dLinitR=dL;
	    }
	  else
	    {
	      xLeft=xM;
	      dLinitL=dL;
	    }
	}
      
      //std::cout <<"mass didd "<< InitMass << " "<<(xLeft+xRight)/2.<<std::endl;

      return (xLeft+xRight)/2.; 
    }
  else
    {
      if(dLinitL>dLinitR)
	{
	  //std::cout <<"mass didd "<< InitMass << " "<<xLeft<<std::endl;
	  return xLeft;
	}
      else
	{
	  //std::cout <<"mass didd "<< InitMass << " "<<xRight<<std::endl;
	  return xRight;
	}
    }
  
}
double ConstraintFit::LikelihoodMass(void)
{
  //--------------
  double p[4][2] = {0.};
  p[0][0] = (*m_parametersInit)(1,1);// (*m_parametersInit)(1,1) * cos((*m_parametersInit)(2,1)) * sin((*m_parametersInit)(3,1));
  p[1][0] = (*m_parametersInit)(2,1);// (*m_parametersInit)(1,1) * sin((*m_parametersInit)(2,1)) * sin((*m_parametersInit)(3,1));
  p[2][0] = (*m_parametersInit)(3,1);//(*m_parametersInit)(1,1) * cos((*m_parametersInit)(3,1));
  p[3][0] = m_objmass[0]*m_objmass[0];
  for(int i=0;i<3;i++)
    p[3][0] += p[i][0]*p[i][0];
  p[3][0] = sqrt(p[3][0]);
  
  p[0][1] = (*m_parametersInit)(4,1);//(*m_parametersInit)(4,1) * cos((*m_parametersInit)(5,1)) * sin((*m_parametersInit)(6,1));
  p[1][1] = (*m_parametersInit)(5,1);//(*m_parametersInit)(4,1) * sin((*m_parametersInit)(5,1)) * sin((*m_parametersInit)(6,1));
  p[2][1] = (*m_parametersInit)(6,1);//(*m_parametersInit)(4,1) * cos((*m_parametersInit)(6,1));
  p[3][1] = m_objmass[1]*m_objmass[1];
  for(int i=0;i<3;i++)
    p[3][1] += p[i][1]*p[i][1];
  p[3][1] = sqrt(p[3][1]);

  double Etot = p[3][0] + p[3][1];
  double Xtot = p[0][0] + p[0][1];
  double Ytot = p[1][0] + p[1][1];
  double Ztot = p[2][0] + p[2][1];
  double InitMass =  Etot*Etot-Xtot*Xtot-Ytot*Ytot-Ztot*Ztot;
  InitMass = sqrt(InitMass);

  //// --- 
  //what is the uncertainty on InitMass?
  HepMatrix JacobianMass(6,1,0);
  JacobianMass(1,1) = (2.*(Etot)*p[0][0]/p[3][0] - 2.* Xtot );
  JacobianMass(2,1) = (2.*(Etot)*p[1][0]/p[3][0] - 2.* Ytot );
  JacobianMass(3,1) = (2.*(Etot)*p[2][0]/p[3][0] - 2.* Ztot );
  JacobianMass(4,1) = (2.*(Etot)*p[0][1]/p[3][1] - 2.* Xtot );
  JacobianMass(5,1) = (2.*(Etot)*p[1][1]/p[3][1] - 2.* Ytot );
  JacobianMass(6,1) = (2.*(Etot)*p[2][1]/p[3][1] - 2.* Ztot );
  double sig = (JacobianMass.T()*(*m_covarianceInit)*JacobianMass)(1,1);
  sig=sqrt(sig);
  sig/=2*InitMass;
  //std::cout<< sig <<std::endl;
  //std::cout << " Sig "<< sig<< " init mass " <<InitMass<<std::endl;
  //double A = 2.*m_conMass+InitMass;
  //double B = 2.*sig*sig+m_conWidth*m_conWidth+m_conMass*m_conMass+2.*InitMass*m_conMass;
  //double C = -(m_conMass*m_conMass*InitMass+m_conWidth*m_conWidth*InitMass+2.*sig*sig*m_conMass);
  //double dx = (m_conMass-InitMass)/100.;
  double xLeft =InitMass;
  double xRight=m_conMass;
  if(m_conMass<InitMass)
    {
      xLeft = m_conMass;
      xRight= InitMass;
    }
  /*xLeft*xLeft*xLeft-A*xLeft*xLeft+B*xLeft+C;/*/
  double dLinitL = (InitMass-xLeft)/sig/sig-2.*(xLeft-m_conMass)/((xLeft-m_conMass)*(xLeft-m_conMass)+m_conWidth*m_conWidth);
  /*xLeft*xLeft*xLeft-A*xLeft*xLeft+B*xLeft+C;/*/
  double dLinitR = (InitMass-xRight)/sig/sig-2.*(xRight-m_conMass)/((xRight-m_conMass)*(xRight-m_conMass)+m_conWidth*m_conWidth);
  if(dLinitL*dLinitR<0.)
    {
      while(xRight-xLeft>1.)//1 MeV
	{
	  double xM = (xRight+xLeft)/2.;
	  /*xM*xM*xM-A*xM*xM+B*xM+C;/*/
	  double dL = (InitMass-xM)/sig/sig-2.*(xM-m_conMass)/((xM-m_conMass)*(xM-m_conMass)+m_conWidth*m_conWidth); 
	  if(dL*dLinitL<0.)
	    {
	      xRight=xM;
	      dLinitR=dL;
	    }
	  else
	    {
	      xLeft=xM;
	      dLinitL=dL;
	    }
	}
      
      //std::cout <<"mass didd "<< InitMass << " "<<(xLeft+xRight)/2.<<std::endl;

      return (xLeft+xRight)/2.; 
    }
  else
    {
      if(dLinitL>dLinitR)
	{
	  //std::cout <<"mass didd "<< InitMass << " "<<xLeft<<std::endl;
	  return xLeft;
	}
      else
	{
	  //std::cout <<"mass didd "<< InitMass << " "<<xRight<<std::endl;
	  return xRight;
	}
    }
  
}

void ConstraintFit::SetResolutionModel(ResolutionModel *resModel) {
  // this method must be called BEFORE running likelihood minimization
  m_resModel = resModel;
}


ConstraintFit::ConstraintFit(void):
  m_conMass(0.),
  m_conWidth(0.),
  m_conHasWidth(false),
  m_parameters(3),
  m_parametersInit(NULL),
  m_covarianceInit(NULL),
  m_parametersFinal(NULL),
  m_covarianceFinal(NULL),
  m_chi2(NULL)
{}
ConstraintFit::ConstraintFit(double mass, bool haswidth, double width):
  m_conMass(mass),
  m_conWidth(width),
  m_conHasWidth(haswidth),
  m_parameters(3),
  m_parametersInit(NULL),
  m_covarianceInit(NULL),
  m_parametersFinal(NULL),
  m_covarianceFinal(NULL),
  m_chi2(NULL)
{}

ConstraintFit::~ConstraintFit(void)
{
  delete m_parametersInit;
  delete m_covarianceInit;
  delete m_parametersFinal;
  delete m_covarianceFinal;
  delete m_objmass;
  delete m_chi2;
}
//double ConstraintFit::MassFitInterface(double (*pin)[4],double (*sigmain)[3][3], int iobj)
void ConstraintFit::MassFitInterface(double (*pin)[4],double sigmain[6][6], int iobj)
{
  // the m_parameters first parameters of pin are the fitted parameters
  // the next one is the mass of the particle
  // by definition 
  // parameter[0] = P /1.e3;
  // parameter[1] = phi;
  // parameter[2] = theta;

  m_obj            = iobj;
  int dimension    = m_parameters*m_obj;
  m_parametersInit = new HepMatrix(dimension,1,0);
  m_covarianceInit = new HepMatrix(dimension,dimension,0);
  m_parametersFinal= new HepMatrix(dimension,1,0);
  m_covarianceFinal= new HepMatrix(dimension,dimension,0);
  m_chi2           = new HepMatrix(1,1,0);
  
  m_objmass        = new double [2];
  for(int i=0;i<m_obj;i++)
    m_objmass[i]  = pin[i][3];
  // 1 p
  // 2 phi
  // 3 theta
  /*
  HepMatrix Jacobian(2*m_parameters,2*m_parameters,0);
  for(int i=0;i<2;i++)
    {
      double p = sqrt(pin[i][0]*pin[i][0]+pin[i][1]*pin[i][1]+pin[i][2]*pin[i][2]);
      double phi = atan2(pin[i][1],pin[i][0]);
      double theta = acos(pin[i][2]/p);
      (*m_parametersInit)(1+3*i,1) = p;
      (*m_parametersInit)(2+3*i,1) = phi;
      (*m_parametersInit)(3+3*i,1) = theta;

      Jacobian(1+3*i,1+3*i)= pin[i][0]/p;
      Jacobian(1+3*i,2+3*i)= 1./(p*cos(theta)*cos(phi));
      Jacobian(1+3*i,3+3*i)= 1./(-p*sin(theta)*sin(phi));
      Jacobian(2+3*i,1+3*i)= pin[i][1]/p;
      Jacobian(2+3*i,2+3*i)= 1./(p*cos(theta)*sin(phi));
      Jacobian(2+3*i,3+3*i)= 1./(p*sin(theta)*cos(phi));
      Jacobian(3+3*i,1+3*i)= pin[i][2]/p;
      Jacobian(3+3*i,2+3*i)= 1./(-p*sin(phi));
      Jacobian(3+3*i,3+3*i)= 0.;
 
    }
  */
  for(int i=0;i<m_parameters;i++)
    {
      (*m_parametersInit)(i+0+1,1) = pin[0][i];
      (*m_parametersInit)(i+3+1,1) = pin[1][i];
    }
  for(int i=0;i<2*m_parameters;i++)
    for(int j=0;j<2*m_parameters;j++)
      (*m_covarianceInit)(i+0+1,j+0+1)= sigmain[i][j];
  //transform the covariance matrices for px,py,pz to p,phi,theta
  //going from d0,z0,phi,theta,P --> d0,z0,px,py,pz
  //(*m_covarianceInit)=Jacobian.T()*(*m_covarianceInit)*Jacobian;
  
}

double ConstraintFit::MassFitRun(double (*pin)[4],double sigmain[6][6])
{
  if(!m_conHasWidth)
    {
      *m_parametersFinal = *m_parametersInit;
      *m_covarianceFinal = *m_covarianceInit;
      double chi2 = MassFit(m_parametersInit,m_covarianceInit,m_conMass,m_parametersFinal,m_covarianceFinal);
      double chi2prob = TMath::Prob(chi2,1);
      double chi22 = CalculateChi2(m_parametersInit,m_covarianceInit,m_conMass);
      for(int i=0;i<m_parameters;i++)
	{
	  pin[0][i] = (*m_parametersFinal)(i+0+1,1);
	  pin[1][i] = (*m_parametersFinal)(i+3+1,1);
	  //(*m_covarianceInit)(i+0+1,i+0+1)= sigmain[0][i][i]*sigmain[0][i][i];
	  //(*m_covarianceInit)(i+3+1,i+3+1)= sigmain[1][i][i]*sigmain[1][i][i];
	}
      return chi2;
    }
  HepMatrix *m_parametersFit=new HepMatrix(*m_parametersInit);
  *m_parametersFit   = *m_parametersInit;
  //std::cout << "--> "<< LikelihoodMass2()<< " "<< LikelihoodMass()<<std::endl;
  double Mass = LikelihoodMass2();
  double chi2 = MassFit(m_parametersInit,m_covarianceInit,Mass,m_parametersFit,m_covarianceFinal);
  *m_parametersFinal = *m_parametersFit;
  *m_covarianceFinal = *m_covarianceInit;
  double chi2prob = TMath::Prob(chi2,1);
  /*
  HepMatrix *m_parametersFit=new HepMatrix(*m_parametersInit);
  int Ntrials = 100;
  double chi2totprob_old = 1.e7;
  for(int i=0;i<Ntrials;i++)
    {
      
      *m_parametersFinal = *m_parametersFit;
      *m_covarianceFinal = *m_covarianceInit;
      *m_parametersFit   = *m_parametersInit;
      double Mass = (m_conMass-InitMass)/(Ntrials-1)*i+InitMass;
      double chi2 = MassFit(m_parametersInit,m_covarianceInit,Mass,m_parametersFit,m_covarianceFinal);
      double bwprob = TMath::BreitWigner(Mass,m_conMass,m_conWidth)/TMath::BreitWigner(m_conMass,m_conMass,m_conWidth);

      const double ndof = 1.;
      double chi2prob = TMath::Prob(chi2,1);

      double chi2totprob =  -TMath::Log(bwprob)-TMath::Log(chi2prob);

      std::cout << "Fitting mass = " << Mass
		<< " Breit-Wigner Prob = "<< bwprob
		<< " Chi2 Prob = " << chi2prob << " Multiplication = " << chi2totprob 
		<< " Chi2 old = " << chi2totprob_old<<std::endl;

 
      if(chi2totprob < chi2totprob_old)
	chi2totprob_old = chi2totprob;
      else
	{
	  std::cout << "found minimum " <<std::endl;
	  //std::cout << (*m_parametersInit) << (*m_parametersFinal) << chi2<<std::endl;
	
	  break;
	}
    }
  */  
  for(int i=0;i<m_parameters;i++)
    {
      pin[0][i] = (*m_parametersFinal)(i+0+1,1);
      pin[1][i] = (*m_parametersFinal)(i+3+1,1);
      //(*m_covarianceInit)(i+0+1,i+0+1)= sigmain[0][i][i]*sigmain[0][i][i];
      //(*m_covarianceInit)(i+3+1,i+3+1)= sigmain[1][i][i]*sigmain[1][i][i];
    }
 
  return chi2;
  // original version below
  //return 0;

}
double ConstraintFit::MassFitRun(double (*pin)[4],double sigmain[6][6],double zresol)
{
  if(!m_conHasWidth)
    {
      *m_parametersFinal = *m_parametersInit;
      *m_covarianceFinal = *m_covarianceInit;
      double chi2 = MassFit(m_parametersInit,m_covarianceInit,m_conMass,m_parametersFinal,m_covarianceFinal);
      double chi2prob = TMath::Prob(chi2,1);
      double chi22 = CalculateChi2(m_parametersInit,m_covarianceInit,m_conMass);
      for(int i=0;i<m_parameters;i++)
	{
	  pin[0][i] = (*m_parametersFinal)(i+0+1,1);
	  pin[1][i] = (*m_parametersFinal)(i+3+1,1);
	  //(*m_covarianceInit)(i+0+1,i+0+1)= sigmain[0][i][i]*sigmain[0][i][i];
	  //(*m_covarianceInit)(i+3+1,i+3+1)= sigmain[1][i][i]*sigmain[1][i][i];
	}
      return chi2;
    }
  HepMatrix *m_parametersFit=new HepMatrix(*m_parametersInit);
  *m_parametersFit   = *m_parametersInit;
  double Mass = LikelihoodMass(zresol);
  double chi2 = MassFit(m_parametersInit,m_covarianceInit,Mass,m_parametersFit,m_covarianceFinal);
  *m_parametersFinal = *m_parametersFit;
  *m_covarianceFinal = *m_covarianceInit;
  double chi2prob = TMath::Prob(chi2,1);
  /*
  HepMatrix *m_parametersFit=new HepMatrix(*m_parametersInit);
  int Ntrials = 100;
  double chi2totprob_old = 1.e7;
  for(int i=0;i<Ntrials;i++)
    {
      
      *m_parametersFinal = *m_parametersFit;
      *m_covarianceFinal = *m_covarianceInit;
      *m_parametersFit   = *m_parametersInit;
      double Mass = (m_conMass-InitMass)/(Ntrials-1)*i+InitMass;
      double chi2 = MassFit(m_parametersInit,m_covarianceInit,Mass,m_parametersFit,m_covarianceFinal);
      double bwprob = TMath::BreitWigner(Mass,m_conMass,m_conWidth)/TMath::BreitWigner(m_conMass,m_conMass,m_conWidth);

      const double ndof = 1.;
      double chi2prob = TMath::Prob(chi2,1);

      double chi2totprob =  -TMath::Log(bwprob)-TMath::Log(chi2prob);

      std::cout << "Fitting mass = " << Mass
		<< " Breit-Wigner Prob = "<< bwprob
		<< " Chi2 Prob = " << chi2prob << " Multiplication = " << chi2totprob 
		<< " Chi2 old = " << chi2totprob_old<<std::endl;

 
      if(chi2totprob < chi2totprob_old)
	chi2totprob_old = chi2totprob;
      else
	{
	  std::cout << "found minimum " <<std::endl;
	  //std::cout << (*m_parametersInit) << (*m_parametersFinal) << chi2<<std::endl;
	
	  break;
	}
    }
  */  
  for(int i=0;i<m_parameters;i++)
    {
      pin[0][i] = (*m_parametersFinal)(i+0+1,1);
      pin[1][i] = (*m_parametersFinal)(i+3+1,1);
      //(*m_covarianceInit)(i+0+1,i+0+1)= sigmain[0][i][i]*sigmain[0][i][i];
      //(*m_covarianceInit)(i+3+1,i+3+1)= sigmain[1][i][i]*sigmain[1][i][i];
    }
 
  return 0;

}
void ConstraintFit::ConstraintCalculation(const HepMatrix*p0,const double Mass,HepMatrix *D,HepMatrix *d)
{
  // the constraint is (E1+E2)^2-(P1x+P2x)^2-(P1y+P2y)^2-(P1z-P2z)^2-m^2 = 0.
  // the partial derivatives of the constraint with respect to px,py,pz are
  // i->particle
  // j->x,y,z
  // dF/dPij = 2 (E1+E2)*Pij/Ei - 2(P1j+P2j)
  //
  // partial derivatives of px,py,pz with respect to the fit parameters
  // dpx/dp = cosPhi *sinTheta
  // dpy/dp = sinPhi *sinTheta
  // dpz/dp = cosTheta
  // --
  // dpx/dtheta = px/tanTheta
  // dpy/dtheta = py/tanTheta
  // dpz/dtheta =-pz*tanTheta
  // --
  // dpx/dphi = -px*tanTheta
  // dpy/dphi =  px/tanTheta
  // dpz/dphi =  0.
  // 
  //calculating the constraint linear expansion at the current point
  //
  double p[4][2] = {0.};
  p[0][0] = (*p0)(1,1) * cos((*p0)(2,1)) * sin((*p0)(3,1));
  p[1][0] = (*p0)(1,1) * sin((*p0)(2,1)) * sin((*p0)(3,1));
  p[2][0] = (*p0)(1,1) * cos((*p0)(3,1));
  for(int i=0;i<3;i++)
    p[3][0] += p[i][0]*p[i][0];
  p[3][0] = sqrt(p[3][0]);
  
  p[0][1] = (*p0)(4,1) * cos((*p0)(5,1)) * sin((*p0)(6,1));
  p[1][1] = (*p0)(4,1) * sin((*p0)(5,1)) * sin((*p0)(6,1));
  p[2][1] = (*p0)(4,1) * cos((*p0)(6,1));
  for(int i=0;i<3;i++)
    p[3][1] += p[i][1]*p[i][1];
  p[3][1] = sqrt(p[3][1]);

  double Etot = p[3][0] + p[3][1];

  double constraintD1 = (2.*(Etot)*p[0][0]/p[3][0] - 2.*(p[0][0]+p[0][1]) ) * (cos((*p0)(2,1)) *sin((*p0)(3,1)))
                      + (2.*(Etot)*p[1][0]/p[3][0] - 2.*(p[1][0]+p[1][1]) ) * (sin((*p0)(2,1)) *sin((*p0)(3,1)))
                      + (2.*(Etot)*p[2][0]/p[3][0] - 2.*(p[2][0]+p[2][1]) ) * (                 cos((*p0)(3,1)));
  double constraintD2 = (2.*(Etot)*p[0][0]/p[3][0] - 2.*(p[0][0]+p[0][1]) ) * (-p[0][0] * tan((*p0)(2,1)))
                      + (2.*(Etot)*p[1][0]/p[3][0] - 2.*(p[1][0]+p[1][1]) ) * ( p[1][0] / tan((*p0)(2,1)))
                      + (2.*(Etot)*p[2][0]/p[3][0] - 2.*(p[2][0]+p[2][1]) ) *  0.;
  double constraintD3 = (2.*(Etot)*p[0][0]/p[3][0] - 2.*(p[0][0]+p[0][1]) ) * ( p[0][0] / tan((*p0)(3,1)))
                      + (2.*(Etot)*p[1][0]/p[3][0] - 2.*(p[1][0]+p[1][1]) ) * ( p[1][0] / tan((*p0)(3,1)))
                      + (2.*(Etot)*p[2][0]/p[3][0] - 2.*(p[2][0]+p[2][1]) ) * (-p[2][0] * tan((*p0)(3,1)));

  double constraintD4 = (2.*(Etot)*p[0][1]/p[3][1] - 2.*(p[0][0]+p[0][1]) ) * (cos((*p0)(5,1)) *sin((*p0)(6,1)))
                      + (2.*(Etot)*p[1][1]/p[3][1] - 2.*(p[1][0]+p[1][1]) ) * (sin((*p0)(5,1)) *sin((*p0)(6,1)))
                      + (2.*(Etot)*p[2][1]/p[3][1] - 2.*(p[2][0]+p[2][1]) ) * (                 cos((*p0)(6,1)));
  double constraintD5 = (2.*(Etot)*p[0][1]/p[3][1] - 2.*(p[0][0]+p[0][1]) ) * (-p[0][1] * tan((*p0)(5,1)))
                      + (2.*(Etot)*p[1][1]/p[3][1] - 2.*(p[1][0]+p[1][1]) ) * ( p[1][1] / tan((*p0)(5,1)))
                      + (2.*(Etot)*p[2][1]/p[3][1] - 2.*(p[2][0]+p[2][1]) ) *  0.;
  double constraintD6 = (2.*(Etot)*p[0][1]/p[3][1] - 2.*(p[0][0]+p[0][1]) ) * ( p[0][1] / tan((*p0)(6,1)))
                      + (2.*(Etot)*p[1][1]/p[3][1] - 2.*(p[1][0]+p[1][1]) ) * ( p[1][1] / tan((*p0)(6,1)))
                      + (2.*(Etot)*p[2][1]/p[3][1] - 2.*(p[2][0]+p[2][1]) ) * (-p[2][1] * tan((*p0)(6,1)));

  
  double constraintd =  Etot*Etot;
  for(int i=0;i<3;i++)
    {
      constraintd = constraintd - (p[i][0]+p[i][1])*(p[i][0]+p[i][1]);
    }
  constraintd = constraintd - Mass*Mass;
  //std::cout << "constraint = " <<constraintd << " mass value " <<Mass<<std::endl;
  D=new HepMatrix(1,6,0);
  (*D)(1,1) = constraintD1;
  (*D)(1,2) = constraintD2;
  (*D)(1,3) = constraintD3;
  (*D)(1,4) = constraintD4;
  (*D)(1,5) = constraintD5;
  (*D)(1,6) = constraintD6;
  d = new HepMatrix(1,1,0);
  (*d)(1,1)=constraintd;

}
double ConstraintFit::CalculateChi2(const HepMatrix* p0, HepMatrix* var,const double Mass)
{
  int ierr;
  double p[4][2] = {0.};
  p[0][0] = (*p0)(1,1);
  p[1][0] = (*p0)(2,1);
  p[2][0] = (*p0)(3,1);
  p[3][0] = m_objmass[0]*m_objmass[0];
  for(int i=0;i<3;i++)
    p[3][0] += p[i][0]*p[i][0];
  p[3][0] = sqrt(p[3][0]);
  
  p[0][1] = (*p0)(4,1);
  p[1][1] = (*p0)(5,1);
  p[2][1] = (*p0)(6,1);
  p[3][1] = m_objmass[1]*m_objmass[1];
  for(int i=0;i<3;i++)
    p[3][1] += p[i][1]*p[i][1];
  p[3][1] = sqrt(p[3][1]);

  double Etot = p[3][0] + p[3][1];
  double Xtot = p[0][0] + p[0][1];
  double Ytot = p[1][0] + p[1][1];
  double Ztot = p[2][0] + p[2][1];
  
  double constraintD1 = (2.*(Etot)*p[0][0]/p[3][0] - 2.* Xtot );
  double constraintD2 = (2.*(Etot)*p[1][0]/p[3][0] - 2.* Ytot );
  double constraintD3 = (2.*(Etot)*p[2][0]/p[3][0] - 2.* Ztot );
  double constraintD4 = (2.*(Etot)*p[0][1]/p[3][1] - 2.* Xtot );
  double constraintD5 = (2.*(Etot)*p[1][1]/p[3][1] - 2.* Ytot );
  double constraintD6 = (2.*(Etot)*p[2][1]/p[3][1] - 2.* Ztot );
  
  double constraintd =  Etot*Etot - Xtot*Xtot - Ytot*Ytot - Ztot*Ztot;
  constraintd = constraintd - Mass*Mass;
  HepMatrix D(1,6,0);
  D(1,1) = constraintD1;
  D(1,2) = constraintD2;
  D(1,3) = constraintD3;
  D(1,4) = constraintD4;
  D(1,5) = constraintD5;
  D(1,6) = constraintD6;
  HepMatrix d(1,1,0);
  d(1,1)=constraintd;
  
  HepMatrix DVD(D*(*var)*D.T());
  DVD.invert(ierr);
  //if(ierr!=0)
  //std::cout << "matrix inversion failed " <<std::endl;
  HepMatrix VD(DVD);
  HepMatrix lambda(VD*d); // afou to anaptygma ginetai panw sto p0 (ie shmeio anaptygmatos kai arxikh timh sympiptoun!)
  DVD.invert(ierr);
  double chi2 = (lambda.T()*DVD*lambda)(1,1);
  //std::cout <<"chi2 "<<chi2<<std::endl;
  return chi2;
}
double ConstraintFit::MassFitCalculation(HepMatrix* p0, HepMatrix* var,const double Mass)
{
  int ierr;
  double p[4][2] = {0.};
  p[0][0] = (*p0)(1,1);//(*p0)(1,1) * cos((*p0)(2,1)) * sin((*p0)(3,1));
  p[1][0] = (*p0)(2,1);//(*p0)(1,1) * sin((*p0)(2,1)) * sin((*p0)(3,1));
  p[2][0] = (*p0)(3,1);//(*p0)(1,1) * cos((*p0)(3,1));
  p[3][0] = m_objmass[0]*m_objmass[0];
  for(int i=0;i<3;i++)
    p[3][0] += p[i][0]*p[i][0];
  p[3][0] = sqrt(p[3][0]);
  
  p[0][1] = (*p0)(4,1);//(*p0)(4,1) * cos((*p0)(5,1)) * sin((*p0)(6,1));
  p[1][1] = (*p0)(5,1);//(*p0)(4,1) * sin((*p0)(5,1)) * sin((*p0)(6,1));
  p[2][1] = (*p0)(6,1);//(*p0)(4,1) * cos((*p0)(6,1));
  p[3][1] = m_objmass[1]*m_objmass[1];
  for(int i=0;i<3;i++)
    p[3][1] += p[i][1]*p[i][1];
  p[3][1] = sqrt(p[3][1]);

  double Etot = p[3][0] + p[3][1];
  double Xtot = p[0][0] + p[0][1];
  double Ytot = p[1][0] + p[1][1];
  double Ztot = p[2][0] + p[2][1];
  
  double constraintD1 = (2.*(Etot)*p[0][0]/p[3][0] - 2.* Xtot );
  double constraintD2 = (2.*(Etot)*p[1][0]/p[3][0] - 2.* Ytot );
  double constraintD3 = (2.*(Etot)*p[2][0]/p[3][0] - 2.* Ztot );
  double constraintD4 = (2.*(Etot)*p[0][1]/p[3][1] - 2.* Xtot );
  double constraintD5 = (2.*(Etot)*p[1][1]/p[3][1] - 2.* Ytot );
  double constraintD6 = (2.*(Etot)*p[2][1]/p[3][1] - 2.* Ztot );
  
  double constraintd =  Etot*Etot - Xtot*Xtot - Ytot*Ytot - Ztot*Ztot;
  /*
  double constraintD1 = (2.*(Etot)*p[0][0]/p[3][0] - 2.* Xtot ) * (cos((*p0)(2,1)) *sin((*p0)(3,1)))
                      + (2.*(Etot)*p[1][0]/p[3][0] - 2.* Ytot ) * (sin((*p0)(2,1)) *sin((*p0)(3,1)))
                      + (2.*(Etot)*p[2][0]/p[3][0] - 2.* Ztot ) * (                 cos((*p0)(3,1)));
  double constraintD2 = (2.*(Etot)*p[0][0]/p[3][0] - 2.* Xtot ) * (-p[0][0] * tan((*p0)(2,1)))
                      + (2.*(Etot)*p[1][0]/p[3][0] - 2.* Ytot ) * ( p[1][0] / tan((*p0)(2,1)))
                      + (2.*(Etot)*p[2][0]/p[3][0] - 2.* Ztot ) *  0.;
  double constraintD3 = (2.*(Etot)*p[0][0]/p[3][0] - 2.* Xtot ) * ( p[0][0] / tan((*p0)(3,1)))
                      + (2.*(Etot)*p[1][0]/p[3][0] - 2.* Ytot ) * ( p[1][0] / tan((*p0)(3,1)))
                      + (2.*(Etot)*p[2][0]/p[3][0] - 2.* Ztot ) * (-p[2][0] * tan((*p0)(3,1)));

  double constraintD4 = (2.*(Etot)*p[0][1]/p[3][1] - 2.* Xtot ) * (cos((*p0)(5,1)) *sin((*p0)(6,1)))
                      + (2.*(Etot)*p[1][1]/p[3][1] - 2.* Ytot ) * (sin((*p0)(5,1)) *sin((*p0)(6,1)))
                      + (2.*(Etot)*p[2][1]/p[3][1] - 2.* Ztot ) * (                 cos((*p0)(6,1)));
  double constraintD5 = (2.*(Etot)*p[0][1]/p[3][1] - 2.* Xtot ) * (-p[0][1] * tan((*p0)(5,1)))
                      + (2.*(Etot)*p[1][1]/p[3][1] - 2.* Ytot ) * ( p[1][1] / tan((*p0)(5,1)))
                      + (2.*(Etot)*p[2][1]/p[3][1] - 2.* Ztot ) *  0.;
  double constraintD6 = (2.*(Etot)*p[0][1]/p[3][1] - 2.* Xtot ) * ( p[0][1] / tan((*p0)(6,1)))
                      + (2.*(Etot)*p[1][1]/p[3][1] - 2.* Ytot ) * ( p[1][1] / tan((*p0)(6,1)))
                      + (2.*(Etot)*p[2][1]/p[3][1] - 2.* Ztot ) * (-p[2][1] * tan((*p0)(6,1)));

  
  double constraintd =  Etot*Etot - Xtot*Xtot - Ytot*Ytot - Ztot*Ztot;
  */
  //std::cout << "mass " << sqrt(constraintd)<<std::endl;
  //for(int i=0;i<3;i++)
  //{
  //    constraintd = constraintd - (p[i][0]+p[i][1])*(p[i][0]+p[i][1]);
  //}
  constraintd = constraintd - Mass*Mass;
  //std::cout << "constraint = " <<constraintd << " mass value " <<Mass<<std::endl;
  HepMatrix D(1,6,0);
  D(1,1) = constraintD1;
  D(1,2) = constraintD2;
  D(1,3) = constraintD3;
  D(1,4) = constraintD4;
  D(1,5) = constraintD5;
  D(1,6) = constraintD6;
  HepMatrix d(1,1,0);
  d(1,1)=constraintd;
  
  HepMatrix DVD(D*(*var)*D.T());
  DVD.invert(ierr);
  //if(ierr!=0)
  //std::cout << "matrix inversion failed " <<std::endl;
  HepMatrix VD(DVD);
  //HepMatrix lambda(VD*(D*(*p0)+d)); // afou to anaptygma ginetai panw sto p0
  HepMatrix lambda(VD*d);
  HepMatrix test((*p0)-(*var)*D.T()*lambda);

  (*p0)(1,1)=test(1,1);
  (*p0)(2,1)=test(2,1);
  (*p0)(3,1)=test(3,1);
  (*p0)(4,1)=test(4,1);
  (*p0)(5,1)=test(5,1);
  (*p0)(6,1)=test(6,1);
  
  double constraintValue = d(1,1);
  return constraintValue;
}
double ConstraintFit::MassFit(const HepMatrix* p0, HepMatrix* var,const double Mass, HepMatrix* pOut, HepMatrix* varout)
{
  HepMatrix ivar(*var);
  int ierr;
  ivar.invert(ierr);
  //if(ierr!=0)
  //std::cout << "Warning matrix inversion failed "<<std::endl;

  bool doIter = true;
  int  maxIterations = 20;
  int  iIter  = 0;
  double constraintValue = -1.e10;
  //std::cout<<(*pOut).T();
  while(doIter)
    {
      HepMatrix p_old(*pOut);
      constraintValue = MassFitCalculation(pOut, var,Mass); //call the fit
      
      // Convergence criteria
      // 1. the parameters should not change too much (<1.e-6 relative change)
      // 2. the constraint should be satisfied very well (<1.e-6 absolute)
      double maxDiff = 5.e15;      
      HepMatrix diff = ((*pOut)-p_old);
      for(int i=0;i<m_parameters*m_obj;i++)
	{
	  diff(i+1,1) = diff(i+1,1)/p_old(i+1,1);
	  if(maxDiff > diff(i,1))
	    maxDiff = diff(i,1);
	}
      if((maxDiff<1.e-6 && TMath::Abs(constraintValue)<1.e-6) || maxIterations<=iIter)
	doIter = false;
      //std::cout << "iteration " << iIter << " " <<constraintValue <<std::endl;
      iIter++;

    }
  //std::cout << "Iterations " <<iIter <<std::endl;
  
  double chi2 = CalculateChi2(p0,var,Mass);
  return chi2;
}

