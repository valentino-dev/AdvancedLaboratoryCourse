#include "math.h"
#include "TMath.h"

double ElecCalib(double e_raw, double pt, double eta, 
		 double phi, double etiso, double eoverp, double mindrjet)
{
  double dummy=pt*eta*phi*etiso*eoverp*mindrjet;
  double energy = e_raw;


  if (fabs(eta)>0.0 && fabs(eta)<0.5) energy = energy * 91.2/91.5;
  else if (fabs(eta)>0.5 && fabs(eta)<1.0) energy = energy * 91.2/91.2;
  else if (fabs(eta)>1.0 && fabs(eta)<1.5) energy = energy * 91.2/90.9;
  else if (fabs(eta)>1.5 && fabs(eta)<2.0) energy = energy * 91.2/92.5;
  else if (fabs(eta)>2.0 && fabs(eta)<2.5) energy = energy * 91.2/89.9;

  if (phi>-0.5*3.2 && phi<-0.25*3.2) energy = energy * 91.2/89.9;
  else if (phi>-0.25*3.2 && phi<-0.*3.2) energy = energy * 91.2/89.8;
  else if (phi>0.0*3.2 && phi<0.25*3.2) energy = energy * 91.2/89.8;
  else if (phi>0.25*3.2 && phi<0.5*3.2) energy = energy * 91.2/89.9;

  else if (phi>0.5*3.2 && phi<0.75*3.2) energy = energy * 91.2/89.9;
  else if (phi>0.75*3.2 && phi<1*3.2) energy = energy * 91.2/89.8;
  else if (phi>-1*3.2 && phi<-0.75*3.2) energy = energy * 91.2/89.9;
  else if (phi>-0.75*3.2 && phi<-0.5*3.2) energy = energy * 91.2/89.8;
  
  if (pt>0 && pt<25) energy = energy * 91.2/89.8;
  else if (pt>25 && pt<50) energy = energy * 91.2/91.1;
  else if (pt>50 && pt<75) energy = energy * 91.2/92.1;
  else if (pt>75 && pt<100) energy = energy * 91.2/91.8;
  else if (pt>100) energy = energy * 91.2/92.1;



  return energy;
} 
