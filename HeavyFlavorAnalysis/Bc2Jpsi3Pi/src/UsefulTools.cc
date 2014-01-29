#include "../interface/UsefulTools.h"

#include "TrackingTools/IPTools/interface/IPTools.h"

#include <stdio.h>
#include <iostream>
#include <TMath.h>

UsefulTools::UsefulTools ()
{
}




//Compute cosine
double UsefulTools::computeCosine (double Lx,
								   double Ly,
								   double Lz,
								   double Px,
								   double Py,
								   double Pz)
{
  double cosine = 0. ;
  double Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  double Pnorm = sqrt(Px*Px + Py*Py + Pz*Pz);
  double LdotP = Lx*Px + Ly*Py + Lz*Pz;
  if ((Lnorm > 0.) && (Pnorm > 0.)) cosine = LdotP / (Lnorm * Pnorm);
  return cosine;
}



double UsefulTools::computeCosineError (double Lx,
										double Ly,
										double Lz,
										double Px,
										double Py,
										double Pz,
										double LxErr2,
										double LyErr2,
										double LzErr2,
										double LxyCov,
										double LxzCov,
										double LyzCov,
										double PxErr2,
										double PyErr2,
										double PzErr2,
										double PxyCov,
										double PxzCov,
										double PyzCov)
{
  double cosine = 0.;
  double cosineError = 0. ;
  double Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  double Pnorm = sqrt(Px*Px + Py*Py + Pz*Pz);
  double LdotW = computeCosine(Lx, Ly, Lz, Px, Py, Pz);
  
  if ((Lnorm > 0.) && (Pnorm > 0.))
  {
    cosine = LdotW / (Lnorm * Pnorm);
    cosineError  = sqrt( ((Lx*Pnorm - LdotW*Px) * (Lx*Pnorm - LdotW*Px) * PxErr2 +
						  (Ly*Pnorm - LdotW*Py) * (Ly*Pnorm - LdotW*Py) * PyErr2 +
						  (Lz*Pnorm - LdotW*Pz) * (Lz*Pnorm - LdotW*Pz) * PzErr2 +
						 
						  (Lx*Pnorm - LdotW*Px) * (Ly*Pnorm - LdotW*Py) * 2.*PxyCov +
						  (Lx*Pnorm - LdotW*Px) * (Lz*Pnorm - LdotW*Pz) * 2.*PxzCov +
						  (Ly*Pnorm - LdotW*Py) * (Lz*Pnorm - LdotW*Pz) * 2.*PyzCov) /
						 (Pnorm*Pnorm*Pnorm*Pnorm) +
						 
						 ((Px*Lnorm - LdotW*Lx) * (Px*Lnorm - LdotW*Lx) * LxErr2 +
						  (Py*Lnorm - LdotW*Ly) * (Py*Lnorm - LdotW*Ly) * LyErr2 +
						  (Pz*Lnorm - LdotW*Lz) * (Pz*Lnorm - LdotW*Lz) * LzErr2 +
						  
						  (Px*Lnorm - LdotW*Lx) * (Py*Lnorm - LdotW*Ly) * 2.*LxyCov +
						  (Px*Lnorm - LdotW*Lx) * (Pz*Lnorm - LdotW*Lz) * 2.*LxzCov +
						  (Py*Lnorm - LdotW*Ly) * (Pz*Lnorm - LdotW*Lz) * 2.*LyzCov) /
						 (Lnorm*Lnorm*Lnorm*Lnorm) ) / (Pnorm*Lnorm);
  }
  else
  {
    cosineError = 0.;
  }
  
  return cosineError;
}


//Compute impact parameter 3D wrt Transient Vtx
std::pair<double,double> UsefulTools::pionImpactParameter(reco::TransientTrack piTT, TransientVertex jpsiVtx)
{
    std::pair<double,double> measure;
	std::pair<bool,Measurement1D>  piIP_pair = IPTools::absoluteImpactParameter3D(piTT, jpsiVtx);
	if (piIP_pair.first)
	{
	  measure.first  = piIP_pair.second.value();
	  measure.second = piIP_pair.second.significance();
	}
    else 
    {
	  measure.first  = 0;
	  measure.second = 0;
    } 
	return measure;
}

//Compute impact parameter 3D wrt Reco Vtx
std::pair<double,double> UsefulTools::pionImpactParameter(reco::TransientTrack piTT, reco::Vertex myVtx)
{
    std::pair<double,double> measure;
	std::pair<bool,Measurement1D>  piIP_pair = IPTools::absoluteImpactParameter3D(piTT, myVtx);
	if (piIP_pair.first)
	{
	  measure.first  = piIP_pair.second.value();
	  measure.second = piIP_pair.second.significance();
	}
    else 
    {
	  measure.first  = 0;
	  measure.second = 0;
    } 
	return measure;
}

//Compute impact parameter 2D wrt beamSpot
std::pair<double,double> UsefulTools::pionIPBeamSpot(reco::TransientTrack piTT, GlobalPoint BsGp)
{
    std::pair<double,double> measureBS;
    TrajectoryStateClosestToPoint pion_BeamSpot = piTT.trajectoryStateClosestToPoint(BsGp);
    if(pion_BeamSpot.isValid())
    {
      measureBS.first = pion_BeamSpot.perigeeParameters().transverseImpactParameter();
      if(pion_BeamSpot.hasError() && !(pion_BeamSpot.hasError()==0)) 
      {
        measureBS.second = measureBS.first/pion_BeamSpot.perigeeError().transverseImpactParameterError();
      }
	}
	else
	{
	measureBS.first  = 99999;
	measureBS.second = 99999;
	}
	
	return measureBS;      

}

