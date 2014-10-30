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
  double cosineError = 0. ;
  double Lnorm = sqrt(Lx*Lx + Ly*Ly + Lz*Lz);
  double Pnorm = sqrt(Px*Px + Py*Py + Pz*Pz);
  double LdotW = computeCosine(Lx, Ly, Lz, Px, Py, Pz);
  
  if ((Lnorm > 0.) && (Pnorm > 0.))
  {
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
    TrajectoryStateClosestToBeamLine pion_BeamSpot = piTT.stateAtBeamLine();
    Measurement1D meas = pion_BeamSpot.transverseImpactParameter();
    measureBS.first  = meas.value();
    measureBS.second = meas.value()/meas.error();
	
    return measureBS;      

}



//Compute longitudinal impact parameter wrt vertex
std::pair<double,double> UsefulTools::LongitudinalIP(reco::TransientTrack transientTrk, reco::Vertex myVtx)
{
    GlobalPoint VtxGP(myVtx.x(), myVtx.y(), myVtx.z()); 
    std::pair<double,double> measureL;
    TrajectoryStateClosestToPoint TSC = transientTrk.trajectoryStateClosestToPoint(VtxGP);
    if(TSC.isValid())
    {
      measureL.first = TSC.perigeeParameters().longitudinalImpactParameter();
      if(TSC.hasError() && !(TSC.hasError()==0)) 
      {
        measureL.second = measureL.first/TSC.perigeeError().longitudinalImpactParameterError();
      }
	}
	else
	{
	measureL.first  = 99999;
	measureL.second = 99999;
	}
	
	return measureL;      

}


//Compute impact parameter 2D significance wrt any vtx
double UsefulTools::piond0IPSignificanceVtx(reco::TransientTrack piTT, GlobalPoint vert)
{
    TrajectoryStateClosestToPoint pion_vGp = piTT.trajectoryStateClosestToPoint(vert);
    double d0 = pion_vGp.perigeeParameters().transverseImpactParameter();
    double d0_error = pion_vGp.perigeeError().transverseImpactParameterError();
    double significance = d0/d0_error;
	
    return significance;      
}

double UsefulTools::piondzIPSignificanceVtx(reco::TransientTrack piTT, GlobalPoint vert)
{
    TrajectoryStateClosestToPoint pion_vGp = piTT.trajectoryStateClosestToPoint(vert);
    double dz = pion_vGp.perigeeParameters().longitudinalImpactParameter();
    double dz_error = pion_vGp.perigeeError().longitudinalImpactParameterError();
    double significance = dz/dz_error;
	
    return significance;      
}




double UsefulTools::computeInvariantMass (double Px1, double Py1, double Pz1, double mass1,
                                          double Px2, double Py2, double Pz2, double mass2,
                                          double Px3, double Py3, double Pz3, double mass3
                                         )
{
  double Energy1 = sqrt(Px1*Px1 + Py1*Py1 + Pz1*Pz1 + mass1*mass1);
  double Energy2 = sqrt(Px2*Px2 + Py2*Py2 + Pz2*Pz2 + mass2*mass2);
  double Energy3 = sqrt(Px3*Px3 + Py3*Py3 + Pz3*Pz3 + mass3*mass3);
  return sqrt((Energy1+Energy2+Energy3) * (Energy1+Energy2+Energy3) - ((Px1+Px2+Px3) * (Px1+Px2+Px3) + (Py1+Py2+Py3) * (Py1+Py2+Py3) + (Pz1+Pz2+Pz3) * (Pz1+Pz2+Pz3)));
}


// double UsefulTools::myDistance(const reco::TransientTrack & transientTrack, reco::Vertex & vertex){
//       AnalyticalImpactPointExtrapolator extrapolator(transientTrack.field());
//       TrajectoryStateOnSurface myTsos = extrapolator.extrapolate(transientTrack.impactPointState(), RecoVertex::convertPos(vertex.position());
//  return 0.44;     
// }      
//       VertexDistance3D dist;
//       return absoluteImpactParameter(extrapolator.extrapolate(transientTrack.impactPointState(), RecoVertex::convertPos(vertex.position())), vertex, dist);
// 
//     std::pair<bool,Measurement1D> absoluteImpactParameter(const TrajectoryStateOnSurface & tsos , const reco::Vertex & vertex, VertexDistance & distanceComputer) {
//         if(!tsos.isValid()) {
//          return pair<bool,Measurement1D>(false,Measurement1D(0.,0.)) ;
//         }
//         GlobalPoint refPoint = tsos.globalPosition();
//         GlobalError refPointErr = tsos.cartesianError().position();
//         GlobalPoint vertexPosition = RecoVertex::convertPos(vertex.position());
//         GlobalError vertexPositionErr = RecoVertex::convertError(vertex.error());
//         return pair<bool,Measurement1D>(true,distanceComputer.distance(VertexState(vertexPosition,vertexPositionErr), VertexState(refPoint, refPointErr)));
//    }
