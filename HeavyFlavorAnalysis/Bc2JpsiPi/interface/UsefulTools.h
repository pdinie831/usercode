#ifndef USEFULTOOLS_H
#define USEFULTOOLS_H

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

class UsefulTools
{


 public:
  
  UsefulTools();
  ~UsefulTools() {};
  

double computeCosine (double Lx, 
					  double Ly,
					  double Lz,
					  double Px,
					  double Py,
					  double Pz);

double computeCosineError   (double Lx,
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
							 double PyzCov);

std::pair<double,double> 	pionImpactParameter(reco::TransientTrack piTT, TransientVertex jpsiVtx	);
std::pair<double,double> 	pionImpactParameter(reco::TransientTrack piTT, reco::Vertex myVtx	);
std::pair<double,double> 	pionIPBeamSpot(reco::TransientTrack piTT, GlobalPoint BsGp				);

  
};

#endif