#ifndef Bc2Jpsi3Pi_h
#define Bc2Jpsi3Pi_h

#include <map>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <vector>

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TMath.h>
#include <TObject.h>
#include <TTree.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TObjArray.h>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Vector.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"

#include "Math/GenVector/PxPyPzE4D.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParametersError.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include <DataFormats/RecoCandidate/interface/RecoChargedCandidate.h>
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "HeavyFlavorAnalysis/Bc2Jpsi3Pi/interface/BcTree.h"


class Bc2Jpsi3Pi : public edm::EDAnalyzer {
   public:
      
      typedef std::vector<std::pair<std::pair<double,double>, std::pair<double,double> > > JTrigDef ;
      typedef std::map<uint, reco::TrackRef  >                                         newTracksDef ;  
 
      explicit	    			Bc2Jpsi3Pi(            const edm::ParameterSet&);
              	   				~Bc2Jpsi3Pi();

   private:
      virtual void  			beginJob() ;
      virtual void  			analyze(                const edm::Event&, 
                                           				const edm::EventSetup&  							);
      void    	    			acquireInputParameters(	void                   								);
      void          			MonteCarloStudies(     	const edm::Event&       							);
      bool          			passPioncuts(			reco::TrackRef pionCand								);
      std::pair<double,double> 	pionImpactParameter(	reco::TransientTrack piTT, TransientVertex jpsiVtx	);
      std::pair<double,double> 	pionIPBeamSpot(			reco::TransientTrack piTT, GlobalPoint BsGp			);
      void          			checkHLT(				const edm::Event& iEvent							);
      
      virtual void  			endJob() ;
      
// ----------member data ---------------------------
	
	TFile * file;
	BcTree *BcTree_;

	// These declarations must be in this order to avoid compiler warnings
	std::vector<std::string> inputDouble_;
	std::vector<std::string> inputString_;
	edm::InputTag 		 thePVs_;
	edm::InputTag 		 thebeamspot_;
	// End of mandatory alphabetic order 

	int runNumber;
	
	std::map<std::string, TH1F* 	 > h1_ ;
	std::map<std::string, TH2F* 	 > h2_ ;
	std::map<std::string, TH1F*      > g1_ ;
	std::map<std::string, double	 > pDouble_ ;
	std::map<std::string, std::string> pString_ ;
	std::map<std::string, TTree*     > t1_ ;
       
	edm::ParameterSet thisConfig_ ;
	unsigned int             nEvents_;

	JTrigDef JPsiTriggered_ ;

};

#endif
