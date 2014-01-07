#ifndef Bc2JpsiPi_h
#define Bc2JpsiPi_h

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

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include <DataFormats/RecoCandidate/interface/RecoChargedCandidate.h>

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "HeavyFlavorAnalysis/Bc2JpsiPi/interface/BcTree.h"
#include "UsefulTools.h"



class Bc2JpsiPi : public edm::EDAnalyzer {

   public:
      
      typedef std::vector<std::pair<std::pair<double,double>, std::pair<double,double> > > JTrigDef ;
      typedef std::map<uint, reco::TrackRef  >                                         newTracksDef ;  
 	  typedef ROOT::Math::SVector<double, 3>    											Vector3 ;
  	  typedef ROOT::Math::SMatrix<double, 3, 3> 										   Matrix33 ;
 
      explicit	    			Bc2JpsiPi(            const edm::ParameterSet&);
              	   				~Bc2JpsiPi();


   private:
      virtual void  			beginJob() ;
      virtual void  			analyze(                const edm::Event&, 
                                           				const edm::EventSetup&  							);
      void    	    			acquireInputParameters(	void                   								);
      void          			MonteCarloStudies(     	const edm::Event&       							);
      bool          			passPioncuts(			reco::TrackRef pionCand								);
      void          			checkHLT(				const edm::Event& iEvent							);
      
      virtual void  			endJob() ;
      
// ----------member data ---------------------------

	  TFile  *file;
      BcTree *BcTree_;
	
      // These declarations must be in this order to avoid compiler warnings
	  bool boolTrg_;  
      std::vector<std::string> inputDouble_;
      std::vector<std::string> inputString_;
      edm::InputTag 	   thePVs_;
      edm::InputTag 	   thebeamspot_;
	  // End of mandatory alphabetic order 
	
	  int runNumber;
	
      std::map<std::string, TH1F*  	> h1_ 		;
      std::map<std::string, TH2F*  	> h2_ 		;
      std::map<std::string, TH1F*  	> g1_ 		;
      std::map<std::string, double 	> pDouble_ 	;
      std::map<std::string, std::string> pString_ ;
      std::map<std::string, TTree* 	> t1_ 		;
       
      edm::ParameterSet thisConfig_ ;
      unsigned int      nEvents_    ;

      JTrigDef JPsiTriggered_ ;
      
      UsefulTools* Utils;

};

#endif
