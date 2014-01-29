// -*- C++ -*-
// Package:    BcAnalysis
// Class:      BcAnalysis
// 
/**\class  BcAnalysis BcAnalysis.cc HeavyFlavor/BcAnalysis/src/BcAnalysis.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
// Original Author:  Silvia Taroni & Sandra Malvezzi
// Put in local cvs
// Added in cvs milano
//         Created:  Thu Apr 23 16:09:09 CEST 2009
// $Id: BcAnalysis.cc,v 1.8 2010/05/26 09:34:53 menasce Exp $
//

// system include files
#include <memory>
#include <boost/regex.hpp> 

#include <typeinfo>

// user include files
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"   
#include "HepMC/GenParticle.h"

#include <DataFormats/Math/interface/deltaR.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CLHEP/Vector/LorentzVector.h"

#include "HeavyFlavorAnalysis/Bc2Jpsi3Pi/interface/Bc2Jpsi3Pi.h"

// constants, enums and typedefs
double c_const  = 0.0299792458;//cm/ps
const ParticleMass muon_mass    = 0.10565837; 
ParticleMass jpsimass_c 	  	= 3.096916; 
ParticleMass pion_mass  	  	= 0.13957018; 
int JpsiId      				= 443; 

TLorentzVector p_muP,p_muM, p_pi1, p_pi2, p_pi3, p_jpsi,p_Bc;

double ELLE       = 0.;
double TAU        = 0.;
double COSINE     = 0.;
double PVVTX[3]   = {0,0,0};
double BCVTX[3]   = {0,0,0};
double JPSIVTX[3] = {0,0,0};
double BSPOT[3]   = {0,0,0};
int PICH[3]       = {0,0,0};
bool boolMC       = false;


using namespace std ;

// constructors and destructor
Bc2Jpsi3Pi::Bc2Jpsi3Pi(const edm::ParameterSet& iConfig) :
   boolTrg_    (iConfig.getParameter<bool>                     ("printTriggers"   )),
   inputDouble_(iConfig.getParameter<std::vector<std::string> >("inputDouble_"    )),
   inputString_(iConfig.getParameter<std::vector<std::string> >("inputString_"    )),
   thePVs_(     iConfig.getParameter<edm::InputTag>	           ("primaryVertexTag")),
   thebeamspot_(iConfig.getParameter<edm::InputTag>	           ("beamSpotTag"     ))
{
  cout << __LINE__ << "]\t" << __PRETTY_FUNCTION__ << "\tCtor" << endl ;
  thisConfig_ = iConfig;
  acquireInputParameters() ;
  Utils = new UsefulTools();
}	

//==============================================================================================
Bc2Jpsi3Pi::~Bc2Jpsi3Pi()
{
  delete Utils;
}

//=============== method called to for each event ===================================================================
void Bc2Jpsi3Pi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  BcTree_  ->BcTree::BcTreeClear();
  runNumber++;
  nEvents_++;

 // Handles to collections 
  edm::Handle<reco::BeamSpot>                    theBeamSpot;
  iEvent.getByLabel(thebeamspot_,                theBeamSpot	);
  reco::BeamSpot bs = *theBeamSpot;

  edm::Handle<reco::MuonCollection>   			 recoMuons      ; 
  iEvent.getByLabel("muons",          			 recoMuons)		;

  edm::Handle<reco::TrackCollection>             tracks 		;
  iEvent.getByLabel(pString_["trackCollection"], tracks)		;

  edm::Handle<reco::VertexCollection>            priVtxs		;
  iEvent.getByLabel(thePVs_,                     priVtxs)		;
  const reco::VertexCollection primvtx 	= *(priVtxs.product())  ;    

  edm::ESHandle<GlobalTrackingGeometry>          theTrackingGeometry;
  edm::ESHandle<TransientTrackBuilder>           theB;

  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  edm::Handle<edm::TriggerResults> 				 hltresults;
  edm::InputTag tag("TriggerResults","","HLT");
  iEvent.getByLabel(tag,           hltresults);

  //Call GEN level analysis
  if (!iEvent.isRealData())
  {
	MonteCarloStudies(iEvent) ;
  }

  BSPOT[0]	= bs.x0(); 
  BSPOT[1]	= bs.y0(); 
  BSPOT[2]	= bs.z0(); 
  GlobalPoint BeamSpotGP(bs.x0(), bs.y0(), bs.z0()); 

  if(!iEvent.isRealData())  
	BcTree_->BcTree::SetBcTreeGENCand( p_Bc 	,
									   p_muP 	,
									   p_muM 	,
									   p_jpsi	,
									   p_pi1 	,
									   p_pi2 	,
									   p_pi3 	,
									   ELLE     ,
									   COSINE   ,
									   TAU      ,
									   PVVTX    ,
									   BCVTX    ,
									   JPSIVTX  ,
									   BSPOT    ,
									   PICH
									   );


  // internal variables 
  int 		counter 	= 0 ;
  int       NumBc   	= 0 ;
  int 		numJpsi 	= 0 ;
  int 		JpsiNoTrig 	= 0 ; 
  int 		JpsiTrig   	= 0 ; 

  float 	Tnpv 		= 0 ;

  double 	jpsiX   	= 0 ;
  double 	jpsiY   	= 0 ;
  double 	jpsiZ   	= 0 ;
  double    DCA			= 0 ;

  double 	BSCovariance[9]		;
  double    HighestPtPVPos[3] 	;
  double    HighestPtPVCov[9]	;

  double 	SVPos[3]			;
  double 	SVCov[9]			;
  double    PointPVPos[3]		;
  double 	PointPVCov[9]		;
  double 	JpsiVtxPos[3]		;

  double 	HltMatch[4]	= {0,0,0,0}	;													     	
  double 	matchMu[2]	= {0,0}		;
  double 	matchPi1[2] = {0,0}		;
  double 	matchPi2[2] = {0,0}		;
  double 	matchPi3[2] = {0,0}		;
  
  float     TrackMuPDR[3];
  float     TrackMuMDR[3];

  const reco::Muon* recomu1   = 0;
  const reco::Muon* recomu2   = 0;
  const reco::Muon* recomuP   = 0;
  const reco::Muon* recomuM   = 0;
  
  std::vector<int>      goodTrig ;

  Matrix33 Covsec, Covprim, Covbs, CovPointPV, CovJpsi;
  Vector3  Deriv;

  //Set event variables
  HighestPtPVPos[0] = primvtx[0].position().x();
  HighestPtPVPos[1] = primvtx[0].position().y();
  HighestPtPVPos[2] = primvtx[0].position().z();
  
  int n = 0;
  for (int i=0; i<3; i++)
  {
	for (int j=0; j<3; j++)
	{ 
	  Covbs(i,j)        = bs.covariance(i,j);  
	  Covprim(i,j)      = primvtx[0].covariance(i,j);  
	  BSCovariance[n]   = Covbs(i,j);
	  HighestPtPVCov[n] = Covprim(i,j);
	  n++;
	}
  }

  //Get info for MC PileUp reweighting 
  if(!iEvent.isRealData())
  {
	edm::Handle<std::vector< PileupSummaryInfo > > PupInfo;
	iEvent.getByLabel("addPileupInfo", PupInfo);
	std::vector<PileupSummaryInfo>::const_iterator PVI;
	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) 
	{
	  int BX = PVI->getBunchCrossing();
	  if(BX == 0) 
	  { 
		Tnpv = PVI->getTrueNumInteractions();
		continue;
	  }
	}
	h1_["TrueNumInteraction"  ]->Fill(Tnpv);
  }

  //Check trigger 
  if (!hltresults.isValid() ) 
  {
	cout << __LINE__ << "\tcould not get the collection!" << endl;
	return;
  }

  edm::TriggerNames triggerNames_ = iEvent.triggerNames(*hltresults);
  for(int itrig = 0; itrig != (int)hltresults->size(); ++itrig)
  { 
     if (boolTrg_) 	std::cout << triggerNames_.triggerName(itrig) << std::endl;
	 if((triggerNames_.triggerName(itrig) == pString_["HLTname1"]) || 
		(triggerNames_.triggerName(itrig) == pString_["HLTname2"]) ||
		(triggerNames_.triggerName(itrig) == pString_["HLTname3"]) ||
		(triggerNames_.triggerName(itrig) == pString_["HLTname4"]) 
	 ) 
	{
	  goodTrig.push_back(itrig) ;
    }
  }
  if ( goodTrig.size() == 0 ) {
	cout << __LINE__ << "\tthe required trigger paths are not present is this dataset!" << endl;
	std::cout << "Looking for " << pString_["HLTname1"] << "," << pString_["HLTname2"] << "," << pString_["HLTname3"] << "," << pString_["HLTname4"]  << std::endl;
	return ;
  }
  bool triggerAccepted = false ;
  for(unsigned int t=0; t<goodTrig.size(); ++t)
  {
	if( hltresults->accept(goodTrig[t]) )   triggerAccepted = true ;
  }

  //Pre-selection cuts on track collection
  newTracksDef qualityTracks;
  for (uint tracksIt =0 ;  tracksIt < tracks->size(); tracksIt++)
  {
     reco::TrackRef checkTrk(tracks,tracksIt) ;										        
	 if (! checkTrk->quality(reco::TrackBase::highPurity))		continue;
     bool flag = false ; 												        
     for (reco::MuonCollection::const_iterator mu =recoMuons->begin(); mu!=recoMuons->end(); mu++)			        
     {														        
       if (mu->track().get()!= 0 && mu->track().get() == checkTrk.get()) { flag=true; break;}
//        if (mu->innerTrack().get()!= 0 && deltaR(mu->innerTrack().get()->eta(), mu->innerTrack().get()->phi(), checkTrk->eta(), checkTrk->phi()) < 0.001)	std::cout <<" can be a muon!!" << std::endl;		        
     }														        
     if (flag)												   											continue;	        
     if (checkTrk->pt()  				  					< pDouble_["cut_Pt_Trk"])			   		continue;	        
     if (fabs(checkTrk->eta())			  					> pDouble_["cut_eta_pi"])			  		continue;	   						  
//      if (checkTrk->numberOfValidHits()	  					< pDouble_["cut_nhits"] )			   		continue;	   					   
//      if (checkTrk->hitPattern().numberOfValidPixelHits() 	< pDouble_["numberOfValidPixelHitsCut"]  )	continue;	        
//      if (checkTrk->normalizedChi2()			 			    > pDouble_["cut_chi2n"])			   		continue; 	        
     qualityTracks[tracksIt]=checkTrk;											   	
  }
  BcTree_->BcTree::SetBcTreeHeader(  nEvents_		       					,
									(int32_t)iEvent.id().run()            	,
									(int32_t)iEvent.id().luminosityBlock()	,
									(int32_t)iEvent.id().event()	       	,
									(int32_t)tracks ->size()              	,
									(int32_t)priVtxs->size()              	,
									qualityTracks.size()					,
									Tnpv                                  	,
									BSPOT                                  	,
									BSCovariance							,
									HighestPtPVPos							,
									HighestPtPVCov							,
									NumBc		                           
							    );

   h1_["filter"      ]->Fill(counter++);//0
   if ( !triggerAccepted ) {
     if (boolMC) t1_["BcTree"]->Fill();
     return ;
   }
  
   this->checkHLT(iEvent) ;  

   h1_["filter"      ]->Fill(counter++);//1
   if (recoMuons->size()<2){
     if (boolMC) t1_["BcTree"]->Fill();
     return ;
   }
 
   h1_["pVertexSize"        ]-> Fill(primvtx.size()			);
   h1_["pVertexTrackSize"	]-> Fill(primvtx[0].tracksSize());
   h1_["EventTrackSize"		]-> Fill(tracks->size() 		);
   h1_["goodTrkSize"		]-> Fill(qualityTracks.size() 	);
   if(qualityTracks.size() < 3 )
   {
	 std::cout << "exit for track size < 3 " << std::endl; 
     if (boolMC) t1_["BcTree"]->Fill();
     return; 
   }

//---------------------------------------------------------------------------------------------------
  //Loop on the first muon

  for(reco::MuonCollection::const_iterator muR1 =recoMuons->begin(); muR1!=recoMuons->end(); ++muR1)
  {
	recomu1 = &*muR1;
	
	if (fabs(muR1->eta()) > pDouble_["cut_eta"]) 							         				continue;
	if (fabs(muR1->pt() ) < pDouble_["cut_Pt_Mu"])							         				continue;
	
    //Loop on the second muon 
    for(reco::MuonCollection::const_iterator muR2 =muR1+1; muR2!=recoMuons->end(); ++muR2)
    {
	  recomu2 = &*muR2;

	  if (muR1->charge() + muR2->charge() != 0)								   						continue; 
	  if (fabs(muR2->eta()) > pDouble_["cut_eta"])								    				continue;
	  if (fabs(muR2->pt() ) < pDouble_["cut_Pt_Mu"])												continue;
	
	  if( !recomu1->track().isNonnull() ||   !recomu2->track().isNonnull())							continue ;

	  recomu1 = &*muR1;
	  recomu2 = &*muR2;

	  reco::TrackRef refMuP, refMuM, innerTrkMuP, innerTrkMuM;
	  if (muR1->charge() == 1){
	    refMuP=recomu1->track(); 
	    refMuM=recomu2->track();
	    innerTrkMuP = recomu1->innerTrack(); 
	    innerTrkMuM = recomu2->innerTrack(); 
	    recomuP = recomu1;
	    recomuM = recomu2;
	  }
	  else{
	    refMuM=recomu1->track(); 
	    refMuP=recomu2->track();
	    innerTrkMuM = recomu1->innerTrack(); 
	    innerTrkMuP = recomu2->innerTrack(); 
	    recomuM = recomu1;
	    recomuP = recomu2;
	  }

      const reco::TransientTrack tTrackmuP=(*theB).build(refMuP.get());
      const reco::TransientTrack tTrackmuM=(*theB).build(refMuM.get());

	  reco::Candidate::LorentzVector MuPp4(refMuP->px(), 
										   refMuP->py(),
										   refMuP->pz(),
										   sqrt(refMuP->p()*refMuP->p() + recomu1->mass()*recomu1->mass()));
	  reco::Candidate::LorentzVector MuMp4(refMuM->px(),
										   refMuM->py(),
										   refMuM->pz(),
										   sqrt(refMuM->p()*refMuM->p() + recomu2->mass()*recomu2->mass()));

      //Trigger requirements - dr
	  if (fabs( (- (refMuP->vx()-bs.x0()) * refMuP->py() + (refMuP->vy()-bs.y0()) * refMuP->px() ) / refMuP->pt() ) > 2.0 ) continue;
	  if (fabs( (- (refMuM->vx()-bs.x0()) * refMuM->py() + (refMuM->vy()-bs.y0()) * refMuM->px() ) / refMuM->pt() ) > 2.0 ) continue;

	  reco::Candidate::LorentzVector jpsiV = MuPp4 + MuMp4;
	  h1_["InvMassJPsi" ]-> Fill(jpsiV.M() );
	  h1_["JpsiPt"      ]-> Fill(jpsiV.pt()); 

      //Trigger requirements - DCA
	  TrajectoryStateClosestToPoint muPTS = tTrackmuP.impactPointTSCP();
	  TrajectoryStateClosestToPoint muMTS = tTrackmuM.impactPointTSCP();
	  if (muPTS.isValid() && muMTS.isValid()) 
	  {
		ClosestApproachInRPhi cApp;
		cApp.calculate(muPTS.theState(), muMTS.theState());
		if (!cApp.status() || cApp.distance() > 0.5) continue;
		DCA = cApp.distance();
	  }
	  h1_["filter"      ]->Fill(counter++);//2

	  //Fit to the dimuon vertex
	  KalmanVertexFitter 			fitterMuMu;
	  TransientVertex    			mumuVertex;
	  vector<reco::TransientTrack> 	tTrMu;
	  tTrMu.push_back(tTrackmuP);										  
	  tTrMu.push_back(tTrackmuM);
	  mumuVertex = fitterMuMu.vertex(tTrMu);
	  if(!mumuVertex.isValid()) 													continue;
	 
	  //Evaluate Jpsi L/Sigma and cosine
	  math::XYZVector pperp(refMuP->px() + refMuM->px(), refMuP->py() + refMuM->py(), 0.);
	  GlobalError JpsiVertexError = mumuVertex.positionError();
	  GlobalPoint displacementFromBeamspot(-1*((bs.x0() - mumuVertex.position().x()) + (mumuVertex.position().z() - bs.z0()) * bs.dxdz()),
										   -1*((bs.y0() - mumuVertex.position().y()) + (mumuVertex.position().z() - bs.z0()) * bs.dydz()), 0);
	 
	  double LxyJpsi 	    = displacementFromBeamspot.perp();
	  double LxyErrJpsi 	= sqrt(JpsiVertexError.rerr(displacementFromBeamspot));
	  reco::Vertex::Point     vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
	  double cosJpsiXY 		= vperp.Dot(pperp)/(vperp.R()*pperp.R());
	  double elsigJpsi  	= 0 ;
	  try	      { elsigJpsi= LxyJpsi / LxyErrJpsi ; }
	  catch (...) {cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t" << "Floating divide by zero!" << endl ;}

	  h1_["ElsigJpsi"   ]-> Fill(elsigJpsi); 
	  h1_["filter"      ]-> Fill(counter++);//3

      //Trigger requirements: CL and pT
	  float ClJpsiVtx(TMath::Prob(mumuVertex.totalChiSquared(),(int)(mumuVertex.degreesOfFreedom())));
	  if (ClJpsiVtx  <     pDouble_["cut_cl_Jpsi"])            						continue;
	  if (jpsiV.pt() <=    pDouble_["ptJpsi"])										continue;
	  if (fabs(jpsiV.M() - pDouble_["JPsiMassPDG"]) > pDouble_["JPsiMassWindow"]) 	continue;

	  h1_["InvMassJPsiCuts"   ]->Fill(jpsiV.M());

      //Trigger Matching
	  JpsiNoTrig++ ; 
	  bool   MuTrigMatch  = false;		
	  for(JTrigDef::iterator it=JPsiTriggered_.begin(); it!=JPsiTriggered_.end(); ++it)								     	
	  {																		     	
		double etaMu1 = (*it).first.first   ;														     	
		double phiMu1 = (*it).first.second  ;														     	
		double etaMu2 = (*it).second.first  ;														     	
		double phiMu2 = (*it).second.second ;														     	
		if( ((reco::deltaR(refMuP->eta(), refMuP->phi(), etaMu1, phiMu1 ) < pDouble_["HLTMatch"]) 	&&								     	
			 (reco::deltaR(refMuM->eta(), refMuM->phi(), etaMu2, phiMu2 ) < pDouble_["HLTMatch"])) 	||								     	
			((reco::deltaR(refMuP->eta(), refMuP->phi(), etaMu2, phiMu2 ) < pDouble_["HLTMatch"]) 	&&								     	
			 (reco::deltaR(refMuM->eta(), refMuM->phi(), etaMu1, phiMu1 ) < pDouble_["HLTMatch"]))									     	
		  ) 																		     	
		{																		     	
		  JpsiTrig++; 
		  MuTrigMatch = true;	
		  HltMatch[0] = reco::deltaR(refMuP->eta(), refMuP->phi(), etaMu1, phiMu1 );
		  HltMatch[1] = reco::deltaR(refMuM->eta(), refMuM->phi(), etaMu2, phiMu2 );
		  HltMatch[2] = reco::deltaR(refMuP->eta(), refMuP->phi(), etaMu2, phiMu2 );
		  HltMatch[3] = reco::deltaR(refMuM->eta(), refMuM->phi(), etaMu1, phiMu1 );
		}																		     	
	  }																		     	
  
	  h1_["JpsiBeforeTrigMatch" 		]->Fill(JpsiNoTrig);
	  if(!MuTrigMatch )					continue;
	  h1_["JpsiAfterTrigMatch"  		]->Fill(JpsiTrig  );
	  h1_["InvMassJpsiPassingTrigMatch" ]->Fill(jpsiV.M()  );

	  TLorentzVector *jpsiTV   = new TLorentzVector(jpsiV.px(), jpsiV.py(), jpsiV.pz(), jpsiV.energy());
	  JpsiVtxPos[0] = mumuVertex.position().x() ;
	  JpsiVtxPos[1] = mumuVertex.position().y() ;
	  JpsiVtxPos[2] = mumuVertex.position().z() ;

	  //Loop on tracks
	  for (newTracksDef::const_iterator it_pi1=qualityTracks.begin(); it_pi1!=qualityTracks.end(); ++it_pi1)
	  {
		reco::TrackRef firstPionTrackCand((*it_pi1).second) ;
		const reco::TransientTrack pi1TrackCand = (*theB).build(firstPionTrackCand);
		
		TrackMuPDR[0] = deltaR(innerTrkMuP->eta(), innerTrkMuP->phi(), firstPionTrackCand->eta(), firstPionTrackCand->phi());
		TrackMuMDR[0] = deltaR(innerTrkMuM->eta(), innerTrkMuM->phi(), firstPionTrackCand->eta(), firstPionTrackCand->phi());
		
		pair<double,double>  pion1_IPPair = Utils->pionImpactParameter(pi1TrackCand,mumuVertex);
		double Pi1_3DIp                   = pion1_IPPair.first;
		double Pi1_3DIpSignificance       = pion1_IPPair.second;
		h1_["Pion1_ImpactParameter"]->Fill(Pi1_3DIp);
		h1_["Pion1_IPsignificance" ]->Fill(Pi1_3DIpSignificance);
		
		pair<double,double>  pion1_BSPair = Utils->pionIPBeamSpot(pi1TrackCand,BeamSpotGP);
		double Pi1_IPbs                   = pion1_BSPair.first;
		double Pi1_IPbsSignificance       = pion1_BSPair.second;
		h1_["Pion1_TransverseImpactParameter_BS" ]-> Fill(Pi1_IPbs);
		h1_["Pion1_TIParameterBSsignificance"    ]-> Fill(Pi1_IPbsSignificance);

		pair<double,double>  pion1_IPPairPV = Utils->pionImpactParameter(pi1TrackCand,primvtx[0]);
		double Pi1_3DIpPV                   = pion1_IPPairPV.first;
		double Pi1_3DIpPVSignificance       = pion1_IPPairPV.second;
		//Require small impact parameter significance wrt Jpsi vertex
		if(Pi1_3DIpPVSignificance > 10) 								continue;

		//Require trk in a cone wrt Jpsi track
		double DRPi1 =  deltaR(firstPionTrackCand->eta(),firstPionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
		h1_["deltaRPi1"     ]-> Fill(DRPi1);
		if (DRPi1 > 5) continue;

		//Loop on second pion
		for (newTracksDef::const_iterator it_pi2=it_pi1; it_pi2!=qualityTracks.end(); ++it_pi2)
		{
		  if( (*it_pi1).first == (*it_pi2).first ) 						continue ;
		  reco::TrackRef secondPionTrackCand((*it_pi2).second) ;
		  const reco::TransientTrack pi2TrackCand = (*theB).build(secondPionTrackCand);

		  TrackMuPDR[1] = deltaR(innerTrkMuP->eta(), innerTrkMuP->phi(), secondPionTrackCand->eta(), secondPionTrackCand->phi());
		  TrackMuMDR[1] = deltaR(innerTrkMuM->eta(), innerTrkMuM->phi(), secondPionTrackCand->eta(), secondPionTrackCand->phi());

		  pair<double,double>  pion2_IPPair = Utils->pionImpactParameter(pi2TrackCand,mumuVertex);
		  double Pi2_3DIp                   = pion2_IPPair.first;
		  double Pi2_3DIpSignificance       = pion2_IPPair.second;
		  h1_["Pion2_ImpactParameter"]-> Fill(Pi2_3DIp);
		  h1_["Pion2_IPsignificance" ]-> Fill(Pi2_3DIpSignificance);
	  
		  pair<double,double>  pion2_BSPair = Utils->pionIPBeamSpot(pi2TrackCand,BeamSpotGP);
		  double Pi2_IPbs                   = pion2_BSPair.first;
		  double Pi2_IPbsSignificance       = pion2_BSPair.second;
		  h1_["Pion2_TransverseImpactParameter_BS"]-> Fill(Pi2_IPbs);
		  h1_["Pion2_TIParameterBSsignificance"   ]-> Fill(Pi2_IPbsSignificance);

		  pair<double,double>  pion2_IPPairPV = Utils->pionImpactParameter(pi2TrackCand,primvtx[0]);
		  double Pi2_3DIpPV                   = pion2_IPPairPV.first;
		  double Pi2_3DIpPVSignificance       = pion2_IPPairPV.second;
		  //Require small impact parameter significance wrt Jpsi vertex
		  if(Pi2_3DIpPVSignificance > 10) 								continue;
		  
		  //Require trk in a cone wrt Jpsi track
		  double DRPi2 =  deltaR(secondPionTrackCand->eta(),secondPionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
		  h1_["deltaRPi2"     ]-> Fill(DRPi2);
		  if (DRPi2 > 5) 												continue;

		  //Loop on third pion
		  for (newTracksDef::const_iterator it_pi3=it_pi2; it_pi3!=qualityTracks.end(); ++it_pi3)
		  {
			if( (*it_pi3).first == (*it_pi2).first ) 					continue ;
			reco::TrackRef thirdPionTrackCand((*it_pi3).second) ;
			const reco::TransientTrack pi3TrackCand = (*theB).build(thirdPionTrackCand);

			TrackMuPDR[2] = deltaR(innerTrkMuP->eta(), innerTrkMuP->phi(), thirdPionTrackCand->eta(), thirdPionTrackCand->phi());
			TrackMuMDR[2] = deltaR(innerTrkMuM->eta(), innerTrkMuM->phi(), thirdPionTrackCand->eta(), thirdPionTrackCand->phi());

			pair<double,double>  pion3_IPPair = Utils->pionImpactParameter(pi3TrackCand,mumuVertex);
			double Pi3_3DIp                   = pion3_IPPair.first;
			double Pi3_3DIpSignificance       = pion3_IPPair.second;
			h1_["Pion3_ImpactParameter"]-> Fill(Pi3_3DIp);
			h1_["Pion3_IPsignificance" ]-> Fill(Pi3_3DIpSignificance);

			pair<double,double>  pion3_BSPair = Utils->pionIPBeamSpot(pi3TrackCand,BeamSpotGP);
			double Pi3_IPbs                   = pion3_BSPair.first;
			double Pi3_IPbsSignificance       = pion3_BSPair.second;
			h1_["Pion3_TransverseImpactParameter_BS"]-> Fill(Pi3_IPbs);
			h1_["Pion3_TIParameterBSsignificance"   ]-> Fill(Pi3_IPbsSignificance);

			pair<double,double>  pion3_IPPairPV = Utils->pionImpactParameter(pi3TrackCand,primvtx[0]);
			double Pi3_3DIpPV                   = pion3_IPPairPV.first;
			double Pi3_3DIpPVSignificance       = pion3_IPPairPV.second;
			//Require small impact parameter significance wrt Jpsi vertex
			if(Pi3_3DIpPVSignificance> 10) 									continue;

			//Require trk in a cone wrt Jpsi track
			double DRPi3 =  deltaR(thirdPionTrackCand->eta(),thirdPionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
			h1_["deltaRPi3"     ]-> Fill(DRPi3);
			if (DRPi3 > 5) continue;

			if(iEvent.isRealData())
			{  	  
			 if (fabs(firstPionTrackCand ->charge()+
				   secondPionTrackCand->charge()+
				   thirdPionTrackCand ->charge()) == 1 )					continue ;
			 }
			if(!iEvent.isRealData())
			{  	  
			 if (firstPionTrackCand ->charge()+
				   secondPionTrackCand->charge()+
				   thirdPionTrackCand ->charge() != 1 )					    continue ;
			}
	  
			KinematicParticleFactoryFromTransientTrack pFactory;
			float muon_sigma 		= muon_mass*1.e-6;
			float pion_sigma        = pion_mass*1.e-6; 
			float chi 		 		= 0.;
			float ndf 		 		= 0.;

			std::vector<RefCountedKinematicParticle> XParticles_Bc;
			XParticles_Bc.push_back(pFactory.particle(tTrackmuP,muon_mass,chi,ndf,muon_sigma));
			XParticles_Bc.push_back(pFactory.particle(tTrackmuM,muon_mass,chi,ndf,muon_sigma));
			XParticles_Bc.push_back(pFactory.particle(pi1TrackCand,pion_mass,chi,ndf,pion_sigma)); 
			XParticles_Bc.push_back(pFactory.particle(pi2TrackCand,pion_mass,chi,ndf,pion_sigma)); 
			XParticles_Bc.push_back(pFactory.particle(pi3TrackCand,pion_mass,chi,ndf,pion_sigma)); 

			MultiTrackKinematicConstraint *  jpsiconstraint = new  TwoTrackMassKinematicConstraint(jpsimass_c);
			KinematicConstrainedVertexFitter kcvFitterBc;
			RefCountedKinematicTree BcVertexFitTree = kcvFitterBc.fit(XParticles_Bc, jpsiconstraint); 

			if (!BcVertexFitTree->isValid()) 
			{
			   h1_["BcVtxFitTreeNonValid"      ]->Fill(2);
			   continue;
			} 

			BcVertexFitTree->movePointerToTheTop();
			RefCountedKinematicParticle BcCand = BcVertexFitTree->currentParticle();
			RefCountedKinematicVertex BcVertex = BcVertexFitTree->currentDecayVertex();
			if (!BcVertex->vertexIsValid())
			{
			   h1_["BcVtxNonValid"             ]->Fill(2);
			   continue;
			}
   
			double ClBcVtx = ChiSquaredProbability((double)(BcVertex->chiSquared()),(double)(BcVertex->degreesOfFreedom()));
			h1_["Bcvtxcl"		      ]-> Fill(ClBcVtx);
			if (ClBcVtx <= pDouble_["cut_cl_Bc"])                                                           continue ;

			BcVertexFitTree->movePointerToTheFirstChild();
			RefCountedKinematicParticle muP =  BcVertexFitTree->currentParticle();
			BcVertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle muM =  BcVertexFitTree->currentParticle();
			BcVertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle pi1 =  BcVertexFitTree->currentParticle();
			BcVertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle pi2 =  BcVertexFitTree->currentParticle();
			BcVertexFitTree->movePointerToTheNextChild();
			RefCountedKinematicParticle pi3 =  BcVertexFitTree->currentParticle();

			reco::Candidate::LorentzVector MuPp4R( muP->currentState().globalMomentum().x(),muP->currentState().globalMomentum().y(), 
												   muP->currentState().globalMomentum().z(),muP->currentState().kinematicParameters().energy());
			reco::Candidate::LorentzVector MuMp4R( muM->currentState().globalMomentum().x(),muM->currentState().globalMomentum().y(), 
												   muM->currentState().globalMomentum().z(),muM->currentState().kinematicParameters().energy());
			reco::Candidate::LorentzVector Pi1p4R( pi1->currentState().globalMomentum().x(),pi1->currentState().globalMomentum().y(), 
												   pi1->currentState().globalMomentum().z(),pi1->currentState().kinematicParameters().energy());
			reco::Candidate::LorentzVector Pi2p4R( pi2->currentState().globalMomentum().x(),pi2->currentState().globalMomentum().y(), 
												   pi2->currentState().globalMomentum().z(),pi2->currentState().kinematicParameters().energy());
			reco::Candidate::LorentzVector Pi3p4R( pi3->currentState().globalMomentum().x(),pi3->currentState().globalMomentum().y(), 
												   pi3->currentState().globalMomentum().z(),pi3->currentState().kinematicParameters().energy());
			reco::Candidate::LorentzVector Bc_3Pi( BcCand->currentState().globalMomentum().x(),BcCand->currentState().globalMomentum().y(), 
												   BcCand->currentState().globalMomentum().z(),BcCand->currentState().kinematicParameters().energy());

			VertexDistance3D vertTool;
			Measurement1D BcJpsiDistance = vertTool.distance(BcVertex->vertexState(),mumuVertex.vertexState());
			double BcJpsiSignificance    = BcJpsiDistance.significance();
			double BcJpsiDistanceV       = BcJpsiDistance.value();
			h1_["BcJpsiSignificance" ]-> Fill(BcJpsiSignificance);
			h1_["BcJpsiDistance"     ]-> Fill(BcJpsiDistanceV);

			SVPos[0] 	= BcVertex->position().x();
			SVPos[1] 	= BcVertex->position().y();
			SVPos[2] 	= BcVertex->position().z(); 
			Covsec(0,0) = BcVertex->error().cxx();
			Covsec(1,0) = BcVertex->error().cyx();
			Covsec(2,0) = BcVertex->error().czx();
			Covsec(0,1) = Covsec(1,0);
			Covsec(1,1) = BcVertex->error().cyy();
			Covsec(2,1) = BcVertex->error().czy();
			Covsec(0,2) = Covsec(2,0);
			Covsec(1,2) = Covsec(2,1);
			Covsec(2,2) = BcVertex->error().czz();

		   //MC match Muons
		   if(!iEvent.isRealData())
		   {		
			  matchMu[0] = deltaR(MuPp4R.eta(),MuPp4R.phi(),p_muP.Eta(),p_muP.Phi());
			  matchMu[1] = deltaR(MuMp4R.eta(),MuMp4R.phi(),p_muM.Eta(),p_muM.Phi());

			//MC match Pions
			if(firstPionTrackCand->charge() == PICH[0])
			{
			  matchPi1[0] = deltaR(Pi1p4R.eta(),Pi1p4R.phi(),p_pi1.Eta(),p_pi1.Phi());
			  matchPi1[1] = deltaR(Pi1p4R.eta(),Pi2p4R.phi(),p_pi2.Eta(),p_pi2.Phi());
			}
			else
			{
			  matchPi1[0] = deltaR(Pi1p4R.eta(),Pi1p4R.phi(),p_pi3.Eta(),p_pi3.Phi());
			  matchPi1[1] = 10000;
			}
	  
			if(secondPionTrackCand->charge() == PICH[0])
			{
			  matchPi2[0] = deltaR(Pi2p4R.eta(),Pi2p4R.phi(),p_pi1.Eta(),p_pi1.Phi());
			  matchPi2[1] = deltaR(Pi2p4R.eta(),Pi2p4R.phi(),p_pi2.Eta(),p_pi2.Phi());
			}
			else
			{
			  matchPi2[0] = deltaR(Pi2p4R.eta(),Pi2p4R.phi(),p_pi3.Eta(),p_pi3.Phi());
			  matchPi2[1] = 10000;
			}

			if(thirdPionTrackCand->charge() == PICH[0])
			{
			  matchPi3[0] = deltaR(Pi3p4R.eta(),Pi3p4R.phi(),p_pi1.Eta(),p_pi1.Phi());
			  matchPi3[1] = deltaR(Pi3p4R.eta(),Pi3p4R.phi(),p_pi2.Eta(),p_pi2.Phi());
			}
			else
			{
			  matchPi3[0] = deltaR(Pi3p4R.eta(),Pi3p4R.phi(),p_pi3.Eta(),p_pi3.Phi());
			  matchPi3[1] = 10000;
			}
		   }//end if MC

		   //Evaluate best pointing PV
		   float pointPVCl =   0 ;
		   double Cosine 	= -10 ;
		   TransientVertex pointPVtrn;
		   reco::Vertex pointPV;
		   for (reco::VertexCollection::const_iterator pvIt = primvtx.begin(); pvIt!=primvtx.end(); pvIt++)		
		   {
			 reco::Vertex iPV = *pvIt;
			 reco::TrackCollection bcLess;
			 bcLess.reserve(tracks->size());
			 reco::Vertex::trackRef_iterator trki;
			 for (trki  = iPV.tracks_begin(); trki != iPV.tracks_end(); ++trki) 
			 {
			   reco::TrackRef trk_now(tracks, (*trki).key());
			   if(!(trk_now != innerTrkMuP && trk_now != innerTrkMuM ))													continue ;
			   if(!(trk_now != firstPionTrackCand && trk_now != secondPionTrackCand && trk_now != thirdPionTrackCand))	continue ;
			   bcLess.push_back(*trk_now);
			 }
			 h1_["PointPVDelTrks"	] -> Fill(iPV.tracksSize() - bcLess.size());

			 std::vector<reco::TransientTrack> primaryTks_tks;
			 if(bcLess.size() < 2)                            				continue;
			 for (reco::TrackCollection::const_iterator iitrk=bcLess.begin(); iitrk!=bcLess.end();iitrk++)
			 {
			  primaryTks_tks.push_back((*theB).build(*iitrk));
			 }

			 AdaptiveVertexFitter checkFitter;
			 TransientVertex checkPvs = checkFitter.vertex(primaryTks_tks);
			 if (!checkPvs.isValid()) 							  			continue;
			 float iPVCl(TMath::Prob(checkPvs.totalChiSquared(),(int)(checkPvs.degreesOfFreedom())));
			 if (iPVCl <= pDouble_["cut_cl_PV"]) 			  				continue; 
			 iPV = checkPvs;	
			 if(!(iPV.isValid())) 							  				continue;

			 double idx  = SVPos[0] - iPV.x();
			 double idy  = SVPos[1] - iPV.y();
			 double idz  = SVPos[2] - iPV.z();
			 double icos = Utils->computeCosine(idx, idy, idz, Bc_3Pi.px(), Bc_3Pi.py(), Bc_3Pi.pz());
			 if (icos > Cosine)
			 {
				pointPV 	= iPV;
				Cosine		= icos;
				pointPVCl 	= iPVCl;
				pointPVtrn  = checkPvs;
			 }
		   }	        

		   //Set vertices values
		   PointPVPos[0] = pointPV.position().x();
		   PointPVPos[1] = pointPV.position().y();
		   PointPVPos[2] = pointPV.position().z();
		   int n=0;
		   for (int i=0; i<3; i++)
		   {
			 for (int j=0; j<3; j++)
			 { 
			   SVCov[n] 		= Covsec(i,j);
			   CovPointPV(i,j) = pointPV.covariance(i,j);  
			   PointPVCov[n] 	= CovPointPV(i,j);
			   n++;
			 }
		   }

		   //Evaluate L and L/sigma and cosine for the Bc cand		
		   math::XYZVector BcPerp(Bc_3Pi.px(), Bc_3Pi.py(), 0.);
		   GlobalError BcVertexError 	= BcVertex->error();
		   GlobalPoint BcDisplacementFromBeamspot(-1*((bs.x0() -  BcVertex->position().x()) + ( BcVertex->position().z() - bs.z0()) * bs.dxdz()),
												  -1*((bs.y0() -  BcVertex->position().y()) + ( BcVertex->position().z() - bs.z0()) * bs.dydz()), 0);
	 
		   double LBcBSXY 		= BcDisplacementFromBeamspot.perp();
		   double LBcBSErr		= sqrt(BcVertexError.rerr(BcDisplacementFromBeamspot));
		   reco::Vertex::Point    vperpBc(BcDisplacementFromBeamspot.x(),BcDisplacementFromBeamspot.y(),0.);
		   double CosBcBsXY   	= vperpBc.Dot(BcPerp)/(vperpBc.R()*BcPerp.R());
		   double ElsBcBSXY	= 0 ;
		   try	    { ElsBcBSXY	= LBcBSXY / LBcBSErr ; }
		   catch   (...) {cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t" << "Floating divide by zero!" << endl ;}

		   VertexDistance3D calc3DDistance;
		   //Wrt best pointing PV
		   Measurement1D BcPointPVDistance = calc3DDistance.distance(BcVertex->vertexState(),pointPVtrn.vertexState());
		   double LBcPointPV       = BcPointPVDistance.value();
		   double LBcPointPVErr    = BcPointPVDistance.error();
		   double ElsBcPointPV     = BcPointPVDistance.significance();

		   //ntupla.fill    
		   TLorentzVector *recoBc  = new TLorentzVector(Bc_3Pi.px(),Bc_3Pi.py(),Bc_3Pi.pz(),Bc_3Pi.energy());
		   TLorentzVector *recoMuP = new TLorentzVector(MuPp4R.px(),MuPp4R.py(),MuPp4R.pz(),MuPp4R.energy());
		   TLorentzVector *recoMuM = new TLorentzVector(MuMp4R.px(),MuMp4R.py(),MuMp4R.pz(),MuMp4R.energy());
		   TLorentzVector *recoPi1 = new TLorentzVector(Pi1p4R.px(),Pi1p4R.py(),Pi1p4R.pz(),Pi1p4R.energy());
		   TLorentzVector *recoPi2 = new TLorentzVector(Pi2p4R.px(),Pi2p4R.py(),Pi2p4R.pz(),Pi2p4R.energy());
		   TLorentzVector *recoPi3 = new TLorentzVector(Pi3p4R.px(),Pi3p4R.py(),Pi3p4R.pz(),Pi3p4R.energy());
		   reco::Candidate::LorentzVector JJpsi = MuPp4R + MuMp4R;
		   TLorentzVector *recoJpsi  = new TLorentzVector(JJpsi.px(),JJpsi.py(),JJpsi.pz(),JJpsi.energy());
		   TLorentzVector *recoJpsiV = new TLorentzVector(jpsiV.px(),jpsiV.py(),jpsiV.pz(),jpsiV.energy());

		   if (jpsiV.px()!=jpsiX && jpsiV.py()!=jpsiY && jpsiV.pz()!=jpsiZ)
		   {
			 jpsiX = jpsiV.px();
			 jpsiY = jpsiV.py();
			 jpsiZ = jpsiV.pz();
			 numJpsi++;
		   }
	
		   NumBc++;

		   //Fill Ntuple
		   BcTreeCand theNewCand;

		   //Bc cand
		   theNewCand.SetBcCand 			(*recoBc															);
		   //Jpsi cand
		   theNewCand.SetJpsi 				(*recoJpsi 															);
		   theNewCand.SetJpsiV				(*recoJpsiV 														);
		   theNewCand.SetDCA  				(DCA	 															);
		   //Muons
		   theNewCand.SetMuP  				(*recoMuP 															);
		   theNewCand.SetMuM  				(*recoMuM 															);
		   theNewCand.SetMuPisGlobal   		(recomuP->isGlobalMuon() 											);
		   theNewCand.SetMuPisTracker  		(recomuP->isTrackerMuon()  											);
		   theNewCand.SetMuPisPFlow    		(recomuP->isPFMuon() 												);
		   theNewCand.SetMuPTMOST      		(muon::isGoodMuon(*recomuP, muon::TMOneStationTight) 				);
		   theNewCand.SetMuPTrkLMeas   		(recomuP->innerTrack()->hitPattern().trackerLayersWithMeasurement() );
		   theNewCand.SetMuPPixLMeas   		(recomuP->innerTrack()->hitPattern().pixelLayersWithMeasurement() 	);
		   theNewCand.SetMuPPixHits    		(recomuP->innerTrack()->hitPattern().numberOfValidPixelHits()  		);
		   theNewCand.SetMuPTrkHits    		(recomuP->innerTrack()->hitPattern().numberOfValidTrackerHits() 	);
		   theNewCand.SetMuPMatchedStat		(recomuP->numberOfMatchedStations()  								);
		   theNewCand.SetMuPNormChi2   		(recomuP->innerTrack()->normalizedChi2() 							);
		   theNewCand.SetMuPDxy   			(recomuP->innerTrack()->dxy(primvtx[0].position() )					);
		   theNewCand.SetMuPDz   			(recomuP->innerTrack()->dz(primvtx[0].position() ) 					);
		   theNewCand.SetMuMisGlobal   		(recomuM->isGlobalMuon() 											);
		   theNewCand.SetMuMisTracker  		(recomuM->isTrackerMuon()  											);
		   theNewCand.SetMuMisPFlow    		(recomuM->isPFMuon() 												);
		   theNewCand.SetMuMTMOST      		(muon::isGoodMuon(*recomuM, muon::TMOneStationTight) 				);
		   theNewCand.SetMuMTrkLMeas   		(recomuM->innerTrack()->hitPattern().trackerLayersWithMeasurement() );
		   theNewCand.SetMuMPixLMeas   		(recomuM->innerTrack()->hitPattern().pixelLayersWithMeasurement() 	);
		   theNewCand.SetMuMPixHits    		(recomuM->innerTrack()->hitPattern().numberOfValidPixelHits()  		);
		   theNewCand.SetMuMTrkHits    		(recomuM->innerTrack()->hitPattern().numberOfValidTrackerHits() 	);
		   theNewCand.SetMuMMatchedStat		(recomuM->numberOfMatchedStations()  								);
		   theNewCand.SetMuMNormChi2   		(recomuM->innerTrack()->normalizedChi2() 							);
		   theNewCand.SetMuMDxy   			(recomuM->innerTrack()->dxy(primvtx[0].position() )					);
		   theNewCand.SetMuMDz   			(recomuM->innerTrack()->dz(primvtx[0].position() ) 					);
		   // Track 1
		   theNewCand.SetPi1            	(*recoPi1															);
		   theNewCand.SetPi1Ch   			(firstPionTrackCand->charge()										);
		   theNewCand.SetTrk1PixLMeas   	(firstPionTrackCand->hitPattern().pixelLayersWithMeasurement()		);
		   theNewCand.SetTrk1TrkLMeas   	(firstPionTrackCand->hitPattern().trackerLayersWithMeasurement()	);
		   theNewCand.SetTrk1PixHits   		(firstPionTrackCand->hitPattern().numberOfValidPixelHits()  		);
		   theNewCand.SetTrk1TrkHits    	(firstPionTrackCand->hitPattern().numberOfValidTrackerHits() 		);
		   theNewCand.SetTrk1NormChi2   	(firstPionTrackCand->normalizedChi2() 								);
		   theNewCand.SetTrk1IP3DJpsi	    (Pi1_3DIp															);
		   theNewCand.SetTrk1IP3DJpsiSign 	(Pi1_3DIpSignificance												);
		   theNewCand.SetTrk1IP2DBS	    	(Pi1_IPbs 															);
		   theNewCand.SetTrk1IP2DBSSign   	(Pi1_IPbsSignificance												);
		   theNewCand.SetTrk1IP3DPV	    	(Pi1_3DIpPV 														);
		   theNewCand.SetTrk1IP3DPVSign  	(Pi1_3DIpPVSignificance												);
		   theNewCand.SetTrk1DeltaR	  		(DRPi1 																);
		   // Track 2
		   theNewCand.SetPi2            	(*recoPi2															);
		   theNewCand.SetPi2Ch   			(secondPionTrackCand->charge()										);
		   theNewCand.SetTrk2PixLMeas   	(secondPionTrackCand->hitPattern().pixelLayersWithMeasurement()		);
		   theNewCand.SetTrk2TrkLMeas   	(secondPionTrackCand->hitPattern().trackerLayersWithMeasurement()	);
		   theNewCand.SetTrk2PixHits   		(secondPionTrackCand->hitPattern().numberOfValidPixelHits()  		);
		   theNewCand.SetTrk2TrkHits    	(secondPionTrackCand->hitPattern().numberOfValidTrackerHits() 		);
		   theNewCand.SetTrk2NormChi2   	(secondPionTrackCand->normalizedChi2() 								);
		   theNewCand.SetTrk2IP3DJpsi	    (Pi2_3DIp															);
		   theNewCand.SetTrk2IP3DJpsiSign 	(Pi2_3DIpSignificance												);
		   theNewCand.SetTrk2IP2DBS	    	(Pi2_IPbs 															);
		   theNewCand.SetTrk2IP2DBSSign   	(Pi2_IPbsSignificance												);
		   theNewCand.SetTrk2IP3DPV	    	(Pi2_3DIpPV 														);
		   theNewCand.SetTrk2IP3DPVSign  	(Pi2_3DIpPVSignificance												);
		   theNewCand.SetTrk2DeltaR	  		(DRPi2 																);
		   // Track 3
		   theNewCand.SetPi3            	(*recoPi3															);
		   theNewCand.SetPi3Ch   			(thirdPionTrackCand->charge()										);
		   theNewCand.SetTrk3PixLMeas   	(thirdPionTrackCand->hitPattern().pixelLayersWithMeasurement()		);
		   theNewCand.SetTrk3TrkLMeas   	(thirdPionTrackCand->hitPattern().trackerLayersWithMeasurement()	);
		   theNewCand.SetTrk3PixHits   		(thirdPionTrackCand->hitPattern().numberOfValidPixelHits()  		);
		   theNewCand.SetTrk3TrkHits    	(thirdPionTrackCand->hitPattern().numberOfValidTrackerHits() 		);
		   theNewCand.SetTrk3NormChi2   	(thirdPionTrackCand->normalizedChi2() 								);
		   theNewCand.SetTrk3IP3DJpsi	    (Pi3_3DIp															);
		   theNewCand.SetTrk3IP3DJpsiSign 	(Pi3_3DIpSignificance												);
		   theNewCand.SetTrk3IP2DBS	    	(Pi3_IPbs 															);
		   theNewCand.SetTrk3IP2DBSSign   	(Pi3_IPbsSignificance												);
		   theNewCand.SetTrk3IP3DPV	    	(Pi3_3DIpPV 														);
		   theNewCand.SetTrk3IP3DPVSign  	(Pi3_3DIpPVSignificance												);
		   theNewCand.SetTrk3DeltaR	  		(DRPi3 																);
		   for (int i=0; i < 3; i++){
		     theNewCand.SetTrackMuPDR			(i, TrackMuPDR[i]												);
		     theNewCand.SetTrackMuMDR			(i, TrackMuMDR[i]												);
		   }
		   // Bc vertex 
		   theNewCand.SetClS		    	(ClBcVtx															);
		   theNewCand.SetEl2DWrtBS     		(LBcBSXY															);
		   theNewCand.SetEls2DWrtBS    		(ElsBcBSXY															);
		   theNewCand.SetSigma2DWrtBS  		(LBcBSErr															);
		   theNewCand.SetCos2DWrtBS    		(CosBcBsXY															);
		   theNewCand.SetEl3DWrtPV     		(LBcPointPV															);
		   theNewCand.SetEls3DWrtPV    		(ElsBcPointPV														);
		   theNewCand.SetSigma3DWrtPV  		(LBcPointPVErr														);
		   theNewCand.SetCos3DWrtPV    		(Cosine																);
		   //Vertices
		   for (int i=0; i < 3; i++){
			 theNewCand.SetBcVtxPosition 	(i, SVPos[i]														);
			 theNewCand.SetJpsiVtxPosition 	(i, JpsiVtxPos[i]													);
			 theNewCand.SetPointPVPosition 	(i, PointPVPos[i]													);
		   }
		   for (int i=0; i < 9; i++){
			 theNewCand.SetBcVtxCovariance	(i, SVCov[i]														);
			 theNewCand.SetPointPVCovariance(i, PointPVCov[i]													);
		   }
		   theNewCand.SetClJpsi	    		(ClJpsiVtx															);
		   theNewCand.SetElsigJpsi	    	(elsigJpsi															);		
		   theNewCand.SetCosJpsi	    	(cosJpsiXY															);
		   theNewCand.SetPointPVCl		    (pointPVCl		   												    );
		   // MC matching 
		   theNewCand.SetMatchMuP	    	(matchMu[0]															);
		   theNewCand.SetMatchMuM	    	(matchMu[1]															);
		   for (int i=0; i < 2; i++){
			 theNewCand.SetMatchPi1	    	(i, matchPi1[i]														);
			 theNewCand.SetMatchPi2	    	(i, matchPi2[i]														);
			 theNewCand.SetMatchPi3	    	(i, matchPi3[i]														);
		   }
		   // HLT matching 
		   for (int i=0; i < 4; i++){
			 theNewCand.SetHltMatch 		(i,HltMatch[i]														);
		   }

		 }//pi3
		 
	   }//pi2
     }//pi1 
     
   }//mu2			      
  }//mu1

  h1_["NumJPSI"      	] ->Fill(numJpsi);
  t1_["BcTree"]->Fill();

}

//=====================================================================================================================
void Bc2Jpsi3Pi::beginJob()
{
	runNumber=0;
	nEvents_ = 0 ;
 
	file    = new TFile(pString_["filename"].c_str(),"recreate");
	const bool oldAddDir = TH1::AddDirectoryStatus();
	TH1::AddDirectory(true);
	
	TDirectory * BcTreeDir      = file      ->mkdir("BcTreeDir"  ) ;
	TDirectory * JPsiDir 		= file  	->mkdir("JPsi"       ) ;
	TDirectory * MonteCarloDir  = file	    ->mkdir("MonteCarlo" ) ;
	TDirectory * VertexingDir	= file  	->mkdir("Vertexing"  ) ;
 
	file->cd() ;
	BcTreeDir->cd() ;
	t1_["BcTree"] = new TTree("BcTree","BcTree");
	BcTree_ = new BcTree();
	t1_["BcTree"] -> Branch("BcBr","BcTree",&BcTree_,64000,2);
	t1_["BcTree"] -> SetAutoSave(500000000);

    file->cd();
  
    //Histograms for JPsi
    JPsiDir->cd();
    
    h1_["ElsigJpsi"		    			] = new TH1F("ElsigJpsi",   				"L/sigma Jpsi",					 80, 0,  20 );
    h1_["JpsiPt"		    			] = new TH1F("JpsiPt",      				"JpsiPt",						100, 0, 100 );
    h1_["JpsiPtAfterTrigger"  			] = new TH1F("JpsiPtAfterTrigger", 			"JpsiPtAfterTrigger",			100, 0, 100 );
    h1_["InvMassJPsi"               	] = new TH1F("InvMassJPsi",		   			"Invariant Mass JPsi",          500, 0,  10 );
    h1_["InvMassJPsiCuts" 	    		] = new TH1F("InvMassJPsiCuts", 			"Inv Mass JPsi after trig cuts",500, 0,  10 ); 
    h1_["JpsiBeforeTrigMatch"       	] = new TH1F("JpsiBeforeTrigMatch",			"Jpsi Before TrigMatch", 		 10, 0,  10 );
    h1_["JpsiAfterTrigMatch"        	] = new TH1F("JpsiAfterTrigMatch",			"Jpsi After TrigMatch", 		 10, 0,  10 );
    h1_["InvMassJpsiPassingTrigMatch" 	] = new TH1F("InvMassJpsiPassingTrigMatch", "InvMassJpsiPassingTrigMatch",	200, 1,   5 ); 
   
    file->cd();
 
//------Histograms for MC--------------------------------------------  CHECKED
   MonteCarloDir->cd();   
   h1_["PYID"        	  	] = new TH1F("PYID",	     "PYTHIA ID",2000,-1000,1000) ;		     
   h1_["EVTA"        		] = new TH1F("EVTA",	     "Size of the generated event",1500,0,1500) ;    
   h1_["ETAMU"        		] = new TH1F("ETAMU",	     "Eta distribution of generated mumu",100,-10,10) ;    
   h1_["ETAPI1"        		] = new TH1F("ETAPI1",	     "Eta distribution of generated pi1",100,-10,10) ;    
   h1_["ETAPI2"        		] = new TH1F("ETAPI2",	     "Eta distribution of generated pi2",100,-10,10) ;    
   h1_["ETAPI3"        		] = new TH1F("ETAPI3",	     "Eta distribution of generated pi3",100,-10,10) ;    
   h1_["CHPI1"        		] = new TH1F("CHPI1",	     "charge of pi1",10,-5,5) ;    
   h1_["CHPI2"        		] = new TH1F("CHPI2",	     "charge of pi2",10,-5,5) ;    
   h1_["CHPI3"        		] = new TH1F("CHPI3",	     "charge of pi3",10,-5,5) ;    

   h1_["DAUID"       		] = new TH1F("DAUID",	     "Bc daughters ID",2000,-1000,1000) ;
   h1_["DAUJPSIID"   		] = new TH1F("DAUJPSIID",    "JPSI daughters ID",1000,-500,500) ;
   h1_["XPRIMBC"     		] = new TH1F("XPRIMBC",		 "X primary vertex Bc", 200, -1., 1.) ;
   h1_["YPRIMBC"     		] = new TH1F("YPRIMBC",		 "Y primary vertex Bc", 200, -1., 1.) ;
   h1_["ZPRIMBC"     		] = new TH1F("ZPRIMBC",		 "Z primary vertex Bc", 300, -15.0,  15.0)  ;
   h1_["XSECBC"      		] = new TH1F("XSECBC", 		 "X secondary vertex Bc", 200, -1., 1.) ;
   h1_["YSECBC"      		] = new TH1F("YSECBC", 		 "Y secondary vertex Bc", 200, -1., 1.) ;
   h1_["ZSECBC"      		] = new TH1F("ZSECBC", 		 "Z secondary vertex Bc", 300, -15.0, 15.0) ;   
   h1_["ELLEBC"      		] = new TH1F("ELLEBC", 		 "Bc decay lenght",	 400,  -2.,   2.0)  ;
   h1_["TAUBC"			    ] = new TH1F("TAUBC", 		 "Bc decay time", 2000,  -10,   10)  ;
   h1_["PTOTMUP"     		] = new TH1F("PTOTMUP",      "Momentum of 1st muon from Jpsi", 100, 0., 100. );
   h1_["PTOTMUM"     		] = new TH1F("PTOTMUM",      "Momentum of 2nd muon from Jpsi", 100, 0., 100. );
   h1_["PTOTPI1"      		] = new TH1F("PTOTPI1",      "Momentum of Pi1 from Bc+", 100, 0., 100. );
   h1_["PTOTPI2"      		] = new TH1F("PTOTPI2",      "Momentum of Pi2 from Bc+", 100, 0., 100. );
   h1_["PTOTPI3"      		] = new TH1F("PTOTPI3",      "Momentum of Pi3 from Bc+", 100, 0., 100. );
   h1_["PTJPSI"      		] = new TH1F("PTJPSI",       "Transverse Momentum of Jpsi", 400, 0., 40. );
   h1_["PTMUP"	     		] = new TH1F("PTMUP",	     "Transverse momentum of 1st muon from Jpsi", 250, 0., 25. ); 
   h1_["PTMUM"	     		] = new TH1F("PTMUM",	     "Transverse momentum of 2nd muon from Jpsi", 250, 0., 25. );
   h1_["PTPI1"	     		] = new TH1F("PTPI1",	     "Transverse momentum of Pi1 from Bc+", 250, 0., 25. );
   h1_["PTPI2"	     		] = new TH1F("PTPI2",	     "Transverse momentum of Pi2 from Bc+", 250, 0., 25. );
   h1_["PTPI3"	     		] = new TH1F("PTPI3",	     "Transverse momentum of Pi3 from Bc+", 250, 0., 25. );
   h1_["BCPT"               ] = new TH1F("BCPT",          "Transverse momentum of Bc"    , 400, 0., 100. );
   h1_["COSINE"   			] = new TH1F("COSINE",	     "Cosine at gen level",  240, -1.2, 1.2) ;
   h2_["PMUPMUM"  		    ] = new TH2F("PMUPMUM",      "MuM momentum vs MuP momentum", 40, 0., 40., 40, 0., 40.);
   h2_["PTMUPMUM"    		] = new TH2F("PTMUPMUM", 	 "MuM pt vs MuP pt", 40, 0., 20., 40, 0., 20.);
   file->cd();				


//-----Histograms for Vertex Stuff-------------------------------------------- CHECKED
   VertexingDir->cd();   

   h1_["pVertexSize"					    ] = new TH1F("pVertexSize", 					  "pVertexSize after trigger requirement ",		 20,0	  ,20	);
   h1_["pVertexTrackSize"				    ] = new TH1F("pVertexTrackSize",				  "pVertexTrackSize", 	 800,0     ,800   );
   h1_["EventTrackSize" 				    ] = new TH1F("EventTrackSize",  				  "EventTrackSize",		 800,0     ,800   );
   h1_["goodTrkSize"        			    ] = new TH1F("goodTrkSize", 	 				  "number of tracks passing selections", 				800 , 0,800 );
   h1_["Bcvtxcl"		    			    ] = new TH1F("Bcvtxcl", 	 					  "Bc vertex CL before cut", 	 10000,0.	, 1.  );
   h1_["Pion1_ImpactParameter"		   		] = new TH1F("Pion1_ImpactParameter",			  "Pion1_ImpactParameter",   10000,0   ,50   );
   h1_["Pion1_IPsignificance"              	] = new TH1F("Pion1_IPsignificance" ,			  "Pion1_IPsignificance" ,   10000,0   ,100  );
   h1_["Pion1_TransverseImpactParameter_BS"	] = new TH1F("Pion1_TransverseImpactParameter_BS","Pion1_TransverseImpactParameter_BS",  5000, 0,10   );
   h1_["Pion1_TIParameterBSsignificance"   	] = new TH1F("Pion1_TIParameterBSsignificance" ,  "Pion1_TIParameterBSsignificance"     , 10000, 0,50   );
   h1_["Pion2_ImpactParameter"		   		] = new TH1F("Pion2_ImpactParameter",			  "Pion2_ImpactParameter",   10000,0   ,50    );
   h1_["Pion2_IPsignificance"              	] = new TH1F("Pion2_IPsignificance" ,			  "Pion2_IPsignificance" ,   10000,0   ,100   );
   h1_["Pion2_TransverseImpactParameter_BS"	] = new TH1F("Pion2_TransverseImpactParameter_BS","Pion2_TransverseImpactParameter_BS",  5000, 0, 10  );
   h1_["Pion2_TIParameterBSsignificance"   	] = new TH1F("Pion2_TIParameterBSsignificance" ,  "Pion2_TIParameterBSsignificance"     , 10000, 0, 50  );
   h1_["Pion3_ImpactParameter"		   		] = new TH1F("Pion3_ImpactParameter",			  "Pion3_ImpactParameter",  10000,0   ,50   );
   h1_["Pion3_IPsignificance"              	] = new TH1F("Pion3_IPsignificance" ,			  "Pion3_IPsignificance" ,  10000,0   ,100  );
   h1_["Pion3_TransverseImpactParameter_BS"	] = new TH1F("Pion3_TransverseImpactParameter_BS","Pion3_TransverseImpactParameter_BS",  5000, 0 ,10 );
   h1_["Pion3_TIParameterBSsignificance"   	] = new TH1F("Pion3_TIParameterBSsignificance" ,  "Pion3_TIParameterBSsignificance" ,     10000, 0, 50  );
   h1_["BcVtxFitTreeNonValid"   		    ] = new TH1F("BcVtxFitTreeNonValid",			  "BcVtxFitTreeNonValid", 			4, 0,  4 );
   h1_["BcVtxNonValid"		    		    ] = new TH1F("BcVtxNonValid",					  "BcVtxNonValid",				4, 0,  4 );
   h1_["BcJpsiSignificance"	    		    ] = new TH1F("BcJpsiSignificance",  			  "BcJpsiSignificance",  	 100, 0.  ,50.0 ) ;
   h1_["BcJpsiDistance"	        		    ] = new TH1F("BcJpsiDistance",  	 			  "BcJpsiDistance",  		 10000, 0.  ,5.0 ) ;
   h1_["PointPVDelTrks"	        			] = new TH1F("PointPVDelTrks",					  "Number of deleted trks from PointPV"   ,       5,   0,   5.0 );
   h1_["deltaRPi1"	            			] = new TH1F("deltaRPi1",	   		 			  "deltaRPi1",		   1000, 0.  ,10.0  ) ;
   h1_["deltaRPi2"	            			] = new TH1F("deltaRPi2",	   					  "deltaRPi2",		   1000, 0.  ,10.0  ) ;
   h1_["deltaRPi3"	            			] = new TH1F("deltaRPi3",	   		 			  "deltaRPi3",		   1000, 0.  ,10.0  ) ;

   file->cd();
   h1_["filter"      		      	] = new TH1F("filter",	   "Binned filter counter",100,0,100); // muon size
   h1_["TrueNumInteraction"      	] = new TH1F("TrueNumInteraction", "TrueNumInteraction",50,0,50); // muon size
   h1_["NumJPSI"      				] = new TH1F("NumJPSI",			               "NumJPSI", 50,0, 50); 

   TH1::AddDirectory(oldAddDir);
}
//============================================================================================================
// ------------ method called once each job just after ending the event loop  ------------
void Bc2Jpsi3Pi::endJob() {
  file->Write();
  file->Close();
  cout << endl
       << endl
       << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t"
       << "======================= Saving histograms to " 
       << pString_["filename"]
       << " ======================="
       << endl 
       << endl ;
}
//===============================================================================================
void Bc2Jpsi3Pi::MonteCarloStudies(const edm::Event& iEvent)
{
   //HepMC study
   	if( iEvent.isRealData() ) {return ;}
   	int id;
   	int JpsiId      =  443;
   	int BcplusId    =  541;
    int piplus      =  211;
    int piminus     = -211;
    double BcMassMC = 6.275;

    boolMC        = false;
	bool boolJpsi = false;
	bool boolpi1  = false;
	bool boolpi2  = false;
	bool boolpi3  = false;
	
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
		
	int quanteBc=0;
   	ELLE = -1.;
    for ( size_t i=0; i< genParticles->size(); ++i) 
    { 
	  const reco::GenParticle &p = (*genParticles)[i];
	  id = p.pdgId();
	  h1_["PYID"]->Fill(id);
	  if (fabs(id) != BcplusId ) continue;
	  quanteBc++;
	  if (quanteBc!=2) continue;

	  // return all daughters of Bc      	
	  int quantipiplus = 0;    
	  int quantipiminus = 0;	       
	  for ( size_t ides=0; ides < p.numberOfDaughters(); ++ides ) 
	  {
		const reco::Candidate *des = p.daughter(ides);
		int dauId = des->pdgId();
		h1_["DAUID"]->Fill(dauId);
		if( dauId == JpsiId ) 
		{
		  BCVTX[0] = des->vertex().x();
		  BCVTX[1] = des->vertex().y();
		  BCVTX[2] = des->vertex().z();
		  boolJpsi = true;

		// return all daughters of J/psi    
		  if (des->numberOfDaughters()!=2)  		continue;  
		  for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
		  {
			const reco::Candidate *mumu = des->daughter(imumu);
			if	( mumu->pdgId() == -13 ) 	
			{
				p_muP.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
			    JPSIVTX[0] =  mumu->vertex().x();
			    JPSIVTX[1] =  mumu->vertex().y();
			    JPSIVTX[2] =  mumu->vertex().z();
			} 
			else if ( mumu->pdgId() == 13 )
			{
				p_muM.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
			}
			else 									continue;		 
			h1_["DAUJPSIID"] -> Fill(mumu->pdgId());    		
			h1_["ETAMU"    ] -> Fill(mumu->momentum().Eta());
		  }
	    } //end if dauId=jpsiId
	      
		if( dauId == piplus && quantipiplus == 0) 
		{
		  quantipiplus++;
		  p_pi1.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
		  PICH[0]=des->charge();
		  boolpi1 = true;
		  h1_["ETAPI1"]->Fill(p_pi1.Eta());
		  h1_["CHPI1" ]->Fill(des->charge());
		}	      
		else if( dauId == piplus && quantipiplus == 1 ) 
		{
		  quantipiplus++; 
		  p_pi2.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
		  PICH[1]=des->charge();
		  boolpi2 = true;
		  h1_["ETAPI2"]->Fill(p_pi2.Eta());
		  h1_["CHPI2" ]->Fill(des->charge());
		}	      
	    if( dauId == piminus && quantipiminus == 0) 	
	    {
		  quantipiminus++;
		  p_pi3.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
		  PICH[2] = des->charge();
		  boolpi3    = true;
	      h1_["ETAPI3"]->Fill(p_pi3.Eta());
		  h1_["CHPI3" ]->Fill(des->charge());
	   } 	           
     } // end for des
      
	  PVVTX[0] = p.vertex().x();
	  PVVTX[1] = p.vertex().y();
	  PVVTX[2] = p.vertex().z();
	  
	  double dx = (BCVTX[0] - PVVTX[0]);
	  double dy = (BCVTX[1] - PVVTX[1]);
	  double dz = (BCVTX[2] - PVVTX[2]);
	  ELLE = sqrt( dx*dx + dy*dy + dz*dz );      
	  p_Bc.SetPxPyPzE(p.px(), p.py(), p.pz(), p.energy());
	  TAU=ELLE*BcMassMC/c_const/sqrt((p.px()*p.px())+ (p.py()*p.py())+ (p.pz()*p.pz()));
      COSINE = Utils->computeCosine(dx, dy, dz, p.px(), p.py(), p.pz());

      if(boolJpsi && boolpi1 && boolpi2 && boolpi3) boolMC = true;
    } // end genParticle loop

    double ptotmup = sqrt(p_muP.Px()*p_muP.Px() + p_muP.Py()*p_muP.Py() + p_muP.Pz()*p_muP.Pz());
    double ptotmum = sqrt(p_muM.Px()*p_muM.Px() + p_muM.Py()*p_muM.Py() + p_muM.Pz()*p_muM.Pz());
    double ptotpi1 = sqrt(p_pi1.Px()*p_pi1.Px()+p_pi1.Py()*p_pi1.Py()+p_pi1.Pz()*p_pi1.Pz());
    double ptotpi2 = sqrt(p_pi2.Px()*p_pi2.Px()+p_pi2.Py()*p_pi2.Py()+p_pi2.Pz()*p_pi2.Pz());
    double ptotpi3 = sqrt(p_pi3.Px()*p_pi3.Px()+p_pi3.Py()*p_pi3.Py()+p_pi3.Pz()*p_pi3.Pz());
    p_jpsi = p_muP+p_muM;

    h1_["PTOTMUP"	] -> Fill(ptotmup);
    h1_["PTOTMUM"	] -> Fill(ptotmum);
    h1_["PTOTPI1"	] -> Fill(ptotpi1);
    h1_["PTOTPI2"	] -> Fill(ptotpi2);
    h1_["PTOTPI3"	] -> Fill(ptotpi3);    
    h1_["PTMUP"  	] -> Fill(p_muP.Perp());
    h1_["PTMUM"  	] -> Fill(p_muM.Perp());
    h1_["PTPI1"		] -> Fill(p_pi1.Perp());
    h1_["PTPI2"		] -> Fill(p_pi2.Perp());
    h1_["PTPI3"		] -> Fill(p_pi3.Perp());
    h1_["PTJPSI"    ] -> Fill(p_jpsi.Perp());		
    h1_["TAUBC"		] -> Fill(TAU);
    h1_["COSINE"	] -> Fill(COSINE);			    
	h1_["XPRIMBC"   ] -> Fill(PVVTX[0]);
	h1_["YPRIMBC"   ] -> Fill(PVVTX[1]);
	h1_["ZPRIMBC"   ] -> Fill(PVVTX[2]);
	h1_["XSECBC"    ] -> Fill(BCVTX[0]);
	h1_["YSECBC"    ] -> Fill(BCVTX[1]);
	h1_["ZSECBC"    ] -> Fill(BCVTX[2]);
	h1_["ELLEBC"    ] -> Fill(ELLE);
    h2_["PMUPMUM"   ] -> Fill(ptotmup,ptotmum);
    h2_["PTMUPMUM"	] -> Fill(p_muP.Perp(),p_muM.Perp());			    

}

//============================================================================================== 
//Acquire Input Parameters

void Bc2Jpsi3Pi::acquireInputParameters(void)
{
   std::cout << "\n\n"
             << __LINE__ << "]\t" 
             << __PRETTY_FUNCTION__  
	     << "\t=================== Parameters read from python configuration file ============================="
	     << endl << endl ;
   
   map<string, vector<string> > containers ;
   
   containers["double"] = inputDouble_ ;
   containers["string"] = inputString_ ;
   for(map<string, vector<string> >::iterator typeIt =containers.begin();
                                              typeIt!=containers.end()  ;
					      ++typeIt)
   {
     for(std::vector<std::string>::iterator it =typeIt->second.begin();
     					    it!=typeIt->second.end();
     					    ++it)
     {
       static const boost::regex exp("(.+)?=(.+)", boost::regex::perl);
       if( typeIt->first == "double" )
         {
           boost::cmatch what ;
           if(boost::regex_match((*it).c_str(), what, exp))
     	     { 
     	       pDouble_[what[1]] = atof(what[2].first) ;
     	       std::cout << what[1]
     	   		 << "\t= "
     	   		 << what[2].first
     	   		 << endl ;//stampa a schermo 
     	     }
           else
     	     {
     	       std::cout << __LINE__ << "]\t" << __PRETTY_FUNCTION__  
     	   		 << "\tInvalid parameter syntax for " 
     	   		 << *it
     	   		 << endl ;
     	       assert(0) ;
     	     }
          }
	else
	  {
     	    boost::smatch what ;
     	    static const boost::regex exp("(.+)?=(.+)");
     	    if(boost::regex_match(*it, what, exp, boost::match_extra)) 
     	      { 
     	    	pString_[what[1]] = what[2] ;
     	    	std::cout << what[1]
	    		  << "\t= "
	    		  << what[2]
	    		  << endl ;
     	      }
     	    else
     	      {
     	    	std::cout << __LINE__ << "]\t" << __PRETTY_FUNCTION__  
	    		  << "\tInvalid parameter syntax for " 
	    		  << *it
	    		  << endl ;
	    	assert(0) ;
     	      }
	  }
     }
   }
   
   std::cout << "\n\n"
             << __LINE__ << "]\t" 
             << __PRETTY_FUNCTION__  
	     << "\t================================================================================================"
	     << endl << endl ;
}
//fine acquisizione parametri input
//========================================================================================
void Bc2Jpsi3Pi::checkHLT(const edm::Event& iEvent)
{
   std::string theName   = pString_["HLTMatchName"];
   std::string trigpart1 = pString_["HLTMatchModule1"];   
   std::string trigpart2 = pString_["HLTMatchModule2"];

   edm::Handle<trigger::TriggerEvent> theTriggerEvent;
   iEvent.getByLabel(theName,theTriggerEvent);
   std::vector<reco::Particle::LorentzVector>       JpsiSandra;
   std::vector<std::pair< double,double> >          MuonTrig;
   JpsiSandra.clear();
   MuonTrig.clear();
   JPsiTriggered_.clear() ;
 
   const trigger::TriggerObjectCollection& theTriggerObjects(theTriggerEvent->getObjects());

   int count = 0 ;
   for(size_t trigFilter = 0; trigFilter < theTriggerEvent->sizeFilters(); trigFilter++)
   {
     std::string triggerName = theTriggerEvent->filterTag(trigFilter).encode();
     if(triggerName.find(trigpart1) != std::string::npos ) count++ ; 
     if(triggerName.find(trigpart2) != std::string::npos ) count++ ;
   }
   
   for(size_t trigFilter = 0; trigFilter < theTriggerEvent->sizeFilters(); trigFilter++)
   {
     std::string triggerName = theTriggerEvent->filterTag(trigFilter).encode();
     if( count==2 &&  (triggerName.find(trigpart1) != std::string::npos ||
                       triggerName.find(trigpart2) != std::string::npos))
     {
	   const trigger::Keys &keys = theTriggerEvent->filterKeys(trigFilter);
	   reco::Particle::LorentzVector JpsiSandra_temp;
	   for(trigger::Keys::const_iterator kIt = keys.begin(); kIt != keys.end(); ++kIt)
	   {
		 MuonTrig.push_back(std::make_pair(theTriggerObjects[*kIt].particle().eta(), theTriggerObjects[*kIt].particle().phi())) ;
	   }
	   if(MuonTrig.size()>0)
	   {
         JPsiTriggered_.push_back(std::make_pair(
	                          std::make_pair(MuonTrig[0].first,MuonTrig[0].second),
	                          std::make_pair(MuonTrig[1].first,MuonTrig[1].second))) ;
       }
     }
   }
}
//===============================================================================================================
//define this as a plug-in
DEFINE_FWK_MODULE(Bc2Jpsi3Pi);
