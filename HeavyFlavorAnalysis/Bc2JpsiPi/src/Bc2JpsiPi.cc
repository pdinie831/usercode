// -*- C++ -*-
// Package:    Bc2JpsiPi
// Class:      Bc2JpsiPi
// 
/**\class  Bc2JpsiPi Bc2JpsiPi.cc HeavyFlavor/Bc2JpsiPi/src/Bc2JpsiPi.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Silvia Taroni & Sandra Malvezzi
// Put in local cvs
// Added in cvs milano
//         Created:  Thu Apr 23 16:09:09 CEST 2009


// system include files
#include <memory>

// user include files
#include <DataFormats/Math/interface/deltaR.h>

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"   
#include "HepMC/GenParticle.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "CLHEP/Vector/LorentzVector.h"

#include <boost/regex.hpp> 
#include <typeinfo>

#include "HeavyFlavorAnalysis/Bc2JpsiPi/interface/Bc2JpsiPi.h"

// constants, enums and typedefs
   double c_const  				= 0.0299792458;//cm/ps
   const ParticleMass muon_mass = 0.10565837; 
   ParticleMass jpsimass_c 	  	= 3.096916; 
   ParticleMass pion_mass  	  	= 0.13957018; 
   ParticleMass kaon_mass      	= 0.493677;
   int JpsiId      				= 443; 

// static data member definitions
   TLorentzVector p_muP,p_muM, p_pi, p_jpsi,p_Bc;

   double ELLE     = 0.;
   double TAU        = 0.;
   double COSINE     = 0.;
   double PVVTX[3]   = {0,0,0};
   double BCVTX[3]   = {0,0,0};
   double JPSIVTX[3] = {0,0,0};
   double BSPOT[3]   = {0,0,0};
   int    PICH       = 0;
   bool   boolMC     = false;

// constructors and destructor
using namespace std ;

Bc2JpsiPi::Bc2JpsiPi(const edm::ParameterSet& iConfig) :
   boolTrg_    (iConfig.getParameter<bool>                    ("printTriggers")),
   inputDouble_(iConfig.getParameter<std::vector<std::string> >("inputDouble_")),
   inputString_(iConfig.getParameter<std::vector<std::string> >("inputString_")),
   thePVs_(     iConfig.getParameter<edm::InputTag>	       ("primaryVertexTag")),
   thebeamspot_(iConfig.getParameter<edm::InputTag>	       ("beamSpotTag"     ))
{
  cout << __LINE__ << "]\t" << __PRETTY_FUNCTION__ << "\tCtor" << endl ;
  thisConfig_ = iConfig;
  acquireInputParameters() ;
  Utils = new UsefulTools();
}	

//==============================================================================================
Bc2JpsiPi::~Bc2JpsiPi()
{
  delete Utils;
}

//=============== method called to for each event ===================================================================
void Bc2JpsiPi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  BcTree_    ->BcTree::BcTreeClear();
  runNumber++;
  nEvents_++ ;

  // Handles to collections 
  edm::Handle<reco::BeamSpot>                    theBeamSpot	;
  iEvent.getByLabel(thebeamspot_,                theBeamSpot)	;
  reco::BeamSpot bs = *theBeamSpot;

  edm::Handle<reco::MuonCollection>   			 recoMuons      ; 
  iEvent.getByLabel("muons",          			 recoMuons)		;

  edm::Handle<reco::TrackCollection>  			 muonTracks     ;
  iEvent.getByLabel("globalMuons",    			 muonTracks)	;

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
									   p_pi 	,
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
  double 	matchPi 	= 10000		;
 	
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
		(triggerNames_.triggerName(itrig) == pString_["HLTname3"]) 
	 ) 
	{
	  goodTrig.push_back(itrig) ;
    }
  }
  if ( goodTrig.size() == 0 ) {
	cout << __LINE__ << "\tthe required trigger paths are not present is this dataset!" << endl;
	std::cout << "Looking for " << pString_["HLTname1"] << "," << pString_["HLTname2"] << "," << pString_["HLTname3"]  << std::endl;
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
	 }														        
	 if (flag)												   	continue;	        
																   
	 if (checkTrk->pt()  				  					< pDouble_["cut_Pt_Trk"])			   			continue;	        
	 if (fabs(checkTrk->eta())			  					> pDouble_["cut_eta_pi"])			  			continue;	   						  
// 	 if (checkTrk->numberOfValidHits()	  					< pDouble_["cut_nhits"] )			   			continue;	   					   
// 	 if (checkTrk->hitPattern().numberOfValidPixelHits() 	< pDouble_["numberOfValidPixelHitsCut"]  )	    continue;	        
// 	 if (checkTrk->normalizedChi2()			 				> pDouble_["cut_chi2n"])			   		    continue; 	        
	 if (!iEvent.isRealData() && checkTrk->charge() !=1)  													continue;   //per Bc
	 qualityTracks[tracksIt]=checkTrk;											   	
  }
 
  //set BcTreeHeader   BcNum=0
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

  h1_["filter"]->Fill(counter); //0
  if ( !triggerAccepted ) {
    if (boolMC) t1_["BcTree"]->Fill();
    return ;
  }
  this->checkHLT(iEvent) ;  //to be checked !!!!!!!!!!!!!!!!

  h1_["filter"      ]->Fill(counter++);//1
  if (recoMuons->size()<2) {
	if (boolMC) t1_["BcTree"]->Fill();
	return ;
  }
    
  h1_["pVertexSize"         ]-> Fill( primvtx.size()		 );
  h1_["pVertexTrackSize"	]-> Fill( primvtx[0].tracksSize());
  h1_["EventTrackSize"		]-> Fill( tracks->size() 		 );
  h1_["goodTrkSize"			]-> Fill( qualityTracks.size() 	 );

  if(qualityTracks.size() < 1 )
  { 
	std::cout << "exit for track size < 1 " << std::endl; 
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

	//POG SOFT Muon Selection 2011 - da togliere!
// 	if (!muon::isGoodMuon(*recomu1,muon::TMOneStationTight)                    )				 	continue;
// 	if (recomu1->innerTrack()->hitPattern().numberOfValidTrackerHits()   <= 10 )					continue;
// 	if (recomu1->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 1  )					continue;
// 	if (recomu1->innerTrack()->normalizedChi2() 						 	>= 1.8)	 				continue;
// 	if (!(fabs(recomu1->innerTrack()->dxy(primvtx[0].position())) < 3. && 
// 		  fabs(recomu1->innerTrack()->dz(primvtx[0].position()))  < 30. )) 							continue;
	
    //Loop on the second muon 
    for(reco::MuonCollection::const_iterator muR2 =muR1+1; muR2!=recoMuons->end(); ++muR2)
    {
	  recomu2 = &*muR2;

	  if (muR1->charge() + muR2->charge() != 0)								   						continue; 
	  if (fabs(muR2->eta()) > pDouble_["cut_eta"])								    				continue;
	  if (fabs(muR2->pt() ) < pDouble_["cut_Pt_Mu"])												continue;
	
      //POG SOFT Muon Selection - da togliere!!
// 	  if (!muon::isGoodMuon(*recomu2,muon::TMOneStationTight)                    )				 	continue;
// 	  if (recomu2->innerTrack()->hitPattern().numberOfValidTrackerHits()   <= 10 )					continue;
// 	  if (recomu2->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 1  )					continue;
// 	  if (recomu2->innerTrack()->normalizedChi2() 						   >= 1.8)	 				continue;
// 	  if (!(fabs(recomu2->innerTrack()->dxy(primvtx[0].position())) < 3. && 
// 			fabs(recomu2->innerTrack()->dz(primvtx[0].position()))  < 30. )) 						continue;
 
	  if( !recomu1->track().isNonnull() || !recomu2->track().isNonnull())
	  {
		assert(0) ;
	  }

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

      //Trigger requirements: CL, L/Sigma, Cosine and pT
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

// 	  TLorentzVector *jpsiMuP, *jpsiMuM;
// 	  jpsiMuP  = new TLorentzVector(MuPp4.px(), MuPp4.py(), MuPp4.pz(), MuPp4.energy());
// 	  jpsiMuM  = new TLorentzVector(MuMp4.px(), MuMp4.py(), MuMp4.pz(), MuMp4.energy());
	  TLorentzVector *jpsiTV   = new TLorentzVector(jpsiV.px(), jpsiV.py(), jpsiV.pz(), jpsiV.energy());
	  JpsiVtxPos[0] = mumuVertex.position().x() ;
	  JpsiVtxPos[1] = mumuVertex.position().y() ;
	  JpsiVtxPos[2] = mumuVertex.position().z() ;

      //MC match Muons
// 	  double matchMuJpsi[2];
// 	  if(!iEvent.isRealData())
// 	  {
// 		  matchMuJpsi[0] = deltaR(jpsiMuP->Eta(),jpsiMuP->Phi(),p_muP.Eta(),p_muP.Phi());
// 		  matchMuJpsi[1] = deltaR(jpsiMuM->Eta(),jpsiMuM->Phi(),p_muM.Eta(),p_muM.Phi());
// 	  }

      //set JpsiCand
// 	  BcTree_->BcTree::AddJpsiCand(  
// 					  *jpsiMuP	    ,
// 					  *jpsiMuM	    ,
// 					  *jpsiTV		,
// 					  elsigJpsi     ,
// 					  cosJpsiXY     ,
// 					  ClJpsiVtx     ,
// 					  JpsiVtx       ,
// 					  matchMuJpsi   
// 					  );

//----------------------------------------------------------------------------
      //Loop on the tracks that have passed pre-selection cuts
// 	  unsigned int trkI = -1 ;
	  for (newTracksDef::const_iterator it_pi=qualityTracks.begin(); it_pi!=qualityTracks.end(); ++it_pi)
	  {
		// trkI=(*it_pi).first ;
		reco::TrackRef PionTrackCand((*it_pi).second) ;
		const reco::TransientTrack piTrackCand = (*theB).build(PionTrackCand);
		
		pair<double,double>  pion_IPPair = Utils->pionImpactParameter(piTrackCand,mumuVertex);
		double pi_3DIp                   = pion_IPPair.first;
		double pi_3DIpSignificance       = pion_IPPair.second;
		h1_["Pion_ImpactParameter"]->Fill(pi_3DIp);
		h1_["Pion_IPsignificance" ]->Fill(pi_3DIpSignificance);
		
		pair<double,double>  Pion_BSPair = Utils->pionIPBeamSpot(piTrackCand,BeamSpotGP);
		double pi_IPbs                   = Pion_BSPair.first;
		double pi_IPbsSignificance       = Pion_BSPair.second;
		h1_["Pion_TransverseImpactParameter_BS" ]-> Fill(pi_IPbs);
		h1_["Pion_TIParameterBSsignificance"    ]-> Fill(pi_IPbsSignificance);

		pair<double,double>  pion_IPPairPV = Utils->pionImpactParameter(piTrackCand,primvtx[0]);
		double pi_3DIpPV                   = pion_IPPairPV.first;
		double pi_3DIpPVSignificance       = pion_IPPairPV.second;

		double DRPi =  deltaR(PionTrackCand->eta(),PionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
		h1_["deltaRPi"     ]-> Fill(DRPi);
        if (DRPi > 5) 								continue;

        //Kinematic vertex fit with Jpsi constraint for Bc->JpsiPi	
		KinematicParticleFactoryFromTransientTrack pFactory;
		float muon_sigma 			  = muon_mass*1.e-6;
		float pion_sigma        	  = pion_mass*1.e-6; 
		float kaon_sigma              = kaon_mass*1.e-6; 
		float chi 		 			  = 0.;
		float ndf 		 			  = 0.;

		std::vector<RefCountedKinematicParticle> XParticles_Bc;
		XParticles_Bc.push_back(pFactory.particle(tTrackmuP,  muon_mass,chi,ndf,muon_sigma));
		XParticles_Bc.push_back(pFactory.particle(tTrackmuM,  muon_mass,chi,ndf,muon_sigma));
		XParticles_Bc.push_back(pFactory.particle(piTrackCand,pion_mass,chi,ndf,pion_sigma)); 

		MultiTrackKinematicConstraint * jpsiconstraint = new TwoTrackMassKinematicConstraint(jpsimass_c);
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
		h1_["Bcvtxcl"]-> Fill(ClBcVtx);
		if (ClBcVtx <= pDouble_["cut_cl_Bc"])                                       continue ;

		BcVertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle muP =  BcVertexFitTree->currentParticle();
		BcVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle muM =  BcVertexFitTree->currentParticle();
		BcVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle pi =  BcVertexFitTree->currentParticle();
		BcVertexFitTree->movePointerToTheNextChild();

		reco::Candidate::LorentzVector MuPp4R( muP->currentState().globalMomentum().x(),muP->currentState().globalMomentum().y(), 
											   muP->currentState().globalMomentum().z(),muP->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector MuMp4R( muM->currentState().globalMomentum().x(),muM->currentState().globalMomentum().y(), 
											   muM->currentState().globalMomentum().z(),muM->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector Pip4R(  pi->currentState().globalMomentum().x(),pi->currentState().globalMomentum().y(), 
											   pi->currentState().globalMomentum().z(),pi->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector Bc_Pi(  BcCand->currentState().globalMomentum().x(),BcCand->currentState().globalMomentum().y(), 
											   BcCand->currentState().globalMomentum().z(),BcCand->currentState().kinematicParameters().energy());

        //Kinematic vertex fit with Jpsi constrain for B+->JpsiK
		chi=0;
		ndf=0;
		pion_sigma = pion_mass*1.e-6; 
		kaon_sigma = kaon_mass*1.e-6; 

		std::vector<RefCountedKinematicParticle> XParticles_Bp;
		XParticles_Bp.push_back(pFactory.particle(tTrackmuP,muon_mass,chi,ndf,muon_sigma));
		XParticles_Bp.push_back(pFactory.particle(tTrackmuM,muon_mass,chi,ndf,muon_sigma));
		XParticles_Bp.push_back(pFactory.particle(piTrackCand,kaon_mass,chi,ndf,kaon_sigma)); 
		KinematicConstrainedVertexFitter kcvFitterBp;
		RefCountedKinematicTree BpVertexFitTree = kcvFitterBp.fit(XParticles_Bp, jpsiconstraint); 

		if (!BpVertexFitTree->isValid())  		continue;
		
		BpVertexFitTree->movePointerToTheTop();
		RefCountedKinematicParticle BpCand = BpVertexFitTree->currentParticle();
		RefCountedKinematicVertex BpVertex = BpVertexFitTree->currentDecayVertex();
		if (!BpVertex->vertexIsValid())			continue;
		
		double ClBpVtx=ChiSquaredProbability((double)(BpVertex->chiSquared()),(double)(BpVertex->degreesOfFreedom()));
		h1_["Bpvtxcl"  ]->Fill(ClBpVtx);
		if ( ClBpVtx<pDouble_["cut_cl_Bc"] ) 	continue;
	   
		BpVertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle muP_2 =  BpVertexFitTree->currentParticle();
		BpVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle muM_2 = BpVertexFitTree->currentParticle();
		BpVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle kaon1 = BpVertexFitTree->currentParticle();

		reco::Candidate::LorentzVector kp4R( kaon1->currentState().globalMomentum().x(),kaon1->currentState().globalMomentum().y(), 
											 kaon1->currentState().globalMomentum().z(),kaon1->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector Bc_K( BpCand->currentState().globalMomentum().x(),BpCand->currentState().globalMomentum().y(), 
											 BpCand->currentState().globalMomentum().z(),BpCand->currentState().kinematicParameters().energy());

		VertexDistance3D vertTool;
		Measurement1D BcJpsiDistance = vertTool.distance(BcVertex->vertexState(), mumuVertex.vertexState());
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

        //MC match 
		if(!iEvent.isRealData())
		{
		  matchMu[0] = deltaR(MuPp4R.eta(),MuPp4R.phi(),p_muP.Eta(),p_muP.Phi());
		  matchMu[1] = deltaR(MuMp4R.eta(),MuMp4R.phi(),p_muM.Eta(),p_muM.Phi());
		  if(PionTrackCand->charge() == PICH) 	matchPi = deltaR(Pip4R.eta(), Pip4R.phi(), p_pi.Eta(), p_pi.Phi()); 
		}

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
			if(!(trk_now != innerTrkMuP && trk_now != innerTrkMuM && trk_now != PionTrackCand))		  	continue ;
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
		  double icos = Utils->computeCosine(idx, idy, idz, Bc_Pi.px(), Bc_Pi.py(), Bc_Pi.pz());
		  if (icos > Cosine)
		  {
			 pointPV 	= iPV;
			 Cosine		= icos;
			 pointPVCl 	= iPVCl;
			 pointPVtrn = checkPvs;
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

		//Evaluate L/Sigma and cosine for the Bc candidate
        math::XYZVector BcPerp(Bc_Pi.px(), Bc_Pi.py(), 0.);
//         GlobalPoint BcVertexPoint   = BcVertex->position();
        GlobalError BcVertexError 	= BcVertex->error();
        GlobalPoint BcDisplacementFromBeamspot(-1*((bs.x0() -  BcVertex->position().x()) + ( BcVertex->position().z() - bs.z0()) * bs.dxdz()),
 	          			     				   -1*((bs.y0() -  BcVertex->position().y()) + ( BcVertex->position().z() - bs.z0()) * bs.dydz()), 0);
       
        double LBcBSXY 		= BcDisplacementFromBeamspot.perp();
        double LBcBSErr		= sqrt(BcVertexError.rerr(BcDisplacementFromBeamspot));
        reco::Vertex::Point   vperpBc(BcDisplacementFromBeamspot.x(),BcDisplacementFromBeamspot.y(),0.);
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
        
        //Set TLorentz vectors for Ntuples    
		TLorentzVector *recoBc    = new TLorentzVector(Bc_Pi.px() , Bc_Pi.py() , Bc_Pi.pz() , Bc_Pi.energy() );
		TLorentzVector *recoB     = new TLorentzVector(Bc_K.px()  , Bc_K.py()  , Bc_K.pz()  , Bc_K.energy()  );
		TLorentzVector *recoMuP   = new TLorentzVector(MuPp4R.px(), MuPp4R.py(), MuPp4R.pz(), MuPp4R.energy());
		TLorentzVector *recoMuM   = new TLorentzVector(MuMp4R.px(), MuMp4R.py(), MuMp4R.pz(), MuMp4R.energy());
		TLorentzVector *recoPi    = new TLorentzVector(Pip4R.px() , Pip4R.py() , Pip4R.pz() , Pip4R.energy() );

		reco::Candidate::LorentzVector JJpsi = MuPp4R+MuMp4R;
		TLorentzVector *recoJpsi  = new TLorentzVector(JJpsi.px() , JJpsi.py() , JJpsi.pz() , JJpsi.energy());
		TLorentzVector *recoJpsiV = new TLorentzVector(jpsiV.px() , jpsiV.py() , jpsiV.pz() , jpsiV.energy());

		if (jpsiV.px()!=jpsiX && jpsiV.py()!=jpsiY &&jpsiV.pz()!=jpsiZ)
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
        theNewCand.SetBcPi 				(*recoBc															);
        theNewCand.SetBcK  				(*recoB 															);
        //Jpsi cand
        theNewCand.SetJpsi 				(*recoJpsi 															);
        theNewCand.SetJpsiV				(*recoJpsiV 														);
        theNewCand.SetDCA  				(DCA	 															);
        //Muons
        theNewCand.SetMuP  				(*recoMuP 															);
        theNewCand.SetMuM  				(*recoMuM 															);
		theNewCand.SetMuPisGlobal   	(recomuP->isGlobalMuon() 											);
		theNewCand.SetMuPisTracker  	(recomuP->isTrackerMuon()  											);
		theNewCand.SetMuPisPFlow    	(recomuP->isPFMuon() 												);
		theNewCand.SetMuPTMOST      	(muon::isGoodMuon(*recomuP, muon::TMOneStationTight) 				);
		theNewCand.SetMuPTrkLMeas   	(recomuP->innerTrack()->hitPattern().trackerLayersWithMeasurement() );
		theNewCand.SetMuPPixLMeas   	(recomuP->innerTrack()->hitPattern().pixelLayersWithMeasurement() 	);
		theNewCand.SetMuPPixHits    	(recomuP->innerTrack()->hitPattern().numberOfValidPixelHits()  		);
		theNewCand.SetMuPTrkHits    	(recomuP->innerTrack()->hitPattern().numberOfValidTrackerHits() 	);
		theNewCand.SetMuPMatchedStat	(recomuP->numberOfMatchedStations()  								);
		theNewCand.SetMuPNormChi2   	(recomuP->innerTrack()->normalizedChi2() 							);
		theNewCand.SetMuPDxy   			(recomuP->innerTrack()->dxy(primvtx[0].position() )					);
		theNewCand.SetMuPDz   			(recomuP->innerTrack()->dz(primvtx[0].position() ) 					);
		theNewCand.SetMuMisGlobal   	(recomuM->isGlobalMuon() 											);
		theNewCand.SetMuMisTracker  	(recomuM->isTrackerMuon()  											);
		theNewCand.SetMuMisPFlow    	(recomuM->isPFMuon() 												);
		theNewCand.SetMuMTMOST      	(muon::isGoodMuon(*recomuM, muon::TMOneStationTight) 				);
		theNewCand.SetMuMTrkLMeas   	(recomuM->innerTrack()->hitPattern().trackerLayersWithMeasurement() );
		theNewCand.SetMuMPixLMeas   	(recomuM->innerTrack()->hitPattern().pixelLayersWithMeasurement() 	);
		theNewCand.SetMuMPixHits    	(recomuM->innerTrack()->hitPattern().numberOfValidPixelHits()  		);
		theNewCand.SetMuMTrkHits    	(recomuM->innerTrack()->hitPattern().numberOfValidTrackerHits() 	);
		theNewCand.SetMuMMatchedStat	(recomuM->numberOfMatchedStations()  								);
		theNewCand.SetMuMNormChi2   	(recomuM->innerTrack()->normalizedChi2() 							);
		theNewCand.SetMuMDxy   			(recomuM->innerTrack()->dxy(primvtx[0].position() )					);
		theNewCand.SetMuMDz   			(recomuM->innerTrack()->dz(primvtx[0].position() ) 					);
		// Track 
		theNewCand.SetPi            	(*recoPi															);
		theNewCand.SetPiCh   			(PionTrackCand->charge()											);
		theNewCand.SetTrkPixLMeas   	(PionTrackCand->hitPattern().pixelLayersWithMeasurement()			);
		theNewCand.SetTrkTrkLMeas   	(PionTrackCand->hitPattern().trackerLayersWithMeasurement()			);
		theNewCand.SetTrkPixHits   		(PionTrackCand->hitPattern().numberOfValidPixelHits()  				);
		theNewCand.SetTrkTrkHits    	(PionTrackCand->hitPattern().numberOfValidTrackerHits() 			);
		theNewCand.SetTrkNormChi2   	(PionTrackCand->normalizedChi2() 									);
		theNewCand.SetIP3DJpsi	    	(pi_3DIp															);
		theNewCand.SetIP3DJpsiSign  	(pi_3DIpSignificance												);
		theNewCand.SetIP2DBS	    	(pi_IPbs 															);
		theNewCand.SetIP2DBSSign    	(pi_IPbsSignificance												);
		theNewCand.SetIP3DPV	    	(pi_3DIpPV 															);
		theNewCand.SetIP3DPVSign  		(pi_3DIpPVSignificance												);
		theNewCand.SetDeltaR	  		(DRPi 																);
		// Bc vertex 
		theNewCand.SetClS		    	(ClBcVtx		);
		theNewCand.SetEl2DWrtBS     	(LBcBSXY		);
		theNewCand.SetEls2DWrtBS    	(ElsBcBSXY		);
		theNewCand.SetSigma2DWrtBS  	(LBcBSErr		);
		theNewCand.SetCos2DWrtBS    	(CosBcBsXY		);
		theNewCand.SetEl3DWrtPV     	(LBcPointPV		);
		theNewCand.SetEls3DWrtPV    	(ElsBcPointPV	);
		theNewCand.SetSigma3DWrtPV  	(LBcPointPVErr	);
		theNewCand.SetCos3DWrtPV    	(Cosine			);
		//Vertices
		for (int i=0; i < 3; i++){
		  theNewCand.SetBcVtxPosition 		(i, SVPos[i]		);
		  theNewCand.SetJpsiVtxPosition 	(i, JpsiVtxPos[i]	);
		  theNewCand.SetPointPVPosition 	(i, PointPVPos[i]	);
		}
		for (int i=0; i < 9; i++){
		  theNewCand.SetBcVtxCovariance		(i, SVCov[i]		);
		  theNewCand.SetPointPVCovariance	(i, PointPVCov[i]	);
		}
  		theNewCand.SetClJpsi	    		(ClJpsiVtx			);
  		theNewCand.SetElsigJpsi	    		(elsigJpsi			);		
  		theNewCand.SetCosJpsi	    		(cosJpsiXY			);
  		theNewCand.SetPointPVCl		        (pointPVCl		    );
		// MC matching 
  		theNewCand.SetMatchMuP	    		(matchMu[0]			);
  		theNewCand.SetMatchMuM	    		(matchMu[1]			);
  		theNewCand.SetMatchPi	    		(matchPi			);
		// HLT matching 
		for (int i=0; i < 4; i++){
  		  theNewCand.SetHltMatch 			(i,HltMatch[i]		);
		}

 		BcTree_->BcTree::AddBcTreeCand( theNewCand );
       }//pi 
     }//mu2			      
   }//mu1

   h1_["NumJPSI"      	] ->Fill(numJpsi);
   t1_["BcTree"]->Fill();
}

//=====================================================================================================================
void Bc2JpsiPi::beginJob()
{
	runNumber = 0 ;
	nEvents_  = 0 ;
 
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

 //------Histograms for JPsi--------------------------------------------  
 
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
 
//------Histograms for MC--------------------------------------------  
	MonteCarloDir->cd();   
	h1_["PYID"        	    ] = new TH1F("PYID",	     "PYTHIA ID",2000,-1000,1000) ;		     
	h1_["EVTA"        	    ] = new TH1F("EVTA",	     "Size of the generated event",1500,0,1500) ;    
	h1_["ETAMU"        	    ] = new TH1F("ETAMU",	     "Eta distbution oof generated mumu",100,-10,10) ;    
	h1_["ETAPI"        	    ] = new TH1F("ETAPI",	     "Eta distbution oof generated PI",100,-10,10) ;    
	h1_["DAUID"       	    ] = new TH1F("DAUID",	     "Bc daughters ID",2000,-1000,1000) ;
	h1_["DAUJPSIID"   	    ] = new TH1F("DAUJPSIID",    "JPSI daughters ID",1000,-500,500) ;
	h1_["XPRIMBC"     	    ] = new TH1F("XPRIMBC",		 "X primary vertex Bc", 200, -1., 1.) ;
	h1_["YPRIMBC"     	    ] = new TH1F("YPRIMBC",		 "Y primary vertex Bc", 200, -1., 1.) ;
	h1_["ZPRIMBC"     	    ] = new TH1F("ZPRIMBC",		 "Z primary vertex Bc", 300, -15.0,  15.0)  ;
	h1_["XSECBC"      	    ] = new TH1F("XSECBC", 		 "X secondary vertex Bc", 200, -1., 1.) ;
	h1_["YSECBC"      	    ] = new TH1F("YSECBC", 		 "Y secondary vertex Bc", 200, -1., 1.) ;
	h1_["ZSECBC"      	    ] = new TH1F("ZSECBC", 		 "Z secondary vertex Bc", 300, -15.0, 15.0) ;   
	h1_["ELLEBC"      	    ] = new TH1F("ELLEBC", 		 "Bc decay lenght",	 400,  -2.,   2.0)  ;
	h1_["TAUBC"				] = new TH1F("TAUBC", 		 "Bc decay time", 2000,  -10,   10)  ;
	h1_["PTOTMUP"     	    ] = new TH1F("PTOTMUP",      "Momentum of 1st muon from Jpsi", 10000, 0., 100. );
	h1_["PTOTMUM"     	    ] = new TH1F("PTOTMUM",      "Momentum of 2nd muon from Jpsi", 10000, 0., 100. );
	h1_["PTOTPI"      	    ] = new TH1F("PTOTPI",       "Momentum of PI from Bc+", 10000, 0., 100. );
	h1_["PTJPSI"      	    ] = new TH1F("PTJPSI",       "Transverse Momentum of Jpsi", 400, 0., 40. );
	h1_["PTMUP"	     	    ] = new TH1F("PTMUP",	     "Transverse momentum of 1st muon from Jpsi", 2500, 0., 25. ); 
	h1_["PTMUM"	     	    ] = new TH1F("PTMUM",	     "Transverse momentum of 2nd muon from Jpsi", 2500, 0., 25. );
	h1_["PTPI"	     	    ] = new TH1F("PTPI",	     "Transverse momentum of PI from Bc+", 2500, 0., 25. );
	h1_["COSINE"   			] = new TH1F("COSINE",	     "Cosine at gen level",  240, -1.2, 1.2) ;
	h2_["PMUPMUM"  		    ] = new TH2F("PMUPMUM",      "MuM momentum vs MuP momentum", 40, 0., 40., 40, 0., 40.);
	h2_["PTMUPMUM"    		] = new TH2F("PTMUPMUM", 	 "MuM pt vs MuP pt", 40, 0., 20., 40, 0., 20.);

    file->cd();				
    VertexingDir->cd();   

	h1_["pVertexSize"						] = new TH1F("pVertexSize"						,"pVertexSize after trigger requirement ", 		20,   0,  20   );
	h1_["pVertexTrackSize"					] = new TH1F("pVertexTrackSize"					,"pVertexTrackSize"					 	 , 	   800,   0, 800   );
	h1_["EventTrackSize" 					] = new TH1F("EventTrackSize"					,"EventTrackSize"					  	 ,	   800,   0, 800   );
	h1_["goodTrkSize"       				] = new TH1F("goodTrkSize"						,"number of tracks passing selections"	 ,     800,   0, 800   );
	h1_["Bcvtxcl"		   				 	] = new TH1F("Bcvtxcl"							,"Bc vertex CL before cut"				 ,   10000,   0,   1.  );
	h1_["Bpvtxcl"           				] = new TH1F("Bpvtxcl"							,"Bp vertex CL after cut"				 ,   10000,   0,   1.  );
	h1_["Pion_ImpactParameter"		   		] = new TH1F("Pion_ImpactParameter"				,"Pion_ImpactParameter"					 ,   10000,   0,  50   );
	h1_["Pion_IPsignificance"              	] = new TH1F("Pion_IPsignificance" 				,"Pion_IPsignificance" 					 ,   10000,   0, 100   );
	h1_["Pion_TransverseImpactParameter_BS"	] = new TH1F("Pion_TransverseImpactParameter_BS","Pion_TransverseImpactParameter_BS"	 ,    5000,   0,  25   );
	h1_["Pion_TIParameterBSsignificance"   	] = new TH1F("Pion_TIParameterBSsignificance"   ,"Pion_TIParameterBSsignificance"        ,   10000,   0,  50   );
	h1_["BcVtxFitTreeNonValid"   			] = new TH1F("BcVtxFitTreeNonValid"				,"BcVtxFitTreeNonValid"					 ,	     4,   0,   4   );
	h1_["BcVtxNonValid"		    			] = new TH1F("BcVtxNonValid"					,"BcVtxNonValid"						 ,	     4,   0,   4   );
	h1_["BcJpsiSignificance"	    		] = new TH1F("BcJpsiSignificance"				,"BcJpsiSignificance"					 ,	   100,   0,  50.0 );
	h1_["BcJpsiDistance"	        		] = new TH1F("BcJpsiDistance"					,"BcJpsiDistance"						 ,   10000,   0,   5.0 );
	h1_["PointPVDelTrks"	        		] = new TH1F("PointPVDelTrks"					,"Number of deleted trks from PointPV"   ,       5,   0,   5.0 );
	h1_["deltaRPi"	        				] = new TH1F("deltaRPi"							,"deltaRPi"								 ,    1000,   0,  10.0 );

	file->cd();
	h1_["filter"      		      	] = new TH1F("filter",	   		 "Binned filter counter", 50,0, 50); 
	h1_["TrueNumInteraction"      	] = new TH1F("TrueNumInteraction",  "TrueNumInteraction", 50,0, 50); 
	h1_["NumJPSI"      				] = new TH1F("NumJPSI",			               "NumJPSI", 50,0, 50); 

    TH1::AddDirectory(oldAddDir);
}

//============================================================================================================
void Bc2JpsiPi::endJob() {
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
void Bc2JpsiPi::MonteCarloStudies(const edm::Event& iEvent)
{
   	if( iEvent.isRealData() ) {return ;}
   	bool isBc = false ;
   	bool isBp = false ;
   	int  id;
   	int BcplusId      = 0 ;
   	int piplus        = 0 ;
   	double BcMassMC   = 0.;

   	std::string mcpart = pString_["MCID"]; 
   	std::string strBc("Bc");
   	std::string strBp("Bp");
   	if (mcpart.find(strBc) != std::string::npos) 
   	{
   	  BcplusId =   541;
   	  piplus   =   211;
   	  BcMassMC = 6.275;
   	  isBc     =  true;
   	}
   	else if (mcpart.find(strBp) != std::string::npos) 
   	{
   	  BcplusId =    521;
   	  piplus   =    321;
   	  BcMassMC = 5.2789;
   	  isBp     =   true;
   	}
    boolMC        = false;
	bool boolJpsi = false;
	bool boolpi   = false;
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
//    	h1_["EVTA"]->Fill(genParticles.size());
		
	int quanteBc=0;
   	ELLE = -1.;
    for ( size_t i=0; i< genParticles->size(); ++i) 
    { 
	  const reco::GenParticle &p = (*genParticles)[i];
	  id = p.pdgId();
	  h1_["PYID"]->Fill(id);
	  if (isBc)
	  {
		if (id != BcplusId ) 		continue;
	  }
	  else if (isBp)
	  {
		if(fabs(id) != BcplusId ) 	continue; 
	  }        
	  quanteBc++;
	  if (isBc)
	  {
		if (quanteBc!=2) 			continue;
	  }

	  // return all daughters of Bc 
	  if (p.numberOfDaughters()!=2 ) continue;  
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
		
		bool dauPi = false;
		if (isBc && dauId == piplus) 		 dauPi= true;
		if (isBp && fabs(dauId) == piplus)   dauPi= true;
		if( dauPi)  
		{
		  p_pi.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
		  PICH = des->charge();
		  h1_["ETAPI"]->Fill(p_pi.Eta());
		  boolpi = true;
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

	  if(boolJpsi && boolpi) boolMC = true;
    } // end MonteCarlo loop

    double ptotmup = sqrt(p_muP.Px()*p_muP.Px() + p_muP.Py()*p_muP.Py() + p_muP.Pz()*p_muP.Pz());
    double ptotmum = sqrt(p_muM.Px()*p_muM.Px() + p_muM.Py()*p_muM.Py() + p_muM.Pz()*p_muM.Pz());
    double ptotpi  = sqrt(p_pi.Px()*p_pi.Px()   + p_pi.Py()*p_pi.Py()   + p_pi.Pz()*p_pi.Pz());
    p_jpsi         = p_muP + p_muM;

    h1_["PTOTMUP"	] -> Fill(ptotmup);
    h1_["PTOTMUM"	] -> Fill(ptotmum);
    h1_["PTOTPI" 	] -> Fill(ptotpi);
    h1_["PTMUP"  	] -> Fill(p_muP.Perp());
    h1_["PTMUM"  	] -> Fill(p_muM.Perp());
    h1_["PTPI"   	] -> Fill(p_pi.Perp());
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

//=============================================================================================================================================== 
bool Bc2JpsiPi::passPioncuts(reco::TrackRef pionCand)
{
	 if (pionCand->pt()<pDouble_["cut_Pt_Trk"])								 return false ;
     if (pionCand->p()<pDouble_["pCandMomCut"])								 return false ;
     if (fabs(pionCand->eta())>pDouble_["cut_eta_pi"])						 return false ;
	 int nhits=pionCand->numberOfValidHits();
  	 if( nhits <pDouble_["cut_nhits"] )										 return false ;
 	 if (pionCand->hitPattern().numberOfValidPixelHits() < pDouble_["numberOfValidPixelHitsCut"]  )	         return false ;
     if( pionCand->normalizedChi2() > pDouble_["cut_chi2n"])				 return false ;
         
	 return true;
}
//============================================================================================== 
//Acquire Input Parameters
void Bc2JpsiPi::acquireInputParameters(void)
{
   std::cout << "\n\n"
             << __LINE__ << "]\t" 
             << __PRETTY_FUNCTION__  
	     << "\t=================== Parameters read from python configuration file ============================="
	     << endl << endl ;
   
   map<string, vector<string> > containers ;
   
   containers["double"] = inputDouble_ ;
   containers["string"] = inputString_ ;
   for(map<string, vector<string> >::iterator typeIt =containers.begin(); typeIt!=containers.end(); ++typeIt)
   {
     for(std::vector<std::string>::iterator it =typeIt->second.begin(); it!=typeIt->second.end(); ++it)
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
//end input parameter reading
//========================================================================================
void Bc2JpsiPi::checkHLT(const edm::Event& iEvent)
{
   std::string theName   = pString_["HLTMatchName"];
   std::string trigpart1 = pString_["HLTMatchModule1"];   
   std::string trigpart2 = pString_["HLTMatchModule2"];

   edm::Handle<trigger::TriggerEvent> theTriggerEvent;
   iEvent.getByLabel(theName,theTriggerEvent);
   std::vector<reco::Particle::LorentzVector>       JpsiReco;
   std::vector<std::pair< double,double> >          MuonTrig;
   JpsiReco.clear();
   MuonTrig.clear();
   JPsiTriggered_.clear() ;
 
   const trigger::TriggerObjectCollection& theTriggerObjects(theTriggerEvent->getObjects());

   int count = 0 ;
   for(size_t trigFilter = 0; trigFilter < theTriggerEvent->sizeFilters(); trigFilter++)
   {
     std::string triggerName = theTriggerEvent->filterTag(trigFilter).encode();
     if (boolTrg_) 	std::cout << "filter:" << triggerName << std::endl;
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
	   reco::Particle::LorentzVector JpsiReco_temp;
   
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
     else
     {
     }
   }
}
//===============================================================================================================
//define this as a plug-in
DEFINE_FWK_MODULE(Bc2JpsiPi);
