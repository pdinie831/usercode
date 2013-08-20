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
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HeavyFlavorAnalysis/Bc2JpsiPi/interface/Bc2JpsiPi.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h" //sara

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"   
#include "HepMC/GenParticle.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 
#include "CLHEP/Vector/LorentzVector.h"

#include "TMath.h"

#include <boost/regex.hpp> 

#include <math.h>
#include <typeinfo>
#include <DataFormats/Math/interface/deltaR.h>

// constants, enums and typedefs
   double c_const  = 0.0299792458;//cm/ps

// static data member definitions
   double xprimevt = 0.;
   double yprimevt = 0.;
   double zprimevt = 0.;
   double xprimBc  = 0.;
   double yprimBc  = 0.;
   double zprimBc  = 0.;
   double xsecBc   = 0.;
   double ysecBc   = 0.;
   double zsecBc   = 0.;
   double elleBc   = 0.;
   double TAU      = 0.;
   double COSENO   = 0.;
   TLorentzVector p_mu1,p_mu2, p_pi, p_jpsi,p_Bc;

   int PiGENch   = 0;
   double PVMC[3]   = {0,0,0};
   double SVMC[3]   = {0,0,0};
   int    MuGENch[2]   = {0,0};
   bool boolMC      = false;

// constructors and destructor
using namespace std ;

Bc2JpsiPi::Bc2JpsiPi(const edm::ParameterSet& iConfig) :
   inputDouble_(iConfig.getParameter<std::vector<std::string> >("inputDouble_")),
   inputString_(iConfig.getParameter<std::vector<std::string> >("inputString_")),
   thePVs_(     iConfig.getParameter<edm::InputTag>	       ("primaryVertexTag")),
   thebeamspot_(iConfig.getParameter<edm::InputTag>	       ("beamSpotTag"     ))
{
  cout << __LINE__ << "]\t" << __PRETTY_FUNCTION__ << "\tCtor" << endl ;
  thisConfig_ = iConfig;
  acquireInputParameters() ;
}	

//==============================================================================================
Bc2JpsiPi::~Bc2JpsiPi()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)  
}

//=============== method called to for each event ===================================================================
void Bc2JpsiPi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  runNumber++;
  nEvents_++;
  cout << __LINE__ << "]\t" << __PRETTY_FUNCTION__ << "\trunNumber: " << runNumber << endl ;
  BcTree_    ->BcTree::BcTreeClear();

  edm::Handle<reco::BeamSpot>                    theBeamSpot;
  iEvent.getByLabel(thebeamspot_,                theBeamSpot	);
  reco::BeamSpot bs = *theBeamSpot;

  if (!iEvent.isRealData())
  {
  MonteCarloStudies(iEvent) ;
  }

  double       BeamSpot[3]={0,0,0};
  BeamSpot[0]	= bs.x0(); 
  BeamSpot[1]	= bs.y0(); 
  BeamSpot[2]	= bs.z0(); 

  if(!iEvent.isRealData())
  {
	BcTree_->BcTree::SetBcTreeGENCand(  
					p_Bc 		,
					p_mu1 		,
					p_mu2 		,
					p_jpsi		,
					p_pi 		,
					elleBc      ,
					COSENO      ,
					SVMC        ,
				    BeamSpot    ,
					PiGENch     ,
					MuGENch
					);
  }

  int counter = 0 ;
  h1_["filter"]->Fill(counter);//0

  const reco::Muon* recomu1   = 0;
  const reco::Muon* recomu2   = 0;

  typedef ROOT::Math::SVector<double, 3>    Vector3;
  typedef ROOT::Math::SMatrix<double, 3, 3> Matrix33;
  Matrix33 Covsec;
  Matrix33 Covprim;
  Vector3  Deriv;
 	
  int                 NumBc       = 0;// remove header from list
  double			  bestVtxS[3] = {0,0,0};
  int                 MuCh[2]     = {0,0};
  int                 triggerBit  = 0;

  int numJPSI    = 0 ;
  double jpsiX   = 0 ;
  double jpsiY   = 0 ;
  double jpsiZ   = 0 ;
  int JpsiGood   = 0 ;
  int JpsiNoTrig = 0 ; 
  int JpsiTrig   = 0 ; 

  if(runNumber == 1)
  {
    if( iEvent.isRealData())
	{
	  cout << "\n\n"
	  << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t"
	  << "These are real data" 
	  << endl ;
	}
    else
	{
	  cout << "\n\n"
	  << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t"
	  << "These are MC data" 
	  << endl ;
	}
  } 
  if( (runNumber % 100) == 0 ) 
  {
    cout << __LINE__ << "]\t" << __PRETTY_FUNCTION__ << "\t" 
    << "Event " << iEvent.id() << " Run Number " << runNumber << endl;
  }

  edm::Handle<reco::MuonCollection>   recoMuons     ; 
  edm::Handle<reco::TrackCollection>  muonTracks    ;
  iEvent.getByLabel("muons",          recoMuons )	;
  iEvent.getByLabel("globalMuons",    muonTracks)	;

  edm::ESHandle<GlobalTrackingGeometry>          theTrackingGeometry;
  edm::ESHandle<TransientTrackBuilder>           theB;
  edm::Handle<reco::TrackCollection>             tracks ;
  edm::Handle<reco::TrackCollection>             pvtracks;   
  edm::Handle<reco::VertexCollection>            priVtxs;

  iEvent.getByLabel(pString_["trackCollection"], tracks	);
  iEvent.getByLabel(thePVs_,                     priVtxs	);

  iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

  //Get info for MC PileUp reweighting 
  float Tnpv = -1;
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
   

  //Trigger   
  edm::Handle<edm::TriggerResults> hltresults;
  edm::InputTag tag("TriggerResults","","HLT");//HLT REDIGI39X REDIGI38XPU
  iEvent.getByLabel(tag,           hltresults);

  if (!hltresults.isValid() ) 
  {
	cout << __LINE__ << "\tcould not get the collection!" << endl;
	if (boolMC) t1_["BcTree"]->Fill();
	return;
  }
  edm::TriggerNames triggerNames_ = iEvent.triggerNames(*hltresults);
  int ntrigs=hltresults->size();
  vector<string> triggernames = triggerNames_.triggerNames();

  bool triggerFound = false ;
  std::vector<int>         goodTrig ;
  
  for(int itrig = 0; itrig != ntrigs; ++itrig)
  { 
//   std::cout << triggerNames_.triggerName(itrig) << std::endl;
	 if((triggerNames_.triggerName(itrig) == pString_["HLTname1"]) || 
		(triggerNames_.triggerName(itrig) == pString_["HLTname2"]) ||
		(triggerNames_.triggerName(itrig) == pString_["HLTname3"]) 
	 ) 
	   {
		triggerFound = true ;
	 goodTrig.push_back(itrig) ;
   }
  }
    
  if ( !triggerFound ) {
	if (boolMC) t1_["BcTree"]->Fill();
	return ;
  }
  
  bool triggerAccepted = false ;
  for(unsigned int t=0; t<goodTrig.size(); ++t)
  {
	if( hltresults->accept(goodTrig[t]) ) 
	{
	  triggerAccepted = true ;
	}
  }

  if ( triggerAccepted ) triggerBit = 1;

  //Pre-selection cuts on track collection
  newTracksDef qualityTracks;
  for (uint tracksIt =0 ;  tracksIt < tracks->size(); tracksIt++)
  {
	 reco::TrackRef checkTrk(tracks,tracksIt) ;										        
	 bool flag = false ; 												        
	 for (reco::MuonCollection::const_iterator mu =recoMuons->begin(); mu!=recoMuons->end(); mu++)			        
	 {														        
	   if (mu->track().get()!= 0 && mu->track().get() == checkTrk.get()) { flag=true; break;}			        
	 }														        
	 if (flag)												   continue;	        
																   
	 if (checkTrk->pt()  				  					< pDouble_["cut_Pt_Trk"])			   			continue;	        
	 if (fabs(checkTrk->eta())			  					> pDouble_["cut_eta_pi"])			  			continue;	   						  
	 if (checkTrk->numberOfValidHits()	  					< pDouble_["cut_nhits"] )			   			continue;	   					   
	 if (checkTrk->hitPattern().numberOfValidPixelHits() 	< pDouble_["numberOfValidPixelHitsCut"]  )	    continue;	        
	 if (checkTrk->normalizedChi2()			 				> pDouble_["cut_chi2n"])			   		    continue; 	        
	 if (!iEvent.isRealData() && checkTrk->charge() !=1)  													continue;   //per Bc
	 qualityTracks[tracksIt]=checkTrk;											   	
  }
 
  //set BcTreeHeader   BcNum=0
  BcTree_->BcTree::SetBcTreeHeader(  nEvents_		       ,
				   (int32_t)iEvent.id().run()            ,
				   (int32_t)iEvent.id().luminosityBlock(),
				   (int32_t)iEvent.id().event()	       ,
					NumBc		                           ,
				   (int32_t)tracks ->size()              ,
				   (int32_t)priVtxs->size()              ,
				   Tnpv                                  ,
				   triggerBit                            , 
				   qualityTracks.size()
				   );

  if ( !triggerAccepted ) {
    if (boolMC) t1_["BcTree"]->Fill();
    return ;
  }
  this->checkHLT(iEvent) ;  

  h1_["filter"      ]->Fill(counter++);//1
  if (recoMuons->size()<2) {
	if (boolMC) t1_["BcTree"]->Fill();
	return ;
  }
    
  const reco::VertexCollection primvtx = *(priVtxs.product() );    
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
	reco::TrackRef refmu1t=recomu1->track(); 
	
	if (fabs(muR1->eta())>pDouble_["cut_eta"]) 							         					continue;
	if (fabs(muR1->pt() )<pDouble_["cut_Pt_Mu"])							         				continue;

	//POG SOFT Muon Selection 2011
	if (!muon::isGoodMuon(*recomu1,muon::TMOneStationTight)                    )				 	continue;
	if (recomu1->innerTrack()->hitPattern().numberOfValidTrackerHits()   <= 10 )					continue;
	if (recomu1->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 1  )					continue;
	if (recomu1->innerTrack()->normalizedChi2() 						 	>= 1.8)	 				continue;
	if (!(fabs(recomu1->innerTrack()->dxy(primvtx[0].position())) < 3. && 
		  fabs(recomu1->innerTrack()->dz(primvtx[0].position()))  < 30. )) 							continue;
	
    //Loop on the second muon 
    for(reco::MuonCollection::const_iterator muR2 =muR1+1; muR2!=recoMuons->end(); ++muR2)
    {
	  recomu2 = &*muR2;
	  reco::TrackRef refmu2t=recomu2->track();

	  if (muR1->charge()+muR2->charge()!=0)									   						continue; // Consider oppsite sign pairs
	  if (fabs(muR2->eta())>pDouble_["cut_eta"])								    				continue;
	  if (fabs(muR2->pt() )<pDouble_["cut_Pt_Mu"])													continue;
	
      //POG SOFT Muon Selection
	  if (!muon::isGoodMuon(*recomu2,muon::TMOneStationTight)                    )				 	continue;
	  if (recomu2->innerTrack()->hitPattern().numberOfValidTrackerHits()   <= 10 )					continue;
	  if (recomu2->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 1  )					continue;
	  if (recomu2->innerTrack()->normalizedChi2() 						   >= 1.8)	 				continue;
	  if (!(fabs(recomu2->innerTrack()->dxy(primvtx[0].position())) < 3. && 
			fabs(recomu2->innerTrack()->dz(primvtx[0].position()))  < 30. )) 						continue;
 
	  if( !recomu1->track().isNonnull() || !recomu2->track().isNonnull())
	  {
		assert(0) ;
	  }

	  recomu1 = &*muR1;
	  recomu2 = &*muR2;
	  const reco::TransientTrack tTrackmu1=(*theB).build(refmu1t.get());
	  const reco::TransientTrack tTrackmu2=(*theB).build(refmu2t.get());

	  reco::Candidate::LorentzVector tr1p4(refmu1t->px(), 
										   refmu1t->py(),
										   refmu1t->pz(),
										   sqrt(refmu1t->p()*refmu1t->p() + recomu1->mass()*recomu1->mass()));
	  reco::Candidate::LorentzVector tr2p4(refmu2t->px(),
										   refmu2t->py(),
										   refmu2t->pz(),
										   sqrt(refmu2t->p()*refmu2t->p() + recomu2->mass()*recomu2->mass()));

      //Trigger requirements - DCA
	  TrajectoryStateClosestToPoint mu1TS = tTrackmu1.impactPointTSCP();
	  TrajectoryStateClosestToPoint mu2TS = tTrackmu2.impactPointTSCP();
	  if (mu1TS.isValid() && mu2TS.isValid()) 
	  {
		ClosestApproachInRPhi cApp;
		cApp.calculate(mu1TS.theState(), mu2TS.theState());
		if (!cApp.status() || cApp.distance() > 0.5) continue;
	  }
	  h1_["filter"      ]->Fill(counter++);//2

	  reco::Candidate::LorentzVector jpsV = tr1p4 + tr2p4;
	  double JPsiMass= jpsV.M();
	  double ptJpsi = jpsV.pt() ;

	  //Fit to the dimuon vertex
	  KalmanVertexFitter fitter1;
	  TransientVertex    mumuVertex;
	  vector<reco::TransientTrack> tTrMu;
	  tTrMu.push_back(tTrackmu1);										  
	  tTrMu.push_back(tTrackmu2);
	  mumuVertex=fitter1.vertex(tTrMu);
	  if(!mumuVertex.isValid()) 													continue;
	 
	  //Evaluate Jpsi L/Sigma and cosine
	  math::XYZVector pperp(refmu1t->px() + refmu2t->px(), refmu1t->py() + refmu2t->py(), 0.);
	  GlobalPoint JpsiVertex      = mumuVertex.position();
	  GlobalError JpsiVertexError = mumuVertex.positionError();
	  GlobalPoint displacementFromBeamspot(-1*((bs.x0() - mumuVertex.position().x()) + (mumuVertex.position().z() - bs.z0()) * bs.dxdz()),
										   -1*((bs.y0() - mumuVertex.position().y()) + (mumuVertex.position().z() - bs.z0()) * bs.dydz()), 0);
	 
	  double LPromptsecXY 	= displacementFromBeamspot.perp();
	  double lxyerr 	 	= sqrt(JpsiVertexError.rerr(displacementFromBeamspot));
	  reco::Vertex::Point     vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
	  double cosJpsiXY 		= vperp.Dot(pperp)/(vperp.R()*pperp.R());
	  double elsigJpsi  	= 0 ;
	  try	      { elsigJpsi= LPromptsecXY / lxyerr ; }
	  catch (...) {cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t" << "Floating divide by zero!" << endl ;}

	  h1_["filter"      ]-> Fill(counter++);//3
	  h1_["InvMassJPsi" ]-> Fill(jpsV.M() );
	  h1_["JpsiPt"      ]-> Fill(ptJpsi   ); 
	  h1_["ElsigJpsi"   ]-> Fill(elsigJpsi); 

      //Trigger requirements: CL, L/Sigma, Cosine and pT
	  float vProbJpsi(TMath::Prob(mumuVertex.totalChiSquared(),(int)(mumuVertex.degreesOfFreedom())));
	  if (vProbJpsi <  pDouble_["cut_cl_Jpsi"])            							continue;
	  h1_["InvMassJPsicl"       ]->Fill(jpsV.M());

	  if (elsigJpsi <= 3   )   														continue;
	  if (cosJpsiXY <= 0.9 )   														continue;
	  if (ptJpsi    <=  6.9)														continue;
	  if( fabs(JPsiMass-pDouble_["JPsiMassPDG"])> pDouble_["JPsiMassWindow"]) 		continue;
	  JpsiGood++;

      //Trigger Matching
	  JpsiNoTrig++ ; 
	  double hltMatchDR = 0.5; 															     	
	  bool MuTrigMatch  = false;		
	  double hltMatch[4];													     	
																				  
	  for(JTrigDef::iterator it=JPsiTriggered_.begin(); it!=JPsiTriggered_.end(); ++it)								     	
	  {																		     	
		double etaMu1 = (*it).first.first   ;														     	
		double phiMu1 = (*it).first.second  ;														     	
		double etaMu2 = (*it).second.first  ;														     	
		double phiMu2 = (*it).second.second ;														     	
		if( ((reco::deltaR(refmu1t->eta(), refmu1t->phi(), etaMu1, phiMu1 ) < hltMatchDR) &&								     	
			 (reco::deltaR(refmu2t->eta(), refmu2t->phi(), etaMu2, phiMu2 ) < hltMatchDR)) ||								     	
			((reco::deltaR(refmu1t->eta(), refmu1t->phi(), etaMu2, phiMu2 ) < hltMatchDR) &&								     	
			 (reco::deltaR(refmu2t->eta(), refmu2t->phi(), etaMu1, phiMu1 ) < hltMatchDR))									     	
		  ) 																		     	
		{																		     	
		  JpsiTrig++; 
		  MuTrigMatch = true;	
		  if (MuTrigMatch)
		  {		
			hltMatch[0] = reco::deltaR(refmu1t->eta(), refmu1t->phi(), etaMu1, phiMu1 );
			hltMatch[1] = reco::deltaR(refmu2t->eta(), refmu2t->phi(), etaMu2, phiMu2 );
			hltMatch[2] = reco::deltaR(refmu1t->eta(), refmu1t->phi(), etaMu2, phiMu2 );
			hltMatch[3] = reco::deltaR(refmu2t->eta(), refmu2t->phi(), etaMu1, phiMu1 );
		  }
		}																		     	
	  }																		     	
  
	  h1_["JpsiBeforeTrigMatch" 		]->Fill(JpsiNoTrig);
	  if(!MuTrigMatch )					continue;
	  h1_["JpsiAfterTrigMatch"  		]->Fill(JpsiTrig  );
	  h1_["InvMassJpsiPassingTrigMatch" ]->Fill(jpsV.M()  );

	  MuCh[0] = muR1->charge();
	  MuCh[1] = muR2->charge();

	  TLorentzVector *jpsiMu1  = new TLorentzVector(tr1p4.px(), tr1p4.py(), tr1p4.pz(), tr1p4.energy());
	  TLorentzVector *jpsiMu2  = new TLorentzVector(tr2p4.px(), tr2p4.py(), tr2p4.pz(), tr2p4.energy());
	  TLorentzVector *jpsiTV   = new TLorentzVector(jpsV.px(), jpsV.py(), jpsV.pz(), jpsV.energy());
	  double JpsiVtx[3];
	  JpsiVtx[0] = mumuVertex.position().x() ;
	  JpsiVtx[1] = mumuVertex.position().y() ;
	  JpsiVtx[2] = mumuVertex.position().z() ;

      //MC match Muons
	  double matchMuJpsi[2];
	  if(!iEvent.isRealData())
	  {
		if(muR1->charge() == MuGENch[0])
		{
		  matchMuJpsi[0] = deltaR(tr1p4.eta(),tr1p4.phi(),p_mu1.Eta(),p_mu1.Phi());
		  matchMuJpsi[1] = deltaR(tr2p4.eta(),tr2p4.phi(),p_mu2.Eta(),p_mu2.Phi());
		}
		else 
		{
		  matchMuJpsi[0] = deltaR(tr2p4.eta(),tr2p4.phi(),p_mu1.Eta(),p_mu1.Phi());
		  matchMuJpsi[1] = deltaR(tr1p4.eta(),tr1p4.phi(),p_mu2.Eta(),p_mu2.Phi());
		}	
	  }

      //set JpsiCand
	  BcTree_->BcTree::AddJpsiCand(  
					  *jpsiMu1	    ,
					  *jpsiMu2	    ,
					  *jpsiTV		,
					  elsigJpsi     ,
					  cosJpsiXY     ,
					  vProbJpsi     ,
					  JpsiVtx       ,
					  matchMuJpsi   ,
					  MuCh
					  );

//----------------------------------------------------------------------------
      //Loop on the tracks that have passed pre-selection cuts
	  unsigned int trkI = -1 ;
	  for (newTracksDef::const_iterator it_pi=qualityTracks.begin(); it_pi!=qualityTracks.end(); ++it_pi)
	  {
		trkI=(*it_pi).first ;
		reco::TrackRef PionTrackCand((*it_pi).second) ;
		const reco::TransientTrack piTrackCand = (*theB).build(PionTrackCand);
		
		pair<double,double>  pion_IPPair = pionImpactParameter(piTrackCand,mumuVertex);
		double pi_3DIp                   = pion_IPPair.first;
		double pi_3DIpSignificance       = pion_IPPair.second;
		h1_["Pion_ImpactParameter"]->Fill(pi_3DIp);
		h1_["Pion_IPsignificance" ]->Fill(pi_3DIpSignificance);
		
		GlobalPoint BeamSpotGP(bs.x0(), bs.y0(), bs.z0()); 
		pair<double,double>  Pion_BSPair = pionIPBeamSpot(piTrackCand,BeamSpotGP);
		double pi_IPbs                   = Pion_BSPair.first;
		double pi_IPbsSignificance       = Pion_BSPair.second;
		h1_["Pion_TransverseImpactParameter_BS" ]-> Fill(pi_IPbs);
		h1_["Pion_TIParameterBSsignificance"    ]-> Fill(pi_IPbsSignificance);

		double DRPi =  deltaR(PionTrackCand->eta(),PionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
		h1_["deltaRPi"     ]-> Fill(DRPi);
//         if (DRPi > 3) continue;

		int trkNHits     = PionTrackCand->numberOfValidHits();
		int trkPixelHits = PionTrackCand->hitPattern().numberOfValidPixelHits();
		float trkChi2    = PionTrackCand->normalizedChi2();

        //Kinematic vertex fit with Jpsi constrain for Bc->JpsiPi	
		KinematicParticleFactoryFromTransientTrack pFactory;
		const ParticleMass muon_mass  = 0.10565837; 
		ParticleMass jpsimass_c 	  = 3.096916; 
		ParticleMass pion_mass  	  = 0.13957018; 
		ParticleMass kaon_mass        = 0.493677;
		float muon_sigma 			  = muon_mass*1.e-6;
		float pion_sigma        	  = pion_mass*1.e-6; 
		float kaon_sigma              = kaon_mass*1.e-6; 
		float chi 		 			  = 0.;
		float ndf 		 			  = 0.;

		std::vector<RefCountedKinematicParticle> XParticles_Bc;
		XParticles_Bc.push_back(pFactory.particle(tTrackmu1,muon_mass,chi,ndf,muon_sigma));
		XParticles_Bc.push_back(pFactory.particle(tTrackmu2,muon_mass,chi,ndf,muon_sigma));
		XParticles_Bc.push_back(pFactory.particle(piTrackCand,pion_mass,chi,ndf,pion_sigma)); 

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
	 
		double vProbBc = ChiSquaredProbability((double)(BcVertex->chiSquared()),(double)(BcVertex->degreesOfFreedom()));
		h1_["Bcvtxcl"]-> Fill(vProbBc);
		if (vProbBc <= pDouble_["cut_cl_Bc"])                                       continue ;

		BcVertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle muP =  BcVertexFitTree->currentParticle();
		BcVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle muN =  BcVertexFitTree->currentParticle();
		BcVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle pi =  BcVertexFitTree->currentParticle();
		BcVertexFitTree->movePointerToTheNextChild();

		reco::Candidate::LorentzVector mu1p4R( muP->currentState().globalMomentum().x(),muP->currentState().globalMomentum().y(), 
											   muP->currentState().globalMomentum().z(),muP->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector mu2p4R( muN->currentState().globalMomentum().x(),muN->currentState().globalMomentum().y(), 
											   muN->currentState().globalMomentum().z(),muN->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector pip4R(  pi->currentState().globalMomentum().x(),pi->currentState().globalMomentum().y(), 
											   pi->currentState().globalMomentum().z(),pi->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector Bc_Pi(  BcCand->currentState().globalMomentum().x(),BcCand->currentState().globalMomentum().y(), 
											   BcCand->currentState().globalMomentum().z(),BcCand->currentState().kinematicParameters().energy());


        //Kinematic vertex fit with Jpsi constrain for B+->JpsiK
		chi=0;
		ndf=0;
		pion_sigma = pion_mass*1.e-6; 
		kaon_sigma = kaon_mass*1.e-6; 

		std::vector<RefCountedKinematicParticle> XParticles_Bp;
		XParticles_Bp.push_back(pFactory.particle(tTrackmu1,muon_mass,chi,ndf,muon_sigma));
		XParticles_Bp.push_back(pFactory.particle(tTrackmu2,muon_mass,chi,ndf,muon_sigma));
		XParticles_Bp.push_back(pFactory.particle(piTrackCand,kaon_mass,chi,ndf,kaon_sigma)); 
		KinematicConstrainedVertexFitter kcvFitterBp;
		RefCountedKinematicTree BpVertexFitTree = kcvFitterBp.fit(XParticles_Bp, jpsiconstraint); 

		if (!BpVertexFitTree->isValid())  		continue;
		
		BpVertexFitTree->movePointerToTheTop();
		RefCountedKinematicParticle BpCand = BpVertexFitTree->currentParticle();
		RefCountedKinematicVertex BpVertex = BpVertexFitTree->currentDecayVertex();
		if (!BpVertex->vertexIsValid())			continue;
		
		double vProbBp=ChiSquaredProbability((double)(BpVertex->chiSquared()),(double)(BpVertex->degreesOfFreedom()));
	     
		h1_["Bpvtxcl"  ]->Fill(vProbBp);
		if ( vProbBp<pDouble_["cut_cl_Bc"] ) 	continue;
	   
		BpVertexFitTree->movePointerToTheFirstChild();
		RefCountedKinematicParticle muP_2 =  BpVertexFitTree->currentParticle();
		BpVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle muN_2 = BpVertexFitTree->currentParticle();
		BpVertexFitTree->movePointerToTheNextChild();
		RefCountedKinematicParticle kaon1 = BpVertexFitTree->currentParticle();

		reco::Candidate::LorentzVector kp4R( kaon1->currentState().globalMomentum().x(),kaon1->currentState().globalMomentum().y(), 
											 kaon1->currentState().globalMomentum().z(),kaon1->currentState().kinematicParameters().energy());
		reco::Candidate::LorentzVector Bc_K( BpCand->currentState().globalMomentum().x(),BpCand->currentState().globalMomentum().y(), 
											 BpCand->currentState().globalMomentum().z(),BpCand->currentState().kinematicParameters().energy());

		VertexDistance3D vertTool;
		Measurement1D BcJpsiDistance = vertTool.distance(BcVertex->vertexState(),mumuVertex.vertexState());
		double BcJpsiSignificance    = BcJpsiDistance.significance();
		double BcJpsiDistanceV       = BcJpsiDistance.value();
		h1_["BcJpsiSignificance" ]-> Fill(BcJpsiSignificance);
		h1_["BcJpsiDistance"     ]-> Fill(BcJpsiDistanceV);

        //MC match Muons ----------------------------------------------------------------------
		double matchMu[2];
		double matchPi = 10000;
		if(!iEvent.isRealData())
		{
		  if(muR1->charge() == MuGENch[0])
		  {
			matchMu[0] = deltaR(mu1p4R.eta(),mu1p4R.phi(),p_mu1.Eta(),p_mu1.Phi());
			matchMu[1] = deltaR(mu2p4R.eta(),mu2p4R.phi(),p_mu2.Eta(),p_mu2.Phi());
		  }
		  else 
		  {
			matchMu[0] = deltaR(mu2p4R.eta(),mu2p4R.phi(),p_mu1.Eta(),p_mu1.Phi());
			matchMu[1] = deltaR(mu1p4R.eta(),mu1p4R.phi(),p_mu2.Eta(),p_mu2.Phi());
		  }

          //MC match Pions ---------------------------------------------------------------------
		  if(PionTrackCand->charge() == PiGENch) 						    
		  {																    
			matchPi = deltaR(pip4R.eta(),pip4R.phi(),p_pi.Eta(),p_pi.Phi()); 
		  }																    
		  else															    
		  {																    
			matchPi = 10000; 											    
		  }																    
		}

		//Evaluate L/Sigma and cosine for the Bc candidate
        math::XYZVector BcPerp(Bc_Pi.px(), Bc_Pi.py(), 0.);
        GlobalPoint BcVertexPoint   = BcVertex->position();
        GlobalError BcVertexError 	= BcVertex->error();
        GlobalPoint BcDisplacementFromBeamspot(-1*((bs.x0() -  BcVertex->position().x()) + ( BcVertex->position().z() - bs.z0()) * bs.dxdz()),
 	          			     				   -1*((bs.y0() -  BcVertex->position().y()) + ( BcVertex->position().z() - bs.z0()) * bs.dydz()), 0);
       
        double LBcXY 		= BcDisplacementFromBeamspot.perp();
        double LBcErr 		= sqrt(BcVertexError.rerr(BcDisplacementFromBeamspot));
        reco::Vertex::Point   vperpBc(BcDisplacementFromBeamspot.x(),BcDisplacementFromBeamspot.y(),0.);
        double CosBcXY   	= vperpBc.Dot(BcPerp)/(vperpBc.R()*BcPerp.R());
    	double elsig   		= 0 ;

        try	    { elsig		= LBcXY / LBcErr ; }
		catch (...) {cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t" << "Floating divide by zero!" << endl ;}

		double xsec = BcVertex->position().x();
		double ysec = BcVertex->position().y();
		double zsec = BcVertex->position().z(); 
		Covsec(0,0) = BcVertex->error().cxx();
		Covsec(1,0) = BcVertex->error().cyx();
		Covsec(2,0) = BcVertex->error().czx();
		Covsec(0,1) = Covsec(1,0);
		Covsec(1,1) = BcVertex->error().cyy();
		Covsec(2,1) = BcVertex->error().czy();
		Covsec(0,2) = Covsec(2,0);
		Covsec(1,2) = Covsec(2,1);
		Covsec(2,2) = BcVertex->error().czz();

		double BSCovariance[9];
		double SVCovariance[9];
		int n=0;
		for (int i=0; i<3; i++)
		{
		  for (int j=0; j<3; j++)
		  { 
			SVCovariance[n] = Covsec(i,j);
		    Covprim(i,j)    = bs.covariance(i,j);  
		    BSCovariance[n] = Covprim(i,j);
			n++;
		  }
		}

        //Set Ntuple variables    
		TLorentzVector *recoBc  = new TLorentzVector(Bc_Pi.px() , Bc_Pi.py() , Bc_Pi.pz() , Bc_Pi.energy() );
		TLorentzVector *recoB   = new TLorentzVector(Bc_K.px()  , Bc_K.py()  , Bc_K.pz()  , Bc_K.energy()  );
		TLorentzVector *recoMu1 = new TLorentzVector(mu1p4R.px(), mu1p4R.py(), mu1p4R.pz(), mu1p4R.energy());
		TLorentzVector *recoMu2 = new TLorentzVector(mu2p4R.px(), mu2p4R.py(), mu2p4R.pz(), mu2p4R.energy());
		TLorentzVector *recoPi  = new TLorentzVector(pip4R.px() , pip4R.py() , pip4R.pz() , pip4R.energy() );

		reco::Candidate::LorentzVector JJpsi = mu1p4R+mu2p4R;
		TLorentzVector *recoJpsi  = new TLorentzVector(JJpsi.px(),JJpsi.py(),JJpsi.pz(),JJpsi.energy());
		TLorentzVector *recoJpsiV = new TLorentzVector(jpsV.px(),jpsV.py(),jpsV.pz(),jpsV.energy());

		if (jpsV.px()!=jpsiX && jpsV.py()!=jpsiY &&jpsV.pz()!=jpsiZ)
		{
		  jpsiX = jpsV.px();
		  jpsiY = jpsV.py();
		  jpsiZ = jpsV.pz();
		  numJPSI++;
		}
	  
		bestVtxS[0]	= xsec  ;
		bestVtxS[1]	= ysec  ;
		bestVtxS[2]	= zsec  ;


		NumBc++;
		//Fill Ntuple
		BcTree_->BcTree::AddBcTreeCand(  
						 *recoBc    			,
						 *recoB    				,
						 *recoMu1   			,
						 *recoMu2   			,
						 *recoJpsi  			,
						 *recoPi    			,
						 *recoJpsiV 			,
						 elsig	  				,
						 CosBcXY	  			,
						 vProbBc    			,
						 vProbJpsi  			,
						 BeamSpot   			,
						 bestVtxS  				,
						 BSCovariance			,
						 SVCovariance			,
						 pi_3DIp     			,
						 pi_IPbs     			,
						 pi_3DIpSignificance    ,
						 pi_IPbsSignificance    ,
						 matchMu    			,
						 matchPi    			,
						 PionTrackCand->charge(),
						 MuCh       			,
						 elsigJpsi  			,           
						 cosJpsiXY  			,
						 trkNHits				,
						 trkPixelHits			,
						 trkChi2				,
						 hltMatch
						 );
								   
								   
       }//pi 
     }//mu2			      
   }//mu1

   h1_["NumJPSI"      	] ->Fill(numJPSI);
   t1_["BcTree"]->Fill();
}

//=====================================================================================================================
// ------------ method called once each job just before starting event loop  ------------
void Bc2JpsiPi::beginJob()
{
	runNumber=0;
 
	file    = new TFile(pString_["filename"].c_str(),"recreate");
 
	nEvents_ = 0 ;
 
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
  
  
 //------Histograms for JPsi--------------------------------------------   CHECKED
 
    JPsiDir->cd();
    
    h1_["ElsigJpsi"		    			] = new TH1F("ElsigJpsi",   				"L/sigma Jpsi",					 80, 0,  20 );
    h1_["JpsiPt"		    			] = new TH1F("JpsiPt",      				"JpsiPt",						100, 0, 100 );
    h1_["JpsiPtAfterTrigger"  			] = new TH1F("JpsiPtAfterTrigger", 			"JpsiPtAfterTrigger",			100, 0, 100 );
    h1_["InvMassJPsi"               	] = new TH1F("InvMassJPsi",		   			"Invariant Mass JPsi",          500, 0,  10 );
    h1_["InvMassJPsicl" 	    		] = new TH1F("InvMassJPsicl", 				"Invariant Mass JPsi cl",		500, 0,  10 ); 
    h1_["JpsiBeforeTrigMatch"       	] = new TH1F("JpsiBeforeTrigMatch",			"Jpsi Before TrigMatch", 		 10, 0,  10 );
    h1_["JpsiAfterTrigMatch"        	] = new TH1F("JpsiAfterTrigMatch",			"Jpsi After TrigMatch", 		 10, 0,  10 );
    h1_["InvMassJpsiPassingTrigMatch" 	] = new TH1F("InvMassJpsiPassingTrigMatch", "InvMassJpsiPassingTrigMatch",	200, 1,   5 ); 
    
    file->cd();
 
//------Histograms for MC--------------------------------------------  CHECKED
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
   h1_["DIFXPRIM"    	    ] = new TH1F("DIFXPRIM", 	 "Xprim evt - Xprim Bc", 2000, -1., 1.) ;
   h1_["DIFYPRIM"    	    ] = new TH1F("DIFYPRIM", 	 "Yprim evt - Yprim Bc", 2000, -1., 1.) ;
   h1_["DIFZPRIM"    	    ] = new TH1F("DIFZPRIM", 	 "Zprim evt - Zprim Bc", 2000, -1., 1.) ;
   h1_["PTOTMU1"     	    ] = new TH1F("PTOTMU1",      "Momentum of 1st muon from Jpsi", 10000, 0., 100. );
   h1_["PTOTMU2"     	    ] = new TH1F("PTOTMU2",      "Momentum of 2nd muon from Jpsi", 10000, 0., 100. );
   h1_["PTOTPI"      	    ] = new TH1F("PTOTPI",       "Momentum of PI from Bc+", 10000, 0., 100. );
   h1_["PTJPSI"      	    ] = new TH1F("PTJPSI",       "Transverse Momentum of Jpsi", 400, 0., 40. );
   h1_["PTMU1"	     	    ] = new TH1F("PTMU1",	     "Transverse momentum of 1st muon from Jpsi", 25000, 0., 25. ); 
   h1_["PTMU2"	     	    ] = new TH1F("PTMU2",	     "Transverse momentum of 2nd muon from Jpsi", 25000, 0., 25. );
   h1_["PTPI"	     	    ] = new TH1F("PTPI",	     "Transverse momentum of PI from Bc+", 25000, 0., 25. );
   h1_["COSENO_MC"   		] = new TH1F("COSENO_MC",	 "COSENO_MC",  300, -1.5, 1.5) ;
   h1_["COSENO_MC2"  		] = new TH1F("COSENO_MC2",	 "COSENO_MC2",  300, -1.5, 1.5) ;
   h2_["PTOTMU1MU2"  		] = new TH2F("PTOTMU1MU2",   "Mu2 momentum vs Mu1 momentum", 40, 0., 40., 40, 0., 40.);
   h2_["PTMU1MU2"    		] = new TH2F("PTMU1MU2", 	 "Mu2 pt vs Mu1 pt", 40, 0., 20., 40, 0., 20.);

   file->cd();				
   VertexingDir->cd();   

   h1_["pVertexSize"						] = new TH1F("pVertexSize"						,"pVertexSize after trigger requirement ", 		20,   0,  20   );
   h1_["pVertexTrackSize"					] = new TH1F("pVertexTrackSize"					,"pVertexTrackSize"					 	 , 	   800,   0, 800   );
   h1_["EventTrackSize" 					] = new TH1F("EventTrackSize"					,"EventTrackSize"					  	 ,	   800,   0, 800   );
   h1_["goodTrkSize"       				 	] = new TH1F("goodTrkSize"						,"number of tracks passing selections"	 ,     800,   0, 800   );
   h1_["d0"			     					] = new TH1F("d0 "        						,"d0"									 ,     100,   0,   0.1 );
   h1_["dz"			     					] = new TH1F("dz "        						,"dz"									 ,   20000, -10,  10   );
   h1_["Bcvtxcl"		   				 	] = new TH1F("Bcvtxcl"							,"Bc vertex CL before cut"				 ,   10000,   0,   1.  );
   h1_["Bpvtxcl"           				 	] = new TH1F("Bpvtxcl"							,"Bp vertex CL after cut"				 ,   10000,   0,   1.  );
   h1_["bcLessPVTrackSize"					] = new TH1F("bcLessPVTrackSize"				,"bcLessPVTrackSize"					 ,	   800,   0, 800   );
   h1_["Pion_ImpactParameter"		   		] = new TH1F("Pion_ImpactParameter"				,"Pion_ImpactParameter"					 ,   10000,   0,  50   );
   h1_["Pion_IPsignificance"              	] = new TH1F("Pion_IPsignificance" 				,"Pion_IPsignificance" 					 ,   10000,   0, 100   );
   h1_["Pion_TransverseImpactParameter_BS"	] = new TH1F("Pion_TransverseImpactParameter_BS","Pion_TransverseImpactParameter_BS"	 ,    5000,   0,  25   );
   h1_["Pion_TIParameterBSsignificance"   	] = new TH1F("Pion_TIParameterBSsignificance"   ,"Pion_TIParameterBSsignificance"        ,   10000,   0,  50   );
   h1_["Pion_PVsImpactParameter"		    ] = new TH1F("Pion_PVsImpactParameter"			,"Pion_PVsImpactParameter"				 ,   10000,   0,  50   );
   h1_["BcVtxFitTreeNonValid"   			] = new TH1F("BcVtxFitTreeNonValid"				,"BcVtxFitTreeNonValid"					 ,	     4,   0,   4   );
   h1_["BcVtxNonValid"		    			] = new TH1F("BcVtxNonValid"					,"BcVtxNonValid"						 ,	     4,   0,   4   );
   h1_["BcJpsiSignificance"	    			] = new TH1F("BcJpsiSignificance"				,"BcJpsiSignificance"					 ,	   100,   0,  50.0 );
   h1_["BcJpsiDistance"	        			] = new TH1F("BcJpsiDistance"					,"BcJpsiDistance"						 ,   10000,   0,   5.0 );
   h1_["PVJpsiSignificance"	    			] = new TH1F("PVJpsiSignificance"				,"PVJpsiSignificance"					 ,     100,   0,  50.0 );
   h1_["PVJpsiDistance"	        			] = new TH1F("PVJpsiDistance"					,"PVJpsiDistance"						 ,   10000,   0,   5.0 );
   h1_["deltaRPi"	        				] = new TH1F("deltaRPi"							,"deltaRPi"								 ,    1000,   0,  10.0 );

   file->cd();
   h1_["filter"      		      	] = new TH1F("filter",	   		 "Binned filter counter", 50,0, 50); 
   h1_["TrueNumInteraction"      	] = new TH1F("TrueNumInteraction",  "TrueNumInteraction", 50,0, 50); 
   h1_["NumJPSI"      				] = new TH1F("NumJPSI",			               "NumJPSI", 50,0, 50); 


   TH1::AddDirectory(oldAddDir);
}

//============================================================================================================
// ------------ method called once each job just after ending the event loop  ------------
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
   	int id;
   	bool isBc = false ;
   	bool isBp = false ;
   	int JpsiId      = 443;
   	int BcplusId    = 0;
   	int piplus      = 0;
   	double BcMassMC = 0.;
   	std::string mcpart = pString_["MCID"]; 
   	std::string strBc("Bc");
   	std::string strBp("Bp");
   	if (mcpart.find(strBc) != std::string::npos) 
   	{
   	  BcplusId =  541;
   	  piplus   =  211;
   	  BcMassMC = 6.286;
   	  isBc     = true;
   	}
   	else if (mcpart.find(strBp) != std::string::npos) 
   	{
   	  BcplusId =  521;
   	  piplus   =  321;
   	  BcMassMC = 5.2789;
   	  isBp     = true;
   	}
    boolMC        = false;
	bool boolJpsi = false;
	bool boolpi   = false;
	
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
   	//h1_["EVTA"]->Fill(myGenEvent->particles_size());
		
	int quanteBc=0;
   	elleBc = -1.;
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
		  h1_["PTJPSI"]->Fill(des->momentum().Rho());		
		  xsecBc  = des->vertex().x()/10.;
		  SVMC[0] = xsecBc;
		  ysecBc  = des->vertex().y()/10.;
		  SVMC[1] = ysecBc;
		  zsecBc  = des->vertex().z()/10.; 
		  SVMC[2] = zsecBc;
		  boolJpsi = true;

		  // return all daughters of J/psi    
		  int quantemu = 0.;
		  if (des->numberOfDaughters()!=2)  		continue;  
		  for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
		  {
			const reco::Candidate *mumu = des->daughter(imumu);
			quantemu++;
			if (fabs(mumu->pdgId())!= 13) 		continue; 
			h1_["DAUJPSIID"]->Fill(mumu->pdgId());    		
			h1_["ETAMU"]->Fill(mumu->momentum().Eta());
			if( quantemu == 1 ) 	
			{
				MuGENch[0]=mumu->charge();
				p_mu1.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
			} 
			else 
			{
				MuGENch[1]=mumu->charge();
				p_mu2.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
			}		 
		  }
		} //end if dauId=jpsiId
		
		bool dauPi = false;
		if (isBc && dauId == piplus) 		dauPi= true;
		if (isBp && fabs(dauId) == piplus)   dauPi= true;
		if( dauPi)  
		{
		  p_pi.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
		  PiGENch=des->charge();
		  boolpi = true;
		  h1_["ETAPI"]->Fill(p_pi.Eta());
		}	      
	  } // end for des
	
	  xprimBc = p.vertex().x();
	  PVMC[0]=xprimBc;
	  yprimBc = p.vertex().y();
	  PVMC[1]=yprimBc;
	  zprimBc = p.vertex().z();
	  PVMC[2]=zprimBc;

	  double dx = (xsecBc - xprimBc);
	  double dy = (ysecBc - yprimBc);
	  double dz = (zsecBc - zprimBc);
	  elleBc = sqrt( dx*dx + dy*dy + dz*dz );      
	  h1_["XPRIMBC" ]->Fill(xprimBc);
	  h1_["YPRIMBC" ]->Fill(yprimBc);
	  h1_["ZPRIMBC" ]->Fill(zprimBc);
	  h1_["XSECBC"  ]->Fill(xsecBc);
	  h1_["YSECBC"  ]->Fill(ysecBc);
	  h1_["ZSECBC"  ]->Fill(zsecBc);
	  h1_["ELLEBC"  ]->Fill(elleBc);
	  p_Bc.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());
	  TAU=elleBc*BcMassMC/c_const/sqrt((p.px()*p.px())+ (p.py()*p.py())+ (p.pz()*p.pz()));
	  h1_["TAUBC"]->Fill(TAU);

	  double cosMC =((xsecBc - xprimBc)*p.px()+ (ysecBc - yprimBc)*p.py()+ (zsecBc - zprimBc)*p.pz())/
					(elleBc*sqrt(p.px()*p.px()+ p.py()*p.py()+ p.pz()*p.pz()));
	 
	  h1_["COSENO_MC2"]->Fill(cosMC);	    
	  
	  if(boolJpsi && boolpi) boolMC = true;
    } // end MonteCarlo loop

    double ptotmu1 = sqrt(p_mu1.Px()*p_mu1.Px()+p_mu1.Py()*p_mu1.Py()+p_mu1.Pz()*p_mu1.Pz());
    double ptotmu2 = sqrt(p_mu2.Px()*p_mu2.Px()+p_mu2.Py()*p_mu2.Py()+p_mu2.Pz()*p_mu2.Pz());
    double ptotpi  = sqrt(p_pi.Px()*p_pi.Px()+p_pi.Py()*p_pi.Py()+p_pi.Pz()*p_pi.Pz());
    p_jpsi         = p_mu1+p_mu2;

    h1_["PTOTMU1"]	  -> Fill(ptotmu1);
    h1_["PTOTMU2"]	  -> Fill(ptotmu2);
    h1_["PTOTPI" ]	  -> Fill(ptotpi);
    h1_["PTMU1"  ]	  -> Fill(p_mu1.Perp());
    h1_["PTMU2"  ]	  -> Fill(p_mu2.Perp());
    h1_["PTPI"   ]	  -> Fill(p_pi.Perp());
    h2_["PTOTMU1MU2"] -> Fill(ptotmu1,ptotmu2);
    h2_["PTMU1MU2"]   -> Fill(p_mu1.Perp(),p_mu2.Perp());			    
    
    COSENO =   ((xsecBc - xprimBc) * (p_mu1.Px()+p_mu2.Px()+p_pi.Px())+
	      		(ysecBc - yprimBc) * (p_mu1.Py()+p_mu2.Py()+p_pi.Py())+
	      		(zsecBc - zprimBc) * (p_mu1.Pz()+p_mu2.Pz()+p_pi.Pz()))/
	      	   ((elleBc)*
	   	   sqrt((p_mu1.Px()+p_mu2.Px()+p_pi.Px())*(p_mu1.Px()+p_mu2.Px()+p_pi.Px())+
	   			(p_mu1.Py()+p_mu2.Py()+p_pi.Py())*(p_mu1.Py()+p_mu2.Py()+p_pi.Py())+
	   			(p_mu1.Pz()+p_mu2.Pz()+p_pi.Pz())*(p_mu1.Pz()+p_mu2.Pz()+p_pi.Pz())));
    
    h1_["COSENO_MC"]->Fill(COSENO);			    

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
pair<double,double> Bc2JpsiPi::pionImpactParameter(reco::TransientTrack piTT, TransientVertex jpsiVtx)
{
    pair<double,double> measure;
	pair<bool,Measurement1D>  piIP_pair = IPTools::absoluteImpactParameter3D(piTT, jpsiVtx);
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
//============================================================================================== 
pair<double,double> Bc2JpsiPi::pionIPBeamSpot(reco::TransientTrack piTT, GlobalPoint BsGp)
{
    pair<double,double> measureBS;
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
