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
#include <math.h>
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"   
#include "HepMC/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" 
#include <DataFormats/Math/interface/deltaR.h>
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TMath.h"

#include "HeavyFlavorAnalysis/Bc2Jpsi3Pi/interface/Bc2Jpsi3Pi.h"

// constants, enums and typedefs
double c_const  = 0.0299792458;//cm/ps
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
TLorentzVector p_mu1,p_mu2, p_pi1, p_pi2, p_pi3, p_jpsi,p_Bc;
double dRGEN[3] = {0,0,0};
double GENres[5]= {0,0,0,0,0};
double PVMC[3]  = {0,0,0};
double SVMC[3]  = {0,0,0};
int PiGENch[3]  = {0,0,0};
int Mu1GENch    = 0;
int Mu2GENch    = 0;
bool boolMC     = false;


using namespace std ;

Bc2Jpsi3Pi::Bc2Jpsi3Pi(const edm::ParameterSet& iConfig) :
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
Bc2Jpsi3Pi::~Bc2Jpsi3Pi()
{
//    for(std::map<std::string, TH1F*>::iterator it=h_.begin(); it!=h_.end(); ++it)
//    {
//     if( it->second ) delete it->second ;
//    }
//    if( file ) delete file ;
}

//=============== method called to for each event ===================================================================
void Bc2Jpsi3Pi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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

   int MuGENch[2] ;
   MuGENch[0] = Mu1GENch; 
   MuGENch[1] = Mu2GENch; 
   
   BcTree_->BcTree::SetBcTreeGENCand(  
					 p_Bc 		    ,
					 p_mu1 		    ,
					 p_mu2 		    ,
					 p_jpsi		    ,
					 p_pi1 		    ,
					 p_pi2 		    ,
					 p_pi3 		    ,
					 elleBc         ,
					 COSENO         ,
					 BeamSpot       ,
					 SVMC           ,
					 GENres			,
					 PiGENch        ,
					 MuGENch
	                 );

   int counter = 0 ;
   h1_["filter"]->Fill(counter);//0

   const reco::Muon* recomu1   = 0;
   const reco::Muon* recomu2   = 0;
 
   typedef ROOT::Math::SVector<double, 3>    Vector3;
   typedef ROOT::Math::SMatrix<double, 3, 3> Matrix33;
   Matrix33 Covsec;
   Matrix33 Covprim;
   Vector3  Deriv;
 	

   int      NumBc       	= 0;// remove header from list
   double   qBestRis[5] 	={0,0,0,0,0};   
   double	bestVtxS[3] 	={0,0,0};
   double	bestIP3D[3] 	= {0,0,0};
   double	bestIPbs[3] 	= {0,0,0};
   double	bestsIP3D[3]	= {0,0,0};
   double	bestsIPbs[3]	= {0,0,0};
   int      bestPiCh[3] 	= {0,0,0};
   int 		MuCh[2]     	= {0,0};
   int      trkNHits[3] 	= {0,0,0};
   int      trkPixelHits[3]	= {0,0,0};
   float	trkChi2[3] 		= {0,0,0};
   int      triggerBit  	= 0;

   reco::Candidate::LorentzVector qCandPsi1   ;
   reco::Candidate::LorentzVector qCandPsi2   ;
   reco::Candidate::LorentzVector qCanda1     ;
   reco::Candidate::LorentzVector qCandRho1   ;
   reco::Candidate::LorentzVector qCandRho2   ;

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

   edm::Handle<reco::MuonCollection>   recoMuons    ; 
   edm::Handle<reco::TrackCollection>  muonTracks   ;
   iEvent.getByLabel("muons",          recoMuons )	;
   iEvent.getByLabel("globalMuons",    muonTracks)	;

   edm::ESHandle<GlobalTrackingGeometry>          theTrackingGeometry;
   edm::ESHandle<TransientTrackBuilder>           theB;
   edm::Handle<reco::TrackCollection>             tracks ;
   edm::Handle<reco::TrackCollection>             pvtracks;   
   edm::Handle<reco::VertexCollection>            priVtxs;
   edm::Handle<reco::BeamSpot>	                  pvbeamspot; 

   iEvent.getByLabel(pString_["trackCollection"], tracks	);
   iEvent.getByLabel(thePVs_,                     priVtxs	);

   iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

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
   
   edm::Handle<edm::TriggerResults> hltresults;
   edm::InputTag tag("TriggerResults","","HLT");
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
  //    std::cout << triggerNames_.triggerName(itrig) << std::endl;
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
  
   //pion tracks selection 
  newTracksDef qualityTracks;
  for (uint tracksIt =0 ;  tracksIt < tracks->size(); tracksIt++)
  {
     reco::TrackRef checkTrk(tracks,tracksIt) ;										        
     bool flag = false ; 												        
     for (reco::MuonCollection::const_iterator mu =recoMuons->begin(); mu!=recoMuons->end(); mu++)			        
     {														        
       if (mu->track().get()!= 0 && mu->track().get() == checkTrk.get()) { flag=true; break;}			        
     }														        
     if (flag)												   											continue;	        
     if (checkTrk->pt()  				  					< pDouble_["cut_Pt_Trk"])			   		continue;	        
     if (fabs(checkTrk->eta())			  					> pDouble_["cut_eta_pi"])			  		continue;	   						  
     if (checkTrk->numberOfValidHits()	  					< pDouble_["cut_nhits"] )			   		continue;	   					   
     if (checkTrk->hitPattern().numberOfValidPixelHits() 	< pDouble_["numberOfValidPixelHitsCut"]  )	continue;	        
     if (checkTrk->normalizedChi2()			 			    > pDouble_["cut_chi2n"])			   		continue; 	        
     qualityTracks[tracksIt]=checkTrk;											   	
  }
     
  
  BcTree_->BcTree::SetBcTreeHeader(  nEvents_		       ,
 				     (int32_t)iEvent.id().run()            ,
 				     (int32_t)iEvent.id().luminosityBlock(),
 				     (int32_t)iEvent.id().event()	       ,
 				      NumBc		                           ,
 				     (int32_t)tracks ->size()              ,
 				     (int32_t)priVtxs->size()              ,
 					 Tnpv       						   ,
 					 triggerBit                            , 
 					 qualityTracks.size()
                     );

   if ( !triggerAccepted ) {
     if (boolMC) t1_["BcTree"]->Fill();
     return ;
   }
  
   this->checkHLT(iEvent) ;  
   h1_["filter"      ]->Fill(counter++);//1
   if (recoMuons->size()<2) 
   {
     if (boolMC) t1_["BcTree"]->Fill();
     return ;
   }
 
   const reco::VertexCollection primvtx = *(priVtxs.product());    
   h1_["pVertexSize"        ]-> Fill(primvtx.size()			);
   h1_["pVertexTrackSize"	]-> Fill(primvtx[0].tracksSize());
   h1_["EventTrackSize"		]-> Fill(tracks->size() 		);
   h1_["goodTrkSize"		]-> Fill(qualityTracks.size() 	);
   if(qualityTracks.size() < 3 )
   {
     if (boolMC) t1_["BcTree"]->Fill();
     return; 
   }

  //#################################################################################################################################
   int JpsiNoTrig = 0; 
   int JpsiTrig = 0; 

   for(reco::MuonCollection::const_iterator muR1 =recoMuons->begin(); muR1!=recoMuons->end(); ++muR1)
   {
    recomu1 = &*muR1;
    reco::TrackRef refmu1t=recomu1->track(); 
     
    if (fabs(muR1->eta())>pDouble_["cut_eta"]) 							         						continue;
    if (fabs(muR1->pt() )<pDouble_["cut_Pt_Mu"])							         					continue;

    //SOFT Muon Selection for Mu1
    if (!muon::isGoodMuon(*recomu1,muon::TMOneStationTight)                    )				 		continue;
    if (recomu1->innerTrack()->hitPattern().numberOfValidTrackerHits()   <= 10 )						continue;
    if (recomu1->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 1  )						continue;
    if (recomu1->innerTrack()->normalizedChi2() 						 >= 1.8)	 					continue;
    if (!(fabs(recomu1->innerTrack()->dxy(primvtx[0].position())) < 3. && 
          fabs(recomu1->innerTrack()->dz(primvtx[0].position()))  < 30. )) 								continue;

    for(reco::MuonCollection::const_iterator muR2 =muR1+1; muR2!=recoMuons->end(); ++muR2)
    {
        recomu2 = &*muR2;
        reco::TrackRef refmu2t=recomu2->track();

        if (muR1->charge()+muR2->charge()!=0)									   						continue; 
        MuCh[0] = muR1->charge();
        MuCh[1] = muR2->charge();
        if (fabs(muR2->eta())>pDouble_["cut_eta"])								    					continue;
        if (fabs(muR2->pt() )<pDouble_["cut_Pt_Mu"])													continue;

        //SOFT Muon Selection for Mu2
    	if (!muon::isGoodMuon(*recomu2,muon::TMOneStationTight)                    )				 	continue;
    	if (recomu2->innerTrack()->hitPattern().numberOfValidTrackerHits()   <= 10 )					continue;
    	if (recomu2->innerTrack()->hitPattern().pixelLayersWithMeasurement() <= 1  )					continue;
    	if (recomu2->innerTrack()->normalizedChi2() 						 >= 1.8)	 				continue;
    	if (!(fabs(recomu2->innerTrack()->dxy(primvtx[0].position())) < 3. && 
              fabs(recomu2->innerTrack()->dz(primvtx[0].position()))  < 30. )) 							continue;
     

  	    if( !recomu1->track().isNonnull() ||   !recomu2->track().isNonnull())							assert(0) ;

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

        h1_["filter"      ]->Fill(counter++);//2

        reco::Candidate::LorentzVector jpsV = tr1p4 + tr2p4;
        double JPsiMass = jpsV.M() ;
        double ptJpsi   = jpsV.pt();


        vector<reco::TransientTrack> tTrMu;
        tTrMu.push_back(tTrackmu1);										  
        tTrMu.push_back(tTrackmu2);

        KalmanVertexFitter fitter1;
        TransientVertex    mumuVertex;
        mumuVertex = fitter1.vertex(tTrMu);
        if(!mumuVertex.isValid()) 																	continue;
        h1_["JpsiPt_noTrigger"                ]->Fill(jpsV.pt());
        if(triggerBit) h1_["JpsiPt_TriggerBit"]->Fill(jpsV.pt());

        h1_["filter"      ]->Fill(counter++);//3
        h1_["InvMassJPsi" ]->Fill(jpsV.M());

        math::XYZVector pperp(refmu1t->px() + refmu2t->px(), refmu1t->py() + refmu2t->py(), 0.);

        //Trigger requirements - DCA
	    TrajectoryStateClosestToPoint mu1TS = tTrackmu1.impactPointTSCP();
	    TrajectoryStateClosestToPoint mu2TS = tTrackmu2.impactPointTSCP();
	    if (mu1TS.isValid() && mu2TS.isValid()) 
	    {
          ClosestApproachInRPhi cApp;
          cApp.calculate(mu1TS.theState(), mu2TS.theState());
          if (!cApp.status() || cApp.distance() > 0.5) continue;
	    }

        GlobalPoint JpsiVertex = mumuVertex.position();
        GlobalError JpsiVertexError = mumuVertex.positionError();
        GlobalPoint displacementFromBeamspot(-1*((bs.x0() - mumuVertex.position().x()) + (mumuVertex.position().z() - bs.z0()) * bs.dxdz()),
 	          			     				 -1*((bs.y0() - mumuVertex.position().y()) + (mumuVertex.position().z() - bs.z0()) * bs.dydz()), 0);
       
        double LPromptsecXY = displacementFromBeamspot.perp();
        double lxyerr 		= sqrt(JpsiVertexError.rerr(displacementFromBeamspot));
        reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
        double cosJpsiXY 	= vperp.Dot(pperp)/(vperp.R()*pperp.R());
    	double elsigJpsi    = 0 ;

        try	    { elsigJpsi= LPromptsecXY / lxyerr ; }
		catch (...) {cout << __LINE__ << " " << __PRETTY_FUNCTION__ << "\t" << "Floating divide by zero!" << endl ;}
  	
	    h1_["elsigJpsi"]-> Fill(elsigJpsi); 
	    h1_["JpsiPt"   ]-> Fill(ptJpsi); 

        //Trigger requirements
        float vProbJpsi(TMath::Prob(mumuVertex.totalChiSquared(),(int)(mumuVertex.degreesOfFreedom())));
        if (vProbJpsi <  pDouble_["cut_cl_Jpsi"])            										continue;
        h1_["filter"      ]->Fill(counter++);
        h1_["InvMassJPsicl"         ]->Fill(jpsV.M());

        if (elsigJpsi <= 3)    																		continue;
        h1_["filter"      ]->Fill(counter++);

        if (cosJpsiXY <= 0.9)   																	continue;
        h1_["filter"      ]->Fill(counter++);
       
        if (ptJpsi <=  6.9)																			continue;
        h1_["filter"      ]->Fill(counter++);

        if( fabs(JPsiMass-pDouble_["JPsiMassPDG"])> pDouble_["JPsiMassWindow"]) 					continue;
        h1_["filter"      ]->Fill(counter++);

        //Trigger Matching
        JpsiNoTrig++ ; 
        double hltMatchDR =0.5; 															     	
        bool MuTrigMatch = false;															     	
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
    
		if(triggerBit) h1_["JpsiPt_TriggerBitAndTriggerCuts"]->Fill(jpsV.pt());

		h1_["JpsiBeforeTrigMatch"  		 ]-> Fill(JpsiNoTrig);       
		if( !MuTrigMatch ) 				 continue ;
		h1_["JpsiAfterTrigMatch"   		 ]-> Fill(JpsiTrig  );
		h1_["JpsiPtAfterTrigger"   		 ]-> Fill(ptJpsi	); 
		h1_["InvMassJpsiPassingTrigMatch"]-> Fill(jpsV.M()  );

	    TLorentzVector *jpsiMu1  = new TLorentzVector(tr1p4.px(), tr1p4.py(), tr1p4.pz(), tr1p4.energy());
	    TLorentzVector *jpsiMu2  = new TLorentzVector(tr2p4.px(), tr2p4.py(), tr2p4.pz(), tr2p4.energy());
	    TLorentzVector *jpsiTV   = new TLorentzVector(jpsV.px(),   jpsV.py(),  jpsV.pz(),  jpsV.energy());
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
						*jpsiTV			,
						elsigJpsi      	,
						cosJpsiXY      	,
						vProbJpsi      	,
						JpsiVtx        	,
						matchMuJpsi    	,
						MuCh			);

		//Loop on tracks
        unsigned int trk1I = -1 ;
        for (newTracksDef::const_iterator it_pi1=qualityTracks.begin(); it_pi1!=qualityTracks.end(); ++it_pi1)
        {
          trk1I=(*it_pi1).first ;
   	      reco::TrackRef firstPionTrackCand((*it_pi1).second) ;
          const reco::TransientTrack pi1TrackCand = (*theB).build(firstPionTrackCand);
          
          pair<double,double>  pion1_IPPair = pionImpactParameter(pi1TrackCand,mumuVertex);
          double pi1_3DIp                   = pion1_IPPair.first;
          double pi1_3DIpSignificance       = pion1_IPPair.second;
          h1_["Pion1_ImpactParameter"]->Fill(pi1_3DIp);
          h1_["Pion1_IPsignificance" ]->Fill(pi1_3DIpSignificance);
          
          GlobalPoint BeamSpotGP(bs.x0(), bs.y0(), bs.z0()); 

          pair<double,double>  pion1_BSPair = pionIPBeamSpot(pi1TrackCand,BeamSpotGP);
          double pi1_IPbs                   = pion1_BSPair.first;
          double pi1_IPbsSignificance       = pion1_BSPair.second;
	      h1_["Pion1_TransverseImpactParameter_BS" ]-> Fill(pi1_IPbs);
	      h1_["Pion1_TIParameterBSsignificance"    ]-> Fill(pi1_IPbsSignificance);

		  //Require small impact parameter significance wrt Jpsi vertex
          if(pi1_3DIpSignificance > 6) 									continue;

		  //Require trk in a cone wrt Jpsi track
          double DRPi1 =  deltaR(firstPionTrackCand->eta(),firstPionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
          h1_["deltaRPi"     ]-> Fill(DRPi1);
          if (DRPi1 > 2.5) continue;

		  //Loop on second pion
          unsigned int trk2I = -1;   
	      for (newTracksDef::const_iterator it_pi2=it_pi1; it_pi2!=qualityTracks.end(); ++it_pi2)
          {
	        if( (*it_pi1).first == (*it_pi2).first ) 					continue ;
	        trk2I = (*it_pi2).first;  ;
	        reco::TrackRef secondPionTrackCand((*it_pi2).second) ;
   	        const reco::TransientTrack pi2TrackCand = (*theB).build(secondPionTrackCand);

            pair<double,double>  pion2_IPPair = pionImpactParameter(pi2TrackCand,mumuVertex);
            double pi2_3DIp                   = pion2_IPPair.first;
            double pi2_3DIpSignificance       = pion2_IPPair.second;
			h1_["Pion2_ImpactParameter"]-> Fill(pi2_3DIp);
	        h1_["Pion2_IPsignificance" ]-> Fill(pi2_3DIpSignificance);
	    
            pair<double,double>  pion2_BSPair = pionIPBeamSpot(pi2TrackCand,BeamSpotGP);
            double pi2_IPbs                   = pion2_BSPair.first;
            double pi2_IPbsSignificance       = pion2_BSPair.second;
	        h1_["Pion2_TransverseImpactParameter_BS"]-> Fill(pi2_IPbs);
	        h1_["Pion2_TIParameterBSsignificance"   ]-> Fill(pi2_IPbsSignificance);

			//Require small impact parameter significance wrt Jpsi vertex
			if(pi2_3DIpSignificance > 6) 								continue;
			
		    //Require trk in a cone wrt Jpsi track
			double DRPi2 =  deltaR(secondPionTrackCand->eta(),secondPionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
			h1_["deltaRPi"     ]-> Fill(DRPi2);
			if (DRPi2 > 2.5) 				continue;

			//Loop on third pion
			unsigned int trk3I = -1 ;	   
			for (newTracksDef::const_iterator it_pi3=it_pi2; it_pi3!=qualityTracks.end(); ++it_pi3)
			{
			  if( (*it_pi3).first == (*it_pi2).first ) 					continue ;
			  trk3I = (*it_pi3).first;  ;
			  reco::TrackRef thirdPionTrackCand((*it_pi3).second) ;
			  const reco::TransientTrack pi3TrackCand = (*theB).build(thirdPionTrackCand);

			  pair<double,double>  pion3_IPPair = pionImpactParameter(pi3TrackCand,mumuVertex);
			  double pi3_3DIp                   = pion3_IPPair.first;
			  double pi3_3DIpSignificance       = pion3_IPPair.second;
			  h1_["Pion3_ImpactParameter"]-> Fill(pi3_3DIp);
			  h1_["Pion3_IPsignificance" ]-> Fill(pi3_3DIpSignificance);

			  pair<double,double>  pion3_BSPair = pionIPBeamSpot(pi3TrackCand,BeamSpotGP);
			  double pi3_IPbs                   = pion3_BSPair.first;
			  double pi3_IPbsSignificance       = pion3_BSPair.second;
			  h1_["Pion3_TransverseImpactParameter_BS"]-> Fill(pi3_IPbs);
			  h1_["Pion3_TIParameterBSsignificance"   ]-> Fill(pi3_IPbsSignificance);

			  //Require small impact parameter significance wrt Jpsi vertex
			  if(pi3_3DIpSignificance> 6) 									  continue;

			  //Require trk in a cone wrt Jpsi track
			  double DRPi3 =  deltaR(thirdPionTrackCand->eta(),thirdPionTrackCand->phi(),jpsiTV->Eta(),jpsiTV->Phi());
			  h1_["deltaRPi"     ]-> Fill(DRPi3);
			  if (DRPi3 > 2.5) continue;

			  double piMassSquare = pow(0.13957018,2); 
			  double energy1 = sqrt(firstPionTrackCand->px()*firstPionTrackCand->px()   + firstPionTrackCand->py()*firstPionTrackCand->py()   + firstPionTrackCand->pz()*firstPionTrackCand->pz()   + piMassSquare);
			  double energy2 = sqrt(secondPionTrackCand->px()*secondPionTrackCand->px() + secondPionTrackCand->py()*secondPionTrackCand->py() + secondPionTrackCand->pz()*secondPionTrackCand->pz() + piMassSquare);
			  double energy3 = sqrt(thirdPionTrackCand->px()*thirdPionTrackCand->px()   + thirdPionTrackCand->py()*thirdPionTrackCand->py()   + thirdPionTrackCand->pz()*thirdPionTrackCand->pz()   + piMassSquare);

			  TLorentzVector *origPi1 = new TLorentzVector(firstPionTrackCand->px(), firstPionTrackCand->py(), firstPionTrackCand->pz(), energy1);
			  TLorentzVector *origPi2 = new TLorentzVector(secondPionTrackCand->px(),secondPionTrackCand->py(),secondPionTrackCand->pz(),energy2);
			  TLorentzVector *origPi3 = new TLorentzVector(thirdPionTrackCand->px(), thirdPionTrackCand->py(), thirdPionTrackCand->pz(), energy3);

			  if(iEvent.isRealData())
			  {  	  
			   if (fabs(firstPionTrackCand ->charge()+
					 secondPionTrackCand->charge()+
					 thirdPionTrackCand ->charge()) == 1 )					  continue ;
			   }
			
			  if(!iEvent.isRealData())
			  {  	  
			   if (firstPionTrackCand ->charge()+
					 secondPionTrackCand->charge()+
					 thirdPionTrackCand ->charge() != 1 )					  continue ;
			  }
		
	
			  KinematicParticleFactoryFromTransientTrack pFactory;
			  const ParticleMass muon_mass 	= 0.10565837; 
			  float muon_sigma 			  	= muon_mass*1.e-6;
			  ParticleMass jpsimass_c 	  	= 3.096916; 
			  ParticleMass pion_mass  	  	= 0.13957018; 
			  float pion_sigma        	  	= pion_mass*1.e-6; 
			  float chi 		 			= 0.;
			  float ndf 		 			= 0.;

			  std::vector<RefCountedKinematicParticle> XParticles_Bc;
			  XParticles_Bc.push_back(pFactory.particle(tTrackmu1,muon_mass,chi,ndf,muon_sigma));
			  XParticles_Bc.push_back(pFactory.particle(tTrackmu2,muon_mass,chi,ndf,muon_sigma));
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
	 
			  double vProbBc = ChiSquaredProbability((double)(BcVertex->chiSquared()),(double)(BcVertex->degreesOfFreedom()));
			  h1_["Bcvtxcl"		      ]-> Fill(vProbBc);
			  if (vProbBc <= pDouble_["cut_cl_Bc"])                                                           continue ;

			  BcVertexFitTree->movePointerToTheFirstChild();
			  RefCountedKinematicParticle muP =  BcVertexFitTree->currentParticle();
			  BcVertexFitTree->movePointerToTheNextChild();
			  RefCountedKinematicParticle muN =  BcVertexFitTree->currentParticle();
			  BcVertexFitTree->movePointerToTheNextChild();
			  RefCountedKinematicParticle pi1 =  BcVertexFitTree->currentParticle();
			  BcVertexFitTree->movePointerToTheNextChild();
			  RefCountedKinematicParticle pi2 =  BcVertexFitTree->currentParticle();
			  BcVertexFitTree->movePointerToTheNextChild();
			  RefCountedKinematicParticle pi3 =  BcVertexFitTree->currentParticle();

			  reco::Candidate::LorentzVector mu1p4R( muP->currentState().globalMomentum().x(),muP->currentState().globalMomentum().y(), 
													 muP->currentState().globalMomentum().z(),muP->currentState().kinematicParameters().energy());
			  reco::Candidate::LorentzVector mu2p4R( muN->currentState().globalMomentum().x(),muN->currentState().globalMomentum().y(), 
													 muN->currentState().globalMomentum().z(),muN->currentState().kinematicParameters().energy());
			  reco::Candidate::LorentzVector pi1p4R( pi1->currentState().globalMomentum().x(),pi1->currentState().globalMomentum().y(), 
													 pi1->currentState().globalMomentum().z(),pi1->currentState().kinematicParameters().energy());
			  reco::Candidate::LorentzVector pi2p4R( pi2->currentState().globalMomentum().x(),pi2->currentState().globalMomentum().y(), 
													 pi2->currentState().globalMomentum().z(),pi2->currentState().kinematicParameters().energy());
			  reco::Candidate::LorentzVector pi3p4R( pi3->currentState().globalMomentum().x(),pi3->currentState().globalMomentum().y(), 
													 pi3->currentState().globalMomentum().z(),pi3->currentState().kinematicParameters().energy());
			  reco::Candidate::LorentzVector Bc_3Pi( BcCand->currentState().globalMomentum().x(),BcCand->currentState().globalMomentum().y(), 
													 BcCand->currentState().globalMomentum().z(),BcCand->currentState().kinematicParameters().energy());

			  VertexDistance3D vertTool;
			  Measurement1D BcJpsiDistance = vertTool.distance(BcVertex->vertexState(),mumuVertex.vertexState());
			  double BcJpsiSignificance    = BcJpsiDistance.significance();
			  double BcJpsiDistanceV       = BcJpsiDistance.value();
			  h1_["BcJpsiSignificance" ]-> Fill(BcJpsiSignificance);
			  h1_["BcJpsiDistance"     ]-> Fill(BcJpsiDistanceV);

			 //MC match Muons
			 double matchMu[2];
			 double matchPi1[2];
			 double matchPi2[2];
			 double matchPi3[2];

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

			  //MC match Pions
			  if(firstPionTrackCand->charge() == PiGENch[0])
			  {
				matchPi1[0] = deltaR(pi1p4R.eta(),pi1p4R.phi(),p_pi1.Eta(),p_pi1.Phi());
				matchPi1[1] = deltaR(pi1p4R.eta(),pi1p4R.phi(),p_pi2.Eta(),p_pi2.Phi());
			  }
			  else
			  {
				matchPi1[0] = deltaR(pi1p4R.eta(),pi1p4R.phi(),p_pi3.Eta(),p_pi3.Phi());
				matchPi1[1] = 0;
			  }
		
			  if(secondPionTrackCand->charge() == PiGENch[0])
			  {
				matchPi2[0] = deltaR(pi2p4R.eta(),pi2p4R.phi(),p_pi1.Eta(),p_pi1.Phi());
				matchPi2[1] = deltaR(pi2p4R.eta(),pi2p4R.phi(),p_pi2.Eta(),p_pi2.Phi());
			  }
			  else
			  {
				matchPi2[0] = deltaR(pi2p4R.eta(),pi2p4R.phi(),p_pi3.Eta(),p_pi3.Phi());
				matchPi2[1] = 0;
			  }

			  if(thirdPionTrackCand->charge() == PiGENch[0])
			  {
				matchPi3[0] = deltaR(pi3p4R.eta(),pi3p4R.phi(),p_pi1.Eta(),p_pi1.Phi());
				matchPi3[1] = deltaR(pi3p4R.eta(),pi3p4R.phi(),p_pi2.Eta(),p_pi2.Phi());
			  }
			  else
			  {
				matchPi3[0] = deltaR(pi3p4R.eta(),pi3p4R.phi(),p_pi3.Eta(),p_pi3.Phi());
				matchPi3[1] = 0;
			  }
			 }//end if MC

			 //L and L/sigma 		
			 math::XYZVector BcPerp(Bc_3Pi.px(), Bc_3Pi.py(), 0.);

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
				 Covprim(i,j) 	 = bs.covariance(i,j);  
				 SVCovariance[n] = Covsec(i,j);
				 BSCovariance[n] = Covprim(i,j);
				 n++;
			   }
			 }

			 bool flag1= false;    
			 bool flag2= false; 
			 bool flag3= false; 
			 if ((firstPionTrackCand->charge() + secondPionTrackCand->charge())==0)  flag1=true;
			 if ((firstPionTrackCand->charge() + thirdPionTrackCand->charge()) ==0)  flag2=true;
			 if ((secondPionTrackCand->charge() + thirdPionTrackCand->charge())==0)  flag3=true;
	   
			 if (flag1)
			 {
				qCandPsi1= pi1p4R+pi2p4R+jpsV;
				qCandRho1= pi1p4R+pi2p4R;
				if (flag2)		
				{
					qCandPsi2= pi1p4R+pi3p4R+jpsV;
					qCandRho2= pi1p4R+pi3p4R;
				}
				else  
				{
					qCandPsi2= pi3p4R+pi2p4R+jpsV;	
					qCandRho2= pi3p4R+pi2p4R;	
				}	
			}
			else    
			{
				qCandPsi1= pi1p4R+pi3p4R+jpsV;
				qCandPsi2= pi2p4R+pi3p4R+jpsV;
				qCandRho1= pi1p4R+pi3p4R;	       
				qCandRho2= pi2p4R+pi3p4R;
			}
	   
			qCanda1=pi1p4R+pi2p4R+pi3p4R;

			//ntupla.fill    
			TLorentzVector *recoBc  = new TLorentzVector(Bc_3Pi.px(),Bc_3Pi.py(),Bc_3Pi.pz(),Bc_3Pi.energy());
			TLorentzVector *recoMu1 = new TLorentzVector(mu1p4R.px(),mu1p4R.py(),mu1p4R.pz(),mu1p4R.energy());
			TLorentzVector *recoMu2 = new TLorentzVector(mu2p4R.px(),mu2p4R.py(),mu2p4R.pz(),mu2p4R.energy());
			TLorentzVector *recoPi1 = new TLorentzVector(pi1p4R.px(),pi1p4R.py(),pi1p4R.pz(),pi1p4R.energy());
			TLorentzVector *recoPi2 = new TLorentzVector(pi2p4R.px(),pi2p4R.py(),pi2p4R.pz(),pi2p4R.energy());
			TLorentzVector *recoPi3 = new TLorentzVector(pi3p4R.px(),pi3p4R.py(),pi3p4R.pz(),pi3p4R.energy());
			reco::Candidate::LorentzVector JJpsi = mu1p4R+mu2p4R;
			TLorentzVector *recoJpsi  = new TLorentzVector(JJpsi.px(),JJpsi.py(),JJpsi.pz(),JJpsi.energy());
			TLorentzVector *recoJpsiV = new TLorentzVector(jpsV.px(),jpsV.py(),jpsV.pz(),jpsV.energy());

			qBestRis[0]   	  = qCandPsi1.M();
			qBestRis[1]   	  = qCandPsi2.M();	
			qBestRis[2]   	  = qCandRho1.M();
			qBestRis[3]   	  = qCandRho2.M();	
			qBestRis[4]   	  = qCanda1.M()  ;

			bestVtxS[0]   	  = xsec ;
			bestVtxS[1]   	  = ysec ;
			bestVtxS[2]   	  = zsec ;
			bestIP3D[0]   	  = pi1_3DIp;
			bestIP3D[1]   	  = pi2_3DIp;
			bestIP3D[2]   	  = pi3_3DIp;
			bestIPbs[0]   	  = pi1_IPbs;
			bestIPbs[1]   	  = pi2_IPbs;
			bestIPbs[2]   	  = pi3_IPbs;
			bestsIP3D[0]  	  = pi1_3DIpSignificance;
			bestsIP3D[1]  	  = pi2_3DIpSignificance;
			bestsIP3D[2]  	  = pi3_3DIpSignificance;
			bestsIPbs[0]  	  = pi1_IPbsSignificance;
			bestsIPbs[1]  	  = pi2_IPbsSignificance;
			bestsIPbs[2]  	  = pi3_IPbsSignificance;

			bestPiCh[0]   	  = firstPionTrackCand  ->charge();
			bestPiCh[1]   	  = secondPionTrackCand ->charge();
			bestPiCh[2]   	  = thirdPionTrackCand  ->charge();
			trkNHits[0]   	  = firstPionTrackCand  ->numberOfValidHits();
			trkNHits[1]   	  = secondPionTrackCand ->numberOfValidHits();
			trkNHits[2]   	  = thirdPionTrackCand  ->numberOfValidHits();
			trkPixelHits[0]   = firstPionTrackCand  ->hitPattern().numberOfValidPixelHits();
			trkPixelHits[1]   = secondPionTrackCand ->hitPattern().numberOfValidPixelHits();
			trkPixelHits[2]   = thirdPionTrackCand  ->hitPattern().numberOfValidPixelHits();
			trkChi2[0] 		  = firstPionTrackCand  ->normalizedChi2();
			trkChi2[1] 		  = secondPionTrackCand ->normalizedChi2();
			trkChi2[2] 		  = thirdPionTrackCand  ->normalizedChi2();

			NumBc++;
			//BcTree_->fBcTreeHeader.SetNumBc(NumBc);
			BcTree_->BcTree::AddBcTreeCand(  
							 *recoBc  	  			,
							 *recoMu1 	  			,
							 *recoMu2 	  			,
							 *recoJpsi	  			,
							 *recoPi1 	  			,
							 *recoPi2 	  			,
							 *recoPi3 	  			,
							 *origPi1 	  			,
							 *origPi2 	  			,
							 *origPi3 	  			,
							 *recoJpsiV   			,
							 qBestRis     			,
							 LBcXY 					, 
							 elsig    	  			,
							 CosBcXY   	  			,
							 vProbBc      			,
							 vProbJpsi    			,
							 BeamSpot     			,
							 bestVtxS     			,
							 BSCovariance 			,
							 SVCovariance 			,
							 bestIP3D     			,
							 bestIPbs     			,
							 bestsIP3D    			,
							 bestsIPbs    			,
							 matchMu      			,
							 matchPi1     			,
							 matchPi2     			,
							 matchPi3     			,
							 bestPiCh     			,
							 MuCh         			,
							 elsigJpsi    			,
							 cosJpsiXY  			,
							 trkNHits				,
							 trkPixelHits			,
							 trkChi2				,
							 hltMatch				);
          }//pi3
         }//pi2
       }//pi1 
     }//mu2			      
   }//mu1

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
    
    h1_["elsigJpsi"		    		] = new TH1F("elsigJpsi",   "elsigJpsi",80,0,20);
    h1_["JpsiPt"		    		] = new TH1F("JpsiPt",      "JpsiPt",100,0,100);
    h1_["JpsiPtAfterTrigger"  		] = new TH1F("JpsiPtAfterTrigger", "JpsiPtAfterTrigger",100,0,100);
    h1_["InvMassJPsi"               ] = new TH1F("InvMassJPsi",		   	"Invariant Mass JPsi",              500 , 0, 10 );
    h1_["InvMassJPsicl" 	    	] = new TH1F("InvMassJPsicl", "Invariant Mass JPsi cl",500,0,10); 
    h1_["JpsiBeforeTrigMatch"       ] = new TH1F("JpsiBeforeTrigMatch","Jpsi Before TrigMatch", 10, 0, 10) ;
    h1_["JpsiAfterTrigMatch"        ] = new TH1F("JpsiAfterTrigMatch","Jpsi After TrigMatch", 10, 0, 10) ;
    h1_["JpsiPt_noTrigger"    		] = new TH1F("JpsiPt_noTrigger", "JpsiPt_noTrigger",100,0,100);
    h1_["JpsiPt_TriggerBit"    		] = new TH1F("JpsiPt_TriggerBit","JpsiPt_TriggerBit",100,0,100);
    h1_["JpsiPt_TriggerBitAndTriggerCuts" ] = new TH1F("JpsiPt_TriggerBitAndTriggerCuts","JpsiPt_TriggerBitAndTriggerCuts",100,0,100);
    h1_["InvMassJpsiPassingTrigMatch" 	] = new TH1F("InvMassJpsiPassingTrigMatch", "InvMassJpsiPassingTrigMatch Mass JPsi cl",200,1,5); 
   
    file->cd();
 
//------Histograms for MC--------------------------------------------  CHECKED
   MonteCarloDir->cd();   
   h1_["PYID"        	  	        ] = new TH1F("PYID",	     "PYTHIA ID",2000,-1000,1000) ;		     
   h1_["EVTA"        		        ] = new TH1F("EVTA",	     "Size of the generated event",1500,0,1500) ;    
   h1_["ETAMU"        		        ] = new TH1F("ETAMU",	     "Eta distribution of generated mumu",100,-10,10) ;    
   h1_["ETAPI1"        		        ] = new TH1F("ETAPI1",	     "Eta distribution of generated pi1",100,-10,10) ;    
   h1_["ETAPI2"        		        ] = new TH1F("ETAPI2",	     "Eta distribution of generated pi2",100,-10,10) ;    
   h1_["ETAPI3"        		        ] = new TH1F("ETAPI3",	     "Eta distribution of generated pi3",100,-10,10) ;    
   h1_["CHPI1"        		        ] = new TH1F("CHPI1",	     "charge of pi1",10,-5,5) ;    
   h1_["CHPI2"        		        ] = new TH1F("CHPI2",	     "charge of pi2",10,-5,5) ;    
   h1_["CHPI3"        		        ] = new TH1F("CHPI3",	     "charge of pi3",10,-5,5) ;    
   h1_["CHMU1"        		        ] = new TH1F("CHMU1",	     "charge of MU1",10,-5,5) ;    
   h1_["CHMU2"        		        ] = new TH1F("CHMU2",	     "charge of MU2",10,-5,5) ;    

   h1_["PRIMVTXNUMBER"				] = new TH1F("PRIMVTXNUMBER","Number of generated primary vertex",50,0,50);
   h1_["GOODEVT"					] = new TH1F("GOODEVT",	     "Event with geometric acceptance in muons",1500,0,1500) ; 
   h1_["DAUID"       		        ] = new TH1F("DAUID",	     "Bc daughters ID",2000,-1000,1000) ;
   h1_["DAUJPSIID"   		        ] = new TH1F("DAUJPSIID",    "JPSI daughters ID",1000,-500,500) ;
   h1_["XPRIMBC"     		        ] = new TH1F("XPRIMBC"," X primary vertex Bc", 200, -1., 1.) ;
   h1_["YPRIMBC"     		        ] = new TH1F("YPRIMBC"," Y primary vertex Bc", 200, -1., 1.) ;
   h1_["ZPRIMBC"     		        ] = new TH1F("ZPRIMBC"," Z primary vertex Bc", 300, -15.0,  15.0)  ;
   h1_["XSECBC"      		        ] = new TH1F("XSECBC", " X secondary vertex Bc", 200, -1., 1.) ;
   h1_["YSECBC"      		        ] = new TH1F("YSECBC", " Y secondary vertex Bc", 200, -1., 1.) ;
   h1_["ZSECBC"      		        ] = new TH1F("ZSECBC", " Z secondary vertex Bc", 300, -15.0, 15.0) ;   
   h1_["ELLEBC"      		        ] = new TH1F("ELLEBC", " Bc decay lenght",	 400,  -2.,   2.0)  ;
   h1_["TAUBC"			        	] = new TH1F("TAUBC", " Bc decay time", 2000,  -10,   10)  ;
   h1_["DIFXPRIM"    		        ] = new TH1F("DIFXPRIM", "Xprim evt - Xprim Bc", 2000, -1., 1.) ;
   h1_["DIFYPRIM"    		        ] = new TH1F("DIFYPRIM", "Yprim evt - Yprim Bc", 2000, -1., 1.) ;
   h1_["DIFZPRIM"    		        ] = new TH1F("DIFZPRIM", "Zprim evt - Zprim Bc", 2000, -1., 1.) ;
   h1_["PTOTMU1"     		        ] = new TH1F("PTOTMU1",      "Momentum of 1st muon from Jpsi", 100, 0., 100. );
   h1_["PTOTMU2"     		        ] = new TH1F("PTOTMU2",      "Momentum of 2nd muon from Jpsi", 100, 0., 100. );
   h1_["PTOTPI1"      		        ] = new TH1F("PTOTPI1",       "Momentum of Pi1 from Bc+", 100, 0., 100. );
   h1_["PTOTPI2"      		        ] = new TH1F("PTOTPI2",       "Momentum of Pi2 from Bc+", 100, 0., 100. );
   h1_["PTOTPI3"      		        ] = new TH1F("PTOTPI3",       "Momentum of Pi3 from Bc+", 100, 0., 100. );
   h1_["PTJPSI"      		        ] = new TH1F("PTJPSI",       "Transverse Momentum of Jpsi", 400, 0., 40. );
   h1_["PTMU1"	     		        ] = new TH1F("PTMU1",	     "Transverse momentum of 1st muon from Jpsi", 250, 0., 25. ); 
   h1_["PTMU2"	     		        ] = new TH1F("PTMU2",	     "Transverse momentum of 2nd muon from Jpsi", 250, 0., 25. );
   h1_["PTPI1"	     		        ] = new TH1F("PTPI1",	     "Transverse momentum of Pi1 from Bc+", 250, 0., 25. );
   h1_["PTPI2"	     		        ] = new TH1F("PTPI2",	     "Transverse momentum of Pi2 from Bc+", 250, 0., 25. );
   h1_["PTPI3"	     		        ] = new TH1F("PTPI3",	     "Transverse momentum of Pi3 from Bc+", 250, 0., 25. );
   h1_["BCPT"                       ] = new TH1F("BCPT",          "Transverse momentum of Bc"    , 400, 0., 100. );
   h1_["COSENO_MC"   				] = new TH1F("COSENO_MC","COSENO_MC",  300, -1.5, 1.5) ;
   h1_["COSENO_MC2"  				] = new TH1F("COSENO_MC2","COSENO_MC2",  300, -1.5, 1.5) ;
   h2_["PTOTMU1MU2"  				] = new TH2F("PTOTMU1MU2",   "Mu2 momentum vs Mu1 momentum", 40, 0., 40., 40, 0., 40.);
   h2_["PTMU1MU2"    				] = new TH2F("PTMU1MU2", 	  "Mu2 pt vs Mu1 pt", 40, 0., 20., 40, 0., 20.);
   file->cd();				


//-----Histograms for Vertex Stuff-------------------------------------------- CHECKED
   VertexingDir->cd();   

   h1_["pVertexSize"					    ] = new TH1F("pVertexSize", 	"pVertexSize after trigger requirement ",		 20,0	  ,20	);
   h1_["pVertexTrackSize"				    ] = new TH1F("pVertexTrackSize","pVertexTrackSize", 	 800,0     ,800   );
   h1_["EventTrackSize" 				    ] = new TH1F("EventTrackSize",  "EventTrackSize",		 800,0     ,800   );
   h1_["goodTrkSize"        			    ] = new TH1F("goodTrkSize", 	 "number of tracks passing selections", 				800 , 0,800 );
   h1_["Bcvtxcl"		    			    ] = new TH1F("Bcvtxcl", 	 "Bc vertex CL before cut", 	 10000,0.	, 1.  );
   h1_["Pion1_ImpactParameter"		   		] = new TH1F("Pion1_ImpactParameter","Pion1_ImpactParameter",   10000,0   ,50   );
   h1_["Pion1_IPsignificance"              	] = new TH1F("Pion1_IPsignificance" ,"Pion1_IPsignificance" ,   10000,0   ,100  );
   h1_["Pion1_TransverseImpactParameter_BS"	] = new TH1F("Pion1_TransverseImpactParameter_BS","Pion1_TransverseImpactParameter_BS",  5000, 0,10   );
   h1_["Pion1_TIParameterBSsignificance"   	] = new TH1F("Pion1_TIParameterBSsignificance" ,"Pion1_TIParameterBSsignificance"     , 10000, 0,50   );
   h1_["Pion2_ImpactParameter"		   		] = new TH1F("Pion2_ImpactParameter","Pion2_ImpactParameter",   10000,0   ,50    );
   h1_["Pion2_IPsignificance"              	] = new TH1F("Pion2_IPsignificance" ,"Pion2_IPsignificance" ,   10000,0   ,100   );
   h1_["Pion2_TransverseImpactParameter_BS"	] = new TH1F("Pion2_TransverseImpactParameter_BS","Pion2_TransverseImpactParameter_BS",  5000, 0, 10  );
   h1_["Pion2_TIParameterBSsignificance"   	] = new TH1F("Pion2_TIParameterBSsignificance" ,"Pion2_TIParameterBSsignificance"     , 10000, 0, 50  );
   h1_["Pion3_ImpactParameter"		   		] = new TH1F("Pion3_ImpactParameter","Pion3_ImpactParameter",  10000,0   ,50   );
   h1_["Pion3_IPsignificance"              	] = new TH1F("Pion3_IPsignificance" ,"Pion3_IPsignificance" ,  10000,0   ,100  );
   h1_["Pion3_TransverseImpactParameter_BS"	] = new TH1F("Pion3_TransverseImpactParameter_BS","Pion3_TransverseImpactParameter_BS",  5000, 0 ,10 );
   h1_["Pion3_TIParameterBSsignificance"   	] = new TH1F("Pion3_TIParameterBSsignificance" ,"Pion3_TIParameterBSsignificance" ,     10000, 0, 50  );
   h1_["BcVtxFitTreeNonValid"   		    ] = new TH1F("BcVtxFitTreeNonValid",		"BcVtxFitTreeNonValid", 			4, 0,  4 );
   h1_["BcVtxNonValid"		    		    ] = new TH1F("BcVtxNonValid",		"BcVtxNonValid",				4, 0,  4 );
   h1_["BcJpsiSignificance"	    		    ] = new TH1F("BcJpsiSignificance",  	 "BcJpsiSignificance",  	 100, 0.  ,50.0 ) ;
   h1_["BcJpsiDistance"	        		    ] = new TH1F("BcJpsiDistance",  	 "BcJpsiDistance",  		 10000, 0.  ,5.0 ) ;
   h1_["deltaRPi"	            ] = new TH1F("deltaRPi",	   		 "deltaRPi",		   1000, 0.  ,10.0  ) ;

   file->cd();
   h1_["filter"      		      	] = new TH1F("filter",	   "Binned filter counter",100,0,100); // muon size
   h1_["TrueNumInteraction"      	] = new TH1F("TrueNumInteraction", "TrueNumInteraction",50,0,50); // muon size

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
   	int JpsiId   =  443;
   	int BcplusId =  541;
    int piplus   =  211;
    int piminus  = -211;
//    int PrimVtxNumber=0;
    double BcMassMC = 6.286;
    bool GoodEvt = true;
    boolMC = false;
	bool boolJpsi = false;
	bool boolpi1 = false;
	bool boolpi2 = false;
	bool boolpi3 = false;
    TLorentzVector GENCandPsi1, GENCandPsi2, GENCanda1, GENCandRho1, GENCandRho2;
	
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel("genParticles", genParticles);
	
/*	const HepMC::GenEvent * myGenEvent = evt->GetEvent();
   	h1_["EVTA"]->Fill(myGenEvent->particles_size());
	for ( 	HepMC::GenEvent::vertex_const_iterator vert = myGenEvent->vertices_begin();
 	 		vert != myGenEvent->vertices_end(); ++vert ) PrimVtxNumber++;
	h1_["PRIMVTXNUMBER"]->Fill(PrimVtxNumber);
    for ( HepMC::GenEvent::vertex_const_iterator vert = myGenEvent->vertices_begin(); vert != myGenEvent->vertices_end(); ++vert ) 
    {
    	xprimevt = (*vert)->position().x()/10.;
    	yprimevt = (*vert)->position().y()/10.;
    	zprimevt = (*vert)->position().z()/10.;
	    break;
    }*/
	
	int quanteBc=0;
   	elleBc = -1.;
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
	      		h1_["PTJPSI"]->Fill(des->momentum().Rho());		
				xsecBc  = des->vertex().x();
				SVMC[0] = xsecBc;
				ysecBc  = des->vertex().y();
				SVMC[1] = ysecBc;
		 	 	zsecBc  = des->vertex().z(); 
				SVMC[2] = zsecBc;
	       		boolJpsi = true;

		// return all daughters of J/psi    
		 		int quantemu = 0.;
	   		    for ( size_t imumu=0; imumu < des->numberOfDaughters(); ++imumu ) 
	            {
					const reco::Candidate *mumu = des->daughter(imumu);
					quantemu++;
					h1_["DAUJPSIID"]->Fill(mumu->pdgId());    		
					h1_["ETAMU"]->Fill(mumu->momentum().Eta());
					if ( fabs(mumu->momentum().Eta()) > 2.5 ) GoodEvt=false;
					if( quantemu == 1 ) 	
					{
					    Mu1GENch=mumu->charge();
					    p_mu1.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
					} 
					else 
					{
					    Mu2GENch=mumu->charge();
					    p_mu2.SetPxPyPzE(mumu->px(),mumu->py(),mumu->pz(),mumu->energy());
					}		 
		 		}
	       	} //end if dauId=jpsiId
	      
			if( dauId == piplus && quantipiplus == 0) 
			{
				quantipiplus++;
			    p_pi1.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
				PiGENch[0]=des->charge();
	       		boolpi1 = true;
				h1_["ETAPI1"]->Fill(p_pi1.Eta());
				h1_["CHPI1" ]->Fill(des->charge());
			}	      

			else if( dauId == piplus && quantipiplus == 1 ) 
			{
				quantipiplus++; 
			    p_pi2.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
				PiGENch[1]=des->charge();
	       		boolpi2 = true;
				h1_["ETAPI2"]->Fill(p_pi2.Eta());
				h1_["CHPI2" ]->Fill(des->charge());
			}	      
	  
			if( dauId == piminus && quantipiminus == 0) 	
			{
				quantipiminus++;
			    p_pi3.SetPxPyPzE(des->px(),des->py(),des->pz(),des->energy());
				PiGENch[2] = des->charge();
	       		boolpi3    = true;
				h1_["ETAPI3"]->Fill(p_pi3.Eta());
				h1_["CHPI3" ]->Fill(des->charge());
			} 	           
	    } // end for des
      
		if (GoodEvt) h1_["GOODEVT"]->Fill(1);      

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
        h1_["BCPT"    ]->Fill(p.pt());
	    p_Bc.SetPxPyPzE(p.px(),p.py(),p.pz(),p.energy());
        TAU=elleBc*BcMassMC/c_const/
				   sqrt((p.px()*p.px())+ (p.py()*p.py())+ (p.pz()*p.pz()));
        h1_["TAUBC"]->Fill(TAU);

        double cosMC =((xsecBc - xprimBc)*p.px()+ (ysecBc - yprimBc)*p.py()+ (zsecBc - zprimBc)*p.pz())/
	             	  (elleBc*sqrt(p.px()*p.px()+ p.py()*p.py()+ p.pz()*p.pz()));
    	h1_["COSENO_MC2"]->Fill(cosMC);	    
	    
        if(boolJpsi && boolpi1 && boolpi2 && boolpi3) boolMC = true;
    } // end genParticle loop

    double ptotmu1 = sqrt(p_mu1.Px()*p_mu1.Px()+p_mu1.Py()*p_mu1.Py()+p_mu1.Pz()*p_mu1.Pz());
    double ptotmu2 = sqrt(p_mu2.Px()*p_mu2.Px()+p_mu2.Py()*p_mu2.Py()+p_mu2.Pz()*p_mu2.Pz());
    double ptotpi1 = sqrt(p_pi1.Px()*p_pi1.Px()+p_pi1.Py()*p_pi1.Py()+p_pi1.Pz()*p_pi1.Pz());
    double ptotpi2 = sqrt(p_pi2.Px()*p_pi2.Px()+p_pi2.Py()*p_pi2.Py()+p_pi2.Pz()*p_pi2.Pz());
    double ptotpi3 = sqrt(p_pi3.Px()*p_pi3.Px()+p_pi3.Py()*p_pi3.Py()+p_pi3.Pz()*p_pi3.Pz());
    p_jpsi = p_mu1+p_mu2;

    h1_["PTOTMU1"]->Fill(ptotmu1);
    h1_["PTOTMU2"]->Fill(ptotmu2);
    h1_["PTOTPI1"]->Fill(ptotpi1);
    h1_["PTOTPI2"]->Fill(ptotpi2);
    h1_["PTOTPI3"]->Fill(ptotpi3);    
    h1_["PTMU1"]->Fill(p_mu1.Perp());
    h1_["PTMU2"]->Fill(p_mu2.Perp());
    h1_["PTPI1"]->Fill(p_pi1.Perp());
    h1_["PTPI2"]->Fill(p_pi2.Perp());
    h1_["PTPI3"]->Fill(p_pi3.Perp());
    h2_["PTOTMU1MU2"]->Fill(ptotmu1,ptotmu2);
    h2_["PTMU1MU2"]->Fill(p_mu1.Perp(),p_mu2.Perp());			    
    
    COSENO =   ((xsecBc - xprimBc) * (p_mu1.Px()+p_mu2.Px()+p_pi1.Px()+p_pi2.Px()+p_pi3.Px())+
	      		(ysecBc - yprimBc) * (p_mu1.Py()+p_mu2.Py()+p_pi1.Py()+p_pi2.Py()+p_pi3.Py())+
	      		(zsecBc - zprimBc) * (p_mu1.Pz()+p_mu2.Pz()+p_pi1.Pz()+p_pi2.Pz()+p_pi3.Pz()))/
	      	   ((elleBc)*
	   	   sqrt((p_mu1.Px()+p_mu2.Px()+p_pi1.Px()+p_pi2.Px()+p_pi3.Px())*(p_mu1.Px()+p_mu2.Px()+p_pi1.Px()+p_pi2.Px()+p_pi3.Px())+
	   			(p_mu1.Py()+p_mu2.Py()+p_pi1.Py()+p_pi2.Py()+p_pi3.Py())*(p_mu1.Py()+p_mu2.Py()+p_pi1.Py()+p_pi2.Py()+p_pi3.Py())+
	   			(p_mu1.Pz()+p_mu2.Pz()+p_pi1.Pz()+p_pi2.Pz()+p_pi3.Pz())*(p_mu1.Pz()+p_mu2.Pz()+p_pi1.Pz()+p_pi2.Pz()+p_pi3.Pz())));
    
    h1_["COSENO_MC"]->Fill(COSENO);			    

	dRGEN[0] = deltaR(p_jpsi.Eta(),p_jpsi.Phi(),p_pi1.Eta(),p_pi1.Phi());
	dRGEN[1] = deltaR(p_jpsi.Eta(),p_jpsi.Phi(),p_pi2.Eta(),p_pi2.Phi());
	dRGEN[2] = deltaR(p_jpsi.Eta(),p_jpsi.Phi(),p_pi3.Eta(),p_pi3.Phi());

//resonances-GEN------------------------------------------------------------------
	GENCandPsi1 = p_pi1+p_pi3+p_jpsi;
	GENCandPsi2 = p_pi2+p_pi3+p_jpsi;
	GENCandRho1 = p_pi1+p_pi3;
	GENCandRho2 = p_pi2+p_pi3;	
	GENCanda1   = p_pi1+p_pi2+p_pi3;

	GENres[0] = GENCandPsi1.M();
	GENres[1] = GENCandPsi2.M();
	GENres[2] = GENCandRho1.M();
	GENres[3] = GENCandRho2.M();
	GENres[4] = GENCanda1.M();
}

//=============================================================================================================================================== 
/*bool Bc2Jpsi3Pi::passPioncuts(reco::TrackRef pionCand)
{
	 if (pionCand->pt()<pDouble_["cut_Pt_Trk"])								 return false ;
     if (pionCand->p()<pDouble_["pCandMomCut"])								 return false ;
     if (fabs(pionCand->eta())>pDouble_["cut_eta_pi"])						 return false ;
	 int nhits=pionCand->numberOfValidHits();
  	 if( nhits <pDouble_["cut_nhits"] )										 return false ;
 	 if (pionCand->hitPattern().numberOfValidPixelHits() < pDouble_["numberOfValidPixelHitsCut"]  )	         return false ;
     if( pionCand->normalizedChi2() > pDouble_["cut_chi2n"])				 return false ;
         
	 return true;
}*/
//============================================================================================== 
pair<double,double> Bc2Jpsi3Pi::pionImpactParameter(reco::TransientTrack piTT, TransientVertex jpsiVtx)
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
pair<double,double> Bc2Jpsi3Pi::pionIPBeamSpot(reco::TransientTrack piTT, GlobalPoint BsGp)
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
