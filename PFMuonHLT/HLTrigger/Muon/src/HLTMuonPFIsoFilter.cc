/** \class HLTMuonIsoFilter
 *
 * See header file for documentation
 *
 *  \author J. Alcaraz
 *
 */   

#include "HLTrigger/Muon/interface/HLTMuonPFIsoFilter.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/HLTReco/interface/TriggerRefsCollections.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
// #include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

// #include "RecoMuon/MuonIsolation/interface/MuonIsolatorFactory.h"

#include <iostream>
//
// constructors and destructor
//
HLTMuonPFIsoFilter::HLTMuonPFIsoFilter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig),
   candTag_ (iConfig.getParameter< edm::InputTag > ("CandTag") ),
   previousCandTag_ (iConfig.getParameter<edm::InputTag > ("PreviousCandTag")),
   depTag_  (iConfig.getParameter< std::vector< edm::InputTag > >("DepTag" ) ),
   maxIso_  (iConfig.getParameter<double>("MaxIso" ) ),
//    theDepositIsolator(0),
   min_N_   (iConfig.getParameter<int> ("MinN")),
   doDeltaB_(iConfig.getParameter<bool> ("doDeltaBeta"))
{
  std::stringstream tags;
  for (unsigned int i=0;i!=depTag_.size();++i)
    tags<<" PFIsoTag["<<i<<"] : "<<depTag_[i].encode()<<" \n";
  LogDebug("HLTMuonPFIsoFilter") << " candTag : " << candTag_.encode()
				<< "\n" << tags 
				<< "  MinN : " << min_N_;

//    edm::ParameterSet isolatorPSet = iConfig.getParameter<edm::ParameterSet>("IsolatorPSet");
//    if (isolatorPSet.empty()) {
//      theDepositIsolator=0;
//        }else{
//      std::string type = isolatorPSet.getParameter<std::string>("ComponentName");
//      theDepositIsolator = MuonIsolatorFactory::get()->create(type, isolatorPSet);
//    }
   
//    if (theDepositIsolator) 
     produces<edm::ValueMap<bool> >();
}

HLTMuonPFIsoFilter::~HLTMuonPFIsoFilter()
{
}

//
// member functions
//

// ------------ method called to produce the data  ------------
 bool
 HLTMuonPFIsoFilter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup, trigger::TriggerFilterObjectWithRefs & filterproduct)
 {
    using namespace std;
    using namespace edm;
    using namespace trigger;
    using namespace reco;
 
    // All HLT filters must create and fill an HLT filter object,
    // recording any reconstructed physics objects satisfying (or not)
    // this HLT filter, and place it in the Event.
 
    //the decision map
    std::auto_ptr<edm::ValueMap<bool> > PFisoMap( new edm::ValueMap<bool> ());
 
    // get hold of trks
    Handle<RecoChargedCandidateCollection> mucands;
    if (saveTags()) filterproduct.addCollectionTag(candTag_);
    iEvent.getByLabel (candTag_,mucands);
    Handle<TriggerFilterObjectWithRefs> previousLevelCands;
    iEvent.getByLabel (previousCandTag_,previousLevelCands);
    vector<RecoChargedCandidateRef> vcands;
    previousLevelCands->getObjects(TriggerMuon,vcands);
    
    //get hold of energy deposition
    unsigned int nDep=depTag_.size();
    std::vector< Handle<edm::ValueMap<double> > > depMap(nDep);
    Handle<edm::ValueMap<bool> > decisionMap;
 
    for (unsigned int i=0;i!=nDep;++i) iEvent.getByLabel (depTag_[i],depMap[i]);

    // look at all mucands,  check cuts and add to filter object
    int nIsolatedMu = 0;
    unsigned int nMu=mucands->size();
    std::vector<bool> isos(nMu, false);
//     std::vector<double> MuonDeposits(nDep, 0);
    double MuonDeposits = 0.;
    unsigned int iMu=0;
    for (; iMu<nMu; iMu++) 
    {
      RecoChargedCandidateRef candref(mucands,iMu);
      
      LogDebug("HLTMuonPFIsoFilter") << "candref isNonnull " << candref.isNonnull(); 

     //did this candidate triggered at previous stage.
      if (!triggerdByPreviousLevel(candref,vcands)) continue;

     //reference to the track
      TrackRef tk = candref->get<TrackRef>();
      LogDebug("HLTMuonPFIsoFilter") << "tk isNonNull " << tk.isNonnull();

       //get the deposits and evaluate relIso if noDeltaBeta correction is applied
      if (!doDeltaB_)
      {
		for(unsigned int iDep=0;iDep!=nDep;++iDep)
		{
		  const edm::ValueMap<double> ::value_type & muonDeposit = (*(depMap[iDep]))[candref];
		  LogDebug("HLTMuonPFIsoFilter") << " Muon with q*pt= " << tk->charge()*tk->pt() << " (" << candref->charge()*candref->pt() << ") " << ", eta= " << tk->eta() << " (" << candref->eta() << ") " << "; has deposit["<<iDep<<"]: " << muonDeposit;

		  MuonDeposits += muonDeposit; 
		}
  	    MuonDeposits = MuonDeposits/tk->pt();
  	  }
  	  
  	  else if (doDeltaB_)
  	  {
  	    double neutralDeposits = 0.;
        for(unsigned int iDep=0;iDep!=nDep;++iDep)
		{
		  const edm::ValueMap<double> ::value_type & muonDeposit = (*(depMap[iDep]))[candref];
		  LogDebug("HLTMuonPFIsoFilter") << " Muon with q*pt= " << tk->charge()*tk->pt() << " (" << candref->charge()*candref->pt() << ") " << ", eta= " << tk->eta() << " (" << candref->eta() << ") " << "; has deposit["<<iDep<<"]: " << muonDeposit;

  	      std::size_t foundCharged = depTag_[iDep].label().find("Charged");
   	      if (foundCharged!=std::string::npos)  MuonDeposits += muonDeposit; 
  	      
  	      std::size_t foundGamma = depTag_[iDep].label().find("Gamma");
   	      if (foundGamma!=std::string::npos) neutralDeposits += muonDeposit;

  	      std::size_t foundNeutral = depTag_[iDep].label().find("Neutral");
   	      if (foundNeutral!=std::string::npos) neutralDeposits += muonDeposit;

  	      std::size_t foundPU = depTag_[iDep].label().find("PU");
   	      if (foundPU!=std::string::npos) neutralDeposits = neutralDeposits - 0.5*muonDeposit;
        }
        
        if (neutralDeposits > 0) MuonDeposits = (MuonDeposits + neutralDeposits)/tk->pt();
        else MuonDeposits = MuonDeposits/tk->pt();
//         MuonDeposits = MuonDeposits + max(neutralDeposits, 0);
        std::cout << "Muon deposits is " << MuonDeposits << std::endl;
  	  }
      
      
      
// RelIso = (pfIsolationR04().sumChargedHadronPt +
//              + max(pfIsolationR04().sumNeutralHadronPt +
//              + pfIsolationR04().sumPhotonEt
//              - pfIsolationR04().sumPUPt / 2.0, 0.0))/ Pt
       
      //get the selection
      if (MuonDeposits < maxIso_) isos[iMu] = true;

      LogDebug("HLTMuonPFIsoFilter") << " Muon with q*pt= " << tk->charge()*tk->pt() << ", eta= " << tk->eta() << "; "<<(isos[iMu]?"Is an isolated muon.":"Is NOT an isolated muon.");
       
      if (!isos[iMu]) continue;

      nIsolatedMu++;
      filterproduct.addObject(TriggerMuon,candref);
    }//for iMu

    // filter decision
    const bool accept (nIsolatedMu >= min_N_);

     //put the decision map
    if (nMu!=0)
    {
      edm::ValueMap<bool> ::Filler isoFiller(*PFisoMap);     

      // get a track ref
      TrackRef aRef = mucands->front().get<TrackRef>();

      // get the corresponding handle
      edm::Handle<reco::TrackCollection> HandleToTrackRef;
      iEvent.get(aRef.id(), HandleToTrackRef);
      isoFiller.insert(HandleToTrackRef, isos.begin(), isos.end());
      isoFiller.fill();
    }

    iEvent.put(PFisoMap);

    LogDebug("HLTMuonPFIsoFilter") << " >>>>> Result of HLTMuonIsoFilter is " << accept << ", number of muons passing isolation cuts= " << nIsolatedMu; 

    return accept;
 }
 

bool HLTMuonPFIsoFilter::triggerdByPreviousLevel(const reco::RecoChargedCandidateRef & candref, const std::vector<reco::RecoChargedCandidateRef>& vcands){
  bool ok=false;
  unsigned int i=0;
  unsigned int i_max=vcands.size();
  for (;i!=i_max;++i){
    if (candref == vcands[i]) { ok=true; break;}
  }

  return ok;
}
																						       
