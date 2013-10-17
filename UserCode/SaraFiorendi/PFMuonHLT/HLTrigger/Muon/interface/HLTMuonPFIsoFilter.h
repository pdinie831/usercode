#ifndef HLTMuonPFIsoFilter_h
#define HLTMuonPFIsoFilter_h

/** \class HLTMuonPFIsoFilter
 *
 *  
 *  This class is an HLTFilter (-> EDFilter) implementing
 *  the PF isolation filtering for HLT muons
 *
 *  \author 
 *
 */

#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
//#include "RecoMuon/MuonIsolation/interface/MuIsoBaseIsolator.h"
#include "PhysicsTools/IsolationAlgos/plugins/CandIsolatorFromDeposits.h"

class HLTMuonPFIsoFilter : public HLTFilter {

   public:
      explicit HLTMuonPFIsoFilter(const edm::ParameterSet&);
      ~HLTMuonPFIsoFilter();
      virtual bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs & filterproduct);
      bool triggerdByPreviousLevel(const reco::RecoChargedCandidateRef &, const std::vector<reco::RecoChargedCandidateRef> &);
   private:
      edm::InputTag candTag_; // input tag identifying muon container
      edm::InputTag previousCandTag_;  // input tag identifying product contains muons passing the previous level
      std::vector<edm::InputTag> depTag_;  // input tag identifying deposit maps

//       const muonisolation::MuIsoBaseIsolator * theDepositIsolator;
      double maxIso_  ;         // max PF iso deposit allowed
      int    min_N_   ;         // minimum number of muons to fire the trigger
      bool   doDeltaB_; 		// if true, apply deltaBeta correction
};

#endif //HLTMuonPFIsoFilter_h
