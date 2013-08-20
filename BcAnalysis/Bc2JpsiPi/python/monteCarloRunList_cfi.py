import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles  = cms.untracked.vstring() 

readFiles.extend( 

 [
# 'file:/gwteray/users/fiorendi/dataReplica/MCBplus44HLTBPH/CAEA7A57-09D7-E111-ADA3-00215E222256.root'
 'file:/gwteray/users/fiorendi/dataReplica/MCjpsiPi44HLTBPH/FC18DF50-9BDC-E111-992A-002618943973.root'
# 'file:/gwteray/users/fiorendi/dataReplica/B0ToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgenRECO/F86E103F-72EE-E111-8009-E41F13181AB4.root'
# 'file:/gwteray/users/fiorendi/dataReplica/B0ToPsiMuMu_2MuPEtaFilter_Tight_7TeV-pythia6-evtgenRECO/F86E103F-72EE-E111-8009-E41F13181AB4.root'
 ] );

secFiles.extend( [  ] )
