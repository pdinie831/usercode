import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('Bc2Jpsi3Pi'
    ,tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
)
