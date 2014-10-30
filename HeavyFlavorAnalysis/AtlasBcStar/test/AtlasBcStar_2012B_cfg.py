# Load the input parameters parsing component (ARGV)
from HeavyFlavorAnalysis.AtlasBcStar.parseArgv_cfi import *

# Load the necessary ancillary configuration files
from HeavyFlavorAnalysis.AtlasBcStar.loadInit_cfi import *

# Instantiate a MessageLogger set of control parameters 
#service = MessageLogger { }

process.maxEvents = cms.untracked.PSet( input    = cms.untracked.int32(-1))
# process.maxEvents = cms.untracked.PSet( input    = cms.untracked.int32(int(maxEvents)))
process.options   = cms.untracked.PSet( 
                                        fileMode 	 = cms.untracked.string ( 'NOMERGE'         ),
					                    IgnoreCompletely = cms.untracked.vstring( 'myProducerLabel' ),
					                    wantSummary = cms.untracked.bool(True) 
				      )

# Limit the amount of output produced
if verbose == 0 :
	process.MessageLogger.cerr.FwkReport.reportEvery =  200 # print out 10 times only overall

import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

from HLTrigger.special.hltPhysicsDeclared_cfi import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed

if dataType == "data":
	print lineno(), " "
	print lineno(), "Acquiring input DATA files list"
	process.GlobalTag.globaltag = 'FT_53_V21_AN4::All'
	from HeavyFlavorAnalysis.AtlasBcStar.dataRunList_cfi       import *
else:
	print lineno(), " "
	print lineno(), "Acquiring input MonteCarlo files list"
        from Configuration.AlCa.autoCond import autoCond
        process.GlobalTag.globaltag = autoCond['mc']
	from HeavyFlavorAnalysis.AtlasBcStar.monteCarloRunList_cfi import *

# Specify list of input files
process.source = cms.Source ("PoolSource", 
                             fileNames          = readFiles, 
			                 secondaryFileNames = secFiles,
# 			                 skipEvents=cms.untracked.uint32(2400)
			                 )

process.hltPhysicsDeclared.L1GtReadoutRecordTag 	     	  = 'gtDigis'
seqPhysDeclBitSelection                         	     	  = cms.Sequence(hltPhysicsDeclared)

##filter on lumisections
if dataType == "data" and useJSON == 1:
	print lineno(), "Specifying lumi sections for DATA"
	from HeavyFlavorAnalysis.AtlasBcStar.JSON_cfi import *
	process.source.lumisToProcess = lumisToProcess   

    
process.AtlasBcStar = cms.EDAnalyzer('AtlasBcStar',
							  printTriggers    = cms.bool(False),
							  cut_pt_Mu        = cms.double(  0),
							  cut_pt_Trk       = cms.double(2.0),
							  cut_pt_Pi        = cms.double(0.4),
							  cut_pt_Jpsi      = cms.double(7.9),
							  cut_cl_PV        = cms.double(0.001),
							  cut_cl_Jpsi      = cms.double(0.005),
							  cut_cl_Bc        = cms.double(0.02),
							  cut_eta_mu       = cms.double(2.1),
							  cut_eta_pi       = cms.double(2.4),
							  cut_jpsi_Mass    = cms.double(0.120),
							  cut_HLT_match    = cms.double(0.5),
							  cut_chi2_trk     = cms.double(3.),
							  cut_nPix_hits    = cms.uint32(1),
							  cut_nTrk_hits    = cms.uint32(6),
							  MCID             = cms.string("Bc"),
							  HLTname1         = cms.string("HLT_Dimuon8_Jpsi_v4"),
							  HLTname2         = cms.string("HLT_Dimuon8_Jpsi_v5"),
							  HLTname3         = cms.string("HLT_Dimuon8_Jpsi_v6"),
							  HLTname4         = cms.string("HLT_Dimuon8_Jpsi_v7"),
							  HLTnameReference1= cms.string("HLT_Dimuon0_Jpsi_v1"),
							  HLTnameReference2= cms.string("HLT_Dimuon0_Jpsi_v3"),
							  HLTMatchModule1  = cms.string("hltDimuon8JpsiL3Filtered"),
							  HLTMatchModule2  = cms.string("hltVertexmumuFilterDimuon8Jpsi"),
							  HLTMatchName     = cms.string("hltTriggerSummaryAOD"),
							  doGenMC          = cms.bool(False),
							  filename         = cms.string(output),
							  trackCollection  = cms.InputTag("generalTracks"),
                              primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
                              beamSpotTag      = cms.InputTag("offlineBeamSpot")    
			  )
			     
if dataType == "data":
	print lineno(), "Specifying path for DATA"
	process.p = cms.Path(
			     process.AtlasBcStar
			    )
else:
	print lineno(), "Specifying path for MonteCarlo"
	process.p = cms.Path(
			     process.AtlasBcStar
			    )

print lineno(), "\n\n"

