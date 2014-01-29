# Load the input parameters parsing component (ARGV)
from HeavyFlavorAnalysis.Bc2Jpsi3Pi.parseArgv_cfi import *

# Load the necessary ancillary configuration files
from HeavyFlavorAnalysis.Bc2Jpsi3Pi.loadInit_cfi import *

# Instantiate a MessageLogger set of control parameters 
#service = MessageLogger { }

process.maxEvents = cms.untracked.PSet( input    = cms.untracked.int32(int(maxEvents)))
process.options   = cms.untracked.PSet( 
                                        fileMode 	 = cms.untracked.string ( 'NOMERGE'         ),
					                    IgnoreCompletely = cms.untracked.vstring( 'myProducerLabel' ),
					                    wantSummary = cms.untracked.bool(True) 
				      )

# Limit the amount of output produced
if verbose == 0 :
	process.MessageLogger.cerr.FwkReport.reportEvery =  200 

import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

from HLTrigger.special.hltPhysicsDeclared_cfi import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed

if dataType == "data":
	print lineno(), " "
	print lineno(), "Acquiring input DATA files list"
	process.GlobalTag.globaltag = 'FT_53_V21_AN4::All'
	from HeavyFlavorAnalysis.Bc2Jpsi3Pi.dataRunList_cfi       import *
else:
	print lineno(), " "
	print lineno(), "Acquiring input MonteCarlo files list"
        from Configuration.AlCa.autoCond import autoCond
        process.GlobalTag.globaltag = autoCond['mc']
	from HeavyFlavorAnalysis.Bc2Jpsi3Pi.monteCarloRunList_cfi import *

# Specify list of input files
process.source = cms.Source ("PoolSource", 
                             fileNames          = readFiles, 
			                 secondaryFileNames = secFiles
			                 )

process.hltPhysicsDeclared.L1GtReadoutRecordTag 	     	  = 'gtDigis'
seqPhysDeclBitSelection                         	     	  = cms.Sequence(hltPhysicsDeclared)


##filter on lumisections
if dataType == "data" and useJSON == 1:
	print lineno(), "Specifying lumi sections for DATA"
	from HeavyFlavorAnalysis.Bc2Jpsi3Pi.JSON_cfi import *
	process.source.lumisToProcess = lumisToProcess

process.Bc2Jpsi3Pi = cms.EDAnalyzer('Bc2Jpsi3Pi',
							  printTriggers    = cms.bool(False),
                              inputDouble_     = cms.vstring ("cut_Pt_Mu=0",
            		        	        	      "cut_Pt_Trk=0.5",
            		      		        	      "cut_chi2n=3",
												  "cut_cl_PV=0.001", 
												  "cut_cl_Jpsi=0.005", 
  						      					  "ptJpsi=7.9",
							    			      "cut_cl_Bc=0.001",
    	    		      		        	      "cut_eta=2.2",
    	    		      		        	      "cut_eta_pi=2.4",
    	    		      		        	      "cut_nhits=6",
							    				  "JPsiMassPDG=3.096916",       
							    				  "JPsiMassWindow=0.120",       
							    				  "pCandMomCut=0.",
					        		 			  "numberOfValidPixelHitsCut=2",
				        	                      "HLTMatch=0.5"
					        	     ),
                               inputString_ 	= cms.vstring ("filename=" + output,
												  "HLTname1=HLT_Dimuon8_Jpsi_v4",
												  "HLTname2=HLT_Dimuon8_Jpsi_v5",
												  "HLTname3=HLT_Dimuon8_Jpsi_v6",
												  "HLTname4=HLT_Dimuon8_Jpsi_v7",
							      				  "trackCollection=generalTracks",
												  "HLTnameReference1=HLT_Dimuon0_Jpsi_v1",
												  "HLTnameReference2=HLT_Dimuon0_Jpsi_v3",
												  "HLTMatchModule1=hltDimuon8JpsiL3Filtered",
												  "HLTMatchModule2=hltVertexmumuFilterDimuon8Jpsi",
							      				  "HLTMatchName=hltTriggerSummaryAOD",
							      				  "dump=false"
			                         ),
                              primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
                              beamSpotTag      = cms.InputTag("offlineBeamSpot"),    
			  )
			     
if dataType == "data":
	print lineno(), "Specifying path for DATA"
	process.p = cms.Path(
			     process.Bc2Jpsi3Pi
			    )
else:
	print lineno(), "Specifying path for MonteCarlo"
	process.p = cms.Path(
			     process.Bc2Jpsi3Pi
			    )

print lineno(), "\n\n"

