# Load the input parameters parsing component (ARGV)
from HeavyFlavorAnalysis.Bc2JpsiPi.parseArgv_cfi import *

# Load the necessary ancillary configuration files
from HeavyFlavorAnalysis.Bc2JpsiPi.loadInit_cfi import *

# Instantiate a MessageLogger set of control parameters 
#service = MessageLogger { }

process.maxEvents = cms.untracked.PSet( input    	 = cms.untracked.int32  ( int(maxEvents))   )
process.options   = cms.untracked.PSet( 
                                        fileMode 	 = cms.untracked.string ( 'NOMERGE'         ),
					IgnoreCompletely = cms.untracked.vstring( 'myProducerLabel' ) 
				      )

# Limit the amount of output produced
if verbose == 0 :
	process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # print out 10 times only overall

import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt

from HLTrigger.special.hltPhysicsDeclared_cfi import *
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed

if dataType == "data":
	print lineno(), " "
	print lineno(), "Acquiring input DATA files list"
#For ReReco 
	process.GlobalTag.globaltag = 'FT_R_44_V11::All'

#For Prompt Reco	
#	process.GlobalTag.globaltag = 'GR_P_V22A::All'
# For reprocessing 4th November  
#	process.GlobalTag.globaltag = 'FT_R_39X_V4A::All'
	from HeavyFlavorAnalysis.Bc2JpsiPi.dataRunList_cfi       import *
else:
	print lineno(), " "
	print lineno(), "Acquiring input MonteCarlo files list"
        from Configuration.AlCa.autoCond import autoCond
        process.GlobalTag.globaltag = autoCond['mc']
	from HeavyFlavorAnalysis.Bc2JpsiPi.monteCarloRunList_cfi import *

# Specify list of input files
process.source = cms.Source ("PoolSource", 
                             fileNames          = readFiles, 
			     secondaryFileNames = secFiles)

process.hltPhysicsDeclared.L1GtReadoutRecordTag 	     	  = 'gtDigis'
seqPhysDeclBitSelection                         	     	  = cms.Sequence(hltPhysicsDeclared)

# accept if all path succeed (explicit AND)
#process.filter_all_explicit = hlt.triggerResultsFilter.clone(
#							     triggerConditions = cms.vstring( 'HLT_*', ),
#							     l1tResults = cms.InputTag(''),
#							     throw = cms.bool(True)
#							    )

# Add summary line at the end of analysis stream
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# bit 40+41 no BH  selection
process.bit4041 = hltLevel1GTSeed.clone(L1TechTriggerSeeding      = cms.bool(True), 
                                        L1SeedsLogicalExpression  = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39)'))

#bit 0 selection: BPTX_AND
process.bptxAnd = hltLevel1GTSeed.clone(L1TechTriggerSeeding      = cms.bool(True), 
                                        L1SeedsLogicalExpression  = cms.string('0'))

# one PV
process.oneGoodVertexFilter 	   = cms.EDFilter("VertexSelector",
				   		  src	          = cms.InputTag("offlinePrimaryVertices"),
				   		  cut	          = cms.string("!isFake && ndof >= 5 && abs(z) <= 24 && position.Rho <= 2"),
				   		  filter          = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
				   		 )

# against scraping events: high purity tracks
process.noScraping	    	   = cms.EDFilter("FilterOutScraping",
			    	  		  applyfilter     = cms.untracked.bool(True),
			    	  		  debugOn         = cms.untracked.bool(False), ## Or 'True' to get some per-event info
			    	  		  numtrack        = cms.untracked.uint32(10),
			    	  		  thresh          = cms.untracked.double(0.2)
			    	  		 )
					   
# one good displaced PV
process.goodDisplacedVertexFilter  = cms.EDFilter("VertexSelector",
				   		   src  	  = cms.InputTag("offlinePrimaryVertices"),
				   		   cut  	  = cms.string("!isFake && ndof >= 5" ),
				   		   filter	  = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
				   		 )


process.allTracksGenParticlesMatch = cms.EDFilter("MCTruthDeltaRMatcherNew",
						  src	      	  = cms.InputTag("generalTracks"),
						  distMin     	  = cms.double(0.15),
						  matchPDGId  	  = cms.vint32(),
						  matched     	  = cms.InputTag("genParticles")
						 )

process.offlinePrimaryVertices=cms.EDProducer("PrimaryVertexProducer",
					      PVSelParameters 	  = cms.PSet(
					        		  	     maxDistanceToBeam  = cms.double(0.05), ## 500 microns
					        		  	     minVertexFitProb   = cms.double(0.01)  ## 1% vertex fit probability
					        		  	    ),
					      verbose	      	  = cms.untracked.bool(False),
					      algorithm       	  = cms.string('AdaptiveVertexFitter'),
					      TkFilterParameters  = cms.PSet(
					        		  	      maxNormalizedChi2 = cms.double(5.0),
					        		  	      minSiliconHits	= cms.int32(7),     ## hits > 7
					        		  	      maxD0Significance = cms.double(5.0),  ## keep most primary tracks
					        		  	      minPt		= cms.double(0.0),  ## better for softish events
					        		  	      minPixelHits	= cms.int32(2)      ## hits > 2
					        		  	    ),
					      beamSpotLabel	  = cms.InputTag("offlineBeamSpot"),
					      TrackLabel	  = cms.InputTag("generalTracks"),                  # label of tracks to be used
					      useBeamConstraint   = cms.bool(False),
					      VtxFinderParameters = cms.PSet(
					        			     minTrackCompatibilityToOtherVertex = cms.double(0.01), ## 1%
					        			     minTrackCompatibilityToMainVertex  = cms.double(0.05), ## 5%
					        			     maxNbVertices			= cms.int32(0)      ## search all vertices in each cluster
					        			    ),
					      TkClusParameters = cms.PSet(zSeparation = cms.double(0.1))                            ## 1 mm max separation betw. clusters
					     )    
##filter on lumisections
if dataType == "data" and useJSON == 1:
	print lineno(), "Specifying lumi sections for DATA"
	from HeavyFlavorAnalysis.Bc2JpsiPi.JSON_cfi import *
	process.source.lumisToProcess = lumisToProcess

#process.filter = cms.Sequence(process.eIDRobustLoose+process.eIDRobustTight+process.eIDRobustHighEnergy+process.eIDLoose+process.eIDTight+process.allTracksGenParticlesMatch)
process.filter = cms.Sequence(process.allTracksGenParticlesMatch)
    
process.Bc2JpsiPi = cms.EDAnalyzer('Bc2JpsiPi',
                              inputDouble_     = cms.vstring ("cut_Pt_Mu=4.0",##tagli che posso cambiare
			      		                  	      "cut_Pt_Trk=0.9", 	    #1 not used 
            		      		        	      "cut_chi2n=3",
    	    		      		        	      "cut_cl_Jpsi=0.15",   	     #0.01 
    	    		      		        	      "cut_cl=0.001",   	     #0.01 
							                      "cut_cl_Bc=0.001",
    	    		      		        	      "cut_eta=2.2",
     	    		      		        	      "cut_eta_pi=2.4",
    	    		      		        	      "cut_nhits=6",
    	    		      		        	      "lSigmaCutmin=1", 	    #
    	    		      		        	      "lSigmaCutmax=12",	    #
							                      "JPsiMassPDG=3.096916",       #
							                      "JPsiMassWindow=.120",        # this is in GeV
							                      "pCandMomCut=0.",    #1
               					        	      "numberOfValidPixelHitsCut=2" #2
              					        	     ),
                              inputString_     = cms.vstring ("filename=" + output,
				                 			      "HLTname1=HLT_DoubleMu4_Jpsi_Displaced_v4",
							                      "HLTname2=HLT_DoubleMu4_Jpsi_Displaced_v1",
							                      "HLTname3=HLT_DoubleMu4_Jpsi_Displaced_v5",
                                                  "HLTnameReference1=HLT_Dimuon0_Jpsi_v6",
                                                  "HLTnameReference2=HLT_Dimuon0_Jpsi_v6",
                                                  "HLTMatchName=hltTriggerSummaryAOD",
 							                      "HLTMatchModule1=hltDisplacedmumuFilterDoubleMu4Jpsi",
							                      "HLTMatchModule2=hltDoubleMu4JpsiDisplacedL3Filtered",
			                        	          "trackCollection=generalTracks",
                                                  "MCID=Bc",
							                      "dump=false"
			                        	         ),
                              primaryVertexTag = cms.InputTag("offlinePrimaryVertices"),
                              beamSpotTag      = cms.InputTag("offlineBeamSpot"),    
    			              TkFilterParameters = cms.PSet(
			        			    maxNormalizedChi2	     = cms.double(5.0),
			        			    minSiliconLayersWithHits = cms.int32(5),	# >=5  (TDR > 7 hits)
			        			    minPixelLayersWithHits   = cms.int32(2),	# >=2  (TDR > 2 hits)
			        			    maxD0Significance	     = cms.double(5.0), # keep most primary tracks
			        			    minPt		     = cms.double(0.0), # better for softish events
			        			    trackQuality	     = cms.string("any")
    			        			   )
			  )
			     

if dataType == "data":
	print lineno(), "Specifying path for DATA"
	process.p = cms.Path(
##	                     seqPhysDeclBitSelection 	 *
##	                     process.bit4041         	 *
##			     process.bptxAnd         	 *
##			     process.filter_all_explicit *
##			     process.oneGoodVertexFilter *
##			     process.noScraping          *
			     process.Bc2JpsiPi       
			       			    )
else:
	print lineno(), "Specifying path for MonteCarlo"
	process.p = cms.Path(
##	                     seqPhysDeclBitSelection     *
##			     process.bit4041             *
#			     process.oneGoodVertexFilter *
#			     process.noScraping          *
			     process.Bc2JpsiPi
			    )


print lineno(), "\n\n"

