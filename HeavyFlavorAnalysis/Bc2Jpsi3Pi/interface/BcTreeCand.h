#ifndef BcTreeCand_h
#define BcTreeCand_h
#include "TObject.h"
#include "TLorentzVector.h"

//################################################
//  BcTreeCand Class
//################################################

class BcTreeCand : public TObject {

private:

    // Bc cand
 	TLorentzVector 	BcCand_ 	   			;

	// Jpsi cand
 	TLorentzVector  Jpsi_  	     			;
 	TLorentzVector  JpsiV_  	    		;
 	Float_t         DCA_					;

	// Muons
 	TLorentzVector 	MuP_   	     			;
 	TLorentzVector 	MuM_   	     			;

 	bool 			MuPGlobal_				;
 	bool 			MuPTracker_				;
 	bool 			MuPPFlow_				;
 	bool 			MuPTMOST_				;
	Int_t 			MuPTrkLayerMeas_		;
	Int_t 			MuPPixLayerMeas_		;
	Int_t 			MuPPixHits_				;
	Int_t 			MuPTrkHits_				;
	Int_t 			MuPMatchedStations_		;
	Float_t 		MuPNormChi2_			;
    Float_t 		MuPDxy_					;
    Float_t 		MuPDz_					;

 	bool 			MuMGlobal_				;
 	bool 			MuMTracker_				;
 	bool 			MuMPFlow_				;
 	bool 			MuMTMOST_				;
	Int_t 			MuMTrkLayerMeas_		;
	Int_t 			MuMPixLayerMeas_		;
	Int_t 			MuMPixHits_				;
	Int_t 			MuMTrkHits_				;
	Int_t 			MuMMatchedStations_		;
	Float_t 		MuMNormChi2_			;
    Float_t 		MuMDxy_					;
    Float_t 		MuMDz_					;
 	
	// Tracks 
 	TLorentzVector 	Pi1_   	     			;
 	TLorentzVector 	Pi2_   	     			;
 	TLorentzVector 	Pi3_   	     			;

 	Int_t	     	Pi1Ch_           		;  	 
	Int_t 			Trk1PixLayerMeas_		;
	Int_t 			Trk1TrkLayerMeas_		;
	Int_t 			Trk1PixHits_			;
	Int_t 			Trk1TrkHits_			;
 	Float_t		 	Trk1NormChi2_       	;
	Double_t     	Trk1IP3DWrtJpsi_ 	    ;
	Double_t     	Trk1IP3DWrtJpsiSign_ 	;
	Double_t     	Trk1IP2DWrtBS_ 	        ;
	Double_t     	Trk1IP2DWrtBSSign_      ;
	Double_t     	Trk1IP3DWrtPV_ 	    	;
	Double_t     	Trk1IP3DWrtPVSign_	 	;
	Double_t 		Trk1DeltaR_				;

 	Int_t	     	Pi2Ch_           		;  	 
	Int_t 			Trk2PixLayerMeas_		;
	Int_t 			Trk2TrkLayerMeas_		;
	Int_t 			Trk2PixHits_			;
	Int_t 			Trk2TrkHits_			;
 	Float_t		 	Trk2NormChi2_       	;
	Double_t     	Trk2IP3DWrtJpsi_ 	    ;
	Double_t     	Trk2IP3DWrtJpsiSign_ 	;
	Double_t     	Trk2IP2DWrtBS_ 	        ;
	Double_t     	Trk2IP2DWrtBSSign_      ;
	Double_t     	Trk2IP3DWrtPV_ 	    	;
	Double_t     	Trk2IP3DWrtPVSign_	 	;
	Double_t 		Trk2DeltaR_				;

 	Int_t	     	Pi3Ch_           		;  	 
	Int_t 			Trk3PixLayerMeas_		;
	Int_t 			Trk3TrkLayerMeas_		;
	Int_t 			Trk3PixHits_			;
	Int_t 			Trk3TrkHits_			;
 	Float_t		 	Trk3NormChi2_       	;
	Double_t     	Trk3IP3DWrtJpsi_ 	    ;
	Double_t     	Trk3IP3DWrtJpsiSign_ 	;
	Double_t     	Trk3IP2DWrtBS_ 	        ;
	Double_t     	Trk3IP2DWrtBSSign_      ;
	Double_t     	Trk3IP3DWrtPV_ 	    	;
	Double_t     	Trk3IP3DWrtPVSign_	 	;
	Double_t 		Trk3DeltaR_				;

	Float_t 		TrkMuPDR_[3]			;
	Float_t 		TrkMuMDR_[3]			;

	// Bc vertex 
 	Double_t     	ClS_     	     		;
 	Double_t     	BcVtxPosition_[3] 	    ;
  	Double_t     	BcVtxCovariance_[9]		;
 	Double_t     	El2DWrtBS_	     		;
 	Double_t     	Els2DWrtBS_	     		;
 	Double_t     	Sigma2DWrtBS_     		;
 	Double_t     	Cos2DWrtBS_	     		;
 	Double_t     	El3DWrtPV_	     		;
 	Double_t     	Els3DWrtPV_	     		;
 	Double_t     	Sigma3DWrtPV_	   		;
 	Double_t     	Cos3DWrtPV_	     		;

	// Jpsi vertex 
 	Double_t     	ClJpsi_  	     		;
    Double_t     	ElsigJpsi_      		;
	Double_t     	CosJpsi_        		;
 	Double_t     	JpsiVtxPosition_[3] 	;

	// Pointing primary vertex 
 	Double_t     	PointPVPosition_[3] 	;
  	Double_t     	PointPVCovariance_[9]	;
  	Double_t     	PointPVCl_				;

	// MC matching 
 	Double_t	 MatchMuP_                  ;   
 	Double_t	 MatchMuM_       			;   
 	Double_t     MatchPi1_[2]      			;   
 	Double_t     MatchPi2_[2]      			;   
 	Double_t     MatchPi3_[2]      			;   

	// HLT matching 
  	Double_t     HltMatch_[4]	 			;



public:

    BcTreeCand(){;}
    ~BcTreeCand() {}


// ------- Set information into ntuple ----------------------

    // Bc cand
  	void   SetBcCand        (TLorentzVector n) 			{   BcCand_   				= n ;}

	// Jpsi cand
  	void   SetJpsi          (TLorentzVector n) 			{   Jpsi_    				= n ;}
  	void   SetJpsiV         (TLorentzVector n) 			{   JpsiV_    				= n ;}
  	void   SetDCA           (Float_t	    n) 			{   DCA_    				= n ;}

	// Muons	
  	void   SetMuP           (TLorentzVector n) 			{   MuP_  					= n ;}
  	void   SetMuM           (TLorentzVector n) 			{   MuM_    				= n ;}
	
  	void   SetMuPisGlobal   (bool			n) 			{   MuPGlobal_ 				= n ;}
  	void   SetMuPisTracker  (bool			n) 			{   MuPTracker_ 			= n ;}
  	void   SetMuPisPFlow    (bool			n) 			{   MuPPFlow_ 				= n ;}
  	void   SetMuPTMOST      (bool			n) 			{   MuPTMOST_ 				= n ;}
  	void   SetMuPTrkLMeas   (Int_t			n) 			{   MuPTrkLayerMeas_		= n ;}
  	void   SetMuPPixLMeas   (Int_t			n) 			{   MuPPixLayerMeas_		= n ;}
  	void   SetMuPPixHits    (Int_t			n) 			{   MuPPixHits_ 			= n ;}
  	void   SetMuPTrkHits    (Int_t			n) 			{   MuPTrkHits_ 			= n ;}
  	void   SetMuPMatchedStat(Int_t			n) 			{   MuPMatchedStations_ 	= n ;}
  	void   SetMuPNormChi2   (Float_t		n) 			{   MuPNormChi2_ 			= n ;}
  	void   SetMuPDxy   		(Float_t		n) 			{   MuPDxy_ 				= n ;}
  	void   SetMuPDz   		(Float_t		n) 			{   MuPDz_ 					= n ;}

  	void   SetMuMisGlobal   (bool			n) 			{   MuMGlobal_ 				= n ;}
  	void   SetMuMisTracker  (bool			n) 			{   MuMTracker_ 			= n ;}
  	void   SetMuMisPFlow    (bool			n) 			{   MuMPFlow_ 				= n ;}
  	void   SetMuMTMOST      (bool			n) 			{   MuMTMOST_ 				= n ;}
  	void   SetMuMTrkLMeas   (Int_t			n) 			{   MuMTrkLayerMeas_		= n ;}
  	void   SetMuMPixLMeas   (Int_t			n) 			{   MuMPixLayerMeas_		= n ;}
  	void   SetMuMPixHits    (Int_t			n) 			{   MuMPixHits_ 			= n ;}
  	void   SetMuMTrkHits    (Int_t			n) 			{   MuMTrkHits_ 			= n ;}
  	void   SetMuMMatchedStat(Int_t			n) 			{   MuMMatchedStations_ 	= n ;}
  	void   SetMuMNormChi2   (Float_t		n) 			{   MuMNormChi2_ 			= n ;}
  	void   SetMuMDxy   		(Float_t		n) 			{   MuMDxy_ 				= n ;}
  	void   SetMuMDz   		(Float_t		n) 			{   MuMDz_ 					= n ;}

	// Track 
  	void   SetPi1           	(TLorentzVector n) 		{   Pi1_  					= n ;}
  	void   SetPi1Ch   			(Int_t			n) 		{   Pi1Ch_					= n ;}
  	void   SetTrk1PixLMeas  	(Int_t			n) 		{   Trk1PixLayerMeas_		= n ;}
  	void   SetTrk1TrkLMeas  	(Int_t			n) 		{   Trk1TrkLayerMeas_		= n ;}
  	void   SetTrk1PixHits   	(Int_t			n) 		{   Trk1PixHits_ 			= n ;}
  	void   SetTrk1TrkHits   	(Int_t			n) 		{   Trk1TrkHits_ 			= n ;}
  	void   SetTrk1NormChi2  	(Float_t		n) 		{   Trk1NormChi2_ 			= n ;}
  	void   SetTrk1IP3DJpsi	    (Double_t		n) 		{   Trk1IP3DWrtJpsi_ 		= n ;}
  	void   SetTrk1IP3DJpsiSign  (Double_t		n) 		{   Trk1IP3DWrtJpsiSign_	= n ;}
  	void   SetTrk1IP2DBS	    (Double_t		n) 		{   Trk1IP2DWrtBS_ 			= n ;}
  	void   SetTrk1IP2DBSSign    (Double_t		n) 		{   Trk1IP2DWrtBSSign_		= n ;}
  	void   SetTrk1IP3DPV	    (Double_t		n) 		{   Trk1IP3DWrtPV_ 			= n ;}
  	void   SetTrk1IP3DPVSign  	(Double_t		n) 		{   Trk1IP3DWrtPVSign_		= n ;}
  	void   SetTrk1DeltaR	  	(Double_t		n) 		{   Trk1DeltaR_				= n ;}

  	void   SetPi2           	(TLorentzVector n) 		{   Pi2_  					= n ;}
  	void   SetPi2Ch   			(Int_t			n) 		{   Pi2Ch_					= n ;}
  	void   SetTrk2PixLMeas  	(Int_t			n) 		{   Trk2PixLayerMeas_		= n ;}
  	void   SetTrk2TrkLMeas  	(Int_t			n) 		{   Trk2TrkLayerMeas_		= n ;}
  	void   SetTrk2PixHits   	(Int_t			n) 		{   Trk2PixHits_ 			= n ;}
  	void   SetTrk2TrkHits   	(Int_t			n) 		{   Trk2TrkHits_ 			= n ;}
  	void   SetTrk2NormChi2  	(Float_t		n) 		{   Trk2NormChi2_ 			= n ;}
  	void   SetTrk2IP3DJpsi	    (Double_t		n) 		{   Trk2IP3DWrtJpsi_ 		= n ;}
  	void   SetTrk2IP3DJpsiSign  (Double_t		n) 		{   Trk2IP3DWrtJpsiSign_	= n ;}
  	void   SetTrk2IP2DBS	    (Double_t		n) 		{   Trk2IP2DWrtBS_ 			= n ;}
  	void   SetTrk2IP2DBSSign    (Double_t		n) 		{   Trk2IP2DWrtBSSign_		= n ;}
  	void   SetTrk2IP3DPV	    (Double_t		n) 		{   Trk2IP3DWrtPV_ 			= n ;}
  	void   SetTrk2IP3DPVSign  	(Double_t		n) 		{   Trk2IP3DWrtPVSign_		= n ;}
  	void   SetTrk2DeltaR	  	(Double_t		n) 		{   Trk2DeltaR_				= n ;}

  	void   SetPi3           	(TLorentzVector n) 		{   Pi3_  					= n ;}
  	void   SetPi3Ch   			(Int_t			n) 		{   Pi3Ch_					= n ;}
  	void   SetTrk3PixLMeas  	(Int_t			n) 		{   Trk3PixLayerMeas_		= n ;}
  	void   SetTrk3TrkLMeas  	(Int_t			n) 		{   Trk3TrkLayerMeas_		= n ;}
  	void   SetTrk3PixHits   	(Int_t			n) 		{   Trk3PixHits_ 			= n ;}
  	void   SetTrk3TrkHits   	(Int_t			n) 		{   Trk3TrkHits_ 			= n ;}
  	void   SetTrk3NormChi2  	(Float_t		n) 		{   Trk3NormChi2_ 			= n ;}
  	void   SetTrk3IP3DJpsi	    (Double_t		n) 		{   Trk3IP3DWrtJpsi_ 		= n ;}
  	void   SetTrk3IP3DJpsiSign  (Double_t		n) 		{   Trk3IP3DWrtJpsiSign_	= n ;}
  	void   SetTrk3IP2DBS	    (Double_t		n) 		{   Trk3IP2DWrtBS_ 			= n ;}
  	void   SetTrk3IP2DBSSign    (Double_t		n) 		{   Trk3IP2DWrtBSSign_		= n ;}
  	void   SetTrk3IP3DPV	    (Double_t		n) 		{   Trk3IP3DWrtPV_ 			= n ;}
  	void   SetTrk3IP3DPVSign  	(Double_t		n) 		{   Trk3IP3DWrtPVSign_		= n ;}
  	void   SetTrk3DeltaR	  	(Double_t		n) 		{   Trk3DeltaR_				= n ;}

  	void   SetTrackMuPDR	  	(int id, Double_t n)	{   TrkMuPDR_[id]			= n ;}
  	void   SetTrackMuMDR	  	(int id, Double_t n)	{   TrkMuMDR_[id]			= n ;}

	// Bc vertex 
  	void   SetClS		    (Double_t		n) 			{   ClS_	 				= n ;}
  	void   SetEl2DWrtBS     (Double_t		n) 			{   El2DWrtBS_	 			= n ;}
  	void   SetEls2DWrtBS    (Double_t		n) 			{   Els2DWrtBS_	 			= n ;}
  	void   SetSigma2DWrtBS  (Double_t		n) 			{   Sigma2DWrtBS_			= n ;}
  	void   SetCos2DWrtBS    (Double_t		n) 			{   Cos2DWrtBS_	 			= n ;}
  	void   SetEl3DWrtPV     (Double_t		n) 			{   El3DWrtPV_	 			= n ;}
  	void   SetEls3DWrtPV    (Double_t		n) 			{   Els3DWrtPV_	 			= n ;}
  	void   SetSigma3DWrtPV  (Double_t		n) 			{   Sigma3DWrtPV_	 		= n ;}
  	void   SetCos3DWrtPV    (Double_t		n) 			{   Cos3DWrtPV_	 			= n ;}

  	void   SetBcVtxPosition   (int id, Double_t	n) 	    {   BcVtxPosition_[id]		= n ;}
  	void   SetBcVtxCovariance (int id, Double_t	n)		{   BcVtxCovariance_[id]	= n ;}

	// Jpsi vertex 
  	void   SetClJpsi	    (Double_t		n) 			{   ClJpsi_	 				= n ;}
  	void   SetElsigJpsi	    (Double_t		n) 			{   ElsigJpsi_	 			= n ;}
  	void   SetCosJpsi	    (Double_t		n) 			{   CosJpsi_	 			= n ;}
  	
  	void   SetJpsiVtxPosition   (int id, Double_t	n)	{   JpsiVtxPosition_[id]	= n ;}

	// Pointing primary vertex 
  	void   SetPointPVPosition   (int id, Double_t	n)	{   PointPVPosition_[id]	= n ;}
  	void   SetPointPVCovariance (int id, Double_t	n)  {   PointPVCovariance_[id]	= n ;}
  	void   SetPointPVCl		    (Double_t			n) 	{   PointPVCl_	 			= n ;}

	// MC matching 
  	void   SetMatchMuP	    (Double_t		n) 			{   MatchMuP_	 			= n ;}
  	void   SetMatchMuM	    (Double_t		n) 			{   MatchMuM_	 			= n ;}
  	void   SetMatchPi1	    (int id, Double_t n) 		{   MatchPi1_[id] 			= n ;}
  	void   SetMatchPi2	    (int id, Double_t n) 		{   MatchPi2_[id]			= n ;}
  	void   SetMatchPi3	    (int id, Double_t n) 		{   MatchPi3_[id]			= n ;}

	// HLT matching 
  	void   SetHltMatch 		(int id, Double_t	n)		{   HltMatch_[id]			= n ;}



// ------- Retrieve information from ntuple ----------------------

    // Bc cand
  	const TLorentzVector   *GetBcCand()          		{   return  &BcCand_   			;}

	// Jpsi cand
  	const TLorentzVector   *GetJpsi()      				{   return  &Jpsi_   			;}
  	const TLorentzVector   *GetJpsiV()    				{   return  &JpsiV_  			;}
  	Float_t   				GetDCA() 					{   return  DCA_ 				;}

	// Muons
  	const TLorentzVector   *GetMuP()           			{   return  &MuP_  				;}
  	const TLorentzVector   *GetMuM()           			{   return  &MuM_   			;}

  	bool   					GetMuPisGlobal() 			{   return  MuPGlobal_ 			;}
  	bool   					GetMuPisTracker() 			{   return  MuPTracker_ 		;}
  	bool   					GetMuPisPFlow() 			{   return  MuPPFlow_ 			;}
  	bool   					GetMuPTMOST() 				{   return  MuPTMOST_ 			;}
  	Int_t   				GetMuPTrkLMeas() 			{   return  MuPTrkLayerMeas_	;}
  	Int_t   				GetMuPPixLMeas() 			{   return  MuPPixLayerMeas_	;}
  	Int_t   				GetMuPPixHits() 			{   return  MuPPixHits_ 		;}
  	Int_t   				GetMuPTrkHits() 			{   return  MuPTrkHits_ 		;}
  	Int_t   				GetMuPMatchedStat() 		{   return  MuPMatchedStations_ ;}
  	Float_t   				GetMuPNormChi2() 			{   return  MuPNormChi2_ 		;}
  	Float_t   				GetMuPDxy() 				{   return  MuPDxy_ 			;}
  	Float_t   				GetMuPDz() 					{   return  MuPDz_ 				;}
	
  	bool   					GetMuMisGlobal() 			{   return  MuMGlobal_ 		 	;}
  	bool   					GetMuMisTracker() 			{   return  MuMTracker_ 		;}
  	bool   					GetMuMisPFlow() 			{   return  MuMPFlow_ 			;}
  	bool   					GetMuMTMOST() 				{   return  MuMTMOST_ 			;}
  	Int_t  					GetMuMTrkLMeas() 			{   return  MuMTrkLayerMeas_	;}
  	Int_t  					GetMuMPixLMeas() 			{   return  MuMPixLayerMeas_	;}
  	Int_t  					GetMuMPixHits() 			{   return  MuMPixHits_ 		;}
  	Int_t  					GetMuMTrkHits() 			{   return  MuMTrkHits_ 		;}
  	Int_t  					GetMuMMatchedStat() 		{   return  MuMMatchedStations_ ;}
  	Float_t 				GetMuMNormChi2() 			{   return  MuMNormChi2_ 	 	;}
  	Float_t 				GetMuMDxy() 				{   return  MuMDxy_ 		 	;}
  	Float_t 				GetMuMDz() 					{   return  MuMDz_ 			 	;}
	
	// Track 
  	const TLorentzVector   *GetPi1() 					{   return  &Pi1_  				;}
  	Int_t   				GetPi1Ch() 					{   return  Pi1Ch_				;}
  	Int_t   				GetTrk1TrkPixLMeas() 		{   return  Trk1PixLayerMeas_	;}
  	Int_t   				GetTrk1TrkTrkLMeas() 		{   return  Trk1TrkLayerMeas_	;}
  	Int_t   				GetTrk1TrkPixHits() 		{   return  Trk1PixHits_ 		;}
  	Int_t   				GetTrk1TrkTrkHits() 		{   return  Trk1TrkHits_ 		;}
  	Float_t   				GetTrk1TrkNormChi2() 		{   return  Trk1NormChi2_ 		;}
  	Double_t   				GetTrk1IP3DJpsi() 			{   return  Trk1IP3DWrtJpsi_ 	;}
  	Double_t   				GetTrk1IP3DJpsiSign() 		{   return  Trk1IP3DWrtJpsiSign_;}
  	Double_t   				GetTrk1IP2DBS() 			{   return  Trk1IP2DWrtBS_ 		;}
  	Double_t   				GetTrk1IP2DBSSign() 		{   return  Trk1IP2DWrtBSSign_	;}
  	Double_t   				GetTrk1IP3DPV() 			{   return  Trk1IP3DWrtPV_ 		;}
  	Double_t   				GetTrk1IP3DPVSign() 		{   return  Trk1IP3DWrtPVSign_	;}
  	Double_t   				GetTrk1DeltaR() 			{   return  Trk1DeltaR_			;}

  	const TLorentzVector   *GetPi2() 					{   return  &Pi2_  				;}
  	Int_t   				GetPi2Ch() 					{   return  Pi2Ch_				;}
  	Int_t   				GetTrk2TrkPixLMeas() 		{   return  Trk2PixLayerMeas_	;}
  	Int_t   				GetTrk2TrkTrkLMeas() 		{   return  Trk2TrkLayerMeas_	;}
  	Int_t   				GetTrk2TrkPixHits() 		{   return  Trk2PixHits_ 		;}
  	Int_t   				GetTrk2TrkTrkHits() 		{   return  Trk2TrkHits_ 		;}
  	Float_t   				GetTrk2TrkNormChi2() 		{   return  Trk2NormChi2_ 		;}
  	Double_t   				GetTrk2IP3DJpsi() 			{   return  Trk2IP3DWrtJpsi_ 	;}
  	Double_t   				GetTrk2IP3DJpsiSign() 		{   return  Trk2IP3DWrtJpsiSign_;}
  	Double_t   				GetTrk2IP2DBS() 			{   return  Trk2IP2DWrtBS_ 		;}
  	Double_t   				GetTrk2IP2DBSSign() 		{   return  Trk2IP2DWrtBSSign_	;}
  	Double_t   				GetTrk2IP3DPV() 			{   return  Trk2IP3DWrtPV_ 		;}
  	Double_t   				GetTrk2IP3DPVSign() 		{   return  Trk2IP3DWrtPVSign_	;}
  	Double_t   				GetTrk2DeltaR() 			{   return  Trk2DeltaR_			;}

  	const TLorentzVector   *GetPi3() 					{   return  &Pi3_  				;}
  	Int_t   				GetPi3Ch() 					{   return  Pi3Ch_				;}
  	Int_t   				GetTrk3TrkPixLMeas() 		{   return  Trk3PixLayerMeas_	;}
  	Int_t   				GetTrk3TrkTrkLMeas() 		{   return  Trk3TrkLayerMeas_	;}
  	Int_t   				GetTrk3TrkPixHits() 		{   return  Trk3PixHits_ 		;}
  	Int_t   				GetTrk3TrkTrkHits() 		{   return  Trk3TrkHits_ 		;}
  	Float_t   				GetTrk3TrkNormChi2() 		{   return  Trk3NormChi2_ 		;}
  	Double_t   				GetTrk3IP3DJpsi() 			{   return  Trk3IP3DWrtJpsi_ 	;}
  	Double_t   				GetTrk3IP3DJpsiSign() 		{   return  Trk3IP3DWrtJpsiSign_;}
  	Double_t   				GetTrk3IP2DBS() 			{   return  Trk3IP2DWrtBS_ 		;}
  	Double_t   				GetTrk3IP2DBSSign() 		{   return  Trk3IP2DWrtBSSign_	;}
  	Double_t   				GetTrk3IP3DPV() 			{   return  Trk3IP3DWrtPV_ 		;}
  	Double_t   				GetTrk3IP3DPVSign() 		{   return  Trk3IP3DWrtPVSign_	;}
  	Double_t   				GetTrk3DeltaR() 			{   return  Trk3DeltaR_			;}

  	Float_t   				GetTrackMuPDR(int i)		{  return   TrkMuPDR_[i]		;}
  	Float_t   				GetTrackMuMDR(int i)		{  return   TrkMuMDR_[i]		;}

	// Bc vertex 
  	Double_t   				GetClS() 					{   return  ClS_	 		    ;}
  	Double_t   				GetEl2DWrtBS() 				{   return  El2DWrtBS_	 		;}
  	Double_t   				GetEls2DWrtBS() 			{   return  Els2DWrtBS_	 		;}
  	Double_t   				GetSigma2DWrtBS() 			{   return  Sigma2DWrtBS_		;}
  	Double_t   				GetCos2DWrtBS() 			{   return  Cos2DWrtBS_	 		;}
  	Double_t   				GetEl3DWrtPV() 				{   return  El3DWrtPV_	 		;}
  	Double_t   				GetEls3DWrtPV() 			{   return  Els3DWrtPV_	 		;}
  	Double_t   				GetSigma3DWrtPV() 			{   return  Sigma3DWrtPV_	 	;}
  	Double_t   				GetCos3DWrtPV() 			{   return  Cos3DWrtPV_	 		;}
  	Double_t   				GetBcVtxPosition(int i) 	{   return  BcVtxPosition_[i]   ;}
  	Double_t   				GetBcVtxCovariance(int i)	{   return  BcVtxCovariance_[i] ;}

	// Jpsi vertex 
  	Double_t   				GetClJpsi() 				{   return  ClJpsi_	 			 ;}
  	Double_t   				GetElsigJpsi() 				{   return  ElsigJpsi_	 		 ;}
  	Double_t   				GetCosJpsi() 				{   return  CosJpsi_	 	 	 ;}
  	Double_t   				GetJpsiVtxPosition(int i)	{   return  JpsiVtxPosition_[i]  ;}

	// Pointing primary vertex 
  	Double_t   				GetPointPVPosition(int i) 	{ 	return  PointPVPosition_[i]  ;}
  	Double_t   				GetPointPVCovariance(int i)	{	return  PointPVCovariance_[i];}
  	Double_t   				GetPointPVCl() 				{   return  PointPVCl_	 		;}

	// MC matching 
  	Double_t   				GetMatchMuP() 				{   return  MatchMuP_	 	    ;}
  	Double_t   				GetMatchMuM() 				{   return  MatchMuM_	 		;}
  	Double_t   				GetMatchPi1(int i) 			{   return  MatchPi1_[i] 		;}
  	Double_t   				GetMatchPi2(int i) 			{   return  MatchPi2_[i] 		;}
  	Double_t   				GetMatchPi3(int i) 			{   return  MatchPi3_[i] 		;}

	// HLT matching 
  	Double_t   				GetHltMatch(int i)	 		{	return  HltMatch_[i]		;}




// ------- Copy a BcCand into another ----------------------

BcTreeCand(const BcTreeCand &orig) : TObject(orig)
{
    // Bc cand
  	BcCand_      			=  orig.BcCand_          		;

	// Jpsi cand
 	Jpsi_      				=  orig.Jpsi_          			;
 	JpsiV_      			=  orig.JpsiV_         			;
 	DCA_      				=  orig.DCA_         			;

	// Muons
 	MuP_       				=  orig.MuP_           			;
 	MuM_       				=  orig.MuM_           			;

 	MuPGlobal_        		=  orig.MuPGlobal_            	;
 	MuPTracker_        		=  orig.MuPTracker_            	;
 	MuPPFlow_        		=  orig.MuPPFlow_            	;
 	MuPTMOST_        		=  orig.MuPTMOST_            	;
 	MuPTrkLayerMeas_        =  orig.MuPTrkLayerMeas_        ;
 	MuPPixLayerMeas_        =  orig.MuPPixLayerMeas_        ;
 	MuPPixHits_        		=  orig.MuPPixHits_            	;
 	MuPTrkHits_        		=  orig.MuPTrkHits_            	;
 	MuPMatchedStations_     =  orig.MuPMatchedStations_     ;
 	MuPNormChi2_        	=  orig.MuPNormChi2_            ;
 	MuPDxy_        			=  orig.MuPDxy_            		;
 	MuPDz_        			=  orig.MuPDz_            		;

 	MuMGlobal_        		=  orig.MuMGlobal_            	;
 	MuMTracker_        		=  orig.MuMTracker_            	;
 	MuMPFlow_        		=  orig.MuMPFlow_            	;
 	MuMTMOST_        		=  orig.MuMTMOST_            	;
 	MuMTrkLayerMeas_        =  orig.MuMTrkLayerMeas_        ;
 	MuMPixLayerMeas_        =  orig.MuMPixLayerMeas_        ;
 	MuMPixHits_        		=  orig.MuMPixHits_            	;
 	MuMTrkHits_        		=  orig.MuMTrkHits_            	;
 	MuMMatchedStations_     =  orig.MuMMatchedStations_     ;
 	MuMNormChi2_        	=  orig.MuMNormChi2_            ;
 	MuMDxy_        			=  orig.MuMDxy_            		;
 	MuMDz_        			=  orig.MuMDz_            		;

	// Track 
 	Pi1_        			=  orig.Pi1_            		;
 	Pi1Ch_        		 	=  orig.Pi1Ch_        		  	;
 	Trk1PixLayerMeas_	    =  orig.Trk1PixLayerMeas_	  	;
 	Trk1TrkLayerMeas_	    =  orig.Trk1TrkLayerMeas_	  	;
 	Trk1PixHits_			=  orig.Trk1PixHits_			;
 	Trk1TrkHits_			=  orig.Trk1TrkHits_			;
 	Trk1NormChi2_         	=  orig.Trk1NormChi2_          	;
 	Trk1IP3DWrtJpsi_ 	    =  orig.Trk1IP3DWrtJpsi_ 	    ;
 	Trk1IP3DWrtJpsiSign_    =  orig.Trk1IP3DWrtJpsiSign_    ;
 	Trk1IP2DWrtBS_ 	        =  orig.Trk1IP2DWrtBS_ 	        ;
 	Trk1IP2DWrtBSSign_      =  orig.Trk1IP2DWrtBSSign_      ;
 	Trk1IP3DWrtPV_ 	        =  orig.Trk1IP3DWrtPV_ 	        ;
 	Trk1IP3DWrtPVSign_	    =  orig.Trk1IP3DWrtPVSign_	    ;
 	Trk1DeltaR_			    =  orig.Trk1DeltaR_			    ;

 	Pi2_        			=  orig.Pi2_            		;
 	Pi2Ch_        		 	=  orig.Pi2Ch_        		  	;
 	Trk2PixLayerMeas_	    =  orig.Trk2PixLayerMeas_	  	;
 	Trk2TrkLayerMeas_	    =  orig.Trk2TrkLayerMeas_	  	;
 	Trk2PixHits_			=  orig.Trk2PixHits_			;
 	Trk2TrkHits_			=  orig.Trk2TrkHits_			;
 	Trk2NormChi2_         	=  orig.Trk2NormChi2_          	;
 	Trk2IP3DWrtJpsi_ 	    =  orig.Trk2IP3DWrtJpsi_ 	    ;
 	Trk2IP3DWrtJpsiSign_    =  orig.Trk2IP3DWrtJpsiSign_    ;
 	Trk2IP2DWrtBS_ 	        =  orig.Trk2IP2DWrtBS_ 	        ;
 	Trk2IP2DWrtBSSign_      =  orig.Trk2IP2DWrtBSSign_      ;
 	Trk2IP3DWrtPV_ 	        =  orig.Trk2IP3DWrtPV_ 	        ;
 	Trk2IP3DWrtPVSign_	    =  orig.Trk2IP3DWrtPVSign_	    ;
 	Trk2DeltaR_			    =  orig.Trk2DeltaR_			    ;

 	Pi3_        			=  orig.Pi3_            		;
 	Pi3Ch_        		 	=  orig.Pi3Ch_        		  	;
 	Trk3PixLayerMeas_	    =  orig.Trk3PixLayerMeas_	  	;
 	Trk3TrkLayerMeas_	    =  orig.Trk3TrkLayerMeas_	  	;
 	Trk3PixHits_			=  orig.Trk3PixHits_			;
 	Trk3TrkHits_			=  orig.Trk3TrkHits_			;
 	Trk3NormChi2_         	=  orig.Trk3NormChi2_          	;
 	Trk3IP3DWrtJpsi_ 	    =  orig.Trk3IP3DWrtJpsi_ 	    ;
 	Trk3IP3DWrtJpsiSign_    =  orig.Trk3IP3DWrtJpsiSign_    ;
 	Trk3IP2DWrtBS_ 	        =  orig.Trk3IP2DWrtBS_ 	        ;
 	Trk3IP2DWrtBSSign_      =  orig.Trk3IP2DWrtBSSign_      ;
 	Trk3IP3DWrtPV_ 	        =  orig.Trk3IP3DWrtPV_ 	        ;
 	Trk3IP3DWrtPVSign_	    =  orig.Trk3IP3DWrtPVSign_	    ;
 	Trk3DeltaR_			    =  orig.Trk3DeltaR_			    ;
 	
 	for(int j=0; j<3; j++){
 	  TrkMuPDR_[j] 			= orig.TrkMuPDR_[j]				;
 	  TrkMuMDR_[j] 			= orig.TrkMuMDR_[j]				;
	}
	// Bc vertex 
 	ClS_     	     	 	=  orig.ClS_			      	;
 	El2DWrtBS_	     	 	=  orig.El2DWrtBS_	    	   	;
 	Els2DWrtBS_	     	 	=  orig.Els2DWrtBS_	    	    ;
 	Sigma2DWrtBS_     	 	=  orig.Sigma2DWrtBS_       	;
 	Cos2DWrtBS_	     	 	=  orig.Cos2DWrtBS_	          	;
 	El3DWrtPV_	     	 	=  orig.El3DWrtPV_	          	;
 	Els3DWrtPV_	     	 	=  orig.Els3DWrtPV_	          	;
 	Sigma3DWrtPV_	   	 	=  orig.Sigma3DWrtPV_	      	;
 	Cos3DWrtPV_	     	 	=  orig.Cos3DWrtPV_	          	;
 	for(int j=0; j<3; j++){
 	  BcVtxPosition_[j] 	=  orig.BcVtxPosition_[j]		;
 	}
 	for(int j=0; j<9; j++){
 	  BcVtxCovariance_[j]	=  orig.BcVtxCovariance_[j]		;
 	}

	// Jpsi vertex 
 	ClJpsi_     	     	=  orig.ClS_			      	;
 	ElsigJpsi_	     	 	=  orig.ElsigJpsi_	    	   	;
 	CosJpsi_	     	 	=  orig.CosJpsi_	    	    ;
 	for(int j=0; j<3; j++){
 	  JpsiVtxPosition_[j] 	=  orig.JpsiVtxPosition_[j]		;
 	}

	// Pointing primary vertex 
 	for(int j=0; j<3; j++){
 	  PointPVPosition_[j] 	=  orig.PointPVPosition_[j]		;
 	}
 	for(int j=0; j<9; j++){
 	  PointPVCovariance_[j]	=  orig.PointPVCovariance_[j]	;
 	}
 	PointPVCl_			 	=  orig.PointPVCl_				;

	// MC matching 
 	MatchMuP_     	     	=  orig.MatchMuP_			   	;
 	MatchMuM_	     	 	=  orig.MatchMuM_	    	   	;
 	for(int i=0; i<9; i++){
	  MatchPi1_[i]	     	=  orig.MatchPi1_[i]    	    ;
	  MatchPi2_[i]	     	=  orig.MatchPi2_[i]    	    ;
	  MatchPi3_[i]	     	=  orig.MatchPi3_[i]    	    ;
    }
	// HLT matching 
 	for(int j=0; j<4; j++){
 	  HltMatch_[j]	 		=  orig.HltMatch_[j]			;
 	}

}

    ClassDef(BcTreeCand,1) 
};
#endif