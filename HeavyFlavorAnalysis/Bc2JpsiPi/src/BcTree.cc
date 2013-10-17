#include "TLorentzVector.h"
#include "TProcessID.h"
#include "TDirectory.h"
#include "../interface/BcTree.h"
//
ClassImp(BcTreeHeader)
ClassImp(BcTree)
ClassImp(BcTreeCand)
ClassImp(BcTreeGENCand)
ClassImp(JpsiCand)

TClonesArray *BcTree::fgBcTreeArrayCand     = 0;
TClonesArray *BcTree::fgBcTreeJpsiArrayCand = 0;

//---------------------------------------------------------------------------------

BcTree::BcTree() : jIsValid(kFALSE)
{
    if (!fgBcTreeArrayCand) 
	{
      std::cout<<"Create fBcTreeArrayCand"<<std::endl;
      fgBcTreeArrayCand = new TClonesArray("BcTreeCand",1000);
    }  

    if (!fgBcTreeJpsiArrayCand) 
	{
      std::cout<<"Create fBcTreeJpsiArrayCand"<<std::endl;
      fgBcTreeJpsiArrayCand = new TClonesArray("JpsiCand",100);
    }  

    fBcTreeJpsiArrayCand = fgBcTreeJpsiArrayCand;
    NJpsiCand	= 0 ;

    fBcTreeArrayCand 	 = fgBcTreeArrayCand;
    fMaxPtCand	= new TRefArray;
    NBcTreeCand	= 0 ;
    NJpsiCand	= 0 ;
    MaxPt		= 0.;
    JPsiPxOld	= 0.;
    JPsiPyOld	= 0.;
    JPsiPzOld	= 0.;
    NumJPsi  	= 0 ;
 
//   TProcessID::SetObjectCount(ObjectNumber);
}

BcTree::~BcTree()
{
   BcTreeClear();
}
//
//======================================================================
//
void BcTree::BcTreeClear(Option_t * /*option*/)
{
       fBcTreeArrayCand->Clear("C");
       fBcTreeJpsiArrayCand->Clear("C");
	   ClearBcTreeHeader();
	   ClearBcTreeGENCand();
       fMaxPtCand->Delete();
       NBcTreeCand  = 0;
       NJpsiCand  = 0;
       MaxPt	    = 0.;
       JPsiPxOld	= 0.;
       JPsiPyOld	= 0.;
       JPsiPzOld	= 0.;
       NumJPsi  	= 0 ;
 }
//
//======================================================================
//

void BcTree::SetBcTreeHeader( UInt_t NEve  ,
 			      UInt_t  Run   ,
 			      UInt_t  LBlk  ,
			      UInt_t  Event ,
			      UInt_t  NumBc ,
			      Int_t   Ntrk  ,
			      Int_t   Nprim ,
				  Float_t TrueNI,
				  Int_t TriggerBit,
				  Int_t GoodTrkSize
				  )
	
{
       NBcTreeCand = 0;
       NJpsiCand   = 0;
       MaxPt	   = 0.;
       JPsiPxOld   = 0.;
       JPsiPyOld   = 0.;
       JPsiPzOld   = 0.;
       fBcTreeHeader.SetNEve(  NEve  );
       fBcTreeHeader.SetRun(   Run   );
       fBcTreeHeader.SetLBlk(  LBlk  );
       fBcTreeHeader.SetEvent( Event );
       fBcTreeHeader.SetNumBc( NumBc );
       fBcTreeHeader.SetNtrk ( Ntrk  );
       fBcTreeHeader.SetNprim( Nprim );
       fBcTreeHeader.SetTrueNI(TrueNI);
       fBcTreeHeader.SetTriggerBit(TriggerBit);
       fBcTreeHeader.SetGoodTrkSize(GoodTrkSize);

	}
//
//======================================================================
//
void BcTree::ClearBcTreeHeader()
{
       fBcTreeHeader.SetNEve(  0);
       fBcTreeHeader.SetRun(   0);
       fBcTreeHeader.SetLBlk(  0);
       fBcTreeHeader.SetEvent( 0);
       fBcTreeHeader.SetNumBc( 0);
       fBcTreeHeader.SetNtrk ( 0);
       fBcTreeHeader.SetNprim( 0);
       fBcTreeHeader.SetTrueNI(0);
       fBcTreeHeader.SetTriggerBit(0);
       fBcTreeHeader.SetGoodTrkSize(0);

}
//
//======================================================================
//
void BcTree::ClearBcTreeGENCand()
{
    TLorentzVector *empty = new TLorentzVector (0,0,0,0);
    fBcTreeGENCand.SetBcPi    (*empty)  ;
    fBcTreeGENCand.SetMu1	  (*empty)  ;
	fBcTreeGENCand.SetMu2	  (*empty)  ;
	fBcTreeGENCand.SetJPsi    (*empty)  ;
	fBcTreeGENCand.SetPi	  (*empty)  ;
 	fBcTreeGENCand.SetEl	  (0 )  ;
 	fBcTreeGENCand.SetCos  	  (0 )  ;
 	fBcTreeGENCand.SetVtxSx	  (0 )  ;
 	fBcTreeGENCand.SetVtxSy	  (0 )  ;
 	fBcTreeGENCand.SetVtxSz	  (0 )  ;
 	fBcTreeGENCand.SetBSx	  (0 )  ;
 	fBcTreeGENCand.SetBSy	  (0 )  ;
 	fBcTreeGENCand.SetBSz	  (0 )  ;
	fBcTreeGENCand.SetPiCh	  (0 )  ;
	fBcTreeGENCand.SetMuCh1	  (0 )  ;
	fBcTreeGENCand.SetMuCh2	  (0 )  ;
}
//======================================================================
//
BcTreeCand::BcTreeCand(const BcTreeCand &orig) : TObject(orig)
{
  		BcPi_      		=  orig.BcPi_          ;
  		BcK_       		=  orig.BcK_           ;
 		Mu1_       		=  orig.Mu1_           ;
 		Mu2_       		=  orig.Mu2_           ;
 		Jpsi_      		=  orig.Jpsi_          ;
 		Pi_        		=  orig.Pi_            ;
 		JpsiV_      	=  orig.JpsiV_         ;
 		Els_	   		=  orig.Els_           ;
 		Cos_	   		=  orig.Cos_           ;
 		ClS_	   		=  orig.ClS_           ;
 		ClJpsi_	   		=  orig.ClJpsi_        ;
		IP3D_      		=  orig.IP3D_          ;
		IPbs_      		=  orig.IPbs_          ;
		sIP3D_     		=  orig.sIP3D_	       ;
  		sIPbs_     		=  orig.sIPbs_	       ;
  		MatchMu_[0]		=  orig.MatchMu_[0]    ;
 		MatchMu_[1]		=  orig.MatchMu_[1]    ;
 		MatchPi_   		=  orig.MatchPi_       ;
 		PiCh_      		=  orig.PiCh_          ;
 		MuCh_[0]   		=  orig.MuCh_[0]       ;
 		MuCh_[1]   		=  orig.MuCh_[1]       ;
        ElsigJpsi_ 		=  orig.ElsigJpsi_     ;
        CosJpsi_   		=  orig.CosJpsi_       ;
 		TrkNHits_    	=  orig.TrkNHits_      ; 
 		TrkPixelHits_   =  orig.TrkPixelHits_  ; 
 		TrkChi2_    	=  orig.TrkChi2_       ; 
 		
 		for (int i=0; i < 3; i++)
 		{
 		  BS_[i]		= 	orig.BS_[i]			;
 		  VtxS_[i]		= 	orig.VtxS_[i]		;
 		}  
 		 
 		for(int j=0; j<9; j++)
 		{
 		  BSCovariance_[j] 	= orig.BSCovariance_[j]   ;
 		  SVCovariance_[j] 	= orig.SVCovariance_[j]   ;
 		}

 		for (int k=0; k < 4; k++)
 		{
 		  HltMatch_[k]		= 	orig.HltMatch_[k]		;
 		}  
 }
//======================================================================
//
BcTreeCand::BcTreeCand( 
 	        const TLorentzVector&     BcPi  	   ,
 	        const TLorentzVector&     BcK   	   ,
            const TLorentzVector&     Mu1   	   ,
            const TLorentzVector&     Mu2   	   ,
            const TLorentzVector&     Jpsi  	   ,
            const TLorentzVector&     Pi    	   ,
            const TLorentzVector&     JpsiV  	   ,
 	        Double_t     			  Els		   ,
 	        Double_t     			  Cos		   ,
 	        Double_t     			  ClS		   ,
 	        Double_t     			  ClJpsi       ,
 	        Double_t     			  BS[3]		   ,
 	        Double_t     			  VtxS[3]	   ,
            Double_t     			  BSCovariance[9],
 	        Double_t     			  SVCovariance[9],
			Double_t     			  IP3D  	   ,
			Double_t     			  IPbs  	   ,
			Double_t     			  sIP3D 	   ,
			Double_t     			  sIPbs 	   ,	 
			Double_t     			  MatchMu[2]   ,
 			Double_t     			  MatchPi	   ,
 	        Int_t	     			  PiCh         ,
	        Int_t	     			  MuCh[2]      ,
 	        Double_t                  ElsigJpsi    ,  
 	        Double_t                  CosJpsi      ,  
  	        Int_t                     TrkNHits     ,  
  	        Int_t                     TrkPixelHits ,  
  	        Float_t                    TrkChi2      ,              
			Double_t     			  HltMatch[4]
                      )
 {
	        BcPi_       =  BcPi 	  		 ;
	        BcK_        =  BcK  	  		 ;
 			Mu1_        =  Mu1  	  		 ;
 			Mu2_        =  Mu2  	  		 ;
 			Jpsi_       =  Jpsi 	  		 ;
 			Pi_         =  Pi		  		 ;
 			JpsiV_      =  JpsiV 	  		 ;
	        Els_	    =  Els  	  		 ;
	        Cos_	    =  Cos  	  		 ;
	        ClS_	    =  ClS  	  		 ;
	        ClJpsi_	    =  ClJpsi  	  		 ;
			IP3D_       =  IP3D 	  		 ;
			IPbs_       =  IPbs 	  		 ;
			sIP3D_      =  sIP3D	  		 ;
			sIPbs_      =  sIPbs	  		 ;
 		 	MatchMu_[0] =  MatchMu[0] 		 ;
 		 	MatchMu_[1] =  MatchMu[1] 		 ;
 		 	MatchPi_    =  MatchPi    		 ;
			PiCh_       =  PiCh	      		 ;
			MuCh_[0]    =  MuCh[0]    		 ;
			MuCh_[1]    =  MuCh[1] 	  		 ;
			ElsigJpsi_  =  ElsigJpsi  		 ; 
			CosJpsi_    =  CosJpsi    		 ; 
			TrkNHits_   =  TrkNHits     	 ;      
			TrkPixelHits_ =  TrkPixelHits    ;      
			TrkChi2_    =  TrkChi2     		 ;     
			
			for (int i=0; i < 3; i++)
			{
			  BS_[i]		= 	BS[i]		;
			  VtxS_[i]		= 	VtxS[i]		;
			} 
			
			for(int j=0; j<9; j++)
			{
			  BSCovariance_[j] 		= BSCovariance[j]  ;
			  SVCovariance_[j] 		= SVCovariance[j]  ;
			}

 			for (int k=0; k < 4; k++)
 			{
 		  	  HltMatch_[k]		= 	HltMatch[k]		;
 			}  
}		        

//======================================================================
//
BcTreeCand::~BcTreeCand() {
}

//======================================================================
//
BcTreeCand *BcTree::AddBcTreeCand(    
   		    const TLorentzVector&	BcPi	       ,
   		    const TLorentzVector&	BcK 	       ,
            const TLorentzVector&	Mu1 	       ,
            const TLorentzVector&	Mu2 	       ,
            const TLorentzVector&	Jpsi	       ,
            const TLorentzVector&	Pi		       ,
            const TLorentzVector&	JpsiV	       ,
	        Double_t     			Els 	       ,
	        Double_t     			Cos 	       ,
	        Double_t     			ClS 	       ,
	        Double_t     			ClJpsi 	       ,
	        Double_t     			BS[3]          ,
	        Double_t     			VtxS[3]        ,
 	        Double_t     			BSCovariance[9],
 	        Double_t     			SVCovariance[9],
		    Double_t	 			IP3D	       ,
		    Double_t	 			IPbs	       ,
		    Double_t	 			sIP3D	       ,
		    Double_t	 			sIPbs	       ,
		    Double_t	 			MatchMu[2]     ,
		    Double_t	 			MatchPi        ,
	        Int_t	     			PiCh           ,
	        Int_t	     			MuCh[2]        ,
	        Double_t 				ElsigJpsi      ,    
	        Double_t 				CosJpsi        ,   
     	    Int_t	     			TrkNHits       ,   
     	    Int_t	     			TrkPixelHits   , 
     	    Float_t	     			TrkChi2		   ,
     	    Double_t				HltMatch[4]      
				   ){
   // Add a new track to the list of tracks for this event.
   // To avoid calling the very time consuming operator new for each track,
   // the standard but not well know C++ operator "new with placement"
   // is called. If tracks[i] is 0, a new Track object will be created
   // otherwise the previous Track[i] will be overwritten.

   TClonesArray &BcTreeArrayCand = *fBcTreeArrayCand ;
   BcTreeCand *_BcTreeCand = new(BcTreeArrayCand[BcTree::NBcTreeCand++])  
   BcTreeCand(  BcPi		,
   				BcK     	,
   				Mu1     	,
				Mu2     	,
				Jpsi    	,
			    Pi      	,
				JpsiV   	,
 	          	Els			,    	   
 				Cos			,	     
 				ClS			,	     
 				ClJpsi  	,	     
 		    	BS 			,	  
 		    	VtxS		,	  
 		    	BSCovariance,	  
 		    	SVCovariance,	  
				IP3D    	,
				IPbs    	,
				sIP3D   	,
				sIPbs   	,
 				MatchMu 	,
 				MatchPi 	,
	        	PiCh		,
	        	MuCh		,
	        	ElsigJpsi   ,      
				CosJpsi     ,
	        	TrkNHits    ,
	        	TrkPixelHits,
				TrkChi2	    ,
				HltMatch
	        	);  
 
    fIsValid = kTRUE;

//
//Save reference to last BcTreeCand in the collection of BcTreeCand
//
    Double_t mxBcPt = BcPi.Pt();
    if( BcPi.Pt() > MaxPt )
	{
      MaxPt      = mxBcPt ;
      fMaxPtCand -> Add( _BcTreeCand );
    }
    
	Double_t JPsiPx = JpsiV.Px();
    Double_t JPsiPy = JpsiV.Py();
    Double_t JPsiPz = JpsiV.Pz();
    if(JPsiPx!=JPsiPxOld&&JPsiPy!=JPsiPyOld&&JPsiPz!=JPsiPzOld)
	{
	  JPsiPxOld = JPsiPx;
      JPsiPyOld = JPsiPy;
      JPsiPzOld = JPsiPz;
      NumJPsi++;
    } 

   return _BcTreeCand;
}

//======================================================================
//GENPARTICLES
//======================================================================

void BcTree::SetBcTreeGENCand(    
	        TLorentzVector  BcPi    ,
	        TLorentzVector  Mu1     ,
	        TLorentzVector  Mu2     ,
	        TLorentzVector  JPsi    ,
	        TLorentzVector  Pi      ,
	        Double_t     	El      ,
	        Double_t     	Cos     ,
	        Double_t     	VtxS[3] ,
	        Double_t     	BS[3]   ,
	        Int_t	     	PiCh    ,
	        Int_t	     	MuCh[2] 
				   ){

    fBcTreeGENCand.SetBcPi    (BcPi    )  ;
    fBcTreeGENCand.SetMu1	  (Mu1     )  ;
	fBcTreeGENCand.SetMu2	  (Mu2     )  ;
	fBcTreeGENCand.SetJPsi    (JPsi    )  ;
	fBcTreeGENCand.SetPi	  (Pi      )  ;
 	fBcTreeGENCand.SetEl	  (El	   )  ;
 	fBcTreeGENCand.SetCos  	  (Cos     )  ;
 	fBcTreeGENCand.SetVtxSx	  (VtxS[0] )  ;
 	fBcTreeGENCand.SetVtxSy	  (VtxS[1] )  ;
 	fBcTreeGENCand.SetVtxSz	  (VtxS[2] )  ;
 	fBcTreeGENCand.SetBSx	  (BS[0]   )  ;
 	fBcTreeGENCand.SetBSy	  (BS[1]   )  ;
 	fBcTreeGENCand.SetBSz	  (BS[2]   )  ;
	fBcTreeGENCand.SetPiCh	  (PiCh    )  ;
	fBcTreeGENCand.SetMuCh1	  (MuCh[0] )  ;
	fBcTreeGENCand.SetMuCh2	  (MuCh[1] )  ;
	 
}

//======================================================================
//JpsiCand 
//======================================================================

JpsiCand::JpsiCand(const JpsiCand &orig) : TObject(orig)
{
 		Mu1_       		=  orig.Mu1_           ;
 		Mu2_       		=  orig.Mu2_           ;
 		Jpsi_      		=  orig.Jpsi_          ;
 		Els_	   		=  orig.Els_           ;
 		Cos_	   		=  orig.Cos_           ;
 		ClS_	   		=  orig.ClS_           ;
  		MatchMu_[0]		=  orig.MatchMu_[0]    ;
 		MatchMu_[1]		=  orig.MatchMu_[1]    ;
 		MuCh_[0]   		=  orig.MuCh_[0]       ;
 		MuCh_[1]   		=  orig.MuCh_[1]       ;
 		
 		for (int i=0; i < 3; i++)
 		{
 		  VtxJ_[i]		= 	orig.VtxJ_[i]		;
 		}  
 		 
 }
//
//======================================================================

JpsiCand::JpsiCand( 
            const TLorentzVector&     Mu1   	   ,
            const TLorentzVector&     Mu2   	   ,
            const TLorentzVector&     Jpsi  	   ,
 	        Double_t     			  Els		   ,
 	        Double_t     			  Cos		   ,
 	        Double_t     			  ClS		   ,
 	        Double_t     			  VtxJ[3]	   ,
			Double_t     			  MatchMu[2]   ,
	        Int_t	     			  MuCh[2]      
                      )
 {
 			Mu1_        =  Mu1  	  		 ;
 			Mu2_        =  Mu2  	  		 ;
 			Jpsi_       =  Jpsi 	  		 ;
	        Els_	    =  Els  	  		 ;
	        Cos_	    =  Cos  	  		 ;
	        ClS_	    =  ClS  	  		 ;
 		 	MatchMu_[0] =  MatchMu[0] 		 ;
 		 	MatchMu_[1] =  MatchMu[1] 		 ;
			MuCh_[0]    =  MuCh[0]    		 ;
			MuCh_[1]    =  MuCh[1] 	  		 ;

			for (int i=0; i < 3; i++)
			{
			  VtxJ_[i]		= 	VtxJ[i]		;
			} 
}		        
//
//======================================================================

JpsiCand::~JpsiCand() {
}
//======================================================================
//
JpsiCand *BcTree::AddJpsiCand(    
            const TLorentzVector&	Mu1 	       ,
            const TLorentzVector&	Mu2 	       ,
            const TLorentzVector&	Jpsi	       ,
	        Double_t     			Els 	       ,
	        Double_t     			Cos 	       ,
	        Double_t     			ClS 	       ,
	        Double_t     			VtxJ[3]        ,
		    Double_t	 			MatchMu[2]     ,
	        Int_t	     			MuCh[2]        
				   ){

   TClonesArray &BcTreeJpsiArrayCand = *fBcTreeJpsiArrayCand ;

   JpsiCand *_JpsiCand = new(BcTreeJpsiArrayCand[BcTree::NJpsiCand++])  
   JpsiCand(    Mu1     	,
				Mu2     	,
				Jpsi    	,
 	          	Els			,    	   
 				Cos			,	     
 				ClS			,	     
 		    	VtxJ		,	  
 				MatchMu 	,
	        	MuCh		
	        	);  
 
    jIsValid = kTRUE;

   return _JpsiCand;
}
