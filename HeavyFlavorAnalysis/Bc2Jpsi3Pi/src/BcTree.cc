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
      fgBcTreeJpsiArrayCand = new TClonesArray("JpsiCand",1000);
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
       NBcTreeCand  = 0 ;
       NJpsiCand    = 0;
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
//       BcTreeHeader  fBcTreeHeader;
       NBcTreeCand = 0;
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
    fBcTreeGENCand.SetBc3Pi   (*empty)  ;
    fBcTreeGENCand.SetMu1	  (*empty)  ;
	fBcTreeGENCand.SetMu2	  (*empty)  ;
	fBcTreeGENCand.SetJPsi    (*empty)  ;
	fBcTreeGENCand.SetPi1	  (*empty)  ;
	fBcTreeGENCand.SetPi2	  (*empty)  ;
	fBcTreeGENCand.SetPi3	  (*empty)  ;
 	fBcTreeGENCand.SetEl	  (0 )  ;
 	fBcTreeGENCand.SetCos  	  (0 )  ;
	fBcTreeGENCand.SetBSx	  (0 )  ;  
 	fBcTreeGENCand.SetBSy	  (0 )  ; 
 	fBcTreeGENCand.SetBSz	  (0 )  ;
 	fBcTreeGENCand.SetVtxSx	  (0 )  ;
 	fBcTreeGENCand.SetVtxSy	  (0 )  ;
 	fBcTreeGENCand.SetVtxSz	  (0 )  ;
	fBcTreeGENCand.SetRis0 	  (0 )  ;
	fBcTreeGENCand.SetRis1 	  (0 )  ;
	fBcTreeGENCand.SetRis2 	  (0 )  ;
	fBcTreeGENCand.SetRis3 	  (0 )  ;
	fBcTreeGENCand.SetRis4 	  (0 )  ;
	fBcTreeGENCand.SetPiCh1	  (0 )  ;
	fBcTreeGENCand.SetPiCh2	  (0 )  ;
	fBcTreeGENCand.SetPiCh3	  (0 )  ;
	fBcTreeGENCand.SetMuCh1	  (0 )  ;
	fBcTreeGENCand.SetMuCh2	  (0 )  ;
	 
}
//
//======================================================================
//
BcTreeCand::BcTreeCand(const BcTreeCand &orig) : TObject(orig)
{
  		Bc3Pi_    		   =  orig.Bc3Pi_      		   ;
 		Mu1_      		   =  orig.Mu1_        		   ;
 		Mu2_      		   =  orig.Mu2_        		   ;
 		Jpsi_     		   =  orig.Jpsi_       		   ;
 		Pi1_      		   =  orig.Pi1_        		   ;
 		Pi2_      		   =  orig.Pi2_        		   ;
 		Pi3_      		   =  orig.Pi3_        		   ;

 		tPi1_      		   =  orig.tPi1_       		   ;
 		tPi2_      		   =  orig.tPi2_       		   ;
 		tPi3_      		   =  orig.tPi3_       		   ;

 		JpsiV_    		   =  orig.JpsiV_      		   ;
		Ris_[0]   		   =  orig.Ris_[0]     		   ;
		Ris_[1]   		   =  orig.Ris_[1]     		   ;
		Ris_[2]   		   =  orig.Ris_[2]     		   ;
		Ris_[3]   		   =  orig.Ris_[3]     		   ;
		Ris_[4]   		   =  orig.Ris_[4]     		   ;
 		Els_	  		   =  orig.Els_        		   ;
 		El_		  		   =  orig.El_        		   ;
 		Cos_	  		   =  orig.Cos_        		   ;
 		ClS_	  		   =  orig.ClS_        		   ;
 		ClJpsi_	  		   =  orig.ClJpsi_     		   ;
		IP3D_[0]  		   =  orig.IP3D_[0]    		   ;
		IP3D_[1]  		   =  orig.IP3D_[1]    		   ;
		IP3D_[2]  		   =  orig.IP3D_[2]    		   ;
		IPbs_[0]  		   =  orig.IPbs_[0]    		   ;
		IPbs_[1]  		   =  orig.IPbs_[1]    		   ;
		IPbs_[2]  		   =  orig.IPbs_[2]    		   ;
		sIP3D_[0] 		   =  orig.sIP3D_[0]   		   ;
		sIP3D_[1] 		   =  orig.sIP3D_[1]   		   ;
		sIP3D_[2] 		   =  orig.sIP3D_[2]   		   ;
		sIPbs_[0] 		   =  orig.sIPbs_[0]   		   ;
		sIPbs_[1] 		   =  orig.sIPbs_[1]   		   ;
		sIPbs_[2] 		   =  orig.sIPbs_[2]   		   ;
 		MatchMu_[0]		   =  orig.MatchMu_[0] 		   ;
 		MatchMu_[1]		   =  orig.MatchMu_[1] 		   ;
 		MatchPi1_[0]	   =  orig.MatchPi1_[0]		   ;
 		MatchPi1_[1]	   =  orig.MatchPi1_[1]		   ;
 		MatchPi2_[0]	   =  orig.MatchPi2_[0]		   ;
 		MatchPi2_[1]	   =  orig.MatchPi2_[1]		   ;
 		MatchPi3_[0]	   =  orig.MatchPi3_[0]		   ;
 		MatchPi3_[1]	   =  orig.MatchPi3_[1]		   ;
 		MuCh_[0]    	   =  orig.MuCh_[0]    		   ;
 		MuCh_[1]   		   =  orig.MuCh_[1]    		   ;
        ElsigJpsi_ 		   =  orig.ElsigJpsi_     	   ;
        CosJpsi_   		   =  orig.CosJpsi_            ;
        
        
   		for (int i=0; i < 3; i++)
 		{
 		  BS_[i]		    =  orig.BS_[i]		      ;
 		  VtxS_[i]		    =  orig.VtxS_[i]		  ;
 		  PiCh_[i]          =  orig.PiCh_[i]          ;
 		  TrkNHits_[i] 	    =  orig.TrkNHits_[i]      ; 
   		  TrkPixelHits_[i]  =  orig.TrkPixelHits_[i]  ; 
 		  TrkChi2_[i]    	=  orig.TrkChi2_[i]       ; 
 		  
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
//
//======================================================================
//
BcTreeCand::BcTreeCand( 
 	        const TLorentzVector&     Bc3Pi ,
            const TLorentzVector&     Mu1   ,
            const TLorentzVector&     Mu2   ,
            const TLorentzVector&     Jpsi  ,
            const TLorentzVector&     Pi1   ,
            const TLorentzVector&     Pi2   ,
            const TLorentzVector&     Pi3   ,
            const TLorentzVector&     tPi1  ,
            const TLorentzVector&     tPi2  ,
            const TLorentzVector&     tPi3  ,
            const TLorentzVector&     JpsiV ,
		    Double_t     Ris[5]  	 		,
 	        Double_t     El     	 		,
 	        Double_t     Els     	 		,
 	        Double_t     Cos     	 		,
 	        Double_t     ClS     	 		,
 	        Double_t     ClJpsi    	 		,
 	        Double_t     BS[3] 	 		    ,
 	        Double_t     VtxS[3] 	 		,
            Double_t     BSCovariance[9]    ,
 	        Double_t     SVCovariance[9]    ,
			Double_t     IP3D[3] 	 		,
			Double_t     IPbs[3] 	 		,
			Double_t     sIP3D[3]	 		,
			Double_t     sIPbs[3]	 		,
			Double_t     MatchMu[2]  		,
 			Double_t     MatchPi1[2] 		,
 			Double_t     MatchPi2[2] 		,
 			Double_t     MatchPi3[2] 		,
	        Int_t	     PiCh[3]     		,
	        Int_t	     MuCh[2]     		,
 	        Double_t     ElsigJpsi   		, 
 	        Double_t     CosJpsi     		, 
  	        Int_t        TrkNHits[3]        ,  
  	        Int_t        TrkPixelHits[3]    ,  
  	        Float_t      TrkChi2[3]         ,              
			Double_t     HltMatch[4]
                       )
 {
	        Bc3Pi_     			 =  Bc3Pi 	           ;
 			Mu1_       			 =  Mu1         	   ;
 			Mu2_       			 =  Mu2         	   ;
 			Jpsi_      			 =  Jpsi        	   ;
 			Pi1_       			 =  Pi1         	   ;
 			Pi2_       			 =  Pi2         	   ;
 			Pi3_       			 =  Pi3         	   ;
 			tPi1_       		 =  tPi1         	   ;
 			tPi2_       		 =  tPi2         	   ;
 			tPi3_       		 =  tPi3         	   ;
 			JpsiV_     			 =  JpsiV 	  		   ;
			Ris_[0]    			 =  Ris[0]	  		   ;
			Ris_[1]    			 =  Ris[1]	  		   ;
			Ris_[2]    			 =  Ris[2]	  		   ;
			Ris_[3]    			 =  Ris[3]	  		   ;
			Ris_[4]    			 =  Ris[4]	  		   ;  
	        Els_	   			 =  Els		  		   ;
	        El_		   			 =  El		  		   ;
	        Cos_	   			 =  Cos		  		   ;
	        ClS_	   			 =  ClS		  		   ;
	        ClJpsi_	   			 =  ClJpsi	  		   ;
 		 	MatchMu_[0] 		 =  MatchMu[0]  	   ;
 		 	MatchMu_[1] 		 =  MatchMu[1]  	   ;
 		 	MatchPi1_[0]		 =  MatchPi1[0] 	   ;
 		 	MatchPi1_[1]		 =  MatchPi1[1] 	   ;
 		 	MatchPi2_[0]		 =  MatchPi2[0] 	   ;
 		 	MatchPi2_[1]		 =  MatchPi2[1] 	   ;
 		 	MatchPi3_[0]		 =  MatchPi3[0] 	   ;
 		 	MatchPi3_[1]		 =  MatchPi3[1] 	   ;
			MuCh_[0]    		 =  MuCh[0]	    	   ;
			MuCh_[1]    		 =  MuCh[1]	    	   ;
			ElsigJpsi_  		 =  ElsigJpsi  		   ; 
			CosJpsi_    		 =  CosJpsi    		   ; 
			
			for (int i=0; i < 3; i++)
			{
			  BS_[i]		   = 	BS[i]		       ;
			  VtxS_[i]		   = 	VtxS[i]		       ;
			  TrkNHits_[i]     =    TrkNHits[i]        ;      
			  TrkPixelHits_[i] =    TrkPixelHits[i]    ;      
			  TrkChi2_[i]      =    TrkChi2[i]     	   ;     
			  PiCh_[i]         =    PiCh[i]	           ;
			  IP3D_[i]   	   =    IP3D[i]	  		   ;
			  IPbs_[i]   	   =    IPbs[i]	  		   ;
			  sIP3D_[i]  	   =    sIP3D[i]    	   ;
			  sIPbs_[i]  	   =    sIPbs[i]   		   ;
			} 
			
			for(int j=0; j<9; j++)
			{
			  BSCovariance_[j] 	= 	BSCovariance[j]    ;
			  SVCovariance_[j] 	= 	SVCovariance[j]    ;
			}

 			for (int k=0; k < 4; k++)
 			{
 		  	  HltMatch_[k]		= 	HltMatch[k]		   ;
 			}  
			
}		        
//
//======================================================================
//
BcTreeCand::~BcTreeCand() {
}
//
//======================================================================
//
BcTreeCand *BcTree::AddBcTreeCand(    
   		    const TLorentzVector&	  Bc3Pi   ,
            const TLorentzVector&	  Mu1	  ,
            const TLorentzVector&	  Mu2	  ,
            const TLorentzVector&	  Jpsi    ,
            const TLorentzVector&	  Pi1	  ,
            const TLorentzVector&	  Pi2	  ,
            const TLorentzVector&	  Pi3	  ,
            const TLorentzVector&	  tPi1	  ,
            const TLorentzVector&	  tPi2	  ,
            const TLorentzVector&	  tPi3	  ,
            const TLorentzVector&	  JpsiV	  ,
	        Double_t     Ris[5]     		  ,
	        Double_t     El        		  	  ,
	        Double_t     Els	        	  ,
	        Double_t     Cos        		  ,
	        Double_t     ClS        		  ,
	        Double_t     ClJpsi     		  ,
	        Double_t     BS[3]    		      ,
	        Double_t     VtxS[3]    		  ,
 	        Double_t     BSCovariance[9]      ,
 	        Double_t     SVCovariance[9]      ,
		    Double_t	 IP3D[3]    		  ,
		    Double_t	 IPbs[3]    		  ,
		    Double_t	 sIP3D[3]   		  ,
		    Double_t	 sIPbs[3]   		  ,
		    Double_t	 MatchMu[2] 		  ,
		    Double_t	 MatchPi1[2]		  ,
		    Double_t	 MatchPi2[2]		  ,
		    Double_t	 MatchPi3[2]		  ,
	        Int_t	     PiCh[3]    		  ,
	        Int_t	     MuCh[2]    		  ,
	        Double_t 	 ElsigJpsi      	  ,    
	        Double_t 	 CosJpsi        	  ,   
     	    Int_t	     TrkNHits[3]          ,   
     	    Int_t	     TrkPixelHits[3]      , 
     	    Float_t	     TrkChi2[3]		      ,
     	    Double_t	 HltMatch[4]      

				   ){
   // Add a new track to the list of tracks for this event.
   // To avoid calling the very time consuming operator new for each track,
   // the standard but not well know C++ operator "new with placement"
   // is called. If tracks[i] is 0, a new Track object will be created
   // otherwise the previous Track[i] will be overwritten.

   TClonesArray &BcTreeArrayCand = *fBcTreeArrayCand ;
//   std::cout<<" NBcTreeCand = "<<NBcTreeCand<<std::endl;

   BcTreeCand *_BcTreeCand = new(BcTreeArrayCand[BcTree::NBcTreeCand++])  
   BcTreeCand(  Bc3Pi				,
   				Mu1     			,
				Mu2     			,
				Jpsi    			,
			    Pi1     			,
				Pi2     			,
				Pi3     			,
			    tPi1     			,
				tPi2     			,
				tPi3     			,
				JpsiV   			,
                Ris     			,    	     
 	          	El					,    	   
 	          	Els					,    	   
 	          	Cos					,    	   
 				ClS					,	     
 				ClJpsi				,	     
 		    	BS					,	  
 		    	VtxS				,	  
 		    	BSCovariance		,	  
 		    	SVCovariance		,	  
				IP3D    			,
				IPbs    			,
				sIP3D   			,
				sIPbs   			,
 				MatchMu 			,
 				MatchPi1 			,
 				MatchPi2 			,
 				MatchPi3 			,
	        	PiCh	 			,
	        	MuCh	 			,
	        	ElsigJpsi   		,      
				CosJpsi     		,
	        	TrkNHits    		,
	        	TrkPixelHits		,
				TrkChi2	    		,
				HltMatch
				);
 
    fIsValid = kTRUE;

//
//Save reference to last BcTreeCand in the collection of BcTreeCand
//
    if( Bc3Pi.Pt() > MaxPt )
    {
      MaxPt      = Bc3Pi.Pt() ;
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

//=================================================================================
//GENPARTICLES
//=================================================================================
//
void BcTree::SetBcTreeGENCand(    
	        TLorentzVector  Bc3Pi   ,
	        TLorentzVector  Mu1     ,
	        TLorentzVector  Mu2     ,
	        TLorentzVector  JPsi    ,
	        TLorentzVector  Pi1     ,
	        TLorentzVector  Pi2     ,
	        TLorentzVector  Pi3     ,
	        Double_t     	El      ,
	        Double_t     	Cos     ,
	        Double_t     	BS[3] ,
	        Double_t     	VtxS[3] ,
	        Double_t     	Ris[5]  ,
	        Int_t	     	PiCh[3] ,
	        Int_t	     	MuCh[2] 
				   ){

    fBcTreeGENCand.SetBc3Pi   (Bc3Pi   )  ;
    fBcTreeGENCand.SetMu1	  (Mu1     )  ;
	fBcTreeGENCand.SetMu2	  (Mu2     )  ;
	fBcTreeGENCand.SetJPsi    (JPsi    )  ;
	fBcTreeGENCand.SetPi1	  (Pi1     )  ;
	fBcTreeGENCand.SetPi2	  (Pi2     )  ;
	fBcTreeGENCand.SetPi3	  (Pi3     )  ;
 	fBcTreeGENCand.SetEl	  (El	   )  ;
 	fBcTreeGENCand.SetCos  	  (Cos     )  ;
 	fBcTreeGENCand.SetBSx	  (BS[0] )  ;  
 	fBcTreeGENCand.SetBSy	  (BS[1] )  ; 
 	fBcTreeGENCand.SetBSz	  (BS[2] )  ;
 	fBcTreeGENCand.SetVtxSx	  (VtxS[0] )  ;
 	fBcTreeGENCand.SetVtxSy	  (VtxS[1] )  ;
 	fBcTreeGENCand.SetVtxSz	  (VtxS[2] )  ;
	fBcTreeGENCand.SetRis0 	  (Ris[0]  )  ;
	fBcTreeGENCand.SetRis1 	  (Ris[1]  )  ;
	fBcTreeGENCand.SetRis2 	  (Ris[2]  )  ;
	fBcTreeGENCand.SetRis3 	  (Ris[3]  )  ;
	fBcTreeGENCand.SetRis4 	  (Ris[4]  )  ;
	fBcTreeGENCand.SetPiCh1	  (PiCh[0] )  ;
	fBcTreeGENCand.SetPiCh2	  (PiCh[1] )  ;
	fBcTreeGENCand.SetPiCh3	  (PiCh[2] )  ;
	fBcTreeGENCand.SetMuCh1	  (MuCh[0] )  ;
	fBcTreeGENCand.SetMuCh2	  (MuCh[1] )  ;
	 
}
//======================================================================
//JpsiCand 
//======================================================================
//
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
//
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
//
JpsiCand::~JpsiCand() {
}
//
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
