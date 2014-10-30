#include "TLorentzVector.h"
#include "TProcessID.h"
#include "TDirectory.h"
#include "../interface/BcStarTree.h"

ClassImp(BcStarTreeHeader)
ClassImp(BcStarTree)
ClassImp(BcStarCand)
ClassImp(BcCand)
ClassImp(BcStarTreeGENCand)
ClassImp(JpsiCand)

TClonesArray *BcStarTree::fgBcStarArrayCand         = 0;
TClonesArray *BcStarTree::fgBcArrayCand             = 0;
TClonesArray *BcStarTree::fgBcStarTreeJpsiArrayCand = 0;

//---------------------------------------------------------------------------------

BcStarTree::BcStarTree() : jIsValid(kFALSE)
{
    if (!fgBcArrayCand) 
	{
      std::cout<<"Create fBcArrayCand"<<std::endl;
      fgBcArrayCand = new TClonesArray("BcCand",1000);
    }  

    if (!fgBcStarArrayCand) 
	{
      std::cout<<"Create fBcStarArrayCand"<<std::endl;
      fgBcStarArrayCand = new TClonesArray("BcStarCand",1000);
    }  

    if (!fgBcStarTreeJpsiArrayCand) 
	{
      std::cout<<"Create fBcStarTreeJpsiArrayCand"<<std::endl;
      fgBcStarTreeJpsiArrayCand = new TClonesArray("JpsiCand",100);
    }  

    fBcStarTreeJpsiArrayCand = fgBcStarTreeJpsiArrayCand;
    NJpsiCand	= 0 ;

    fBcArrayCand = fgBcArrayCand;
    fMaxPtCand	 = new TRefArray;
    NBcCand	= 0 ;

    fBcStarArrayCand = fgBcStarArrayCand;
    NBcStarCand	= 0 ;

    MaxPt		= 0.;
    JPsiPxOld	= 0.;
    JPsiPyOld	= 0.;
    JPsiPzOld	= 0.;
    NumJPsi  	= 0 ;
 
//   TProcessID::SetObjectCount(ObjectNumber);
}

BcStarTree::~BcStarTree()
{
   BcStarTreeClear();
}

// ----------  Delete all the ntuple contents -----------------------

void BcStarTree::BcStarTreeClear(Option_t * /*option*/)
{
       fBcStarArrayCand->Clear("C");
       fBcArrayCand->Clear("C");
       fBcStarTreeJpsiArrayCand->Clear("C");

	   ClearBcStarTreeHeader()	;
	   ClearBcStarTreeGENCand()	;
       fMaxPtCand->Delete()	;

       NBcCand      = 0 ;
       NBcStarCand  = 0 ;
       NJpsiCand  	= 0 ;
       NumJPsi  	= 0 ;
       MaxPt	    = 0 ;
       JPsiPxOld	= 0 ;
       JPsiPyOld	= 0 ;
       JPsiPzOld	= 0 ;
 }

//Clear Header
void BcStarTree::ClearBcStarTreeHeader()
{
    fBcStarTreeHeader.SetNEve(    0	);
    fBcStarTreeHeader.SetRun(     0	);
    fBcStarTreeHeader.SetLBlk(    0	);
    fBcStarTreeHeader.SetEvent(   0	);
    fBcStarTreeHeader.SetNtrk(    0	);
    fBcStarTreeHeader.SetNprim(   0	);
    fBcStarTreeHeader.SetNGoodTrk(0	);
    fBcStarTreeHeader.SetTrueNI(  0 );

	for (int i = 0; i < 3; i++){
 	  fBcStarTreeHeader.SetBS( 			i, 0);
 	  fBcStarTreeHeader.SetHighestPtPV( i, 0);
	}
	for (int i = 0; i < 9; i++){
 	  fBcStarTreeHeader.SetBSCovariance(  i, 0);
 	  fBcStarTreeHeader.SetHighestPtPVCov(i, 0);
	}
    fBcStarTreeHeader.SetNumBc(   0 );

}

//Clear GEN Cand
void BcStarTree::ClearBcStarTreeGENCand()
{
    TLorentzVector *empty = new TLorentzVector (0,0,0,0);
    fBcStarTreeGENCand.SetBcPi    (*empty)  ;
    fBcStarTreeGENCand.SetMuP	  (*empty)  ;
	fBcStarTreeGENCand.SetMuM	  (*empty)  ;
	fBcStarTreeGENCand.SetJPsi    (*empty)  ;
	fBcStarTreeGENCand.SetPi	  (*empty)  ;
 	fBcStarTreeGENCand.SetEl	  ( 0 )  	;
 	fBcStarTreeGENCand.SetCos  	  ( 0 )  	;
 	fBcStarTreeGENCand.SetTau  	  ( 0 )  	;
 	fBcStarTreeGENCand.SetPiCh 	  ( 0 )  	;

	for (int i = 0; i < 3; i++){
 	  fBcStarTreeGENCand.SetPVtx	( i, 0 );
 	  fBcStarTreeGENCand.SetBcVtx	( i, 0 );
 	  fBcStarTreeGENCand.SetJpsiVtx	( i, 0 );
 	  fBcStarTreeGENCand.SetBS		( i, 0 );
	}
}

//------------------ Set Header -------------------------------

void BcStarTree::SetBcStarTreeHeader( 
   				  UInt_t    NEve   			 ,
 				  UInt_t    Run    			 , 
 				  UInt_t    LBlk   			 ,
 				  UInt_t    Event  			 ,
				  Int_t     Ntrk   			 ,
				  Int_t     Nprim  			 ,
				  Int_t   	NGoodTrk		 ,
				  Float_t   TrueNI 			 ,
				  Double_t  BS[3]			 , 
				  Double_t  BSCovariance[9]	 , 
				  Double_t  HighestPtPV[3]	 , 
				  Double_t  HighestPtPVCov[9], 
 				  UInt_t    NumBc  
				  )
{
    fBcStarTreeHeader.SetNEve    ( NEve		);
    fBcStarTreeHeader.SetRun	 ( Run		);
    fBcStarTreeHeader.SetLBlk	 ( LBlk		);
    fBcStarTreeHeader.SetEvent	 ( Event	);
    fBcStarTreeHeader.SetNtrk	 ( Ntrk		);
    fBcStarTreeHeader.SetNprim	 ( Nprim	);
    fBcStarTreeHeader.SetNGoodTrk( NGoodTrk	);
    fBcStarTreeHeader.SetTrueNI	 ( TrueNI 	);

	for (int i = 0; i < 3; i++){
 	  fBcStarTreeHeader.SetBS		  ( i, BS[i]		 );
 	  fBcStarTreeHeader.SetHighestPtPV( i, HighestPtPV[i]);
	}
	for (int i = 0; i < 9; i++){
 	  fBcStarTreeHeader.SetBSCovariance  ( i, BSCovariance[i]  );
 	  fBcStarTreeHeader.SetHighestPtPVCov( i, HighestPtPVCov[i]);
	}
    fBcStarTreeHeader.SetNumBc(   NumBc);
}

//------------------ Set GEN Cand -------------------------------
void BcStarTree::SetBcStarTreeGENCand(    
 				  TLorentzVector  BcPi     	,
 				  TLorentzVector  MuP     	,
 				  TLorentzVector  MuM     	,
 				  TLorentzVector  JPsi    	,
 				  TLorentzVector  Pi	  	,
 				  Double_t        El	  	,
 				  Double_t        Cos	  	,
 				  Double_t        Tau	  	,
 				  Double_t        PVtx[3] 	,
 				  Double_t        BcVtx[3] 	,
 				  Double_t        JpsiVtx[3],
 				  Double_t        BS[3]     ,
 				  Int_t			  PiCh
				   )
{
    fBcStarTreeGENCand.SetBcPi    (BcPi    )  ;
    fBcStarTreeGENCand.SetMuP	  (MuP     )  ;
	fBcStarTreeGENCand.SetMuM	  (MuM     )  ;
	fBcStarTreeGENCand.SetJPsi    (JPsi    )  ;
	fBcStarTreeGENCand.SetPi	  (Pi      )  ;
 	fBcStarTreeGENCand.SetEl	  (El	   )  ;
 	fBcStarTreeGENCand.SetCos  	  (Cos     )  ;
 	fBcStarTreeGENCand.SetTau  	  (Tau     )  ;
	fBcStarTreeGENCand.SetPiCh	  (PiCh    )  ;

	for (int i=0; i < 3; i++)
	{
   	  fBcStarTreeGENCand.SetPVtx	  ( i, PVtx[i] 		);
   	  fBcStarTreeGENCand.SetBcVtx	  ( i, BcVtx[i] 	);
   	  fBcStarTreeGENCand.SetJpsiVtx	  ( i, JpsiVtx[i] 	);
   	  fBcStarTreeGENCand.SetBS	  	  ( i, BS[i] 		);
	}  
	 
}

//======================================================================
//
BcCand *BcStarTree::AddBcCand(    
   		    const BcCand&	myBcCand  ){

   TClonesArray &BcArrayCand = *fBcArrayCand ;
   BcCand *_BcCand = new(BcArrayCand[BcStarTree::NBcCand++])  
   BcCand( myBcCand);
 
   fIsValid = kTRUE;
   return _BcCand;
}

//======================================================================

BcStarCand *BcStarTree::AddBcStarCand(    
   		    const BcStarCand&	myBcStarCand  ){

   TClonesArray &BcStarArrayCand = *fBcStarArrayCand ;
   BcStarCand *_BcStarCand = new(BcStarArrayCand[BcStarTree::NBcStarCand++])  
   BcStarCand( myBcStarCand);
 
   gIsValid = kTRUE;
   return _BcStarCand;
}

//======================================================================
JpsiCand *BcStarTree::AddJpsiCand(    
            	  const TLorentzVector&	 	MuP   		    ,
            	  const TLorentzVector&	 	MuM   		    ,
            	  const TLorentzVector&	 	Jpsi  		    ,
 				  Double_t     				Els     	    ,
 				  Double_t    				Cos     		,
 				  Double_t    				ClS     		,
 				  Double_t    				VtxJ[3] 		,
				  Double_t    				MatchMu[2]      
				   ){

   TClonesArray &BcStarTreeJpsiArrayCand = *fBcStarTreeJpsiArrayCand ;

   JpsiCand *_JpsiCand = new(BcStarTreeJpsiArrayCand[BcStarTree::NJpsiCand++])  
   JpsiCand(    MuP     	,
				MuM     	,
				Jpsi    	,
 	          	Els			,    	   
 				Cos			,	     
 				ClS			,	     
 		    	VtxJ		,	  
 				MatchMu 	
	        	);  
 
    jIsValid = kTRUE;

   return _JpsiCand;
}
