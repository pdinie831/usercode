#include "TLorentzVector.h"
#include "TProcessID.h"
#include "TDirectory.h"
#include "../interface/BcTree.h"

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

// ----------  Delete all the ntuple contents -----------------------

void BcTree::BcTreeClear(Option_t * /*option*/)
{
       fBcTreeArrayCand->Clear("C");
       fBcTreeJpsiArrayCand->Clear("C");

	   ClearBcTreeHeader()	;
	   ClearBcTreeGENCand()	;
       fMaxPtCand->Delete()	;

       NBcTreeCand  = 0 ;
       NJpsiCand  	= 0 ;
       NumJPsi  	= 0 ;
       MaxPt	    = 0 ;
       JPsiPxOld	= 0 ;
       JPsiPyOld	= 0 ;
       JPsiPzOld	= 0 ;
 }

//Clear Header
void BcTree::ClearBcTreeHeader()
{
    fBcTreeHeader.SetNEve(    0	);
    fBcTreeHeader.SetRun(     0	);
    fBcTreeHeader.SetLBlk(    0	);
    fBcTreeHeader.SetEvent(   0	);
    fBcTreeHeader.SetNtrk(    0	);
    fBcTreeHeader.SetNprim(   0	);
    fBcTreeHeader.SetNGoodTrk(0	);
    fBcTreeHeader.SetTrueNI(  0 );

	for (int i = 0; i < 3; i++){
 	  fBcTreeHeader.SetBS( 			i, 0);
 	  fBcTreeHeader.SetHighestPtPV( i, 0);
	}
	for (int i = 0; i < 9; i++){
 	  fBcTreeHeader.SetBSCovariance(  i, 0);
 	  fBcTreeHeader.SetHighestPtPVCov(i, 0);
	}
    fBcTreeHeader.SetNumBc(   0 );

}

//Clear GEN Cand
void BcTree::ClearBcTreeGENCand()
{
    TLorentzVector *empty = new TLorentzVector (0,0,0,0);
    fBcTreeGENCand.SetBcPi    (*empty)  ;
    fBcTreeGENCand.SetMuP	  (*empty)  ;
	fBcTreeGENCand.SetMuM	  (*empty)  ;
	fBcTreeGENCand.SetJPsi    (*empty)  ;
	fBcTreeGENCand.SetPi	  (*empty)  ;
 	fBcTreeGENCand.SetEl	  ( 0 )  	;
 	fBcTreeGENCand.SetCos  	  ( 0 )  	;
 	fBcTreeGENCand.SetTau  	  ( 0 )  	;
 	fBcTreeGENCand.SetPiCh 	  ( 0 )  	;

	for (int i = 0; i < 3; i++){
 	  fBcTreeGENCand.SetPVtx	( i, 0 );
 	  fBcTreeGENCand.SetBcVtx	( i, 0 );
 	  fBcTreeGENCand.SetJpsiVtx	( i, 0 );
 	  fBcTreeGENCand.SetBS		( i, 0 );
	}
}

//------------------ Set Header -------------------------------

void BcTree::SetBcTreeHeader( 
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
    fBcTreeHeader.SetNEve    ( NEve		);
    fBcTreeHeader.SetRun	 ( Run		);
    fBcTreeHeader.SetLBlk	 ( LBlk		);
    fBcTreeHeader.SetEvent	 ( Event	);
    fBcTreeHeader.SetNtrk	 ( Ntrk		);
    fBcTreeHeader.SetNprim	 ( Nprim	);
    fBcTreeHeader.SetNGoodTrk( NGoodTrk	);
    fBcTreeHeader.SetTrueNI	 ( TrueNI 	);

	for (int i = 0; i < 3; i++){
 	  fBcTreeHeader.SetBS		  ( i, BS[i]		 );
 	  fBcTreeHeader.SetHighestPtPV( i, HighestPtPV[i]);
	}
	for (int i = 0; i < 9; i++){
 	  fBcTreeHeader.SetBSCovariance  ( i, BSCovariance[i]  );
 	  fBcTreeHeader.SetHighestPtPVCov( i, HighestPtPVCov[i]);
	}
    fBcTreeHeader.SetNumBc(   NumBc);
}

//------------------ Set GEN Cand -------------------------------
void BcTree::SetBcTreeGENCand(    
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
    fBcTreeGENCand.SetBcPi    (BcPi    )  ;
    fBcTreeGENCand.SetMuP	  (MuP     )  ;
	fBcTreeGENCand.SetMuM	  (MuM     )  ;
	fBcTreeGENCand.SetJPsi    (JPsi    )  ;
	fBcTreeGENCand.SetPi	  (Pi      )  ;
 	fBcTreeGENCand.SetEl	  (El	   )  ;
 	fBcTreeGENCand.SetCos  	  (Cos     )  ;
 	fBcTreeGENCand.SetTau  	  (Tau     )  ;
	fBcTreeGENCand.SetPiCh	  (PiCh    )  ;

	for (int i=0; i < 3; i++)
	{
   	  fBcTreeGENCand.SetPVtx	  ( i, PVtx[i] 		);
   	  fBcTreeGENCand.SetBcVtx	  ( i, BcVtx[i] 	);
   	  fBcTreeGENCand.SetJpsiVtx	  ( i, JpsiVtx[i] 	);
   	  fBcTreeGENCand.SetBS	  	  ( i, BS[i] 		);
	}  
	 
}

//======================================================================
//
BcTreeCand *BcTree::AddBcTreeCand(    
   		    const BcTreeCand&	myBcTreeCand  ){

   TClonesArray &BcTreeArrayCand = *fBcTreeArrayCand ;
   BcTreeCand *_BcTreeCand = new(BcTreeArrayCand[BcTree::NBcTreeCand++])  
   BcTreeCand( myBcTreeCand);
 
    fIsValid = kTRUE;

//
//Save reference to last BcTreeCand in the collection of BcTreeCand
//
//     Double_t mxBcPt = *myBcTreeCand->GetBcPi()->Pt();
//     if( BcPi.Pt() > MaxPt )
// 	{
//       MaxPt      = mxBcPt ;
//       fMaxPtCand -> Add( _BcTreeCand );
//     }
//     
// 	Double_t JPsiPx = JpsiV.Px();
//     Double_t JPsiPy = JpsiV.Py();
//     Double_t JPsiPz = JpsiV.Pz();
//     if(JPsiPx!=JPsiPxOld&&JPsiPy!=JPsiPyOld&&JPsiPz!=JPsiPzOld)
// 	{
// 	  JPsiPxOld = JPsiPx;
//       JPsiPyOld = JPsiPy;
//       JPsiPzOld = JPsiPz;
//       NumJPsi++;
//     } 

   return _BcTreeCand;
}

//======================================================================
//GENPARTICLES
//======================================================================


//
//======================================================================
JpsiCand *BcTree::AddJpsiCand(    
            	  const TLorentzVector&	 	MuP   		    ,
            	  const TLorentzVector&	 	MuM   		    ,
            	  const TLorentzVector&	 	Jpsi  		    ,
 				  Double_t     				Els     	    ,
 				  Double_t    				Cos     		,
 				  Double_t    				ClS     		,
 				  Double_t    				VtxJ[3] 		,
				  Double_t    				MatchMu[2]      
				   ){

   TClonesArray &BcTreeJpsiArrayCand = *fBcTreeJpsiArrayCand ;

   JpsiCand *_JpsiCand = new(BcTreeJpsiArrayCand[BcTree::NJpsiCand++])  
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
