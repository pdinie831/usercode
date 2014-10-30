#ifndef BcStarTree_h
#define BcStarTree_h
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "BcCand.h"
#include "BcStarCand.h"
#include "BcStarTreeJpsiCand.h"
#include "BcStarTreeGENCand.h"
#include "BcStarTreeHeader.h"

#include <TChain.h>
#include <iostream>
#include <TMath.h>

//####################################################################
//  BcStarTree Class
//####################################################################

class BcStarTree : public TObject {

private:
   TClonesArray  *fBcArrayCand;            
   static TClonesArray *fgBcArrayCand;

   TClonesArray  *fBcStarArrayCand;            
   static TClonesArray *fgBcStarArrayCand;

   TClonesArray  *fBcStarTreeJpsiArrayCand;            
   static TClonesArray *fgBcStarTreeJpsiArrayCand;

   BcStarTreeHeader   fBcStarTreeHeader		;
   BcStarTreeGENCand  fBcStarTreeGENCand	;

   Int_t          NBcCand           ;         
   Int_t          NBcStarCand       ;        
   Int_t          NJpsiCand			;          
   Int_t          NumJPsi  			;        
   Double_t       MaxPt				;
   Double_t       JPsiPxOld			;          
   Double_t       JPsiPyOld			;          
   Double_t       JPsiPzOld			;          

   TRefArray     *fMaxPtCand		;            //array of High Pt Cand Only

   Bool_t         fIsValid			;           
   Bool_t         gIsValid			;           
   Bool_t         jIsValid			;           

   void           ClearBcStarTreeHeader();
   void           ClearBcStarTreeGENCand();
   void           ClearJpsiCand();
   void           ClearBcCand();
   void           ClearBcStarCand();
  

public:
   BcStarTree();
   virtual       ~BcStarTree();

   void           BcStarTreeClear(Option_t *option ="");

   TClonesArray  *GetBcArrayCand() 		        const {	return fBcArrayCand	            ;}
   TClonesArray  *GetBcStarArrayCand() 		    const {	return fBcStarArrayCand		    ;}
   TClonesArray  *GetBcStarTreeJpsiArrayCand() 	const {	return fBcStarTreeJpsiArrayCand	;}
   BcStarTreeHeader  *GetBcStarTreeHeader()           {	return &fBcStarTreeHeader  		;}
   BcStarTreeGENCand *GetBcStarTreeGENCand()          {	return &fBcStarTreeGENCand 		;}

   TRefArray     *GetMaxPtCand()       		const {	return fMaxPtCand      		;}
   Int_t          GetNBcCand()              const {	return NBcCand     	     	;}
   Int_t          GetNBcStarCand()          const {	return NBcStarCand     	    ;}
   Int_t          GetNJpsiCand()     		const {	return NJpsiCand       		;}
   Int_t          GetNumJPsi()         		const {	return NumJPsi         		;}

   BcCand        *AddBcCand     (  const BcCand&     );
   BcStarCand    *AddBcStarCand (  const BcStarCand& );

   void SetBcStarTreeGENCand(  
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
 				  Double_t        BS[3]		,   
 				  Int_t        	  PiCh   
				);

   void SetBcStarTreeHeader(	  
   				  UInt_t    NEve   			 ,
 				  UInt_t    Run    			 , 
 				  UInt_t    LBlk   			 ,
 				  UInt_t    Event  			 ,
				  Int_t     Ntrk   			 ,
				  Int_t     Nprim  			 ,
				  Int_t   	NGoodTrk		 ,
				  Float_t   TrueNI 			 ,
				  Double_t  BS[3]			 , 
				  Double_t  BSCovariance[3]	 , 
				  Double_t  HighestPtPV[3]	 , 
				  Double_t  HighestPtPVCov[3], 
 				  UInt_t    NumBc  
				);


   JpsiCand    	 *AddJpsiCand(  
            	  const TLorentzVector&	 	MuP   		    ,
            	  const TLorentzVector&	 	MuM   		    ,
            	  const TLorentzVector&	 	Jpsi  		    ,
 				  Double_t     				Els     	    ,
 				  Double_t    				Cos     		,
 				  Double_t    				ClS     		,
 				  Double_t    				VtxJ[3] 		,
				  Double_t    				MatchMu[2]      
				);

   ClassDef(BcStarTree,1)  //Event structure
};
#endif
