#ifndef BcTree_h
#define BcTree_h
#include "TObject.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "BcTreeCand.h"
#include "BcTreeJpsiCand.h"
#include "BcTreeGENCand.h"
#include "BcTreeHeader.h"

#include <TChain.h>
#include <iostream>
#include <TMath.h>

//####################################################################
//  BcTree Class
//####################################################################

class BcTree : public TObject {

private:
   TClonesArray  *fBcTreeArrayCand;            
   static TClonesArray *fgBcTreeArrayCand;

   TClonesArray  *fBcTreeJpsiArrayCand;            
   static TClonesArray *fgBcTreeJpsiArrayCand;

   BcTreeHeader   fBcTreeHeader		;
   BcTreeGENCand  fBcTreeGENCand	;

   Int_t          NBcTreeCand		;        
   Int_t          NJpsiCand			;          
   Int_t          NumJPsi  			;        
   Double_t       MaxPt				;
   Double_t       JPsiPxOld			;          
   Double_t       JPsiPyOld			;          
   Double_t       JPsiPzOld			;          

   TRefArray     *fMaxPtCand		;            //array of High Pt Cand Only

   Bool_t         fIsValid			;           
   Bool_t         jIsValid			;           

   void           ClearBcTreeHeader();
   void           ClearBcTreeGENCand();
   void           ClearJpsiCand();
   void           ClearBcCand();


public:
   BcTree();
   virtual       ~BcTree();

   void           BcTreeClear(Option_t *option ="");

   TClonesArray  *GetBcTreeArrayCand() 		const {	return fBcTreeArrayCand		;}
   TClonesArray  *GetBcTreeJpsiArrayCand() 	const {	return fBcTreeJpsiArrayCand	;}
   BcTreeHeader  *GetBcTreeHeader()          	  {	return &fBcTreeHeader  		;}
   BcTreeGENCand *GetBcTreeGENCand()         	  {	return &fBcTreeGENCand 		;}

   TRefArray     *GetMaxPtCand()       		const {	return fMaxPtCand      		;}
   Int_t          GetNBcTreeCand()     		const {	return NBcTreeCand     		;}
   Int_t          GetNJpsiCand()     		const {	return NJpsiCand       		;}
   Int_t          GetNumJPsi()         		const {	return NumJPsi         		;}

   BcTreeCand    *AddBcTreeCand(  const BcTreeCand&);

   void SetBcTreeGENCand(  
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

   void SetBcTreeHeader(	  
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

   ClassDef(BcTree,1)  //Event structure
};
#endif
