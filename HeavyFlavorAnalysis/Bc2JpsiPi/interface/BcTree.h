#ifndef BcTree_h
#define BcTree_h
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <TMath.h>
#include "TRefArray.h"
#include "TRef.h"

#include <TChain.h>
#include <iostream>

//################################################
//  BcTreeCand Class
//################################################


class BcTreeCand : public TObject {

private:

 	TLorentzVector BcPi_ 	     ;
 	TLorentzVector BcK_ 	     ;
 	TLorentzVector Mu1_   	     ;
 	TLorentzVector Mu2_   	     ;
 	TLorentzVector Jpsi_  	     ;
 	TLorentzVector Pi_   	     ;
 	TLorentzVector JpsiV_  	     ;
 	Double_t     Els_     	     ;
 	Double_t     Cos_     	     ;
 	Double_t     ClS_     	     ;
 	Double_t     ClJpsi_  	     ;
 	Double_t     BS_[3] 	     ;
 	Double_t     VtxS_[3] 	     ;
    Double_t     BSCovariance_[9];
  	Double_t     SVCovariance_[9];
	Double_t     IP3D_ 	         ;
	Double_t     IPbs_ 	         ;
	Double_t     sIP3D_	         ;
	Double_t     sIPbs_	         ;
 	Double_t	 MatchMu_[2]     ;   
 	Double_t     MatchPi_        ;   
 	Int_t	     PiCh_           ;   
 	Int_t	     MuCh_[2]        ;   
    Double_t     ElsigJpsi_      ;
	Double_t     CosJpsi_        ;
 	Int_t 		 TrkNHits_       ;
 	Int_t  		 TrkPixelHits_   ;
 	Float_t		 TrkChi2_        ;
  	Double_t     HltMatch_[4]	 ;

 
//      
public:
//
//
        BcTreeCand(){;}
        BcTreeCand(const BcTreeCand& orig);
        BcTreeCand( 	
            const TLorentzVector&     BcPi  ,
            const TLorentzVector&     BcK   ,
            const TLorentzVector&     Mu1   ,
            const TLorentzVector&     Mu2   ,
            const TLorentzVector&     Jpsi  ,
            const TLorentzVector&     Pi    ,
            const TLorentzVector&     JpsiV ,
 			Double_t     Els            ,
 			Double_t     Cos            ,
 			Double_t     ClS            ,
 			Double_t     ClJpsi         ,
 			Double_t     BS[3]          ,
 			Double_t     VtxS[3]        ,
 		 	Double_t     BSCovariance[9],
 		 	Double_t     SVCovariance[9],
			Double_t     IP3D           ,
			Double_t     IPbs           ,
			Double_t     sIP3D          ,
			Double_t     sIPbs          ,
			Double_t     MatchMu[2]     ,
			Double_t     MatchPi        ,
 			Int_t	     PiCh           ,
 			Int_t	     MuCh[2]        ,
          	Double_t     ElsigJpsi      ,
          	Double_t     CosJpsi        ,
 		 	Int_t        TrkNHits       , 
 		 	Int_t        TrkPixelHits   ,  
 		 	Float_t      TrkChi2 		,     
 		 	Double_t     HltMatch[4]      
                   );

	virtual ~BcTreeCand();

 	const TLorentzVector   *GetBcPi()	     	 	 	{return  &BcPi_				;}
 	const TLorentzVector   *GetBcK()	     	 		{return  &BcK_		    	;}
 	const TLorentzVector   *GetMu1()	     	 	 	{return  &Mu1_				;}
 	const TLorentzVector   *GetMu2()	     	 	 	{return  &Mu2_				;}
 	const TLorentzVector   *GetJpsi()	     	 	 	{return  &Jpsi_				;}
 	const TLorentzVector   *GetPi()	         	 	 	{return  &Pi_ 				;}
 	const TLorentzVector   *GetJpsiV()	     	 	 	{return  &JpsiV_			;}
 	Double_t    			GetClS()			 	 	{return  ClS_ 				;}
 	Double_t    			GetEls()			 	 	{return  Els_ 				;}
 	Double_t    			GetCos()			 	 	{return  Cos_ 				;}
 	Double_t    			GetClJpsi() 		 	 	{return  ClJpsi_			;}
 	Double_t    			GetBS(Int_t Id_)   	 		{return  BS_[Id_]			;}
 	Double_t    			GetVtxS(Int_t Id_)   	 	{return  VtxS_[Id_]			;}
 	Double_t				GetBSCovariance(Int_t Id_)	{return  BSCovariance_[Id_] ;}
 	Double_t				GetSVCovariance(Int_t Id_)	{return  SVCovariance_[Id_] ;}
	Double_t    			GetIP3D()			 	 	{return  IP3D_				;}
	Double_t    			GetIPbs()			 		{return  IPbs_		 		;}
	Double_t    			GetsIP3D()  		 	 	{return  sIP3D_				;}
	Double_t    			GetsIPbs()  		 	 	{return  sIPbs_				;}
	Double_t    			GetMatchMu(Int_t Id_)	 	{return  MatchMu_[Id_] 		;}
	Double_t    			GetMatchPi()		 	 	{return  MatchPi_ 			;}
	Int_t	    			GetPiCh()			 	 	{return  PiCh_				;}  	
	Int_t	   				GetMuCh(Int_t Id_)   	 	{return  MuCh_[Id_]    		;}	   
	Double_t    			GetElsigJpsi()       	 	{return  ElsigJpsi_    		;}
	Double_t    			GetCosJpsi()         	 	{return  CosJpsi_      		;}
	Int_t	    			GetTrkNHits()		 	 	{return  TrkNHits_			;}  	
	Int_t	    			GetTrkPixelHIts()	 	 	{return  TrkPixelHits_	 	;}  	
	Float_t	    			GetTrkChi2()		 	 	{return  TrkChi2_			;}  	
	Double_t    			GetHltMatch(Int_t Id_)	 	{return  HltMatch_[Id_] 	;}

        ClassDef(BcTreeCand,1) 
};


//##################################################################################################################
//  BcTreeGENCand Class
//##################################################################################################################


class BcTreeGENCand : public TObject {

private:

	TLorentzVector     BcPi_    ;
 	TLorentzVector     Mu1_     ;
 	TLorentzVector     Mu2_     ;
 	TLorentzVector     JPsi_    ;
 	TLorentzVector     Pi_      ;
  	Double_t     	   El_  	;
 	Double_t     	   Cos_ 	;
 	Double_t     	   VtxS_[3] ;
 	Double_t     	   BS_[3]   ;
 	Int_t	     	   PiCh_ 	;	
 	Int_t	    	   MuCh_[2] ;   
 

public:

    BcTreeGENCand(){;}
    virtual ~BcTreeGENCand(){ ;}


  	void   SetBcPi            (TLorentzVector n) {   BcPi_    = n ;}
    void   SetMu1	      	  (TLorentzVector n) {   Mu1_	  = n ;}
	void   SetMu2	      	  (TLorentzVector n) {   Mu2_	  = n ;}
	void   SetJPsi	          (TLorentzVector n) {   JPsi_    = n ;}
	void   SetPi	      	  (TLorentzVector n) {   Pi_	  = n ;}
 	void   SetEl	      	  (Double_t n)  	 {   El_	  = n ;}
 	void   SetCos	      	  (Double_t n)  	 {   Cos_	  = n ;}
 	void   SetVtxSx           (Double_t n)  	 {   VtxS_[0] = n ;}
 	void   SetVtxSy           (Double_t n)  	 {   VtxS_[1] = n ;}
 	void   SetVtxSz           (Double_t n)  	 {   VtxS_[2] = n ;}
 	void   SetBSx             (Double_t n)  	 {   BS_[0]   = n ;}
 	void   SetBSy             (Double_t n)  	 {   BS_[1]   = n ;}
 	void   SetBSz             (Double_t n)  	 {   BS_[2]   = n ;}
	void   SetPiCh    		  (Int_t n)	    	 {   PiCh_    = n ;}	
	void   SetMuCh1    		  (Int_t n)	         {   MuCh_[0] = n ;}	
	void   SetMuCh2    	  	  (Int_t n)	         {   MuCh_[1] = n ;}	

 	TLorentzVector   GetBcPi()      	{return  BcPi_  	;}
    TLorentzVector   GetMu1()	    	{return  Mu1_		;}
	TLorentzVector   GetMu2()	    	{return  Mu2_		;}
	TLorentzVector   GetJPsi()	    	{return  JPsi_  	;}
	TLorentzVector   GetPi()	    	{return  Pi_		;}
 	Double_t    	 GetEl()			{return  El_		;}
 	Double_t    	 GetCos()			{return  Cos_		;}
 	Double_t    	 GetVtxS(Int_t Id_) {return  VtxS_[Id_] ;}
 	Double_t    	 GetBS(Int_t Id_)   {return  BS_[Id_]   ;}
	Int_t	    	 GetPiCh()  		{return  PiCh_  	;}     
	Int_t	         GetMuCh(Int_t Id_) {return  MuCh_[Id_] ;}	   

    ClassDef(BcTreeGENCand,1) 
};


//##################################################################################################################
//  JpsiCand Class
//##################################################################################################################


class JpsiCand : public TObject {

private:

 	TLorentzVector     Mu1_     	;
 	TLorentzVector     Mu2_     	;
 	TLorentzVector     Jpsi_    	;
  	Double_t     	   Els_  		;
 	Double_t     	   Cos_ 		;
 	Double_t     	   ClS_ 		;
 	Double_t     	   VtxJ_[3] 	;
 	Double_t	 	   MatchMu_[2]  ;   
 	Int_t	    	   MuCh_[2] 	;   
 

public:

    JpsiCand(){;}
    JpsiCand(const JpsiCand& orig);
    JpsiCand( 	
            const TLorentzVector&     Mu1   ,
            const TLorentzVector&     Mu2   ,
            const TLorentzVector&     Jpsi  ,
 			Double_t     Els            ,
 			Double_t     Cos            ,
 			Double_t     ClS            ,
 			Double_t     VtxJ[3]        ,
			Double_t     MatchMu[2]     ,
 			Int_t	     MuCh[2]        
                   );

    virtual ~JpsiCand();

    const TLorentzVector   *GetMu1()	    		{return  &Mu1_		;}
	const TLorentzVector   *GetMu2()	    		{return  &Mu2_		;}
	const TLorentzVector   *GetJPsi()	    		{return  &Jpsi_  	;}
 	Double_t    	 		GetEls()				{return  Els_		;}
 	Double_t    	 		GetCos()				{return  Cos_		;}
 	Double_t    	 		GetClS()				{return  ClS_		;}
 	Double_t    	 		GetVtx(Int_t Id_)  		{return  VtxJ_[Id_] ;}
 	Double_t    	 		GetMatchMu(Int_t Id_)  	{return  MatchMu_[Id_] ;}
	Int_t	         		GetMuCh(Int_t Id_) 		{return  MuCh_[Id_] ;}	   

    ClassDef(JpsiCand,1) 
};


//################################################################################################
//  BcTreeHeader Class
//################################################################################################

class BcTreeHeader {

private:
   UInt_t       NEve_  ;
   UInt_t       Run_   ;
   UInt_t       LBlk_  ;
   UInt_t       Event_ ;
   UInt_t       NumBc_ ;
   Int_t        Ntrk_  ;
   Int_t        Nprim_ ;
   Float_t      TrueNI_;
   Int_t        TriggerBit_ ;
   Int_t        GoodTrkSize_ ;

public:
   BcTreeHeader() { ;}
   virtual ~BcTreeHeader()  { ;}
   void	   SetNEve       (UInt_t    n){ NEve_   = n;}
   void	   SetRun        (UInt_t    n){ Run_    = n;}
   void	   SetLBlk       (UInt_t    n){ LBlk_   = n;}  
   void	   SetEvent      (UInt_t    n){ Event_  = n;}
   void	   SetNumBc      (UInt_t    n){ NumBc_  = n;}	
   void    SetNtrk       (Int_t     n){ Ntrk_   = n;}
   void    SetNprim      (Int_t     n){ Nprim_  = n;}
   void    SetTrueNI     (Float_t   n){ TrueNI_ = n;}
   void    SetTriggerBit (Int_t     n){ TriggerBit_  = n;}
   void    SetGoodTrkSize(Int_t     n){ GoodTrkSize_ = n;}
   UInt_t  GetNEve()   {return  NEve_  ;}
   UInt_t  GetRun()    {return  Run_   ;}
   UInt_t  GetLBlk()   {return  LBlk_  ;}
   UInt_t  GetEvent()  {return  Event_ ;}
   UInt_t  GetNumBc()  {return  NumBc_ ;}
   Int_t   GetNtrk()   {return  Ntrk_  ;}
   Int_t   GetNprim()  {return  Nprim_ ;}
   Float_t GetTrueNI() {return  TrueNI_;}
   Int_t   GetTriggerBit()  {return  TriggerBit_ ;}
   Int_t   GetGoodTrkSize() {return  GoodTrkSize_;}
 
   ClassDef(BcTreeHeader,1)  //Event Header
};

//################################################################################################
//  BcTree Class
//################################################################################################

class BcTree : public TObject {

private:
   TClonesArray  *fBcTreeArrayCand;            
   static TClonesArray *fgBcTreeArrayCand;

   TClonesArray  *fBcTreeJpsiArrayCand;            
   static TClonesArray *fgBcTreeJpsiArrayCand;

   BcTreeHeader   fBcTreeHeader;
   BcTreeGENCand  fBcTreeGENCand;
   Int_t          NBcTreeCand;        //Number of Bc Candidates
   Int_t          NJpsiCand;          //Number of Jpsi Candidates
   Double_t       MaxPt;
   Double_t       JPsiPxOld;          // Count JPsi 
   Double_t       JPsiPyOld;          // ...
   Double_t       JPsiPzOld;          // ...
   int            NumJPsi  ;          // Count JPsi

   TRefArray     *fMaxPtCand;            //array of High Pt Cand Only

   Bool_t         fIsValid;           //
   Bool_t         jIsValid;           //
   void           ClearBcTreeHeader();
   void           ClearBcTreeGENCand();
   void           ClearJpsiCand();


public:
   BcTree();
   virtual       ~BcTree();
   void           BcTreeClear(Option_t *option ="");
   TClonesArray  *GetBcTreeArrayCand() 		const {return fBcTreeArrayCand;}
   TClonesArray  *GetBcTreeJpsiArrayCand() 	const {return fBcTreeJpsiArrayCand;}
   TRefArray     *GetMaxPtCand()       		const {return fMaxPtCand      ;}
   Int_t          GetNBcTreeCand()     		const {return NBcTreeCand     ;}
   Int_t          GetNJpsiCand()     		const {return NJpsiCand       ;}
   Int_t          GetNumJPsi()         		const {return NumJPsi         ;}
   BcTreeHeader  *GetBcTreeHeader()          	  {return &fBcTreeHeader  ;}
   BcTreeGENCand *GetBcTreeGENCand()         	  {return &fBcTreeGENCand ;}

   BcTreeCand    *AddBcTreeCand(  
   				  const TLorentzVector&  	BcPi  		    ,
   				  const TLorentzVector&  	BcK   		    ,
            	  const TLorentzVector&	 	Mu1   		    ,
            	  const TLorentzVector&	 	Mu2   		    ,
            	  const TLorentzVector&	 	Jpsi  		    ,
            	  const TLorentzVector&	 	Pi    		    ,
            	  const TLorentzVector&	 	JpsiV 		    ,
 				  Double_t     				Els     	    ,
 				  Double_t    				Cos     		,
 				  Double_t    				ClS     		,
 				  Double_t    				ClJpsi  		,
 				  Double_t    				BS[3] 		    ,
 				  Double_t    				VtxS[3] 		,
				  Double_t    				BSCovariance[9] ,
				  Double_t    				SVCovariance[9] ,
				  Double_t    				IP3D            ,
				  Double_t    				IPbs            ,
				  Double_t    				sIP3D  		    ,
				  Double_t    				sIPbs   		,
				  Double_t    				MatchMu[2]      ,
				  Double_t    				MatchPi         ,
 				  Int_t       				PiCh    		, 
 				  Int_t       				MuCh[2] 		,
 				  Double_t    				ElsigJpsi	    ,    
 				  Double_t    				CosJpsi 	  	,     
 				  Int_t       				TrkNHits		,  
 				  Int_t       				TrkPixelHits    , 
 				  Float_t       			TrkChi2         ,
				  Double_t    				HltMatch[4]      
				);

   void SetBcTreeGENCand(  
 				  TLorentzVector  BcPi    ,
 				  TLorentzVector  Mu1     ,
 				  TLorentzVector  Mu2     ,
 				  TLorentzVector  JPsi    ,
 				  TLorentzVector  Pi	  ,
 				  Double_t        El	  ,
 				  Double_t        Cos	  ,
 				  Double_t        VtxS[3] ,
 				  Double_t        BS[3]   ,
 				  Int_t           PiCh    ,
 				  Int_t           MuCh[2] 
				);

   void SetBcTreeHeader(	  
   				  UInt_t     NEve   ,
 				  UInt_t     Run    ,
 				  UInt_t     LBlk   ,
 				  UInt_t     Event  ,
 				  UInt_t     NumBc  ,
				  Int_t      Ntrk   ,
				  Int_t      Nprim  ,
				  Float_t    TrueNI ,
				  Int_t   TriggerBit,
				  Int_t   GoodTrkSize
				  );

   JpsiCand    	 *AddJpsiCand(  
            	  const TLorentzVector&	 	Mu1   		    ,
            	  const TLorentzVector&	 	Mu2   		    ,
            	  const TLorentzVector&	 	Jpsi  		    ,
 				  Double_t     				Els     	    ,
 				  Double_t    				Cos     		,
 				  Double_t    				ClS     		,
 				  Double_t    				VtxJ[3] 		,
				  Double_t    				MatchMu[2]      ,
 				  Int_t       				MuCh[2] 		
				);

   ClassDef(BcTree,1)  //Event structure
};
#endif
