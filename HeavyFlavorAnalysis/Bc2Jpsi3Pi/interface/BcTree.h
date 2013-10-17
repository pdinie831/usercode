#ifndef BcTree_h
#define BcTree_h
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TRefArray.h"
#include "TRef.h"
#include <TChain.h>
#include <iostream>

//################################################
//  BcTreeCand Class
//################################################


class BcTreeCand : public TObject {

private:

 	TLorentzVector Bc3Pi_ 	 ;
 	TLorentzVector Mu1_   	 ;
 	TLorentzVector Mu2_   	 ;
 	TLorentzVector Jpsi_  	 ;
 	TLorentzVector Pi1_   	 ;
 	TLorentzVector Pi2_   	 ;
 	TLorentzVector Pi3_   	 ;
 	TLorentzVector tPi1_   	 ;
 	TLorentzVector tPi2_   	 ;
 	TLorentzVector tPi3_   	 ;
 	TLorentzVector JpsiV_  	 ;
 	Double_t     Ris_[5]  	 ;
 	Double_t     El_     	 ;
 	Double_t     Els_     	 ;
 	Double_t     Cos_     	 ;
 	Double_t     ClS_     	 ;
 	Double_t     ClJpsi_  	 ;
 	Double_t     BS_[3] 	 ;
 	Double_t     VtxS_[3] 	 ;
    Double_t     BSCovariance_[9];
  	Double_t     SVCovariance_[9];
	Double_t     IP3D_[3] 	 ;
	Double_t     IPbs_[3] 	 ;
	Double_t     sIP3D_[3]	 ;
	Double_t     sIPbs_[3]	 ;
 	Double_t	 MatchMu_[2] ;   
 	Double_t     MatchPi1_[2];   
 	Double_t     MatchPi2_[2];   
 	Double_t     MatchPi3_[2];   
 	Int_t	     PiCh_[3]    ;   
 	Int_t	     MuCh_[2]    ;   
    Double_t     ElsigJpsi_  ;
	Double_t     CosJpsi_    ;
 	Int_t 		 TrkNHits_[3]       ;
 	Int_t  		 TrkPixelHits_[3]   ;
 	Float_t		 TrkChi2_[3]        ;
  	Double_t     HltMatch_[4]	 	;

 
//      
public:
//
//
        BcTreeCand(){;}
        BcTreeCand(const BcTreeCand& orig);
        BcTreeCand( 	
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
			Double_t     Ris[5]     ,
 			Double_t     El         ,
 			Double_t     Els        ,
 			Double_t     Cos        ,
 			Double_t     ClS        ,
 			Double_t     ClJpsi     ,
 			Double_t     BS[3] 	    ,
 			Double_t     VtxS[3]    ,
 		 	Double_t     BSCovariance[9],
 		 	Double_t     SVCovariance[9],
			Double_t     IP3D[3]    ,
			Double_t     IPbs[3]    ,
			Double_t     sIP3D[3]   ,
			Double_t     sIPbs[3]   ,
			Double_t     MatchMu[2] ,
			Double_t     MatchPi1[2],
			Double_t     MatchPi2[2],
			Double_t     MatchPi3[2],
 			Int_t	     PiCh[3]    ,
 			Int_t	     MuCh[2]    ,
          	Double_t     ElsigJpsi  ,
          	Double_t     CosJpsi    ,
 		 	Int_t        TrkNHits[3]    , 
 		 	Int_t        TrkPixelHits[3],  
 		 	Float_t      TrkChi2[3] 	,     
 		 	Double_t     HltMatch[4]      
                   );

	virtual ~BcTreeCand();
 
  	const TLorentzVector   *GetBc3Pi()	     {return  &Bc3Pi_	   		;}
 	const TLorentzVector   *GetMu1()	     {return  &Mu1_  	   		;}
 	const TLorentzVector   *GetMu2()	     {return  &Mu2_  	   		;}
 	const TLorentzVector   *GetJpsi()	     {return  &Jpsi_ 	   		;}
 	const TLorentzVector   *GetPi1()	     {return  &Pi1_  	   		;}
 	const TLorentzVector   *GetPi2()	     {return  &Pi2_  	   		;}
 	const TLorentzVector   *GetPi3()	     {return  &Pi3_  	   		;}
 	const TLorentzVector   *GettPi1()	     {return  &tPi1_  	   		;}
 	const TLorentzVector   *GettPi2()	     {return  &tPi2_  	   		;}
 	const TLorentzVector   *GettPi3()	     {return  &tPi3_  	   		;}
 	const TLorentzVector   *GetJpsiV()	     {return  &JpsiV_	   		;}
	Double_t   GetRis(Int_t Id_)     		 {return  Ris_[Id_]    		;}
 	Double_t   GetEl()	             		 {return  El_	  	   		;}
 	Double_t   GetEls()	             		 {return  Els_  	   		;}
 	Double_t   GetCos()	             		 {return  Cos_  	   		;}
 	Double_t   GetClS()	             		 {return  ClS_  	   		;}
 	Double_t   GetClJpsi()	         		 {return  ClJpsi_  	   		;}
 	Double_t   GetBS(Int_t Id_)    			 {return  BS_[Id_]   		;}
 	Double_t   GetVtxS(Int_t Id_)    		 {return  VtxS_[Id_]   		;}
 	Double_t   GetBSCovariance(Int_t Id_)	 {return  BSCovariance_[Id_];}
 	Double_t   GetSVCovariance(Int_t Id_)	 {return  SVCovariance_[Id_];}
	Double_t   GetIP3D(Int_t Id_)    		 {return  IP3D_[Id_]   		;}
	Double_t   GetIPbs(Int_t Id_)    		 {return  IPbs_[Id_]   		;}
	Double_t   GetsIP3D(Int_t Id_)   		 {return  sIP3D_[Id_]  		;}
	Double_t   GetsIPbs(Int_t Id_)   		 {return  sIPbs_[Id_]  		;}
	Double_t   GetMatchMu(Int_t Id_) 		 {return  MatchMu_[Id_]		;}
	Double_t   GetMatchPi1(Int_t Id_)		 {return  MatchPi1_[Id_]	;}
	Double_t   GetMatchPi2(Int_t Id_)		 {return  MatchPi2_[Id_]	;}
	Double_t   GetMatchPi3(Int_t Id_) 		 {return  MatchPi3_[Id_]	;}
	Int_t	   GetPiCh(Int_t Id_)    		 {return  PiCh_[Id_]   		;}	   
	Int_t	   GetMuCh(Int_t Id_)    		 {return  MuCh_[Id_]   		;}	   
	Double_t   GetElsigJpsi()		 		 {return  ElsigJpsi_   		;}
	Double_t   GetCosJpsi() 		 		 {return  CosJpsi_     		;}
	Int_t	   GetTrkNHits(Int_t Id_)		 {return  TrkNHits_[Id_]	;}  	
	Int_t	   GetTrkPixelHIts(Int_t Id_)	 {return  TrkPixelHits_[Id_];}  	
	Float_t	   GetTrkChi2(Int_t Id_)		 {return  TrkChi2_[Id_]		;}  	
	Double_t   GetHltMatch(Int_t Id_)		 {return  HltMatch_[Id_] 	;}

        ClassDef(BcTreeCand,1) 
};


//##################################################################################################################
//  BcTreeGENCand Class
//##################################################################################################################


class BcTreeGENCand : public TObject {

private:

	TLorentzVector  Bc3Pi_  ;
 	TLorentzVector  Mu1_    ;
 	TLorentzVector  Mu2_    ;
 	TLorentzVector  JPsi_   ;
 	TLorentzVector  Pi1_    ;
 	TLorentzVector  Pi2_    ;
 	TLorentzVector  Pi3_    ;
 	Double_t     	El_     ;
 	Double_t     	Cos_    ;
 	Double_t     	BS_[3]  ;
 	Double_t     	VtxS_[3];
	Double_t     	Ris_[5] ;
 	Int_t	     	PiCh_[3];   
 	Int_t	     	MuCh_[2];   
 
//      
public:
//
//
    BcTreeGENCand(){;}
    virtual ~BcTreeGENCand(){ ;}


  	void   SetBc3Pi           (TLorentzVector n)  	{   Bc3Pi_  	= n ;}
    void   SetMu1	      	  (TLorentzVector n)  	{   Mu1_	  	= n ;}
	void   SetMu2	      	  (TLorentzVector n)  	{   Mu2_	  	= n ;}
	void   SetJPsi	          (TLorentzVector n)  	{   JPsi_   	= n ;}
	void   SetPi1	      	  (TLorentzVector n)  	{   Pi1_	  	= n ;}
	void   SetPi2	      	  (TLorentzVector n)  	{   Pi2_	  	= n ;}
	void   SetPi3	      	  (TLorentzVector n)  	{   Pi3_	  	= n ;}
 	void   SetEl	      	  (Double_t n)  		{   El_	      	= n ;}
 	void   SetCos	      	  (Double_t n)  		{   Cos_	 	= n ;}
 	void   SetBSx             (Double_t n)  		{   BS_[0]  	= n ;}
 	void   SetBSy             (Double_t n)  		{   BS_[1]  	= n ;}
 	void   SetBSz             (Double_t n)  		{   BS_[2]  	= n ;}
 	void   SetVtxSx           (Double_t n)  		{   VtxS_[0]  	= n ;}
 	void   SetVtxSy           (Double_t n)  		{   VtxS_[1]  	= n ;}
 	void   SetVtxSz           (Double_t n)  		{   VtxS_[2]  	= n ;}
	void   SetRis0			  (Double_t n)  		{   Ris_[0]   	= n ;}
	void   SetRis1  		  (Double_t n)  		{   Ris_[1]   	= n ;}
	void   SetRis2  		  (Double_t n)  		{   Ris_[2]   	= n ;}
	void   SetRis3  		  (Double_t n)  		{   Ris_[3]   	= n ;}
	void   SetRis4  		  (Double_t n)  		{   Ris_[4]   	= n ;}
	void   SetPiCh1    		  (Int_t n)	    		{   PiCh_[0]    = n ;}	
	void   SetPiCh2    	  	  (Int_t n)	    		{   PiCh_[1]    = n ;}	
	void   SetPiCh3     	  (Int_t n)	    		{   PiCh_[2]    = n ;}	
	void   SetMuCh1    		  (Int_t n)	    		{   MuCh_[0]    = n ;}	
	void   SetMuCh2    	  	  (Int_t n)	    		{   MuCh_[1]    = n ;}	

 	TLorentzVector   GetBc3Pi()    {return  Bc3Pi_     ;}
    TLorentzVector   GetMu1()	   {return  Mu1_	   ;}
	TLorentzVector   GetMu2()	   {return  Mu2_	   ;}
	TLorentzVector   GetJPsi()	   {return  JPsi_      ;}
	TLorentzVector   GetPi1()	   {return  Pi1_	   ;}
	TLorentzVector   GetPi2()	   {return  Pi2_	   ;}
	TLorentzVector   GetPi3()	   {return  Pi3_	   ;}
 	Double_t   GetEl()	           {return  El_	       ;}
 	Double_t   GetCos()	           {return  Cos_	   ;}
 	Double_t   GetBS(Int_t Id_)    {return  BS_[Id_] ;}
 	Double_t   GetVtxS(Int_t Id_)  {return  VtxS_[Id_] ;}
	Double_t   GetRis(Int_t Id_)   {return  Ris_[Id_]  ;}
	Int_t	   GetPiCh(Int_t Id_)  {return  PiCh_[Id_] ;}	   
	Int_t	   GetMuCh(Int_t Id_)  {return  MuCh_[Id_] ;}	   

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
 
//      
public:
//
//
    JpsiCand(){;}
    JpsiCand(const JpsiCand& orig);
    JpsiCand( 	
            const TLorentzVector&   Mu1   		,
            const TLorentzVector&   Mu2   		,
            const TLorentzVector&   Jpsi  		,
 			Double_t     			Els         ,
 			Double_t     			Cos         ,
 			Double_t     			ClS         ,
 			Double_t     			VtxJ[3]     ,
			Double_t     			MatchMu[2]  ,
 			Int_t	     			MuCh[2]        
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
   UInt_t       NEve_  		;
   UInt_t       Run_   		;
   UInt_t       LBlk_  		;
   UInt_t       Event_ 		;
   UInt_t       NumBc_ 		;
   Int_t        Ntrk_  		;
   Int_t        Nprim_ 		;
   Int_t        TriggerBit_ ;
   Int_t        GoodTrkSize_;
   Float_t      TrueNI_		;

public:
   BcTreeHeader() { ;}
   virtual ~BcTreeHeader()  { ;}
   void	   SetNEve     		(UInt_t    n){ NEve_   = n;}
   void	   SetRun      		(UInt_t    n){ Run_    = n;}
   void	   SetLBlk     		(UInt_t    n){ LBlk_   = n;}  
   void	   SetEvent    		(UInt_t    n){ Event_  = n;}
   void	   SetNumBc    		(UInt_t    n){ NumBc_  = n;}	
   void    SetNtrk     		(Int_t     n){ Ntrk_   = n;}
   void    SetNprim    		(Int_t     n){ Nprim_  = n;}
   void    SetTrueNI   		(Float_t   n){ TrueNI_ = n;}
   void    SetTriggerBit 	(Int_t     n){ TriggerBit_  = n;}
   void    SetGoodTrkSize	(Int_t     n){ GoodTrkSize_ = n;}

   UInt_t  GetNEve()   		{return  NEve_  		;}
   UInt_t  GetRun()    		{return  Run_   		;}
   UInt_t  GetLBlk()   		{return  LBlk_  		;}
   UInt_t  GetEvent()  		{return  Event_ 		;}
   UInt_t  GetNumBc()  		{return  NumBc_ 		;}
   Int_t   GetNtrk()   		{return  Ntrk_  		;}
   Int_t   GetNprim()  		{return  Nprim_ 		;}
   Float_t GetTrueNI() 		{return  TrueNI_		;}
   Int_t   GetTriggerBit()  {return  TriggerBit_ 	;}
   Int_t   GetGoodTrkSize() {return  GoodTrkSize_	;}

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
   TClonesArray  *GetBcTreeJpsiArrayCand() 	const {return fBcTreeJpsiArrayCand	;}
   TClonesArray  *GetBcTreeArrayCand() 		const {return fBcTreeArrayCand		;}
   TRefArray     *GetMaxPtCand()       		const {return fMaxPtCand      		;}
   Int_t          GetNBcTreeCand()     		const {return NBcTreeCand     		;}
   Int_t          GetNJpsiCand()     		const {return NJpsiCand       		;}
   Int_t          GetNumJPsi()         		const {return NumJPsi         		;}
   BcTreeHeader *GetBcTreeHeader()          	  {return &fBcTreeHeader  		;}
   BcTreeGENCand*GetBcTreeGENCand()         	  {return &fBcTreeGENCand 		;}
//
   BcTreeCand    *AddBcTreeCand(  
   				  const TLorentzVector& Bc3Pi 			,
            	  const TLorentzVector&	Mu1   			,
            	  const TLorentzVector&	Mu2   			,
            	  const TLorentzVector&	Jpsi  			,
            	  const TLorentzVector&	Pi1   			,
            	  const TLorentzVector&	Pi2   			,
            	  const TLorentzVector&	Pi3   			,
            	  const TLorentzVector&	tPi   			,	
            	  const TLorentzVector&	tPi2  			,
            	  const TLorentzVector&	tPi3  			,
            	  const TLorentzVector&	JpsiV			,
   				  	    Double_t     	Ris[5]  	  	,	
 				  		Double_t     	El      	  	,
 				  		Double_t     	Els     	  	,
 				  		Double_t     	Cos     	  	,
 				  		Double_t     	ClS     	  	,
 				  		Double_t     	ClJpsi  	  	,
 				  		Double_t     	BS[3]    		,
 				  		Double_t     	VtxS[3]    		,
				  		Double_t     	BSCovariance[9] ,
				  		Double_t     	SVCovariance[9] ,
				  		Double_t     	IP3D[3]    		,
				  		Double_t     	IPbs[3]    		,
				  		Double_t     	sIP3D[3]   		,
				  		Double_t     	sIPbs[3]   		,
				  		Double_t     	MatchMu[2] 		,
				  		Double_t     	MatchPi1[2]		,
				  		Double_t     	MatchPi2[2]		,
				  		Double_t     	MatchPi3[2]		,
 				  		Int_t        	PiCh[3]   		,
 				  		Int_t        	MuCh[2]    		,
 				  		Double_t     	ElsigJpsi  		,     
 				  		Double_t     	CosJpsi    		,
 				  		Int_t        	TrkNHits[3]		,  
 				  		Int_t        	TrkPixelHits[3] , 
 				  		Float_t      	TrkChi2[3]      ,
				  		Double_t     	HltMatch[4]      
 				  
 				  				);
//
   void SetBcTreeGENCand(  
 				  TLorentzVector    Bc3Pi 	,
 				  TLorentzVector    Mu1   	,
 				  TLorentzVector    Mu2   	,
 				  TLorentzVector    JPsi  	,
 				  TLorentzVector    Pi1   	,
 				  TLorentzVector    Pi2   	,
 				  TLorentzVector    Pi3   	,
 				  Double_t     		El      ,
 				  Double_t     		Cos     ,
 				  Double_t     		BS[3] 	,
 				  Double_t     		VtxS[3] ,
   				  Double_t     		Ris[5]  ,	
 				  Int_t        		PiCh[3] ,
 				  Int_t        		MuCh[2] 
				);
//
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
//
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
