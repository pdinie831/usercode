#ifndef BcStarTreeGENCand_h
#define BcStarTreeGENCand_h
#include "TObject.h"
#include "TLorentzVector.h"

//##############################################################################
//  BcStarTreeGENCand Class
//##############################################################################


class BcStarTreeGENCand : public TObject {

private:

	TLorentzVector     BcPi_    	;
 	TLorentzVector     MuP_     	;
 	TLorentzVector     MuM_     	;
 	TLorentzVector     JPsi_    	;
 	TLorentzVector     Pi_      	;
  	Double_t     	   El_  		;
 	Double_t     	   Cos_ 		;
  	Double_t     	   Tau_  		;
 	Double_t     	   PVtx_[3] 	;
 	Double_t     	   BcVtx_[3]	;
 	Double_t     	   JpsiVtx_[3]	;
 	Double_t     	   BS_[3]   	;
  	Int_t	     	   PiCh_  		;
 

public:

    BcStarTreeGENCand(){;}
    virtual ~BcStarTreeGENCand(){ ;}

// ------- Set information into ntuple ----------------------

   	void   SetBcPi          (TLorentzVector n) 		{   BcPi_    	= n ;}
    void   SetMuP	      	(TLorentzVector n) 		{   MuP_	  	= n ;}
	void   SetMuM	      	(TLorentzVector n) 		{   MuM_	  	= n ;}
	void   SetJPsi	        (TLorentzVector n) 		{   JPsi_    	= n ;}
	void   SetPi	      	(TLorentzVector n) 		{   Pi_	  		= n ;}
 	void   SetEl	      	(Double_t n)  	 		{   El_	  		= n ;}
 	void   SetCos	      	(Double_t n)  	 		{   Cos_	  	= n ;}
 	void   SetTau	      	(Double_t n)  	 		{   Tau_	  	= n ;}
 	void   SetPiCh	      	(Int_t n)  	 			{   PiCh_	  	= n ;}

  	void   SetBcVtx 		(int id, Double_t	n)	{   BcVtx_[id]	= n ;}
  	void   SetPVtx 			(int id, Double_t	n)	{   PVtx_[id]	= n ;}
  	void   SetJpsiVtx 		(int id, Double_t	n)	{   JpsiVtx_[id]= n ;}
  	void   SetBS	 		(int id, Double_t	n)	{   BS_[id]		= n ;}


// ------- Retrieve information from ntuple ----------------------

   	const TLorentzVector  *GetBcPi() 			{   return &BcPi_     	;}
    const TLorentzVector  *GetMuP() 			{   return &MuP_	   	;}
	const TLorentzVector  *GetMuM() 			{   return &MuM_	   	;}
	const TLorentzVector  *GetJPsi() 			{   return &JPsi_     	;}
	const TLorentzVector  *GetPi() 				{   return &Pi_	  		;}
 	Double_t   	GetEl()  	 					{   return El_	  		;}
 	Double_t   	GetCos()  	 					{   return Cos_	  	 	;}
 	Double_t   	GetTau()  	 					{   return Tau_	  		;}
 	Double_t   	GetPiCh()  	 					{   return PiCh_  		;}
  	Double_t   	GetBcVtx  (int id)				{   return BcVtx_[id]	;}
  	Double_t   	GetPVtx   (int id)				{   return PVtx_[id]	;}
  	Double_t   	GetJpsiVtx(int id)				{   return JpsiVtx_[id]	;}
  	Double_t   	GetBS     (int id)				{   return BS_[id]		;}

    ClassDef(BcStarTreeGENCand,1) 
};


#endif
