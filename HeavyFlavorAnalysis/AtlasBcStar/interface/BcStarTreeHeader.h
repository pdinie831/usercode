#ifndef BcStarTreeHeader_h
#define BcStarTreeHeader_h
#include "TObject.h"
#include "TLorentzVector.h"


//##########################################################################
//  BcStarTreeHeader Class
//##########################################################################

class BcStarTreeHeader {

private:

   	UInt_t       			NEve_  				;
   	UInt_t       			Run_   				;
   	UInt_t       			LBlk_  				;
   	UInt_t      	 		Event_ 				;
   	Int_t        			Ntrk_  				;
   	Int_t        			Nprim_ 				;
   	Int_t        			NGoodTrk_	 		;
   	Float_t      			TrueNI_				;
 	Double_t     			BS_[3] 	     		; 	
    Double_t     			BSCovariance_[9]	;
 	Double_t     			HighestPtPV_[3] 	; 	
    Double_t     			HighestPtPVCov_[9]	;

   	Int_t	      		 	NumBc_ 				;//controllare se si riesce a settare alla fine


public:
    BcStarTreeHeader() { ;}
    virtual ~BcStarTreeHeader()  { ;}


// ------- Set information into ntuple ----------------------

   	void	SetNEve       		(UInt_t    n)			{ 	NEve_  					= n	;}
   	void	SetRun        		(UInt_t    n)			{ 	Run_   			 		= n	;}
   	void	SetLBlk       		(UInt_t    n)			{ 	LBlk_   				= n	;}  
   	void	SetEvent      		(UInt_t    n)			{ 	Event_  				= n	;}
   	void    SetNtrk       		(Int_t     n)			{ 	Ntrk_   				= n	;}
   	void    SetNprim      		(Int_t     n)			{ 	Nprim_  				= n	;}
   	void    SetNGoodTrk			(Int_t     n)			{ 	NGoodTrk_ 				= n	;}
   	void    SetTrueNI     		(Float_t   n)			{ 	TrueNI_ 				= n	;}
  	void   	SetBS   			(int id, Double_t	n) 	{   BS_[id]					= n ;}
  	void   	SetBSCovariance   	(int id, Double_t	n) 	{   BSCovariance_[id]		= n ;}
  	void   	SetHighestPtPV   	(int id, Double_t	n)  {   HighestPtPV_[id]		= n ;}
  	void   	SetHighestPtPVCov   (int id, Double_t	n)  {   HighestPtPVCov_[id]		= n ;}
    void	SetNumBc      		(UInt_t    n)			{ 	NumBc_  				= n	;}	



 // ------- Retrieve information from ntuple ----------------------

  	UInt_t		GetNEve()								{ 	return 	NEve_  				;}
   	UInt_t		GetRun()								{ 	return 	Run_   			 	;}
   	UInt_t		GetLBlk()								{ 	return 	LBlk_   			;}  
   	UInt_t		GetEvent()								{ 	return 	Event_  			;}
   	Int_t    	GetNtrk()								{ 	return 	Ntrk_   			;}
   	Int_t    	GetNprim()								{ 	return 	Nprim_  			;}
   	Int_t    	GetNGoodTrk()							{ 	return 	NGoodTrk_	 		;}
   	Float_t    	GetTrueNI()								{ 	return 	TrueNI_ 			;}
  	Double_t   	GetBS   			(int id) 			{ 	return 	BS_[id]				;}
  	Double_t   	GetBSCovariance   	(int id) 			{   return 	BSCovariance_[id]	;}
  	Double_t   	GetHighestPtPV   	(int id)  			{  	return  HighestPtPV_[id]	;}
  	Double_t   	GetHighestPtPVCov   (int id)  			{ 	return  HighestPtPVCov_[id]	;}
    UInt_t	    GetNumBc()								{ 	return 	NumBc_  			;}	


   ClassDef(BcStarTreeHeader,1)  //Event Header
};

#endif
