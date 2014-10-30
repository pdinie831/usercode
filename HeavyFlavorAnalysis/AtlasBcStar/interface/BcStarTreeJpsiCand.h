#ifndef BcStarTreeJpsiCand_h
#define BcStarTreeJpsiCand_h
#include "TObject.h"
#include "TLorentzVector.h"


//############################################################################
//  JpsiCand Class
//############################################################################

//serve solo per la trigger turn on???

class JpsiCand : public TObject {

private:

 	TLorentzVector     MuP_     	;
 	TLorentzVector     MuM_     	;
 	TLorentzVector     Jpsi_    	;
  	Double_t     	   Els_  		;
 	Double_t     	   Cos_ 		;
 	Double_t     	   ClS_ 		;
 	Double_t     	   VtxJ_[3] 	;
 	Double_t	 	   MatchMu_[2]  ;   
 

public:

    JpsiCand(){;}
    JpsiCand( 	
            const TLorentzVector&     MuP   ,
            const TLorentzVector&     MuM   ,
            const TLorentzVector&     Jpsi  ,
 			Double_t     Els            	,
 			Double_t     Cos            	,
 			Double_t     ClS            	,
 			Double_t     VtxJ[3]        	,
			Double_t     MatchMu[2]     
            )
        {
       	  MuP_       	=  MuP           	;
 		  MuM_       	=  MuM           	;
 		  Jpsi_      	=  Jpsi          	;
 		  Els_	   		=  Els           	;
 		  Cos_	   		=  Cos           	;
 		  ClS_	   		=  ClS           	;
  		  MatchMu_[0]	=  MatchMu[0]    	;
 		  MatchMu_[1]	=  MatchMu[1]    	;
 		  for (int i=0; i < 3; i++)
 		  {
 		    VtxJ_[i]	= 	VtxJ[i]		;
 		  }  
        }

    ~JpsiCand(){};

    const TLorentzVector   *GetMuP()	    		{	return      &MuP_			;}
	const TLorentzVector   *GetMuM()	    		{	return      &MuM_			;}
	const TLorentzVector   *GetJPsi()	    		{	return      &Jpsi_  		;}
 	Double_t    	 		GetEls()				{	return  	Els_			;}
 	Double_t    	 		GetCos()				{	return  	Cos_			;}
 	Double_t    	 		GetClS()				{	return  	ClS_			;}
 	Double_t    	 		GetVtx(Int_t Id_)  		{	return  	VtxJ_[Id_] 		;}
 	Double_t    	 		GetMatchMu(Int_t Id_)  	{	return  	MatchMu_[Id_] 	;}


// ------- Copy a JpsiCand into another ----------------------

	JpsiCand(const JpsiCand &orig) : TObject(orig)
	{
 		MuP_       		=  orig.MuP_           ;
 		MuM_       		=  orig.MuM_           ;
 		Jpsi_      		=  orig.Jpsi_          ;
 		Els_	   		=  orig.Els_           ;
 		Cos_	   		=  orig.Cos_           ;
 		ClS_	   		=  orig.ClS_           ;
  		MatchMu_[0]		=  orig.MatchMu_[0]    ;
 		MatchMu_[1]		=  orig.MatchMu_[1]    ;
 		for (int i=0; i < 3; i++)
 		{
 		  VtxJ_[i]		= 	orig.VtxJ_[i]		;
 		}  
	}



    ClassDef(JpsiCand,1) 
};



#endif
