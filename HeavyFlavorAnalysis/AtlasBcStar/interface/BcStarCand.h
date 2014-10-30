#ifndef BcStarCand_h
#define BcStarCand_h
#include "TObject.h"
#include "TLorentzVector.h"


//################################################
//  BcStarCand Class
//################################################


class BcStarCand : public TObject {

private:

    // Bc cand
    TLorentzVector     BcStar_                 ;
    TLorentzVector     Bc_                     ;
    TLorentzVector     BcPrefit_               ;
    TLorentzVector     Pi1_                    ;
    TLorentzVector     Pi2_                    ;

    Int_t              Pi1Ch_                  ;   
    Int_t              Pi1PixHits_             ;
    Int_t              Pi1TrkHits_             ;
    Float_t            Pi1NormChi2_            ;
    Double_t           Pi1DeltaR_              ;

    Int_t              Pi2Ch_                  ;   
    Int_t              Pi2PixHits_             ;
    Int_t              Pi2TrkHits_             ;
    Float_t            Pi2NormChi2_            ;
    Double_t           Pi2DeltaR_              ;

    // Bc star vertex 
    Double_t           Cl_                     ;
    Double_t           Bc_BcStarDistance_      ;
    Double_t           Bc_BcStarSignificance_  ;
    Double_t           Chi2BcVertices_         ;

    Double_t           BcVtxPosition_[3]       ;
    Double_t           BcVtxCovariance_[9]     ;
    Double_t           Cosine_                 ;



public:

    BcStarCand(){;}
    ~BcStarCand() {}

// ------- Set information into ntuple ----------------------
    // Bc cand
    void   SetBcStar                 (TLorentzVector n)          {   BcStar_               = n ;}
    void   SetBc                     (TLorentzVector n)          {   Bc_                   = n ;}
    void   SetBcPrefit               (TLorentzVector n)          {   BcPrefit_             = n ;}
    void   SetPi1                    (TLorentzVector n)          {   Pi1_                  = n ;}
    void   SetPi2                    (TLorentzVector n)          {   Pi2_                  = n ;}

    // Track 
    void   SetPi1Ch                  (Int_t          n)          {   Pi1Ch_                   = n ;}
    void   SetPi1PixHits             (Int_t          n)          {   Pi1PixHits_              = n ;}
    void   SetPi1TrkHits             (Int_t          n)          {   Pi1TrkHits_              = n ;}
    void   SetPi1NormChi2            (Float_t        n)          {   Pi1NormChi2_             = n ;}
    void   SetPi1DeltaR              (Double_t       n)          {   Pi1DeltaR_               = n ;}

    void   SetPi2Ch                  (Int_t          n)          {   Pi2Ch_                   = n ;}
    void   SetPi2PixHits             (Int_t          n)          {   Pi2PixHits_              = n ;}
    void   SetPi2TrkHits             (Int_t          n)          {   Pi2TrkHits_              = n ;}
    void   SetPi2NormChi2            (Float_t        n)          {   Pi2NormChi2_             = n ;}
    void   SetPi2DeltaR              (Double_t       n)          {   Pi2DeltaR_               = n ;}

    // Bc vertex 
    void   SetCl                     (Double_t       n)          {   Cl_                      = n ;}
    void   SetBcBcStarDistance       (Double_t       n)          {   Bc_BcStarDistance_       = n ;}
    void   SetBcBcStarSignificance   (Double_t       n)          {   Bc_BcStarSignificance_   = n ;}
    void   SetChi2BcVertices         (Double_t       n)          {   Chi2BcVertices_          = n ;}
    void   SetCosine                 (Double_t       n)          {   Cosine_                  = n ;}

    void   SetBcVtxPosition       (int id, Double_t   n)         {   BcVtxPosition_[id]       = n ;}
    void   SetBcVtxCovariance     (int id, Double_t   n)         {   BcVtxCovariance_[id]     = n ;}

// ------- Retrieve information from ntuple ----------------------

    // Bc cand
    const TLorentzVector   *GetBcStar()                     {   return  &BcStar_                ;}
    const TLorentzVector   *GetBc()                         {   return  &Bc_                    ;}
    const TLorentzVector   *GetBcPrefit()                   {   return  &BcPrefit_              ;}
    const TLorentzVector   *GetPi1()                        {   return  &Pi1_                   ;}
    const TLorentzVector   *GetPi2()                        {   return  &Pi2_                   ;}

    // Track  
    Int_t                   GetPi1Ch()                      {   return  Pi1Ch_                  ;}
    Int_t                   GetPi1PixHits()                 {   return  Pi1PixHits_             ;}
    Int_t                   GetPi1TrkHits()                 {   return  Pi1TrkHits_             ;}
    Float_t                 GetPi1NormChi2()                {   return  Pi1NormChi2_            ;}
    Double_t                GetPi1DeltaR()                  {   return  Pi1DeltaR_              ;}
 
    Int_t                   GetPi2Ch()                      {   return  Pi2Ch_                  ;}
    Int_t                   GetPi2PixHits()                 {   return  Pi2PixHits_             ;}
    Int_t                   GetPi2TrkHits()                 {   return  Pi2TrkHits_             ;}
    Float_t                 GetPi2NormChi2()                {   return  Pi2NormChi2_            ;}
    Double_t                GetPi2DeltaR()                  {   return  Pi2DeltaR_              ;}
 
    // Bc vertex 
    Double_t                GetCl()                         {   return  Cl_                     ;}
    Double_t                GetBcBcStarDistance()           {   return  Bc_BcStarDistance_      ;}
    Double_t                GetBcBcStarSignificance()       {   return  Bc_BcStarSignificance_  ;}
    Double_t                GetChi2BcVertices()             {   return  Chi2BcVertices_         ;}
    Double_t                GetCosine()                     {   return  Cosine_                 ;}
    Double_t                GetBcVtxPosition(int i)         {   return  BcVtxPosition_[i]       ;}
    Double_t                GetBcVtxCovariance(int i)       {   return  BcVtxCovariance_[i]     ;}

// ------- Copy a BcCand into another ----------------------

BcStarCand(const BcStarCand &orig) : TObject(orig)
{
    // Bc cand 
    BcStar_                 =  orig.BcStar_                 ;
    Bc_                     =  orig.Bc_                     ;
    BcPrefit_               =  orig.BcPrefit_               ;
    Pi1_                    =  orig.Pi1_                    ;
    Pi2_                    =  orig.Pi2_                    ;


    // Track  ;
    Pi1Ch_                  =  orig.Pi1Ch_                  ;
    Pi1PixHits_             =  orig.Pi1PixHits_             ;
    Pi1TrkHits_             =  orig.Pi1TrkHits_             ;
    Pi1NormChi2_            =  orig.Pi1NormChi2_            ;
    Pi1DeltaR_              =  orig.Pi1DeltaR_              ;

    Pi2Ch_                  =  orig.Pi2Ch_                  ;
    Pi2PixHits_             =  orig.Pi2PixHits_             ;
    Pi2TrkHits_             =  orig.Pi2TrkHits_             ;
    Pi2NormChi2_            =  orig.Pi2NormChi2_            ;
    Pi2DeltaR_              =  orig.Pi2DeltaR_              ;
    
    // Bc vertex  ;
    Cl_                     =  orig.Cl_                     ;
    Bc_BcStarDistance_      =  orig.Bc_BcStarDistance_      ;
    Bc_BcStarSignificance_  =  orig.Bc_BcStarSignificance_  ;
    Chi2BcVertices_         =  orig.Chi2BcVertices_         ;
    Cosine_                 =  orig.Cosine_                 ;

    for(int j=0; j<3; j++){
      BcVtxPosition_[j]     =  orig.BcVtxPosition_[j]       ;
    }
    for(int j=0; j<9; j++){
      BcVtxCovariance_[j]   =  orig.BcVtxCovariance_[j]     ;
    }

}

    ClassDef(BcStarCand,1) 
};
#endif
