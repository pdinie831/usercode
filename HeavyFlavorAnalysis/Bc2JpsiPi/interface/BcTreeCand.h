#ifndef BcTreeCand_h
#define BcTreeCand_h
#include "TObject.h"
#include "TLorentzVector.h"


//################################################
//  BcTreeCand Class
//################################################


class BcTreeCand : public TObject {

private:

    // Bc cand
    TLorentzVector     BcPi_                   ;
    TLorentzVector     BcK_                    ;

    // Jpsi cand
    TLorentzVector     Jpsi_                   ;
    TLorentzVector     JpsiV_                  ;
    Float_t            DCA_                    ;

    // Muons
    TLorentzVector     MuP_                    ;
    TLorentzVector     MuM_                    ;

    bool               MuPGlobal_              ;
    bool               MuPTracker_             ;
    bool               MuPPFlow_               ;
    bool               MuPTMOST_               ;
    Int_t              MuPTrkLayerMeas_        ;
    Int_t              MuPPixLayerMeas_        ;
    Int_t              MuPPixHits_             ;
    Int_t              MuPTrkHits_             ;
    Int_t              MuPMatchedStations_     ;
    Float_t            MuPNormChi2_            ;
    Float_t            MuPDxy_                 ;
    Float_t            MuPDz_                  ;

    bool               MuMGlobal_              ;
    bool               MuMTracker_             ;
    bool               MuMPFlow_               ;
    bool               MuMTMOST_               ;
    Int_t              MuMTrkLayerMeas_        ;
    Int_t              MuMPixLayerMeas_        ;
    Int_t              MuMPixHits_             ;
    Int_t              MuMTrkHits_             ;
    Int_t              MuMMatchedStations_     ;
    Float_t            MuMNormChi2_            ;
    Float_t            MuMDxy_                 ;
    Float_t            MuMDz_                  ;
     
    // Track 
    TLorentzVector     Pi_                     ;
    Int_t              PiCh_                   ;   
    Int_t              TrkPixLayerMeas_        ;
    Int_t              TrkTrkLayerMeas_        ;
    Int_t              TrkPixHits_             ;
    Int_t              TrkTrkHits_             ;

    Float_t            TrkNormChi2_            ;
    Double_t           IP3DWrtJpsi_            ;
    Double_t           IP3DWrtJpsiSign_        ;
    Double_t           IP2DWrtBS_              ;
    Double_t           IP2DWrtBSSign_          ;
    Double_t           IP3DWrtLongPV_          ;
    Double_t           IP3DWrtLongPVSign_      ;
    Double_t           IP3DWrtPointPV_         ;
    Double_t           IP3DWrtPointPVSign_     ;
    Double_t           DeltaR_                 ;
    std::vector<float> TrackDxyPVs_            ;
    std::vector<float> TrackDzPVs_             ;

    // Bc vertex 
    Double_t           ClS_                    ;
    Double_t           BcVtxPosition_[3]       ;
    Double_t           BcVtxCovariance_[9]     ;
    Double_t           El2DWrtBS_              ;
    Double_t           Els2DWrtBS_             ;
    Double_t           Sigma2DWrtBS_           ;
    Double_t           Cos2DWrtBS_             ;
    Double_t           El3DWrtPointPV_         ;
    Double_t           Els3DWrtPointPV_        ;
    Double_t           Sigma3DWrtPointPV_      ;
    Double_t           Cos3DWrtPointPV_        ;
    Double_t           El3DWrtLongPV_          ;
    Double_t           Els3DWrtLongPV_         ;
    Double_t           Sigma3DWrtLongPV_       ;
    Double_t           Cos3DWrtLongPV_         ;

    // Jpsi vertex 
    Double_t           ClJpsi_                 ;
    Double_t           ElsigJpsi_              ;
    Double_t           CosJpsi_                ;
    Double_t           JpsiVtxPosition_[3]     ;

    // Pointing primary vertex 
    Double_t           PointPVPosition_[3]     ;
    Double_t           PointPVCovariance_[9]   ;
    Double_t           PointPVCl_              ;

    // Longitudinal primary vertex 
    Double_t           LongPVPosition_[3]      ;
    Double_t           LongPVCovariance_[9]    ;
    Double_t           LongPVCl_               ;

    // MC matching 
    Double_t           MatchMuP_               ;   
    Double_t           MatchMuM_               ;   
    Double_t           MatchPi_                ;   
         
    // HLT matching          
    Double_t           HltMatch_[4]            ;



public:

    BcTreeCand(){;}
//     BcTreeCand(const BcTreeCand& orig);
    ~BcTreeCand() {}
//     virtual ~BcTreeCand();

// ------- Set information into ntuple ----------------------

    // Bc cand
    void   SetBcPi                (TLorentzVector n)          {   BcPi_                   = n ;}
    void   SetBcK                 (TLorentzVector n)          {   BcK_                    = n ;}
      
    // Jpsi cand      
    void   SetJpsi                (TLorentzVector n)          {   Jpsi_                   = n ;}
    void   SetJpsiV               (TLorentzVector n)          {   JpsiV_                  = n ;}
    void   SetDCA                 (Float_t        n)          {   DCA_                    = n ;}
      
    // Muons          
    void   SetMuP                 (TLorentzVector n)          {   MuP_                    = n ;}
    void   SetMuM                 (TLorentzVector n)          {   MuM_                    = n ;}
   
    void   SetMuPisGlobal         (bool           n)          {   MuPGlobal_              = n ;}
    void   SetMuPisTracker        (bool           n)          {   MuPTracker_             = n ;}
    void   SetMuPisPFlow          (bool           n)          {   MuPPFlow_               = n ;}
    void   SetMuPTMOST            (bool           n)          {   MuPTMOST_               = n ;}
    void   SetMuPTrkLMeas         (Int_t          n)          {   MuPTrkLayerMeas_        = n ;}
    void   SetMuPPixLMeas         (Int_t          n)          {   MuPPixLayerMeas_        = n ;}
    void   SetMuPPixHits          (Int_t          n)          {   MuPPixHits_             = n ;}
    void   SetMuPTrkHits          (Int_t          n)          {   MuPTrkHits_             = n ;}
    void   SetMuPMatchedStat      (Int_t          n)          {   MuPMatchedStations_     = n ;}
    void   SetMuPNormChi2         (Float_t        n)          {   MuPNormChi2_            = n ;}
    void   SetMuPDxy              (Float_t        n)          {   MuPDxy_                 = n ;}
    void   SetMuPDz               (Float_t        n)          {   MuPDz_                  = n ;}

    void   SetMuMisGlobal         (bool           n)          {   MuMGlobal_              = n ;}
    void   SetMuMisTracker        (bool           n)          {   MuMTracker_             = n ;}
    void   SetMuMisPFlow          (bool           n)          {   MuMPFlow_               = n ;}
    void   SetMuMTMOST            (bool           n)          {   MuMTMOST_               = n ;}
    void   SetMuMTrkLMeas         (Int_t          n)          {   MuMTrkLayerMeas_        = n ;}
    void   SetMuMPixLMeas         (Int_t          n)          {   MuMPixLayerMeas_        = n ;}
    void   SetMuMPixHits          (Int_t          n)          {   MuMPixHits_             = n ;}
    void   SetMuMTrkHits          (Int_t          n)          {   MuMTrkHits_             = n ;}
    void   SetMuMMatchedStat      (Int_t          n)          {   MuMMatchedStations_     = n ;}
    void   SetMuMNormChi2         (Float_t        n)          {   MuMNormChi2_            = n ;}
    void   SetMuMDxy              (Float_t        n)          {   MuMDxy_                 = n ;}
    void   SetMuMDz               (Float_t        n)          {   MuMDz_                  = n ;}

    // Track 
    void   SetPi                  (TLorentzVector n)          {   Pi_                     = n ;}
    void   SetPiCh                (Int_t          n)          {   PiCh_                   = n ;}
    void   SetTrkPixLMeas         (Int_t          n)          {   TrkPixLayerMeas_        = n ;}
    void   SetTrkTrkLMeas         (Int_t          n)          {   TrkTrkLayerMeas_        = n ;}
    void   SetTrkPixHits          (Int_t          n)          {   TrkPixHits_             = n ;}
    void   SetTrkTrkHits          (Int_t          n)          {   TrkTrkHits_             = n ;}
    void   SetTrkNormChi2         (Float_t        n)          {   TrkNormChi2_            = n ;}
    void   SetIP3DJpsi            (Double_t       n)          {   IP3DWrtJpsi_            = n ;}
    void   SetIP3DJpsiSign        (Double_t       n)          {   IP3DWrtJpsiSign_        = n ;}
    void   SetIP2DBS              (Double_t       n)          {   IP2DWrtBS_              = n ;}
    void   SetIP2DBSSign          (Double_t       n)          {   IP2DWrtBSSign_          = n ;}
    void   SetIP3DLongPV          (Double_t       n)          {   IP3DWrtLongPV_          = n ;}
    void   SetIP3DLongPVSign      (Double_t       n)          {   IP3DWrtLongPVSign_      = n ;}
    void   SetIP3DPointPV         (Double_t       n)          {   IP3DWrtPointPV_         = n ;}
    void   SetIP3DPointPVSign     (Double_t       n)          {   IP3DWrtPointPVSign_     = n ;}
    void   SetDeltaR              (Double_t       n)          {   DeltaR_                 = n ;}
    void   SetTrackDxyPVs         (std::vector<float> n)      {   TrackDxyPVs_            = n ;}
    void   SetTrackDzPVs          (std::vector<float> n)      {   TrackDzPVs_             = n ;}


    // Bc vertex 
    void   SetClS                 (Double_t       n)          {   ClS_                    = n ;}
    void   SetEl2DWrtBS           (Double_t       n)          {   El2DWrtBS_              = n ;}
    void   SetEls2DWrtBS          (Double_t       n)          {   Els2DWrtBS_             = n ;}
    void   SetSigma2DWrtBS        (Double_t       n)          {   Sigma2DWrtBS_           = n ;}
    void   SetCos2DWrtBS          (Double_t       n)          {   Cos2DWrtBS_             = n ;}
    void   SetEl3DWrtPointPV      (Double_t       n)          {   El3DWrtPointPV_         = n ;}
    void   SetEls3DWrtPointPV     (Double_t       n)          {   Els3DWrtPointPV_        = n ;}
    void   SetSigma3DWrtPointPV   (Double_t       n)          {   Sigma3DWrtPointPV_      = n ;}
    void   SetCos3DWrtPointPV     (Double_t       n)          {   Cos3DWrtPointPV_        = n ;}
    void   SetEl3DWrtLongPV       (Double_t       n)          {   El3DWrtLongPV_          = n ;}
    void   SetEls3DWrtLongPV      (Double_t       n)          {   Els3DWrtLongPV_         = n ;}
    void   SetSigma3DWrtLongPV    (Double_t       n)          {   Sigma3DWrtLongPV_       = n ;}
    void   SetCos3DWrtLongPV      (Double_t       n)          {   Cos3DWrtLongPV_         = n ;}

    void   SetBcVtxPosition       (int id, Double_t   n)      {   BcVtxPosition_[id]      = n ;}
    void   SetBcVtxCovariance     (int id, Double_t   n)      {   BcVtxCovariance_[id]    = n ;}

    // Jpsi vertex 
    void   SetClJpsi              (Double_t       n)          {   ClJpsi_                 = n ;}
    void   SetElsigJpsi           (Double_t       n)          {   ElsigJpsi_              = n ;}
    void   SetCosJpsi             (Double_t       n)          {   CosJpsi_                = n ;}
    
    void   SetJpsiVtxPosition     (int id, Double_t   n)      {   JpsiVtxPosition_[id]    = n ;}

    // Pointing primary vertex 
    void   SetPointPVPosition     (int id, Double_t   n)      {   PointPVPosition_[id]    = n ;}
    void   SetPointPVCovariance   (int id, Double_t   n)      {   PointPVCovariance_[id]  = n ;}
    void   SetPointPVCl           (Double_t           n)      {   PointPVCl_              = n ;}

    // Longitudinal primary vertex 
    void   SetLongPVPosition      (int id, Double_t   n)      {   LongPVPosition_[id]     = n ;}
    void   SetLongPVCovariance    (int id, Double_t   n)      {   LongPVCovariance_[id]   = n ;}
    void   SetLongPVCl            (Double_t           n)      {   LongPVCl_               = n ;}

    // MC matching 
    void   SetMatchMuP            (Double_t       n)          {   MatchMuP_               = n ;}
    void   SetMatchMuM            (Double_t       n)          {   MatchMuM_               = n ;}
    void   SetMatchPi             (Double_t       n)          {   MatchPi_                = n ;}

    // HLT matching 
    void   SetHltMatch            (int id, Double_t   n)      {   HltMatch_[id]           = n ;}




// ------- Retrieve information from ntuple ----------------------

    // Bc cand
    const TLorentzVector   *GetBcPi()                    {   return  &BcPi_               ;}
    const TLorentzVector   *GetBcK()                     {   return  &BcK_                ;}

    // Jpsi cand
    const TLorentzVector   *GetJpsi()                    {   return  &Jpsi_               ;}
    const TLorentzVector   *GetJpsiV()                   {   return  &JpsiV_              ;}
    Float_t                 GetDCA()                     {   return  DCA_                 ;}

    // Muons
    const TLorentzVector   *GetMuP()                     {   return  &MuP_                ;}
    const TLorentzVector   *GetMuM()                     {   return  &MuM_                ;}

    bool                    GetMuPisGlobal()             {   return  MuPGlobal_           ;}
    bool                    GetMuPisTracker()            {   return  MuPTracker_          ;}
    bool                    GetMuPisPFlow()              {   return  MuPPFlow_            ;}
    bool                    GetMuPTMOST()                {   return  MuPTMOST_            ;}
    Int_t                   GetMuPTrkLMeas()             {   return  MuPTrkLayerMeas_     ;}
    Int_t                   GetMuPPixLMeas()             {   return  MuPPixLayerMeas_     ;}
    Int_t                   GetMuPPixHits()              {   return  MuPPixHits_          ;}
    Int_t                   GetMuPTrkHits()              {   return  MuPTrkHits_          ;}
    Int_t                   GetMuPMatchedStat()          {   return  MuPMatchedStations_  ;}
    Float_t                 GetMuPNormChi2()             {   return  MuPNormChi2_         ;}
    Float_t                 GetMuPDxy()                  {   return  MuPDxy_              ;}
    Float_t                 GetMuPDz()                   {   return  MuPDz_               ;}
    
    bool                    GetMuMisGlobal()             {   return  MuMGlobal_           ;}
    bool                    GetMuMisTracker()            {   return  MuMTracker_          ;}
    bool                    GetMuMisPFlow()              {   return  MuMPFlow_            ;}
    bool                    GetMuMTMOST()                {   return  MuMTMOST_            ;}
    Int_t                   GetMuMTrkLMeas()             {   return  MuMTrkLayerMeas_     ;}
    Int_t                   GetMuMPixLMeas()             {   return  MuMPixLayerMeas_     ;}
    Int_t                   GetMuMPixHits()              {   return  MuMPixHits_          ;}
    Int_t                   GetMuMTrkHits()              {   return  MuMTrkHits_          ;}
    Int_t                   GetMuMMatchedStat()          {   return  MuMMatchedStations_  ;}
    Float_t                 GetMuMNormChi2()             {   return  MuMNormChi2_         ;}
    Float_t                 GetMuMDxy()                  {   return  MuMDxy_              ;}
    Float_t                 GetMuMDz()                   {   return  MuMDz_               ;}
    
    // Track 
    const TLorentzVector   *GetPi()                      {   return  &Pi_                 ;}
    Int_t                   GetPiCh()                    {   return  PiCh_                ;}
    Int_t                   GetTrkPixLMeas()             {   return  TrkPixLayerMeas_     ;}
    Int_t                   GetTrkTrkLMeas()             {   return  TrkTrkLayerMeas_     ;}
    Int_t                   GetTrkPixHits()              {   return  TrkPixHits_          ;}
    Int_t                   GetTrkTrkHits()              {   return  TrkTrkHits_          ;}
    Float_t                 GetTrkNormChi2()             {   return  TrkNormChi2_         ;}
    Double_t                GetIP3DJpsi()                {   return  IP3DWrtJpsi_         ;}
    Double_t                GetIP3DJpsiSign()            {   return  IP3DWrtJpsiSign_     ;}
    Double_t                GetIP2DBS()                  {   return  IP2DWrtBS_           ;}
    Double_t                GetIP2DBSSign()              {   return  IP2DWrtBSSign_       ;}
    Double_t                GetIP3DLongPV()              {   return  IP3DWrtLongPV_       ;}
    Double_t                GetIP3DLongPVSign()          {   return  IP3DWrtLongPVSign_   ;}
    Double_t                GetIP3DPointPV()             {   return  IP3DWrtPointPV_      ;}
    Double_t                GetIP3DPointPVSign()         {   return  IP3DWrtPointPVSign_  ;}
    Double_t                GetDeltaR()                  {   return  DeltaR_              ;}
    std::vector<float>      GetTrackDxyPVs()             {   return  TrackDxyPVs_         ;}
    std::vector<float>      GetTrackDzPVs()              {   return  TrackDzPVs_          ;}

    // Bc vertex 
    Double_t                GetClS()                     {   return  ClS_                 ;}
    Double_t                GetEl2DWrtBS()               {   return  El2DWrtBS_           ;}
    Double_t                GetEls2DWrtBS()              {   return  Els2DWrtBS_          ;}
    Double_t                GetSigma2DWrtBS()            {   return  Sigma2DWrtBS_        ;}
    Double_t                GetCos2DWrtBS()              {   return  Cos2DWrtBS_          ;}
    Double_t                GetEl3DWrtPointPV()          {   return  El3DWrtPointPV_      ;}
    Double_t                GetEls3DWrtPointPV()         {   return  Els3DWrtPointPV_     ;}
    Double_t                GetSigma3DWrtPointPV()       {   return  Sigma3DWrtPointPV_   ;}
    Double_t                GetCos3DWrtPointPV()         {   return  Cos3DWrtPointPV_     ;}
    Double_t                GetEl3DWrtLongPV()           {   return  El3DWrtLongPV_       ;}
    Double_t                GetEls3DWrtLongPV()          {   return  Els3DWrtLongPV_      ;}
    Double_t                GetSigma3DWrtLongPV()        {   return  Sigma3DWrtLongPV_    ;}
    Double_t                GetCos3DWrtLongPV()          {   return  Cos3DWrtLongPV_      ;}
    Double_t                GetBcVtxPosition(int i)      {   return  BcVtxPosition_[i]    ;}
    Double_t                GetBcVtxCovariance(int i)    {   return  BcVtxCovariance_[i]  ;}

    // Jpsi vertex 
    Double_t                GetClJpsi()                  {   return  ClJpsi_              ;}
    Double_t                GetElsigJpsi()               {   return  ElsigJpsi_           ;}
    Double_t                GetCosJpsi()                 {   return  CosJpsi_             ;}
    Double_t                GetJpsiVtxPosition(int i)    {   return  JpsiVtxPosition_[i]  ;}

    // Pointing primary vertex 
    Double_t                GetPointPVPosition(int i)    {   return  PointPVPosition_[i]  ;}
    Double_t                GetPointPVCovariance(int i)  {   return  PointPVCovariance_[i];}
    Double_t                GetPointPVCl()               {   return  PointPVCl_           ;}

    // Longitudinal primary vertex 
    Double_t                GetLongPVPosition(int i)     {   return  LongPVPosition_[i]   ;}
    Double_t                GetLongPVCovariance(int i)   {   return  LongPVCovariance_[i] ;}
    Double_t                GetLongPVCl()                {   return  LongPVCl_            ;}

    // MC matching 
    Double_t                GetMatchMuP()                {   return  MatchMuP_            ;}
    Double_t                GetMatchMuM()                {   return  MatchMuM_            ;}
    Double_t                GetMatchPi()                 {   return  MatchPi_             ;}

    // HLT matching 
    Double_t                GetHltMatch(int i)           {    return  HltMatch_[i]        ;}



// ------- Copy a BcCand into another ----------------------

BcTreeCand(const BcTreeCand &orig) : TObject(orig)
{
    // Bc cand
    BcPi_                   =  orig.BcPi_                   ;
    BcK_                    =  orig.BcK_                    ;

    // Jpsi cand
    Jpsi_                   =  orig.Jpsi_                   ;
    JpsiV_                  =  orig.JpsiV_                  ;
    DCA_                    =  orig.DCA_                    ;

    // Muons
    MuP_                    =  orig.MuP_                    ;
    MuM_                    =  orig.MuM_                    ;

    MuPGlobal_              =  orig.MuPGlobal_              ;
    MuPTracker_             =  orig.MuPTracker_             ;
    MuPPFlow_               =  orig.MuPPFlow_               ;
    MuPTMOST_               =  orig.MuPTMOST_               ;
    MuPTrkLayerMeas_        =  orig.MuPTrkLayerMeas_        ;
    MuPPixLayerMeas_        =  orig.MuPPixLayerMeas_        ;
    MuPPixHits_             =  orig.MuPPixHits_             ;
    MuPTrkHits_             =  orig.MuPTrkHits_             ;
    MuPMatchedStations_     =  orig.MuPMatchedStations_     ;
    MuPNormChi2_            =  orig.MuPNormChi2_            ;
    MuPDxy_                 =  orig.MuPDxy_                 ;
    MuPDz_                  =  orig.MuPDz_                  ;

    MuMGlobal_              =  orig.MuMGlobal_              ;
    MuMTracker_             =  orig.MuMTracker_             ;
    MuMPFlow_               =  orig.MuMPFlow_               ;
    MuMTMOST_               =  orig.MuMTMOST_               ;
    MuMTrkLayerMeas_        =  orig.MuMTrkLayerMeas_        ;
    MuMPixLayerMeas_        =  orig.MuMPixLayerMeas_        ;
    MuMPixHits_             =  orig.MuMPixHits_             ;
    MuMTrkHits_             =  orig.MuMTrkHits_             ;
    MuMMatchedStations_     =  orig.MuMMatchedStations_     ;
    MuMNormChi2_            =  orig.MuMNormChi2_            ;
    MuMDxy_                 =  orig.MuMDxy_                 ;
    MuMDz_                  =  orig.MuMDz_                  ;

    // Track  ;
    Pi_                     =  orig.Pi_                     ;
    PiCh_                   =  orig.PiCh_                   ;
    TrkPixLayerMeas_        =  orig.TrkPixLayerMeas_        ;
    TrkTrkLayerMeas_        =  orig.TrkTrkLayerMeas_        ;
    TrkPixHits_             =  orig.TrkPixHits_             ;
    TrkTrkHits_             =  orig.TrkTrkHits_             ;
    TrkNormChi2_            =  orig.TrkNormChi2_            ;
    IP3DWrtJpsi_            =  orig.IP3DWrtJpsi_            ;
    IP3DWrtJpsiSign_        =  orig.IP3DWrtJpsiSign_        ;
    IP2DWrtBS_              =  orig.IP2DWrtBS_              ;
    IP2DWrtBSSign_          =  orig.IP2DWrtBSSign_          ;
    IP3DWrtLongPV_          =  orig.IP3DWrtLongPV_          ;
    IP3DWrtLongPVSign_      =  orig.IP3DWrtLongPVSign_      ;
    IP3DWrtPointPV_         =  orig.IP3DWrtPointPV_         ;
    IP3DWrtPointPVSign_     =  orig.IP3DWrtPointPVSign_     ;
    DeltaR_                 =  orig.DeltaR_                 ;
    TrackDxyPVs_            =  orig.TrackDxyPVs_            ;
    TrackDzPVs_             =  orig.TrackDzPVs_             ;
    
    // Bc vertex  ;
    ClS_                    =  orig.ClS_                    ;
    El2DWrtBS_              =  orig.El2DWrtBS_              ;
    Els2DWrtBS_             =  orig.Els2DWrtBS_             ;
    Sigma2DWrtBS_           =  orig.Sigma2DWrtBS_           ;
    Cos2DWrtBS_             =  orig.Cos2DWrtBS_             ;
    El3DWrtPointPV_         =  orig.El3DWrtPointPV_         ;
    Els3DWrtPointPV_        =  orig.Els3DWrtPointPV_        ;
    Sigma3DWrtPointPV_      =  orig.Sigma3DWrtPointPV_      ;
    Cos3DWrtPointPV_        =  orig.Cos3DWrtPointPV_        ;
    El3DWrtLongPV_          =  orig.El3DWrtLongPV_          ;
    Els3DWrtLongPV_         =  orig.Els3DWrtLongPV_         ;
    Sigma3DWrtLongPV_       =  orig.Sigma3DWrtLongPV_       ;
    Cos3DWrtLongPV_         =  orig.Cos3DWrtLongPV_         ;

    // Jpsi vertex 
    ClJpsi_                 =  orig.ClS_                    ;
    ElsigJpsi_              =  orig.ElsigJpsi_              ;
    CosJpsi_                =  orig.CosJpsi_                ;

    //Primary vertices
    PointPVCl_              =  orig.PointPVCl_              ;
    LongPVCl_               =  orig.LongPVCl_               ;

    // MC matching 
    MatchMuP_               =  orig.MatchMuP_               ;
    MatchMuM_               =  orig.MatchMuM_               ;
    MatchPi_                =  orig.MatchPi_                ;

    for(int j=0; j<3; j++){
      BcVtxPosition_[j]     =  orig.BcVtxPosition_[j]       ;
      JpsiVtxPosition_[j]   =  orig.JpsiVtxPosition_[j]     ;
      PointPVPosition_[j]   =  orig.PointPVPosition_[j]     ;
      LongPVPosition_[j]    =  orig.LongPVPosition_[j]      ;
    }
    for(int j=0; j<9; j++){
      BcVtxCovariance_[j]   =  orig.BcVtxCovariance_[j]     ;
      PointPVCovariance_[j] =  orig.PointPVCovariance_[j]   ;
      LongPVCovariance_[j]  =  orig.LongPVCovariance_[j]    ;
    }

    // HLT matching 
    for(int j=0; j<4; j++){
      HltMatch_[j]          =  orig.HltMatch_[j]            ;
    }

}


    ClassDef(BcTreeCand,1) 
};
#endif
