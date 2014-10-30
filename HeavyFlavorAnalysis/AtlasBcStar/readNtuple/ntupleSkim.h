// Include file for Skim ntuples to plain tree
// 
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include "Riostream.h"
#include <map>
#include <string>
#include <vector>
#include <math.h>
#include <TCint.h>
#include <TMatrix.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TApplication.h>
#include <TFile.h>
#include <TStyle.h> 

#include "TBranch.h"
#include "TLorentzVector.h"
#include "interface/BcStarTree.h"

//constants
double JpsiMass = 3.096916;

// variables
char * InputFile    = new char [200];
char * OutputFile   = new char [200];
char * TreeFileName = new char [200];

std::stringstream ss_, ssW ;

BcStarTree 			*newBcTree_;
BcStarTreeHeader 	*fBcTreeHeader;

std::map<std::string, TTree* 	> t1_ ;


// Ntupla header 
UInt_t              NEve, Run, LBlk, Event              ;   
Int_t               Ntrk, Nprim, NGoodTrk               ;
Float_t             TrueNI                              ;
Double_t            BS[3], BSCovariance[9]              ;
Double_t            HighestPtPV[3] , HighestPtPVCov[9]  ;

// Ntupla Cand
TLorentzVector      BcPi, BcK, MuP, MuM, Pi, Jpsi, JpsiV;
TLorentzVector      BcStar, BcPrefit, BcPostfit, Pi1, Pi2;


// Muons
Float_t             DCA ;
bool                MuPGlobal, MuPTracker, MuPPFlow, MuPTMOST, MuMGlobal, MuMTracker, MuMPFlow, MuMTMOST;
Int_t               MuPTrkLayerMeas, MuPPixLayerMeas, MuPPixHits, MuPTrkHits, MuPMatchedStations;
Int_t               MuMTrkLayerMeas, MuMPixLayerMeas, MuMPixHits, MuMTrkHits, MuMMatchedStations;
Float_t             MuPNormChi2, MuPDxy, MuPDz, MuMNormChi2, MuMDxy, MuMDz;

// Track 
Int_t               PiCh, TrkPixLayerMeas, TrkTrkLayerMeas, TrkPixHits, TrkTrkHits;   
Float_t             TrkNormChi2;
Double_t            IP3DWrtJpsi, IP3DWrtJpsiSign, IP2DWrtBS, IP2DWrtBSSign, IP3DWrtPV, IP3DWrtPVSign, DeltaR;

// Bc vertex 
Double_t            ClS, BcVtxPosition[3], BcVtxCovariance[9];
Double_t            El2DWrtBS,      Els2DWrtBS,      Sigma2DWrtBS,      Cos2DWrtBS;
Double_t            El3DWrtPointPV, Els3DWrtPointPV, Sigma3DWrtPointPV, Cos3DWrtPointPV;
Double_t            El3DWrtLongPV,  Els3DWrtLongPV,  Sigma3DWrtLongPV,  Cos3DWrtLongPV;

// Jpsi vertex 
Double_t            ClJpsi, ElsigJpsi, CosJpsi, JpsiVtxPosition[3];

// Primary vertices 
Double_t            PointPVPosition[3], PointPVCovariance[9], PointPVCl;
Double_t            LongPVPosition[3],  LongPVCovariance[9],  LongPVCl;

// MC matching 
Double_t            MatchMuP, MatchMuM, MatchPi;   

// HLT matching 
Double_t            HltMatch[4];

//BcStar Cand
Int_t               Pi1Ch, Pi1PixHits, Pi1TrkHits, Pi2Ch, Pi2PixHits, Pi2TrkHits;
Float_t             Pi1NormChi2, Pi2NormChi2;
Double_t            Pi1DeltaR, Pi2DeltaR;
Double_t            ClBcStar, Bc_BcStarDistance, Bc_BcStarSignificance, BcStarVtxPosition[3], BcStarVtxCovariance[9], Cosine;


//Variables for plain tree
Double_t    BcPiM                  = -9999  ; 
Double_t    BcKM                   = -9999  ; 
Double_t    BcPiPt                 = -9999  ; 
Double_t    BcPiEta                = -9999  ; 
Double_t    BcPiPhi                = -9999  ; 
Double_t    BcY                    = -9999  ;  
Double_t    BpY                    = -9999  ;  

Double_t    PiPt                   = -9999  ; 
Double_t    PiEta                  = -9999  ; 
Double_t    PiPhi                  = -9999  ; 

Double_t    MuPPt                  = -9999  ;  
Double_t    MuMPt                  = -9999  ;  
Double_t    MuPEta                 = -9999  ;  
Double_t    MuMEta                 = -9999  ;  

Double_t    JpsiM                  = -9999  ;  
Double_t    JpsiEta                = -9999  ;  
Double_t    JpsiPhi                = -9999  ;  
Double_t    JpsiPt                 = -9999  ;  

Double_t    tau                    = -9999  ;

Double_t    BcStarM                = -9999  ;      
Double_t    BcStarPt               = -9999  ;      
Double_t    BcStarY                = -9999  ;      
Double_t    Pi1Pt                  = -9999  ; 
Double_t    Pi1Eta                 = -9999  ; 
Double_t    Pi2Pt                  = -9999  ; 
Double_t    Pi2Eta                 = -9999  ; 


int CountBc                       = 0      ;
int CountBcStar                   = 0      ;
int CountBcDR                     = 0      ;
bool isNewJpsi                    = false  ;
Double_t MaxBcPtArbitration       = -9999  ;                                         
Int_t    fMaxBcPtArbitration      = -9999  ;                                         

