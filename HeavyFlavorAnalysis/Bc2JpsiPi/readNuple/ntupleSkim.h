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
#include "interface/BcTree.h"

//constants
double JpsiMass = 3.096916;

// variables
char * InputFile    = new char [200];
char * OutputFile   = new char [200];
char * TreeFileName = new char [200];

std::stringstream ss_, ssW ;

BcTree              *newBcTree_;
BcTreeHeader        *fBcTreeHeader;

std::map<std::string, TTree*> t1_ ;

// Ntupla header 
UInt_t              NEve, Run, LBlk, Event              ;   
Int_t               Ntrk, Nprim, NGoodTrk               ;
Float_t             TrueNI                              ;
Double_t            BS[3], BSCovariance[9]              ;
Double_t            HighestPtPV[3] , HighestPtPVCov[9]  ;

// Ntupla Cand
TLorentzVector      BcPi, BcK, MuP, MuM, Pi, Jpsi, JpsiV;

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
Double_t            ClS, BcVtxChi2, BcVtxNdof, BcVtxPosition[3], BcVtxCovariance[9];
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

//Variables for plain tree
Double_t    BcPiM                  = -9999  ; 
Double_t    BcKM                   = -9999  ; 
Double_t    BcPiPt                 = -9999  ; 
Double_t    BcPiP                  = -9999  ; 
Double_t    BcPiEta                = -9999  ; 
Double_t    BcPiPhi                = -9999  ; 
Double_t    BcY                    = -9999  ;  
Double_t    BpY                    = -9999  ;  

Double_t    PiPt                   = -9999  ; 
Double_t    PiPx                   = -9999  ; 
Double_t    PiPy                   = -9999  ; 
Double_t    PiPz                   = -9999  ; 
Double_t    PiEta                  = -9999  ; 
Double_t    PiPhi                  = -9999  ; 

Double_t    MuPPt                  = -9999  ;  
Double_t    MuPPx                  = -9999  ;  
Double_t    MuPPy                  = -9999  ;  
Double_t    MuPPz                  = -9999  ;  
Double_t    MuMPt                  = -9999  ;  
Double_t    MuMPx                  = -9999  ;  
Double_t    MuMPy                  = -9999  ;  
Double_t    MuMPz                  = -9999  ;  
Double_t    MuPEta                 = -9999  ;  
Double_t    MuMEta                 = -9999  ;  

Double_t    JpsiM                  = -9999  ;  
Double_t    JpsiEta                = -9999  ;  
Double_t    JpsiPhi                = -9999  ;  
Double_t    JpsiPt                 = -9999  ;  

Double_t    tau                    = -9999  ;

int CountBc                        = 0      ;
int CountBcDR                      = 0      ;
bool isNewJpsi                     = false  ;
Double_t MaxBcPtArbitration        = -9999  ;                                         
Int_t    fMaxBcPtArbitration       = -9999  ;                                         

