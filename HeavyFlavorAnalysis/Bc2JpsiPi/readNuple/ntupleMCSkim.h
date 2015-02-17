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

BcTree              *newBcTree_ ;
BcTreeHeader     *fBcTreeHeader ;
BcTreeGENCand    *fBcTreeGenCand;

std::map<std::string, TTree*     > t1_ ;


// Ntupla header 
UInt_t              NEve, Run, LBlk, Event              ;   
Int_t               Ntrk, Nprim, NGoodTrk               ;
Float_t             TrueNI                              ;
Double_t            w                                   ;
Double_t            BS[3], BSCovariance[9]              ;
Double_t            HighestPtPV[3] , HighestPtPVCov[9]  ;

// Ntupla Cand
TLorentzVector      BcPi, BcK, MuP, MuM, Pi, Jpsi, JpsiV, BcGen;

// Muons - not saved in plain ntuple
Float_t             DCA ;
bool                MuPGlobal, MuPTracker, MuPPFlow, MuPTMOST, MuMGlobal, MuMTracker, MuMPFlow, MuMTMOST;
Int_t               MuPTrkLayerMeas, MuPPixLayerMeas, MuPPixHits, MuPTrkHits, MuPMatchedStations;
Int_t               MuMTrkLayerMeas, MuMPixLayerMeas, MuMPixHits, MuMTrkHits, MuMMatchedStations;
Float_t             MuPNormChi2, MuPDxy, MuPDz, MuMNormChi2, MuMDxy, MuMDz;

// Bc vertex 
Double_t           BcVtxPosition[3], BcVtxCovariance[9];
Double_t           Sigma2DWrtBS     ;
Double_t           Sigma3DWrtPointPV;
Double_t           Sigma3DWrtLongPV ;

// Jpsi vertex 
Double_t            JpsiVtxPosition[3];

// Primary vertices 
Double_t            PointPVPosition[3], PointPVCovariance[9];
Double_t            LongPVPosition[3],  LongPVCovariance[9] ;

// HLT matching 
Double_t            HltMatch[4];

//GEN quantities
Double_t            xGENPV, yGENPV, zGENPV, xGENSV, yGENSV, zGENSV  ;
Double_t            TauGen, BcGenPt, BcGenY;



//Variables for plain tree
std::vector<Double_t>    T_BcPiM, T_BcKM,    T_BcPt,    T_BcP,   T_BcEta, T_BcPhi, T_BcY, T_BpY; 
std::vector<Double_t>    T_PiPt,  T_PiPx,    T_PiPy,    T_PiPz,  T_PiEta,   T_PiPhi; 
std::vector<Double_t>    T_MuPPt, T_MuPEta,  T_MuMPt,   T_MuMEta; 
std::vector<Double_t>    T_MuPPx, T_MuPPy,   T_MuPPz,   T_MuMPx, T_MuMPy,   T_MuMPz; 
std::vector<Double_t>    T_JpsiM, T_JpsiEta, T_JpsiPhi, T_JpsiPt; 
std::vector<Double_t>    T_Tau; 

std::vector<Int_t>       T_PiCh, T_TrkPixLayerMeas, T_TrkTrkLayerMeas, T_TrkPixHits, T_TrkTrkHits; 
std::vector<Float_t>     T_TrkNormChi2; 
std::vector<Double_t>    T_IP3DWrtJpsi, T_IP3DWrtJpsiSign, T_IP2DWrtBS, T_IP2DWrtBSSign, T_IP3DWrtPV, T_IP3DWrtPVSign, T_DeltaR;

// Bc vertex 
std::vector<Double_t>     T_ClS,            T_BcVtxChi2,        T_BcVtxNdof;// BcVtxPosition[3], BcVtxCovariance[9];
std::vector<Double_t>     T_El2DWrtBS,      T_Els2DWrtBS,       T_Cos2DWrtBS;
std::vector<Double_t>     T_El3DWrtPointPV, T_Els3DWrtPointPV,  T_Cos3DWrtPointPV;
std::vector<Double_t>     T_El3DWrtLongPV,  T_Els3DWrtLongPV,   T_Cos3DWrtLongPV;

std::vector<Double_t>     T_ClJpsi, T_ElsigJpsi, T_CosJpsi;//, JpsiVtxPosition[3];
std::vector<Double_t>     T_PointPVCl, T_LongPVCl;

// MC matching 
std::vector<Double_t>     T_MatchMuP, T_MatchMuM, T_MatchPi;   

std::vector<Double_t>     T_PointPVPosition_x, T_PointPVPosition_y, T_PointPVPosition_z; 
std::vector<Double_t>     T_LongPVPosition_x,  T_LongPVPosition_y,  T_LongPVPosition_z; 

int CountBc                       = 0      ;
int CountBcDR                     = 0      ;
bool isNewJpsi                    = false  ;
Double_t MaxBcPtArbitration       = -9999  ;                                         
Int_t    fMaxBcPtArbitration      = -9999  ;                                         



