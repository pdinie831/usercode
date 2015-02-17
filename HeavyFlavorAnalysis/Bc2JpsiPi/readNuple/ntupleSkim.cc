/* 
g++ -o ntupleSkim ntupleSkim.cc dictionary/BcTreeDictionary.cc BcTree.cc -I interface/ -I dictionary/ `root-config --cflags --libs`
*/
// Skim ntuples to plain tree
//

#include "ntupleSkim.h"

using namespace std ;

void Analysis(int argc, char** argv) ;
void ClearVariables();
void ClearEventVariables();
std::map<std::string, std::string>  ReadNamelist(int argc, char** argv) ;
char * fileName;
//=============================================================
int main ( int argc, char** argv ) {
  
  fileName=argv[1];
  TApplication app("App",&argc, argv);
  Analysis(argc, argv) ;
  cout <<"End of ntupleSkim" <<endl;
  return 0;
}

//=========================================================================================
// Analysis Routine
void Analysis(int argc, char** argv) {

  std::map<std::string, std::string> mappa = ReadNamelist(argc, argv);
  gROOT ->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1);

  Char_t    ctemp[40] ;
  
  //Input and output files
  strcpy (  InputFile, (mappa["InputFile"]).c_str()  );
  if(mappa.find("OutputFile") != mappa.end()){
   strcpy (  OutputFile, (mappa["OutputFile"]).c_str()  );
  }else{
   strcpy (  OutputFile, "AnalysisBc.root"  );
  } 
  if(mappa.find("TreeFile") != mappa.end()){
   strcpy (  TreeFileName, (mappa["TreeFile"]).c_str()  );
  }else{
   strcpy (  TreeFileName, "treeFile.root"  );
  } 
  sprintf(ctemp,"rm %s",OutputFile);
  cout<<ctemp<<endl;
  gSystem->Exec(ctemp);
  TFile*OutPutFile = new TFile(OutputFile,"RECREATE","ROOT file");

  // Plot ------------------------------------------------------------------------------
  const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);
  TH1::SetDefaultSumw2() ;

  TDirectory * ControlPlots     = OutPutFile  ->mkdir("ControlPlots") ;

  OutPutFile->cd();
  TH1F*HNumBc           = new TH1F("HNumBc"            ,"Num. Cand Bc."                ,   2000,   0.  , 2000. );
  TH1F*HNumBcDR         = new TH1F("HNumBcDR"          ,"Num. Cand Bc.DR"              ,   2000,   0.  , 2000. );
  TH1F*HNtrk            = new TH1F("HNtrk"             ,"Ntrk"                         ,   1500,   0.  , 1500. );
  TH1F*HNprim           = new TH1F("HNprim"            ,"Nprim"                        ,     60,   0.  ,  60.  );
  TH1F*HNumJPsi         = new TH1F("HNumJPsi"          ,"NumJPsi"                      ,     10,   0.  ,  10.  );

  ControlPlots->cd(); 
  TH1F*HJpsiMass        = new TH1F("HJpsiMass"         ,"Jpsi Mass"                    ,   1000,  2.5  , 3.5   );
  TH1F*HBcPiMass        = new TH1F("HBcPiMass"         ,"BcPi Mass"                    ,   7000,   3.  ,  10.  );
  TH1F*HBcKMass         = new TH1F("HBcKMass"          ,"BcK Mass"                     ,   7000,   3.  ,  10.  );

  TH1F*HJpsiMassTMost   = new TH1F("HJpsiMassTMost"    ,"Jpsi Mass if TMOST"           ,   1000,  2.5  , 3.5   );
  TH1F*HJpsiMassTrkMeas = new TH1F("HJpsiMassTrkMeas"  ,"Jpsi Mass if trk layers"      ,   1000,  2.5  , 3.5   );
  TH1F*HJpsiMassPixMeas = new TH1F("HJpsiMassPixMeas"  ,"Jpsi Mass if pix layers"      ,   1000,  2.5  , 3.5   );
  TH1F*HJpsiMassChi2    = new TH1F("HJpsiMassChi2"     ,"Jpsi Mass if chi2"            ,   1000,  2.5  , 3.5   );
  TH1F*HJpsiMassDxyz    = new TH1F("HJpsiMassDxyz"     ,"Jpsi Mass if dxy and dz"      ,   1000,  2.5  , 3.5   );
  TH1F*HJpsiMassEta     = new TH1F("HJpsiMassEta"      ,"Jpsi Mass if eta mu < 2.1"    ,   1000,  2.5  , 3.5   );

  TH1F*HTrackPixL       = new TH1F("HTrackPixL"        ,"Number of Pix layers"         ,     20,    0  , 20    );
  TH1F*HTrackTrkL       = new TH1F("HTrackTrkL"        ,"Number of Trk layers"         ,     20,    0  , 20    );
  TH1F*HTrackPixHits    = new TH1F("HTrackPixHits"     ,"Number of Pix hits"           ,     20,    0  , 20    );
  TH1F*HTrackTrkHits    = new TH1F("HTrackTrkHits"     ,"Number of Trk hits"           ,     20,    0  , 20    );
  TH1F*HTrackNormChi2   = new TH1F("HTrackNormChi2"    ,"Trk chi2"                     ,    100,    0  , 10    );
  TH1F*HTrackDeltaR     = new TH1F("HTrackDeltaR"      ,"Trk DR"                       ,    100,    0  , 10    );
  TH1F*HTrackIP3DJpsi   = new TH1F("HTrackIP3DJpsi"    ,"Trk IP sign from Jpsi,3D"     ,    100,    0  , 10    );
  TH1F*HTrackIP2DBS     = new TH1F("HTrackIP2DBS"      ,"Trk IP sign from BS,2D"       ,    100,    0  , 10    );
  TH1F*HTrackIP3DPV     = new TH1F("HTrackIP3DPV"      ,"Trk IP sign from PV,3D"       ,    100,    0  , 10    );

  TH1F*HBcVtxClS        = new TH1F("HBcVtxClS"         ,"CL for Bc vtx"                ,   1000,    0  , 1     );
  TH1F*HBcPt            = new TH1F("HBcPt"             ,"Pt of the Bc cand"            ,   1200,    0  , 120   );
  TH1F*HMuPPt           = new TH1F("HMuPPt"            ,"Pt of the MuP cand"           ,   1200,    0  , 120   );
  TH1F*HMuMPt           = new TH1F("HMuMPt"            ,"Pt of the MuM cand"           ,   1200,    0  , 120   );
  TH1F*HTrackPt         = new TH1F("HTrackPt"          ,"Pt of the pion cand"          ,   1200,    0  , 120   );

  OutPutFile->cd(); 
  TH1::AddDirectory(oldAddDir);

//end plots --------------------------------------------------------------------------

  TFile * FileBc = TFile::Open(InputFile,"READ");
  if(!FileBc){
    std::cout << "File: "<< InputFile << " not found!!!" << std::endl;
    return;
  }

  FileBc->cd("BcTreeDir");
  TTree *TTBcTree = (TTree*)gROOT->FindObject("BcTree;1");
  if (TTBcTree == 0){
     std::cout << "TTree: not found!!!" << std::endl;
     exit(1);
  } 
  BcTree *BcTree_ = new BcTree();
  TTBcTree->SetBranchAddress("BcBr", &BcTree_);

  TFile *TreeFile   = new TFile(TreeFileName,"RECREATE", "ROOT file");
  TTree *signalTree = new TTree("signalTree","Bc mass tree");

  signalTree->Branch("BcPiM"           ,    &BcPiM           );   
  signalTree->Branch("BcKM"            ,    &BcKM            );   
  signalTree->Branch("BcPiPt"          ,    &BcPiPt          );
  signalTree->Branch("BcPiP"           ,    &BcPiP           );
  signalTree->Branch("BcY"             ,    &BcY             );
  signalTree->Branch("BpY"             ,    &BpY             );

  signalTree->Branch("PiPt"            ,    &PiPt            );
  signalTree->Branch("PiPx"            ,    &PiPx            );
  signalTree->Branch("PiPy"            ,    &PiPy            );
  signalTree->Branch("PiPz"            ,    &PiPz            );
  signalTree->Branch("PiEta"           ,    &PiEta           );

  signalTree->Branch("MuPPt"           ,    &MuPPt           );
  signalTree->Branch("MuPPx"           ,    &MuPPx           );
  signalTree->Branch("MuPPy"           ,    &MuPPy           );
  signalTree->Branch("MuPPz"           ,    &MuPPz           );
  signalTree->Branch("MuMPt"           ,    &MuMPt           );
  signalTree->Branch("MuMPx"           ,    &MuMPx           );
  signalTree->Branch("MuMPy"           ,    &MuMPy           );
  signalTree->Branch("MuMPz"           ,    &MuMPz           );
  signalTree->Branch("MuPEta"          ,    &MuPEta          );
  signalTree->Branch("MuMEta"          ,    &MuMEta          );

  signalTree->Branch("JpsiM"           ,    &JpsiM           );
  signalTree->Branch("JpsiEta"         ,    &JpsiEta         );
  signalTree->Branch("JpsiPhi"         ,    &JpsiPhi         );
  signalTree->Branch("JpsiPt"          ,    &JpsiPt          );

  signalTree->Branch("Els2DWrtBS"      ,    &Els2DWrtBS      );
  signalTree->Branch("El2DWrtBS"       ,    &El2DWrtBS       );
  signalTree->Branch("Els3DWrtPointPV" ,    &Els3DWrtPointPV );
  signalTree->Branch("El3DWrtPointPV"  ,    &El3DWrtPointPV  );
  signalTree->Branch("Els3DWrtLongPV"  ,    &Els3DWrtLongPV  );
  signalTree->Branch("El3DWrtLongPV"   ,    &El3DWrtLongPV   );
  signalTree->Branch("Cos2DWrtBS"      ,    &Cos2DWrtBS      );
  signalTree->Branch("Cos3DWrtPointPV" ,    &Cos3DWrtPointPV );
  signalTree->Branch("Cos3DWrtLongPV"  ,    &Cos3DWrtLongPV  );
  signalTree->Branch("ClS"             ,    &ClS             );
  signalTree->Branch("BcVtxChi2"       ,    &BcVtxChi2       );
  signalTree->Branch("BcVtxNdof"       ,    &BcVtxNdof       );
  signalTree->Branch("PointPVCl"       ,    &PointPVCl       );
  signalTree->Branch("LongPVCl"        ,    &LongPVCl        );
  signalTree->Branch("DeltaR"          ,    &DeltaR          );

  signalTree->Branch("IP3DJpsiSign"    ,    &IP3DWrtJpsiSign );
  signalTree->Branch("IP2DBSSign"      ,    &IP2DWrtBSSign   );
  signalTree->Branch("IP3DPVSign"      ,    &IP3DWrtPVSign   );

  signalTree->Branch("MatchMuP"        ,    &MatchMuP        );
  signalTree->Branch("MatchMuM"        ,    &MatchMuM        );
  signalTree->Branch("MatchPi"         ,    &MatchPi         );

  signalTree->Branch("ClJpsi"          ,    &ClJpsi          );
  signalTree->Branch("ElsigJpsi"       ,    &ElsigJpsi       );
  signalTree->Branch("CosJpsi"         ,    &CosJpsi         );
  signalTree->Branch("PiCh"            ,    &PiCh            );
  signalTree->Branch("Ntrk"            ,    &Ntrk            );
  signalTree->Branch("Nprim"           ,    &Nprim           );
  signalTree->Branch("NEve"            ,    &NEve            );
  signalTree->Branch("Run"             ,    &Run             );
  signalTree->Branch("LBlk"            ,    &LBlk            );
  signalTree->Branch("Event"           ,    &Event           );
  signalTree->Branch("TrueNI"          ,    &TrueNI          );

  signalTree->Branch("TrkPixLayerMeas" ,    &TrkPixLayerMeas );
  signalTree->Branch("TrackTrkL"       ,    &TrkTrkLayerMeas );
  signalTree->Branch("TrackPixHits"    ,    &TrkPixHits      );
  signalTree->Branch("TrackTrkHits"    ,    &TrkTrkHits      );
  signalTree->Branch("TrackNormChi2"   ,    &TrkNormChi2     );


  Int_t NBcTreeCand, NumJPsi;
  Int_t nentries = Int_t(TTBcTree->GetEntries());
  std::cout << "Found " << nentries << " Entries" << std::endl;

  Double_t JpsiPxOld, JpsiPyOld, JpsiPzOld;
  TClonesArray *fBcTreeArrayCand = new TClonesArray("BcTreeCand",10000);
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    ClearEventVariables();
    Int_t IgetEvent   = TTBcTree->GetEvent(eventNo);
    fBcTreeHeader     = BcTree_ ->GetBcTreeHeader();
    fBcTreeArrayCand  = BcTree_ ->GetBcTreeArrayCand();
    NBcTreeCand       = BcTree_ ->GetNBcTreeCand();
    TClonesArray &BcTreeArrayCand = *fBcTreeArrayCand;
    Ntrk     =  fBcTreeHeader->GetNtrk(  );
    NEve     =  fBcTreeHeader->GetNEve(  );
    Run      =  fBcTreeHeader->GetRun(   );
    LBlk     =  fBcTreeHeader->GetLBlk(  );
    Event    =  fBcTreeHeader->GetEvent( );
    TrueNI   =  fBcTreeHeader->GetTrueNI();
    Nprim    =  fBcTreeHeader->GetNprim( );
    HNtrk    -> Fill( Ntrk      );
    HNprim   -> Fill( Nprim     );

//=================================================================================
    for (Int_t countBcCand=0; countBcCand < NBcTreeCand; countBcCand++)            
    {
      ClearVariables();
      CountBc++;
      BcTreeCand *_BcTreeCand = (BcTreeCand*) BcTreeArrayCand.At(countBcCand);        
      BcPi        =  *_BcTreeCand->GetBcPi()    ;
      BcK         =  *_BcTreeCand->GetBcK()     ;
      JpsiV       =  *_BcTreeCand->GetJpsiV()   ;
      MuP         =  *_BcTreeCand->GetMuP()     ;
      MuM         =  *_BcTreeCand->GetMuM()     ;
      Pi          =  *_BcTreeCand->GetPi()      ;

      //Check if this is a new Jpsi or the same of last cand       
      if(JpsiV.Px() != JpsiPxOld && JpsiV.Py()!=JpsiPyOld && JpsiV.Pz()!=JpsiPzOld)
      {
        JpsiPxOld = JpsiV.Px();
        JpsiPyOld = JpsiV.Py();
        JpsiPzOld = JpsiV.Pz();
        NumJPsi++;
        isNewJpsi = true;
      } 

      if (isNewJpsi)    HJpsiMass -> Fill(JpsiV.M());
      HBcPiMass -> Fill(BcPi.M());
      HBcKMass  -> Fill(BcK.M());
       
      //New Soft muon requirements
      if (!(_BcTreeCand->GetMuPTMOST()) || !(_BcTreeCand->GetMuMTMOST()) )                        continue;
      if (isNewJpsi)     HJpsiMassTMost      -> Fill(JpsiV.M());
      if (!(_BcTreeCand->GetMuPTrkLMeas() > 5 ) || !(_BcTreeCand->GetMuMTrkLMeas() > 5) )         continue;
      if (isNewJpsi)     HJpsiMassTrkMeas    -> Fill(JpsiV.M());
      if (!(_BcTreeCand->GetMuPPixLMeas() > 0 ) || !(_BcTreeCand->GetMuMPixLMeas() > 0) )         continue;
      if (isNewJpsi)     HJpsiMassPixMeas    -> Fill(JpsiV.M());
      if (!(_BcTreeCand->GetMuPNormChi2() < 3 ) || !(_BcTreeCand->GetMuMNormChi2() < 3) )         continue;
      if (isNewJpsi)     HJpsiMassChi2       -> Fill(JpsiV.M());
      if (!(_BcTreeCand->GetMuPDxy() < 0.3 ) || !(_BcTreeCand->GetMuPDz() < 20) )                 continue;
      if (!(_BcTreeCand->GetMuMDxy() < 0.3 ) || !(_BcTreeCand->GetMuMDz() < 20) )                 continue;
      if (isNewJpsi)     HJpsiMassDxyz     -> Fill(JpsiV.M());
      //L1 requires eta(mu) < 2.1        
      if (!(fabs(MuP.Eta()) < 2.1 ) || !(fabs(MuM.Eta()) < 2.1 ) )                                continue;
      if (isNewJpsi)     HJpsiMassEta     -> Fill(JpsiV.M());


      HTrackPixL        ->Fill(_BcTreeCand->GetTrkPixLMeas()    );
      HTrackTrkL        ->Fill(_BcTreeCand->GetTrkTrkLMeas()    );
      HTrackPixHits     ->Fill(_BcTreeCand->GetTrkPixHits()     );
      HTrackTrkHits     ->Fill(_BcTreeCand->GetTrkTrkHits()     );
      HTrackNormChi2    ->Fill(_BcTreeCand->GetTrkNormChi2()    );
      HTrackDeltaR      ->Fill(_BcTreeCand->GetDeltaR()         );
      HTrackIP3DJpsi    ->Fill(_BcTreeCand->GetIP3DJpsiSign()   );
      HTrackIP2DBS      ->Fill(_BcTreeCand->GetIP2DBSSign()     );
      HTrackIP3DPV      ->Fill(_BcTreeCand->GetIP3DLongPVSign() );
      
      HBcVtxClS         ->Fill(_BcTreeCand->GetClS()            );
      HBcPt             ->Fill(BcPi.Pt()                        );
      HTrackPt          ->Fill(Pi.Pt()                          );
      if (isNewJpsi)    HMuPPt    ->Fill(MuP.Pt()               );
      if (isNewJpsi)    HMuMPt    ->Fill(MuM.Pt()               );

      if (_BcTreeCand->GetTrkTrkHits()  <= 5){  continue;}
      if (_BcTreeCand->GetTrkPixHits()  <= 0){  continue;}
      if (_BcTreeCand->GetTrkNormChi2() >  3){  continue;}
      
      if(BcPi.Pt() < 10 ) continue;
      if(BcPi.M()  < 4.5) continue;
      if(BcPi.M()  > 7.5) continue;
      if(MuP.Pt()  < 3.5) continue;
      if(MuM.Pt()  < 3.5) continue;
       
      Jpsi             =  *_BcTreeCand->GetJpsi()     ;

      Cos2DWrtBS       =  _BcTreeCand->GetCos2DWrtBS()        ;     
      Els2DWrtBS       =  _BcTreeCand->GetEls2DWrtBS()        ;     
      El2DWrtBS        =  _BcTreeCand->GetEl2DWrtBS()         ;     
      Cos3DWrtPointPV  =  _BcTreeCand->GetCos3DWrtPointPV()   ;     
      Els3DWrtPointPV  =  _BcTreeCand->GetEls3DWrtPointPV()   ;     
      El3DWrtPointPV   =  _BcTreeCand->GetEl3DWrtPointPV()    ;     
      Cos3DWrtLongPV   =  _BcTreeCand->GetCos3DWrtLongPV()    ;     
      Els3DWrtLongPV   =  _BcTreeCand->GetEls3DWrtLongPV()    ;     
      El3DWrtLongPV    =  _BcTreeCand->GetEl3DWrtLongPV()     ;     

      ClS              =  _BcTreeCand->GetClS()               ;          
      BcVtxChi2        =  _BcTreeCand->GetBcVtxChi2()         ;          
      BcVtxNdof        =  _BcTreeCand->GetBcVtxNdof()         ;          
      PointPVCl        =  _BcTreeCand->GetPointPVCl()         ;          
      LongPVCl         =  _BcTreeCand->GetLongPVCl()          ;          
      DeltaR           =  _BcTreeCand->GetDeltaR()            ;     

      BcPiM            =  BcPi.M()                            ;
      BcKM             =  BcK.M()                             ;                           
      BcPiP            =  BcPi.P()                            ;                             
      BcPiPt           =  BcPi.Pt()                           ;                             
      BcY              =  BcPi.Rapidity()                     ;                                   
      BpY              =  BcK.Rapidity()                      ;                                  
                                                            
      PiPt             =  Pi.Pt()                             ;                           
      PiPx             =  Pi.Px()                             ;                           
      PiPy             =  Pi.Py()                             ;                           
      PiPz             =  Pi.Pz()                             ;                           
      PiEta            =  Pi.Eta()                            ;                            
      PiPhi            =  Pi.Phi()                            ;                            
      MuPPt            =  MuP.Pt()                            ;
      MuPPx            =  MuP.Px()                            ;
      MuPPy            =  MuP.Py()                            ;
      MuPPz            =  MuP.Pz()                            ;
      MuMPt            =  MuM.Pt()                            ;                            
      MuMPx            =  MuM.Px()                            ;                            
      MuMPy            =  MuM.Py()                            ;                            
      MuMPz            =  MuM.Pz()                            ;                            
      MuPEta           =  MuP.Eta()                           ;                             
      MuMEta           =  MuM.Eta()                           ;                             
      JpsiEta          =  Jpsi.Eta()                          ;                              
      JpsiPhi          =  Jpsi.Phi()                          ;  
      JpsiM            =  JpsiV.M()                           ;  
      JpsiPt           =  Jpsi.Pt()                           ;
        
      PiCh             =  _BcTreeCand->GetPiCh()              ;
      MatchMuP         =  _BcTreeCand->GetMatchMuP()          ;
      MatchMuM         =  _BcTreeCand->GetMatchMuM()          ;
      MatchPi          =  _BcTreeCand->GetMatchPi()           ;

      ElsigJpsi        =  _BcTreeCand->GetElsigJpsi()         ;     
      ClJpsi           =  _BcTreeCand->GetClJpsi()            ;     
      CosJpsi          =  _BcTreeCand->GetCosJpsi()           ;     
    
      TrkPixLayerMeas  =    _BcTreeCand->GetTrkPixLMeas()     ;
      TrkTrkLayerMeas  =    _BcTreeCand->GetTrkTrkLMeas()     ;
      TrkPixHits       =    _BcTreeCand->GetTrkPixHits()      ;
      TrkTrkHits       =    _BcTreeCand->GetTrkTrkHits()      ;
      TrkNormChi2      =    _BcTreeCand->GetTrkNormChi2()     ;
          
      IP3DWrtJpsiSign  =   _BcTreeCand->GetIP3DJpsiSign()     ;
      IP2DWrtBSSign    =   _BcTreeCand->GetIP2DBSSign()       ;
      IP3DWrtPVSign    =   _BcTreeCand->GetIP3DLongPVSign()   ;

//       BcVtxPosition[0]  =   _BcTreeCand->GetBcVtxPosition(0)  ;
//       BcVtxPosition[1]  =   _BcTreeCand->GetBcVtxPosition(1)  ;
//       BcVtxPosition[2]  =   _BcTreeCand->GetBcVtxPosition(2)  ;
// 
//       PointPVPosition[0]=   _BcTreeCand->GetPointPVPosition(0)  ;
//       PointPVPosition[1]=   _BcTreeCand->GetPointPVPosition(1)  ;
//       PointPVPosition[2]=   _BcTreeCand->GetPointPVPosition(2)  ;
//       
//       BS[0]              =   fBcTreeHeader->GetBS(0)  ;
//       BS[1]              =   fBcTreeHeader->GetBS(1)  ;
//       BS[2]              =   fBcTreeHeader->GetBS(2)  ;
    
        signalTree->Fill();
    } //end loop candidates

    HNumJPsi -> Fill(NumJPsi);
    HNumBc   -> Fill(CountBc);
    HNumBcDR -> Fill(CountBcDR);

  }//end loop events

  TreeFile   ->cd();
  signalTree ->Write();    
  TreeFile   ->Write();
  OutPutFile ->cd();
  OutPutFile ->Write();
  cout<<"Write OutPutFile : "<< OutputFile << ", " << TreeFileName <<endl;
  exit(1);
}
//=========================================================================================
void ClearVariables() 
{ 
  BcPiM                  = -9999 ;
  BcKM                   = -9999 ;
  BcPiPt                 = -9999 ;
  BcPiP                  = -9999 ;
  BcPiEta                = -9999 ;
  BcPiPhi                = -9999 ;
  BcY                    = -9999 ;
  BpY                    = -9999 ;
  PiPt                   = -9999 ;
  PiPx                   = -9999 ;
  PiPy                   = -9999 ;
  PiPz                   = -9999 ;
  PiEta                  = -9999 ;
  PiPhi                  = -9999 ;
  MuPPt                  = -9999 ;
  MuPPx                  = -9999 ;
  MuPPy                  = -9999 ;
  MuPPz                  = -9999 ;
  MuMPt                  = -9999 ;
  MuMPx                  = -9999 ;
  MuMPy                  = -9999 ;
  MuMPz                  = -9999 ;
  MuPEta                 = -9999 ;
  MuMEta                 = -9999 ;
  JpsiM                  = -9999 ;
  JpsiEta                = -9999 ;
  JpsiPhi                = -9999 ;
  JpsiPt                 = -9999 ;
  tau                    = -9999 ;

  BcPi .SetPxPyPzE(0,0,0,0); 
  BcK  .SetPxPyPzE(0,0,0,0); 
  MuP  .SetPxPyPzE(0,0,0,0); 
  MuM  .SetPxPyPzE(0,0,0,0); 
  Jpsi .SetPxPyPzE(0,0,0,0); 
  JpsiV.SetPxPyPzE(0,0,0,0); 
  Pi   .SetPxPyPzE(0,0,0,0); 

  MuPGlobal              = false ;
  MuPTracker             = false ;
  MuPPFlow               = false ;
  MuPTMOST               = false ;
  MuPTrkLayerMeas        = -9999 ;
  MuPPixLayerMeas        = -9999 ;
  MuPPixHits             = -9999 ;
  MuPTrkHits             = -9999 ;
  MuPMatchedStations     = -9999 ;
  MuPNormChi2            = -9999 ;
  MuPDxy                 = -9999 ;
  MuPDz                  = -9999 ;

  MuMGlobal              = false ;
  MuMTracker             = false ;
  MuMPFlow               = false ;
  MuMTMOST               = false ;
  MuMTrkLayerMeas        = -9999 ;
  MuMPixLayerMeas        = -9999 ;
  MuMPixHits             = -9999 ;
  MuMTrkHits             = -9999 ;
  MuMMatchedStations     = -9999 ;
  MuMNormChi2            = -9999 ;
  MuMDxy                 = -9999 ;
  MuMDz                  = -9999 ;

  DCA                    = -9999 ;

  PiCh                   = -9999 ;
  TrkPixLayerMeas        = -9999 ;
  TrkTrkLayerMeas        = -9999 ;
  TrkPixHits             = -9999 ;
  TrkTrkHits             = -9999 ;
  TrkNormChi2            = -9999 ;
  IP3DWrtJpsi            = -9999 ;
  IP3DWrtJpsiSign        = -9999 ;
  IP2DWrtBS              = -9999 ;
  IP2DWrtBSSign          = -9999 ;
  IP3DWrtPV              = -9999 ;
  IP3DWrtPVSign          = -9999 ;
  DeltaR                 = -9999 ;

  ClS                    = -9999 ;
  BcVtxChi2              = -9999 ;
  BcVtxNdof              = -9999 ;
  El2DWrtBS              = -9999 ;
  Els2DWrtBS             = -9999 ;
  Sigma2DWrtBS           = -9999 ;
  Cos2DWrtBS             = -9999 ;
  El3DWrtPointPV         = -9999 ;
  Els3DWrtPointPV        = -9999 ;
  Sigma3DWrtPointPV      = -9999 ;
  Cos3DWrtPointPV        = -9999 ;
  El3DWrtLongPV          = -9999 ;
  Els3DWrtLongPV         = -9999 ;
  Sigma3DWrtLongPV       = -9999 ;
  Cos3DWrtLongPV         = -9999 ;
  for(int j=0; j<3; j++){
    BcVtxPosition[j]     = -9999 ;
    JpsiVtxPosition[j]   = -9999 ;
    PointPVPosition[j]   = -9999 ;
    LongPVPosition[j]    = -9999 ;
  }
  for(int j=0; j<9; j++){
    PointPVCovariance[j] = -9999 ;
    LongPVCovariance[j]  = -9999 ;
    BcVtxCovariance[j]   = -9999 ;
  }

  ClJpsi                 = -9999 ;
  ElsigJpsi              = -9999 ;
  CosJpsi                = -9999 ;
  PointPVCl              = -9999 ;
  LongPVCl               = -9999 ;

  MatchMuP               = -9999 ;
  MatchMuM               = -9999 ;
  MatchPi                = -9999 ;

  for(int j=0; j<4; j++){
    HltMatch[j]          = -9999 ;
  }
  
  isNewJpsi              = false ;
}
//=========================================================================================
void ClearEventVariables() 
{ 

  NEve                    =  9999 ; 
  Run                     =  9999 ; 
  LBlk                    =  9999 ; 
  Event                   =  9999 ; 
  Ntrk                    = -9999 ;
  Nprim                   = -9999 ;
  NGoodTrk                = -9999 ;
  TrueNI                  = -9999 ;
  for(int j=0; j<3; j++){
    BS[j]                 = -9999 ;
    HighestPtPV[j]        = -9999 ;
  }
  for(int j=0; j<9; j++){
    BSCovariance[j]       = -9999 ;
    HighestPtPVCov[j]     = -9999 ;
  }
  
  CountBc                 = 0;
  CountBcDR               = 0;
  MaxBcPtArbitration      = -9999.;                                         
  fMaxBcPtArbitration     = -9999 ;                                         

}

//=========================================================================================
vector<string> split( char *str, char c = ' ')
{
  vector<string> result;
  while(1)
  {
    char *begin = str;
    while(*str != c && *str)
      str++;
    result.push_back(string(begin, str));
    if(0 == *str++)     break;
  }
  return result;
}
//=========================================================================================
// Read parameters from namelist file
std::map<std::string, std::string> ReadNamelist(int argc, char** argv){
  vector<string> split( char *str, char c = ' ');
  ifstream indata;
  std::map<std::string, std::string> mappa;
  std::string line;
  vector<string>vstring ;
  indata.open(fileName);

  if(!indata) { // file couldn't be opened
   std::cout << "ntupleSkim Error: file could not be opened" << std::endl;
   exit(1);
  }
  while(std::getline(indata, line)) {
    line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), ' ' ), line.end());
    char *cstr = new char [line.size()+1];
    strcpy (cstr, line.c_str());
    vstring = split(cstr,'=');
    mappa.insert( std::pair<string,string>(vstring[0],vstring[1]) );
  }
  indata.close();    
  return mappa ;
}
