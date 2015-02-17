/* 
g++ -o ntupleMCSkim ntupleMCSkim.cc dictionary/BcTreeDictionary.cc BcTree.cc -I interface/ -I dictionary/ `root-config --cflags --libs`
*/
// Skim ntuples to plain tree
//

#include "ntupleMCSkim.h"

using namespace std ;

void Analysis(int argc, char** argv) ;
void ClearVariables();
void ClearEventVariables();
void endEvent();
void resizeVectors(int);
void setNullVectors(int);
std::map<std::string, std::string>  ReadNamelist(int argc, char** argv) ;
char * fileName;
//=============================================================
int main ( int argc, char** argv ) {
  
  fileName=argv[1];
  TApplication app("App",&argc, argv);
  Analysis(argc, argv) ;
  cout <<"End of ntupleMCSkim" <<endl;
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
  
  std::cout << " starting " << std::endl;
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
  TH1F*HTrueNI          = new TH1F("HTrueNI"           ,"HTrueNI"                      ,    600,   0.  ,   60. );
  TH1F*HNtrk            = new TH1F("HNtrk"             ,"Ntrk"                         ,   1500,   0.  , 1500. );
  TH1F*HNprim           = new TH1F("HNprim"            ,"Nprim"                        ,     60,   0.  ,  60.  );
  TH1F*HTrueNIW         = new TH1F("HTrueNIW"          ,"HTrueNIW"                     ,    600,   0.  ,   60. );
  TH1F*HNtrkW           = new TH1F("HNtrkW"            ,"NtrkW"                        ,   1500,   0.  , 1500. );
  TH1F*HNprimW          = new TH1F("HNprimW"           ,"NprimW"                       ,     60,   0.  ,  60.  );
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
  TFile * FilePU = TFile::Open("/gwpool/users/fiorendi/Dottorato/2012/CMSSW_5_3_7_patch5/src/HeavyFlavorAnalysis/Bc2JpsiPi/test/weightHisto.root","READ");
  if(!FilePU){
    std::cout << "PU File not found!!!" << std::endl;
    return;
  }
  TH1F* hPU = (TH1F*) FilePU->Get("ComputedWeights");
  if(!hPU){
    std::cout << "PU Histogram not found!!!" << std::endl;
    return;
  }


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

  signalTree->Branch("Ntrk"            ,    &Ntrk                );
  signalTree->Branch("Nprim"           ,    &Nprim               );
  signalTree->Branch("NEve"            ,    &NEve                );
  signalTree->Branch("Run"             ,    &Run                 );
  signalTree->Branch("LBlk"            ,    &LBlk                );
  signalTree->Branch("Event"           ,    &Event               );
  signalTree->Branch("TrueNI"          ,    &TrueNI              );
  signalTree->Branch("w"               ,    &w                   );
  
  signalTree->Branch("TauGen"          ,    &TauGen              );
  signalTree->Branch("BcGenPt"         ,    &BcGenPt             );
  signalTree->Branch("BcGenY"          ,    &BcGenY              );
  signalTree->Branch("xGENPV"          ,    &xGENPV              );
  signalTree->Branch("yGENPV"          ,    &yGENPV              );
  signalTree->Branch("zGENPV"          ,    &zGENPV              );
  signalTree->Branch("xGENSV"          ,    &xGENSV              );
  signalTree->Branch("yGENSV"          ,    &yGENSV              );
  signalTree->Branch("zGENSV"          ,    &zGENSV              );
  
  signalTree->Branch("BcPiM"           ,    &T_BcPiM             );   
  signalTree->Branch("BcKM"            ,    &T_BcKM              );   
  signalTree->Branch("BcPt"            ,    &T_BcPt              );
  signalTree->Branch("BcP"             ,    &T_BcP               );
  signalTree->Branch("BcY"             ,    &T_BcY               );
  signalTree->Branch("BpY"             ,    &T_BpY               );
  
  signalTree->Branch("PiPt"            ,    &T_PiPt              );
  signalTree->Branch("PiPx"            ,    &T_PiPx              );
  signalTree->Branch("PiPy"            ,    &T_PiPy              );
  signalTree->Branch("PiPz"            ,    &T_PiPz              );
  signalTree->Branch("PiEta"           ,    &T_PiEta             );
  signalTree->Branch("PiCh"            ,    &T_PiCh              );
  signalTree->Branch("TrkPixLayerMeas" ,    &T_TrkPixLayerMeas   );
  signalTree->Branch("TrackTrkL"       ,    &T_TrkTrkLayerMeas   );
  signalTree->Branch("TrackPixHits"    ,    &T_TrkPixHits        );
  signalTree->Branch("TrackTrkHits"    ,    &T_TrkTrkHits        );
  signalTree->Branch("TrackNormChi2"   ,    &T_TrkNormChi2       );

  signalTree->Branch("MuPPt"           ,    &T_MuPPt             );
  signalTree->Branch("MuPPx"           ,    &T_MuPPx             );
  signalTree->Branch("MuPPy"           ,    &T_MuPPy             );
  signalTree->Branch("MuPPz"           ,    &T_MuPPz             );
  signalTree->Branch("MuMPt"           ,    &T_MuMPt             );
  signalTree->Branch("MuMPx"           ,    &T_MuMPx             );
  signalTree->Branch("MuMPy"           ,    &T_MuMPy             );
  signalTree->Branch("MuMPz"           ,    &T_MuMPz             );
  signalTree->Branch("MuPEta"          ,    &T_MuPEta            );
  signalTree->Branch("MuMEta"          ,    &T_MuMEta            );
  
  signalTree->Branch("JpsiM"           ,    &T_JpsiM             );
  signalTree->Branch("JpsiEta"         ,    &T_JpsiEta           );
  signalTree->Branch("JpsiPhi"         ,    &T_JpsiPhi           );
  signalTree->Branch("JpsiPt"          ,    &T_JpsiPt            );
  signalTree->Branch("ClJpsi"          ,    &T_ClJpsi            );
  signalTree->Branch("ElsigJpsi"       ,    &T_ElsigJpsi         );
  signalTree->Branch("CosJpsi"         ,    &T_CosJpsi           );
  
  signalTree->Branch("Els2DWrtBS"      ,    &T_Els2DWrtBS        );
  signalTree->Branch("El2DWrtBS"       ,    &T_El2DWrtBS         );
  signalTree->Branch("Els3DWrtPointPV" ,    &T_Els3DWrtPointPV   );
  signalTree->Branch("El3DWrtPointPV"  ,    &T_El3DWrtPointPV    );
  signalTree->Branch("Els3DWrtLongPV"  ,    &T_Els3DWrtLongPV    );
  signalTree->Branch("El3DWrtLongPV"   ,    &T_El3DWrtLongPV     );
  signalTree->Branch("Cos2DWrtBS"      ,    &T_Cos2DWrtBS        );
  signalTree->Branch("Cos3DWrtPointPV" ,    &T_Cos3DWrtPointPV   );
  signalTree->Branch("Cos3DWrtLongPV"  ,    &T_Cos3DWrtLongPV    );
  signalTree->Branch("ClS"             ,    &T_ClS               );
  signalTree->Branch("BcVtxChi2"       ,    &T_BcVtxChi2         );
  signalTree->Branch("BcVtxNdof"       ,    &T_BcVtxNdof         );
  signalTree->Branch("PointPVCl"       ,    &T_PointPVCl         );
  signalTree->Branch("LongPVCl"        ,    &T_LongPVCl          );
  signalTree->Branch("DeltaR"          ,    &T_DeltaR            );
  
  signalTree->Branch("IP3DJpsiSign"    ,    &T_IP3DWrtJpsiSign   );
  signalTree->Branch("IP2DBSSign"      ,    &T_IP2DWrtBSSign     );
  signalTree->Branch("IP3DPVSign"      ,    &T_IP3DWrtPVSign     );
  
  signalTree->Branch("MatchMuP"        ,    &T_MatchMuP          );
  signalTree->Branch("MatchMuM"        ,    &T_MatchMuM          );
  signalTree->Branch("MatchPi"         ,    &T_MatchPi           );
  
  signalTree->Branch("xPointPV"        ,    &T_PointPVPosition_x );
  signalTree->Branch("yPointPV"        ,    &T_PointPVPosition_y );
  signalTree->Branch("zPointPV"        ,    &T_PointPVPosition_z );
  signalTree->Branch("xLongPV"         ,    &T_LongPVPosition_x  );
  signalTree->Branch("yLongPV"         ,    &T_LongPVPosition_y  );
  signalTree->Branch("zLongPV"         ,    &T_LongPVPosition_z  );

//   signalTree->Branch("xxPointPV"       ,    &PointPVCovariance[0]);
//   signalTree->Branch("yyPointPV"       ,    &PointPVCovariance[4]);
//   signalTree->Branch("zzPointPV"       ,    &PointPVCovariance[8]);
//   signalTree->Branch("xxLongPV"        ,    &LongPVCovariance[0] );
//   signalTree->Branch("yyLongPV"        ,    &LongPVCovariance[4] );
//   signalTree->Branch("zzLongPV"        ,    &LongPVCovariance[8] );

  Int_t NBcTreeCand, NumJPsi;
  Int_t nentries = Int_t(TTBcTree->GetEntries());
  std::cout << "Found " << nentries << " Entries" << std::endl;

  Double_t JpsiPxOld, JpsiPyOld, JpsiPzOld;
  TClonesArray *fBcTreeArrayCand = new TClonesArray("BcTreeCand",10000);
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    ClearEventVariables();
//     std::cout << "event:" << eventNo << std::endl;
    Int_t IgetEvent   = TTBcTree->GetEvent(eventNo);
    fBcTreeHeader     = BcTree_ ->GetBcTreeHeader();
    fBcTreeArrayCand  = BcTree_ ->GetBcTreeArrayCand();
    fBcTreeGenCand    = BcTree_ ->GetBcTreeGENCand();
    NBcTreeCand       = BcTree_ ->GetNBcTreeCand();
    TClonesArray      &BcTreeArrayCand = *fBcTreeArrayCand;
    Ntrk              =  fBcTreeHeader->GetNtrk(  );
    NEve              =  fBcTreeHeader->GetNEve(  );
    Run               =  fBcTreeHeader->GetRun(   );
    LBlk              =  fBcTreeHeader->GetLBlk(  );
    Event             =  fBcTreeHeader->GetEvent( );
    TrueNI            =  fBcTreeHeader->GetTrueNI();
    Nprim             =  fBcTreeHeader->GetNprim( );
    
    HTrueNI   -> Fill( TrueNI       );
    HNtrk     -> Fill( Ntrk         );
    HNprim    -> Fill( Nprim        );

    //Get PU weight
    int puBin = hPU->FindBin(TrueNI);
    w = hPU -> GetBinContent(puBin); 

    HTrueNIW  -> Fill( TrueNI, w    );
    HNtrkW    -> Fill( Ntrk  , w    );
    HNprimW   -> Fill( Nprim , w    );
    
    BcGen   = *fBcTreeGenCand -> GetBcPi( );
    BcGenPt = BcGen.Pt();
    BcGenY  = BcGen.Rapidity();

    TauGen  = fBcTreeGenCand -> GetTau(   );
    xGENPV  = fBcTreeGenCand -> GetPVtx( 0);
    yGENPV  = fBcTreeGenCand -> GetPVtx( 1);
    zGENPV  = fBcTreeGenCand -> GetPVtx( 2);
    xGENSV  = fBcTreeGenCand -> GetBcVtx(0);
    yGENSV  = fBcTreeGenCand -> GetBcVtx(1);
    zGENSV  = fBcTreeGenCand -> GetBcVtx(2);
    
    if (NBcTreeCand > 0) {
      resizeVectors(NBcTreeCand);

      for (Int_t countBcCand=0; countBcCand < NBcTreeCand; countBcCand++)            
      {            
        ClearVariables();
        CountBc++;
        BcTreeCand *_BcTreeCand = (BcTreeCand*) BcTreeArrayCand.At(countBcCand);        
        BcPi        =  *_BcTreeCand->GetBcPi()   ;
        BcK         =  *_BcTreeCand->GetBcK()    ;
        JpsiV       =  *_BcTreeCand->GetJpsiV()  ;
        Jpsi        =  *_BcTreeCand->GetJpsi()   ;
        MuP         =  *_BcTreeCand->GetMuP()    ;
        MuM         =  *_BcTreeCand->GetMuM()    ;
        Pi          =  *_BcTreeCand->GetPi()     ;

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
        if (!(_BcTreeCand->GetMuPTMOST()) || !(_BcTreeCand->GetMuMTMOST()) )                        { setNullVectors(countBcCand);  continue;}
        if (isNewJpsi)     HJpsiMassTMost      -> Fill(JpsiV.M());
        if (!(_BcTreeCand->GetMuPTrkLMeas() > 5 ) || !(_BcTreeCand->GetMuMTrkLMeas() > 5) )         { setNullVectors(countBcCand);  continue;}
        if (isNewJpsi)     HJpsiMassTrkMeas    -> Fill(JpsiV.M());
        if (!(_BcTreeCand->GetMuPPixLMeas() > 0 ) || !(_BcTreeCand->GetMuMPixLMeas() > 0) )         { setNullVectors(countBcCand);  continue;}
        if (isNewJpsi)     HJpsiMassPixMeas    -> Fill(JpsiV.M());
        if (!(_BcTreeCand->GetMuPNormChi2() < 3 ) || !(_BcTreeCand->GetMuMNormChi2() < 3) )         { setNullVectors(countBcCand);  continue;}
        if (isNewJpsi)     HJpsiMassChi2       -> Fill(JpsiV.M());
        if (!(_BcTreeCand->GetMuPDxy() < 0.3 ) || !(_BcTreeCand->GetMuPDz() < 20) )                 { setNullVectors(countBcCand);  continue;}
        if (!(_BcTreeCand->GetMuMDxy() < 0.3 ) || !(_BcTreeCand->GetMuMDz() < 20) )                 { setNullVectors(countBcCand);  continue;}
        if (isNewJpsi)     HJpsiMassDxyz     -> Fill(JpsiV.M());
        //L1 requires eta(mu) < 2.1        
        if (!(fabs(MuP.Eta()) < 2.1 ) || !(fabs(MuM.Eta()) < 2.1 ) )                                { setNullVectors(countBcCand);  continue;}
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
        
        
        if (_BcTreeCand->GetTrkTrkHits()  <= 5){ setNullVectors(countBcCand);  continue;}
        if (_BcTreeCand->GetTrkPixHits()  <= 0){ setNullVectors(countBcCand);  continue;}
        if (_BcTreeCand->GetTrkNormChi2() >  3){ setNullVectors(countBcCand);  continue;}
        if (MuP.Pt()                     < 3.5){ setNullVectors(countBcCand);  continue;}
        if (MuM.Pt()                     < 3.5){ setNullVectors(countBcCand);  continue;}

        T_BcPiM            [countBcCand] = BcPi.M()                            ;
        T_BcKM             [countBcCand] = BcK.M()                             ;                           
        T_BcP              [countBcCand] = BcPi.P()                            ;                             
        T_BcPt             [countBcCand] = BcPi.Pt()                           ;                             
        T_BcY              [countBcCand] = BcPi.Rapidity()                     ;                                   
        T_BpY              [countBcCand] = BcK.Rapidity()                      ;                                  
                                                                            
        T_PiPt             [countBcCand] = Pi.Pt()                             ;                           
        T_PiPx             [countBcCand] = Pi.Px()                             ;                           
        T_PiPy             [countBcCand] = Pi.Py()                             ;                           
        T_PiPz             [countBcCand] = Pi.Pz()                             ;                           
        T_PiEta            [countBcCand] = Pi.Eta()                            ;                            
        T_PiPhi            [countBcCand] = Pi.Phi()                            ;                            
        T_PiCh             [countBcCand] = _BcTreeCand->GetPiCh()              ;
        T_TrkPixLayerMeas  [countBcCand] = _BcTreeCand->GetTrkPixLMeas()       ;
        T_TrkTrkLayerMeas  [countBcCand] = _BcTreeCand->GetTrkTrkLMeas()       ;
        T_TrkPixHits       [countBcCand] = _BcTreeCand->GetTrkPixHits()        ;
        T_TrkTrkHits       [countBcCand] = _BcTreeCand->GetTrkTrkHits()        ;
        T_TrkNormChi2      [countBcCand] = _BcTreeCand->GetTrkNormChi2()       ;

        T_MuPPt            [countBcCand] = MuP.Pt()                            ;
        T_MuPPx            [countBcCand] = MuP.Px()                            ;
        T_MuPPy            [countBcCand] = MuP.Py()                            ;
        T_MuPPz            [countBcCand] = MuP.Pz()                            ;
        T_MuMPt            [countBcCand] = MuM.Pt()                            ;                            
        T_MuMPx            [countBcCand] = MuM.Px()                            ;                            
        T_MuMPy            [countBcCand] = MuM.Py()                            ;                            
        T_MuMPz            [countBcCand] = MuM.Pz()                            ;                            
        T_MuPEta           [countBcCand] = MuP.Eta()                           ;                             
        T_MuMEta           [countBcCand] = MuM.Eta()                           ;                             

        T_JpsiEta          [countBcCand] = Jpsi.Eta()                          ;                              
        T_JpsiPhi          [countBcCand] = Jpsi.Phi()                          ;  
        T_JpsiM            [countBcCand] = Jpsi.M()                            ;  
        T_JpsiPt           [countBcCand] = Jpsi.Pt()                           ;
        T_ElsigJpsi        [countBcCand] = _BcTreeCand->GetElsigJpsi()         ;     
        T_ClJpsi           [countBcCand] = _BcTreeCand->GetClJpsi()            ;     
        T_CosJpsi          [countBcCand] = _BcTreeCand->GetCosJpsi()           ;     

        T_MatchMuP         [countBcCand] = _BcTreeCand->GetMatchMuP()          ;
        T_MatchMuM         [countBcCand] = _BcTreeCand->GetMatchMuM()          ;
        T_MatchPi          [countBcCand] = _BcTreeCand->GetMatchPi()           ;

        T_IP3DWrtJpsiSign  [countBcCand] = _BcTreeCand->GetIP3DJpsiSign()      ;
        T_IP2DWrtBSSign    [countBcCand] = _BcTreeCand->GetIP2DBSSign()        ;
        T_IP3DWrtPVSign    [countBcCand] = _BcTreeCand->GetIP3DLongPVSign()    ;

        T_Cos2DWrtBS       [countBcCand] = _BcTreeCand->GetCos2DWrtBS()        ;     
        T_Els2DWrtBS       [countBcCand] = _BcTreeCand->GetEls2DWrtBS()        ;     
        T_El2DWrtBS        [countBcCand] = _BcTreeCand->GetEl2DWrtBS()         ;     
        T_Cos3DWrtPointPV  [countBcCand] = _BcTreeCand->GetCos3DWrtPointPV()   ;     
        T_Els3DWrtPointPV  [countBcCand] = _BcTreeCand->GetEls3DWrtPointPV()   ;     
        T_El3DWrtPointPV   [countBcCand] = _BcTreeCand->GetEl3DWrtPointPV()    ;     
        T_Cos3DWrtLongPV   [countBcCand] = _BcTreeCand->GetCos3DWrtLongPV()    ;     
        T_Els3DWrtLongPV   [countBcCand] = _BcTreeCand->GetEls3DWrtLongPV()    ;     
        T_El3DWrtLongPV    [countBcCand] = _BcTreeCand->GetEl3DWrtLongPV()     ;     
        T_ClS              [countBcCand] = _BcTreeCand->GetClS()               ;          
        T_BcVtxChi2        [countBcCand] = _BcTreeCand->GetBcVtxChi2()         ;          
        T_BcVtxNdof        [countBcCand] = _BcTreeCand->GetBcVtxNdof()         ;          
        T_PointPVCl        [countBcCand] = _BcTreeCand->GetPointPVCl()         ;          
        T_LongPVCl         [countBcCand] = _BcTreeCand->GetLongPVCl()          ;          
        T_DeltaR           [countBcCand] = _BcTreeCand->GetDeltaR()            ;     
  
        T_PointPVPosition_x[countBcCand] = _BcTreeCand->GetPointPVPosition(0)  ;
        T_PointPVPosition_y[countBcCand] = _BcTreeCand->GetPointPVPosition(1)  ;
        T_PointPVPosition_z[countBcCand] = _BcTreeCand->GetPointPVPosition(2)  ;

        T_LongPVPosition_x [countBcCand] = _BcTreeCand->GetLongPVPosition(0)   ;
        T_LongPVPosition_y [countBcCand] = _BcTreeCand->GetLongPVPosition(1)   ;
        T_LongPVPosition_z [countBcCand] = _BcTreeCand->GetLongPVPosition(2)   ;
        
       }// end loop on Bc cand
     } //end if NCand > 0
     else{
       resizeVectors(1);
       setNullVectors(0);
     }
     
    signalTree->Fill();
    endEvent();
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

  BcPi.SetPxPyPzE(0,0,0,0); 
  BcK.SetPxPyPzE(0,0,0,0); 
  MuP.SetPxPyPzE(0,0,0,0); 
  MuM.SetPxPyPzE(0,0,0,0); 
  Jpsi.SetPxPyPzE(0,0,0,0); 
  JpsiV.SetPxPyPzE(0,0,0,0); 
  Pi.SetPxPyPzE(0,0,0,0); 

  MuPGlobal                = false ;
  MuPTracker               = false ;
  MuPPFlow                 = false ;
  MuPTMOST                 = false ;
  MuPTrkLayerMeas          = -9999 ;
  MuPPixLayerMeas          = -9999 ;
  MuPPixHits               = -9999 ;
  MuPTrkHits               = -9999 ;
  MuPMatchedStations       = -9999 ;
  MuPNormChi2              = -9999 ;
  MuPDxy                   = -9999 ;
  MuPDz                    = -9999 ;
   
  MuMGlobal                = false ;
  MuMTracker               = false ;
  MuMPFlow                 = false ;
  MuMTMOST                 = false ;
  MuMTrkLayerMeas          = -9999 ;
  MuMPixLayerMeas          = -9999 ;
  MuMPixHits               = -9999 ;
  MuMTrkHits               = -9999 ;
  MuMMatchedStations       = -9999 ;
  MuMNormChi2              = -9999 ;
  MuMDxy                   = -9999 ;
  MuMDz                    = -9999 ;
   
  DCA                      = -9999 ;
   
  Sigma2DWrtBS             = -9999 ;
  Sigma3DWrtPointPV        = -9999 ;
  Sigma3DWrtLongPV         = -9999 ;

//   for(int j=0; j<3; j++){
    // BcVtxPosition[j]     = -9999 ;
//     JpsiVtxPosition[j]   = -9999 ;
//   }
  for(int j=0; j<9; j++){
    PointPVCovariance[j] = -9999 ;
    LongPVCovariance[j]  = -9999 ;
    BcVtxCovariance[j]   = -9999 ;
  }

//   for(int j=0; j<4; j++){
//     HltMatch[j]          = -9999 ;
//   }
  
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
  TauGen                  = -9999 ;

  for(int j=0; j<3; j++){
    BS[j]                 = -9999 ;
    HighestPtPV[j]        = -9999 ;
  }
  for(int j=0; j<9; j++){
    BSCovariance[j]       = -9999 ;
    HighestPtPVCov[j]     = -9999 ;
  }
  
  CountBc                 =     0 ;
  CountBcDR               =     0 ;
  MaxBcPtArbitration      = -9999.;                                         
  fMaxBcPtArbitration     = -9999 ;  
  
  vector<double>().swap(T_BcPiM              ); 
  vector<double>().swap(T_BcKM               ); 
  vector<double>().swap(T_BcPt               ); 
  vector<double>().swap(T_BcEta              ); 
  vector<double>().swap(T_BcPhi              ); 
  vector<double>().swap(T_BcY                ); 
  vector<double>().swap(T_BpY                ); 
  vector<double>().swap(T_PiPt               ); 
  vector<double>().swap(T_PiPx               ); 
  vector<double>().swap(T_PiPy               ); 
  vector<double>().swap(T_PiPz               ); 
  vector<double>().swap(T_PiEta              ); 
  vector<double>().swap(T_PiPhi              ); 
  vector<double>().swap(T_MuPPt              ); 
  vector<double>().swap(T_MuPPx              ); 
  vector<double>().swap(T_MuPPy              ); 
  vector<double>().swap(T_MuPPz              ); 
  vector<double>().swap(T_MuMPt              ); 
  vector<double>().swap(T_MuMPx              ); 
  vector<double>().swap(T_MuMPy              ); 
  vector<double>().swap(T_MuMPz              ); 
  vector<double>().swap(T_MuPEta             ); 
  vector<double>().swap(T_MuMEta             ); 
  vector<double>().swap(T_JpsiM              ); 
  vector<double>().swap(T_JpsiEta            ); 
  vector<double>().swap(T_JpsiPhi            ); 
  vector<double>().swap(T_JpsiPt             ); 
  vector<double>().swap(T_Tau                ); 

  vector<int>().swap(T_PiCh                  ); 
  vector<int>().swap(T_TrkPixLayerMeas       ); 
  vector<int>().swap(T_TrkTrkLayerMeas       ); 
  vector<int>().swap(T_TrkPixHits            ); 
  vector<int>().swap(T_TrkTrkHits            ); 
  vector<float>().swap(T_TrkNormChi2         ); 

  vector<double>().swap(T_IP3DWrtJpsi        ); 
  vector<double>().swap(T_IP3DWrtJpsiSign    ); 
  vector<double>().swap(T_IP2DWrtBS          ); 
  vector<double>().swap(T_IP2DWrtBSSign      ); 
  vector<double>().swap(T_IP3DWrtPV          ); 
  vector<double>().swap(T_IP3DWrtPVSign      ); 
  vector<double>().swap(T_DeltaR             ); 

  vector<double>().swap(T_ClS                ); 
  vector<double>().swap(T_BcVtxChi2          ); 
  vector<double>().swap(T_BcVtxNdof          ); 
  vector<double>().swap(T_El2DWrtBS          ); 
  vector<double>().swap(T_Els2DWrtBS         ); 
  vector<double>().swap(T_Cos2DWrtBS         ); 
  vector<double>().swap(T_El3DWrtPointPV     ); 
  vector<double>().swap(T_Els3DWrtPointPV    ); 
  vector<double>().swap(T_Cos3DWrtPointPV    ); 

  vector<double>().swap(T_El3DWrtLongPV      ); 
  vector<double>().swap(T_Els3DWrtLongPV     ); 
  vector<double>().swap(T_Cos3DWrtLongPV     ); 
  vector<double>().swap(T_LongPVPosition_x   ); 
  vector<double>().swap(T_LongPVPosition_y   ); 
  vector<double>().swap(T_LongPVPosition_z   ); 
  vector<double>().swap(T_PointPVPosition_x  ); 
  vector<double>().swap(T_PointPVPosition_y  ); 
  vector<double>().swap(T_PointPVPosition_z  ); 

  vector<double>().swap(T_ClJpsi             ); 
  vector<double>().swap(T_ElsigJpsi          ); 
  vector<double>().swap(T_CosJpsi            ); 
  vector<double>().swap(T_PointPVCl          ); 
  vector<double>().swap(T_LongPVCl           ); 

  vector<double>().swap(T_MatchMuP           ); 
  vector<double>().swap(T_MatchMuM           ); 
  vector<double>().swap(T_MatchPi            ); 
  
}
//=========================================================================================
void endEvent() 
{ 
}
//=========================================================================================
void resizeVectors(int NBcTreeCand){
T_Cos2DWrtBS       .resize(NBcTreeCand);
T_Els2DWrtBS       .resize(NBcTreeCand);
T_El2DWrtBS        .resize(NBcTreeCand);
T_Cos3DWrtPointPV  .resize(NBcTreeCand);
T_Els3DWrtPointPV  .resize(NBcTreeCand);
T_El3DWrtPointPV   .resize(NBcTreeCand);
T_Cos3DWrtLongPV   .resize(NBcTreeCand);
T_Els3DWrtLongPV   .resize(NBcTreeCand);
T_El3DWrtLongPV    .resize(NBcTreeCand);
T_ClS              .resize(NBcTreeCand);
T_BcVtxChi2        .resize(NBcTreeCand);
T_BcVtxNdof        .resize(NBcTreeCand);
T_PointPVCl        .resize(NBcTreeCand);
T_LongPVCl         .resize(NBcTreeCand);
T_DeltaR           .resize(NBcTreeCand);

T_BcPiM            .resize(NBcTreeCand);
T_BcKM             .resize(NBcTreeCand);
T_BcP              .resize(NBcTreeCand);
T_BcPt             .resize(NBcTreeCand);
T_BcY              .resize(NBcTreeCand);
T_BpY              .resize(NBcTreeCand);

T_PiPt             .resize(NBcTreeCand);
T_PiPx             .resize(NBcTreeCand);
T_PiPy             .resize(NBcTreeCand);
T_PiPz             .resize(NBcTreeCand);
T_PiEta            .resize(NBcTreeCand);
T_PiPhi            .resize(NBcTreeCand);
T_MuPPt            .resize(NBcTreeCand);
T_MuPPx            .resize(NBcTreeCand);
T_MuPPy            .resize(NBcTreeCand);
T_MuPPz            .resize(NBcTreeCand);
T_MuMPt            .resize(NBcTreeCand);
T_MuMPx            .resize(NBcTreeCand);
T_MuMPy            .resize(NBcTreeCand);
T_MuMPz            .resize(NBcTreeCand);
T_MuPEta           .resize(NBcTreeCand);
T_MuMEta           .resize(NBcTreeCand);
T_JpsiEta          .resize(NBcTreeCand);
T_JpsiPhi          .resize(NBcTreeCand);
T_JpsiM            .resize(NBcTreeCand);
T_JpsiPt           .resize(NBcTreeCand);

T_PiCh             .resize(NBcTreeCand);
T_MatchMuP         .resize(NBcTreeCand);
T_MatchMuM         .resize(NBcTreeCand);
T_MatchPi          .resize(NBcTreeCand);

T_ElsigJpsi        .resize(NBcTreeCand);
T_ClJpsi           .resize(NBcTreeCand);
T_CosJpsi          .resize(NBcTreeCand);

T_TrkPixLayerMeas  .resize(NBcTreeCand);
T_TrkTrkLayerMeas  .resize(NBcTreeCand);
T_TrkPixHits       .resize(NBcTreeCand);
T_TrkTrkHits       .resize(NBcTreeCand);
T_TrkNormChi2      .resize(NBcTreeCand);

T_IP3DWrtJpsiSign  .resize(NBcTreeCand);
T_IP2DWrtBSSign    .resize(NBcTreeCand);
T_IP3DWrtPVSign    .resize(NBcTreeCand);
T_PointPVPosition_x.resize(NBcTreeCand);
T_PointPVPosition_y.resize(NBcTreeCand);
T_PointPVPosition_z.resize(NBcTreeCand);

T_LongPVPosition_x .resize(NBcTreeCand);
T_LongPVPosition_y .resize(NBcTreeCand);
T_LongPVPosition_z .resize(NBcTreeCand);

}
//==============================================
void setNullVectors(int n){
      T_BcPiM            [n] = -9999                               ;
      T_BcKM             [n] = -9999                               ;                           
      T_BcP              [n] = -9999                               ;                             
      T_BcPt             [n] = -9999                               ;                             
      T_BcY              [n] = -9999                               ;                                   
      T_BpY              [n] = -9999                               ;                                  

      T_PiPt             [n] = -9999                               ;                           
      T_PiPx             [n] = -9999                               ;                           
      T_PiPy             [n] = -9999                               ;                           
      T_PiPz             [n] = -9999                               ;                           
      T_PiEta            [n] = -9999                               ;                            
      T_PiPhi            [n] = -9999                               ;                            
      T_PiCh             [n] = -9999                               ;
      T_TrkPixLayerMeas  [n] = -9999                               ;
      T_TrkTrkLayerMeas  [n] = -9999                               ;
      T_TrkPixHits       [n] = -9999                               ;
      T_TrkTrkHits       [n] = -9999                               ;
      T_TrkNormChi2      [n] = -9999                               ;

      T_MuPPt            [n] = -9999                               ;
      T_MuPPx            [n] = -9999                               ;
      T_MuPPy            [n] = -9999                               ;
      T_MuPPz            [n] = -9999                               ;
      T_MuMPt            [n] = -9999                               ;                            
      T_MuMPx            [n] = -9999                               ;                            
      T_MuMPy            [n] = -9999                               ;                            
      T_MuMPz            [n] = -9999                               ;                            
      T_MuPEta           [n] = -9999                               ;                             
      T_MuMEta           [n] = -9999                               ;                             

      T_JpsiEta          [n] = -9999                               ;                              
      T_JpsiPhi          [n] = -9999                               ;  
      T_JpsiM            [n] = -9999                               ;  
      T_JpsiPt           [n] = -9999                               ;
      T_ElsigJpsi        [n] = -9999                               ;     
      T_ClJpsi           [n] = -9999                               ;     
      T_CosJpsi          [n] = -9999                               ;     

      T_MatchMuP         [n] = 9999                                ;
      T_MatchMuM         [n] = 9999                                ;
      T_MatchPi          [n] = 9999                                ;

      T_IP3DWrtJpsiSign  [n] = 9999                                ;
      T_IP2DWrtBSSign    [n] = 9999                                ;
      T_IP3DWrtPVSign    [n] = 9999                                ;

      T_Cos2DWrtBS       [n] = -9999                               ;     
      T_Els2DWrtBS       [n] = -9999                               ;     
      T_El2DWrtBS        [n] = -9999                               ;     
      T_Cos3DWrtPointPV  [n] = -9999                               ;     
      T_Els3DWrtPointPV  [n] = -9999                               ;     
      T_El3DWrtPointPV   [n] = -9999                               ;     
      T_Cos3DWrtLongPV   [n] = -9999                               ;     
      T_Els3DWrtLongPV   [n] = -9999                               ;     
      T_El3DWrtLongPV    [n] = -9999                               ;     
      T_ClS              [n] = -9999                               ;          
      T_BcVtxChi2        [n] = -9999                               ;          
      T_BcVtxNdof        [n] = -9999                               ;          
      T_PointPVCl        [n] = -9999                               ;          
      T_LongPVCl         [n] = -9999                               ;          
      T_DeltaR           [n] =  9999                               ;     

      T_PointPVPosition_x[n] = -9999                               ;
      T_PointPVPosition_y[n] = -9999                               ;
      T_PointPVPosition_z[n] = -9999                               ;

      T_LongPVPosition_x [n] = -9999                               ;
      T_LongPVPosition_y [n] = -9999                               ;
      T_LongPVPosition_z [n] = -9999                               ;
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
std::map<std::string, std::string> ReadNamelist(int argc, char** argv)
{
  vector<string> split( char *str, char c = ' ');
  ifstream indata;
  std::map<std::string, std::string> mappa;
  std::string line;
  vector<string>vstring ;
  indata.open(fileName);

  if(!indata) { // file couldn't be opened
   std::cout << "ntupleMCSkim Error: file could not be opened" << std::endl;
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
