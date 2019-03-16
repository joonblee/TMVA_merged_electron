#include <iostream> // Stream declarations
#include <vector>

#include "TChain.h"
#include "TCut.h"
#include "TH1.h"
#include "TLeaf.h"
#include "TLegend.h"
#include "TFile.h"


#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h" //required to load dataset
#include "TMVA/Reader.h"

using namespace std;

using namespace TMVA;

void ANN(){

  std::string factoryOptions( "!V:!Silent:Transformations=I;D;P;U;G,D" );
  TChain *chain_sig = new TChain("signal_tree");
  chain_sig->Add("/scratch/fs2114/Scripts/Merged_Electrons_Update/ROOTFile_sig_2000_1.root");
  //chain_sig->Add("/Users/nanamaan/Documents/ROOT/Files/Scripts/ROOTFile_gsf_sig_200_1_tracks.root");  
  TChain *chain_bkg = new TChain("background_tree");
  chain_bkg->Add("/scratch/fs2114/Scripts/Background_Sample/ROOTFile_background.root");
  //chain_bkg->Add("/Users/nanamaan/Documents/ROOT/Files/Scripts/ROOTFile_gsf_bkg_tracks.root");


  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  TString outfileName( "TMVA_ANN_output.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );


  TMVA::Factory *factory = new TMVA::Factory( "ANN_MC_signal", outputFile, factoryOptions );
  TMVA::DataLoader *dataloader=new TMVA::DataLoader("ANN_dataset");

  dataloader->AddVariable( "Electron_eOverP", "Electron_eOverP", "", 'D' );
  dataloader->AddVariable( "Electron_HoverE", "Electron_HoverE", "", 'D' );
  dataloader->AddVariable( "Electron_fbrem", "Electron_fbrem", "", 'D' );
  dataloader->AddVariable( "Electron_Full5x5_SigmaIEtaIEta", "Electron_Full5x5_SigmaIEtaIEta", "", 'D' );
  dataloader->AddVariable( "Electron_InvEminusInvP", "Electron_InvEminusInvP", "", 'D' );
  dataloader->AddVariable( "Electron_etaWidth", "Electron_etaWidth", "", 'D' );
  dataloader->AddVariable( "Electron_phiWidth", "Electron_phiWidth", "", 'D' );
  dataloader->AddVariable( "Electron_r9", "Electron_r9", "", 'D' );
  dataloader->AddVariable( "Electron_NTTracks_Iso", "Electron_NTTracks_Iso", "", 'D' );
  dataloader->AddVariable( "Electron_NGTracks_Iso", "Electron_NGTracks_Iso", "", 'D' );
   

  dataloader->AddSignalTree    ( chain_sig ,     signalWeight  );
  dataloader->AddBackgroundTree( chain_bkg , backgroundWeight  );

  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "SplitMode=Random:NormMode=NumEvents:!V" );


  factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "!H:!V:NeuronType=radial:VarTransform=G,D,N:NCycles=500:HiddenLayers=10,5:TestRate=10:UseRegulator=True" );
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();

  delete factory;
  delete dataloader;



  //reader obj
  TMVA::Reader *reader= new TMVA::Reader("!Silent");

  Float_t Electron_eOverP, Electron_HoverE, Electron_InvEminusInvP, Electron_Full5x5_SigmaIEtaIEta,Electron_fbrem,Electron_etaWidth,Electron_phiWidth,Electron_r9,Electron_NGTracks_Iso,Electron_NTTracks_Iso;
  auto className="Signal";
  
  reader->AddVariable("Electron_eOverP", &Electron_eOverP);
  reader->AddVariable("Electron_HoverE", &Electron_HoverE);
  reader->AddVariable("Electron_fbrem", &Electron_fbrem);
  reader->AddVariable("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta);
  reader->AddVariable("Electron_InvEminusInvP", &Electron_InvEminusInvP);
  reader->AddVariable("Electron_etaWidth", &Electron_etaWidth);
  reader->AddVariable("Electron_phiWidth", &Electron_phiWidth);
  reader->AddVariable("Electron_r9", &Electron_r9);
  reader->AddVariable( "Electron_NTTracks_Iso", &Electron_NTTracks_Iso );
  reader->AddVariable( "Electron_NGTracks_Iso", &Electron_NGTracks_Iso );


  reader->BookMVA("ANN method" , "ANN_dataset/weights/ANN_MC_signal_MLP.weights.xml");
  
  TFile *input = TFile::Open("./TMVA_ANN_output.root");
  
  TTree *theTree = (TTree*)input->Get("ANN_dataset/TestTree");
  
  
  theTree->SetBranchAddress("Electron_eOverP", &Electron_eOverP);
  theTree->SetBranchAddress("Electron_HoverE", &Electron_HoverE);
  theTree->SetBranchAddress("Electron_fbrem", &Electron_fbrem);
  theTree->SetBranchAddress("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta);
  theTree->SetBranchAddress("Electron_InvEminusInvP", &Electron_InvEminusInvP);
  theTree->SetBranchAddress("Electron_etaWidth", &Electron_etaWidth);
  theTree->SetBranchAddress("Electron_phiWidth", &Electron_phiWidth);
  theTree->SetBranchAddress("Electron_r9", &Electron_r9);
  theTree->SetBranchAddress( "Electron_NGTracks_Iso", &Electron_NGTracks_Iso );
  theTree->SetBranchAddress( "Electron_NTTracks_Iso", &Electron_NTTracks_Iso );
  //theTree->SetBranchAddress("className", &className);




  TH1F *bg = new TH1F("bgh", "NN output", 70, -0.2, 1.1);
  TH1F *sig = new TH1F("sigh", "NN output", 70, -0.2, 1.1);
  bg->SetDirectory(0);
  sig->SetDirectory(0);
  for (Int_t i = 0; i < theTree->GetEntries(); i++) {
    theTree->GetEntry(i);
    TString event = theTree->GetLeaf("className")->GetValue(0);
    if (event == "S"){
      sig->Fill( reader->EvaluateMVA( "ANN method" ) );
    }
    else {
      bg->Fill( reader->EvaluateMVA( "ANN method" ) );
    }
  } 

  bg->SetLineColor(kBlue);
  bg->SetFillStyle(3008);   bg->SetFillColor(kBlue);
  sig->SetLineColor(kRed);
  sig->SetFillStyle(3003); sig->SetFillColor(kRed);
  bg->SetStats(0);
  sig->SetStats(0);
  bg->DrawNormalized();
  sig->DrawNormalized("same");
  TLegend *legend = new TLegend(.75, .80, .95, .95);
  legend->AddEntry(bg, "Background");
  legend->AddEntry(sig, "Signal");
  legend->Draw();
   

}
