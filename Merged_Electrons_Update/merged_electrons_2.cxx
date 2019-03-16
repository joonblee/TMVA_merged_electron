#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TLegend.h>
#include <THStack.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TAttMarker.h>
#include <TF1.h>
#include <TStyle.h>
#include <TEfficiency.h>
#include <vector>
#include "TMath.h"
#include "TVector2.h"
#include "NtupleHandle.h"
#include "Object.h"



Bool_t CompareObject( Object lep1, Object lep2 )
{
  return lep1.Pt > lep2.Pt;
}



static inline void loadBar(int x, int n, int r, int w)
{
   // Only update r times.
   if( x == n )
     cout << endl;

   if ( x % (n/r +1) != 0 ) return;

 
   // Calculuate the ratio of complete-to-incomplete.
   float ratio = x/(float)n;
   int   c     = ratio * w;
 
   // Show the percentage complete.
   printf("%3d%% [", (int)(ratio*100) );
 
   // Show the load bar.
   for (int x=0; x<c; x++) cout << "=";
 
   for (int x=c; x<w; x++) cout << " ";
 
   // ANSI Control codes to go back to the
   // previous line and clear it.
   cout << "]\r" << flush;
}

void merged_electrons_2()
{
  TFile *f_output = TFile::Open("ROOTFile_sig_2000_1-test.root", "RECREATE");
  // -- make chain -- //
  TChain *chain = new TChain("recoTree/DYTree");
  //chain->Add("/Users/nanamaan/Documents/ROOT/Files/Signal_Sample/ntuple_skim_200_1.root");
  //chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_HToAATo4e/H_2000GeV_A_1GeV/ntuple_*.root");
  chain->Add("/u/user/joonblee/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/v01_HAA4e/H_2000GeV_A_1GeV/ntuple_skim.root");

  TTree *sig_tree = new TTree("signal_tree", " updated signal tree");

  Double_t Electron_NTTracks_Iso;
  Double_t Electron_NGTracks_Iso;
  Double_t Electron_etaWidth;
  Double_t Electron_phiWidth;   
  Double_t Electron_Full5x5_SigmaIEtaIEta;
  Double_t Electron_HoverE;
  Double_t Electron_fbrem;
  Double_t Electron_eOverP;
  Double_t Electron_InvEminusInvP;
  Double_t Electron_r9;   
   
  sig_tree->Branch("Electron_etaWidth", &Electron_etaWidth,"Electron_etaWidth/D");
  sig_tree->Branch("Electron_phiWidth", &Electron_phiWidth,"Electron_phiWidth/D");
  sig_tree->Branch("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta,"Electron_Full5x5_SigmaIEtaIEta/D");
  sig_tree->Branch("Electron_HoverE", &Electron_HoverE, "Electron_HoverE/D");
  sig_tree->Branch("Electron_fbrem", &Electron_fbrem, "Electron_fbrem/D");
  sig_tree->Branch("Electron_eOverP", &Electron_eOverP, "Electron_eOverP/D");
  sig_tree->Branch("Electron_InvEminusInvP", &Electron_InvEminusInvP, "Electron_InvEminusInvP/D");
  sig_tree->Branch("Electron_r9", &Electron_r9, "Electron_r9/D");
  sig_tree->Branch("Electron_NTTracks_Iso", &Electron_NTTracks_Iso, "Electron_NTTracks_Iso/D");
  sig_tree->Branch("Electron_NGTracks_Iso", &Electron_NGTracks_Iso, "Electron_NGTracks_Iso/D");
   

  NtupleHandle *ntuple = new NtupleHandle( chain );
  ntuple->TurnOnBranches_GenLepton();
  ntuple->TurnOnBranches_Electron();
  ntuple->TurnOnBranches_TTrack();
  ntuple->TurnOnBranches_GTrack();
   
  Int_t nEvent = chain->GetEntries();
  nEvent = 1000;
  cout << "\t[Total Events: " << nEvent << "]" << endl;
   
  for(Int_t i=0; i<nEvent; i++)
  {

    loadBar(i+1, nEvent, 100, 100);
    ntuple->GetEvent(i);
      
    //store merged electrons during each event
    vector<Electron> merged_electrons;
      
    //collect all the relevant electrons from generator level
    vector<GenLepton> GenLeptonCollection_A;   //collect electrons with mother A
    vector<GenLepton> GenLeptonCollection_e; //collect electrons from FSR
    Int_t NGenLeptons = ntuple->gnpair;
      
    for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
    {
      GenLepton genlep;
      genlep.FillFromNtuple(ntuple, i_gen);
      if( fabs(genlep.Mother) == 36 && genlep.fromHardProcessFinalState == 1 && fabs(genlep.ID) == 11) {
        GenLepton elec;
        elec.FillFromNtuple(ntuple,i_gen);
        GenLeptonCollection_A.push_back( elec );
      }
      if( fabs(genlep.Mother) == 11 && genlep.fromHardProcessFinalState == 1 && fabs(genlep.ID) ==11) {
        GenLepton elec;
        elec.FillFromNtuple(ntuple,i_gen);
        GenLeptonCollection_e.push_back( elec );
      }
    }
      
      
    vector<GenLepton> ElectronCollection;
      
    //compute mother_pt (A) for all electrons in the generator level:
    //1) case 4e from 2A:
    if( (Int_t)GenLeptonCollection_A.size() ==4 )
    {   
      ElectronCollection.push_back( GenLeptonCollection_A[0] );
      ElectronCollection.push_back( GenLeptonCollection_A[1] );
      ElectronCollection.push_back( GenLeptonCollection_A[2] );
      ElectronCollection.push_back( GenLeptonCollection_A[3] );
    }
      
      
    // case of: 2e from A + 2e from FSR / reconstruct mother_pt 
    if( (Int_t)GenLeptonCollection_A.size() ==2 && (Int_t)GenLeptonCollection_e.size() == 2 )
    {   
      ElectronCollection.push_back( GenLeptonCollection_A[0] );   //include e's from A and FSR at even/odd positions
      ElectronCollection.push_back( GenLeptonCollection_e[0] );
      ElectronCollection.push_back( GenLeptonCollection_A[1] );
      ElectronCollection.push_back( GenLeptonCollection_e[1] );
    }
      
    //case of 4e from FSR:
      
    if( (Int_t)GenLeptonCollection_e.size() == 4 )
    {         
      ElectronCollection.push_back( GenLeptonCollection_e[0] );
      ElectronCollection.push_back( GenLeptonCollection_e[1] );
      ElectronCollection.push_back( GenLeptonCollection_e[2] );
      ElectronCollection.push_back( GenLeptonCollection_e[3] );
    }      
      
    Bool_t Flag_Is4Gen_e = kFALSE;
         
    if ((Int_t)ElectronCollection.size() == 4 ) Flag_Is4Gen_e = kTRUE;
      
    if ( Flag_Is4Gen_e ) {
      //collect all electrons at reconstruction level
      vector<Electron> RecoElectronCollection;
      Int_t nElec = ntuple->Nelectrons;
      for (Int_t i_elec=0; i_elec<nElec; i_elec++) {
        Electron elec;
        elec.FillFromNtuple(ntuple, i_elec);
        RecoElectronCollection.push_back( elec );
      }

      //collect all particles for tracker track in this event
      vector<TTrack> TTrackCollection;
      Int_t nTTrack = ntuple->NTT;
      for (Int_t i_ttrack=0; i_ttrack<nTTrack; i_ttrack++) {
        TTrack particle;
        particle.FillFromNtuple(ntuple, i_ttrack);
        TTrackCollection.push_back( particle );
      }
            
      //collect all particles for general track in this event
      vector<GTrack> GTrackCollection;
      Int_t nGTrack = ntuple->nGTrack;
      for (Int_t i_gtrack=0; i_gtrack<nGTrack; i_gtrack++) {
        GTrack particle;
        particle.FillFromNtuple(ntuple, i_gtrack);
        GTrackCollection.push_back( particle );
      }
                  
            
      // fill in merged_electrons vectors will all electrons falling into DeltaR<0.3 cone
         

      TLorentzVector vec_1,vec_2,vec_3,vec_4;
      vec_1.SetPtEtaPhiM(ElectronCollection[0].Pt,ElectronCollection[0].eta,ElectronCollection[0].phi,ElectronCollection[0].Mass);
      vec_2.SetPtEtaPhiM(ElectronCollection[1].Pt,ElectronCollection[1].eta,ElectronCollection[1].phi,ElectronCollection[1].Mass);
      vec_3.SetPtEtaPhiM(ElectronCollection[2].Pt,ElectronCollection[2].eta,ElectronCollection[2].phi,ElectronCollection[2].Mass);
      vec_4.SetPtEtaPhiM(ElectronCollection[3].Pt,ElectronCollection[3].eta,ElectronCollection[3].phi,ElectronCollection[3].Mass);      

      //make sure each reco-e contains 2gen-e within DeltaR<0.3
      for (Int_t k=0; k<(Int_t)RecoElectronCollection.size();k++){
        
        Double_t dEt_ratio,Et_gen1,Et_gen2,Et_gen3,Et_gen4;
        Et_gen1 = vec_1.Et();
        Et_gen2 = vec_2.Et();
        Et_gen3 = vec_3.Et();
        Et_gen4 = vec_4.Et();         
         
        Double_t DeltaR1,DeltaR2,DeltaR3,DeltaR4;
        TLorentzVector vec_reco;
        vec_reco.SetPtEtaPhiE(RecoElectronCollection[k].Pt,RecoElectronCollection[k].eta,RecoElectronCollection[k].phi,RecoElectronCollection[k].Energy);
        DeltaR1 = vec_1.DeltaR(vec_reco);
        DeltaR2 = vec_2.DeltaR(vec_reco);
        DeltaR3 = vec_3.DeltaR(vec_reco);
        DeltaR4 = vec_4.DeltaR(vec_reco);
        //collect merged electrons with missing transverse energy <10%
        if ( DeltaR1 < 0.3 ) {
          if ( DeltaR2 < 0.3  ) {
            dEt_ratio = fabs( RecoElectronCollection[k].etSC - Et_gen1 - Et_gen2 )/(RecoElectronCollection[k].etSC);
            if ( dEt_ratio <0.1 ) {
              merged_electrons.push_back( RecoElectronCollection[k] );
            }
          }
          if ( DeltaR3 < 0.3  ) {
            dEt_ratio = fabs( RecoElectronCollection[k].etSC - Et_gen1 - Et_gen3 )/(RecoElectronCollection[k].etSC);
            if ( dEt_ratio <0.1 ) {
              merged_electrons.push_back( RecoElectronCollection[k] );
            }
          }
          if ( DeltaR4 < 0.3  ) {
            dEt_ratio = fabs( RecoElectronCollection[k].etSC - Et_gen1 - Et_gen4 )/(RecoElectronCollection[k].etSC);
            if ( dEt_ratio <0.1 ) {
              merged_electrons.push_back( RecoElectronCollection[k] );
            }
          }
        }
        if ( DeltaR2 < 0.3 ) {
          if ( DeltaR3 < 0.3  ) {
            dEt_ratio = fabs( RecoElectronCollection[k].etSC - Et_gen2 - Et_gen3 )/(RecoElectronCollection[k].etSC);
            if ( dEt_ratio <0.1 ) {
              merged_electrons.push_back( RecoElectronCollection[k] );
            }
          }
          if ( DeltaR4 < 0.3  ) {
            dEt_ratio = fabs( RecoElectronCollection[k].etSC - Et_gen2 - Et_gen4 )/(RecoElectronCollection[k].etSC);
            if ( dEt_ratio <0.1 ) {
              merged_electrons.push_back( RecoElectronCollection[k] );
            }
          }
        }
        if ( DeltaR3 < 0.3 ) {
          if ( DeltaR4 < 0.3  ) {
            dEt_ratio = fabs( RecoElectronCollection[k].etSC - Et_gen3 - Et_gen4 )/(RecoElectronCollection[k].etSC);
            if ( dEt_ratio <0.1 ) {
              merged_electrons.push_back( RecoElectronCollection[k] );
            }
          }   
        }         
      }
            
         
                  
      //cout<<(Int_t)chain_bkg->GetEntries()<<endl;
      for ( Int_t k=0; k<(Int_t)merged_electrons.size();k++) {
        Double_t gsf_tracks = 0;
        Double_t gsf_pT = 0;
        Double_t DeltaR_TTrack,eta_TTrack,phi_TTrack;
            
        //collect all gsf tracks within merged electrons
        for (Int_t j=0; j<(Int_t)TTrackCollection.size();j++) {
          if ( TTrackCollection[j].Pt>30 ) {
            eta_TTrack = TTrackCollection[j].eta;
            phi_TTrack = TTrackCollection[j].phi;
            DeltaR_TTrack = TMath::Sqrt( (eta_TTrack - merged_electrons[k].gsfEta)*(eta_TTrack - merged_electrons[k].gsfEta) + TVector2::Phi_mpi_pi( phi_TTrack - merged_electrons[k].gsfPhi )*TVector2::Phi_mpi_pi( phi_TTrack - merged_electrons[k].gsfPhi ) );
            if ( DeltaR_TTrack < 0.3  ) {
              gsf_tracks++;
              gsf_pT = gsf_pT + TTrackCollection[j].Pt;
            }
          }
        }
            
        //collect all general tracks within merged electrons
        Double_t gen_tracks = 0;
        Double_t gen_pT = 0;
        Double_t DeltaR_GTrack,eta_GTrack,phi_GTrack;
            
        for (Int_t j=0; j<(Int_t)GTrackCollection.size();j++) {
          if ( GTrackCollection[j].Pt>30 ) {
            eta_GTrack = GTrackCollection[j].eta;
            phi_GTrack = GTrackCollection[j].phi;
            DeltaR_GTrack = TMath::Sqrt( (eta_GTrack - merged_electrons[k].gsfEta)*(eta_GTrack - merged_electrons[k].gsfEta) + TVector2::Phi_mpi_pi( phi_GTrack - merged_electrons[k].gsfPhi )*TVector2::Phi_mpi_pi( phi_GTrack - merged_electrons[k].gsfPhi ) );
            if ( DeltaR_GTrack < 0.3  ) {
              gen_tracks++;
              gen_pT = gen_pT + GTrackCollection[j].Pt;
            }
          }
        }            
            
        Electron_etaWidth= merged_electrons[k].etaWidth;
        Electron_phiWidth= merged_electrons[k].phiWidth;
        Electron_Full5x5_SigmaIEtaIEta= merged_electrons[k].Full5x5_SigmaIEtaIEta;
        Electron_HoverE= merged_electrons[k].HoverE;
        Electron_fbrem= merged_electrons[k].fbrem;
        Electron_eOverP= merged_electrons[k].eOverP;
        Electron_InvEminusInvP= merged_electrons[k].InvEminusInvP;
        Electron_r9= merged_electrons[k].r9;
        //if ( gsf_tracks > 0 ) {
        //   gsf_tracks = gsf_pT;      
        //}
        Electron_NTTracks_Iso = gsf_tracks;
        //if ( gen_tracks > 0 ) {
        //   gen_tracks = gen_pT;      
        //}
        Electron_NGTracks_Iso = gen_tracks;                     
        sig_tree->Fill();
      }               
    }      
  }
   
  sig_tree->Write();

}
