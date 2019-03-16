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

void background()
{
    TFile *f_output = TFile::Open("ROOTFile_background-test2.root", "RECREATE");
    // -- make chain -- //
    TChain *chain = new TChain("recoTree/DYTree");
    //chain->Add("/Users/nanamaan/Documents/ROOT/Files/Background_Sample/ntuple_*.root");    
    //chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_ZZto4L_AOD/ntuple_*.root");    
    chain->Add("/u/user/joonblee/SE_UserHome/ZZTo4L_13TeV_powheg_pythia8/AODToNtuple/170910_120306/0000/ntuple_skim_*.root");

    TTree *bkg_tree = new TTree("background_tree", " updated background tree");

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
    
    bkg_tree->Branch("Electron_etaWidth", &Electron_etaWidth,"Electron_etaWidth/D");
    bkg_tree->Branch("Electron_phiWidth", &Electron_phiWidth,"Electron_phiWidth/D");
    bkg_tree->Branch("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta,"Electron_Full5x5_SigmaIEtaIEta/D");
    bkg_tree->Branch("Electron_HoverE", &Electron_HoverE, "Electron_HoverE/D");
    bkg_tree->Branch("Electron_fbrem", &Electron_fbrem, "Electron_fbrem/D");
    bkg_tree->Branch("Electron_eOverP", &Electron_eOverP, "Electron_eOverP/D");
    bkg_tree->Branch("Electron_InvEminusInvP", &Electron_InvEminusInvP, "Electron_InvEminusInvP/D");
    bkg_tree->Branch("Electron_r9", &Electron_r9, "Electron_r9/D");
    bkg_tree->Branch("Electron_NTTracks_Iso", &Electron_NTTracks_Iso, "Electron_NTTracks_Iso/D");
    bkg_tree->Branch("Electron_NGTracks_Iso", &Electron_NGTracks_Iso, "Electron_NGTracks_Iso/D");


    NtupleHandle *ntuple = new NtupleHandle( chain );
    ntuple->TurnOnBranches_Electron();
    ntuple->TurnOnBranches_TTrack();
    ntuple->TurnOnBranches_GTrack();
    
    Int_t nEvent = chain->GetEntries();
    nEvent = 15000;
    cout << "\t[Total Events: " << nEvent << "]" << endl;
    
    for(Int_t i=0; i<nEvent; i++)
    {

        loadBar(i+1, nEvent, 100, 100);
        ntuple->GetEvent(i);
        
        Int_t nElec = ntuple->Nelectrons;
        
        vector<Electron> ElectronCollection;


        //collect all particles for gsf track in this event
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

        
        //collect information from GSF tracks
        for (Int_t j=0; j< nElec; j++) {
            Electron elec;
            elec.FillFromNtuple(ntuple,j);
            
            Double_t gsf_tracks = 0;
            Double_t gsf_pT = 0;
            Double_t DeltaR_TTrack,eta_TTrack,phi_TTrack;
            //collect all gsf tracks within merged electrons
            for (Int_t k=0; k<(Int_t)TTrackCollection.size();k++) {
                if ( TTrackCollection[k].Pt>30 ) {
                    eta_TTrack = TTrackCollection[k].eta;
                    phi_TTrack = TTrackCollection[k].phi;
                    DeltaR_TTrack = TMath::Sqrt( (eta_TTrack - elec.gsfEta)*(eta_TTrack - elec.gsfEta) + TVector2::Phi_mpi_pi( phi_TTrack - elec.gsfPhi )*TVector2::Phi_mpi_pi( phi_TTrack - elec.gsfPhi ) );
                    if ( DeltaR_TTrack < 0.3  ) {
                        gsf_tracks++;
                        gsf_pT = gsf_pT + TTrackCollection[k].Pt;
                    }
                }
            }
            
            //collect information from general tracks
            Double_t gen_tracks = 0;
            Double_t gen_pT = 0;
            Double_t DeltaR_GTrack,eta_GTrack,phi_GTrack;
            //collect all gen tracks within merged electrons
            for (Int_t k=0; k<(Int_t)GTrackCollection.size();k++) {
                if ( GTrackCollection[k].Pt>30 ) {
                    eta_GTrack = GTrackCollection[k].eta;
                    phi_GTrack = GTrackCollection[k].phi;
                    DeltaR_GTrack = TMath::Sqrt( (eta_GTrack - elec.gsfEta)*(eta_GTrack - elec.gsfEta) + TVector2::Phi_mpi_pi( phi_GTrack - elec.gsfPhi )*TVector2::Phi_mpi_pi( phi_GTrack - elec.gsfPhi ) );
                    if ( DeltaR_GTrack < 0.3  ) {
                        gen_tracks++;
                        gen_pT = gen_pT + GTrackCollection[k].Pt;
                    }
                }
            }
            
                                            
            Electron_etaWidth= elec.etaWidth;
            Electron_phiWidth= elec.phiWidth;
            Electron_Full5x5_SigmaIEtaIEta= elec.Full5x5_SigmaIEtaIEta;
            Electron_HoverE= elec.HoverE;
            Electron_fbrem= elec.fbrem;
            Electron_eOverP= elec.eOverP;
            Electron_InvEminusInvP= elec.InvEminusInvP;
               Electron_r9= elec.r9;
               //if ( gsf_tracks > 0 ) {
               //    gsf_tracks = gsf_pT;        
            //}
            Electron_NTTracks_Iso = gsf_tracks;
               //if ( gen_tracks > 0 ) {
               //    gen_tracks = gen_pT;        
            //}
            Electron_NGTracks_Iso = gen_tracks;                    
            bkg_tree->Fill();                    
        }                    
        
    }
    
    bkg_tree->Write();

}
