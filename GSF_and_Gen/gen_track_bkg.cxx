
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

void gen_track_bkg()
{
	TFile *f_output = TFile::Open("ROOTFile_general_track.root", "RECREATE");
	// -- make chain -- //
	TChain *chain = new TChain("recoTree/DYTree");
	chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_ZZto4L_AOD/ntuple_*.root");	
	
	TH1D *hist_Ngtrack = new TH1D ( "hist_NGT", "Number of GenTrack within DeltaR<0.3 and Pt>30GeV" , 10 , 0 , 10 );
	TH1D *hist_PtGTrack = new TH1D ( "hist_pT", "Pt distribution for gen-tracks within DeltaR<0.3 and Pt>30GeV" , 100 , 0 , 500 );
	
	NtupleHandle *ntuple = new NtupleHandle( chain );
	ntuple->TurnOnBranches_GenLepton();
	ntuple->TurnOnBranches_Electron();
	ntuple->TurnOnBranches_GTrack();
	
	Int_t nEvent = chain->GetEntries();
	cout << "\t[Total Events: " << nEvent << "]" << endl;
	Int_t n=0;
	for(Int_t i=0; i<nEvent; i++)
	{
		loadBar(i+1, nEvent, 100, 100);
		ntuple->GetEvent(i);
				
			//collect all electrons at reconstruction level
		vector<Electron> RecoElectronCollection;
		Int_t nElec = ntuple->Nelectrons;
		for (Int_t i_elec=0; i_elec<nElec; i_elec++) {
			Electron elec;
			elec.FillFromNtuple(ntuple, i_elec);
			RecoElectronCollection.push_back( elec );
		}
			
			//collect all particles for general track in this event
		vector<GTrack> GTrackElectronCollection;
		Int_t nGTrack = ntuple->nGTrack;
		for (Int_t i_gtrack=0; i_gtrack<nGTrack; i_gtrack++) {
			GTrack particle;
			particle.FillFromNtuple(ntuple, i_gtrack);
			GTrackElectronCollection.push_back( particle );
		}
		
		for (Int_t j=0; j<(Int_t)RecoElectronCollection.size();j++) {
			Int_t genTracks=0;
			Double_t DeltaR_GTrack,eta_GTrack,phi_GTrack ;
			for (Int_t k=0; k<(Int_t)GTrackElectronCollection.size();k++) {
				if ( GTrackElectronCollection[k].Pt > 30 ) {
					TLorentzVector vec_genT;
					eta_GTrack = GTrackElectronCollection[k].eta;
					phi_GTrack = GTrackElectronCollection[k].phi;
					DeltaR_GTrack = TMath::Sqrt( (eta_GTrack - RecoElectronCollection[j].gsfEta)*(eta_GTrack - RecoElectronCollection[j].gsfEta) + TVector2::Phi_mpi_pi( phi_GTrack - RecoElectronCollection[j].gsfPhi )*TVector2::Phi_mpi_pi( phi_GTrack - RecoElectronCollection[j].gsfPhi ) );
					if ( DeltaR_GTrack < 0.3 ) {
						genTracks++;
						hist_PtGTrack->Fill ( GTrackElectronCollection[k].Pt );
					}
				}	
			}
			hist_Ngtrack->Fill( genTracks );
		}
	}	
	TCanvas *myC = new TCanvas("canvas");
	myC->Divide(2,1);
	myC->cd(1);
	hist_PtGTrack->DrawNormalized();
	myC->cd(2);
	hist_Ngtrack->DrawNormalized();		
}
