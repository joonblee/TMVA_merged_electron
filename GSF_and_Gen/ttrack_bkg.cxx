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

void ttrack_bkg()
{
	TFile *f_output = TFile::Open("ROOTFile_ttrack_bkg.root", "RECREATE");
	// -- make chain -- //
	TChain *chain = new TChain("recoTree/DYTree");
	chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_ZZto4L_AOD/ntuple_*.root");	
		
	TH1D *hist_pT = new TH1D ( "hist_pT gsf tracks", " pT distribution of highest pT gsf tracks within DeltaR<0.3 from reco-e" , 200 , 0 , 700 );
	TH1D *hist_tracks = new TH1D ( "hist_tracks", "gsf tracks within DeltaR<0.3 from reco-e" , 10 , 0 , 10 );

	NtupleHandle *ntuple = new NtupleHandle( chain );
	ntuple->TurnOnBranches_Electron();
	ntuple->TurnOnBranches_TTrack();
	
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
		
			
		//collect all particles for gsf track in this event
		vector<TTrack> TTrackElectronCollection;
		Int_t nTTrack = ntuple->NTT;
		for (Int_t i_ttrack=0; i_ttrack<nTTrack; i_ttrack++) {
			TTrack particle;
			particle.FillFromNtuple(ntuple, i_ttrack);
			TTrackElectronCollection.push_back( particle );
		}
			
		for (Int_t j=0; j<(Int_t)RecoElectronCollection.size();j++) {
			Int_t TTracks=0;
			Double_t DeltaR_TTrack,eta_TTrack,phi_TTrack ;
			for (Int_t k=0; k<(Int_t)TTrackElectronCollection.size();k++) {
				eta_TTrack = TTrackElectronCollection[k].eta;
				phi_TTrack = TTrackElectronCollection[k].phi;
				DeltaR_TTrack = TMath::Sqrt( (eta_TTrack - RecoElectronCollection[j].gsfEta)*(eta_TTrack - RecoElectronCollection[j].gsfEta) + TVector2::Phi_mpi_pi( phi_TTrack - RecoElectronCollection[j].gsfPhi )*TVector2::Phi_mpi_pi( phi_TTrack - RecoElectronCollection[j].gsfPhi ) );
				if ( DeltaR_TTrack < 0.3 ) {
					TTracks++;
					hist_pT->Fill ( TTrackElectronCollection[k].Pt );
				}	
			}
			hist_tracks->Fill( TTracks );
		}
	}	
	TCanvas *myC = new TCanvas("canvas");
	myC->Divide(2,1);
	myC->cd(1);
	hist_pT->DrawNormalized();
	myC->cd(2);
	hist_tracks->DrawNormalized();

}
