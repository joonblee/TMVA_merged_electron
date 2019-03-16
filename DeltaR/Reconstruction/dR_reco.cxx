#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
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
#include "TMath.h"
#include <vector>
#include <algorithm>

#include "NtupleHandle.h"
#include "Object.h"

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
 

	cout << "]\r" << flush;
}

void dR_reco()
{
	TFile *f_output = TFile::Open("ROOTFile_reco_dR_200_1_.root", "RECREATE");
	TH1D* h_dR_sig = new TH1D( "h_dR", "", 200, -0.03, 5 );
	TH1F* h_Elec = new TH1F( "h_Elec", "", 10, -0.03, 10 );

	// -- make chain -- //
	TChain *chain = new TChain("recoTree/DYTree");
	// chain->Add("./ntuple_*.root");
	chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_HToAATo4e/H_200GeV_A_1GeV/ntuple_*.root");

	NtupleHandle *ntuple = new NtupleHandle( chain );
	ntuple->TurnOnBranches_Electron();

	Int_t nEvent = chain->GetEntries();
	cout << "\t[Total Events: " << nEvent << "]" << endl;
	for(Int_t i=0; i<nEvent; i++)
	{
	
		loadBar(i+1, nEvent, 100, 100);
		
		ntuple->GetEvent(i);
		
		vector<Electron> ElectronCollection;
		
		// collect all electrons at reconstruction level
		Int_t nElec = ntuple->Nelectrons;
		for(Int_t i_elec=0; i_elec<nElec; i_elec++)
		{
			Electron elec;
			elec.FillFromNtuple( ntuple, i_elec );
			ElectronCollection.push_back( elec );
		}
		
		h_Elec->Fill(ElectronCollection.size());
		
		//compute all possible DeltaR values between reconstructed electrons
		for (Int_t j=0; j<ElectronCollection.size(); j++) {
			TLorentzVector vec_1;
			vec_1.SetPtEtaPhiE(ElectronCollection[j].Pt,ElectronCollection[j].eta,ElectronCollection[j].phi,ElectronCollection[j].Energy);
			for (Int_t k= j+1; k<ElectronCollection.size(); k++ ) {
				TLorentzVector vec_2;
				vec_2.SetPtEtaPhiE(ElectronCollection[k].Pt,ElectronCollection[k].eta,ElectronCollection[k].phi,ElectronCollection[k].Energy);
				h_dR_sig->Fill(vec_2.DeltaR(vec_1));
			}
		}
	}	
	
		
	TCanvas* canvas = new TCanvas("canvas","");
	canvas->Divide(2,1);	
	canvas->cd(1);		
	h_dR_sig->GetXaxis()->SetTitle("|DeltaR|");
	h_dR_sig->SetTitle("DeltaR Distribution, H=200Gev, A=1GeV");			
	h_dR_sig->SetLineColor(kRed);
	h_dR_sig->SetStats(0);
	h_dR_sig->DrawNormalized();
	TLegend *legend = new TLegend(.75, .80, .95, .95);
	legend->AddEntry(h_dR_sig, "Signal");
	legend->Draw();
	canvas->cd(2);
	h_Elec->GetXaxis()->SetTitle("Number of e ");
	h_Elec->SetTitle("Number of reconstructed Electrons per event");	
	h_Elec->DrawNormalized();
	
	f_output->Close();
	
}
