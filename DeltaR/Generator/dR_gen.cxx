// This code implements reconstruction of DeltaR between generated
// electrons using information from the generator. 3 cases are analyzed
// 1) 2e from A + 2e from FSR, in this case we assumed that both e either 
// came from A or from FSR, hence matching two-by-two electrons;
// 2) 4e from 2A, in this case all the combinations of DeltaR were 
// computed and only the smallest values was included in the histograms;
// 3) 4e from FST, similarly to 2), all combinations were evaluated and only
// the smallest DeltaR values were used.
// 

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
#include <vector>
#include "TMath.h"
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

void dR_gen()
{
	TFile *f_output = TFile::Open("ROOTFile_dR_gen.root", "RECREATE");
	TH1D* h_DeltaR_2A2e = new TH1D( "h_DeltaR_2A2e", "", 350, -0.03, 1.04 );
	TH1D* h_DeltaR_4A = new TH1D( "h_DeltaR_4A", "", 350, -0.03, 1.04 );
	TH1D* h_DeltaR_4e = new TH1D( "h_DeltaR_4e", "", 350, -0.03, 1.04 );
	// -- make chain -- //
	TChain *chain = new TChain("recoTree/DYTree");
	chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_HToAATo4e/H_2000GeV_A_50GeV/ntuple_*.root");	//import signal sample

	NtupleHandle *ntuple = new NtupleHandle( chain );
	ntuple->TurnOnBranches_GenLepton();
	ntuple->TurnOnBranches_Electron();

	Int_t nEvent = chain->GetEntries();
	cout << "\t[Total Events: " << nEvent << "]" << endl;
	Double_t Event_4A=0, Event_4e=0, Event_2A2e=0;
	for(Int_t i=0; i<nEvent; i++)
	{
		loadBar(i+1, nEvent, 100, 100);
		
		ntuple->GetEvent(i);
		
		vector<GenLepton> GenLeptonCollection_A;	//collect electrons with mother A
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
		
		     
		if( (Int_t)GenLeptonCollection_A.size() ==2 && (Int_t)GenLeptonCollection_e.size() == 2 )	//if 2e from A and 2e from FSR
		{	
			Event_2A2e++;
			vector<GenLepton> ElectronCollection;
			ElectronCollection.push_back( GenLeptonCollection_A[0] );	//include e's from A and FSR at even/odd positions
			ElectronCollection.push_back( GenLeptonCollection_e[0] );
			ElectronCollection.push_back( GenLeptonCollection_A[1] );
			ElectronCollection.push_back( GenLeptonCollection_e[1] );
			
			for (Int_t j=0; j<2; j++){	//compute DeltaR using TLorentzVector between A or FSR electrons only
					TLorentzVector vec_1,vec_2;
					vec_1.SetPtEtaPhiM(ElectronCollection[j].Pt,ElectronCollection[j].eta,ElectronCollection[j].phi,ElectronCollection[j].Mass);
					vec_2.SetPtEtaPhiM(ElectronCollection[j+2].Pt,ElectronCollection[j+2].eta,ElectronCollection[j+2].phi,ElectronCollection[j+2].Mass);
					h_DeltaR_2A2e->Fill(vec_1.DeltaR(vec_2));		
			}	
		}
		
		//gen-electrons only from dark photon
		if( (Int_t)GenLeptonCollection_A.size() ==4 )
		{	
			Event_4A++;
			vector<GenLepton> ElectronCollection;
			ElectronCollection.push_back( GenLeptonCollection_A[0] );
			ElectronCollection.push_back( GenLeptonCollection_A[1] );
			ElectronCollection.push_back( GenLeptonCollection_A[2] );
			ElectronCollection.push_back( GenLeptonCollection_A[3] );
			vector<Double_t> values;
			
			for(Int_t j=0; j<4;j++){			
				TLorentzVector vec_1;
				vec_1.SetPtEtaPhiM(ElectronCollection[j].Pt,ElectronCollection[j].eta,ElectronCollection[j].phi,ElectronCollection[j].Mass);
				for (Int_t k=j+1;k<4;k++){
					TLorentzVector vec_2;
					vec_2.SetPtEtaPhiM(ElectronCollection[k].Pt,ElectronCollection[k].eta,ElectronCollection[k].phi,ElectronCollection[k].Mass);
					values.push_back(vec_1.DeltaR(vec_2));
				}	
			}
			std::sort (values.begin(),values.end());
			h_DeltaR_4A->Fill(values[0]);
			h_DeltaR_4A->Fill(values[1]);		
		}
		
		//gen-electrons only from FSR
		if( (Int_t)GenLeptonCollection_e.size() ==4 )
		{	
			Event_4e++;
			vector<GenLepton> ElectronCollection;
			ElectronCollection.push_back( GenLeptonCollection_e[0] );
			ElectronCollection.push_back( GenLeptonCollection_e[1] );
			ElectronCollection.push_back( GenLeptonCollection_e[2] );
			ElectronCollection.push_back( GenLeptonCollection_e[3] );
			vector<Double_t> values_e;
			
			for(Int_t j=0; j<4;j++){			
				TLorentzVector vec_1;
				vec_1.SetPtEtaPhiM(ElectronCollection[j].Pt,ElectronCollection[j].eta,ElectronCollection[j].phi,ElectronCollection[j].Mass);
				for (Int_t k=j+1;k<4;k++){
					TLorentzVector vec_2;
					vec_2.SetPtEtaPhiM(ElectronCollection[k].Pt,ElectronCollection[k].eta,ElectronCollection[k].phi,ElectronCollection[k].Mass);
					values_e.push_back(vec_1.DeltaR(vec_2));
				}
			}
			std::sort (values_e.begin(),values_e.end());
			h_DeltaR_4e->Fill(values_e[0]);	
			h_DeltaR_4e->Fill(values_e[1]);	
		}
	}
	
	Double_t f_4A,f_4e,f_2A2e;
	f_4A = Event_4A/(Event_4A + Event_4e + Event_2A2e);
	f_4e = Event_4e/(Event_4A + Event_4e + Event_2A2e);
	f_2A2e = Event_2A2e/(Event_4A + Event_4e + Event_2A2e);
	h_DeltaR_2A2e->SetLineColor(kRed);
	h_DeltaR_4A->SetLineColor(kBlue);
	h_DeltaR_4e->SetLineColor(kGreen);
	h_DeltaR_2A2e->SetStats(0);
	h_DeltaR_4A->SetStats(0);
	h_DeltaR_4e->SetStats(0);
	h_DeltaR_4e->GetXaxis()->SetTitle("DeltaR ");
	h_DeltaR_4e->SetTitle("DeltaR distribution, min DeltaR for H=2000GeV A=50GeV");	
	h_DeltaR_4e->DrawNormalized("",f_4e);
	h_DeltaR_2A2e->DrawNormalized("same",f_2A2e);
	h_DeltaR_4A->DrawNormalized("same",f_4A);
	TLegend *legend = new TLegend(.75, .80, .95, .95);
	legend->AddEntry(h_DeltaR_4e, "mother: 4e ");
	legend->AddEntry(h_DeltaR_4A, "mother: 2A ");
	legend->AddEntry(h_DeltaR_2A2e, " mother: A + 2e");
	legend->Draw();
	
	cout<<f_4A<<", "<<f_4e<<", "<<f_2A2e<<endl;
}
