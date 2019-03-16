// Compute DeltaR distribution using mother_pT variable. all mother's pT
// are found looping on those electrons coming from FSR. Once all mothers
// have been found, their pT are matched in order to compute the DeltaR value 
// between electron pairs.
//
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

void dR_gen_motherPt()
{
	TFile *f_output = TFile::Open("ROOTFile_dR_gen_motherPt.root", "RECREATE");
	//TH1D* h_DeltaR_2A2e = new TH1D( "h_DeltaR_2A2e", "", 350, -0.00, 0.1 );
	//TH1D* h_DeltaR_4A = new TH1D( "h_DeltaR_4A", "", 350, -0.0, 0.5 );
	//TH1D* h_DeltaR_4e = new TH1D( "h_DeltaR_4e", "", 350, -0.0, 0.5 );
	TH1D* h_DeltaR = new TH1D( "h_DeltaR", "", 350, -0.0, 1. );
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
		
		vector<GenLepton> NonFinalLeptonCollection;
		vector<GenLepton> GenLeptonCollection_A;	//collect electrons with mother A
		vector<Double_t> MotherPtCollection_A;
		vector<GenLepton> GenLeptonCollection_e; //collect electrons from FSR
		Int_t NGenLeptons = ntuple->gnpair;
		vector<Double_t> MotherPtCollection_e;
		
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			if( fabs(genlep.Mother) == 36 && genlep.fromHardProcessFinalState == 1 && fabs(genlep.ID) == 11) {
				GenLepton elec;
				elec.FillFromNtuple(ntuple,i_gen);
				GenLeptonCollection_A.push_back( elec );
				MotherPtCollection_A.push_back( elec.Mother_Pt );
			}
			if( fabs(genlep.Mother) == 11 && genlep.fromHardProcessFinalState == 1 && fabs(genlep.ID) ==11) {
				GenLepton elec;
				elec.FillFromNtuple(ntuple,i_gen);
				GenLeptonCollection_e.push_back( elec );
				MotherPtCollection_e.push_back( elec.Mother_Pt );
			}
			if( fabs(genlep.ID) == 11 && genlep.isPromptFinalState == 0 ){ // -- electron (ID==11) & not final state  -- //
				NonFinalLeptonCollection.push_back( genlep );
			}	
		}
		
		
		if( (Int_t)GenLeptonCollection_A.size() ==4 )
		{	
			Event_4A++;
			vector<GenLepton> ElectronCollection;
			
			ElectronCollection.push_back( GenLeptonCollection_A[0] );
			ElectronCollection.push_back( GenLeptonCollection_A[1] );
			ElectronCollection.push_back( GenLeptonCollection_A[2] );
			ElectronCollection.push_back( GenLeptonCollection_A[3] );
			
			vector<Double_t> MotherPtCollection;
			MotherPtCollection.push_back( MotherPtCollection_A[0] );
			MotherPtCollection.push_back( MotherPtCollection_A[1] );
			MotherPtCollection.push_back( MotherPtCollection_A[2] );
			MotherPtCollection.push_back( MotherPtCollection_A[3] );
			
			//compute DeltaR
			for(Int_t j=0; j<ElectronCollection.size();j++){			
				TLorentzVector vec_1;
				vec_1.SetPtEtaPhiM(ElectronCollection[j].Pt,ElectronCollection[j].eta,ElectronCollection[j].phi,ElectronCollection[j].Mass);
				for (Int_t k=j+1;k<4;k++){
					if ( MotherPtCollection[j] == MotherPtCollection[k] ) {
						TLorentzVector vec_2;
						vec_2.SetPtEtaPhiM(ElectronCollection[k].Pt,ElectronCollection[k].eta,ElectronCollection[k].phi,ElectronCollection[k].Mass);
						h_DeltaR->Fill(vec_1.DeltaR(vec_2));
					}
					
				}		
			}
		}
		
		//case for 2e from FSR and 2e from A
		
		if( (Int_t)GenLeptonCollection_A.size() ==2 && (Int_t)GenLeptonCollection_e.size() == 2 )
		{	
			Event_2A2e++;
			
			vector<Int_t> mother_reconstructed;	//keep track of electrons with reconstructed mother (e.g. entry = 1 )
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			
			// loop until all electrons are reconstructed from A particles
			while ( mother_reconstructed[0] != 1 && mother_reconstructed[1] != 1 )
			{
				for (Int_t j=0; j<NonFinalLeptonCollection.size(); j++)
				{
					for (Int_t k=0; k<GenLeptonCollection_e.size();k++)
						{
						if ( NonFinalLeptonCollection[j].Pt == MotherPtCollection_e[k] ){
							MotherPtCollection_e[k] = NonFinalLeptonCollection[j].Mother_Pt;
							if ( fabs(NonFinalLeptonCollection[j].Mother) == 36 ) {
								mother_reconstructed[k] = 1;
							}
						}
					}
				}
			}
			

			vector<GenLepton> ElectronCollection;
			ElectronCollection.push_back( GenLeptonCollection_A[0] );	//include e's from A and FSR at even/odd positions
			ElectronCollection.push_back( GenLeptonCollection_e[0] );
			ElectronCollection.push_back( GenLeptonCollection_A[1] );
			ElectronCollection.push_back( GenLeptonCollection_e[1] );
			
			vector<Double_t> MotherPtCollection;
			MotherPtCollection.push_back( MotherPtCollection_A[0] );
			MotherPtCollection.push_back( MotherPtCollection_e[0] );
			MotherPtCollection.push_back( MotherPtCollection_A[1] );
			MotherPtCollection.push_back( MotherPtCollection_e[1] );
			
			
			for(Int_t j=0; j<ElectronCollection.size();j++){			
				TLorentzVector vec_1;
				vec_1.SetPtEtaPhiM(ElectronCollection[j].Pt,ElectronCollection[j].eta,ElectronCollection[j].phi,ElectronCollection[j].Mass);
				for (Int_t k=j+1;k<4;k++){
					if ( MotherPtCollection[j] == MotherPtCollection[k] ) {
						TLorentzVector vec_2;
						vec_2.SetPtEtaPhiM(ElectronCollection[k].Pt,ElectronCollection[k].eta,ElectronCollection[k].phi,ElectronCollection[k].Mass);
						h_DeltaR->Fill(vec_1.DeltaR(vec_2));
					}			
				}		
			}
		}					       
		

		if( (Int_t)GenLeptonCollection_e.size() == 4 )
		{	
			Event_4e++;
						
			vector<Double_t> MotherPtCollection;
			MotherPtCollection.push_back( MotherPtCollection_e[0] );
			MotherPtCollection.push_back( MotherPtCollection_e[1] );
			MotherPtCollection.push_back( MotherPtCollection_e[2] );
			MotherPtCollection.push_back( MotherPtCollection_e[3] );
			
			vector<Int_t> mother_reconstructed;	//keep track of electrons with reconstructed mother (e.g. entry = 1 )
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			
			
			vector<GenLepton> ElectronCollection;
			ElectronCollection.push_back( GenLeptonCollection_e[0] );
			ElectronCollection.push_back( GenLeptonCollection_e[1] );
			ElectronCollection.push_back( GenLeptonCollection_e[2] );
			ElectronCollection.push_back( GenLeptonCollection_e[3] );
			
			
			// loop until all electrons are reconstructed from A particles
			while ( mother_reconstructed[0] != 1 && mother_reconstructed[1] != 1 && mother_reconstructed[2] != 1 && mother_reconstructed[3] != 1  )
			{
				for (Int_t j=0; j<NonFinalLeptonCollection.size(); j++)
				{
					for (Int_t k=0; k<ElectronCollection.size();k++)
						{
						if ( NonFinalLeptonCollection[j].Pt == MotherPtCollection[k] ){
							MotherPtCollection[k] = NonFinalLeptonCollection[j].Mother_Pt;
							if ( fabs(NonFinalLeptonCollection[j].Mother) == 36 ) {
								mother_reconstructed[k] = 1;
							}
						}
					}
				}
			}
			
			
			for(Int_t j=0; j<ElectronCollection.size();j++){			
				TLorentzVector vec_1;
				vec_1.SetPtEtaPhiM(ElectronCollection[j].Pt,ElectronCollection[j].eta,ElectronCollection[j].phi,ElectronCollection[j].Mass);
				for (Int_t k=j+1;k<4;k++){
					if ( MotherPtCollection[j] == MotherPtCollection[k] ) {
						TLorentzVector vec_2;
						vec_2.SetPtEtaPhiM(ElectronCollection[k].Pt,ElectronCollection[k].eta,ElectronCollection[k].phi,ElectronCollection[k].Mass);
						h_DeltaR->Fill(vec_1.DeltaR(vec_2));
					}				
				}		
			}
			
		}
	}
	//Double_t f_4A,f_4e,f_2A2e;
	//f_4A = Event_4A/(Event_4A + Event_4e + Event_2A2e);
	//f_4e = Event_4e/(Event_4A + Event_4e + Event_2A2e);
	//f_2A2e = Event_2A2e/(Event_4A + Event_4e + Event_2A2e);
	h_DeltaR->SetLineColor(kRed);
	//h_DeltaR_4A->SetLineColor(kBlue);
	//h_DeltaR_4e->SetLineColor(kGreen);
	//h_DeltaR_2A2e->SetStats(0);
	//h_DeltaR_4A->SetStats(0);
	//h_DeltaR->SetStats(0);
	h_DeltaR->GetXaxis()->SetTitle("DeltaR ");
	h_DeltaR->SetTitle("DeltaR distribution using mother_pT, A=50GeV");	
	h_DeltaR->DrawNormalized();
	//h_DeltaR_2A2e->DrawNormalized("same",f_2A2e);
	//h_DeltaR_4A->DrawNormalized("same",f_4A);
	//TLegend *legend = new TLegend(.75, .80, .95, .95);
	//legend->AddEntry(h_DeltaR_4e, "mother: 4e ");
	//legend->AddEntry(h_DeltaR_4A, "mother: 2A ");
	//legend->AddEntry(h_DeltaR_2A2e, " mother: A + 2e");
	//legend->Draw();
	
	//cout<<f_4A<<", "<<f_4e<<", "<<f_2A2e<<endl;
}
