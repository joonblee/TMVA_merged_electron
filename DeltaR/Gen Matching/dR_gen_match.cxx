
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

void dR_gen_match()
{
	TFile *f_output = TFile::Open("ROOTFile_DeltaR_gen_match_dEt.root", "RECREATE");
	TH1D* h_DeltaR_match = new TH1D( "h_DeltaR_match", "", 350, -0.03, 0.5 );
	TH1D* h_dEt = new TH1D( "h_dEt", "", 150, -0.03, 2.5 );
	// -- make chain -- //
	TChain *chain = new TChain("recoTree/DYTree");
	//chain->Add("/Users/nanamaan/Documents/ROOT/Files/Signal_Sample/ntuple_skim_2000_50.root");
	chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_HToAATo4e/H_2000GeV_A_50GeV/ntuple_*.root");	//import signal sample

	NtupleHandle *ntuple = new NtupleHandle( chain );
	ntuple->TurnOnBranches_GenLepton();
	ntuple->TurnOnBranches_Electron();

	Int_t nEvent = chain->GetEntries();
	cout << "\t[Total Events: " << nEvent << "]" << endl;
	
	
	for(Int_t i=0; i<nEvent; i++)
	{
		loadBar(i+1, nEvent, 100, 100);
		
		ntuple->GetEvent(i);
		
		
		//collect all the relevant electrons from generator level
		vector<GenLepton> NonFinalLeptonCollection;
		vector<GenLepton> GenLeptonCollection_A;	//collect electrons with mother A
		vector<Double_t> MotherPtCollection_A;
		vector<GenLepton> GenLeptonCollection_e; //collect electrons from FSR
		vector<Double_t> MotherPtCollection_e;
		vector<GenLepton> ElectronCollection;
		vector<Double_t> MotherPtCollection;	
		Int_t NGenLeptons = ntuple->gnpair;
		
		for(Int_t i_gen=0; i_gen<NGenLeptons; i_gen++)
		{
			GenLepton genlep;
			genlep.FillFromNtuple(ntuple, i_gen);
			// collect gen-e with mother dark photon
			if( fabs(genlep.Mother) == 36 && genlep.fromHardProcessFinalState == 1 && fabs(genlep.ID) == 11) {
				GenLepton elec;
				elec.FillFromNtuple(ntuple,i_gen);
				GenLeptonCollection_A.push_back( elec );
				MotherPtCollection_A.push_back( elec.Mother_Pt );
			}
			//collect gen-e from FSR
			if( fabs(genlep.Mother) == 11 && genlep.fromHardProcessFinalState == 1 && fabs(genlep.ID) ==11) {
				GenLepton elec;
				elec.FillFromNtuple(ntuple,i_gen);
				GenLeptonCollection_e.push_back( elec );
				MotherPtCollection_e.push_back( elec.Mother_Pt );
			}
			//collect gen-e that are possible mothers of FSR electrons
			if( fabs(genlep.ID) == 11 && genlep.isPromptFinalState == 0 ){ // -- electron (ID==11) & not final state  -- //
				NonFinalLeptonCollection.push_back( genlep );
			}	
		}
		
	
		//compute mother_pt (A) for all electrons in the generator level:
		//1) case 4e from 2A:
		
		if( (Int_t)GenLeptonCollection_A.size() ==4 )
		{	
			
			ElectronCollection.push_back( GenLeptonCollection_A[0] );
			ElectronCollection.push_back( GenLeptonCollection_A[1] );
			ElectronCollection.push_back( GenLeptonCollection_A[2] );
			ElectronCollection.push_back( GenLeptonCollection_A[3] );
			
			MotherPtCollection.push_back( MotherPtCollection_A[0] );
			MotherPtCollection.push_back( MotherPtCollection_A[1] );
			MotherPtCollection.push_back( MotherPtCollection_A[2] );
			MotherPtCollection.push_back( MotherPtCollection_A[3] );
			
		}
		
		
		// case of: 2e from A + 2e from FSR / reconstruct mother_pt 
		if( (Int_t)GenLeptonCollection_A.size() ==2 && (Int_t)GenLeptonCollection_e.size() == 2 )
		{	
			vector<Int_t> mother_reconstructed;	//keep track of electrons with reconstructed mother (e.g. entry = 1 )
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			
			// loop until all electrons are reconstructed from A (dark photon) particles
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
			

			ElectronCollection.push_back( GenLeptonCollection_A[0] );	//include e's from A and FSR at even/odd positions
			ElectronCollection.push_back( GenLeptonCollection_e[0] );
			ElectronCollection.push_back( GenLeptonCollection_A[1] );
			ElectronCollection.push_back( GenLeptonCollection_e[1] );
			
			MotherPtCollection.push_back( MotherPtCollection_A[0] );
			MotherPtCollection.push_back( MotherPtCollection_e[0] );
			MotherPtCollection.push_back( MotherPtCollection_A[1] );
			MotherPtCollection.push_back( MotherPtCollection_e[1] );
		}
		
		//case of 4e from FSR:
		
		if( (Int_t)GenLeptonCollection_e.size() == 4 )
		{			
			MotherPtCollection.push_back( MotherPtCollection_e[0] );
			MotherPtCollection.push_back( MotherPtCollection_e[1] );
			MotherPtCollection.push_back( MotherPtCollection_e[2] );
			MotherPtCollection.push_back( MotherPtCollection_e[3] );
			
			vector<Int_t> mother_reconstructed;	//keep track of electrons with reconstructed mother (e.g. entry = 1 )
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			mother_reconstructed.push_back( 0 );
			
			
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
		}		
		
		
		//collect all electrons at reconstruction level
		vector<Electron> RecoElectronCollection;
		Int_t nElec = ntuple->Nelectrons;
		for (Int_t i_elec=0; i_elec<nElec; i_elec++) {
			Electron elec;
			elec.FillFromNtuple(ntuple, i_elec);
			RecoElectronCollection.push_back( elec );
			}	
	
		
		TLorentzVector vec_1,vec_2,vec_3,vec_4;
		vec_1.SetPtEtaPhiM(ElectronCollection[0].Pt,ElectronCollection[0].eta,ElectronCollection[0].phi,ElectronCollection[0].Mass);
		vec_2.SetPtEtaPhiM(ElectronCollection[1].Pt,ElectronCollection[1].eta,ElectronCollection[1].phi,ElectronCollection[1].Mass);
		vec_3.SetPtEtaPhiM(ElectronCollection[2].Pt,ElectronCollection[2].eta,ElectronCollection[2].phi,ElectronCollection[2].Mass);
		vec_4.SetPtEtaPhiM(ElectronCollection[3].Pt,ElectronCollection[3].eta,ElectronCollection[3].phi,ElectronCollection[3].Mass);
		
		
		
			
		// fill in Electron_matched vectors will all electrons falling into DeltaR<0.3 cone
		vector<Electron> Electron_matched_1;
		vector<Electron> Electron_matched_2;
		
		for (Int_t j=0; j<ElectronCollection.size();j++){
			for (Int_t k=0; k<RecoElectronCollection.size();k++){
				TLorentzVector vec_reco_2;
				vec_reco_2.SetPtEtaPhiE(RecoElectronCollection[k].Pt,RecoElectronCollection[k].eta,RecoElectronCollection[k].phi,RecoElectronCollection[k].Energy);
				Double_t DeltaR;
				DeltaR = vec_1.DeltaR(vec_reco_2);
				if (DeltaR < 0.3) {
					// fill in Electron_matched_1 with electrons that are paired up with the gen-electron that is in first position in ElectronCollection
					if ( MotherPtCollection[j] == MotherPtCollection[0] ) {		
						Int_t n_copies=0;
						for (Int_t q=0; q<(Int_t)Electron_matched_1.size();q++){
							if ( RecoElectronCollection[k].Pt == Electron_matched_1[q].Pt ) {
								n_copies++;
							}
						}
						if (n_copies == 0 ) {
							Electron_matched_1.push_back( RecoElectronCollection[k] );
						}
					}
					else {
						Int_t n_copies=0;
						for (Int_t q=0; q<(Int_t)Electron_matched_2.size();q++){
							if ( RecoElectronCollection[k].Pt == Electron_matched_2[q].Pt ) {
								n_copies++;
							}
						}
						if (n_copies == 0 ) {
							Electron_matched_2.push_back( RecoElectronCollection[k] );
						}
					}																		
				}
			}
		}		
		
			
		//sort by pT, only compute DeltaR between highest pT tracks	
		sort( Electron_matched_1.begin(), Electron_matched_1.end(), CompareObject );	
		sort( Electron_matched_2.begin(), Electron_matched_2.end(), CompareObject );	
		
				
		if ( (Int_t)Electron_matched_1.size() >= 2  ) {
			TLorentzVector vec_reco_1,vec_reco_2;
			vec_reco_1.SetPtEtaPhiE(Electron_matched_1[0].Pt,Electron_matched_1[0].eta,Electron_matched_1[0].phi,Electron_matched_1[0].Energy);
			vec_reco_2.SetPtEtaPhiE(Electron_matched_1[1].Pt,Electron_matched_1[1].eta,Electron_matched_1[1].phi,Electron_matched_1[1].Energy);
			h_DeltaR_match->Fill ( vec_reco_1.DeltaR(vec_reco_2) );			
		}
		
		if ( (Int_t)Electron_matched_2.size() >= 2  ) {
			TLorentzVector vec_reco_1,vec_reco_2;
			vec_reco_1.SetPtEtaPhiE(Electron_matched_2[0].Pt,Electron_matched_2[0].eta,Electron_matched_2[0].phi,Electron_matched_2[0].Energy);
			vec_reco_2.SetPtEtaPhiE(Electron_matched_2[1].Pt,Electron_matched_2[1].eta,Electron_matched_2[1].phi,Electron_matched_2[1].Energy);
			h_DeltaR_match->Fill ( vec_reco_1.DeltaR(vec_reco_2) );
		}
		
		
		

		Bool_t Flag_IsLessThan4e = kFALSE;		
		if ((Int_t)ElectronCollection.size() == 4 && (Int_t)Electron_matched_1.size() + (Int_t)Electron_matched_2.size() < 4 ) Flag_IsLessThan4e = kTRUE;	
	
	
	// compute missing transverse energy ratio for events where not all electrons are reconstructed
		if ( Flag_IsLessThan4e ) {
			Double_t dEt,Et_gen1,Et_gen2;
			
			if (MotherPtCollection[1] == MotherPtCollection[0]) {
				for (Int_t j=0; j<Electron_matched_1.size();j++) {
					Et_gen1 = vec_1.Et();
					Et_gen2 = vec_2.Et();
					dEt = fabs( Electron_matched_1[j].etSC - Et_gen1 - Et_gen2 );
					h_dEt->Fill ( dEt/(Electron_matched_1[j].etSC) );
				}

				for (Int_t j=0; j<Electron_matched_2.size();j++) {
					Et_gen1 = vec_3.Et();
					Et_gen2 = vec_4.Et();
					dEt = fabs( Electron_matched_2[j].etSC - Et_gen1 - Et_gen2 );
					h_dEt->Fill ( dEt/(Electron_matched_2[j].etSC) );
				}
			}
			
			if (MotherPtCollection[2] == MotherPtCollection[0]) {
				for (Int_t j=0; j<Electron_matched_1.size();j++) {
					Et_gen1 = vec_1.Et();
					Et_gen2 = vec_3.Et();
					dEt = fabs( Electron_matched_1[j].etSC - Et_gen1 - Et_gen2 );
					h_dEt->Fill ( dEt/(Electron_matched_1[j].etSC) );
				}
				for (Int_t j=0; j<Electron_matched_2.size();j++) {
					Et_gen1 = vec_2.Et();
					Et_gen2 = vec_4.Et();
					dEt = fabs( Electron_matched_2[j].etSC - Et_gen1 - Et_gen2 );
					h_dEt->Fill ( dEt/(Electron_matched_2[j].etSC) );
				}
			}
			else {
				for (Int_t j=0; j<Electron_matched_1.size();j++) {
					Et_gen1 = vec_1.Et();
					Et_gen2 = vec_4.Et();
					dEt = fabs( Electron_matched_1[j].etSC - Et_gen1 - Et_gen2 );
					h_dEt->Fill ( dEt/(Electron_matched_1[j].etSC) );
				}
				for (Int_t j=0; j<Electron_matched_2.size();j++) {
					Et_gen1 = vec_2.Et();
					Et_gen2 = vec_3.Et();
					dEt = fabs( Electron_matched_2[j].etSC - Et_gen1 - Et_gen2 );
					h_dEt->Fill ( dEt/(Electron_matched_2[j].etSC) );
				}
			}			
		}
	}	
	TCanvas* canvas = new TCanvas("canvas","");
	canvas->Divide(2,1);	
	canvas->cd(1);	
	h_DeltaR_match->SetLineColor(kBlue);
	h_DeltaR_match->SetStats(0);
	h_DeltaR_match->GetXaxis()->SetTitle("DeltaR ");
	h_DeltaR_match->SetTitle("DeltaR distribution using gen-matching, H=2000GeV, A=50GeV");	
	h_DeltaR_match->DrawNormalized();
	TLegend *legend = new TLegend(.75, .80, .95, .95);
	legend->AddEntry(h_DeltaR_match, "matching function");
	legend->Draw();
	canvas->cd(2);
	h_dEt->SetTitle(" dEt/Et w/ dEt = abs (Et_reco - Et_gen1 - Et_gen2), A=50GeV");	
	h_dEt->DrawNormalized();
	
	
}
