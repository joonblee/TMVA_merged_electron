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

void ttrack_sig()
{
	TFile *f_output = TFile::Open("ROOTFile_ttrack_sig.root", "RECREATE");
	// -- make chain -- //
	TChain *chain = new TChain("recoTree/DYTree");
	chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_HToAATo4e/H_2000GeV_A_1GeV/ntuple_*.root");	
	
	TH1D *hist_pT_ttrack = new TH1D ( "hist_pT", "Delta pT distribution between gsf tracks and gen-electrons" , 100 , 0 , 1.3 );
	TH1D *hist_pT_mer = new TH1D ( "hist_pT_", " Delta pT distribution between gsf tracks and gen-electrons" , 100 , 0 , 1.3 );
	TH1D *hist_pT = new TH1D ( "hist_pT gsf tracks", " pT distribution of highest pT gsf tracks" , 200 , 0 , 700 );
	TH1D *hist_DeltaR_ttrack = new TH1D ( "hist_DeltaR_", "DeltaR distribution between gsf tracks and gen-electrons" , 100 , 0 , 0.1 );
	TH1D *hist_DeltaR_mer = new TH1D ( "hist_DeltaR", " DeltaR distribution between gsf tracks and gen-electrons" , 100 , 0 , 0.1 );
	TH1D *hist_tracks = new TH1D ( "hist_tracks", "gsf tracks within DeltaR<0.3 from merged-e" , 10 , 0 , 10 );

	NtupleHandle *ntuple = new NtupleHandle( chain );
	ntuple->TurnOnBranches_GenLepton();
	ntuple->TurnOnBranches_Electron();
	ntuple->TurnOnBranches_TTrack();
	
	Int_t nEvent = chain->GetEntries();
	cout << "\t[Total Events: " << nEvent << "]" << endl;
	Int_t n=0;
	for(Int_t i=0; i<nEvent; i++)
	{
		loadBar(i+1, nEvent, 100, 100);
		ntuple->GetEvent(i);
		
		//store merged electrons during each event
		vector<Electron> merged_electrons;
		
		//collect all the relevant electrons from generator level
		vector<GenLepton> NonFinalLeptonCollection;
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
		
		
		vector<GenLepton> ElectronCollection;
		vector<Double_t> MotherPtCollection;
		
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
			ElectronCollection.push_back( GenLeptonCollection_A[0] );	//include e's from A and FSR at even/odd positions
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
		
			
			//collect all particles for general track in this event
			vector<TTrack> TTrackElectronCollection;
			Int_t nTTrack = ntuple->NTT;
			for (Int_t i_ttrack=0; i_ttrack<nTTrack; i_ttrack++) {
				TTrack particle;
				particle.FillFromNtuple(ntuple, i_ttrack);
				TTrackElectronCollection.push_back( particle );
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
				Int_t n_gen_e=0;
				TLorentzVector vec_reco;
				vec_reco.SetPtEtaPhiE(RecoElectronCollection[k].Pt,RecoElectronCollection[k].eta,RecoElectronCollection[k].phi,RecoElectronCollection[k].Energy);
				DeltaR1 = vec_1.DeltaR(vec_reco);
				DeltaR2 = vec_2.DeltaR(vec_reco);
				DeltaR3 = vec_3.DeltaR(vec_reco);
				DeltaR4 = vec_4.DeltaR(vec_reco);
				//collect merged electrons with missing transverse energy ratio less than 10%
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
			
			
			for ( Int_t k=0; k<(Int_t)merged_electrons.size();k++) {
				Int_t gsf_tracks = 0;
				vector<Double_t> DeltaPt_merged;
				vector<Double_t> DeltaPt_ttrack;
				vector<Double_t> DeltaR_merged;
				vector<Double_t> DeltaR_ttrack;
				vector<TTrack> candidate_gsf_tracks;
				Double_t DeltaR_TTrack,eta_TTrack,phi_TTrack,DeltaR;
				//collect all gsf tracks within merged electrons
				for (Int_t j=0; j<(Int_t)TTrackElectronCollection.size();j++) {
					eta_TTrack = TTrackElectronCollection[j].eta;
					phi_TTrack = TTrackElectronCollection[j].phi;
					DeltaR_TTrack = TMath::Sqrt( (eta_TTrack - merged_electrons[k].gsfEta)*(eta_TTrack - merged_electrons[k].gsfEta) + TVector2::Phi_mpi_pi( phi_TTrack - merged_electrons[k].gsfPhi )*TVector2::Phi_mpi_pi( phi_TTrack - merged_electrons[k].gsfPhi ) );
					if ( DeltaR_TTrack < 0.3 ) {
						candidate_gsf_tracks.push_back( TTrackElectronCollection[j] );
						gsf_tracks++;
					}
				}
				
				if ( (Int_t)candidate_gsf_tracks.size() != 0 ) {
					sort( candidate_gsf_tracks.begin(), candidate_gsf_tracks.end(), CompareObject );										
					hist_pT->Fill( candidate_gsf_tracks[0].Pt );
					for (Int_t j=0; j<(Int_t)ElectronCollection.size();j++) {
						DeltaPt_ttrack.push_back( fabs((ElectronCollection[j].Pt - candidate_gsf_tracks[0].Pt)/ElectronCollection[j].Pt ));
						DeltaPt_merged.push_back( fabs((ElectronCollection[j].Pt - merged_electrons[k].gsfpT)/ElectronCollection[j].Pt ));
						DeltaR_ttrack.push_back( TMath::Sqrt( (candidate_gsf_tracks[0].eta - ElectronCollection[j].eta)*(candidate_gsf_tracks[0].eta - ElectronCollection[j].eta) + TVector2::Phi_mpi_pi( candidate_gsf_tracks[0].phi - ElectronCollection[j].phi )*TVector2::Phi_mpi_pi( candidate_gsf_tracks[0].phi - ElectronCollection[j].phi ) ) );
						DeltaR_merged.push_back( TMath::Sqrt( (merged_electrons[k].gsfEta - ElectronCollection[j].eta)*(merged_electrons[k].gsfEta - ElectronCollection[j].eta) + TVector2::Phi_mpi_pi( merged_electrons[k].gsfPhi - ElectronCollection[j].phi )*TVector2::Phi_mpi_pi( merged_electrons[k].gsfPhi - ElectronCollection[j].phi ) ) );
					}
					std::sort ( DeltaPt_ttrack.begin(),DeltaPt_ttrack.end() );
					std::sort ( DeltaPt_merged.begin(),DeltaPt_merged.end() );
					std::sort ( DeltaR_ttrack.begin(),DeltaR_ttrack.end() );
					std::sort ( DeltaR_merged.begin(),DeltaR_merged.end() );
					hist_pT_ttrack->Fill( DeltaPt_ttrack[0] );
					hist_pT_mer->Fill( DeltaPt_merged[0] );
					hist_DeltaR_ttrack->Fill ( DeltaR_ttrack[0] );
					hist_DeltaR_mer->Fill ( DeltaR_merged[0] );
					
				}	
				hist_tracks->Fill( gsf_tracks );					
			}
		}
	}	
	
	TCanvas *myC1 = new TCanvas("my super canvas 1");
	myC1->Divide(2,1);
	myC1->cd(1);
	hist_pT_mer->SetLineColor(kRed);
	hist_pT_ttrack->SetLineColor(kBlue);
	hist_pT_mer->SetStats(0);
	//hist_pT_mer->GetXaxis()->SetTitle(" pT (GeV) " );
	hist_pT_mer->DrawNormalized();	
	hist_pT_ttrack->DrawNormalized("same");
	TLegend *legend = new TLegend(.75, .80, .95, .95);
	legend->AddEntry(hist_pT_mer, "matched electrons");
	legend->AddEntry(hist_pT_ttrack, "gsf track w/ highest pT");
	legend->Draw();
	myC1->cd(2);
	hist_pT->GetXaxis()->SetTitle(" pT (GeV) " );
	hist_pT->SetTitle( "distribution of highest pT gsf tracks (not necessarily matched to gen-electron)");
	hist_pT->DrawNormalized();
	TCanvas *myC2 = new TCanvas("my super canvas 2");
	hist_DeltaR_mer->SetLineColor(kRed);
	hist_DeltaR_ttrack->SetLineColor(kBlue);
	hist_DeltaR_mer->SetStats(0);
	hist_DeltaR_mer->GetXaxis()->SetTitle(" DeltaR " );
	hist_DeltaR_mer->DrawNormalized();	
	hist_DeltaR_ttrack->DrawNormalized("same");
	TLegend *legend_ = new TLegend(.75, .80, .95, .95);
	legend_->AddEntry(hist_DeltaR_mer, "matched electrons");
	legend_->AddEntry(hist_DeltaR_ttrack, "gsf track w/ highest pT ");
	legend_->Draw();
	TCanvas *myC4 = new TCanvas("my super canvas 4");
	hist_tracks->SetTitle( "gsf tracks within DeltaR<0.3 from merged-e");
	hist_tracks->DrawNormalized();



}
