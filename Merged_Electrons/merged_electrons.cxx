
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

void merged_electrons()
{
	TFile *f_output = TFile::Open("ROOTFile_merged_electrons.root", "RECREATE");
	// -- make chain -- //
	TChain *chain = new TChain("recoTree/DYTree");
	//chain->Add("/Users/nanamaan/Documents/ROOT/Files/Signal_Sample/ntuple_skim_2000_1.root");
	chain->Add("/scratch/kplee/DYntuple/80X/DYntuple_v20170728_GeneralTrack_HToAATo4e/H_2000GeV_A_1GeV/ntuple_*.root");	

	NtupleHandle *ntuple = new NtupleHandle( chain );
	ntuple->TurnOnBranches_GenLepton();
	ntuple->TurnOnBranches_Electron();


	TTree *tree_merged_electrons = chain->CloneTree(0);
		
	Int_t nEvent = chain->GetEntries();
	cout << "\t[Total Events: " << nEvent << "]" << endl;
	Int_t n=0;
	for(Int_t i=0; i<nEvent; i++)
	{
		loadBar(i+1, nEvent, 100, 100);
		
		ntuple->GetEvent(i);
		
		vector<Electron> merged_electrons;
		
		//collect all the relevant electrons from generator level
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
				//collect merged electrons with mother near gen_electron1 (necessary for next step)
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
		}	
	
		//cout<<(Int_t)merged_electrons.size()<<endl;
		
		//initialize a new tree with only merged electrons
		
		Int_t num_elec = (Int_t)merged_electrons.size();
		Int_t Nelectrons;
    	Double_t Electron_Energy[(Int_t)merged_electrons.size()];
    	Double_t Electron_pT[(Int_t)merged_electrons.size()];
    	Double_t Electron_eta[(Int_t)merged_electrons.size()];
    	Double_t Electron_phi[(Int_t)merged_electrons.size()];
    	Int_t Electron_charge[(Int_t)merged_electrons.size()];
    	Double_t Electron_gsfpT[(Int_t)merged_electrons.size()];
    	Double_t Electron_gsfPx[(Int_t)merged_electrons.size()];
    	Double_t Electron_gsfPy[(Int_t)merged_electrons.size()];
    	Double_t Electron_gsfPz[(Int_t)merged_electrons.size()];
    	Double_t Electron_gsfEta[(Int_t)merged_electrons.size()];
    	Double_t Electron_gsfPhi[(Int_t)merged_electrons.size()];
    	Double_t Electron_etaSC[(Int_t)merged_electrons.size()];
    	Double_t Electron_phiSC[(Int_t)merged_electrons.size()];
    	Double_t Electron_etaWidth[(Int_t)merged_electrons.size()];
    	Double_t Electron_phiWidth[(Int_t)merged_electrons.size()];
    	Double_t Electron_dEtaIn[(Int_t)merged_electrons.size()];
    	Double_t Electron_dPhiIn[(Int_t)merged_electrons.size()];
    	Double_t Electron_sigmaIEtaIEta[(Int_t)merged_electrons.size()];
    	Double_t Electron_Full5x5_SigmaIEtaIEta[(Int_t)merged_electrons.size()];
    	Double_t Electron_HoverE[(Int_t)merged_electrons.size()];
    	Double_t Electron_fbrem[(Int_t)merged_electrons.size()];
    	Double_t Electron_eOverP[(Int_t)merged_electrons.size()];
    	Double_t Electron_InvEminusInvP[(Int_t)merged_electrons.size()];
    	Double_t Electron_dxyVTX[(Int_t)merged_electrons.size()];
    	Double_t Electron_dzVTX[(Int_t)merged_electrons.size()];
    	Double_t Electron_dxy[(Int_t)merged_electrons.size()];
    	Double_t Electron_dz[(Int_t)merged_electrons.size()];
    	Double_t Electron_dxyBS[(Int_t)merged_electrons.size()];
    	Double_t Electron_dzBS[(Int_t)merged_electrons.size()];
    	Double_t Electron_chIso03[(Int_t)merged_electrons.size()];
		Double_t Electron_nhIso03[(Int_t)merged_electrons.size()];
   		Double_t Electron_phIso03[(Int_t)merged_electrons.size()];
   		Double_t Electron_ChIso03FromPU[(Int_t)merged_electrons.size()];
    	Int_t Electron_mHits[(Int_t)merged_electrons.size()];
    	Double_t Electron_EnergySC[(Int_t)merged_electrons.size()];
		Double_t Electron_preEnergySC[(Int_t)merged_electrons.size()];
   		Double_t Electron_rawEnergySC[(Int_t)merged_electrons.size()];
   		Double_t Electron_etSC[(Int_t)merged_electrons.size()];
   		Double_t Electron_E15[(Int_t)merged_electrons.size()];
   		Double_t Electron_E25[(Int_t)merged_electrons.size()];
   		Double_t Electron_E55[(Int_t)merged_electrons.size()];
    	Double_t Electron_RelPFIso_dBeta[(Int_t)merged_electrons.size()];
    	Double_t Electron_RelPFIso_Rho[(Int_t)merged_electrons.size()];
    	Double_t Electron_r9[(Int_t)merged_electrons.size()];
    	Int_t Electron_ecalDriven[(Int_t)merged_electrons.size()];
    	Bool_t Electron_passMediumID[(Int_t)merged_electrons.size()];
	
	
		for ( Int_t j=0; j<(Int_t)merged_electrons.size();j++) {
			Electron_Energy[j]= merged_electrons[j].Energy;
	   		Electron_pT[j] = merged_electrons[j].Pt;
			Electron_eta[j]= merged_electrons[j].eta;
			Electron_phi[j]= merged_electrons[j].phi;
			Electron_charge[j]= merged_electrons[j].charge;
			Electron_gsfpT[j]= merged_electrons[j].gsfpT;
			Electron_gsfPx[j]= merged_electrons[j].gsfPx;
			Electron_gsfPy[j]= merged_electrons[j].gsfPy;
			Electron_gsfPz[j]= merged_electrons[j].gsfPz;
			Electron_gsfEta[j]= merged_electrons[j].gsfEta;
			Electron_gsfPhi[j]= merged_electrons[j].gsfPhi;
			Electron_etaSC[j]= merged_electrons[j].etaSC;
			Electron_phiSC[j]= merged_electrons[j].phiSC;
			Electron_etaWidth[j]= merged_electrons[j].etaWidth;
			Electron_phiWidth[j]= merged_electrons[j].phiWidth;
			Electron_dEtaIn[j]= merged_electrons[j].dEtaIn;
			Electron_dPhiIn[j]= merged_electrons[j].dPhiIn;
			Electron_sigmaIEtaIEta[j]= merged_electrons[j].sigmaIEtaIEta;
			Electron_Full5x5_SigmaIEtaIEta[j]= merged_electrons[j].Full5x5_SigmaIEtaIEta;
			Electron_HoverE[j]= merged_electrons[j].HoverE;
			Electron_fbrem[j]= merged_electrons[j].fbrem;
			Electron_eOverP[j]= merged_electrons[j].eOverP;
			Electron_InvEminusInvP[j]= merged_electrons[j].InvEminusInvP;
			Electron_dxyVTX[j]= merged_electrons[j].dxyVTX;
			Electron_dzVTX[j]= merged_electrons[j].dzVTX;
			Electron_dxy[j]= merged_electrons[j].dxy;
			Electron_dz[j]= merged_electrons[j].dz;
			Electron_dxyBS[j]= merged_electrons[j].dxyBS;
			Electron_dzBS[j]= merged_electrons[j].dzBS;
			Electron_chIso03[j]= merged_electrons[j].chIso03;
			Electron_nhIso03[j]= merged_electrons[j].nhIso03;
			Electron_phIso03[j]= merged_electrons[j].phIso03;
			Electron_ChIso03FromPU[j]= merged_electrons[j].ChIso03FromPU;
			Electron_mHits[j]= merged_electrons[j].mHits;
			Electron_EnergySC[j]= merged_electrons[j].EnergySC;
	  		Electron_preEnergySC[j]= merged_electrons[j].preEnergySC;
	   		Electron_rawEnergySC[j]= merged_electrons[j].rawEnergySC;
	   		Electron_etSC[j]= merged_electrons[j].etSC;
	   		Electron_E15[j]= merged_electrons[j].E15;
	   		Electron_E25[j]= merged_electrons[j].E25;
	   		Electron_E55[j]= merged_electrons[j].E55;
	   		Electron_RelPFIso_dBeta[j]= merged_electrons[j].RelPFIso_dBeta;
	   		Electron_RelPFIso_Rho[j]= merged_electrons[j].RelPFIso_Rho;
	   		Electron_r9[j]= merged_electrons[j].r9;
	   		Electron_ecalDriven[j]= merged_electrons[j].ecalDriven;
			Electron_passMediumID[j]= merged_electrons[j].passMediumID;
		
		
			//cout<<merged_electrons[j].Energy<<endl;
		}
			    	
	    
	    
	    tree_merged_electrons->SetBranchAddress("Nelectrons", &num_elec);
	   	tree_merged_electrons->SetBranchAddress("Electron_Energy", &Electron_Energy);
	   	tree_merged_electrons->SetBranchAddress("Electron_pT", &Electron_pT);
	   	tree_merged_electrons->SetBranchAddress("Electron_eta", &Electron_eta);
		tree_merged_electrons->SetBranchAddress("Electron_phi", &Electron_phi);
	    tree_merged_electrons->SetBranchAddress("Electron_charge", &Electron_charge);
	    tree_merged_electrons->SetBranchAddress("Electron_gsfpT", &Electron_gsfpT);
	    tree_merged_electrons->SetBranchAddress("Electron_gsfPx", &Electron_gsfPx);
	    tree_merged_electrons->SetBranchAddress("Electron_gsfPy", &Electron_gsfPy);
	    tree_merged_electrons->SetBranchAddress("Electron_gsfPz", &Electron_gsfPz);
	    tree_merged_electrons->SetBranchAddress("Electron_gsfEta", &Electron_gsfEta);
	    tree_merged_electrons->SetBranchAddress("Electron_gsfPhi", &Electron_gsfPhi);
	    tree_merged_electrons->SetBranchAddress("Electron_etaSC", &Electron_etaSC);
	    tree_merged_electrons->SetBranchAddress("Electron_phiSC", &Electron_phiSC);
	    tree_merged_electrons->SetBranchAddress("Electron_etaWidth", &Electron_etaWidth);
	    tree_merged_electrons->SetBranchAddress("Electron_phiWidth", &Electron_phiWidth);
	    tree_merged_electrons->SetBranchAddress("Electron_dEtaIn", &Electron_dEtaIn);
	    tree_merged_electrons->SetBranchAddress("Electron_dPhiIn", &Electron_dPhiIn);
		tree_merged_electrons->SetBranchAddress("Electron_sigmaIEtaIEta", &Electron_sigmaIEtaIEta);
	   	tree_merged_electrons->SetBranchAddress("Electron_Full5x5_SigmaIEtaIEta", &Electron_Full5x5_SigmaIEtaIEta);
    	tree_merged_electrons->SetBranchAddress("Electron_HoverE", &Electron_HoverE);
    	tree_merged_electrons->SetBranchAddress("Electron_fbrem", &Electron_fbrem);
		tree_merged_electrons->SetBranchAddress("Electron_eOverP", &Electron_eOverP);
	   	tree_merged_electrons->SetBranchAddress("Electron_InvEminusInvP", &Electron_InvEminusInvP);
	   	tree_merged_electrons->SetBranchAddress("Electron_dxyVTX", &Electron_dxyVTX);
	   	tree_merged_electrons->SetBranchAddress("Electron_dzVTX", &Electron_dzVTX);
	   	tree_merged_electrons->SetBranchAddress("Electron_dxy", &Electron_dxy);
	   	tree_merged_electrons->SetBranchAddress("Electron_dz", &Electron_dz);
	   	tree_merged_electrons->SetBranchAddress("Electron_dxyBS", &Electron_dxyBS);
	   	tree_merged_electrons->SetBranchAddress("Electron_dzBS", &Electron_dzBS);
	   	tree_merged_electrons->SetBranchAddress("Electron_chIso03", &Electron_chIso03);
	   	tree_merged_electrons->SetBranchAddress("Electron_nhIso03", &Electron_nhIso03);
	   	tree_merged_electrons->SetBranchAddress("Electron_phIso03", &Electron_phIso03);
	   	tree_merged_electrons->SetBranchAddress("Electron_ChIso03FromPU", &Electron_ChIso03FromPU);
		tree_merged_electrons->SetBranchAddress("Electron_mHits", &Electron_mHits);
   		tree_merged_electrons->SetBranchAddress("Electron_EnergySC", &Electron_EnergySC);
   		tree_merged_electrons->SetBranchAddress("Electron_preEnergySC", &Electron_preEnergySC);
	   	tree_merged_electrons->SetBranchAddress("Electron_rawEnergySC", &Electron_rawEnergySC);
	   	tree_merged_electrons->SetBranchAddress("Electron_etSC", &Electron_etSC);
	   	tree_merged_electrons->SetBranchAddress("Electron_E15", &Electron_E15);
	   	tree_merged_electrons->SetBranchAddress("Electron_E25", &Electron_E25);
		tree_merged_electrons->SetBranchAddress("Electron_E55", &Electron_E55);
	   	tree_merged_electrons->SetBranchAddress("Electron_RelPFIso_dBeta", &Electron_RelPFIso_dBeta);
	   	tree_merged_electrons->SetBranchAddress("Electron_RelPFIso_Rho", &Electron_RelPFIso_Rho);
	   	tree_merged_electrons->SetBranchAddress("Electron_r9", &Electron_r9);
	   	tree_merged_electrons->SetBranchAddress("Electron_ecalDriven", &Electron_ecalDriven);
	    tree_merged_electrons->SetBranchAddress("Electron_passMediumID", &Electron_passMediumID);
	    
	    tree_merged_electrons->Fill();
	}
		tree_merged_electrons->Write();
}
