

//get rid of the sideband contribution, by subtracting the correlation
//of the SB, scaled down by the purity.
//To have a correct result, we need to work with the yield-normalized Ck
//Since for the plots we want different normalization, we make a dirty trick:
//We do the step above, but the output correlation function is multiplied 
//by the difference in the normalization of the yiel-normed to the range-normed correlation
void ForPaperProposal_pL(const double& Purity, const double& FracLeft, const int& SL, const int& SR){
	TString InputFolderData = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/NormTotYield/DataSignal/";
	TString InputFolderDataNorm = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/";
	TString InputFolderSB = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/NormTotYield/SideBands/";
	//the wider variation, as it has more statistics (and identical to the other variations)
	TString InputFileSL = TString::Format("CkSidebands09062020_SL%i_pL.root",SL);
	//the narrower variation, as here the different variations are not identical, so we stay close to the peak
	TString InputFileSR = TString::Format("CkSidebands09062020_SR%i_pL.root",SR);
	
	//DIMI PURITY
    TString NameSePurFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Purity_vs_kstar/Test/fOutputSE_0_64.root";
    TString NameMePurFile = "/home/dmihaylov/Dudek_Ubuntu/Work/Kclus/GeneralFemtoStuff/Using_CATS3/Output/pLambda_1/Purity_vs_kstar/Test/fOutputME_0_64.root";
    TString NamePurHist = "hPurity";
	
	TString OutputFolder = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/";
	//TString::Format("/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/Purity%.0f/",Purity*100.);
	const unsigned NumReb = 6;
	const unsigned NumVars = 45;
	
	for(unsigned uVar=0; uVar<NumVars; uVar++){
		TFile* OutputFile = new TFile(OutputFolder+TString::Format("CkSB_pL_L%.0f_SL%i_SR%i_P%.0f_%u.root",FracLeft*100,SL,SR,Purity*100,uVar),"recreate");
		//TDirectory** dOutput = new TDirectory* [NumReb];
		for(unsigned uReb=0; uReb<NumReb; uReb++){
			if(uReb==NumReb-1) uReb = 9;
			
			//original correlation
			TH1F* hCk;
			TH1F* hCk_Norm;
			//background left and right
			TH1F* hCkBL;
			TH1F* hCkBR;
			//background
			TH1F* hCkB;
			TH1F* hCkB_Norm;
			TH1F* hCkB_SillyNorm;
			//signal (corrected correlation)
			TH1F* hCkS;
			TH1F* hCkS_Norm;
			TH1F* hCkS_hCk_ratio;
			
			//purity
			TH1F* hPurSe;
			TH1F* hPurMe;
			TH1F* hPurRatioCkS;
			TH1F* hPurRatioCkB;
			//corrected using DIMI purity method
			TH1F* hCkS_Dimi;
			TH1F* hCkS_Dimi_Norm;
//UP TO HERE !!!!!!!!!!!
			//TH1F* h
			
			double Norm;

			//read the files and copy the histograms
			TFile* fInput;
			TDirectoryFile* dInput;
			TH1D* hInput;
			
			fInput = new TFile(InputFolderData+TString::Format("Ck_pL_%u.root",uVar));
			dInput = (TDirectoryFile*)(fInput->FindObjectAny(TString::Format("Binning_%.0f",float(uReb+1)*4.)));
			hInput = NULL;
			dInput->GetObject("hCkmult_SUM",hInput);
			OutputFile->cd();
			hCk = (TH1F*)hInput->Clone(TString::Format("hCk_%.0fMeV",float(uReb+1)*4.));
//printf("hInput=%p\n",hInput);
			delete fInput;
			//
			fInput = new TFile(InputFolderDataNorm+TString::Format("Ck_pL_%u.root",uVar));
			dInput = (TDirectoryFile*)(fInput->FindObjectAny(TString::Format("Binning_%.0f",float(uReb+1)*4.)));
			hInput = NULL;
			dInput->GetObject("hCkmult_SUM",hInput);
			OutputFile->cd();
			hCk_Norm = (TH1F*)hInput->Clone(TString::Format("hCk_Norm_%.0fMeV",float(uReb+1)*4.));
			delete fInput;
			//
			fInput = new TFile(InputFolderSB+InputFileSL);
			dInput = (TDirectoryFile*)(fInput->FindObjectAny(TString::Format("Binning_%.0f",float(uReb+1)*4.)));
			hInput = NULL;
			dInput->GetObject("hCkmult_SUM",hInput);
			OutputFile->cd();
			hCkBL = (TH1F*)hInput->Clone(TString::Format("hCkBL_%.0fMeV",float(uReb+1)*4.));
			delete fInput;
			//
			fInput = new TFile(InputFolderSB+InputFileSR);
			dInput = (TDirectoryFile*)(fInput->FindObjectAny(TString::Format("Binning_%.0f",float(uReb+1)*4.)));
			hInput = NULL;
			dInput->GetObject("hCkmult_SUM",hInput);
			OutputFile->cd();
			hCkBR = (TH1F*)hInput->Clone(TString::Format("hCkBR_%.0fMeV",float(uReb+1)*4.));
			delete fInput;
			//
			//should be the same in ALL bins, we just take one where we are sure we have some entries
			Norm = hCk_Norm->GetBinContent(5)/hCk->GetBinContent(5);
//printf("Norm = %f\n",Norm);
			//create the background histogram
			hCkB = (TH1F*)hCkBL->Clone(TString::Format("hCkB_%.0fMeV",float(uReb+1)*4.));
			hCkB->Scale(FracLeft);
			hCkB->Add(hCkBR,1.-FracLeft);
			
			//correct the correclation
			hCkS = (TH1F*)hCk->Clone(TString::Format("hCkS_%.0fMeV",float(uReb+1)*4.));
			//build the ratio to the original correlation
			hCkS_hCk_ratio	= (TH1F*)hCk->Clone(TString::Format("hCkS_hCk_ratio_%.0fMeV",float(uReb+1)*4.));
			for(unsigned uBin=0; uBin<hCkS->GetNbinsX(); uBin++){
				if(uBin>=hCkB->GetNbinsX()) break;
				float OriginalValue = hCk->GetBinContent(uBin+1);
				float OriginalError = hCk->GetBinError(uBin+1);
				float BackgroundValue = hCkB->GetBinContent(uBin+1);
				float BackgroundError = hCkB->GetBinError(uBin+1);
				float NewValue = OriginalValue-(1.-Purity)*BackgroundValue;
				float NewError = sqrt(OriginalError*OriginalError+(1.-Purity)*BackgroundError*(1.-Purity)*BackgroundError);
				hCkS->SetBinContent(uBin+1,NewValue);
				hCkS->SetBinError(uBin+1,NewError);	
				hCkS_hCk_ratio->SetBinContent(uBin+1,NewValue/Purity/OriginalValue);
				hCkS_hCk_ratio->SetBinError(uBin+1,NewValue/Purity/OriginalValue*(1.-Purity)*BackgroundError);							
			}
			
			//what we end up with is purity*Ck
			//we would like to have Ck alone, so we divide by the purity
			//in addition, we want to get back to the original normalization,
			//hence we multiply by that factor
			hCkS->Scale(1./Purity);
			
			//renomralize to the original normalization
			hCkS_Norm = (TH1F*)hCkS->Clone(TString::Format("hCkS_Norm_%.0fMeV",float(uReb+1)*4.));
			hCkS_Norm->Scale(Norm);
			hCkB_Norm = (TH1F*)hCkB->Clone(TString::Format("hCkB_Norm_%.0fMeV",float(uReb+1)*4.));
			hCkB_Norm->Scale(Norm);
			hCkB_SillyNorm = (TH1F*)hCkB->Clone(TString::Format("hCkB_SillyNorm_%.0fMeV",float(uReb+1)*4.));
			hCkB_SillyNorm->Scale(	hCkB_SillyNorm->Integral(hCkB_SillyNorm->FindBin(700),hCkB_SillyNorm->FindBin(1000))/
									(hCkB_SillyNorm->FindBin(1000)-hCkB_SillyNorm->FindBin(700)));
			
			//save the output
			OutputFile->cd();
			/*
			dOutput[uReb] = OutputFile->mkdir(TString::Format("Binning_%.0f",float(uReb+1)*4.));
			dOutput[uReb]->cd();
			dOutput[uReb]->Add(hCk);
			dOutput[uReb]->Add(hCkS);
			dOutput[uReb]->Add(hCkB);
			dOutput[uReb]->Add(hCk_Norm);
			dOutput[uReb]->Add(hCkS_Norm);
			dOutput[uReb]->Add(hCkB_Norm);
			OutputFile->cd();
			dOutput[uReb]->Write();
			*/
			hCk->Write();
			hCkS->Write();
			hCkB->Write();
			hCk_Norm->Write();
			hCkS_Norm->Write();
			hCkB_Norm->Write();
			hCkB_SillyNorm->Write();
			hCkS_hCk_ratio->Write();
		}//uReb
		//delete [] dOutput;
		delete OutputFile;
	}//uVar
}


//l4 r6 are the tightest
//l1 r1 were the original choice
void CreateCk_CorrectForSB(){
	/*
	for(float fracL=0; fracL<=1.0001; fracL+=0.1){
		for(float pur=0.94; pur<=0.97; pur+=0.01){
			ForPaperProposal_pL(pur,fracL,4,6);
		}
	}
	*/
	
	//ForPaperProposal_pL(0.96,0.0,4,6);
	//ForPaperProposal_pL(0.96,0.5,4,6);
	//ForPaperProposal_pL(0.96,0.6,4,6);
	//ForPaperProposal_pL(0.96,1.0,4,6);
	
	ForPaperProposal_pL(0.953,0.520,4,6);
	ForPaperProposal_pL(0.953,0.554,4,6);
	ForPaperProposal_pL(0.963,0.529,4,6);
	ForPaperProposal_pL(0.963,0.574,4,6);
	
}
