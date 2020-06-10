TH1F* add_CF(TH1F* hist_CF1, TH1F* hist_CF2,TString HistName)
{
  //Calculate CFs with error weighting
  TH1F* hist_CF_sum = (TH1F*)hist_CF1->Clone(HistName.Data());

  Int_t NBins = hist_CF_sum->GetNbinsX();

  for(Int_t i=0;i<NBins;i++)
  {
    double CF1_val = hist_CF1->GetBinContent(i+1);
    double CF1_err = hist_CF1->GetBinError(i+1);
    double CF2_val = hist_CF2->GetBinContent(i+1);
    double CF2_err = hist_CF2->GetBinError(i+1);

    //average for bin i:
    if(CF1_val != 0. && CF2_val != 0.)
    {
      double CF1_err_weight = 1./TMath::Power(CF1_err,2.);
      double CF2_err_weight = 1./TMath::Power(CF2_err,2.);

      double CF_sum_average = (CF1_err_weight*CF1_val + CF2_err_weight*CF2_val)/(CF1_err_weight+CF2_err_weight);
      double CF_sum_err = 1./TMath::Sqrt(CF1_err_weight+CF2_err_weight);

      hist_CF_sum->SetBinContent(i+1,CF_sum_average);
      hist_CF_sum->SetBinError(i+1,CF_sum_err);
    }
    else if(CF1_val == 0. && CF2_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF2_val);
      hist_CF_sum->SetBinError(i+1,CF2_err);
    }
    else if(CF2_val == 0 && CF1_val != 0.)
    {
      hist_CF_sum->SetBinContent(i+1,CF1_val);
      hist_CF_sum->SetBinError(i+1,CF1_err);
    }

  }
  return hist_CF_sum;
}


TDirectory* WorkHorse(unsigned uReb, TString Correlation, 
				//TString NameOutputFolder,TString NameOutputFile,
				TFile* OutputFile,
				TString NameInputFolder,TString NameInputFile,TString NameList){

	TString NamePP;
	TString NameAPAP;
	if(Correlation=="pp"){
		NamePP = "Particle0_Particle0";
		NameAPAP = "Particle1_Particle1";
	}
	else if(Correlation=="pL"){
		NamePP = "Particle0_Particle2";
		NameAPAP = "Particle1_Particle3";
	}	
	else if(Correlation=="pXim"){
		NamePP = "Particle0_Particle4";
		NameAPAP = "Particle1_Particle5";
	}		
	else if(Correlation=="LL"){
		NamePP = "Particle2_Particle2";
		NameAPAP = "Particle3_Particle3";
	}
	else{
		return NULL;
	}
	
	const double NormMin=200;
	const double NormMax=400;
	
	TH1F* hInSE_PP;
	TH1F* hInME_PP;
	TH1F* hInSE_APAP;
	TH1F* hInME_APAP;
	
	TH2F* hInSEmult_PP;
	TH2F* hInMEmult_PP;
	TH2F* hInSEmult_APAP;
	TH2F* hInMEmult_APAP;
	
	//[0] = raw
	//[1] = norm to yield and width
	//[2] = norm in range
	TH1F** hSE_PP = new TH1F* [3];
	TH1F** hME_PP = new TH1F* [3];
	TH1F** hSE_APAP = new TH1F* [3];
	TH1F** hME_APAP = new TH1F* [3];	
		
	TH1F** hCk_PP = new TH1F* [3];
	TH1F** hCk_APAP = new TH1F* [3];
	
	TH1F* hCk_WEIGHTED;
	
	TH1F* hSE_SUM;
	TH1F* hME_SUM;
	TH1F* hCk_SUM;
	

	TH2F* hSEmult_PP;
	TH2F* hMEmult_PP;
	TH2F* hSEmult_APAP;
	TH2F* hMEmult_APAP;
	TH2F* hSEmult_SUM;
	TH2F* hMEmult_SUM;
	TH1F* hCkmult_PP;
	TH1F* hCkmult_APAP;
	TH1F* hCkmult_SUM;
	
	//read file
    TFile* InputFile = new TFile(NameInputFolder+NameInputFile,"read");
    TDirectoryFile* dInput = (TDirectoryFile*)(InputFile->FindObjectAny(NameList));
    TList* lInput1 = NULL;
    dInput->GetObject(NameList,lInput1);
    TList* lInputPP = (TList*)lInput1->FindObject(NamePP);
    TList* lInputAPAP = (TList*)lInput1->FindObject(NameAPAP);
    hInSE_PP = (TH1F*)lInputPP->FindObject(TString::Format("SEDist_%s",NamePP.Data()));
    hInME_PP = (TH1F*)lInputPP->FindObject(TString::Format("MEDist_%s",NamePP.Data()));
    hInSE_APAP = (TH1F*)lInputAPAP->FindObject(TString::Format("SEDist_%s",NameAPAP.Data()));
    hInME_APAP = (TH1F*)lInputAPAP->FindObject(TString::Format("MEDist_%s",NameAPAP.Data())); 
    hInSEmult_PP = (TH2F*)lInputPP->FindObject(TString::Format("SEMultDist_%s",NamePP.Data()));
    hInMEmult_PP = (TH2F*)lInputPP->FindObject(TString::Format("MEMultDist_%s",NamePP.Data()));
    hInSEmult_APAP = (TH2F*)lInputAPAP->FindObject(TString::Format("SEMultDist_%s",NameAPAP.Data()));
    hInMEmult_APAP = (TH2F*)lInputAPAP->FindObject(TString::Format("MEMultDist_%s",NameAPAP.Data()));
    
    hInSE_PP->Sumw2();
    hInSE_PP->Rebin(uReb);
    hInSE_PP->Scale(1./double(uReb));
    
    hInME_PP->Sumw2();
    hInME_PP->Rebin(uReb);
    hInME_PP->Scale(1./double(uReb));    
    
    hInSE_APAP->Sumw2();
    hInSE_APAP->Rebin(uReb);
    hInSE_APAP->Scale(1./double(uReb));
    
    hInME_APAP->Sumw2();
    hInME_APAP->Rebin(uReb);
    hInME_APAP->Scale(1./double(uReb));   
    
    hInSEmult_PP->Sumw2();
    hInSEmult_PP->Rebin2D(uReb,1);
    hInSEmult_PP->Scale(1./double(uReb));
    
    hInMEmult_PP->Sumw2();
    hInMEmult_PP->Rebin2D(uReb,1);
    hInMEmult_PP->Scale(1./double(uReb));
    
    hInSEmult_APAP->Sumw2();
    hInSEmult_APAP->Rebin2D(uReb,1);
    hInSEmult_APAP->Scale(1./double(uReb));
    
    hInMEmult_APAP->Sumw2();
    hInMEmult_APAP->Rebin2D(uReb,1);
    hInMEmult_APAP->Scale(1./double(uReb));        

	//copy
	hSE_PP[0] = new TH1F(	"hSE_PP_Raw","hSE_PP_Raw",
							hInSE_PP->GetNbinsX(),
							hInSE_PP->GetBinLowEdge(1)*1000.,
							hInSE_PP->GetXaxis()->GetBinUpEdge(hInSE_PP->GetNbinsX())*1000.);
	for(unsigned uBin=0; uBin<hInSE_PP->GetNbinsX(); uBin++)
		hSE_PP[0]->SetBinContent(uBin+1,hInSE_PP->GetBinContent(uBin+1));
		
	hME_PP[0] = new TH1F(	"hME_PP_Raw","hME_PP_Raw",
							hInME_PP->GetNbinsX(),
							hInME_PP->GetBinLowEdge(1)*1000.,
							hInME_PP->GetXaxis()->GetBinUpEdge(hInME_PP->GetNbinsX())*1000.);
	for(unsigned uBin=0; uBin<hInME_PP->GetNbinsX(); uBin++)
		hME_PP[0]->SetBinContent(uBin+1,hInME_PP->GetBinContent(uBin+1));	
				
	hSE_PP[1] = (TH1F*)hSE_PP[0]->Clone("hSE_PP_Yield");
	hSE_PP[2] = (TH1F*)hSE_PP[0]->Clone("hSE_PP_Norm");
	
	hME_PP[1] = (TH1F*)hME_PP[0]->Clone("hME_PP_Yield");
	hME_PP[2] = (TH1F*)hME_PP[0]->Clone("hME_PP_Norm");
	
	hSE_PP[0]->Sumw2();
	hSE_PP[1]->Sumw2();
	hSE_PP[2]->Sumw2();
	hME_PP[0]->Sumw2();
	hME_PP[1]->Sumw2();
	hME_PP[2]->Sumw2();
	
	
	hSE_APAP[0] = new TH1F(	"hSE_APAP_Raw","hSE_APAP_Raw",
							hInSE_APAP->GetNbinsX(),
							hInSE_APAP->GetBinLowEdge(1)*1000.,
							hInSE_APAP->GetXaxis()->GetBinUpEdge(hInSE_APAP->GetNbinsX())*1000.);
	for(unsigned uBin=0; uBin<hInSE_APAP->GetNbinsX(); uBin++)
		hSE_APAP[0]->SetBinContent(uBin+1,hInSE_APAP->GetBinContent(uBin+1));
		
	hME_APAP[0] = new TH1F(	"hME_APAP_Raw","hME_APAP_Raw",
							hInME_APAP->GetNbinsX(),
							hInME_APAP->GetBinLowEdge(1)*1000.,
							hInME_APAP->GetXaxis()->GetBinUpEdge(hInME_APAP->GetNbinsX())*1000.);
	for(unsigned uBin=0; uBin<hInME_APAP->GetNbinsX(); uBin++)
		hME_APAP[0]->SetBinContent(uBin+1,hInME_APAP->GetBinContent(uBin+1));	
				
	hSE_APAP[1] = (TH1F*)hSE_APAP[0]->Clone("hSE_APAP_Yield");
	hSE_APAP[2] = (TH1F*)hSE_APAP[0]->Clone("hSE_APAP_Norm");
	
	hME_APAP[1] = (TH1F*)hME_APAP[0]->Clone("hME_APAP_Yield");
	hME_APAP[2] = (TH1F*)hME_APAP[0]->Clone("hME_APAP_Norm");
	
	hSE_APAP[0]->Sumw2();
	hSE_APAP[1]->Sumw2();
	hSE_APAP[2]->Sumw2();
	hME_APAP[0]->Sumw2();
	hME_APAP[1]->Sumw2();
	hME_APAP[2]->Sumw2();
	
	
	//normalizations
	hSE_PP[1]->Scale(1./hSE_PP[2]->Integral(),"width");
	hME_PP[1]->Scale(1./hME_PP[2]->Integral(),"width");
	
	hME_PP[2]->Scale(	hSE_PP[2]->Integral(hSE_PP[2]->FindBin(NormMin),hSE_PP[2]->FindBin(NormMax))/
						hME_PP[2]->Integral(hME_PP[2]->FindBin(NormMin),hME_PP[2]->FindBin(NormMax)));
	
	hSE_APAP[1]->Scale(1./hSE_APAP[2]->Integral(),"width");
	hME_APAP[1]->Scale(1./hME_APAP[2]->Integral(),"width");
	
	hME_APAP[2]->Scale(	hSE_APAP[2]->Integral(hSE_APAP[2]->FindBin(NormMin),hSE_APAP[2]->FindBin(NormMax))/
						hME_APAP[2]->Integral(hME_APAP[2]->FindBin(NormMin),hME_APAP[2]->FindBin(NormMax)));
	
	
	//PP and APAP correlations
	hCk_PP[0] = (TH1F*)hSE_PP[0]->Clone("hCk_PP_Raw");
	hCk_PP[0]->Divide(hME_PP[0]);
	
	hCk_PP[1] = (TH1F*)hSE_PP[1]->Clone("hCk_PP_Yield");
	hCk_PP[1]->Divide(hME_PP[1]);	
	
	hCk_PP[2] = (TH1F*)hSE_PP[2]->Clone("hCk_PP_Norm");
	hCk_PP[2]->Divide(hME_PP[2]);	
	
	hCk_APAP[0] = (TH1F*)hSE_APAP[0]->Clone("hCk_APAP_Raw");
	hCk_APAP[0]->Divide(hME_APAP[0]);

	hCk_APAP[1] = (TH1F*)hSE_APAP[1]->Clone("hCk_APAP_Yield");
	hCk_APAP[1]->Divide(hME_APAP[1]);	
	
	hCk_APAP[2] = (TH1F*)hSE_APAP[2]->Clone("hCk_APAP_Norm");
	hCk_APAP[2]->Divide(hME_APAP[2]);	
	
	//compute the Ck base on the weighted mean
	hCk_WEIGHTED = new TH1F("hCk_WEIGHTED","hCk_WEIGHTED",
							hCk_PP[2]->GetNbinsX(),
							hCk_PP[2]->GetBinLowEdge(1),
							hCk_PP[2]->GetXaxis()->GetBinUpEdge(hCk_PP[2]->GetNbinsX()));
	for(unsigned uBin=0; uBin<hCk_WEIGHTED->GetNbinsX(); uBin++){
		double Val1 = hCk_PP[2]->GetBinContent(uBin+1);
		double Val2 = hCk_APAP[2]->GetBinContent(uBin+1);
		double Err1 = hCk_PP[2]->GetBinError(uBin+1);
		double Err2 = hCk_APAP[2]->GetBinError(uBin+1);
		if(Err1&&Err2){
			hCk_WEIGHTED->SetBinContent(uBin+1,(Val1/Err1/Err1+Val2/Err2/Err2)/(1./Err1/Err1+1./Err2/Err2));
			hCk_WEIGHTED->SetBinError(uBin+1,1./sqrt(1./Err1/Err1+1./Err2/Err2));		
		}
		else if(Err1){
			hCk_WEIGHTED->SetBinContent(uBin+1,Val1);
			hCk_WEIGHTED->SetBinError(uBin+1,Err1);				
		}
		else if(Err2){
			hCk_WEIGHTED->SetBinContent(uBin+1,Val2);
			hCk_WEIGHTED->SetBinError(uBin+1,Err2);	
		}
		else{
			hCk_WEIGHTED->SetBinContent(uBin+1,1);
			hCk_WEIGHTED->SetBinError(uBin+1,1000);	
		}
	}
	
	//compute the Ck from the direct sum of the SE and ME
	
	double W_PP = hSE_PP[0]->Integral()/(hSE_PP[0]->Integral()+hSE_APAP[0]->Integral());
	
	hSE_SUM = (TH1F*)hSE_PP[1]->Clone("hSE_SUM");
	hSE_SUM->Scale(W_PP);
	hSE_SUM->Add(hSE_APAP[1],1.-W_PP);
	
	hME_SUM = (TH1F*)hME_PP[1]->Clone("hME_SUM");
	hME_SUM->Scale(W_PP);
	hME_SUM->Add(hME_APAP[1],1.-W_PP);	

//
/*
	hSE_SUM = (TH1F*)hSE_PP[1]->Clone("hSE_SUM");
	hME_SUM = (TH1F*)hME_PP[1]->Clone("hME_SUM");
	for(unsigned uBin=0; uBin<hSE_SUM->GetNbinsX(); uBin++){
		double WEIGHT = hSE_PP[1]->GetBinContent(uBin+1);
		WEIGHT /= (hSE_PP[1]->GetBinContent(uBin+1)+hSE_APAP[1]->GetBinContent(uBin+1));
		double Val1 = WEIGHT*hSE_PP[1]->GetBinContent(uBin+1);
		double Val2 = (1.-WEIGHT)*hSE_APAP[1]->GetBinContent(uBin+1);
		double Err1 = WEIGHT*hSE_PP[1]->GetBinError(uBin+1);
		double Err2 = (1.-WEIGHT)*hSE_APAP[1]->GetBinError(uBin+1);
		hSE_SUM->SetBinContent(uBin+1,Val1+Val2);
		hSE_SUM->SetBinError(uBin+1,sqrt(Err1*Err1+Err2*Err2));
	}
	for(unsigned uBin=0; uBin<hME_SUM->GetNbinsX(); uBin++){
		double WEIGHT = hME_PP[1]->GetBinContent(uBin+1);
		WEIGHT /= (hME_PP[1]->GetBinContent(uBin+1)+hME_APAP[1]->GetBinContent(uBin+1));
		double Val1 = WEIGHT*hME_PP[1]->GetBinContent(uBin+1);
		double Val2 = (1.-WEIGHT)*hME_APAP[1]->GetBinContent(uBin+1);
		double Err1 = WEIGHT*hME_PP[1]->GetBinError(uBin+1);
		double Err2 = (1.-WEIGHT)*hME_APAP[1]->GetBinError(uBin+1);
		hME_SUM->SetBinContent(uBin+1,Val1+Val2);
		hME_SUM->SetBinError(uBin+1,sqrt(Err1*Err1+Err2*Err2));
	}
*/
//
	//hSE_SUM = (TH1F*)hSE_PP[0]->Clone("hSE_SUM");
	//hSE_SUM->Add(hSE_APAP[0]);
	
	//hME_SUM = (TH1F*)hME_PP[0]->Clone("hME_SUM");
	//hME_SUM->Add(hME_APAP[0]);
	
	//hCk_SUM = (TH1F*)hSE_PP[1]->Clone("hCk_SUM");
	hCk_SUM = (TH1F*)hSE_SUM->Clone("hCk_SUM");
	hCk_SUM->Scale(hME_SUM->Integral()/hSE_SUM->Integral());
	hCk_SUM->Divide(hME_SUM);

unsigned Bin1 = hCk_SUM->FindBin(NormMin);
unsigned Bin2 = hCk_SUM->FindBin(NormMax);
//hCk_PP[0]->Scale(double(Bin2-Bin1+1)/hCk_PP[0]->Integral(Bin1,Bin2));
//hCk_PP[1]->Scale(double(Bin2-Bin1+1)/hCk_PP[1]->Integral(Bin1,Bin2));
//hCk_PP[2]->Scale(double(Bin2-Bin1+1)/hCk_PP[2]->Integral(Bin1,Bin2));
//hCk_APAP[0]->Scale(double(Bin2-Bin1+1)/hCk_APAP[0]->Integral(Bin1,Bin2));
//hCk_APAP[1]->Scale(double(Bin2-Bin1+1)/hCk_APAP[1]->Integral(Bin1,Bin2));
//hCk_APAP[2]->Scale(double(Bin2-Bin1+1)/hCk_APAP[2]->Integral(Bin1,Bin2));
hCk_WEIGHTED->Scale(double(Bin2-Bin1+1)/hCk_WEIGHTED->Integral(Bin1,Bin2));
hCk_SUM->Scale(double(Bin2-Bin1+1)/hCk_SUM->Integral(Bin1,Bin2));
TH1F* hADD =  add_CF(hCk_PP[2], hCk_APAP[2],"hADD");

/*
	//save the output
	TFile* OutputFile = new TFile(NameOutputFolder+NameOutputFile,"recreate");
	for(unsigned uSce=0; uSce<3; uSce++){
		hSE_PP[uSce]->Write();
		hME_PP[uSce]->Write();
		hCk_PP[uSce]->Write();
		hSE_APAP[uSce]->Write();
		hME_APAP[uSce]->Write();
		hCk_APAP[uSce]->Write();		
	}
	hCk_WEIGHTED->Write();
	hSE_SUM->Write();
	hME_SUM->Write();
	hCk_SUM->Write();
hADD->Write();
delete hADD;
*/
	

	// MULTIPLICITY BASE
	// idea: add all mult. Ck directly (SE and ME, which are normed to yield)
	
	//copy
	//hSEmult_PP = (TH2F*)hInSEmult_PP->Clone("hSEmult_PP");
	//hMEmult_PP = (TH2F*)hInMEmult_PP->Clone("hMEmult_PP");
	//hSEmult_APAP = (TH2F*)hInSEmult_APAP->Clone("hSEmult_APAP");
	//hMEmult_APAP = (TH2F*)hInMEmult_APAP->Clone("hMEmult_APAP");


	hSEmult_PP = new TH2F("hSEmult_PP","hSEmult_PP",
		hInSEmult_PP->GetXaxis()->GetNbins(),
		hInSEmult_PP->GetXaxis()->GetBinLowEdge(1)*1000.,
		hInSEmult_PP->GetXaxis()->GetBinUpEdge(hInSEmult_PP->GetXaxis()->GetNbins()*1000.),
		hInSEmult_PP->GetYaxis()->GetNbins(),
		hInSEmult_PP->GetYaxis()->GetBinLowEdge(1),
		hInSEmult_PP->GetYaxis()->GetBinUpEdge(hInSEmult_PP->GetYaxis()->GetNbins()));

	hMEmult_PP = new TH2F("hMEmult_PP","hMEmult_PP",
		hInMEmult_PP->GetXaxis()->GetNbins(),
		hInMEmult_PP->GetXaxis()->GetBinLowEdge(1)*1000.,
		hInMEmult_PP->GetXaxis()->GetBinUpEdge(hInMEmult_PP->GetXaxis()->GetNbins()*1000.),
		hInMEmult_PP->GetYaxis()->GetNbins(),
		hInMEmult_PP->GetYaxis()->GetBinLowEdge(1),
		hInMEmult_PP->GetYaxis()->GetBinUpEdge(hInMEmult_PP->GetYaxis()->GetNbins()));

	hSEmult_APAP = new TH2F("hSEmult_APAP","hSEmult_APAP",
		hInSEmult_APAP->GetXaxis()->GetNbins(),
		hInSEmult_APAP->GetXaxis()->GetBinLowEdge(1)*1000.,
		hInSEmult_APAP->GetXaxis()->GetBinUpEdge(hInSEmult_APAP->GetXaxis()->GetNbins()*1000.),
		hInSEmult_APAP->GetYaxis()->GetNbins(),
		hInSEmult_APAP->GetYaxis()->GetBinLowEdge(1),
		hInSEmult_APAP->GetYaxis()->GetBinUpEdge(hInSEmult_APAP->GetYaxis()->GetNbins()));

	hMEmult_APAP = new TH2F("hMEmult_APAP","hMEmult_APAP",
		hInMEmult_APAP->GetXaxis()->GetNbins(),
		hInMEmult_APAP->GetXaxis()->GetBinLowEdge(1)*1000.,
		hInMEmult_APAP->GetXaxis()->GetBinUpEdge(hInMEmult_APAP->GetXaxis()->GetNbins()*1000.),
		hInMEmult_APAP->GetYaxis()->GetNbins(),
		hInMEmult_APAP->GetYaxis()->GetBinLowEdge(1),
		hInMEmult_APAP->GetYaxis()->GetBinUpEdge(hInMEmult_APAP->GetYaxis()->GetNbins()));
	
	for(unsigned uMom=0; uMom<hSEmult_PP->GetXaxis()->GetNbins()+1; uMom++){
		for(unsigned uMult=0; uMult<hSEmult_PP->GetYaxis()->GetNbins()+1; uMult++){
			hSEmult_PP->SetBinContent(uMom+1,uMult+1,hInSEmult_PP->GetBinContent(uMom+1,uMult+1));
			hMEmult_PP->SetBinContent(uMom+1,uMult+1,hInMEmult_PP->GetBinContent(uMom+1,uMult+1));
			hSEmult_APAP->SetBinContent(uMom+1,uMult+1,hInSEmult_APAP->GetBinContent(uMom+1,uMult+1));
			hMEmult_APAP->SetBinContent(uMom+1,uMult+1,hInMEmult_APAP->GetBinContent(uMom+1,uMult+1));
		}
	}

	hSEmult_PP->Sumw2();
	hMEmult_PP->Sumw2();
	hSEmult_APAP->Sumw2();
	hMEmult_APAP->Sumw2();

	//norm
	//note, due to a bug, the last real bin in mult is the overflow bin
	unsigned NumMultBins = hSEmult_PP->GetYaxis()->GetNbins();
	double Integral;
	//dummy x-axis
	//this is the weight 
	TH1F* mult_weight_PP = new TH1F("mult_weight_PP","mult_weight_PP",NumMultBins+1,0,NumMultBins+1);
	TH1F* mult_weight_APAP = new TH1F("mult_weight_APAP","mult_weight_APAP",NumMultBins+1,0,NumMultBins+1);
	TH1F* mult_weight_SUM;// = new TH1F("mult_weight_SUM","mult_weight_SUM",NumMultBins+1,0,NumMultBins+1);
	double Weight_PP=0;
	double Weight_APAP=0;
	//norm hSEmult_PP and hSEmult_APAP
	for(unsigned uMult=0; uMult<NumMultBins+1; uMult++){
		Integral = hSEmult_PP->Integral(0,hSEmult_PP->GetXaxis()->GetNbins()+1,uMult+1,uMult+1);
		Weight_PP += Integral;
		mult_weight_PP->SetBinContent(uMult+1,Integral);
		for(unsigned uMom=0; uMom<hSEmult_PP->GetXaxis()->GetNbins()+1; uMom++){
			if(Integral){
				hSEmult_PP->SetBinContent(uMom+1,uMult+1,hSEmult_PP->GetBinContent(uMom+1,uMult+1)/Integral);
				hSEmult_PP->SetBinError(uMom+1,uMult+1,hSEmult_PP->GetBinError(uMom+1,uMult+1)/Integral);
			}
			else{
				hSEmult_PP->SetBinContent(uMom+1,uMult+1,0);
				hSEmult_PP->SetBinError(uMom+1,uMult+1,0);				
			}
		}
		
		Integral = hMEmult_PP->Integral(0,hMEmult_PP->GetXaxis()->GetNbins()+1,uMult+1,uMult+1);
		for(unsigned uMom=0; uMom<hMEmult_PP->GetXaxis()->GetNbins()+1; uMom++){
			if(Integral){
				hMEmult_PP->SetBinContent(uMom+1,uMult+1,hMEmult_PP->GetBinContent(uMom+1,uMult+1)/Integral);
				hMEmult_PP->SetBinError(uMom+1,uMult+1,hMEmult_PP->GetBinError(uMom+1,uMult+1)/Integral);
			}
			else{
				hMEmult_PP->SetBinContent(uMom+1,uMult+1,0);
				hMEmult_PP->SetBinError(uMom+1,uMult+1,0);			
			}
		}	
		
		Integral = hSEmult_APAP->Integral(0,hSEmult_APAP->GetXaxis()->GetNbins()+1,uMult+1,uMult+1);
		Weight_APAP += Integral;
		mult_weight_APAP->SetBinContent(uMult+1,Integral);
		for(unsigned uMom=0; uMom<hSEmult_APAP->GetXaxis()->GetNbins()+1; uMom++){
			if(Integral){
				hSEmult_APAP->SetBinContent(uMom+1,uMult+1,hSEmult_APAP->GetBinContent(uMom+1,uMult+1)/Integral);
				hSEmult_APAP->SetBinError(uMom+1,uMult+1,hSEmult_APAP->GetBinError(uMom+1,uMult+1)/Integral);
			}
			else{
				hSEmult_APAP->SetBinContent(uMom+1,uMult+1,0);
				hSEmult_APAP->SetBinError(uMom+1,uMult+1,0);
			}
		}
		
		Integral = hMEmult_APAP->Integral(0,hMEmult_APAP->GetXaxis()->GetNbins()+1,uMult+1,uMult+1);
		for(unsigned uMom=0; uMom<hMEmult_APAP->GetXaxis()->GetNbins()+1; uMom++){
			if(Integral){
				hMEmult_APAP->SetBinContent(uMom+1,uMult+1,hMEmult_APAP->GetBinContent(uMom+1,uMult+1)/Integral);
				hMEmult_APAP->SetBinError(uMom+1,uMult+1,hMEmult_APAP->GetBinError(uMom+1,uMult+1)/Integral);
			}
			else{
				hMEmult_APAP->SetBinContent(uMom+1,uMult+1,0);
				hMEmult_APAP->SetBinError(uMom+1,uMult+1,0);				
			}
		}
	}

	//first add the yields for PP and APAP, than normalize
	mult_weight_SUM = (TH1F*)mult_weight_PP->Clone("mult_weight_SUM");
	mult_weight_SUM->Add(mult_weight_APAP);
	mult_weight_SUM->Scale(1./mult_weight_SUM->Integral());
	
	mult_weight_PP->Scale(1./mult_weight_PP->Integral());
	mult_weight_APAP->Scale(1./mult_weight_APAP->Integral());
	
	Weight_PP = Weight_PP/(Weight_PP+Weight_APAP);
	Weight_APAP = 1.-Weight_PP;
	
	//sum up PP and APAP (they are normed here)
	hSEmult_SUM = (TH2F*)hSEmult_PP->Clone("hSEmult_SUM");
	hSEmult_SUM->Scale(Weight_PP);
	hSEmult_SUM->Add(hSEmult_APAP,Weight_APAP);
	
	hMEmult_SUM = (TH2F*)hMEmult_PP->Clone("hMEmult_SUM");
	hMEmult_SUM->Scale(Weight_PP);
	hMEmult_SUM->Add(hMEmult_APAP,Weight_APAP);

	//now we compute all correlation functions
	TH1F* hTemp_SE;
	TH1F* hTemp_ME;
	TH1D* hProjection;
	
	//
	hProjection = hSEmult_PP->ProjectionX("hProjection",1,1);
	hTemp_SE = (TH1F*)hProjection->Clone("hTemp_SE");
	hTemp_SE->Scale(mult_weight_PP->GetBinContent(1));
	delete hProjection;
	for(unsigned uMult=1; uMult<NumMultBins+1; uMult++){
		hProjection = hSEmult_PP->ProjectionX("hProjection",uMult+1,uMult+1);
		hTemp_SE->Add(hProjection,mult_weight_PP->GetBinContent(uMult+1));
		delete hProjection;
	}
			
	hProjection = hMEmult_PP->ProjectionX("hProjection",1,1);
	hTemp_ME = (TH1F*)hProjection->Clone("hTemp_ME");
	hTemp_ME->Scale(mult_weight_PP->GetBinContent(1));
	delete hProjection;
	for(unsigned uMult=1; uMult<NumMultBins+1; uMult++){
		hProjection = hMEmult_PP->ProjectionX("hProjection",uMult+1,uMult+1);
		hTemp_ME->Add(hProjection,mult_weight_PP->GetBinContent(uMult+1));
		delete hProjection;
	}
	
	hCkmult_PP = (TH1F*)hTemp_SE->Clone("hCkmult_PP");
	hCkmult_PP->Divide(hTemp_ME);
	delete hTemp_SE;
	delete hTemp_ME;

	//
	hProjection = hSEmult_APAP->ProjectionX("hProjection",1,1);
	hTemp_SE = (TH1F*)hProjection->Clone("hTemp_SE");
	hTemp_SE->Scale(mult_weight_APAP->GetBinContent(1));
	delete hProjection;
	for(unsigned uMult=1; uMult<NumMultBins+1; uMult++){
		hProjection = hSEmult_APAP->ProjectionX("hProjection",uMult+1,uMult+1);
		hTemp_SE->Add(hProjection,mult_weight_APAP->GetBinContent(uMult+1));
		delete hProjection;
	}
		
	hProjection = hMEmult_APAP->ProjectionX("hProjection",1,1);
	hTemp_ME = (TH1F*)hProjection->Clone("hTemp_ME");
	hTemp_ME->Scale(mult_weight_APAP->GetBinContent(1));
	delete hProjection;
	for(unsigned uMult=1; uMult<NumMultBins+1; uMult++){
		hProjection = hMEmult_APAP->ProjectionX("hProjection",uMult+1,uMult+1);
		hTemp_ME->Add(hProjection,mult_weight_APAP->GetBinContent(uMult+1));
		delete hProjection;
	}
	
	hCkmult_APAP = (TH1F*)hTemp_SE->Clone("hCkmult_APAP");
	hCkmult_APAP->Divide(hTemp_ME);
	delete hTemp_SE;
	delete hTemp_ME;
		
	//
	hProjection = hSEmult_SUM->ProjectionX("hProjection",1,1);
	hTemp_SE = (TH1F*)hProjection->Clone("hTemp_SE");
	hTemp_SE->Scale(mult_weight_SUM->GetBinContent(1));
	delete hProjection;
	for(unsigned uMult=1; uMult<NumMultBins+1; uMult++){
		hProjection = hSEmult_SUM->ProjectionX("hProjection",uMult+1,uMult+1);
		hTemp_SE->Add(hProjection,mult_weight_SUM->GetBinContent(uMult+1));
		delete hProjection;
	}
		
	hProjection = hMEmult_SUM->ProjectionX("hProjection",1,1);
	hTemp_ME = (TH1F*)hProjection->Clone("hTemp_ME");
	hTemp_ME->Scale(mult_weight_SUM->GetBinContent(1));
	delete hProjection;
	for(unsigned uMult=1; uMult<NumMultBins+1; uMult++){
		hProjection = hMEmult_SUM->ProjectionX("hProjection",uMult+1,uMult+1);
		hTemp_ME->Add(hProjection,mult_weight_SUM->GetBinContent(uMult+1));
		delete hProjection;
	}
	
	hCkmult_SUM = (TH1F*)hTemp_SE->Clone("hCkmult_SUM");
	hCkmult_SUM->Divide(hTemp_ME);
	delete hTemp_SE;
	delete hTemp_ME;
	
hCkmult_PP->Scale(double(Bin2-Bin1+1)/hCkmult_PP->Integral(Bin1,Bin2));
hCkmult_APAP->Scale(double(Bin2-Bin1+1)/hCkmult_APAP->Integral(Bin1,Bin2));
hCkmult_SUM->Scale(double(Bin2-Bin1+1)/hCkmult_SUM->Integral(Bin1,Bin2));

/*
	TFile* OutputFile = new TFile(NameOutputFolder+NameOutputFile,"recreate");
	for(unsigned uSce=0; uSce<3; uSce++){
		hSE_PP[uSce]->Write();
		hME_PP[uSce]->Write();
		hCk_PP[uSce]->Write();
		hSE_APAP[uSce]->Write();
		hME_APAP[uSce]->Write();
		hCk_APAP[uSce]->Write();		
	}
	hCk_WEIGHTED->Write();
	hSE_SUM->Write();
	hME_SUM->Write();
	hCk_SUM->Write();
hADD->Write();
delete hADD;

	hSEmult_PP->Write();
	hMEmult_PP->Write();
	hSEmult_APAP->Write();
	hMEmult_APAP->Write();
	hSEmult_SUM->Write();
	hMEmult_SUM->Write();
	hCkmult_PP->Write();
	hCkmult_APAP->Write();
	hCkmult_SUM->Write();
*/

	//TDirectoryFile* OutputDir = new TDirectoryFile(
	//TString::Format("Binning_%.0f",hCkmult_SUM->GetBinWidth(1)),
	//TString::Format("Binning_%.0f",hCkmult_SUM->GetBinWidth(1))
	//);
	OutputFile->cd();
	
	TDirectory* OutputDir = OutputFile->mkdir(TString::Format("Binning_%.0f",hCkmult_SUM->GetBinWidth(1)));
	OutputDir->cd();

	for(unsigned uSce=0; uSce<3; uSce++){
		OutputDir->Add(hSE_PP[uSce]);
		OutputDir->Add(hME_PP[uSce]);
		OutputDir->Add(hCk_PP[uSce]);
		OutputDir->Add(hSE_APAP[uSce]);
		OutputDir->Add(hME_APAP[uSce]);
		OutputDir->Add(hCk_APAP[uSce]);		
	}
	OutputDir->Add(hCk_WEIGHTED);
	OutputDir->Add(hSE_SUM);
	OutputDir->Add(hME_SUM);
	OutputDir->Add(hCk_SUM);

	OutputDir->Add(hSEmult_PP);
	OutputDir->Add(hMEmult_PP);
	OutputDir->Add(hSEmult_APAP);
	OutputDir->Add(hMEmult_APAP);
	OutputDir->Add(hSEmult_SUM);
	OutputDir->Add(hMEmult_SUM);
	OutputDir->Add(hCkmult_PP);
	OutputDir->Add(hCkmult_APAP);
	OutputDir->Add(hCkmult_SUM);
/*
	//delete
	for(unsigned uSce=0; uSce<3; uSce++){
		delete hSE_PP[uSce];
		delete hME_PP[uSce];
		delete hSE_APAP[uSce];
		delete hME_APAP[uSce];		
		delete hCk_PP[uSce];
		delete hCk_APAP[uSce];
	}
	delete [] hSE_PP;
	delete [] hME_PP;
	delete [] hSE_APAP;
	delete [] hME_APAP;	
	delete [] hCk_PP;
	delete [] hCk_APAP;
	
	delete hCk_WEIGHTED;
	
	delete hSE_SUM;
	delete hME_SUM;
	delete hCk_SUM;
	
	delete mult_weight_PP;
	delete mult_weight_APAP;
	delete mult_weight_SUM;
	
	delete hSEmult_PP;
	delete hMEmult_PP;
	delete hSEmult_APAP;
	delete hMEmult_APAP;
	delete hSEmult_SUM;
	delete hMEmult_SUM;

	delete InputFile;
	//delete OutputFile;
*/
	return OutputDir;
}

void CreateCk_from_FemtoDeam(){

	//pp, pL, pXim, LL
	//TString Correlation = "pp";
	//TString NameOutputFolder = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/CorrelationsDimi2020/";
	//TString NameOutputFile = "CkDimi2020_pp.root";
	//TString NameInputFolder = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/";
	//TString NameInputFile = "AnalysisResults.root";
	//TString NameList = "HMResults";
/*

	TString Correlation = "pL";
	TString NameOutputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";
	TString NameOutputFile = "CkDimi2020_pL.root";
	TString NameInputFolder = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/";
	TString NameInputFile = "AnalysisResults.root";
	TString NameList = "HMResults";

	TString Correlation = "pL";
	TString NameOutputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";
	TString NameOutputFile = "CkSidebands09062020_SL1_pL.root";
	TString NameInputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";	
	TString NameInputFile = "AnalysisResults_Sidebands09062020.root";
	TString NameList = "HMDimiResultsSL1";

	TString Correlation = "pL";
	TString NameOutputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";
	TString NameOutputFile = "CkSidebands09062020_SL2_pL.root";
	TString NameInputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";	
	TString NameInputFile = "AnalysisResults_Sidebands09062020.root";
	TString NameList = "HMDimiResultsSL2";

	TString Correlation = "pL";
	TString NameOutputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";
	TString NameOutputFile = "CkSidebands09062020_SR1_pL.root";
	TString NameInputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";	
	TString NameInputFile = "AnalysisResults_Sidebands09062020.root";
	TString NameList = "HMDimiResultsSR1";


	TString Correlation = "pL";
	TString NameOutputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";
	TString NameOutputFile = "CkSidebands09062020_SR2_pL.root";
	TString NameInputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";	
	TString NameInputFile = "AnalysisResults_Sidebands09062020.root";
	TString NameList = "HMDimiResultsSR2";
*/
	TString Correlation = "pL";
	TString NameOutputFolder = "/home/dmihaylov/CernBox/CatsFiles_Dimi/pLambda/";
	TString NameOutputFile = "CkDimi2020_pL.root";
	TString NameInputFolder = "/home/dmihaylov/CernBox/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/";
	TString NameInputFile = "AnalysisResults.root";
	TString NameList = "HMResults";
	const unsigned NumReb = 5;
	TFile* OutputFile = new TFile(NameOutputFolder+NameOutputFile,"recreate");
	TDirectory** dOutput = new TDirectory* [NumReb];
	
	for(unsigned uReb=0; uReb<NumReb; uReb++){
		printf("uReb=%u\n",uReb);
		dOutput[uReb] = WorkHorse(uReb+1,Correlation,OutputFile,
						NameInputFolder,NameInputFile,NameList);
		dOutput[uReb]->ls();
		OutputFile->cd();
		dOutput[uReb]->Write();
	}
	
	for(unsigned uReb=0; uReb<NumReb; uReb++){
		delete dOutput[uReb];
	}
	delete [] dOutput;
	delete OutputFile;
	
}
