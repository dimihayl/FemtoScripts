
void Create_pL_Ck(){

  const unsigned Rebin = 3;
  TString FolderName = "/home/dimihayl/CernBox/Sync/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/FemtoIntro/";
  TString FileName_PP_S = FolderName+"hPP_S.root";
  TString FileName_PP_R = FolderName+"hPP_R.root";
  TString FileName_AA_S = FolderName+"hAA_S.root";
  TString FileName_AA_R = FolderName+"hAA_R.root";

  TString HistoName_PP_S = "SEMultDist_Particle0_Particle2";
  TString HistoName_PP_R = "MEMultDist_Particle0_Particle2";
  TString HistoName_AA_S = "SEMultDist_Particle1_Particle3";
  TString HistoName_AA_R = "MEMultDist_Particle1_Particle3";

  TFile* fInput_PP_S = new TFile(FileName_PP_S,"read");
  TH2F* hPP_S = (TH2F*)fInput_PP_S->Get(HistoName_PP_S);
  hPP_S->Rebin2D(Rebin,1);

  TFile* fInput_PP_R = new TFile(FileName_PP_R,"read");
  TH2F* hPP_R = (TH2F*)fInput_PP_R->Get(HistoName_PP_R);
  hPP_R->Rebin2D(Rebin,1);

  TFile* fInput_AA_S = new TFile(FileName_AA_S,"read");
  TH2F* hAA_S = (TH2F*)fInput_AA_S->Get(HistoName_AA_S);
  hAA_S->Rebin2D(Rebin,1);

  TFile* fInput_AA_R = new TFile(FileName_AA_R,"read");
  TH2F* hAA_R = (TH2F*)fInput_AA_R->Get(HistoName_AA_R);
  hAA_R->Rebin2D(Rebin,1);

  const unsigned NumKstarBins = hPP_S->GetXaxis()->GetNbins();
  const unsigned NumMultBins = hPP_S->GetYaxis()->GetNbins();

  //evaluate the total yields of all samples, important to get the
  //normalized weights of each multiplicity bin
  const double TotalYield_PP_S = hPP_S->Integral();
  const double TotalYield_PP_R = hPP_R->Integral();
  const double TotalYield_AA_S = hAA_S->Integral();
  const double TotalYield_AA_R = hAA_R->Integral();
  //the total yield for PP and AA in the SE
  const double TotalYield_S = TotalYield_PP_S+TotalYield_AA_S;

  TH1F* hCk_tot = NULL;

  //without reweighting
  TH1F* hCk_simple = NULL;
  //the make the projection including under and overflow bins
  TH1F* hSE_simple = (TH1F*)hPP_S->ProjectionX("hSE_simple",0,NumMultBins+1);
  TH1F* hME_simple = (TH1F*)hPP_R->ProjectionX("hME_simple",0,NumMultBins+1);
  hME_simple->Scale(hSE_simple->Integral()/hME_simple->Integral());
  hCk_simple = (TH1F*)hSE_simple->Clone("hCk_simple");
  hCk_simple->Divide(hME_simple);

  TFile fOutput("Create_pL_Ck.root","recreate");

  //evaluate the correlation in each mult bin
  //again, both under and overflow are accounted for
  for(unsigned uMult=0; uMult<=NumMultBins+1; uMult++){
    //the projection function picks up the correlation in one mult bin (counts)
    TH1F* hSE_PP = (TH1F*) hPP_S->ProjectionX("hSE_PP",uMult,uMult);
    TH1F* hME_PP = (TH1F*) hPP_R->ProjectionX("hME_PP",uMult,uMult);
    TH1F* hSE_AA = (TH1F*) hAA_S->ProjectionX("hSE_AA",uMult,uMult);
    TH1F* hME_AA = (TH1F*) hAA_R->ProjectionX("hME_AA",uMult,uMult);

    const double Yield_SE_PP = hSE_PP->Integral();
    const double Yield_ME_PP = hME_PP->Integral();
    const double Yield_SE_AA = hSE_AA->Integral();
    const double Yield_ME_AA = hME_AA->Integral();
    //the yield with this mult bin (for both PP and AA)
    const double Yield_SE = Yield_SE_PP+Yield_SE_AA;
    //divide by the yield in all mult bins to get the weight of this mult bin
    const double WeightMult = Yield_SE/TotalYield_S;
    //if we want to put equal weights, i.e. no renormaliation
    const double AvgWeightMult = 1./double(NumMultBins);

    //each of these historam is scaled to its inverse integral
    //this way the total integral becomes 1, to have the properties of a
    //probability density function use the option "width".
    //The latter is optional here, as all histos have identical binning
    hSE_PP->Scale(1./Yield_SE_PP,"width");
    hME_PP->Scale(1./Yield_ME_PP,"width");
    hSE_AA->Scale(1./Yield_SE_AA,"width");
    hME_AA->Scale(1./Yield_ME_AA,"width");

    //divide the SE and ME to get the correlation function
    TH1F* hCk_PP = (TH1F*)hSE_PP->Clone("hCk_PP");
    hCk_PP->Divide(hME_PP);

    //check for division by zero in some bin
    for(unsigned uBin=0; uBin<hCk_PP->GetNbinsX(); uBin++){
      if(isnan(hCk_PP->GetBinContent(uBin+1))){
        hCk_PP->SetBinContent(uBin+1,0);
        hCk_PP->SetBinError(uBin+1,0);
      }
    }

    TH1F* hCk_AA = (TH1F*)hSE_AA->Clone("hCk_AA");
    hCk_AA->Divide(hME_AA);
    //check for division by zero in some bin
    for(unsigned uBin=0; uBin<hCk_AA->GetNbinsX(); uBin++){
      if(isnan(hCk_AA->GetBinContent(uBin+1))){
        hCk_AA->SetBinContent(uBin+1,0);
        hCk_AA->SetBinError(uBin+1,0);
      }
    }

    //add the PP and AA correlations according to their weights
    TH1F* hCk = (TH1F*)hCk_PP->Clone("hCk");
    if(Yield_SE){
      hCk->Scale(Yield_SE_PP/Yield_SE);
      hCk->Add(hCk_AA,Yield_SE_AA/Yield_SE);
    }

    //in the first iter we need to create this object
    if(!hCk_tot){
      hCk_tot = (TH1F*)hCk->Clone("hCk_tot");
      //we scale according to the weight
      hCk_tot->Scale(WeightMult);
    }
    //later on we just add next weighted term
    else{
      hCk_tot->Add(hCk,WeightMult);
    }

    delete hSE_PP;
    delete hME_PP;
    delete hSE_AA;
    delete hME_AA;
    delete hCk;
    delete hCk_PP;
    delete hCk_AA;
  }

  TH1F* hRatio = (TH1F*)hCk_tot->Clone("hRatio");
  hRatio->Divide(hCk_simple);

  fOutput.cd();
  hCk_tot->Write();
  hCk_simple->Write();
  hRatio->Write();

  //at this point Laura will ask you for a plot where the correlation
  //functions is 1 in some region. Lets pretend for us this is just above
  //the cust, so 300-350 MeV. There is a solution for you
  TH1F* hCk_tot_N = (TH1F*)hCk_tot->Clone("hCk_tot_N");
  TH1F* hCk_simple_N = (TH1F*)hCk_simple->Clone("hCk_simple_N");

  TF1* fNorm_tot = new TF1("fNorm_tot","[0]",0.300,0.350);
  fNorm_tot->SetParameter(0,1);
  hCk_tot_N->Fit(fNorm_tot,"S, N, R, M");
  hCk_tot_N->Scale(1./fNorm_tot->GetParameter(0));
  hCk_tot_N->Write();

  TF1* fNorm_simple = new TF1("fNorm_simple","[0]",0.300,0.350);
  fNorm_simple->SetParameter(0,1);
  hCk_simple_N->Fit(fNorm_simple,"S, N, R, M");
  hCk_simple_N->Scale(1./fNorm_simple->GetParameter(0));
  hCk_simple_N->Write();


  delete hCk_tot;
  delete hCk_simple;
  delete hCk_tot_N;
  delete hCk_simple_N;
  delete hRatio;
  delete hSE_simple;
  delete hME_simple;
  delete fInput_PP_S;
  delete fInput_PP_R;
  delete fInput_AA_S;
  delete fInput_AA_R;
}

void FemtoIntro(){
  Create_pL_Ck();
}
