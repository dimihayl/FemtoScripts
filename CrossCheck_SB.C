
//we take two Ck from the grid, one the default, one with supposedly 100% purity and similar CPA
//from these two we should be able to obtain the expected shape of the SB correlation, and compare
//it with the SB we actually have extracted. This was used as a cross check during publication of pLambda
//for cleaner comparission, all is done in natural normalization
//in theory: C_SB = (C_exp - P*C_corrected)/(1-P)

//finding: the amount of fractions has a big influence on the expected correlation
//if there is a k* dependence this could fuck up a lot of stuff!
void CrossCheck_SB(){

  const double Purity = 0.953;
  //the values are for pT = 1.25 GeV -> 0.756 and 0.904
  //at 1.65 ->  0.768 and 0.890
  double PrimFrac_Original = 0.756;
  double PrimFrac_100 = 0.904;//904
  const TString OutputFileName = "./Output/CrossCheck_SB.root";

  const TString InputFile_Original = "$CERN_BOX/CatsFiles/ExpData/ALICE_pp_13TeV_HM/DimiJun20/Norm240_340/DataSignal/CkSB_pL_L55_SL4_SR6_P95_0.root";
  const TString InputHistoCk_Original = "hCk_12MeV";//C_exp
  const TString InputHistoCkNorm_Original = "hCk_Norm_12MeV";//C_exp
  const TString InputHistoCkS_Original = "hCkS_12MeV";//C_S
  const TString InputHistoCkB_Original = "hCkB_12MeV";//C_SB


  //100% purity based on V0 daughter cuts
  //originally I though this would not bias the fractions, but the resutl is weird,
  //so probably it does
  const TString InputFile_100 = "$CERN_BOX/CatsFiles_Dimi/pLambda/V0daughters/Natural/CkS4_pL_DS1040.root";
  const TString InputHistoCk_100 = "hCkmult_SUM";//C_corrected

  //this is the CPA. It does bias the fractions, but I can comppute by how much
  //and take this into account
  //const TString InputFile_100 = "$CERN_BOX/CatsFiles_Dimi/pLambda/CPA_purity/Natural/CkS4_pL_CPA1004.root";
  //const TString InputHistoCk_100 = "hCkmult_SUM";//C_corrected

  TDirectoryFile* dTemp;
  TH1F* hTemp;
  TH1D* hdTemp;
  TH1F* hCkOriginal;
  TH1F* hCkNormOriginal;
  TH1F* hCkSOriginal;
  TH1F* hCkBOriginal;
  TH1F* h100;
  TH1F* hSB;
  //expectation for h100 based on the differences in the lambda parameters,
  //evaluated at pT = 1.25 GeV with the TemplateFit_CPA (templates are the same)
  //the uncertainties are related to the feed-down of Xis
  TH1F* hExpected100;
  TH1F* hExpectedRatio;
  TFile* FileOriginal = new TFile(InputFile_Original,"read");;
  TFile* File100 = new TFile(InputFile_100,"read");;
  TFile* FileOutput = new TFile(OutputFileName,"recreate");;

  FileOriginal->cd();
  hTemp = (TH1F*)FileOriginal->Get(InputHistoCk_Original);
  FileOutput->cd();
  hCkOriginal = (TH1F*)hTemp->Clone("hCkOriginal");

  FileOriginal->cd();
  hTemp = (TH1F*)FileOriginal->Get(InputHistoCkNorm_Original);
  FileOutput->cd();
  hCkNormOriginal = (TH1F*)hTemp->Clone("hCkNormOriginal");

  hTemp = (TH1F*)FileOriginal->Get(InputHistoCkS_Original);
  FileOutput->cd();
  hCkSOriginal = (TH1F*)hTemp->Clone("hCkSOriginal");

  hTemp = (TH1F*)FileOriginal->Get(InputHistoCkB_Original);
  FileOutput->cd();
  hCkBOriginal = (TH1F*)hTemp->Clone("hCkBOriginal");

  File100->cd();
  dTemp = (TDirectoryFile*)File100->Get("Binning_12");
  dTemp->GetObject(InputHistoCk_100,hdTemp);
  FileOutput->cd();
  h100 = (TH1F*)hdTemp->Clone("h100");

  delete FileOriginal; FileOriginal=NULL;
  delete File100; File100=NULL;

  FileOutput->cd();
  hSB = (TH1F*)hCkOriginal->Clone("hSB");
  hExpected100 = (TH1F*)h100->Clone("hExpected100");
  for(unsigned uBin=0; uBin<hCkOriginal->GetNbinsX(); uBin++){
    if(uBin>=h100->GetNbinsX()) break;
    float val1,val2,err1,err2,val,err,mom;
    float ExpMax100,ExpMin100,valexp,errexp;
    mom = hCkOriginal->GetBinCenter(uBin+1);
    //1.25
    if(mom<200){
      PrimFrac_Original = 0.755;//0.756
      PrimFrac_100 = 0.901;//0.904
    }
    //1.75
    else{
      PrimFrac_Original = 0.761;//0.770
      PrimFrac_100 = 0.896;//0.864
    }
  //PrimFrac_Original = 0.756;
  //PrimFrac_100 = 0.904;
  //PrimFrac_Original = 0.770;
  //PrimFrac_100 = 0.864;
  PrimFrac_Original = 0.755;
  PrimFrac_100 = 0.755;
    val1 = hCkOriginal->GetBinContent(uBin+1);
    err1 = hCkOriginal->GetBinError(uBin+1);
    val2 = h100->GetBinContent(uBin+1);
    err2 = h100->GetBinError(uBin+1);
    val = (val1 - Purity*val2)/(1.-Purity);
    err = sqrt(err1*err1+err2*err2*Purity*Purity)/(1.-Purity);
    err=0;
    hSB->SetBinContent(uBin+1,val);
    hSB->SetBinError(uBin+1,err);

    val1 = hCkSOriginal->GetBinContent(uBin+1);

    //flat Cxi
    ExpMax100 = (PrimFrac_100*val1-(PrimFrac_100-PrimFrac_Original)*1.00)/PrimFrac_Original;
    //Cxi = 1.14 (roughly the value before lambda pars)
    ExpMin100 = (PrimFrac_100*val1-(PrimFrac_100-PrimFrac_Original)*1.14)/PrimFrac_Original;

    if(mom>120){
      //double CXI = 1.1-0.2*exp(-pow(mom*0.3/197.,2.));
      double CXI = 1;
      valexp = (PrimFrac_100*val1-(PrimFrac_100-PrimFrac_Original)*CXI)/PrimFrac_Original;
    }
    else if(mom>80){
      valexp = (PrimFrac_100*val1-(PrimFrac_100-PrimFrac_Original)*1.07)/PrimFrac_Original;
    }
    else{
      valexp = ExpMin100;
    }

    //valexp = (ExpMax100+ExpMin100)*0.5;
    //errexp = (ExpMax100-ExpMin100)*0.5;
    errexp=0;
    //errexp = err1;

    //gExpected100->SetPoint(uBin,hCkOriginal->GetBinCenter(uBin+1),valexp);
    //gExpected100->SetPointError(uBin,0,errexp);
    hExpected100->SetBinContent(uBin+1,valexp);
    hExpected100->SetBinError(uBin+1,errexp);
  }
  hExpectedRatio = (TH1F*)h100->Clone("hExpectedRatio");
  hExpectedRatio->Divide(hExpected100);


  double Norm = hCkNormOriginal->GetBinContent(3)/hCkOriginal->GetBinContent(3);

  hCkOriginal->Scale(Norm);
  hCkOriginal->SetLineColor(kRed+1);
  hCkOriginal->SetLineWidth(4);
  hCkOriginal->Write();

  hCkSOriginal->Scale(Norm);
  hCkSOriginal->SetLineColor(kRed+2);
  hCkSOriginal->SetLineWidth(3);
  hCkSOriginal->Write();

  hCkBOriginal->Scale(Norm);
  hCkBOriginal->SetLineColor(kAzure+1);
  hCkBOriginal->SetLineWidth(3);
  hCkBOriginal->Write();

  h100->Scale(Norm);
  h100->SetLineColor(kGreen+1);
  h100->SetLineWidth(3);
  h100->Write();

  hSB->Scale(Norm);
  hSB->SetLineColor(kViolet+1);
  hSB->SetLineWidth(3);
  hSB->Write();

  hExpected100->Scale(Norm);
  hExpected100->Write();
  hExpectedRatio->Write();

}
