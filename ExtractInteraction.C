R__ADD_LIBRARY_PATH("/home/dimihayl/Apps/CATS3/lib/");
R__LOAD_LIBRARY(libCATSbasic.so);
R__LOAD_LIBRARY(libCATSdev.so);
R__LOAD_LIBRARY(libCATSextended.so);
#include "/home/dimihayl/Apps/CATS3/include/DLM_CkModels.h"
#include "/home/dimihayl/Apps/CATS3/include/CATS.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_Source.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_Potentials.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_RootFit.h"

//const double& Momentum, const double* SourcePar, const double* PotPar
double LedniFit(double* x, double* par){
  return par[0]*Lednicky_Identical_Singlet(*x,&par[1],&par[2]);
}
//0 = norm
//1 = radius
//2,3,4,5 -> double gaussian potential
CATS* KITTY_FITTER=NULL;
double CATS_FITTER(double* x, double* par){
  KITTY_FITTER->SetAnaSource(0,par[1],true);
  KITTY_FITTER->SetShortRangePotential(0,0,0,par[2]);
  KITTY_FITTER->SetShortRangePotential(0,0,1,par[3]);
  KITTY_FITTER->SetShortRangePotential(0,0,2,par[4]);
  KITTY_FITTER->SetShortRangePotential(0,0,3,par[5]);
  KITTY_FITTER->KillTheCat();
  return par[0]*KITTY_FITTER->EvalCorrFun(*x);
}

bool Eval_ScattParameters(CATS& Kitty, double& ScatLen, double& EffRan, TH1F*& hFit, TF1*& fitSP, const int& Nterms=2, const bool& Fixf0=false, const bool& Fixd0=false){
  Kitty.KillTheCat();
  double* MomBins = Kitty.CopyMomBin();
  hFit = new TH1F("hFit","hFit",Kitty.GetNumMomBins(),MomBins);
  double LAST_POINT;
  double CURRENT_POINT;
  for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
    CURRENT_POINT = Kitty.GetMomentum(uMom)/tan(Kitty.GetPhaseShift(uMom,0,0));
    if(uMom){
      if(CURRENT_POINT*LAST_POINT<0&&fabs(CURRENT_POINT-LAST_POINT)>1000&&Kitty.GetMomentum(uMom)<120)
      {fitSP=NULL;delete[]MomBins;return false;}
    }
    hFit->SetBinContent(uMom+1,CURRENT_POINT);
    hFit->SetBinError(uMom+1,1.);
    LAST_POINT = CURRENT_POINT;
  }
  TF1* fitSP2 = new TF1("fitSP2","0.5*[1]/197.327*x*x+197.327*[0]", 10, 90);
  TF1* fitSP4 = new TF1("fitSP4","[2]*pow(x,4.)+0.5*[1]/197.327*x*x+197.327*[0]", 10, 90);
  TF1* fitSP6 = new TF1("fitSP6","[3]*pow(x,6.)+[2]*pow(x,4.)+0.5*[1]/197.327*x*x+197.327*[0]", 10, 90);
  double inv_f0 = ScatLen==0?0:1./ScatLen;
  //EffRan = 0;
//printf("ScatLen = %e\n",ScatLen);
//printf("inv_f0 = %f\n",inv_f0);
//printf("EffRan = %f\n",EffRan);
  if(Fixf0) {fitSP2->FixParameter(0,inv_f0);fitSP4->FixParameter(0,inv_f0);fitSP6->FixParameter(0,inv_f0);}
  else {fitSP2->SetParameter(0,inv_f0);fitSP2->SetParLimits(0,-100,100);
        fitSP4->SetParameter(0,inv_f0);fitSP4->SetParLimits(0,-100,100);
        fitSP6->SetParameter(0,inv_f0);fitSP6->SetParLimits(0,-100,100);}
  if(Fixd0) { fitSP2->FixParameter(1,EffRan);
              fitSP4->FixParameter(1,EffRan);
              fitSP6->FixParameter(1,EffRan);}
  else {fitSP2->SetParameter(1,EffRan);fitSP2->SetParLimits(1,0,50);
        fitSP4->SetParameter(1,EffRan);fitSP4->SetParLimits(1,0,50);
        fitSP6->SetParameter(1,EffRan);fitSP6->SetParLimits(1,0,50);}
  fitSP4->SetParameter(2,0);fitSP6->SetParameter(2,0);
  fitSP6->SetParameter(3,0);


  double Chi2_Old = 1e64;

  //hFit->Fit(fitSP2, "Q, S, N, R, M");
  ROOT::Math::MinimizerOptions MinOpt;
  MinOpt.SetMinimizerType("Minuit2");
  MinOpt.SetPrintLevel(0);
  DLM_FitHisto(hFit, fitSP2, "Q, S, N, R, M", "", &MinOpt);
  //printf("f0 %f\n", 1./fitSP2->GetParameter(0));
  ScatLen = 1./fitSP2->GetParameter(0);
  EffRan = fitSP2->GetParameter(1);
  if(Nterms<=2){delete fitSP4; delete fitSP6; fitSP=fitSP2; delete[]MomBins; return true;}

  hFit->Fit(fitSP4, "S, N, R, M");
  DLM_FitHisto(hFit, fitSP4, "Q, S, N, R, M", "", &MinOpt);
  if(fitSP4->GetChisquare()/fitSP2->GetChisquare()>0.8)
  {delete fitSP4; delete fitSP6; fitSP=fitSP2; delete[]MomBins; return true;}
  ScatLen = 1./fitSP4->GetParameter(0);
  EffRan = fitSP4->GetParameter(1);
  if(Nterms<=3){delete fitSP2; delete fitSP6; fitSP=fitSP4; delete[]MomBins; return true;}

  //hFit->Fit(fitSP6, "Q, S, N, R, M");
  DLM_FitHisto(hFit, fitSP6, "Q, S, N, R, M", "", &MinOpt);
  if(fitSP6->GetChisquare()/fitSP4->GetChisquare()>0.8)
  {delete fitSP2; delete fitSP6; fitSP=fitSP4; delete[]MomBins; return true;}
  ScatLen = 1./fitSP6->GetParameter(0);
  EffRan = fitSP6->GetParameter(1);
  delete fitSP2; delete fitSP4; fitSP=fitSP6; delete[]MomBins; return true;
}


//attempt to find scattering parameters or a potential that
//can describe (approx) a Gaussian shaped correlation function
void Rafa_2body_expCk(TString WorkFolder,TString InputFile,TString InputHisto){
  const double kMin = 0;
  const double kMax = 250;
  TFile fInput(InputFile,"read");
  TH1D* hInput  = (TH1D*)fInput.Get(InputHisto);
  TF1* fLedni = new TF1("fLedni",LedniFit,kMin,kMax,4);
  fLedni->SetParameter(0,1);
  fLedni->SetParameter(1,1.2);
  fLedni->SetParameter(2,2.0);
  fLedni->SetParameter(3,0.1);
  hInput->Fit(fLedni,"S, N, R, M");
  TH1F* hFit = new TH1F("hFit","hFit",int(kMax-kMin),kMin,kMax);
  for(unsigned uBin=0; uBin<hFit->GetNbinsX(); uBin++){
    hFit->SetBinContent(uBin+1,fLedni->Eval(hFit->GetBinCenter(uBin+1)));
  }
  TFile fOutput(WorkFolder+"FuckThisShit.root","recreate");
  hInput->Write();
  hFit->Write();
  fLedni->Write();
  delete hInput;
}

//same as above, but with a Gaussian potential
void Rafa_2body_expCk_CATS(TString WorkFolder,TString InputFile,TString InputHisto){
  const double kMin = 0;
  const double kMax = 240;
  CATS Kitty;
  Kitty.SetMomBins(TMath::Nint(kMax/2.),kMin,kMax);
  CATSparameters cPars (CATSparameters::tSource,1,true);
  cPars.SetParameter(0,1.2);
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetExcludeFailedBins(false);
  Kitty.SetQ1Q2(0);
  Kitty.SetQuantumStatistics(false);
  Kitty.SetRedMass( 0.5*938. );

  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0,1);
  Kitty.SetSpin(0,0);
  Kitty.SetChannelWeight(0, 1.);
  CATSparameters pPars(CATSparameters::tPotential,4,true);
  pPars.SetParameter(0,-50);
  pPars.SetParameter(1,1.0);
  pPars.SetParameter(2,0);
  pPars.SetParameter(3,1);
  Kitty.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  Kitty.SetNotifications(CATS::nWarning);
  Kitty.SetMaxRad(128);
  Kitty.SetMaxRho(64);
  Kitty.KillTheCat();
  KITTY_FITTER = &Kitty;

  TFile fInput(InputFile,"read");
  TH1D* hInput  = (TH1D*)fInput.Get(InputHisto);
  TF1* fCATS = new TF1("fCATS",CATS_FITTER,kMin,kMax,6);

  fCATS->SetParameter(0,1);
  fCATS->SetParameter(1,1.2);//r
  fCATS->SetParLimits(1,1.0,1.5);
  fCATS->SetParameter(2,-20);//V0
  fCATS->SetParLimits(2,-100,0);
  fCATS->SetParameter(3,1.0);//mu0
  fCATS->SetParLimits(3,0.5,5.0);
  fCATS->SetParameter(4,100.0);//V1
  fCATS->SetParLimits(4,0,1000);
  fCATS->SetParameter(5,0.1);//mu1
  fCATS->SetParLimits(5,0.05,0.5);

  hInput->Fit(fCATS,"S, N, R, M");
  TH1F* hFitResult = new TH1F("hFitResult","hFitResult",int(kMax-kMin),kMin,kMax);
  for(unsigned uBin=0; uBin<hFitResult->GetNbinsX(); uBin++){
    hFitResult->SetBinContent(uBin+1,fCATS->Eval(hFitResult->GetBinCenter(uBin+1)));
  }

  TH1F* h_kcotd=NULL;
  TF1* f_kcotd=NULL;
  double c_f0,c_d0=1;
  c_f0=0;
  c_d0=1;
  Kitty.SetEpsilonConv(1e-8);
  Kitty.SetEpsilonProp(1e-8);

  Eval_ScattParameters(Kitty,c_f0,c_d0,h_kcotd,f_kcotd,2,false,false);
  printf("f0 = %.2f; d0 = %.2f\n",c_f0,c_d0);
  printf("V0=%.2f; mu=%.2f; V1=%.2f; mu1=%.2f\n",fCATS->GetParameter(2),fCATS->GetParameter(3),fCATS->GetParameter(4),fCATS->GetParameter(5));

  TFile fOutput(WorkFolder+"FuckThisCat.root","recreate");
  hInput->Write();
  hFitResult->Write();
  fCATS->Write();
  h_kcotd->Write();
  f_kcotd->Write();

  delete fCATS;
  delete h_kcotd;
  delete f_kcotd;
  delete hInput;
}



void ExtractInteraction(){
	gSystem->Load("/home/dimihayl/Apps/CATS3/lib/libCATSbasic.so");
	gSystem->Load("/home/dimihayl/Apps/CATS3/lib/libCATSdev.so");
	gSystem->Load("/home/dimihayl/Apps/CATS3/lib/libCATSextended.so");

	//These is the file/histo that you have sent me
	TString WorkFolder = TString::Format("/home/dimihayl/Software/LocalFemto/Output/OtherTasks/Rafa_2body_expCk/");
  TString InputFile = WorkFolder+"femto.root";
  TString InputHisto = "CF2PartNorm";
	Rafa_2body_expCk(WorkFolder,InputFile,InputHisto);
	//Rafa_2body_expCk_CATS(WorkFolder,InputFile,InputHisto);
}
