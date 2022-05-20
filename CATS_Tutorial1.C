R__ADD_LIBRARY_PATH("/home/dimihayl/Apps/CATS3/lib/");
//R__ADD_LIBRARY_PATH("/usr/local/include/gsl/");
//R__ADD_LIBRARY_PATH("/usr/local/lib/");
R__LOAD_LIBRARY(libCATSbasic.so);
R__LOAD_LIBRARY(libCATSdev.so);
R__LOAD_LIBRARY(libCATSextended.so);

#include "/home/dimihayl/Apps/CATS3/include/CATS.h"
#include "/home/dimihayl/Apps/CATS3/include/CATSconstants.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_Source.h"
#include "/home/dimihayl/Apps/CATS3/include/CommonAnaFunctions.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_Histo.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_RootFit.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_Potentials.h"
#include "/usr/local/include/gsl/gsl_sf_dawson.h"


//all about pions
void Test1(){
  const unsigned NumMomBins = 200;
  const double kMin = 0;
  const double kMax = 400;

  CATS Kitty;
  Kitty.SetMomBins(NumMomBins,kMin,kMax);

  Kitty.SetUseAnalyticSource(true);
  CATSparameters source_func(CATSparameters::tSource, 1, true);
  //source_func.SetParameter(0,1.2);
  Kitty.SetAnaSource(GaussSource, source_func);
  Kitty.SetAnaSource(0,1.2*2);
  Kitty.SetAutoNormSource(false);
  Kitty.SetNormalizedSource(true);
  Kitty.SetMomentumDependentSource(false);
  Kitty.SetThetaDependentSource(false);
  Kitty.SetNumChannels(1);
  Kitty.SetSpin(0,0);
  //Kitty.SetSpin(1,1);
  Kitty.SetChannelWeight(0, 1./1.);
  //Kitty.SetChannelWeight(1, 3./4.);
  //Kitty.SetOnlyNumericalPw(0,true);
  Kitty.SetQuantumStatistics(0);
  Kitty.SetQ1Q2(0);
  Kitty.SetRedMass((Mass_pic * Mass_pic) / (Mass_pic + Mass_pic)); // reduced mass

  CATSparameters pot_func(CATSparameters::tPotential, 2, true);
  Kitty.SetNumPW(0,1);
  Kitty.SetShortRangePotential(0,0,SingleGauss,pot_func);
  Kitty.SetShortRangePotential(0,0,0,-500);
  Kitty.SetShortRangePotential(0,0,1,1);
  //Kitty.SetNotifications(CATS::nSilent);
  Kitty.SetEpsilonConv(1e-8);
  Kitty.SetEpsilonProp(1e-8);
  Kitty.KillTheCat();

  TGraph grCf;
  grCf.SetName("grCf");
  grCf.SetLineWidth(3);
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    grCf.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
  }

  TGraph grWfU;
  grWfU.SetName("grWfU");
  grWfU.SetLineWidth(6);
  grWfU.SetLineColor(kRed);

  TGraph grWfA;
  grWfA.SetName("grWfA");
  grWfA.SetLineWidth(4);
  grWfA.SetLineColor(kBlue);

  TGraph grWfRef;
  grWfRef.SetName("grWfRef");
  grWfRef.SetLineWidth(3);
  grWfRef.SetLineColor(kGreen);

  TGraph grWfTot2;
  grWfTot2.SetName("grWfTot2");
  grWfTot2.SetLineWidth(3);
  grWfTot2.SetLineColor(kMagenta);
  unsigned counter=0;
  for(double rad=0.05; rad<200; rad+=0.1){
    const double MOM = 40;
    grWfU.SetPoint(counter,rad,real(Kitty.EvalRadialWaveFunction(Kitty.GetMomBin(MOM),0,0,rad,false)));
    grWfA.SetPoint(counter,rad,real(Kitty.EvalAsymptoticRadialWF(Kitty.GetMomBin(MOM),0,0,rad,false)));
    grWfRef.SetPoint(counter,rad,real(Kitty.EvalReferenceRadialWF(Kitty.GetMomBin(MOM),0,rad,false)));
    grWfTot2.SetPoint(counter,rad,Kitty.EvalWaveFun2(Kitty.GetMomBin(MOM),rad,0));

    counter++;
  }

  TH1F* hPS = new TH1F("hPS","hPS",NumMomBins,kMin,kMax);
  TH1F* hkcotg = new TH1F("hkcotg","hkcotg",NumMomBins,kMin,kMax);
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    double PS = Kitty.GetPhaseShift(uBin,0,0);
    double kcotg = Kitty.GetMomentum(uBin)/tan(PS);
    hPS->SetBinContent(uBin,PS);
    hPS->SetBinError(uBin,1e-5);

    hkcotg->SetBinContent(uBin,kcotg);
    hkcotg->SetBinError(uBin,1e-5);
  }
  TF1* fitkcotg = new TF1("fitkcotg","197.327/[0]+0.5*[1]/197.327*x*x", 10, 90);
  hkcotg->Fit(fitkcotg, "S, N, R, M");

  TFile fOutput("CATS_Tutorial1.root","recreate");
  grCf.Write();
  grWfU.Write();
  grWfA.Write();
  grWfRef.Write();
  grWfTot2.Write();
  hPS->Write();
  hkcotg->Write();
  fitkcotg->Write();

  delete hPS;
}

void CATS_Tutorial1(){
  Test1();
}
