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
  const unsigned NumMomBins = 100;

  CATS Kitty;
  Kitty.SetMomBins(NumMomBins,0,200);

  Kitty.SetUseAnalyticSource(true);
  CATSparameters source_func(CATSparameters::tSource, 1, true);
  //source_func.SetParameter(0,1.2);
  Kitty.SetAnaSource(GaussSource, source_func);
  Kitty.SetAnaSource(0,1.2);
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
  Kitty.SetShortRangePotential(0,0,0,-100);
  Kitty.SetShortRangePotential(0,0,1,1);
  //Kitty.SetNotifications(CATS::nSilent);

  Kitty.KillTheCat();

  TGraph grCf;
  grCf.SetName("grCf");
  grCf.SetLineWidth(3);

  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    grCf.SetPoint(uBin,Kitty.GetMomentum(uBin),Kitty.GetCorrFun(uBin));
  }

  TFile fOutput("CATS_Tutorial1.root","recreate");
  grCf.Write();
}

void CATS_Tutorial1(){
  Test1();
}
