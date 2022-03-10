R__ADD_LIBRARY_PATH("/home/dimihayl/Apps/CATS3/lib/");
R__LOAD_LIBRARY(libCATSbasic.so);
R__LOAD_LIBRARY(libCATSdev.so);
R__LOAD_LIBRARY(libCATSextended.so);

#include "/home/dimihayl/Apps/CATS3/include/CATS.h"
#include "/home/dimihayl/Apps/CATS3/include/CATSconstants.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_Source.h"
#include "/home/dimihayl/Apps/CATS3/include/CommonAnaFunctions.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_Histo.h"


void PlugInWaveFunction(){
  //N.B. the unit convetion for CATS : 'r' in fm, and 'k' in MeV

  //to be used for binning. In this example: 5 MeV bins
  const unsigned NumMomBins = 40;
  //N.B. CATS assumes t
  const double kstar_min = 0;
  const double kstar_max = 200;

  //later to be used for the radial relation of the wave function
  //this is important to be rather fine. The maximum value should be such that
  //the tail of the source function is negligible. As a rule of thumb, for
  //a Gaussian source this corresponds to c.a. 10 times the r0 value.
  const unsigned NumRadBins = 1000;
  const double r_min = 0;
  const double r_max = 50;

  //source size
  const double r0 = 3.0;

  //a container for the source parameters.
  //in this case trivial, as we will use 1D Gauss source (single parameter)
  CATSparameters cPars(CATSparameters::tSource,1,true);
  CATS Kitty;
  Kitty.SetMomBins(NumMomBins,kstar_min,kstar_max);
  //set the source function and source parameters (size)
  //definition of GaussSource in DLM_Source.h
  Kitty.SetAnaSource(GaussSource, cPars);
  Kitty.SetAnaSource(0,r0);
  Kitty.SetUseAnalyticSource(true);
  Kitty.SetMomentumDependentSource(false);

  //charge, q1*q2, for same charged particle 1
  //for opposite charged particle -1, if one particle is neutral 0
  //N.B. For p-Xi change to -1
  Kitty.SetQ1Q2(0);
  //set true to account for QS of identical particles
  Kitty.SetQuantumStatistics(false);
  Kitty.SetRedMass( (Mass_p*Mass_Xim)/(Mass_p+Mass_Xim) );

  //in case of different spin/isospin channels, they can be separately defined here
  //in this example, we only define a single effective channel, that contains
  //a single partial wave (PW), the s-wave, that is deviating from the asymptotic form
  Kitty.SetNumChannels(1);
  Kitty.SetNumPW(0,1);
  //dummy value in this case, only relevant in case of identical particle
  Kitty.SetSpin(0,0);
  //as we have 1 channel, we take it with weight of 1 (100%)
  Kitty.SetChannelWeight(0,1.);

  //the wave function has to be passed using some custom made histograms
  //the histograms should be 2D, with one axis corresponding to momentum k (x)
  //and the other to the distance r (y). The value (z-axis) corresponds to the
  //value of the radial partial wave functions (k,r) for the corresponding channel
  //in this example: single channel and single PW, thus a single wave function is needed
  DLM_Histo<complex<double>> histWF;
  //two dimensions
  histWF.SetUp(2);
  //how many bins/range in each dimension
  histWF.SetUp(0,NumMomBins,kstar_min,kstar_max);
  histWF.SetUp(1,NumRadBins,r_min,r_max);
  histWF.Initialize();

  //we can also provide phase shifts, but this is not needed, so concider
  //this histogram a dummy
  DLM_Histo<complex<double>> histPS(histWF);

  //example how to set the bin values of the wave function histogram
  //in this example, we will simply set it equal to the Bessel function,
  //so we should get a flat correlation as it corresponds to no interaction
  //in your case, you would need to include the wave function from the Lednicky model.
  //N.B. In case of Coulomb, the correction has to be included in your solution,
  //CATS will NOT add it on top of the solution you provide
  unsigned axis[2];
  for(axis[0]=0; axis[0]<histWF.GetNbins(0); axis[0]++){//over kstar
    for(axis[1]=0; axis[1]<histWF.GetNbins(1); axis[1]++){//over radial dependence
      double Momentum = histWF.GetBinCenter(0,axis[0]);
      double Radius = histWF.GetBinCenter(1,axis[1]);
      double Rho = Radius/hbarc*Momentum;
      //the value here needs to be u = r*R, where R is the radial partial wave function
      //and the radial parameter 'r' has to be given in units of fermi
      complex<double> CPW = sin(Rho)/Rho*Radius;
      histWF.SetBinContent(axis,CPW);
    }
  }

  Kitty.SetExternalWaveFunction(0,0,histWF,histPS);
  Kitty.KillTheCat();

  TH1F* h_CorrFunction = new TH1F("h_CorrFunction","h_CorrFunction",NumMomBins,kstar_min,kstar_max);
  for(unsigned uBin=0; uBin<NumMomBins; uBin++){
    h_CorrFunction->SetBinContent(uBin+1,Kitty.GetCorrFun(uBin));
  }
  TFile fCatsExample("fCatsExample.root","recreate");
  h_CorrFunction->Write();
}


void CatsExample(){

  PlugInWaveFunction();

}
