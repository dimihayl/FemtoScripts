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



//SET(GSL_INCLUDE "/usr/local/include/gsl")#where are all GSL related .h files
//SET(GSL_LIB "/usr/local/lib")#where are the GSL .a and .so files

int Get_MMclass(const double& r0, const double& f0, const double& d0){
  //f0>0
  //r0>d0
  //|f0|>2*|d0| (BS cond) -> satisfied for 2 dominant or 3 surpressed
  //r0>|f|/2
  //interesting classes:
  //+/- AB, with A,B = 1,2,3 (r0,f0,d0)
  //A = which one is dominant (e.g. 1 if r0>d0 && r0>|f|/2)
  //B = which one is surpressed (e.g. 3 if d0<r0 && d0<|f0|/2)
  int CLASS = 0;
  if(r0<fabs(d0)&&r0<fabs(f0)*0.5) CLASS = 1;
  else if(fabs(f0)<2.*r0&&fabs(f0)<2.*fabs(d0)) CLASS = 2;
  else if(fabs(d0)<r0&&fabs(d0)<fabs(f0)*0.5) CLASS = 3;

  if(r0>fabs(d0)&&r0>fabs(f0)*0.5) CLASS += 10;
  else if(fabs(f0)>2.*r0&&fabs(f0)>2.*fabs(d0)) CLASS += 20;
  else if(fabs(d0)>r0&&fabs(d0)>fabs(f0)*0.5) CLASS += 30;
  CLASS = f0>0?CLASS:-CLASS;
  return CLASS;
}


//DimiFlag -> any further info. At the moment:
//0 - randomly generated potential
//10 - ManufacturePotential used (random phase)
//11 - ManufacturePotential used (converged phase)
void Fill_ntMM(TNtuple* ntMM, const CATS& Kitty, const double& r0, const double& f0, const double& d0,
                const double& V1, const double& mu1, const double& V2, const double& mu2, const double& s2, const int& DimiFlag){
    //ROOT::Math::MinimizerOptions MinOpt;
    //MinOpt.SetMinimizerType("Minuit2");
    //MinOpt.SetPrintLevel(0);

    const unsigned nMom = Kitty.GetNumMomBins();
    double* MomBins = Kitty.CopyMomBin();
    const double kCorrection = 120;
    TH1F* hMMX = new TH1F("hMMX","hMMX",nMom,MomBins);
    TF1* fMMX = new TF1("fMMX","[0]+[1]*sqrt(1.+[2]*x*x)",Kitty.GetMomBinLowEdge(0),kCorrection);
    fMMX->SetParameter(0,1);
    fMMX->SetParLimits(0,-40,40);
    fMMX->SetParameter(1,0);
    fMMX->SetParLimits(1,-40,40);
    fMMX->FixParameter(2,fabs(f0*d0*FmToNu*FmToNu*0.5));
    //the weight of the f^2 factor with respect to F1F2 (after correction)
    double wf2=0;
    const double if0=f0?1./f0:0;
    Float_t ntBuffer[16];
    for(unsigned uMom=0; uMom<nMom; uMom++){
      double MOM = Kitty.GetMomentum(uMom);
      double F1 = gsl_sf_dawson(2.*MOM*r0*FmToNu)/(2.*MOM*r0*FmToNu);
      double F2 = (1.-exp(-4.*MOM*MOM*r0*r0*FmToNu*FmToNu))/(2.*MOM*r0*FmToNu);
      complex<double> Val_f = pow(if0/FmToNu+0.5*d0*FmToNu*MOM*MOM-i*MOM,-1.);
      double v_f2 = 0.5*pow(abs(Val_f)/(r0*FmToNu),2.);
      double F1F2 = 2*real(Val_f)*F1/(sqrt(Pi)*r0*FmToNu)-imag(Val_f)*F2/(r0*FmToNu);
      //double Ck_Ledni = 1.+v_f2+F1F2;
      double Ck_Cats = Kitty.GetCorrFun(uMom);
      //the correction factor we want to study
      double MMX = (Ck_Cats-1.-F1F2)/v_f2;
      hMMX->SetBinContent(uMom+1,MMX);
      hMMX->SetBinError(uMom+1,0.001);
      wf2 += fabs((Ck_Cats-1.-F1F2)/F1F2);
    }
    wf2 /= double(nMom);
    //hMMX->Fit(fMMX,"Q, S, N, R, M");
    ROOT::Math::MinimizerOptions MinOpt;
    MinOpt.SetMinimizerType("Minuit2");
    MinOpt.SetPrintLevel(0);
    DLM_FitHisto(hMMX, fMMX, "Q, S, N, R, M", "", &MinOpt);
    double csmall = (1-(d0)/(2*sqrt(Pi)*r0));
    ntBuffer[0] = r0;
    ntBuffer[1] = f0;
    ntBuffer[2] = d0;
    ntBuffer[3] = fMMX->GetParameter(0);
    ntBuffer[4] = fMMX->GetParameter(1);
    ntBuffer[5] = hMMX->GetBinContent(1);
    ntBuffer[6] = hMMX->GetBinContent(nMom);
    ntBuffer[7] = csmall;
    ntBuffer[8] = V1;
    ntBuffer[9] = mu1;
    ntBuffer[10] = V2;
    ntBuffer[11] = mu2;
    ntBuffer[12] = s2;
    ntBuffer[13] = wf2;
    ntBuffer[14] = Get_MMclass(r0,f0,d0);
    ntBuffer[15] = DimiFlag;

    ntMM->Fill(ntBuffer);
    delete hMMX; delete fMMX;
    delete [] MomBins;
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
  else {fitSP2->SetParameter(1,EffRan);fitSP2->SetParLimits(1,-50,50);
        fitSP4->SetParameter(1,EffRan);fitSP4->SetParLimits(1,-50,50);
        fitSP6->SetParameter(1,EffRan);fitSP6->SetParLimits(1,-50,50);}
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


//evaluated w/o coulomb and for red mass of p-phi (so around 490 MeV)
void OkayishStartingPars( const TString Potential, const double& f0, const double d0,
                          double& V1, double& mu1, double& V2, double& mu2, double& s2){
  V1 = 0.1;
  mu1 = 1;
  V2 = 0.1;
  mu2 = 0.5;
  if(f0>0&&fabs(f0)>fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){

    }
    else if(Potential=="YukawaDimiCore"){
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){

    }
    else if(Potential=="Yukawa"){

    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }
  }
  else if(f0>0&&fabs(f0)<=fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){
      //f0 = 0.993 fm
      //d0 = 3.97 fm
      V1 = -8.576280e+02;
      mu1 = 6.943480e-01;
      V2 = 7.481622e+03;
      mu2 = 3.903410e-01;
      s2 = 0;
    }
    else if(Potential=="YukawaDimiCore"){
      //f0 = 1.027 fm
      //d0 = 4.19 fm
      V1 = 3.455340e+00;
      mu1 = 3.881920e-01;
      V2 = 2.211957e+02;
      mu2 = 6.202096e-01;
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){
      //f0 = 0.976 fm
      //d0 = 3.77 fm
      V1 = -2.469257e+01;
      mu1 = 1.301668e+00;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -1;
    }
    else if(Potential=="Yukawa"){
      s2 = -2;
    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }
  }
  else if(f0<0&&fabs(f0)>fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){
      //f0 = -3.976 fm
      //d0 = 1.01 fm
      V1 = -7.414734e+01;
      mu1 = 1.111177e+00;
      V2 = -6.127154e+02;
      mu2 = 3.208735e-01;
      s2 = 0;
    }
    else if(Potential=="YukawaDimiCore"){
      //f0 = -4.094 fm
      //d0 = 1.04 fm
      V1 = 1.246311e+01;
      mu1 = 2.290515e-01;
      V2 = -8.227462e+01;
      mu2 = 9.258949e-01;
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){
      //f0 = -3.970 fm
      //d0 = 1.05 fm
      V1 = -1.881959e+02;
      mu1 = 8.635644e-01;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -1;
    }
    else if(Potential=="Yukawa"){
      //f0 = -3.837 fm
      //d0 = 1.05 fm
      //V1 = 8.792615e-01;
      //mu1 = 1.045134e+02;
      //V2 = 0.000000e+00;
      //mu2 = 0.000000e+00;
      s2 = -2;
    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }

  }
  else if(f0<0&&fabs(f0)<=fabs(2.*d0)){
    if(Potential=="DoubleGaussSum"){
      //f0 = -0.948 fm
      //d0 = 0.95 fm
      V1 = -1.418785e+03;
      mu1 = 1.123364e+00;
      V2 = -1.147570e+03;
      mu2 = 1.031657e+00;
      s2 = 0;
    }
    else if(Potential=="YukawaDimiCore"){
      //f0 = -1.048 fm
      //d0 = 1.04 fm
      V1 = 1.945222e+01;
      mu1 = 4.377242e-01;
      V2 = 1.085863e+03;
      mu2 = 1.111706e+00;
      s2 = mu2*0.2;
    }
    else if(Potential=="Gaussian"){
      //f0 = -1.054 fm
      //d0 = 1.02 fm
      V1 = -1.260122e+02;
      mu1 = 1.991708e+00;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -1;
    }
    else if(Potential=="Yukawa"){
      //f0 = -1.042 fm
      //d0 = 1.04 fm
      V1 = 7.584738e+00;
      mu1 = 5.422729e-01;
      V2 = 0.000000e+00;
      mu2 = 0.000000e+00;
      s2 = -2;
    }
    else{
      V1 = 0;
      mu1 = 0;
      V2 = 0;
      mu2 = 0;
      s2 = 0;
    }
  }
  else{
    printf(" OkayishStartingPars says f*** !?\n");
  }




}



void Write_MM_Plots(TFile*fOutput, const CATS& Kitty, const double& r0, const double& f0, const double& d0, const TString& suffix){
    fOutput->cd();
    TGraph gCk_SE;
    gCk_SE.SetName("gCk"+suffix);
    gCk_SE.SetLineWidth(3);
    gCk_SE.SetLineColor(kGreen);
    TGraph gCk_LL;
    gCk_LL.SetName("gCk_LL"+suffix);
    gCk_LL.SetLineWidth(3);
    gCk_LL.SetLineColor(kBlack);
    TGraph gCk_LLX;
    gCk_LLX.SetName("gCk_LLX"+suffix);
    gCk_LLX.SetLineWidth(3);
    gCk_LLX.SetLineColor(kRed+2);
    TGraph g_LLX;
    g_LLX.SetName("g_LLX"+suffix);
    g_LLX.SetLineWidth(2);
    g_LLX.SetLineColor(kRed);
    TGraph g_f2;
    g_f2.SetName("g_f2"+suffix);
    g_f2.SetLineWidth(2);
    g_f2.SetLineColor(kGreen+1);
    TGraph g_f2X;
    g_f2X.SetName("g_f2X"+suffix);
    g_f2X.SetLineWidth(2);
    g_f2X.SetLineColor(kYellow+2);
    TGraph g_F1F2;
    g_F1F2.SetName("g_F1F2"+suffix);
    g_F1F2.SetLineWidth(2);
    g_F1F2.SetLineColor(kBlue);

    TGraph g_F1F2_f2X;
    g_F1F2_f2X.SetName("g_F1F2_f2X"+suffix);
    g_F1F2_f2X.SetLineWidth(2);
    g_F1F2_f2X.SetLineColor(kGray+1);

    TGraph g_PS;
    g_PS.SetName("g_PS"+suffix);
    g_PS.SetLineWidth(4);
    g_PS.SetLineStyle(1);
    g_PS.SetLineColor(kPink+2);

    TGraph g_kcot;
    g_kcot.SetName("g_kcot"+suffix);
    g_kcot.SetLineWidth(4);
    g_kcot.SetLineStyle(2);
    g_kcot.SetLineColor(kPink+2);

    TGraph g_sqrtkcot;
    g_sqrtkcot.SetName("g_sqrtkcot"+suffix);
    g_sqrtkcot.SetLineWidth(4);
    g_sqrtkcot.SetLineStyle(3);
    g_sqrtkcot.SetLineColor(kPink+2);


    //the true correction
    TGraph g_MMX;
    g_MMX.SetName("g_MMX"+suffix);
    g_MMX.SetLineWidth(4);
    g_MMX.SetLineStyle(2);
    g_MMX.SetLineColor(kPink+1);

    TGraph g_MMB;
    g_MMB.SetName("g_MMB"+suffix);
    g_MMB.SetLineWidth(4);
    g_MMB.SetLineColor(kAzure+10);

    TGraph g_MMA;
    g_MMA.SetName("g_MMA"+suffix);
    g_MMA.SetLineWidth(4);
    g_MMA.SetLineColor(kViolet+1);

    unsigned NumPts = 0;
    double if0 = 1./f0;;
    for(unsigned uMom=0; uMom<Kitty.GetNumMomBins(); uMom++){
      double MOM = Kitty.GetMomentum(uMom);
      double F1 = gsl_sf_dawson(2.*MOM*r0*FmToNu)/(2.*MOM*r0*FmToNu);
      double F2 = (1.-exp(-4.*MOM*MOM*r0*r0*FmToNu*FmToNu))/(2.*MOM*r0*FmToNu);
      complex<double> Val_f = pow(if0/FmToNu+0.5*d0*FmToNu*MOM*MOM-i*MOM,-1.);
      double LLX = (1-(d0)/(2*sqrt(Pi)*r0));
      double v_f2 = 0.5*pow(abs(Val_f)/r0/FmToNu,2);
      double F1F2 = 2*real(Val_f)*F1/(sqrt(Pi)*r0*FmToNu)-imag(Val_f)*F2/r0/FmToNu;

      double Ck_LL = 1.+v_f2+F1F2;
      double Ck_LLX = 1.+v_f2*LLX+F1F2;

      gCk_LL.SetPoint(NumPts,MOM,Ck_LL);
      gCk_LLX.SetPoint(NumPts,MOM,Ck_LLX);
      g_LLX.SetPoint(NumPts,MOM,LLX);
      g_f2.SetPoint(NumPts,MOM,v_f2);
      g_f2X.SetPoint(NumPts,MOM,v_f2*LLX);
      g_F1F2.SetPoint(NumPts,MOM,F1F2);
      g_F1F2_f2X.SetPoint(NumPts,MOM,fabs(F1F2/v_f2/LLX));

      gCk_SE.SetPoint(NumPts,MOM,Kitty.GetCorrFun(NumPts));
      g_PS.SetPoint(NumPts,MOM,Kitty.GetPhaseShift(NumPts,0,0));
      g_kcot.SetPoint(NumPts,MOM,MOM/tan(Kitty.GetPhaseShift(NumPts,0,0)));
      g_sqrtkcot.SetPoint(NumPts,MOM,sqrt(fabs(MOM/tan(Kitty.GetPhaseShift(NumPts,0,0)))));

      //the correction as in ledni
      double MMX = (Kitty.GetCorrFun(NumPts)-1.-F1F2)/v_f2;
      g_MMX.SetPoint(NumPts,MOM,MMX);

      NumPts++;
    }


    TGraph g_pot;
    g_pot.SetName("g_pot");
    g_pot.SetLineWidth(4);
    g_pot.SetLineStyle(2);
    g_pot.SetLineColor(kBlue-1);
    unsigned uRad=0;
    for(double RAD=0.1; RAD<5; RAD+=0.025){
      g_pot.SetPoint(uRad++,RAD,Kitty.EvaluateThePotential(0,0,50,RAD));
    }

    gCk_SE.Write();
    gCk_LL.Write();
    gCk_LLX.Write();
    g_LLX.Write();
    g_f2.Write();
    g_f2X.Write();
    g_F1F2.Write();
    g_F1F2_f2X.Write();
    g_PS.Write();
    g_kcot.Write();
    g_sqrtkcot.Write();
    g_MMX.Write();
    g_pot.Write();
}



//Potential==Dynamic, means we try to do the double Gaussian, if it fails we repeat it all with Yukawa
//FineTune < 1 to decrease the intrisic minimum step for the current potential
void ManufacturePotential(const double f0, const double df0,
                                const double d0, const double dd0, const double* Radii, const unsigned NumRad,
                                double& V1, double& mu1, double& V2, double& mu2,
                                const TString Potential, const TString OutputFolder, const int& SEED=11, const double* StartPars=NULL,
                                const double RedMass=500, const double FineTune=1){


  //b = best, l = last, d = difference to desired
  TString CurrentPot = Potential;
  if(Potential=="Dynamic") CurrentPot = "DoubleGaussSum";
  double V_1,bV_1;
  double V_2,bV_2;
  double mu_1,bmu_1;
  double mu_2,bmu_2;
  double s_2;
  double f_0,lf_0,bf_0,df_0,bdf_0;
  double d_0,ld_0,bd_0,dd_0,bdd_0;
  bool Starting = true;
  bool Started = false;
  //dist,f0,d0,V1,mu1,V2,mu2: used when VeryStuck>=2 resets the whole thing
  double FallBack[8];FallBack[0]=1e16;
  f_0 = 0; d_0 = 0;
  bV_1 = 0; bV_2 = 0;
  bmu_1 = 1.5; bmu_2 = 0.5;
  bdf_0 = 1e6; bdd_0 = 1e6;
  bf_0 = f0; bd_0 = d0;
  lf_0=f0; ld_0=d0;
  const unsigned MaxIter = 10000;
  unsigned uIter=0;
  TRandom3 rangen(SEED);
  CATSparameters sPars(CATSparameters::tSource,1,true);
  sPars.SetParameter(0,2.5);
  const double kStepFine=3;
  const double kStepCoarse=10;
  const unsigned nMomFine = 40;
  const unsigned nMomCoarse = 18;
  const unsigned nMom = nMomFine+nMomCoarse;
  double* MomBins = new double [nMom+1];
  for(unsigned uMom=0; uMom<nMomFine; uMom++){
    MomBins[uMom] = kStepFine*double(uMom);
  }
  for(unsigned uMom=0; uMom<=nMomCoarse; uMom++){
    MomBins[nMomFine+uMom] = kStepFine*double(nMomFine)+kStepCoarse*double(uMom);
  }
  CATS Kitty_SE;
  Kitty_SE.SetMomBins(nMom,MomBins);
  Kitty_SE.SetAnaSource(GaussSource, sPars);
  Kitty_SE.SetUseAnalyticSource(true);
  Kitty_SE.SetQ1Q2(0);
  Kitty_SE.SetQuantumStatistics(false);
  Kitty_SE.SetRedMass( RedMass );
  //Kitty_SE.SetRedMass(Mass_L*0.5);
//Kitty_SE.SetRedMass(488.6);//pphi
  Kitty_SE.SetNumChannels(1);
  Kitty_SE.SetNumPW(0,1);
  Kitty_SE.SetSpin(0,0);
  Kitty_SE.SetChannelWeight(0, 1.);
  Kitty_SE.SetEpsilonConv(1e-8);
  Kitty_SE.SetEpsilonProp(1e-8);
  //Kitty_SE.SetMaxRad(96);
  //Kitty_SE.SetMaxRho(32);
  Kitty_SE.SetNotifications(CATS::nSilent);
  //Kitty_SE.SetGridEpsilon(1./512.);
  double RadForGrid;
  Kitty_SE.KillTheCat();
  CATSparameters pPars(CATSparameters::tPotential,5,true);
  if(CurrentPot=="YukawaDimiCore") Kitty_SE.SetShortRangePotential(0,0,YukawaDimiCore,pPars);
  else if(CurrentPot=="Gassian") Kitty_SE.SetShortRangePotential(0,0,Gaussian,pPars);
  else if(CurrentPot=="Yukawa") Kitty_SE.SetShortRangePotential(0,0,YukawaDimiSmooth,pPars);
  else Kitty_SE.SetShortRangePotential(0,0,DoubleGaussSum,pPars);
  bool GoodGoing=false;
  double dist=1e16;
  double bdist=1e16;
  //the scale to which 1fm roughly correponds
  const double dist_unit = sqrt(0.5/df0/df0+0.5/dd0/dd0);

  //Fluctuations: should be the same order as bdist
  //if too large -> decrease Convergence
  //if too small -> increse Convergence
  const unsigned fluctN = 8;
  double fluct[fluctN];
  unsigned ufluct = 0;
  double avg_fluct;
  const unsigned dirN = 32;
  //how many of the last dirN iterations moved towards the true solution,
  //which is defined as +/-45 deg angle on the 1/f0 d0 plane (1/4 of the angle)
  //we expect 4 movements on avg, if its less than that we increase Convergence
  bool dir[dirN];
  unsigned udir = 0;
  double avg_dir;
  //double Adjust = 10;
  double Convergence[5];
  for(unsigned uPar=0; uPar<5; uPar++)Convergence[uPar]=1;
  unsigned Stuck = 0;
  unsigned BadPhaseShifts = 0;
  unsigned UnstuckCounter = 0;
  const unsigned StuckLimit = 64;
  //how many iterations to do by changing only single pars
  const unsigned UnstuckPerPar = 32*2;
  //1 is V1, 2 is mu1, 3 is V2 and 4 is mu2. 0 the normal mode
  unsigned SinglePar = 0;
  for(unsigned uf=0; uf<fluctN; uf++) fluct[uf]=-1000*dist_unit;
  for(unsigned ud=0; ud<dirN; ud++) dir[ud]=true;
//TFile fDump(TString::Format("%s/OtherTasks/ManufactureYukawaPotential/fDump.root",GetFemtoOutputFolder()),"recreate");
const bool DEBUG = false;
  unsigned VeryStuck=0;
  TFile fOutput(TString::Format("%s/fOut_%.2f_%.1f_%s.root",OutputFolder.Data(),f0,d0,Potential.Data()),"recreate");
  TNtuple* ntMM = new TNtuple("ntMM", "ntMM","r0:f0:d0:MMA:MMB:C_X0:C_X150:L_X0:V1:mu1:V2:mu2:s2:wf2:Class:Manufactured");
  TH1F* hClassCount = new TH1F("hClassCount","hClassCount",100,-50.5,49.5);
  TH2F* hif0d0Count = new TH2F("hf0d0Count","hf0d0Count",20,-10,10,20,0,40);
  bool FoundIt = false;
  while(uIter<MaxIter&&!FoundIt){
    if(!DEBUG){
      printf("\r\033[K Goal: %s (f0,d0) [Achieved]: (%.3f,%.3f)+/-(%.3f,%.3f) [%3.0f%%,%3.0f%%], Break (%5u/%5u)",Potential.Data(),f0,d0,fabs(df0),fabs(dd0),fabs(df0/bdf_0*100.),fabs(dd0/bdd_0*100.),uIter,MaxIter);
    }
    else printf("\nProgress: Break (%5u/%5u), Precision f0(%3.0f%%) d0(%3.0f%%)",uIter,MaxIter,fabs(df0/bdf_0*100.),fabs(dd0/bdd_0*100.));
    //printf("Progress %5u from %5u, bd=%.1f(%.1f)\n",uIter,MaxIter,bdist,sqrt(0.5*bdf_0*bdf_0/df0/df0+0.5*bdd_0*bdd_0/dd0/dd0));

    avg_fluct = 0;
    for(unsigned uf=0; uf<fluctN; uf++) avg_fluct+=fluct[uf];
    avg_fluct /= double(fluctN);
    double desired_fluct = bdist;

    avg_dir = 0;
    for(unsigned ud=0; ud<dirN; ud++) avg_dir+=double(dir[ud]);
    avg_dir /= double(dirN);

    if(DEBUG){
      printf("\n B: Stuck=%u(%u); UnstuckCounter=%u; SinglePar=%u",Stuck,VeryStuck,UnstuckCounter,SinglePar);
    }

    unsigned Npars = 4;
    if(CurrentPot=="YukawaDimiCore") Npars=4;
    if(CurrentPot=="Gaussian") Npars=2;
    if(CurrentPot=="Yukawa") Npars=2;
    else Npars=4;

    if(Stuck==StuckLimit&&SinglePar==0){
      VeryStuck=1;
      UnstuckCounter=0;
      Stuck=0;
      SinglePar=1;
      for(unsigned uPar=1;uPar<=Npars;uPar++)Convergence[uPar]*=Convergence[0];
      Convergence[0]=1;
      desired_fluct /= 2.;
    }
    if(BadPhaseShifts==10) {VeryStuck=3;Stuck=100000;UnstuckCounter=UnstuckPerPar*Npars;}
    if(UnstuckCounter==UnstuckPerPar*Npars){
      //if VeryStuck goes two times: reset and start from scratch
      if(Stuck<UnstuckCounter){
        VeryStuck=0;
        UnstuckCounter=0;
        Stuck=0;
        SinglePar=0;
        Convergence[0]=0;
        //make it the max value
        for(unsigned uPar=1;uPar<=Npars;uPar++)
          if(Convergence[uPar]>Convergence[0])Convergence[0]=Convergence[uPar];
        for(unsigned uPar=1;uPar<=Npars;uPar++)Convergence[uPar]/=Convergence[0];
      }
      else if(VeryStuck<3){
        VeryStuck = 1+Stuck/UnstuckCounter;
        desired_fluct /= (2.*VeryStuck);
        UnstuckCounter=0;
        //Stuck=0;
        SinglePar=1;
      }
      else{
        if(Potential=="Dynamic"){
          if(CurrentPot=="DoubleGaussSum") CurrentPot="YukawaDimiCore";
          else CurrentPot="DoubleGaussSum";
        }
        if(!DEBUG){
          //cout<<flush;
          printf("\n\n----- RESET (%s) -----\n",CurrentPot.Data());
          printf("\r\033[K Progress: Break (%5u/%5u)",uIter,MaxIter);
        }
        else{
          printf(" --- RESET (%s) ---\n",CurrentPot.Data());
        }

        //make it the max value
        for(unsigned uPar=0;uPar<=Npars;uPar++)
          Convergence[uPar]=1;
        if(bdist<FallBack[0]){
          FallBack[0] = bdist;
          FallBack[1] = bf_0;
          FallBack[2] = bd_0;
          FallBack[3] = bV_1;
          FallBack[4] = bmu_1;
          FallBack[5] = bV_2;
          FallBack[6] = bmu_2;
          FallBack[7] = s_2;
        }
        VeryStuck=0;
        UnstuckCounter=0;
        Stuck=0;
        SinglePar=0;
        bdist = 1e16;
        f_0 = 0; d_0 = 0;
        bV_1 = 0; bV_2 = 0;
        bmu_1 = 1.5; bmu_2 = 0.5;
        bdf_0 = 1e6; bdd_0 = 1e6;
        bf_0 = f0; bd_0 = d0;
        for(unsigned uf=0; uf<fluctN; uf++) fluct[uf]=-1000*dist_unit;
        for(unsigned ud=0; ud<dirN; ud++) dir[ud]=true;
        Starting = true;
        Started = false;
      }
    }
    if(SinglePar){
      SinglePar = 1+UnstuckCounter/UnstuckPerPar;
      if(Stuck<UnstuckPerPar*Npars) printf(" stuck, but working on it...");
      else printf(" very stuck, Jesus is working on it...");
      desired_fluct = bdist/(2.*VeryStuck);
      //desired_fluct = 10.*dist_unit;//10fm... i.e. almost start over
    }
    else{
      desired_fluct = bdist;
    }
    if(DEBUG) printf("\n A: Stuck=%u(%u); UnstuckCounter=%u; SinglePar=%u",Stuck,VeryStuck,UnstuckCounter,SinglePar);
//desired_fluct = 1;

    //if(CurrentPot=="Yukawa"){
    //  Kitty_SE.SetEpsilonConv(5e-8);
    //  Kitty_SE.SetEpsilonProp(5e-8);
    //  Kitty_SE.SetMaxRad(128);
    //  Kitty_SE.SetMaxRho(64);
    //}

    double MaxConv;
    if(CurrentPot=="YukawaDimiCore") MaxConv = 4;
    else MaxConv = 2;
    if(avg_fluct<0||Starting) Convergence[SinglePar]=1;
    else if(avg_fluct/desired_fluct>1) Convergence[SinglePar]/=1.25;
    else if(Convergence[SinglePar]<0.8) Convergence[SinglePar]*=1.25;
    //not good for the potentials to get too scattered values
    else if(Convergence[SinglePar]<MaxConv) Convergence[SinglePar]*=1.05;
    else Convergence[SinglePar] = MaxConv;

    //this will make min:max value Convergence[SinglePar]:1
    //will go towards one if we have no entries getting towards the solution
    //apply only in normal mode, not when already stuck.
    if(Convergence[SinglePar]<1&&SinglePar==0){
      Convergence[SinglePar] *= (1.+(1./Convergence[SinglePar]-1)/(1.+exp((avg_dir-0.25/2.)/(0.25/8.))));
    }



//ADD NUM ITER WITHOUT IMPROVEMENT TO INCREASE THE Convergence
    if(DEBUG){
      printf("\n fluct = %.2f(%.2f)",avg_fluct,desired_fluct);
      printf("\n avg_dir = %.2f(%.2f)",avg_dir,(1.+(1./Convergence[SinglePar]-1)/(1.+exp((avg_dir-0.25/2.)/(0.25/8.)))));
      printf("\n Convergence = %.3f (%.3f)",Convergence[SinglePar],Convergence[0]);
      printf(" %.2e; %.2f; %.2e %.2e\n",bV_1,bmu_1,bV_2,bmu_2);
      //printf("\n bdf_0=%.2f; bdd_0=%.2f",bdf_0,bdd_0);
      //printf("\n bf_0=%.2f; bd_0=%.2f",bf_0,bd_0);
    }

    if(Starting&&!Started){
      if(!StartPars) OkayishStartingPars(CurrentPot,f0,d0,bV_1,bmu_1,bV_2,bmu_2,s_2);
      else {bV_1=StartPars[0];bmu_1=StartPars[1];bV_2=StartPars[2];bmu_2=StartPars[3];}
      V_1 = bV_1;
      mu_1 = bmu_1;
      V_2 = bV_2;
      mu_2 = bmu_2;
      Started = true;
    }
    else{
      if(CurrentPot=="YukawaDimiCore"){
        if(!SinglePar||SinglePar==1){
          do V_1 = rangen.Gaus(bV_1,0.005+5.*Convergence[1]*Convergence[0]*FineTune);
          while(V_1<0);
        }
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.*Convergence[2]*Convergence[0]*FineTune);}
          while(mu_1<0.||mu_1>3.0);
        }
        if(!SinglePar||SinglePar==3){
          V_2 = rangen.Gaus(bV_2,0.2+200.*Convergence[3]*Convergence[0]*FineTune);
          //do V_2 = rangen.Gaus(bV_2,0.2+200.*Convergence[3]*Convergence[0]);
          //while(V_2<0);
        }
        if(!SinglePar||SinglePar==4){
          do{mu_2 = rangen.Gaus(bmu_2,0.0005+0.5*Convergence[4]*Convergence[0]*FineTune);}
          while(mu_2<0.1||mu_2>1.5);
        }
      }
      else if(CurrentPot=="Gaussian"){
        if(!SinglePar||SinglePar==1){
          V_1 = rangen.Gaus(bV_1,0.05+40.*Convergence[1]*Convergence[0]*FineTune);
        }
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.0*Convergence[2]*Convergence[0]*FineTune);}//----
          while(mu_1<0.||mu_1>3.0);
        }
        V_2 = 0;
        mu_2 = 0;
      }
      else if(CurrentPot=="Yukawa"){
        if(!SinglePar||SinglePar==1){
          do V_1 = rangen.Gaus(bV_1,0.05+20.*Convergence[1]*Convergence[0]*FineTune);
          while(V_1<0);
        }
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.0*Convergence[2]*Convergence[0]*FineTune);}
          while(mu_1<0.||mu_1>3.0);
        }
        V_2 = 0;
        mu_2 = 0;
      }
      else{
        //this is the long range guy, we expect it to be negative for positive f0
        //do {
        if(!SinglePar||SinglePar==1){
          V_1 = rangen.Gaus(bV_1,0.2+500.*Convergence[1]*Convergence[0]*FineTune);
        }
        //}
        //while(f0>0&&V_1>0);
        if(!SinglePar||SinglePar==2){
          do{mu_1 = rangen.Gaus(bmu_1,0.001+1.*Convergence[2]*Convergence[0]*FineTune);}
          while(mu_1<0.||mu_1>2.5);
        }
        //the second gaussian should be short ranged and stronger in amplitude
        //do{
        if(!SinglePar||SinglePar==3){
          V_2 = rangen.Gaus(bV_2,0.2+500.*Convergence[3]*Convergence[0]*FineTune);
        }
        //}
        //while(fabs(V_2)<fabs(V_1));
        if(!SinglePar||SinglePar==4){
          do{mu_2 = rangen.Gaus(bmu_2,0.001+1.*Convergence[4]*Convergence[0]*FineTune);}
          while(mu_2<0.||mu_2>2.5);
        }
      }
    }


    TH1F* hDummy; TF1* fDummy;
    Kitty_SE.SetShortRangePotential(0,0,0,V_1);
    Kitty_SE.SetShortRangePotential(0,0,1,mu_1);
    Kitty_SE.SetShortRangePotential(0,0,2,V_2);
    Kitty_SE.SetShortRangePotential(0,0,3,mu_2);
    if(CurrentPot=="YukawaDimiCore"){
      if(StartPars&&StartPars[4]){
        s_2 = StartPars[4];
      }
      else s_2 = mu_2*0.2;
    }
    else if(CurrentPot=="Gaussian") s_2 = -1;
    else if(CurrentPot=="Yukawa") s_2 = -2;
    else s_2 = 0;
    Kitty_SE.SetShortRangePotential(0,0,4,s_2);
    //Kitty_SE.SetNotifications(CATS::nAll);
    //printf(" KillTheCat\n");
    Kitty_SE.SetAnaSource(0,Radii[0]);
    Kitty_SE.KillTheCat();

    if(DEBUG) printf("\n    V1=%.1f mu1=%.3f V2=%.1f mu2=%.3f",V_1,mu_1,V_2,mu_2);
    if(DEBUG) printf("\n b: V1=%.1f mu1=%.3f V2=%.1f mu2=%.3f",bV_1,bmu_1,bV_2,bmu_2);
//printf("kfgnslkjglfkgnjdkf\n");
    if(!DEBUG){
      printf("\n Current solution: V1=%.4e  mu1=%.4e  V2=%.4e  mu2=%.4e <--> f0=%.3f  d0=%.3f      ",
      bV_1,bmu_1,bV_2,bmu_2,bf_0,bd_0);
    }

    if(Eval_ScattParameters(Kitty_SE,f_0,d_0,hDummy,fDummy)) BadPhaseShifts = 0;
    else BadPhaseShifts++;
    //hDummy->Write();
    //fDummy->Write();

    if(hDummy) delete hDummy; if(fDummy) delete fDummy;
    if(BadPhaseShifts) {
      //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
      //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
      if(!DEBUG){
        //cout << "JULI" <<endl;
        cout << flush;
        cout << "\e[A";
        cout << flush;
      }
      else printf("\n ISSUE");
      continue;
    }
    df_0 = fabs(f_0-f0);
    dd_0 = fabs(d_0-d0);
    dist = sqrt(0.5*df_0*df_0/df0/df0+0.5*dd_0*dd_0/dd0/dd0);
    if(DEBUG) printf("\n    f_0=%.2f df_0=%.2e; d_0=%.2f dd_0=%.2e",f_0,df_0,d_0,dd_0);
    if(DEBUG) printf("\n b: f_0=%.2f df_0=%.2e; d_0=%.2f dd_0=%.2e",bf_0,bdf_0,bd_0,bdd_0);

    //we get a better result, that has the correct sign
    //GoodGoing = ( (df_0<bdf_0||df_0<df0)&&(dd_0<bdd_0||dd_0<dd0)
    //              &&(f_0*bf_0>0)&&(f_0*f0>0)&&(d_0*bd_0>0)&&(d_0*d0>0));
    //we reduce df_0 or dd_0 (obsolete?)(df_0<bdf_0||dd_0<bdd_0)&&
    //we have the correct sign
    //we reduce the distance
    //GoodGoing = (dist<bdist&&(f_0*bf_0>0)&&(f_0*f0>0)&&(d_0*bd_0>0)&&(d_0*d0>0));
    //the above line files for negative effective range
    GoodGoing = (dist<bdist&&(f_0*bf_0>0)&&(f_0*f0>0));
    //if(bdist<1e15) fluct[ufluct] = fabs(dist-bdist);
    //if(bdist<1e15){
//CHANGED TO BEST, SHOULD I???
      //fluct[ufluct] = sqrt(0.5*fabs(f_0-lf_0)*fabs(f_0-lf_0)/df0/df0+0.5*fabs(d_0-ld_0)*fabs(d_0-ld_0)/dd0/dd0);
      fluct[ufluct] = sqrt(0.5*fabs(f_0-bf_0)*fabs(f_0-bf_0)/df0/df0+0.5*fabs(d_0-bd_0)*fabs(d_0-bd_0)/dd0/dd0);
      //imagine 3 pts, B,N,G (best, new, goal), spanning vectors from BG=g and BN=n
      //than the vectors are given as (nx-bx,ny-by) etc., where x in my case is f0 and y is d0
      //than we compute the angle from the dot product rule. We want it |alpha|<45 deg to move in the right direction
      double nx=(f_0-bf_0); double ny=(d_0-bd_0);
      double gx=(f0-bf_0); double gy=(d0-bd_0);
      double CosAlpha;
      if((nx==0&&ny==0)||(gx==0&&gy==0)) CosAlpha = 1;
      else CosAlpha = (nx*gx+ny*gy)/(sqrt(nx*nx+ny*ny)*sqrt(gx*gx+gy*gy));
      //so true, if we fluctuated towards the soltion, without going more than dist in the other direction
//perhaps stuck mode, where we increase fluct, and change the best solution as long as we go to "the other side"
      dir[udir] = (CosAlpha>(sqrt(2.)/2.))&&(fluct[ufluct]<2.*bdist);

      if(DEBUG){
        printf("\n Current iter: GG=%i; CosAlpha=%.2f",GoodGoing,CosAlpha);
      }
    //}
    //else{
    //  fluct[ufluct] = dist;
    //}
    ufluct++;ufluct=ufluct%fluctN;
    udir++;udir=udir%dirN;

    if(GoodGoing){
      bf_0 = f_0;
      bdf_0 = df_0;
      bd_0 = d_0;
      bdd_0 = dd_0;
      bV_1 = V_1;
      bmu_1 = mu_1;
      bV_2 = V_2;
      bmu_2 = mu_2;
      //Distance = sqrt(bf_0*bf_0+bd_0*bd_0);
      bdist = sqrt(0.5*df_0*df_0/df0/df0+0.5*dd_0*dd_0/dd0/dd0);
      Stuck = 0;
      Starting = false;
    }
    lf_0 = f_0;
    ld_0 = d_0;

    FoundIt = (fabs(df0)>=df_0&&fabs(dd0)>=dd_0);
    for(unsigned uRad=0; uRad<NumRad; uRad++){
      if(uRad){
        Kitty_SE.SetAnaSource(0,Radii[uRad]);
        Kitty_SE.KillTheCat();
      }
      int Class = Get_MMclass(Radii[uRad],f_0,d_0);
      Fill_ntMM(ntMM,Kitty_SE,Radii[uRad],f_0,d_0,V_1,mu_1,V_2,mu_2,s_2,10+FoundIt);
      hClassCount->Fill(Class);
      if(uRad==0) hif0d0Count->Fill(1./f_0,d_0);
    }

    if(SinglePar){
      UnstuckCounter++;
    }
    else{
      uIter++;
    }
    if(bdf_0<0.999e6||bdf_0>1.001e6){
      Stuck++;
    }
    if(DEBUG){
      printf("\n FoundIt=(df0>=df_0&&dd0>=dd_0); %i=(%.3f>=%.3f && %.3f>=%.3f)\n-----------------------------",FoundIt,fabs(df0),df_0,fabs(dd0),dd_0);
      usleep(500e3);
    }
//printf("hi\n");
    //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
    //cout<<flush;cout<<"\033[F";cout<<flush;cout<<"\033[F";
    //cout<<flush;cout<<"\033[F";
    if(!DEBUG){
      cout << flush;
      cout << "\e[A";
      cout << flush;
    }

  }
  if(DEBUG){printf("\ndist: %.2f vs %.2f\n",bdist,FallBack[0]);}
  if(bdist>FallBack[0]){
    bdist = FallBack[0];
    bf_0 = FallBack[1];
    bd_0 = FallBack[2];
    bV_1 = FallBack[3];
    bmu_1 = FallBack[4];
    bV_2 = FallBack[5];
    bmu_2 = FallBack[6];
    s_2 = FallBack[7];
    if(fabs(s_2)<0.001) CurrentPot = "DoubleGaussSum";
    else CurrentPot = "YukawaDimiCore";
  }
  //TH1F* hClassCount = new TH1F("hClassCount","hClassCount",100,-50.5,49.5);
  //TH2F* hif0d0Count = new TH2F("hf0d0Count","hf0d0Count",20,-10,10,20,0,40);
  fOutput.cd();
  ntMM->Write();
  hClassCount->Write();
  hif0d0Count->Write();
  for(unsigned uRad=0; uRad<NumRad; uRad++){
    Kitty_SE.SetAnaSource(0,Radii[uRad]);
    Kitty_SE.KillTheCat();
    //printf("\n Cross check the scattering parameters");
//BINNING_ISSUE->JUST LOOK AT THE KCOTG IN BOTH PLOTS
    //double cc_f0=bf_0; double cc_d0=bd_0; TH1F* hCC; TF1* fCC;
//THE BINNING COMING OUT OF HERE
    //bool statusCC = Eval_ScattParameters(Kitty_SE,cc_f0,cc_d0,hCC,fCC);
    //printf("\n status=%i; f0=%.3f; d0=%.3f",statusCC,cc_f0,cc_d0);
    //hCC->Write();
    //fCC->Write();
    //delete hCC; delete fCC;
//IS DIFFERENT THAN OUT OF THERE
    Write_MM_Plots(&fOutput,Kitty_SE,Radii[uRad],bf_0,bd_0,TString::Format("_r%id%i",int(Radii[uRad]*100.)/100,int(Radii[uRad]*100.)%100));
  }

  printf("\n");
  printf("Suitable %s potential found:\n",CurrentPot.Data());
  printf(" f0 = %.3f fm\n", bf_0);
  printf(" d0 = %.2f fm\n", bd_0);
  printf("  V1 = %.6e\n",bV_1);
  printf("  mu1 = %.6e\n",bmu_1);
  printf("  V2 = %.6e\n",bV_2);
  printf("  mu2 = %.6e\n",bmu_2);

  V1 = bV_1;
  V2 = bV_2;
  mu1 = bmu_1;
  mu2 = bmu_2;

  delete hClassCount;
  delete hif0d0Count;
  delete ntMM;
  delete [] MomBins;
}

void Example_ProfS(int flag){
  printf("Make a potential %i\n",flag);
  double V1,V2,mu1,mu2;
  //ManufactureYukawaPotential(0.5,0.0025,1,0.005,V1,mu1,V2,mu2);
  //ManufactureYukawaPotential(0.5,0.0025,2,0.005*2,V1,mu1,V2,mu2);
  //ManufactureYukawaPotential(0.5,0.0025,4,0.005*4,V1,mu1,V2,mu2);
  //ManufactureYukawaPotential(0.5,0.0025,8,0.005*8,V1,mu1,V2,mu2);

//make the dist and fluct refer to d0 and f0 separately, check if the parametes get closer
//to their goal (not only distance, but direction...)

  unsigned NumR=16;
  double Radii[NumR];

  //NumR=8;
  //Radii[0] = 0.9;Radii[1] = 1.1;Radii[2] = 1.3;Radii[3] = 1.5;
  //Radii[4] = 2.0;Radii[5] = 2.5;Radii[6] = 3.0;Radii[7] = 4.0;

  NumR=1;
  Radii[0]=1.08;Radii[1]=1.08*2;Radii[2]=1.08*4;

  double f0,d0;
  double ef0 = 0.005;
  double ed0 = 0.005;
//NumR=1; f0=25.00; d0=25.00; ManufacturePotential(f0,f0*ef0,d0,d0*ed0,Radii,NumR,V1,mu1,V2,mu2,"Dynamic");
//return;

  //ADD-ONS FOR EACH FLAG
  TString OutputFolder = TString::Format("./ManufacturedPotentials");
  //flag = last digit is for the potential 0 - dynamic, 1 - DoubleGaussSum, 2 - YukawaDimiCore, 3 - Gaussian, 4 - Yukawa
  TString Potential = "Dynamic";
  if(flag%10==1) Potential = "DoubleGaussSum";
  else if(flag%10==2) Potential = "YukawaDimiCore";
  else if(flag%10==3) Potential = "Gaussian";
  else if(flag%10==4) Potential = "Yukawa";

//for Prof. S 50X, 51X....
// 50X,51X - doublet potentials
// 52X,53,54X - quarted potentials, if other -> perform all of them
// the 1XXX are for systematics
// e.g. for the flag 50X, systematics are 150X, 160X, 170X, 180X,
// for flag 51X we will have 151X, 161X, 171X, 181X etc.
  ef0 = 0.1/10.;
  ed0 = 0.05/5.;
  NumR = 1;
  Radii[0] = 1.2;
  //OutputFolder += "ProfS/";
  int VAR_FLAG = (flag/10)%10;
  if(flag/10>=150){
    VAR_FLAG += (flag/100)*10;//+150,160,170,180
  }
  //for(unsigned uVar=0; uVar<5; uVar++){
    //if(VAR_FLAG<5&&uVar!=VAR_FLAG) continue;
    switch (VAR_FLAG) {
      case 0: f0=-16.8; d0 = 2.3; break;
      case 150: f0=-21.2; d0 = 2.3; break;
      case 160: f0=-14.4; d0 = 2.3; break;
      case 170: f0=-16.8; d0 = 2.6; break;
      case 180: f0=-16.8; d0 = 2.0; break;
      case 1: f0=-16.3; d0 = 3.2; break;
      case 2: f0=7.6; d0 = 3.6; break;
      case 3: f0=10.8; d0 = 3.8; break;
      case 4: f0=17.3; d0 = 3.6; break;
      default: printf("Weird flags for producing the potentials for Prof. S.\n"); return;
    }
    double StartPars[5];
    double FT = 1;

    if(Potential=="DoubleGaussSum"){
      //this does not work, ones you have result for uVar==1, plug in those here
      if(VAR_FLAG==0){
        //StartPars[0]=-8.999757e+02;
        //StartPars[1]=5.081559e-01;
        //StartPars[2]= -4.021392e+02;
        //StartPars[3]=1.019864e+00;
        //StartPars[4]=0;
        //StartPars[0]=-1.943168e+02;
        //StartPars[1]=1.372106e+00;
        //StartPars[2]= 3.646892e+02;
        //StartPars[3]=9.833930e-01;
        //StartPars[4]=0;
        //StartPars[0]=-1.455808e+02;//V1 (MeV)
        //StartPars[1]=1.165081e+00;//mu1 (fm)
        //StartPars[2]=3.700546e+02;//V2 (MeV)
        //StartPars[3]=6.652901e-01;//mu2 (fm)
        StartPars[0]=1.666386e+02;//V1 (MeV)
        StartPars[1]=9.827926e-01;//mu1 (fm)
        StartPars[2]=-7.677560e+02;//V2 (MeV)
        StartPars[3]=9.310511e-01;//mu2 (fm)
        StartPars[4]=0;
      }
      else if(VAR_FLAG==150){
        StartPars[0]=-1.454633e+02;
        StartPars[1]=1.156062e+00;
        StartPars[2]=3.649164e+02;
        StartPars[3]=6.643602e-01;
        StartPars[4]=0;
        FT = 1./128.;
//        Make a potential 1501
//         Goal: DoubleGaussSum (f0,d0) [Achieved]: (-21.200,2.300)+/-(0.010,0.010) [ 38%,251%], Break (   97/10000) stuck, but working on it...
//        Suitable DoubleGaussSum potential found:.1559e+00  V2=3.6454e+02  mu2=6.6436e-01 <--> f0=-21.174  d0=2.304
//         f0 = -21.197 fm
//         d0 = 2.30 fm
//          V1 = -1.454633e+02
//          mu1 = 1.156062e+00
//          V2 = 3.649164e+02
//          mu2 = 6.643602e-01


      }
      else if(VAR_FLAG==160){
        StartPars[0]=-1.406252e+02;
        StartPars[1]=1.183071e+00;
        StartPars[2]=3.671337e+02;
        StartPars[3]=6.646221e-01;
        StartPars[4]=0;
        FT = 1./1024.;
        // Current solution: V1=-1.4051e+02  mu1=1.1833e+00  V2=3.6563e+02  mu2=6.6525e-01 <--> f0=-14.397  d0=2.289
      }
      else if(VAR_FLAG==170){
        StartPars[0]=-1.442826e+02;
        StartPars[1]=1.251641e+00;
        StartPars[2]=3.755765e+02;
        StartPars[3]=7.397668e-01;
        StartPars[4]=0;
        FT = 1./256.;
        //ef0 = 1.0;

  //      Make a potential 1701
  //       Goal: DoubleGaussSum (f0,d0) [Achieved]: (-16.800,2.600)+/-(0.010,0.010) [ 15%,139%], Break (   89/10000) stuck, but working on it...
  //      Suitable DoubleGaussSum potential found:.2535e+00  V2=3.7566e+02  mu2=7.4164e-01 <--> f0=-16.867  d0=2.593
  //       f0 = -16.801 fm
  //       d0 = 2.59 fm
  //        V1 = -1.446716e+02
  //        mu1 = 1.253497e+00
  //        V2 = 3.756614e+02
  //        mu2 = 7.415382e-01

      }
      else if(VAR_FLAG==180){
        StartPars[0]=-1.420678e+02;
        StartPars[1]=1.081553e+00;
        StartPars[2]=3.579042e+02;
        StartPars[3]=5.812257e-01;
        StartPars[4]=0;
        FT = 1./2048.;
        //ef0 = 1.0;
      }
      else if(VAR_FLAG==1){
        //StartPars[0]=-8.999757e+02;
        //StartPars[1]=5.081559e-01;
        //StartPars[2]= -4.021392e+02;
        //StartPars[3]=1.019864e+00;
        //StartPars[4]=0;
        StartPars[0]=-1.707025e+02;
        StartPars[1]=1.495874e+00;
        StartPars[2]= 3.626808e+02;
        StartPars[3]=9.711152e-01;
        StartPars[4]=0;
      }
      else if(VAR_FLAG==2){
        StartPars[0]=-1.008636e+01;
        StartPars[1]=1.348757e+00;
        StartPars[2]= -1.410447e+01;
        StartPars[3]=2.147635e+00;
        StartPars[4]=0;
      }
      else if(VAR_FLAG==3){
        StartPars[0]=-4.867921e+02;
        StartPars[1]=1.175196e+00;
        StartPars[2]= -8.110269e+02;
        StartPars[3]=1.169148e-01;
        StartPars[4]=0;
      }
      else{
        StartPars[0]=-4.940004e+02;
        StartPars[1]=1.177133e+00;
        StartPars[2]= -8.133919e+02;
        StartPars[3]=1.007169e-01;
        StartPars[4]=0;
      }
      //if(uVar==1) ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars);
      //else ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder);
      FT = 1./10.;
      ef0 = 0.5;
      ed0 = 0.1;
      ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,StartPars,
        (Mass_L*Mass_d)/(Mass_L+Mass_d),FT);

    }
    else{
      ManufacturePotential(f0,ef0,d0,ed0,Radii,NumR,V1,mu1,V2,mu2,Potential,OutputFolder,11,NULL,
        (Mass_L*Mass_d)/(Mass_L+Mass_d),FT);
    }
}


void ScatteringPars(){
  printf("hello\n");
  Example_ProfS(501);
}
