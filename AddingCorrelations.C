

//we assume we have two independent measurements of the correlation function,
//each with a fixed value for the mixed events.
//We assume a fixed C_true
//we evaluate the expected number of SE
//we sample the SE following a Poisson and build up the Ck
//we do this iteratevly many times and see which method reproduces better the mean of C_true,
//as well as what the distribution around C_mean (C_stdv) is
//finally, we repeat all for different values of ME, i.e. to study the effect of getting more statistics

//Conclusions:
//N.B. Glen will NOT work if one of the samples has 0 entires, mine will (as long as at least 1 of them is non-zero)
void Dimi_vs_Glen(){
  const unsigned IterPerConfig = 1000;
  const unsigned IterMe = 100;
  const float C_true = 1.5;
  const float Me_Min = 1;
  const float Me_Max = 1000;

  TRandom3 rangen(11);
  TH2F* hDimi_Cmean = new TH2F("hDimi_Cmean","hDimi_Cmean",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);
  TH2F* hDimi_Cstdv = new TH2F("hDimi_Cstdv","hDimi_Cstdv",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);
  TH2F* hGlen_Cmean = new TH2F("hGlen_Cmean","hGlen_Cmean",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);
  TH2F* hGlen_Cstdv = new TH2F("hGlen_Cstdv","hGlen_Cstdv",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);
  TH2F* hDiff_w1 = new TH2F("hDiff_w1","hDiff_w1",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);
  TH2F* hDiff_w2 = new TH2F("hDiff_w2","hDiff_w2",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);

  TH2F* hDiff_new1 = new TH2F("hDiff_new1","hDiff_new1",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);
  TH2F* hDiff_new2 = new TH2F("hDiff_new2","hDiff_new2",IterMe,Me_Min,Me_Max,IterMe,Me_Min,Me_Max);

  for(unsigned uMe1=0; uMe1<IterMe; uMe1++){
    for(unsigned uMe2=uMe1; uMe2<IterMe; uMe2++){
      const int BinId1 = 1+uMe1;
      const int BinId2 = 1+uMe2;
      const float Me_Val1 = hDimi_Cmean->GetXaxis()->GetBinCenter(BinId1);
      const float Me_Val2 = hDimi_Cmean->GetXaxis()->GetBinCenter(BinId2);
      float DevMean_Dimi=0;
      float DevMean_Glen=0;
      float DevStdv_Dimi=0;
      float DevStdv_Glen=0;
      float DiffMean_w1=0;
      float DiffMean_w2=0;

      for(unsigned uConfig=0; uConfig<IterPerConfig; uConfig++){
        const float Se_Exp1 = C_true*Me_Val1;
        const float Se_Val1 = rangen.Poisson(Se_Exp1);
        const float Se_Err1 = sqrt(Se_Val1);
        const float Ck_Val1 = Se_Val1/Me_Val1;
        const float Ck_Err1 = Se_Err1/Me_Val1;


        const float Se_Exp2 = C_true*Me_Val2;
        const float Se_Val2 = rangen.Poisson(Se_Exp2);
        const float Se_Err2 = sqrt(Se_Val2);
        const float Ck_Val2 = Se_Val2/Me_Val2;
        const float Ck_Err2 = Se_Err2/Me_Val2;

        const float Ck_Dimi = (Se_Val1+Se_Val2)/(Me_Val1+Me_Val2);
        const float w1_Dimi = (Me_Val1)/(Me_Val1+Me_Val2);
        const float w2_Dimi = (Me_Val2)/(Me_Val1+Me_Val2);
        const float Cerr_Dimi = sqrt(Se_Val1+Se_Val2)/(Me_Val1+Me_Val2);
        DevMean_Dimi += (Ck_Dimi-C_true);
        DevStdv_Dimi += pow((Ck_Dimi-C_true),2.);

        const float w1_Glen = (1./Ck_Err1/Ck_Err1)/(1./Ck_Err1/Ck_Err1+1./Ck_Err2/Ck_Err2);
        const float w2_Glen = (1./Ck_Err2/Ck_Err2)/(1./Ck_Err1/Ck_Err1+1./Ck_Err2/Ck_Err2);
        const float Ck_Glen = w1_Glen*Ck_Val1+w2_Glen*Ck_Val2;
        const float Cerr_Glen = 1./sqrt(1./Ck_Err1/Ck_Err1+1./Ck_Err2/Ck_Err2);
        DevMean_Glen += (Ck_Glen-C_true);
        DevStdv_Glen += pow((Ck_Glen-C_true),2.);

        const float w1_New = (Ck_Val1/Ck_Err1/Ck_Err1)/(Ck_Val1/Ck_Err1/Ck_Err1+Ck_Val2/Ck_Err2/Ck_Err2);
        const float w2_New = (Ck_Val2/Ck_Err2/Ck_Err2)/(Ck_Val1/Ck_Err1/Ck_Err1+Ck_Val2/Ck_Err2/Ck_Err2);
        const float Ck_New = w1_New*Ck_Val1+w2_New*Ck_Val2;
        const float Cerr_New = sqrt(w1_New*w1_New*Ck_Err1*Ck_Err1+w2_New*w2_New*Ck_Err2*Ck_Err2);

        printf("Dimi vs Glen vs New\n");
        printf("w1: %.3f vs %.3f vs %.3f\n",w1_Dimi,w1_Glen,w1_New);
        printf("Ck: %.3f vs %.3f vs %.3f\n",Ck_Dimi,Ck_Glen,Ck_New);
        printf("er: %.3f vs %.3f vs %.3f\n",Cerr_Dimi,Cerr_Glen,Cerr_New);
        /*
        printf("  Rx = %f\n",Me_Val1);
        printf("  Ry = %f\n",Me_Val2);
        printf("  wx = %f\n",(Me_Val1)/(Me_Val1+Me_Val2));
        printf("  Sx = %f\n",Se_Val1);
        printf(" dSx = %f\n",sqrt(Se_Val1));
        printf("  Cx = %f\n",Se_Val1/Me_Val1;
        printf(" dCx = %f\n",sqrt(Se_Val1)/Me_Val1);

        printf("  rx = %f\n",1./Se_Val1);
        printf("  sx = %f\n",rx*);

        printf("  Sy = %f\n",Me_Val2);
        printf("  rx = %f\n",);
        printf("  sx = %f\n");
        printf("  wx = %f\n");
        */

        usleep(250e3);

        DiffMean_w1 += (w1_Glen-w1_Dimi);
        DiffMean_w2 += (w2_Glen-w2_Dimi);
      }
      DevMean_Dimi /= float(IterPerConfig);
      DevStdv_Dimi /= float(IterPerConfig);
      DevStdv_Dimi = sqrt(DevStdv_Dimi-DevMean_Dimi*DevMean_Dimi);

      DevMean_Glen /= float(IterPerConfig);
      DevStdv_Glen /= float(IterPerConfig);
      DevStdv_Glen = sqrt(DevStdv_Glen-DevMean_Glen*DevMean_Glen);

      DiffMean_w1 /= float(IterPerConfig);
      DiffMean_w2 /= float(IterPerConfig);

      hDimi_Cmean->SetBinContent(BinId1,BinId2,DevMean_Dimi);
      hDimi_Cstdv->SetBinContent(BinId1,BinId2,DevStdv_Dimi);

      hGlen_Cmean->SetBinContent(BinId1,BinId2,DevMean_Glen);
      hGlen_Cstdv->SetBinContent(BinId1,BinId2,DevStdv_Glen);

      hDiff_w1->SetBinContent(BinId1,BinId2,DiffMean_w1);
      hDiff_w2->SetBinContent(BinId1,BinId2,DiffMean_w2);
    }
  }
  TFile fOutput("./Output/Dimi_vs_Glen.root","recreate");
  hDimi_Cmean->Write();
  hDimi_Cstdv->Write();
  hGlen_Cmean->Write();
  hGlen_Cstdv->Write();
  hDiff_w1->Write();
  hDiff_w2->Write();
}

//0 =
//1 = one ME linear, another uniform
void Dani_Suggestion(int type){

  unsigned Yield1;
  unsigned Yield2;
  unsigned NumIter;
  double kMin;
  double kMax;
  unsigned NumMomBins;
  double k_mean1;

  TF1* fCtrue;

  if(type==0){
    Yield1 = 200*1000;
    Yield2 = 50*1000;
    NumIter = (Yield1+Yield2)*64;
    kMin = 0;
    kMax = 800;
    NumMomBins = 160;
    k_mean1 = 200;
    fCtrue = new TF1("fCtrue","1.+[0]*TMath::Gaus(x,0,[1],false)",kMin,kMax);
    fCtrue->FixParameter(0,1);
    fCtrue->FixParameter(1,60);
  }
  else{
    kMin = 0;
    kMax = 64;
    NumMomBins = 64;
    Yield2 = NumMomBins*64;//the uniform
    Yield1 = Yield2/2;

    fCtrue = new TF1("fCtrue","1.+[0]*TMath::Gaus(x,0,[1],false)",kMin,kMax);
    fCtrue->FixParameter(0,-0.5);
    fCtrue->FixParameter(1,20);
  }

  //const double k_mean2 = 360;

  TRandom rnd(11);

  double Me_kstar;
  TH1F* hMe1 = new TH1F("hMe1","hMe1",NumMomBins,kMin,kMax);
  TH1F* hMe2 = new TH1F("hMe2","hMe2",NumMomBins,kMin,kMax);

  TH1F* hSe1 = new TH1F("hSe1","hSe1",NumMomBins,kMin,kMax);
  TH1F* hSe2 = new TH1F("hSe2","hSe2",NumMomBins,kMin,kMax);

  if(type==0){
    for(unsigned uIter=0; uIter<NumIter; uIter++){
      Me_kstar = sqrt(pow(rnd.Gaus(0,k_mean1),2)+pow(rnd.Gaus(0,k_mean1),2)+pow(rnd.Gaus(0,k_mean1),2));
      hMe1->Fill(Me_kstar);

      Me_kstar = kMin+rnd.Uniform()*(kMax-kMin);
      //Me_kstar = sqrt(pow(rnd.Gaus(0,k_mean1),2)+pow(rnd.Gaus(0,k_mean1),2)+pow(rnd.Gaus(0,k_mean1),2));
      hMe2->Fill(Me_kstar);
    }
  }
  else{
    for(unsigned uMom=0; uMom<NumMomBins; uMom++){
      Me_kstar = hMe1->GetBinCenter(uMom+1);
      hMe1->SetBinContent(uMom+1,Me_kstar+hMe1->GetBinCenter(1));
      hMe1->SetBinError(uMom+1,sqrt(hMe1->GetBinContent(uMom+1)));

      hMe2->SetBinContent(uMom+1,Yield2);
      hMe2->SetBinError(uMom+1,sqrt(hMe2->GetBinContent(uMom+1)));
    }
  }

  hMe1->Scale(Yield1/hMe1->Integral());
  hMe2->Scale(Yield2/hMe2->Integral());

  double SeExp;
  double Mom;

  double SeVal,SeVal1,SeVal2;
  double MeVal,MeVal1,MeVal2;
  double CkVal,CkVal1,CkVal2;
  double CkErr,CkErr1,CkErr2;
  double w1,w2;
  double Ck1,Ck2;

  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    Mom = hMe1->GetBinCenter(uMom+1);

    SeExp = hMe1->GetBinContent(uMom+1)*fCtrue->Eval(Mom);
    SeVal = rnd.Poisson(SeExp);
    hSe1->SetBinContent(uMom+1,SeVal);

    SeExp = hMe2->GetBinContent(uMom+1)*fCtrue->Eval(Mom);
    SeVal = rnd.Poisson(SeExp);
    hSe2->SetBinContent(uMom+1,SeVal);
  }
  hSe1->Sumw2();
  hSe2->Sumw2();


  TH1F* hCk1 = (TH1F*)hSe1->Clone("hCk1");
  hCk1->Divide(hMe1);

  TH1F* hCk2 = (TH1F*)hSe2->Clone("hCk2");
  hCk2->Divide(hMe2);

  hCk1->Fit(fCtrue,"Q,S, N, R, M");
  printf("chi2_1 = %.2f\n",fCtrue->GetChisquare()/double(NumMomBins));
  hCk2->Fit(fCtrue,"Q,S, N, R, M");
  printf("chi2_2 = %.2f\n",fCtrue->GetChisquare()/double(NumMomBins));


  TH1F* hCkGlen = new TH1F("hCkGlen","hCkGlen",NumMomBins,kMin,kMax);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    CkErr1 = hCk1->GetBinError(uMom+1);
    CkErr2 = hCk2->GetBinError(uMom+1);
    CkVal1 = hCk1->GetBinContent(uMom+1);
    CkVal2 = hCk2->GetBinContent(uMom+1);
    w1 = (1./CkErr1/CkErr1)/(1./CkErr1/CkErr1+1./CkErr2/CkErr2);
    w2 = (1./CkErr2/CkErr2)/(1./CkErr1/CkErr1+1./CkErr2/CkErr2);
    CkVal = w1*CkVal1+w2*CkVal2;
    CkErr = 1./sqrt(1./CkErr1/CkErr1+1./CkErr2/CkErr2);
    hCkGlen->SetBinContent(uMom+1,CkVal);
    hCkGlen->SetBinError(uMom+1,CkErr);
  }
  hCkGlen->Fit(fCtrue,"Q,S, N, R, M");
  printf("chi2_glen = %.2f\n",fCtrue->GetChisquare()/double(NumMomBins));

  TH1F* hCkDimi = new TH1F("hCkDimi","hCkDimi",NumMomBins,kMin,kMax);
  for(unsigned uMom=0; uMom<NumMomBins; uMom++){
    SeVal1 = hSe1->GetBinContent(uMom+1);
    SeVal2 = hSe2->GetBinContent(uMom+1);
    MeVal1 = hMe1->GetBinContent(uMom+1);
    MeVal2 = hMe2->GetBinContent(uMom+1);
    CkVal = (SeVal1+SeVal2)/(MeVal1+MeVal2);
    CkErr = sqrt(SeVal1+SeVal2)/(MeVal1+MeVal2);

    hCkDimi->SetBinContent(uMom+1,CkVal);
    hCkDimi->SetBinError(uMom+1,CkErr);
  }
  hCkDimi->Fit(fCtrue,"Q,S, N, R, M");
  printf("chi2_dimi = %.2f\n",fCtrue->GetChisquare()/double(NumMomBins));

  TFile fOutput(TString::Format("./Output/Dani_Suggestion_%i.root",type),"recreate");
  hSe1->Write();
  hMe1->Write();
  hCk1->Write();
  hSe2->Write();
  hMe2->Write();
  hCk2->Write();
  fCtrue->Write();
  hCkGlen->Write();
  hCkDimi->Write();


}

void PrintW1W2(unsigned S1, unsigned R1, unsigned S2, unsigned R2){
  //printf("\n");
  double S_1 = S1;
  double S_2 = S2;
  double R_1 = R1;
  double R_2 = R2;
  double C_1 = S_1/R_1;
  double Ce_1 = sqrt(S_1)/R_1+1e-16;
  double C_2 = S_2/R_2;
  double Ce_2 = sqrt(S_2)/R_2+1e-16;
  printf("C_1 = S_1/R_1 = %u/%u = %.2f +/- %.2f\n",S1,R1,C_1,Ce_1);
  printf("C_2 = S_2/R_2 = %u/%u = %.2f +/- %.2f\n",S2,R2,C_2,Ce_2);

  double wGF_1 = (1./(Ce_1)/Ce_1)/(1./Ce_1/Ce_1+1./Ce_2/Ce_2);
  double wGF_2 = (1./Ce_2/Ce_2)/(1./Ce_1/Ce_1+1./Ce_2/Ce_2);
  //if(Ce_1=0&&Ce_2){wGF_1=1; wGF_2}
  printf("In GentleFemto:\n");
  printf(" w1 = %.0f%%, w2 = %.0f%%\n",wGF_1*100.,wGF_2*100.);
  printf(" C = %.2f +/- %.2f\n",wGF_1*C_1+wGF_2*C_2,1./sqrt(1./Ce_1/Ce_1+1./Ce_2/Ce_2));

  double wD_1 = (R_1)/(R1+R_2);
  double wD_2 = (R_2)/(R1+R_2);;
  printf("Dimi:\n");
  printf(" w1 = %.0f%%, w2 = %.0f%%\n",wD_1*100.,wD_2*100.);
  printf(" C = %.2f +/- %.2f\n",wD_1*C_1+wD_2*C_2,sqrt(S_1+S_2)/(R_1+R_2));
}

//R1 = 12 fixed
//R2 = 3 fixed
//C_true = 1
void VerySimpleExample(unsigned R1, unsigned R2, double Ctrue=1, unsigned NumIter=1000000){
  double R_1 = R1;
  double R_2 = R2;
  double S_1,S_2;
  double C_1,C_2,C_tot;
  double Ce_1,Ce_2,Ce_tot;
  double w_1,w_2;
  TH1F* hCk_Gentle = new TH1F("hCk_Gentle","hCk_Gentle",128,0,Ctrue*2.);
  TH1F* hCk_Dimi = new TH1F("hCk_Dimi","hCk_Dimi",128,0,Ctrue*2.);

  TH1F* hCk_GentleSD = new TH1F("hCk_GentleSD","hCk_GentleSD",128,0,Ctrue*2.);
  TH1F* hCk_DimiSD = new TH1F("hCk_DimiSD","hCk_DimiSD",128,0,Ctrue*2.);

  TH1F* hCk_GentleBD = new TH1F("hCk_GentleBD","hCk_GentleBD",128,0,Ctrue*2.);
  TH1F* hCk_DimiBD = new TH1F("hCk_DimiBD","hCk_DimiBD",128,0,Ctrue*2.);

  TRandom3 rnd(11);
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    S_1 = rnd.Poisson(R1)*Ctrue;
    S_2 = rnd.Poisson(R2)*Ctrue;
    C_1 = S_1/R_1;
    C_2 = S_2/R_2;
    Ce_1 = sqrt(S_1)/R_1+1e-16;
    Ce_2 = sqrt(S_2)/R_2+1e-16;
    //bool BigDiff = fabs(C_2-C_1)>sqrt(Ce_1*Ce_1+Ce_2*Ce_2);
    bool BigDiff = fabs(C_2-C_1)>(Ce_1+Ce_2);

    //GF
    w_1 = (1./(Ce_1)/Ce_1)/(1./Ce_1/Ce_1+1./Ce_2/Ce_2);
    w_2 = (1./Ce_2/Ce_2)/(1./Ce_1/Ce_1+1./Ce_2/Ce_2);
    C_tot = w_1*C_1+w_2*C_2;
    Ce_tot = 1./sqrt(1./Ce_1/Ce_1+1./Ce_2/Ce_2);
    hCk_Gentle->Fill(C_tot);
    if(BigDiff) hCk_GentleBD->Fill(C_tot);
    else hCk_GentleSD->Fill(C_tot);

    //Dimi
    w_1 = (R_1)/(R_1+R_2);
    w_2 = (R_2)/(R_1+R_2);
    C_tot = (S_1+S_2)/(R_1+R_2);
    Ce_tot = sqrt(S_1+S_2)/(R_1+R_2);
    hCk_Dimi->Fill(C_tot);
    if(BigDiff) hCk_DimiBD->Fill(C_tot);
    else hCk_DimiSD->Fill(C_tot);
  }

  TFile fOutput(TString::Format("./Output/VerySimpleExample_%u_%u_%.2f.root",R1,R2,Ctrue),"recreate");
  hCk_Gentle->Write();
  hCk_GentleSD->Write();
  hCk_GentleBD->Write();

  hCk_Dimi->Write();
  hCk_DimiSD->Write();
  hCk_DimiBD->Write();
}


void AddingCorrelations(){
  //Dimi_vs_Glen();
  //Dani_Suggestion(0);
  //Dani_Suggestion(1);
  //VerySimpleExample(12,3);

  //VerySimpleExample(16,12,0.6,10*1000*1000);
  PrintW1W2(14,16,4,12);

/*
  //printf("\nIn all cases, the C_true = 1\n");
  printf("\nExtreme case, one correlation is 0\n");
  //PrintW1W2(0,1,10,9);
  PrintW1W2(0,3,16,12);
  printf("\nAgain extreme, a single entry in one of the same events\n");
  //PrintW1W2(1,3,10,9);
  PrintW1W2(1,3,16,12);
  printf("\n20-30%% difference in yield, c.a. 20 entries\n");
  PrintW1W2(20,19,14,15);
  //printf("20-30% difference in yield, c.a. 100 entries\n");
  //PrintW1W2(100,95,70,75);
  */
}
