

//what is the probability that out of NumBins a single one (or more)
//go beyond the Nsigma value
void ProbOfSingleFluct(int NumBins, double Nsigma){
  printf("Running with NumBins=%i; Nsigma=%.2f\n",NumBins,Nsigma);
  int NumIter = 0;
  int FoundIt = 0;
  int FindN = 32;
  //TH1F* hTest = new TH1F("hTest","hTest",NumBins,0,1);
  TRandom3 rangen(11);
  //for(unsigned uBin=0; uBin<NumBins; uBin++){
  //  hTest->SetBinContent(uBin+1,0);
  //  hTest->SetBinError(uBin+1,1);
  //}
  while(FoundIt<FindN){
    for(unsigned uBin=0; uBin<NumBins; uBin++){
      NumIter++;
      if(fabs(rangen.Gaus(0,1))>0) {FoundIt++;}
      //FoundIt = FoundIt>=FindN;
    }
    if(NumIter%1000000==0){printf(" %iM (%i)\n",NumIter/1000000,FoundIt);}
  }

  double pval;
  printf(" Probability: 1:%.0f\n",1./pval);
}

void SingleToDoubleSource(double r_sin){
  const unsigned NumIter = 1000000;

  TH1F* hSingle = new TH1F("hSingle","hSingle",1024,0,16);
  TH1F* hDouble = new TH1F("hDouble","hDouble",1024,0,16);

  TRandom3 rangen(11);
  double x1,y1,z1;
  double x2,y2,z2;
  double xdiff,ydiff,zdiff;
  double rel,rsin1,rsin2;
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    x1 = rangen.Gaus(0,r_sin);
    y1 = rangen.Gaus(0,r_sin);
    z1 = rangen.Gaus(0,r_sin);
    x2 = rangen.Gaus(0,r_sin);
    y2 = rangen.Gaus(0,r_sin);
    z2 = rangen.Gaus(0,r_sin);

    xdiff = x2-x1;
    ydiff = y2-y1;
    zdiff = z2-z1;
    rsin1 = sqrt(x1*x1+y1*y1+z1*z1);
    rsin2 = sqrt(x2*x2+y2*y2+z2*z2);
    rel = sqrt(xdiff*xdiff+ydiff*ydiff+zdiff*zdiff);

    hSingle->Fill(rsin1);
    hSingle->Fill(rsin2);
    hDouble->Fill(rel);
  }

  TFile fOutput("./Output/SingleToDoubleSource.root","recreate");
  hSingle->Write();
  hDouble->Write();
}


void ComparePoissonGauss(){
  double mean = 5;
  double stdv = 5;
  double effcounts = pow(stdv/mean,-2.);
  effcounts *= 1;

  TF1* fGauss = new TF1("fGauss","[0]*TMath::Gaus(x,[1],[2],1)",-30,30);
  fGauss->FixParameter(0,1);
  fGauss->FixParameter(1,mean);
  fGauss->FixParameter(2,stdv);

  TF1* fPoission = new TF1("fGauss","[0]*TMath::Poisson(x*[1],[2])",-30,30);
  fPoission->FixParameter(0,effcounts/mean);
  fPoission->FixParameter(1,effcounts/mean);
  fPoission->FixParameter(2,effcounts);

  TCanvas* can = new TCanvas("can", "can", 1);
  can->cd(0); can->SetCanvasSize(1920/2, 1080/2); can->SetMargin(0.15,0.05,0.2,0.05);//lrbt
  fPoission->Draw();
  fGauss->Draw("same");

  can->SaveAs("./Output/ComparePoissonGauss.png");
}


void RandomDistos(){
  //ProbOfSingleFluct(1,5);
  //ProbOfSingleFluct(2,5);
  //ProbOfSingleFluct(4,5);
  //ProbOfSingleFluct(8,5);
  //ProbOfSingleFluct(16,5);
  //ProbOfSingleFluct(32,5);
  //SingleToDoubleSource(1.0);
  ComparePoissonGauss();
}
