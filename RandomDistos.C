

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



void RandomDistos(){
  ProbOfSingleFluct(1,5);
  ProbOfSingleFluct(2,5);
  ProbOfSingleFluct(4,5);
  ProbOfSingleFluct(8,5);
  ProbOfSingleFluct(16,5);
  ProbOfSingleFluct(32,5);
}
