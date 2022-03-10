
void ErrorPvalRafa(){
  const double Cth = 1;//2.8,3.01
  const double Cth_error = Cth*0.0;
  const unsigned Nse = 5;//5
  //normalied !!!
  const double Nme = double(Nse)+2.5*sqrt(double(Nse));
  //const double Nme = 3.145+4.77;//3.145+4.77
  //const double Nme = 7.91;
  const double Nme_error = Nme*0.0;

  double Cexp = double(Nse)/Nme;
  double Cexp_error = sqrt(double(Nse))/Nme;

  unsigned NumIter = 0;

  TRandom3 rangen(11);

  unsigned NumOutside = 0;
  while(NumOutside<100000){
    double C_rnd = rangen.Gaus(Cth,Cth_error);
    double Nme_rnd = rangen.Gaus(Nme,Nme_error);
    //expected Nse is C_rnd*Nme_rnd
    double Nse_exp = C_rnd*double(Nme_rnd);
    unsigned Nse_rnd = rangen.Poisson(Nse_exp);
    if(double(Nse)<Nse_exp){
      if(Nse_rnd<Nse) NumOutside++;
    }
    else if(double(Nse)>Nse_exp){
    //else{
      if(Nse_rnd>Nse) NumOutside++;
    }
    NumIter++;
  }

  double pval = 2.*double(NumOutside)/double(NumIter);
  double nsig = sqrt(2)*TMath::ErfcInverse(pval);

  double nsig_gaus = sqrt(pow((Cexp-Cth)/sqrt(pow(Cexp_error,2.)+pow(Cth_error,2.)),2.));

  printf(" pval = %.5f\n",pval);
  printf(" nsig = %.3f (%.3f)\n",nsig,nsig_gaus);
}
