//a script to QA the calculations related to the mean lifetime
//the suspected issue: the effective Gamma = 1/tau needs to be evaluated instead of tau directly

void Lifetime_QA(){
  const hbarc = 197.327;
  const unsigned NumIter = 1000000;
  //each primary vector is a reso
  //each secondary vector contains info on the decay
  std::vector<float> Gamma;
  std::vector<std::vector<float>> Gamma_i;
  //std::vector<std::vector<float>> BrRat_i;

  Gamma.push_back(2*hbarc);//0.5fm/c
  std::vector<float> Gamma_0;
  //std::vector<float> BrRat_0;
  Gamma_0.push_back(Gamma.at(0)*0.5);


  Gamma.push_back(1*hbarc);//1fm/c
  Gamma.push_back(0.5*hbarc);//2fm/c
  Gamma.push_back(0.25*hbarc);//4fm/c


}
