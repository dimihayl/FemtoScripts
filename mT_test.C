void mT_test(){
  double M1 = 938;
  double M2 = 2.*M1;
  double AvgM = (M1+M2)*0.5;
  double RedM = (M1*M2)/(M1+M2);
  TLorentzVector P1;
  P1.SetXYZM(1000,1000,1000,M1);
  TLorentzVector P2;
  P2.SetXYZM(2000,1000,1000,M2);
  printf("Particle1:\n");
  printf(" m = %.2f; pT = %.2f; mT = %.2f\n",P1.M(),P1.Pt(),sqrt(P1.M2()+P1.Pt()*P1.Pt()));
  printf("Particle2:\n");
  printf(" m = %.2f; pT = %.2f; mT = %.2f\n",P2.M(),P2.Pt(),sqrt(P2.M2()+P2.Pt()*P2.Pt()));
  TLorentzVector P12 = P1+P2;
  printf("SUM:\n");
  printf(" AvgM = %.2f; 2*RedM = %.2f; M/2 = %.2f\n",AvgM, 2.*RedM, P12.M()/2.);
  printf(" m = %.2f; pT = %.2f; mT = %.2f\n",P12.M(),P12.Pt(),sqrt(P12.M2()+P12.Pt()*P12.Pt()));
  printf(" mT per particle %.2f\n",0.5*sqrt(P12.M2()+P12.Pt()*P12.Pt()));
  printf(" basic mT = %.2f\n",sqrt(AvgM*AvgM+0.25*P12.Pt()*P12.Pt()));
}
