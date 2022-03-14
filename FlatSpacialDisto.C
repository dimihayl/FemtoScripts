
void FlatSpacialDisto(){
  //sample xyz
  //sample pT eta phi
  const int Mode = 5;
  const unsigned NumIter = 1*1000*1000;
  const double Sigma = 1;
  const double Pi = TMath::Pi();
  TH1F* hX = new TH1F("hX","hX",1024,-5.*Sigma,5.*Sigma);
  TH1F* hY = new TH1F("hY","hY",1024,-5.*Sigma,5.*Sigma);
  TH1F* hZ = new TH1F("hZ","hZ",1024,-5.*Sigma,5.*Sigma);
  TH1F* hPt = new TH1F("hPt","hPt",1024,0,10.*Sigma);
  TH1F* hPtot = new TH1F("hPtot","hPtot",1024,0,10.*Sigma);
  TH1F* hEta = new TH1F("hEta","hEta",1024,-8,8);
  TH1F* hPhi = new TH1F("hPhi","hPhi",1024,-2.*Pi,2.*Pi);
  TH1F* hTheta = new TH1F("hTheta","hTheta",1024,-2.*Pi,2.*Pi);
  TH1F* hCosTheta = new TH1F("hCosTheta","hCosTheta",1024,-1,1);
  TH1F* hSinTheta = new TH1F("hSinTheta","hSinTheta",1024,-1,1);
  double px,py,pz,ptot,pt,eta,phi,theta,cos_th,sin_th,tan_th,cotg_th;
  TRandom3 rangen(11);
  for(unsigned uIter=0; uIter<NumIter; uIter++){
    if(Mode==0){
      //x,y,z from Gauss, default choice
      px = rangen.Gaus(0,Sigma);
      py = rangen.Gaus(0,Sigma);
      pz = rangen.Gaus(0,Sigma);
      //px = rangen.Uniform(-Sigma,Sigma);
      //py = rangen.Uniform(-Sigma,Sigma);
      //pz = rangen.Uniform(-Sigma,Sigma);

      ptot = sqrt(px*px+py*py+pz*pz);
      pt = sqrt(px*px+py*py);
      cos_th = pz/ptot;
      //tan_th = pt/pz;
      //sin_th = pt/ptot;
      //theta = atan(tan_th);
      theta = acos(cos_th);
      sin_th = sin(theta);
      tan_th = tan(theta);
      //printf("QA: pt=%.2f pz*tan_th=%.2f\n",pt,pz*tan_th);
      //usleep(500e3);

      phi = atan2(py,px)+Pi;
      eta = -log(tan(0.5*theta));
      double eta_2 = -0.5*log((1.-cos_th)/(1.+cos_th));
      if( fabs(eta/eta_2 - 1)>1e-9){
        printf("FUCK eta: %f %f\n",eta,eta_2);
      }
      //printf("eta: %f %f\n",eta,eta_2);


    }
    else if(Mode==1){
      //sample pT, cos_th, phi => not good
      cos_th = rangen.Uniform(-1,1);
      sin_th = sqrt(1.-cos_th*cos_th)+1e-128;
      tan_th = sin_th/cos_th;
      phi = rangen.Uniform(0,2.*Pi);
      //ptot = sqrt(pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.));
      pt = sqrt(pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.));
      //pt = ptot*sin_th;

      pz = pt/tan_th;
      ptot = sqrt(pt*pt+pz*pz);
      px = ptot*cos(phi)*sin_th;
      py = ptot*sin(phi)*sin_th;
      eta = -0.5*log((1.-cos_th)/(1.+cos_th));
      theta = acos(cos_th);
    }
    //WORKS!!!
    else if(Mode==2){
      //sample pt,pz,phi => works
      pt = sqrt(pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.));
      pz = rangen.Gaus(0,Sigma);
      tan_th = pt/pz;
      theta = atan(tan_th);
      if(theta<0) theta += Pi;
      cos_th = cos(theta);
      sin_th = sin(theta);
      ptot = pz/cos_th;
      eta = -0.5*log((1.-cos_th)/(1.+cos_th));
      phi = rangen.Uniform(0,2.*Pi);
      px = ptot*cos(phi)*sin_th;
      py = ptot*sin(phi)*sin_th;
    }
    //WORKS!!!
    else if(Mode==3){
      //sample ptot,cos_th,phi => works
      cos_th = rangen.Uniform(-1,1);
      sin_th = sqrt(1.-cos_th*cos_th)+1e-128;
      tan_th = sin_th/cos_th;
      phi = rangen.Uniform(0,2.*Pi);
      ptot = sqrt(pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.));
      //pt = sqrt(pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.));
      pt = ptot*sin_th;
      pz = pt/tan_th;
      //ptot = sqrt(pt*pt+pz*pz);
      px = ptot*cos(phi)*sin_th;
      py = ptot*sin(phi)*sin_th;
      eta = -0.5*log((1.-cos_th)/(1.+cos_th));
      theta = acos(cos_th);
    }
    else if(Mode==4){
      //pt,eta,phi => does not work
      pt = sqrt(pow(rangen.Gaus(0,Sigma),2.)+pow(rangen.Gaus(0,Sigma),2.));
      eta = rangen.Uniform(-1,1);
      //do eta = rangen.Gaus(0,2);
      //while(fabs(eta)>1.5);
      phi = rangen.Uniform(0,2.*Pi);
      sin_th = 2.*exp(-eta)/(1.+exp(-2.*eta));
      cotg_th = (1.-exp(-2.*eta))/(2.*exp(-eta));
      cos_th = (1.-exp(-2.*eta))/(1.+exp(-2.*eta));
      pz = pt*cotg_th;
      ptot = sqrt(pt*pt+pz*pz);
      px = ptot*cos(phi)*sin_th;
      py = ptot*sin(phi)*sin_th;
      theta = acos(cos_th);
    }
    //WORKS
    else if(Mode==5){
      //sample px,py,pz, but only evaluate pT and eta
      //save the pT eta, sample random phi => we restore all correctly
      //this mean, if we have a 2D pT vs eta, we can sample phi randomly
      px = rangen.Gaus(0,Sigma);
      py = rangen.Gaus(0,Sigma);
      pz = rangen.Gaus(0,Sigma);
      ptot = sqrt(px*px+py*py+pz*pz);
      pt = sqrt(px*px+py*py);
      cos_th = pz/ptot;
      theta = acos(cos_th);
      phi = atan2(py,px)+Pi;
      eta = -log(tan(0.5*theta));

      //now we reset the true solution, keeping the info on pT and eta
      //phi is resampled
      phi = rangen.Uniform(0,2.*Pi);
      sin_th = 2.*exp(-eta)/(1.+exp(-2.*eta));
      cotg_th = (1.-exp(-2.*eta))/(2.*exp(-eta));
      cos_th = (1.-exp(-2.*eta))/(1.+exp(-2.*eta));
      pz = pt*cotg_th;
      ptot = sqrt(pt*pt+pz*pz);
      px = ptot*cos(phi)*sin_th;
      py = ptot*sin(phi)*sin_th;
      theta = acos(cos_th);
    }

    hX->Fill(px);
    hY->Fill(py);
    hZ->Fill(pz);
    hPt->Fill(pt);
    hEta->Fill(eta);
    hPhi->Fill(phi);
    hTheta->Fill(theta);
    hCosTheta->Fill(cos_th);
    hSinTheta->Fill(sin_th);
    hPtot->Fill(ptot);
  }

  TFile fOutput(TString::Format("FlatSpacialDisto_%i.root",Mode),"recreate");
  hX->Write();
  hY->Write();
  hZ->Write();
  hPtot->Write();
  hPt->Write();
  hEta->Write();
  hPhi->Write();
  hTheta->Write();
  hCosTheta->Write();
  hSinTheta->Write();
/*
  TCanvas c1; hX->Draw();
  TCanvas c2; hY->Draw();
  TCanvas c3; hZ->Draw();
  TCanvas c4; hPt->Draw();
  TCanvas c5; hEta->Draw();
  TCanvas c6; hPhi->Draw();
  TCanvas c7; hTheta->Draw();
  TCanvas c8; hCosTheta->Draw();
  TCanvas c9; hSinTheta->Draw();
*/


}
