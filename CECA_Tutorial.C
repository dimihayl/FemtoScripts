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
#include "/home/dimihayl/Apps/CATS3/include/CECA.h"
#include "/home/dimihayl/Apps/CATS3/include/TREPNI.h"
#include "/home/dimihayl/Apps/CATS3/include/DLM_RootWrapper.h"
#include "/usr/local/include/gsl/gsl_sf_dawson.h"


void Test1(){

  const double Lambda_width = (0.0789)*pow(10,-15)*FmToNu;

  TREPNI Database(0);
  TreParticle* p_pim = Database.NewParticle("pi-");
  TreParticle* p_pi0 = Database.NewParticle();
  TreParticle* p_proton = Database.NewParticle();
  TreParticle* p_neutron = Database.NewParticle();
  TreParticle* p_Lambda = Database.NewParticle();
  //TreParticle p_pim(Database);
  //TreParticle p_pi0(Database);
  //TreParticle p_proton(Database);
  //TreParticle p_neutron(Database);
  //TreParticle p_Lambda(Database);


  //p_pim->SetName("pi-");
  Database.GetParticle("pi-")->SetMass(1.*Mass_pic);
  p_pim->SetMassLimit(1.*Mass_pic,1.*Mass_pic);
  p_pim->SetWidth(0);
  p_pim->SetWidthLimit(0,0);
  p_pim->SetAbundance(0);
  p_pim->SetAbundanceLimit(0,0);

  p_pi0->SetName("pi0");
  p_pi0->SetMass(1.*Mass_pi0);
  p_pi0->SetMassLimit(1.*Mass_pi0,1.*Mass_pi0);
  p_pi0->SetWidth(0);
  p_pi0->SetWidthLimit(0,0);
  p_pi0->SetAbundance(0);
  p_pi0->SetAbundanceLimit(0,0);

  p_proton->SetName("proton");
  p_proton->SetMass(Mass_p);
  p_proton->SetMassLimit(Mass_p,Mass_p);
  p_proton->SetWidth(0);
  p_proton->SetWidthLimit(0,0);
  p_proton->SetAbundance(0);
  p_proton->SetAbundanceLimit(0,0);

  p_neutron->SetName("neutron");
  p_neutron->SetMass(Mass_n);
  p_neutron->SetMassLimit(Mass_n,Mass_n);
  p_neutron->SetWidth(0);
  p_neutron->SetWidthLimit(0,0);
  p_neutron->SetAbundance(0);
  p_neutron->SetAbundanceLimit(0,0);

  p_Lambda->SetName("Lambda");
  p_Lambda->SetMass(Mass_L);
  p_Lambda->SetMassLimit(Mass_L,Mass_L);
  p_Lambda->SetWidth(Lambda_width);
  p_Lambda->SetWidthLimit(Lambda_width,Lambda_width);
  p_Lambda->SetAbundance(100);
  p_Lambda->SetAbundanceLimit(50,150);

  p_Lambda->NewDecay();
  p_Lambda->GetDecay(0)->AddDaughter(*p_pim);
  p_Lambda->GetDecay(0)->AddDaughter(*p_proton);
  p_Lambda->GetDecay(0)->SetBranching(67);
  //p_Lambda->GetDecay(0)->SetBranchingLimit(65, 65.5);

  p_Lambda->NewDecay();
  p_Lambda->GetDecay(1)->AddDaughter(*p_pi0);
  p_Lambda->GetDecay(1)->AddDaughter(*p_neutron);
  p_Lambda->GetDecay(1)->SetBranching(33);
  //p_Lambda->GetDecay(1)->SetBranchingLimit(33.6, 33.9);
  //
  p_pim->Print();
  p_proton->Print();
  p_Lambda->Print();

  printf("Maxi wants this: %s\n",p_Lambda->GetDecay(0)->GetName().c_str());


  Database.SetTotalYield(100);
  bool database_qa = Database.QA();
  printf("Database QA: %i\n",database_qa);

}

void Test2(){

  std::vector<std::string> ListOfParticles;
  ListOfParticles.push_back("Proton");
  //ListOfParticles.push_back("Proton");
  ListOfParticles.push_back("Lambda");

  printf("Running for %s-%s pairs\n",ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str());

  const double Disp = 1.0;

  const double FracProtonReso = 0.6422*1;
  const double FracLambdaReso = 0.6438*0;

  const double Tau_ProtonReso = 1.65;
  const double Tau_LambdaReso = 4.69;//4.69

  const double Mass_ProtonReso = 1362;
  const double Mass_LambdaReso = 1462;

  TREPNI Database(0);
  Database.SetSeed(11);
  std::vector<TreParticle*> ParticleList;
  ParticleList.push_back(Database.NewParticle("Proton"));
  ParticleList.push_back(Database.NewParticle("Lambda"));
  ParticleList.push_back(Database.NewParticle("Pion"));
  ParticleList.push_back(Database.NewParticle("ProtonReso"));
  ParticleList.push_back(Database.NewParticle("LambdaReso"));

  double ExpectedPP,ExpectedPR,ExpectedRP,ExpectedRR;
  ExpectedPP = (1.-FracProtonReso)*(1.-FracLambdaReso);
  ExpectedPR = (1.-FracProtonReso)*FracLambdaReso;
  ExpectedRP = FracProtonReso*(1.-FracLambdaReso);
  ExpectedRR = FracProtonReso*FracLambdaReso;

  for(TreParticle* prt : ParticleList){
    if(prt->GetName()=="Proton"){
      prt->SetMass(Mass_p);
      prt->SetAbundance((1.-FracProtonReso));
    }
    else if(prt->GetName()=="Lambda"){
      prt->SetMass(Mass_L);
      prt->SetAbundance((1.-FracLambdaReso));
      prt->SetDelayTau(30);
    }
    else if(prt->GetName()=="Pion"){
      prt->SetMass(Mass_pic);
      prt->SetAbundance(0);
    }
    else if(prt->GetName()=="ProtonReso"){
      prt->SetMass(Mass_ProtonReso);
      prt->SetAbundance(FracProtonReso);
      prt->SetWidth(hbarc/Tau_ProtonReso);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Proton"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }
    else if(prt->GetName()=="LambdaReso"){
      prt->SetMass(Mass_LambdaReso);
      prt->SetAbundance(FracLambdaReso);
      prt->SetWidth(hbarc/Tau_LambdaReso);

      prt->NewDecay();
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Lambda"));
      prt->GetDecay(0)->AddDaughter(*Database.GetParticle("Pion"));
      prt->GetDecay(0)->SetBranching(100);
    }

    prt->SetPtPz(prt->GetMass()*0.8,500);
    //prt->SetAcceptance_Eta(-EtaCut,EtaCut);
    //prt->SetAcceptance_pT(500,1e6);
  }


  CECA Ivana(Database,ListOfParticles);
  Ivana.SetDisplacement(Disp);
  Ivana.SetHadronization(0);
  Ivana.SetTau(0);
  Ivana.SetEventMult(2);
  Ivana.SetSourceDim(2);
  Ivana.SetThreadTimeout(30);
  Ivana.SetFemtoRegion(200);
  Ivana.SetPropagateMother(false);
  Ivana.GHETTO_EVENT = true;
  Ivana.GoBabyGo(0);

  TFile fOutput(TString::Format("Putka_%s_%s_d%.2f_fp%.2f_fl%.2f.root",ListOfParticles.at(0).c_str(),ListOfParticles.at(1).c_str(),Disp,FracProtonReso,FracLambdaReso),"recreate");
  Ivana.GhettoFemto_rstar->ComputeError();
  TH1F* h_rstar_Ceca = Convert_DlmHisto_TH1F(Ivana.GhettoFemto_rstar,"h_rstar_Ceca");
  h_rstar_Ceca->Scale(1./h_rstar_Ceca->Integral(),"width");
  fOutput.cd();
  h_rstar_Ceca->Write();

  //TF1* f_reff;
  double reff = Get_reff(h_rstar_Ceca);
  //f_reff->Write();
  printf("reff = %.3f fm\n",reff);

  Ivana.GhettoFemto_mT_rstar->ComputeError();
  TH2F* h_mT_rstar = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rstar,"h_mT_rstar");
  h_mT_rstar->Scale(1./h_mT_rstar->Integral(),"width");
  fOutput.cd();
  h_mT_rstar->Write();

  Ivana.GhettoFemto_mT_rcore->ComputeError();
  TH2F* h_mT_rcore = Convert_DlmHisto_TH2F(Ivana.GhettoFemto_mT_rcore,"h_mT_rcore");
  h_mT_rcore->Scale(1./h_mT_rcore->Integral(),"width");
  fOutput.cd();
  h_mT_rcore->Write();

  Ivana.Ghetto_kstar_rstar->ComputeError();
  TH2F* h_kstar_rstar = Convert_DlmHisto_TH2F(Ivana.Ghetto_kstar_rstar,"h_kstar_rstar");
  fOutput.cd();
  h_kstar_rstar->Write();



}

void CECA_Tutorial(){
  //Test1();
  Test2();
}
