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
void CECA_Tutorial(){
  Test1();
}
