#include <iostream> 
#include <sstream> 
#include <TRandom3.h>;
using namespace std; 

TGraphErrors* gppstat;

TGraphErrors* gppsyst;




//to make the figure:
Bool_t DRAWFIG = kTRUE;


void fitmTnew(){

TFile *_file0 = TFile::Open("RadppvsmT_0.root");
gppstat = (TGraphErrors*) _file0->Get("mTRadiusStat")->Clone("gppstat");gppstat->SetName("gppstat");
gppsyst = (TGraphErrors*) _file0->Get("mTRadiusSyst")->Clone("gppsyst");gppsyst->SetName("gppsyst");


//pp:
 const int  npp = 7;
 Float_t ppmT[npp],pprcore[npp],ppstat[npp],ppsyst[npp];
 Float_t ppex[npp] = {0,0,0,0,0,0,0};

 //fill pprcore with Y, ppstat and ppsyst with the corresponding stat and syst error:
 double xx,yy,eyy;
 for(int ii=0;ii<npp;ii++){
  gppstat->GetPoint(ii,xx,yy);
  eyy = gppstat->GetErrorY(ii);
  pprcore[ii] = yy;
  ppmT[ii] = xx;
  ppstat[ii] = eyy;
  eyy = gppsyst->GetErrorY(ii);
  ppsyst[ii] = eyy;
 }


 //draw them
 gppstat->Draw("");
 gppsyst->Draw("same");


//define some things
TRandom3 r3;
d21();
c1();
const int ntries=1000;
TString sfunction = "[0]*pow(x,[1])+[2]";
Float_t rdraw = .99; //draw only 10% of them

//first try with only proton:
//--------------------------------
//--------------------------------
TGraphErrors* gpp[ntries];
TF1* f1pp[ntries];
gppstat->Draw("ap");
Float_t pointspp[npp][ntries];
Float_t meanpp[npp];
Float_t sigmapp[npp];
Float_t uppp[npp];
Float_t downpp[npp];
for(int ip=0;ip<npp;ip++){
 meanpp[ip]=0.;
 sigmapp[ip]=0.;
}

//do the fits
for(int itry = 0;itry<ntries;itry++){
 Float_t esystpp[npp];
 for(int ip=0;ip<npp;ip++){
  esystpp[ip] = pprcore[ip]+r3.Gaus(0.,ppsyst[ip]);
  //cout<<" itry "<<itry<<" ip "<<ip<<" pprcore[ip]="<<pprcore[ip]<<" ppsyst[i] "<<ppsyst[ip]<<" esystpp[ip] "<<esystpp[ip]<<endl;
 }
 gpp[itry] = new TGraphErrors(npp,ppmT,esystpp,ppex,ppstat);
 //if(r3.Rndm()<rdraw)gp[itry]->Draw("same");
 f1pp[itry] = new TF1(Form("f1pp_%i",itry),sfunction,0.1,100.);
 gpp[itry]->Fit(Form("f1pp_%i",itry),"Q");
 if(r3.Rndm()<rdraw)f1pp[itry]->Draw("same");
 //get the points at every point 
 for(int ip=0;ip<npp;ip++){
  pointspp[ip][itry] = f1pp[itry]->Eval(ppmT[ip]);
  meanpp[ip] = meanpp[ip] + f1pp[itry]->Eval(ppmT[ip]);
 }
}

//create gmean with the mean and gup/gdown with mean+-std deviation
for(int itry = 0;itry<ntries;itry++){
 for(int ip=0;ip<npp;ip++){
  if(itry==0) meanpp[ip] = meanpp[ip] / ntries;
  sigmapp[ip] = sigmapp[ip] +  pow(pointspp[ip][itry] - meanpp[ip] ,2);
 }
}

  
//how many sigmas?:
Float_t nsigmas = 3.;
for(int ip=0;ip<npp;ip++){
 sigmapp[ip] = sqrt(sigmapp[ip]/(ntries-1));
 uppp[ip] = meanpp[ip] + nsigmas * sigmapp[ip];
 downpp[ip] = meanpp[ip] - nsigmas * sigmapp[ip];
}

TGraph* gmeanpp = new TGraphErrors(npp,ppmT,meanpp);
TGraph* guppp = new TGraphErrors(npp,ppmT,uppp);
TGraph* gdownpp = new TGraphErrors(npp,ppmT,downpp);

  
//and fit them!:
TF1* fmeanpp = new TF1("fmeanpp",sfunction,0.1,100.);
TF1* fuppp = new TF1("fuppp",sfunction,0.1,100.);
TF1* fdownpp = new TF1("fdownpp",sfunction,0.1,100.);

//gmeanpp->Fit(fmeanpp,"Q");
//guppp->Fit(fuppp,"Q");
//gdownpp->Fit(fdownpp,"Q");
gmeanpp->Fit(fmeanpp,"");
guppp->Fit(fuppp,"");
gdownpp->Fit(fdownpp,"");

cout<<" MEAN, parameters [0]="<<fmeanpp->GetParameter(0)<<" [1]="<<fmeanpp->GetParameter(1)<<" [2]="<<fmeanpp->GetParameter(2)<<endl;
cout<<" up, parameters [0]="<<fuppp->GetParameter(0)<<" [1]="<<fuppp->GetParameter(1)<<" [2]="<<fuppp->GetParameter(2)<<endl;
cout<<" down, parameters [0]="<<fdownpp->GetParameter(0)<<" [1]="<<fdownpp->GetParameter(1)<<" [2]="<<fdownpp->GetParameter(2)<<endl;


fuppp->SetLineWidth(3);fuppp->SetLineColor(3);fuppp->Draw("same");
fdownpp->SetLineWidth(3);fdownpp->SetLineColor(3);fdownpp->Draw("same");
fmeanpp->SetLineWidth(3);fmeanpp->SetLineColor(3);fmeanpp->Draw("same");
gppstat->SetLineColor(1);gppstat->SetLineWidth(2);
gppsyst->SetLineColor(10);gppsyst->SetFillColor(10);
gppsyst->SetFillStyle(3001);
gppstat->Draw("samee");
gppsyst->Draw("samee2");

//eval radius
Float_t pLambda = 1.55;
Float_t pXi = 1.83;
Float_t pSigma = 2.07;
Float_t pOmega = 2.231;
cout<<" pLambda("<<pLambda<<"):   low "<<fdownpp->Eval(pLambda)<<" mean "<<fmeanpp->Eval(pLambda)<<" up "<<fuppp->Eval(pLambda)<<endl;
cout<<" pXi("<<pXi<<"):   low "<<fdownpp->Eval(pXi)<<" mean "<<fmeanpp->Eval(pXi)<<" up "<<fuppp->Eval(pXi)<<endl;
cout<<" pSigma("<<pSigma<<"):   low "<<fdownpp->Eval(pSigma)<<" mean "<<fmeanpp->Eval(pSigma)<<" up "<<fuppp->Eval(pSigma)<<endl;
cout<<" pOmega("<<pOmega<<"):   low "<<fdownpp->Eval(pOmega)<<" mean "<<fmeanpp->Eval(pOmega)<<" up "<<fuppp->Eval(pOmega)<<endl;


/*
//eval radius for bernie 25/6/2019
Float_t mt1 = 1.21;
Float_t mt2 = 1.29;
Float_t mt3 = 1.9;
Float_t mt4 = 1.54;
Float_t mt5 = 1.74;
Float_t mt6 = 2.25;
cout<<" mt1("<<mt1<<"):   low "<<fdownpp->Eval(mt1)<<" mean "<<fmeanpp->Eval(mt1)<<" up "<<fuppp->Eval(mt1)<<endl;
cout<<" mt2("<<mt2<<"):   low "<<fdownpp->Eval(mt2)<<" mean "<<fmeanpp->Eval(mt2)<<" up "<<fuppp->Eval(mt2)<<endl;
cout<<" mt3("<<mt3<<"):   low "<<fdownpp->Eval(mt3)<<" mean "<<fmeanpp->Eval(mt3)<<" up "<<fuppp->Eval(mt3)<<endl;
cout<<" mt4("<<mt4<<"):   low "<<fdownpp->Eval(mt4)<<" mean "<<fmeanpp->Eval(mt4)<<" up "<<fuppp->Eval(mt4)<<endl;
cout<<" mt5("<<mt5<<"):   low "<<fdownpp->Eval(mt5)<<" mean "<<fmeanpp->Eval(mt5)<<" up "<<fuppp->Eval(mt5)<<endl;
cout<<" mt6("<<mt6<<"):   low "<<fdownpp->Eval(mt6)<<" mean "<<fmeanpp->Eval(mt6)<<" up "<<fuppp->Eval(mt6)<<endl;
*/





//Draw the figure
//_______________
if(DRAWFIG){

    gppsyst->SetFillColorAlpha(kBlue+2,0.3);
    gppsyst->SetLineColor(kBlue+2);
    gppsyst->SetLineWidth(4);
    gppstat->SetFillColorAlpha(kBlue+2,0.3);
    gppstat->SetLineColor(kBlue+2);
    gppstat->SetLineWidth(4);

    TH1F* hAxis = new TH1F("hAxis", "hAxis", 128, 0.9, 2.7);
    hAxis->SetStats(false);
    hAxis->SetTitle("");
    hAxis->GetXaxis()->SetTitle("<m_{T}> (GeV/#it{c}^{2})");
    hAxis->GetXaxis()->SetTitleSize(0.06);
    hAxis->GetXaxis()->SetLabelSize(0.06);
    hAxis->GetXaxis()->CenterTitle();
    hAxis->GetXaxis()->SetTitleOffset(1.3);
    hAxis->GetXaxis()->SetLabelOffset(0.02);
    hAxis->GetYaxis()->SetTitle("Gaussian core (fm)");
    hAxis->GetYaxis()->SetTitleSize(0.06);
    hAxis->GetYaxis()->SetLabelSize(0.06);
    hAxis->GetYaxis()->CenterTitle();
    hAxis->GetYaxis()->SetTitleOffset(0.90);
    hAxis->GetYaxis()->SetRangeUser(0.55, 1.35);

    TCanvas* cmT = new TCanvas("cmT", "cmT", 1);
    cmT->cd(0); cmT->SetCanvasSize(1920/2, 1080/2); cmT->SetMargin(0.15,0.05,0.2,0.05);//lrbt
    hAxis->Draw("axis");
    hAxis->Draw("axis");
    gppsyst->Draw("e2,same");
    gppstat->Draw("e,same");



//create TGraph and then histo with resulting band:
TGraph *gf = new TGraph();
//Int_t npx = fuppp->GetNpx();
Int_t npx = fuppp->GetNpx();
Int_t npoints=0;
//Double_t xmin = 0.9;
//Double_t xmax = 3.0;
//Double_t ymin = 0.4;
//Double_t ymax = 1.9;
  float xmin = 0.89;
  float xmax = 2.71;
  float ymin = 0.75;
  float ymax = 1.42;

Double_t dx = (xmax-xmin)/npx;
Double_t x = xmin+0.5*dx;
while (x <= xmax) {
 //Double_t y = fuppp->Eval(x);
 Double_t y = fuppp->Eval(x);
 if (y < ymin) y = ymin;
 if (y > ymax) y = ymax;
 gf->SetPoint(npoints,x,y);
 npoints++;
 x += dx;
}

x = xmax-0.5*dx;
while (x >= xmin) {
 //Double_t y = fdownpp->Eval(x);
 Double_t y = fdownpp->Eval(x);
 if (y < ymin) y = ymin;
 if (y > ymax) y = ymax;
 gf->SetPoint(npoints,x,y);
 npoints++;
 x -= dx;
}
gf->SetLineColor(3);
gf->SetFillColor(3);
gf->SetFillStyle(3001);
gf->Draw("same");
gf->Draw("f");
    gppsyst->Draw("e2,same");
    gppstat->Draw("e,same");

    TLegend* lLegend = new TLegend(0.60,0.65,0.95,0.95);//lbrt
    lLegend->SetName(TString::Format("lLegend"));
    lLegend->SetTextSize(0.045);
    lLegend->AddEntry(gppsyst,"p#minusp (AV18)");
    lLegend->AddEntry(gf,"r_{core} = a #upoint <m_{T}>^{b} + c");
    lLegend->Draw("same");


//output for preliminary figure:
gppsyst->SetName("gppsyst");
gppstat->SetName("gppstat");
gf->SetName("gf");

/*
TString filename = "";
filename = "MtForFig.root";
TFile fout(filename,"recreate");
gpLNLOsyst->Write();
gpLNLOstat->Write();
gppsyst->Write();
gppstat->Write();
gf->Write();
fout.Close();
*/


}//DRAWFIG





}
