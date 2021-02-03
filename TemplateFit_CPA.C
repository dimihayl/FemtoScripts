TH1F* pri;
TH1F* mat;
TH1F* sec;
TH1F* cont;

TH1F* data;

TGraph* graph_pri = new TGraph();
TGraph* graph_mat = new TGraph();
TGraph* graph_sec = new TGraph();
TGraph* graph_cont = new TGraph();
TGraph* graph_chi = new TGraph();


///change to 8 bins?
const unsigned Num_pT_bins = 8;
float parPri[Num_pT_bins]; //20 pT bins
float parMat[Num_pT_bins]; //20 pT bins
float parSec[Num_pT_bins]; //20 pT bins
float parCont[Num_pT_bins]; //20 pT bins
float dataEntries[Num_pT_bins]; //20 pT bins
float mcEntries[Num_pT_bins]; //20 pT bins
//float SoverL[Num_pT_bins]; //20 pT bins

Double_t ftotal(Double_t *x, Double_t *par) {
   Double_t xx = x[0];
   Int_t bin = data->GetXaxis()->FindBin(xx);
   Double_t apri = par[0]*pri->GetBinContent(bin);
   Double_t amat = par[1]*mat->GetBinContent(bin);
   Double_t asec = par[2]*sec->GetBinContent(bin);
   Double_t acont = par[3]*cont->GetBinContent(bin);
   Double_t out = apri + amat + asec + acont;
   return out;
}


void TemplateFit_CPA(){

gStyle->SetOptStat(0);

//get templates
//"hm" 5% andi analysis, lhc16
///MC file
TFile *_file0 = TFile::Open("$CERN_BOX/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/AnalysisResultsMC.root");
TDirectory *d_0 = (TDirectory*) _file0->Get("MBv0CutsMC5");
TList *l_0 = (TList*) d_0->Get("MBv0CutsMC5");
TList *ll_0 = (TList*) l_0->FindObject("v0MonteCarlo");
TList *lll_0 = (TList*) ll_0->FindObject("CPAPtBinning");
TH2F* DCAPtBinningPri = (TH2F*) lll_0->FindObject("CPAPtBinningPri");
TH2F* DCAPtBinningMat = (TH2F*) lll_0->FindObject("CPAPtBinningMat");
TH2F* DCAPtBinningSec = (TH2F*) lll_0->FindObject("CPAPtBinningSec");
TH2F* DCAPtBinningCont = (TH2F*) lll_0->FindObject("CPAPtBinningCont");

TDirectory *d_a0 = (TDirectory*) _file0->Get("MBAntiv0CutsMC5");
TList *l_a0 = (TList*) d_a0->Get("MBAntiv0CutsMC5");
TList *ll_a0 = (TList*) l_a0->FindObject("v0MonteCarlo");
TList *lll_a0 = (TList*) ll_a0->FindObject("CPAPtBinning");
TH2F* antiDCAPtBinningPri = (TH2F*) lll_a0->FindObject("CPAPtBinningPri");
TH2F* antiDCAPtBinningMat = (TH2F*) lll_a0->FindObject("CPAPtBinningMat");
TH2F* antiDCAPtBinningSec = (TH2F*) lll_a0->FindObject("CPAPtBinningSec");
TH2F* antiDCAPtBinningCont = (TH2F*) lll_0->FindObject("CPAPtBinningCont");

DCAPtBinningPri->Add(antiDCAPtBinningPri);
DCAPtBinningMat->Add(antiDCAPtBinningMat);
DCAPtBinningSec->Add(antiDCAPtBinningSec);
DCAPtBinningCont->Add(antiDCAPtBinningCont);

//get data
//TFile *_file1 = TFile::Open("250251252.root");
///DATA file original
TFile *_file1 = TFile::Open("$CERN_BOX/CatsFiles/ExpData/ALICE_pp_13TeV_HM/Sample10HM/AnalysisResults.root");
//TString cuts = "HM_NEWCUTS2";

//add protons and antiplrotons
  TString folder1 = "HMv0Cuts";
  TString folder2 = "HMAntiv0Cuts";
  //folder1 = cuts + folder1;
  //folder2 = cuts + folder2;

TDirectory *d_1 = (TDirectory*) _file1->Get(folder1);
TList *l_1 = (TList*) d_1->Get(folder1);
TList *ll_1 = (TList*) l_1->FindObject("v0Cuts");
TH2F* DCAXYPtBinningTot_1 = (TH2F*) ll_1->FindObject("CPAPtBinsTot");
TH2F* DCAXYPtMult_0_27_1 = (TH2F*) ll_1->FindObject("CPAPtBinsMult_0_27");
TH2F* DCAXYPtMult_27_55_1 = (TH2F*) ll_1->FindObject("CPAPtBinsMult_27_55");
TH2F* DCAXYPtMult_55_inf_1 = (TH2F*) ll_1->FindObject("CPAPtBinsMult_55_inf");
TList *lll_1 = (TList*) ll_1->FindObject("after");
TH1F* pTDist_after_1 = (TH1F*) lll_1->FindObject("pTDist_after");

TDirectory *d_2 = (TDirectory*) _file1->Get(folder2);
TList *l_2 = (TList*) d_2->Get(folder2);
TList *ll_2 = (TList*) l_2->FindObject("v0Cuts");
TH2F* DCAXYPtBinningTot_2 = (TH2F*) l_2->FindObject("CPAPtBinsTot");
TH2F* DCAXYPtMult_0_27_2 = (TH2F*) l_2->FindObject("CPAPtBinsMult_0_27");
TH2F* DCAXYPtMult_27_55_2 = (TH2F*) l_2->FindObject("CPAPtBinsMult_27_55");
TH2F* DCAXYPtMult_55_inf_2 = (TH2F*) l_2->FindObject("CPAPtBinsMult_55_inf");
TList *lll_2 = (TList*) ll_2->FindObject("after");
TH1F* pTDist_after_2 = (TH1F*) lll_2->FindObject("pTDist_after");

//printf("DCAXYPtBinningTot_2=%p\n",DCAXYPtBinningTot_2);
DCAXYPtBinningTot_1->Add(DCAXYPtBinningTot_2);
DCAXYPtMult_0_27_1->Add(DCAXYPtMult_0_27_2);
DCAXYPtMult_27_55_1->Add(DCAXYPtMult_27_55_2);
DCAXYPtMult_55_inf_1->Add(DCAXYPtMult_55_inf_2);
pTDist_after_1->Add(pTDist_after_2);


TFile* fOut = new TFile("./Output/TemplateFit_CPA.root","recreate");

//which to use?
//tot
TH2F* data2d = (TH2F*) DCAXYPtBinningTot_1->Clone("data2d");
Double_t* BinRanges_pT = new Double_t [data2d->GetXaxis()->GetNbins()+1];
data2d->GetXaxis()->GetLowEdge(BinRanges_pT);
BinRanges_pT[data2d->GetXaxis()->GetNbins()] = data2d->GetXaxis()->GetBinUpEdge(data2d->GetXaxis()->GetNbins());
//TH2F* data2d = DCAXYPtMult_27_55_1->Clone("data2d");
//TH2F* data2d = DCAXYPtMult_55_inf_1->Clone("data2d");

 //build function
TF1 *ftot = new TF1("ftot",ftotal,0.9,1,4);

Int_t npoint = 0;
//loop over the 20 bins
for (int xbin = 1;xbin <= Num_pT_bins; xbin++){
 //project templates
 pri = (TH1F*) DCAPtBinningPri->ProjectionY("pri",xbin,xbin);
 pri->SetName(TString::Format("pri_%i",xbin));
 mat = (TH1F*) DCAPtBinningMat->ProjectionY("mat",xbin,xbin);
 mat->SetName(TString::Format("mat_%i",xbin));
 sec = (TH1F*) DCAPtBinningSec->ProjectionY("sec",xbin,xbin);
 sec->SetName(TString::Format("sec_%i",xbin));
 cont = (TH1F*) DCAPtBinningCont->ProjectionY("cont",xbin,xbin);
 cont->SetName(TString::Format("cont_%i",xbin));

 //get sigma/lambda fraction from MC for the specific pt bin from MC
 //(in DCA integrated, since it is where we fit)
 //SoverL[xbin-1] = secS->Integral(1,secS->GetNbinsX())/secL->Integral(1,secL->GetNbinsX());
 TCanvas *Can1 = new TCanvas(TString::Format("Can1_%i",xbin),TString::Format("Can1_%i",xbin), 0,0,650,550);
 Can1->cd();
 //normalize
 pri->Scale(1./pri->Integral());
 mat->Scale(1./mat->Integral());
 sec->Scale(1./sec->Integral());
 cont->Scale(1./cont->Integral());

 gPad->SetLogy();
 //project data
 data = (TH1F*) data2d->ProjectionY(TString::Format("data_%i",xbin),xbin,xbin);
  Float_t momLOW  = data2d->GetXaxis()->GetBinLowEdge(xbin);
  Float_t momHIGH  = data2d->GetXaxis()->GetBinLowEdge(xbin)+data2d->GetXaxis()->GetBinWidth(xbin);

 //get fractions for weighted total fractions in the cut   trackCuts->SetDCAVtxXY(0.1);
 //float dcacut = 0.1;
 //default = 0.99
 const float CPAcut = 0.9999;
 dataEntries[xbin-1] =  data->Integral(data->FindBin(CPAcut),data->FindBin(1));
 //normalize
 data->Scale(1./data->Integral());

 ftot->SetParameters(0.6,0.05,0.25,0.1);
 ftot->SetParLimits(0,0.,1.);
 ftot->SetParLimits(1,0.,1.);
 ftot->SetParLimits(2,0.,1.);
 ftot->SetParLimits(3,0.,1.);

 data->Fit("ftot","S, N, R, M","",0.9,1);

 data->SetLineColor(kBlack);data->SetLineWidth(2);data->Draw("hist");
 pri->Scale(ftot->GetParameter(0));pri->SetLineColor(kRed+2);pri->SetLineWidth(2);pri->Draw("same");
 mat->Scale(ftot->GetParameter(1));mat->SetLineColor(kGreen+3);mat->SetLineWidth(2);mat->Draw("same");
 sec->Scale(ftot->GetParameter(2));sec->SetLineColor(kBlue+1);sec->SetLineWidth(2);sec->Draw("same");
 cont->Scale(ftot->GetParameter(3));cont->SetLineColor(kOrange+1);cont->SetLineWidth(2);cont->Draw("same");
 TH1F* htot = (TH1F*)pri->Clone(TString::Format("htot_%i",xbin));
 htot->SetName(TString::Format("htot_%i",xbin));
 htot->Add(mat);
 htot->Add(sec);
 htot->Add(cont);
 htot->SetLineWidth(2);htot->SetLineColor(kRed);htot->Draw("same");

	fOut->cd();
	pri->Write();
	mat->Write();
	sec->Write();
	cont->Write();
	htot->Write();
	data->Write();



 mcEntries[xbin-1] =  htot->Integral(htot->FindBin(CPAcut),htot->FindBin(1));
 parPri[xbin-1] = pri->Integral(pri->FindBin(CPAcut),pri->FindBin(1))/mcEntries[xbin-1];
 parMat[xbin-1] = mat->Integral(mat->FindBin(CPAcut),mat->FindBin(1))/mcEntries[xbin-1];
 parSec[xbin-1] = sec->Integral(sec->FindBin(CPAcut),sec->FindBin(1))/mcEntries[xbin-1];
 parCont[xbin-1] = cont->Integral(cont->FindBin(CPAcut),cont->FindBin(1))/mcEntries[xbin-1];

 double TotParticles = htot->Integral(1, htot->GetNbinsX());
 double SelectedParticles = htot->Integral(htot->FindBin(CPAcut),htot->FindBin(1)-1);
 double TotPrim = pri->Integral(1, pri->GetNbinsX());
 double SelectedPrim = pri->Integral(pri->FindBin(CPAcut),pri->FindBin(1)-1);
 printf("tot: %i --> %i\n",1, pri->GetNbinsX());
 printf("pri: %i --> %i\n",pri->FindBin(CPAcut),pri->FindBin(1));
 printf("You have selected %.1f %% of all primaries, with an S/(S+B) of %.1f %%\n",
 SelectedPrim/TotPrim*100,SelectedPrim/SelectedParticles*100);

cout<<"par fit Pri "<<ftot->GetParameter(0)<<"  vs calculation in dca limits = "<<parPri[xbin-1]<<endl;
cout<<"par fit Mat "<<ftot->GetParameter(1)<<"  vs calculation in dca limits = "<<parMat[xbin-1]<<endl;
cout<<"par fit Sec "<<ftot->GetParameter(2)<<"  vs calculation in dca limits = "<<parSec[xbin-1]<<endl;
cout<<"par fit Cont "<<ftot->GetParameter(3)<<"  vs calculation in dca limits = "<<parCont[xbin-1]<<endl;
//cout<<"par fit SecS "<<SoverL[xbin-1]*ftot->GetParameter(2)<<"  vs calculation in dca limits = "<<parSecS[xbin-1]<<endl;

 TString textpri = "prim. ";
 TString textmat = "mat.  ";
 TString textsec = "sec.  ";
 TString textcont = "fakes  ";
 textpri = textpri + Form("%.0f",parPri[xbin-1]*100.) + "%   mom=" + Form("%.1f",momLOW) + " - "  + Form("%.1f",momHIGH) ;
;
 textmat = textmat + Form("%.0f",parMat[xbin-1]*100.) + "%";
 textsec = textsec + Form("%.0f",parSec[xbin-1]*100.) + "%";
 textcont = textcont + Form("%.0f",parCont[xbin-1]*100.) + "%";
 //textsecS = textsecS + Form("%.0f",parSecS[xbin-1]*100.) + "%";
 TLegend *leg2 = new TLegend( 0.55, 0.5, 0.90, 0.85);
 leg2->AddEntry(pri,textpri);
 leg2->AddEntry(mat,textmat);
 leg2->AddEntry(sec,textsec);
 leg2->AddEntry(cont,textcont);
 leg2->SetLineColor(0);
 leg2->Draw();

 gPad->SetLogy();
 TString pout = "./Output/";
 pout = pout + Form("%i",xbin) + "cpa.png";
 gPad->SaveAs(pout);

graph_pri->SetName("graph_pri");
graph_mat->SetName("graph_mat");
graph_sec->SetName("graph_sec");
graph_cont->SetName("graph_cont");
graph_chi->SetName("graph_chi");

graph_pri->SetPoint(npoint,data2d->GetXaxis()->GetBinCenter(xbin),parPri[xbin-1]);
graph_mat->SetPoint(npoint,data2d->GetXaxis()->GetBinCenter(xbin),parMat[xbin-1]);
graph_sec->SetPoint(npoint,data2d->GetXaxis()->GetBinCenter(xbin),parSec[xbin-1]);
graph_cont->SetPoint(npoint,data2d->GetXaxis()->GetBinCenter(xbin),parCont[xbin-1]);
//graph_secS->SetPoint(npoint,data2d->GetXaxis()->GetBinCenter(xbin),parSecS[xbin-1]);
graph_chi->SetPoint(npoint,data2d->GetXaxis()->GetBinCenter(xbin),ftot->GetChisquare()/ftot->GetNDF());
npoint++;
}





 TCanvas *c2 = new TCanvas(TString::Format("c2"),TString::Format("c2"), 0,0,650,550);
 c2->cd();
graph_chi->SetLineColor(1);graph_chi->Draw("");
//calcular final pT weighted result
float m2tot = 0.;
float m2totpri = 0.;
float m2totmat = 0.;
float m2totsec = 0.;
float m2totcont = 0.;
for(int i=0;i<Num_pT_bins;i++){
 m2totpri += dataEntries[i]*parPri[i];
 m2totmat += dataEntries[i]*parMat[i];
 m2totsec += dataEntries[i]*parSec[i];
 m2totcont += dataEntries[i]*parCont[i];
 m2tot += dataEntries[i];
}

TH1F* pT_Weights = new TH1F("pT_Weights","pT_Weights",Num_pT_bins,BinRanges_pT);
for(int i=0;i<Num_pT_bins;i++){
	pT_Weights->SetBinContent(i+1,dataEntries[i]/m2tot);
	printf("Weight %i = %f\n",i,dataEntries[i]/m2tot);
}

//d12();
 TCanvas *d12 = new TCanvas(TString::Format("d12"),TString::Format("d12"), 0,0,650,550);
 d12->cd();
graph_pri->SetLineWidth(2);
graph_mat->SetLineWidth(2);
graph_sec->SetLineWidth(2);
graph_cont->SetLineWidth(2);
graph_chi->SetLineWidth(2);
graph_pri->GetXaxis()->SetLabelSize(0.05);
graph_pri->GetYaxis()->SetLabelSize(0.05);
graph_pri->GetXaxis()->SetTitleSize(0.05);
graph_pri->GetXaxis()->SetTitle("p (GeV)");
graph_chi->GetXaxis()->SetLabelSize(0.05);
graph_chi->GetYaxis()->SetLabelSize(0.05);
graph_chi->GetYaxis()->SetTitle("chi2/ndf");
graph_chi->GetXaxis()->SetTitle("p (GeV)");
graph_chi->GetXaxis()->SetTitleSize(0.05);
graph_chi->GetYaxis()->SetTitleSize(0.05);
graph_pri->SetLineColor(kRed+2);graph_pri->GetYaxis()->SetRangeUser(0.,1.);graph_pri->Draw("");
graph_mat->SetLineColor(kGreen+3);graph_mat->Draw("same");
graph_sec->SetLineColor(kBlue+1);graph_sec->Draw("same");
graph_cont->SetLineColor(kOrange+1);graph_cont->Draw("same");
TF1* primAvg = new TF1("primAvg","[0]",0,5);
primAvg->SetParameter(0,m2totpri/m2tot);
primAvg->SetLineColor(kBlack);
primAvg->SetLineStyle(2);
primAvg->Draw("same");

cout<<endl;
cout<<" pri = "<<m2totpri/m2tot<<endl;
cout<<" mat = "<<m2totmat/m2tot<<endl;
cout<<" sec = "<<m2totsec/m2tot<<endl;
cout<<" fak = "<<m2totcont/m2tot<<endl;

double Renorm = 1.-m2totcont/m2tot-m2totmat/m2tot;
cout<<" frac pri = "<<m2totpri/m2tot/Renorm<<endl;
cout<<" frac sec = "<<m2totsec/m2tot/Renorm<<endl;

cout<<" frac lambda = "<<m2totpri/m2tot/Renorm*0.75<<endl;
cout<<" frac sigma0 = "<<m2totpri/m2tot/Renorm*0.25<<endl;
cout<<" frac xim = "<<m2totsec/m2tot/Renorm*0.5<<endl;
cout<<" frac rest = "<<m2totsec/m2tot/Renorm*0.5<<endl;

const float pT_val = 2.00;
double Dimi_LS = graph_pri->Eval(pT_val);
double Dimi_mat = graph_mat->Eval(pT_val);
double Dimi_Xi = graph_sec->Eval(pT_val);
double Dimi_fake = graph_cont->Eval(pT_val);
double Dimi_Tot = Dimi_LS+Dimi_mat+Dimi_Xi+Dimi_fake;
double Dimi_Renorm = 1./Dimi_Tot;
Dimi_LS /= Dimi_Renorm;
Dimi_mat /= Dimi_Renorm;
Dimi_Xi /= Dimi_Renorm;
Dimi_fake /= Dimi_Renorm;
Dimi_Tot /= Dimi_Renorm;
printf("Estimation based on pT = %.2f\n",pT_val);
printf(" Renrom = %.2f %%\n",Dimi_Renorm*100.);
printf(" Λ+Σ = %.2f (%.2f) [%.2f] %%\n",Dimi_LS*100.,Dimi_LS/(1.-Dimi_fake)*100.,Dimi_LS/(1.-Dimi_fake-Dimi_mat)*100.);
printf(" Ξ = %.2f (%.2f) [%.2f] %%\n",Dimi_Xi*100.,Dimi_Xi/(1.-Dimi_fake)*100.,Dimi_Xi/(1.-Dimi_fake-Dimi_mat)*100.);
printf(" mat = %.2f (%.2f) [%.2f] %%\n",Dimi_mat*100.,Dimi_mat/(1.-Dimi_fake)*100.,0.);
printf(" fake = %.2f (%.2f) [%.2f] %%\n",Dimi_fake*100.,0.,0.);

fOut->cd();
graph_chi->Write();
graph_pri->Write();
graph_mat->Write();
graph_sec->Write();
graph_cont->Write();
primAvg->Write();
pTDist_after_1->Write();
pT_Weights->Write();
}
