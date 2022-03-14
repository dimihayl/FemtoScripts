

void TemplatePlot_CPA(){
	unsigned Num_pT_bins = 8;
	unsigned NumCols = 3;
	unsigned NumRows = Num_pT_bins/NumCols + bool(Num_pT_bins%NumCols);

	double MinRange = 0.9;
	double MaxRange = 1.0;
	double ZoomedMinRange = 0.99;
	double ZoomedMaxRange = 1.0;

	TString InputFileName = "$FEMTO_OUTPUT/Scripts/TemplateFit_CPA.root";

	//in the last bin, we add all of them up
	TH1F** pri = new TH1F* [Num_pT_bins+1];
	TH1F** mat = new TH1F* [Num_pT_bins+1];
	TH1F** sec = new TH1F* [Num_pT_bins+1];
	TH1F** cont = new TH1F* [Num_pT_bins+1];
	TH1F** htot = new TH1F* [Num_pT_bins+1];
	TH1F** data = new TH1F* [Num_pT_bins+1];
	TH1F* pT_Weights;
	TH1F* pT_Dist;
	TGraph* graph_pri;
	TGraph* graph_mat;
	TGraph* graph_sec;
	TGraph* graph_cont;
	TF1* primAvg;


	TFile* InputFile = new TFile(InputFileName,"read");

	for(int ipT=0;ipT<Num_pT_bins;ipT++){
		pri[ipT] = (TH1F*)InputFile->Get(TString::Format("pri_%i",ipT+1));
		mat[ipT] = (TH1F*)InputFile->Get(TString::Format("mat_%i",ipT+1));
		sec[ipT] = (TH1F*)InputFile->Get(TString::Format("sec_%i",ipT+1));
		cont[ipT] = (TH1F*)InputFile->Get(TString::Format("cont_%i",ipT+1));
		htot[ipT] = (TH1F*)InputFile->Get(TString::Format("htot_%i",ipT+1));
		data[ipT] = (TH1F*)InputFile->Get(TString::Format("data_%i",ipT+1));
	}
	pT_Weights = (TH1F*)InputFile->Get(TString::Format("pT_Weights"));
	pT_Dist = (TH1F*)InputFile->Get(TString::Format("pTDist_after"));
	graph_pri = (TGraph*)InputFile->Get(TString::Format("graph_pri"));
	graph_mat = (TGraph*)InputFile->Get(TString::Format("graph_mat"));
	graph_sec = (TGraph*)InputFile->Get(TString::Format("graph_sec"));
	graph_cont = (TGraph*)InputFile->Get(TString::Format("graph_cont"));
	primAvg = (TF1*)InputFile->Get(TString::Format("primAvg"));

	unsigned NumCpaBins = data[0]->GetNbinsX();

	//pri[Num_pT_bins] = new TH1F("pri_tot","pri_tot",pri[0]->GetNbinsX(),pri[0]->GetBinLowEdge(1),pri[0]->GetXaxis()->GetBinUpEdge(pri[0]->GetNbinsX()));
	//mat[Num_pT_bins] = new TH1F("mat_tot","mat_tot",mat[0]->GetNbinsX(),mat[0]->GetBinLowEdge(1),mat[0]->GetXaxis()->GetBinUpEdge(mat[0]->GetNbinsX()));
	//sec[Num_pT_bins] = new TH1F("sec_tot","sec_tot",sec[0]->GetNbinsX(),sec[0]->GetBinLowEdge(1),sec[0]->GetXaxis()->GetBinUpEdge(sec[0]->GetNbinsX()));
	//cont[Num_pT_bins] = new TH1F("cont_tot","cont_tot",cont[0]->GetNbinsX(),cont[0]->GetBinLowEdge(1),cont[0]->GetXaxis()->GetBinUpEdge(cont[0]->GetNbinsX()));
	//data[Num_pT_bins] = new TH1F("data_tot","data_tot",data[0]->GetNbinsX(),data[0]->GetBinLowEdge(1),data[0]->GetXaxis()->GetBinUpEdge(data[0]->GetNbinsX()));

	pri[Num_pT_bins] = (TH1F*)pri[0]->Clone("pri_tot");
	mat[Num_pT_bins] = (TH1F*)mat[0]->Clone("mat_tot");
	sec[Num_pT_bins] = (TH1F*)sec[0]->Clone("sec_tot");
	cont[Num_pT_bins] = (TH1F*)cont[0]->Clone("cont_tot");
	htot[Num_pT_bins] = (TH1F*)htot[0]->Clone("htot_tot");
	data[Num_pT_bins] = (TH1F*)data[0]->Clone("data_tot");

	pri[Num_pT_bins]->Scale(pT_Weights->GetBinContent(1));
	mat[Num_pT_bins]->Scale(pT_Weights->GetBinContent(1));
	sec[Num_pT_bins]->Scale(pT_Weights->GetBinContent(1));
	cont[Num_pT_bins]->Scale(pT_Weights->GetBinContent(1));
	htot[Num_pT_bins]->Scale(pT_Weights->GetBinContent(1));
	data[Num_pT_bins]->Scale(pT_Weights->GetBinContent(1));

	for(int ipT=1;ipT<Num_pT_bins;ipT++){
		pri[Num_pT_bins]->Add(pri[ipT],pT_Weights->GetBinContent(ipT+1));
		mat[Num_pT_bins]->Add(mat[ipT],pT_Weights->GetBinContent(ipT+1));
		sec[Num_pT_bins]->Add(sec[ipT],pT_Weights->GetBinContent(ipT+1));
		cont[Num_pT_bins]->Add(cont[ipT],pT_Weights->GetBinContent(ipT+1));
		htot[Num_pT_bins]->Add(htot[ipT],pT_Weights->GetBinContent(ipT+1));
		data[Num_pT_bins]->Add(data[ipT],pT_Weights->GetBinContent(ipT+1));

	}

	TFile* fOutPlot = new TFile("$FEMTO_OUTPUT/Scripts/fOutDcaPlot.root","recreate");

	TGraph* gpri = new TGraph [Num_pT_bins+1];
	TGraph* gmat = new TGraph [Num_pT_bins+1];
	TGraph* gsec = new TGraph [Num_pT_bins+1];
	TGraph* gcont = new TGraph [Num_pT_bins+1];
	TGraph* gtot = new TGraph [Num_pT_bins+1];
	TGraph* gdata = new TGraph [Num_pT_bins+1];
	unsigned* Num_g_bins_pri = new unsigned[Num_pT_bins+1];
	unsigned* Num_g_bins_mat = new unsigned[Num_pT_bins+1];
	unsigned* Num_g_bins_sec = new unsigned[Num_pT_bins+1];
	unsigned* Num_g_bins_cont = new unsigned[Num_pT_bins+1];
	unsigned* Num_g_bins_tot = new unsigned[Num_pT_bins+1];
	unsigned* Num_g_bins_data = new unsigned[Num_pT_bins+1];

	for(int ipT=0;ipT<=Num_pT_bins;ipT++){
		Num_g_bins_pri[ipT]=0;
		Num_g_bins_mat[ipT]=0;
		Num_g_bins_sec[ipT]=0;
		Num_g_bins_cont[ipT]=0;
		Num_g_bins_tot[ipT]=0;
		Num_g_bins_data[ipT]=0;

		gpri[ipT].SetName(TString::Format("gpri_%i",ipT+1));
		gmat[ipT].SetName(TString::Format("gmat_%i",ipT+1));
		gsec[ipT].SetName(TString::Format("gsec_%i",ipT+1));
		gcont[ipT].SetName(TString::Format("gcont_%i",ipT+1));
		gtot[ipT].SetName(TString::Format("gtot_%i",ipT+1));
		gdata[ipT].SetName(TString::Format("gdata_%i",ipT+1));

		//gpri[ipT].Set(NumCpaBins);
		//gmat[ipT].Set(NumCpaBins);
		//gsec[ipT].Set(NumCpaBins);
		//gcont[ipT].Set(NumCpaBins);
		//gtot[ipT].Set(NumCpaBins);
		//gdata[ipT].Set(NumCpaBins);



		for(unsigned uCpa=0; uCpa<NumCpaBins; uCpa++){
			if(pri[ipT]->GetBinContent(uCpa+1)) {gpri[ipT].SetPoint(Num_g_bins_pri[ipT],pri[ipT]->GetBinCenter(uCpa+1),pri[ipT]->GetBinContent(uCpa+1));Num_g_bins_pri[ipT]++;}
			if(mat[ipT]->GetBinContent(uCpa+1)) {gmat[ipT].SetPoint(Num_g_bins_mat[ipT],mat[ipT]->GetBinCenter(uCpa+1),mat[ipT]->GetBinContent(uCpa+1));Num_g_bins_mat[ipT]++;}
			if(sec[ipT]->GetBinContent(uCpa+1)) {gsec[ipT].SetPoint(Num_g_bins_sec[ipT],sec[ipT]->GetBinCenter(uCpa+1),sec[ipT]->GetBinContent(uCpa+1));Num_g_bins_sec[ipT]++;}
			if(cont[ipT]->GetBinContent(uCpa+1)) {gcont[ipT].SetPoint(Num_g_bins_cont[ipT],cont[ipT]->GetBinCenter(uCpa+1),cont[ipT]->GetBinContent(uCpa+1));Num_g_bins_cont[ipT]++;}
			if(htot[ipT]->GetBinContent(uCpa+1)) {gtot[ipT].SetPoint(Num_g_bins_tot[ipT],htot[ipT]->GetBinCenter(uCpa+1),htot[ipT]->GetBinContent(uCpa+1));Num_g_bins_tot[ipT]++;}
			if(data[ipT]->GetBinContent(uCpa+1)) {gdata[ipT].SetPoint(Num_g_bins_data[ipT],data[ipT]->GetBinCenter(uCpa+1),data[ipT]->GetBinContent(uCpa+1));Num_g_bins_data[ipT]++;}

			if(ipT<Num_pT_bins){
				pri[ipT]->SetBinError(uCpa+1,0);
				mat[ipT]->SetBinError(uCpa+1,0);
				sec[ipT]->SetBinError(uCpa+1,0);
				cont[ipT]->SetBinError(uCpa+1,0);
				htot[ipT]->SetBinError(uCpa+1,0);
				data[ipT]->SetBinError(uCpa+1,0);
			}

		}

		pri[ipT]->SetLineColor(kRed+2);
		mat[ipT]->SetLineColor(kGreen+3);
		sec[ipT]->SetLineColor(kBlue+1);
		cont[ipT]->SetLineColor(kOrange+1);
		htot[ipT]->SetLineColor(kRed);
		data[ipT]->SetLineColor(kBlack);

		pri[ipT]->SetLineWidth(3);
		mat[ipT]->SetLineWidth(3);
		sec[ipT]->SetLineWidth(3);
		cont[ipT]->SetLineWidth(3);
		htot[ipT]->SetLineWidth(3);
		data[ipT]->SetLineWidth(4);

		//pri[ipT]->SetLineStyle(2);
		//mat[ipT]->SetLineStyle(3);
		//sec[ipT]->SetLineStyle(4);
		//cont[ipT]->SetLineStyle(5);
		//htot[ipT]->SetLineStyle(1);
		//data[ipT]->SetLineStyle(1);

		gpri[ipT].SetLineColor(kRed+2);
		gmat[ipT].SetLineColor(kGreen+3);
		gsec[ipT].SetLineColor(kBlue+1);
		gcont[ipT].SetLineColor(kOrange+1);
		gtot[ipT].SetLineColor(kRed);
		gdata[ipT].SetLineColor(kBlack);

		gpri[ipT].SetLineWidth(4);
		gmat[ipT].SetLineWidth(4);
		gsec[ipT].SetLineWidth(4);
		gcont[ipT].SetLineWidth(4);
		gtot[ipT].SetLineWidth(4);
		gdata[ipT].SetLineWidth(6);

		gpri[ipT].SetLineStyle(2);
		gmat[ipT].SetLineStyle(3);
		gsec[ipT].SetLineStyle(4);
		gcont[ipT].SetLineStyle(5);
		gtot[ipT].SetLineStyle(1);
		gdata[ipT].SetLineStyle(1);

		//gpri[ipT].SetMarkerSize(0);
		//gmat[ipT].SetMarkerSize(0);
		//gsec[ipT].SetMarkerSize(0);
		//gcont[ipT].SetMarkerSize(0);
		//gdata[ipT].SetMarkerSize(2);

		//gdata[ipT].SetMarkerStyle(34);

		pri[ipT]->Write();
		mat[ipT]->Write();
		sec[ipT]->Write();
		cont[ipT]->Write();
		htot[ipT]->Write();
		data[ipT]->Write();

		gpri[ipT].Write();
		gmat[ipT].Write();
		gsec[ipT].Write();
		gcont[ipT].Write();
		gtot[ipT].Write();
		gdata[ipT].Write();
	}

	TPaveText** PT1 = new TPaveText*[Num_pT_bins];
	for(int ipT=0;ipT<Num_pT_bins;ipT++){
		PT1[ipT] = new TPaveText(0.2,0.725,0.65,0.875, "blNDC");//lbrt
		PT1[ipT]->SetName(TString::Format("PT1_%i",ipT));
		PT1[ipT]->SetBorderSize(1);
		PT1[ipT]->SetTextSize(0.090);
		PT1[ipT]->SetFillColor(kWhite);
		PT1[ipT]->SetTextFont(22);
		PT1[ipT]->AddText(TString::Format("#it{p}_{T}#in[%.1f, %.1f] GeV/#it{c}",pT_Weights->GetXaxis()->GetBinLowEdge(ipT+1),pT_Weights->GetXaxis()->GetBinUpEdge(ipT+1)));
	}



	TCanvas* cAll = new TCanvas("cAll","cAll",1);
	cAll->cd(0);
	cAll->SetCanvasSize(1920, 1080);
	cAll->Divide(3,3);

	cAll->cd(1);
    TLegend* LegendAll = new TLegend(0.1,0.1,0.9,0.9);//lbrt
    LegendAll->SetName(TString::Format("LegendAll"));
    LegendAll->SetTextSize(0.090);
    LegendAll->AddEntry(pri[0],"Primary (#Lambda or #Sigma^{0}#rightarrow#Lambda)");
    LegendAll->AddEntry(sec[0],"Secondary (#Xi#rightarrow#Lambda)");
    LegendAll->AddEntry(mat[0],"Material #Lambda");
    LegendAll->AddEntry(cont[0],"Fake #Lambda");
    LegendAll->AddEntry(htot[0],"Template fit");
    LegendAll->AddEntry(data[0],"Data");
    LegendAll->Draw("same");

	for(int ipT=0;ipT<Num_pT_bins;ipT++){
		cAll->cd(ipT+2);
		cAll->cd(ipT+2)->SetLogy();
		cAll->cd(ipT+2)->SetMargin(0.15,0.05,0.2,0.05);//lrbt
		mat[ipT]->SetStats(false);
		mat[ipT]->SetTitle("");
		mat[ipT]->GetXaxis()->SetLabelSize(0.090);
		mat[ipT]->GetXaxis()->SetTitle("cos(#alpha)");
		mat[ipT]->GetXaxis()->CenterTitle();
		mat[ipT]->GetXaxis()->SetTitleOffset(1.2);
		mat[ipT]->GetXaxis()->SetLabelOffset(0.02);
		mat[ipT]->GetXaxis()->SetTitleSize(0.090);
		mat[ipT]->GetYaxis()->SetLabelSize(0.090);
		mat[ipT]->GetYaxis()->SetTitle("dN/N");
		mat[ipT]->GetYaxis()->CenterTitle();
		mat[ipT]->GetYaxis()->SetTitleOffset(0.90);
		mat[ipT]->GetYaxis()->SetTitleSize(0.090);

		mat[ipT]->GetYaxis()->SetRangeUser(1e-6,0.5);
		mat[ipT]->GetXaxis()->SetNdivisions(506);
		mat[ipT]->GetXaxis()->SetRangeUser(MinRange, MaxRange);

		mat[ipT]->Draw();
		cont[ipT]->Draw("same");
		sec[ipT]->Draw("same");
		pri[ipT]->Draw("same");
		htot[ipT]->Draw("same");
		data[ipT]->Draw("same");
		PT1[ipT]->Draw("same");
	}

	gStyle->SetLineWidth(2.5);
	cAll->SaveAs("$FEMTO_OUTPUT/Scripts/TemplateFitsCPA.png");
	gStyle->SetLineWidth(1);
	for(int ipT=0;ipT<Num_pT_bins;ipT++){
		pri[ipT]->SetLineWidth(pri[ipT]->GetLineWidth()/2.5);
		mat[ipT]->SetLineWidth(mat[ipT]->GetLineWidth()/2.5);
		sec[ipT]->SetLineWidth(sec[ipT]->GetLineWidth()/2.5);
		cont[ipT]->SetLineWidth(cont[ipT]->GetLineWidth()/2.5);
		htot[ipT]->SetLineWidth(htot[ipT]->GetLineWidth()/2.5);
		data[ipT]->SetLineWidth(data[ipT]->GetLineWidth()/2.5);
	}
	cAll->SaveAs("$FEMTO_OUTPUT/Scripts/TemplateFitsCPA.pdf");
	gStyle->SetLineWidth(2.5);
	for(int ipT=0;ipT<Num_pT_bins;ipT++){
		pri[ipT]->SetLineWidth(pri[ipT]->GetLineWidth()*2.5);
		mat[ipT]->SetLineWidth(mat[ipT]->GetLineWidth()*2.5);
		sec[ipT]->SetLineWidth(sec[ipT]->GetLineWidth()*2.5);
		cont[ipT]->SetLineWidth(cont[ipT]->GetLineWidth()*2.5);
		htot[ipT]->SetLineWidth(htot[ipT]->GetLineWidth()*2.5);
		data[ipT]->SetLineWidth(data[ipT]->GetLineWidth()*2.5);
	}

	//cSource->SetMargin(0.15,0.05,0.2,0.05);//lrbt

	TCanvas* cTot = new TCanvas("cTot","cTot",1);
	cTot->cd(0);
	cTot->SetCanvasSize(1920, 1080);
	cTot->SetMargin(0.15,0.05,0.2,0.05);//lrbt
	cTot->SetLogy();

	pri[Num_pT_bins]->SetLineWidth(pri[Num_pT_bins]->GetLineWidth()*2);
	mat[Num_pT_bins]->SetLineWidth(mat[Num_pT_bins]->GetLineWidth()*2);
	sec[Num_pT_bins]->SetLineWidth(sec[Num_pT_bins]->GetLineWidth()*2);
	cont[Num_pT_bins]->SetLineWidth(cont[Num_pT_bins]->GetLineWidth()*2);
	htot[Num_pT_bins]->SetLineWidth(htot[Num_pT_bins]->GetLineWidth()*4);
	data[Num_pT_bins]->SetLineWidth(data[Num_pT_bins]->GetLineWidth()*2);

	mat[Num_pT_bins]->SetStats(false);
	mat[Num_pT_bins]->SetTitle("");
	mat[Num_pT_bins]->GetXaxis()->SetLabelSize(0.06);
	mat[Num_pT_bins]->GetXaxis()->SetTitle("cos(#alpha)");
	mat[Num_pT_bins]->GetXaxis()->CenterTitle();
	mat[Num_pT_bins]->GetXaxis()->SetTitleOffset(1.2);
	mat[Num_pT_bins]->GetXaxis()->SetLabelOffset(0.02);
	mat[Num_pT_bins]->GetXaxis()->SetTitleSize(0.06);
	mat[Num_pT_bins]->GetYaxis()->SetLabelSize(0.06);
	mat[Num_pT_bins]->GetYaxis()->SetTitle("dN/N");
	mat[Num_pT_bins]->GetYaxis()->CenterTitle();
	mat[Num_pT_bins]->GetYaxis()->SetTitleOffset(1.00);
	mat[Num_pT_bins]->GetYaxis()->SetTitleSize(0.06);

	mat[Num_pT_bins]->GetYaxis()->SetRangeUser(1e-6,0.5);
	mat[Num_pT_bins]->GetXaxis()->SetNdivisions(506);
	mat[Num_pT_bins]->GetXaxis()->SetRangeUser(ZoomedMinRange, ZoomedMaxRange);

    TLegend* LegendTot = new TLegend(0.2,0.65,0.55,0.925);//lbrt
    LegendTot->SetName(TString::Format("LegendTot"));
    LegendTot->SetTextSize(0.035);
    LegendTot->AddEntry(pri[Num_pT_bins],"Primary (#Lambda or #Sigma^{0}#rightarrow#Lambda)");
    LegendTot->AddEntry(sec[Num_pT_bins],"Secondary (#Xi#rightarrow#Lambda)");
    LegendTot->AddEntry(mat[Num_pT_bins],"Material #Lambda");
    LegendTot->AddEntry(cont[Num_pT_bins],"Fake #Lambda");
    LegendTot->AddEntry(htot[Num_pT_bins],"Template fit");
    LegendTot->AddEntry(data[Num_pT_bins],"Data");

	mat[Num_pT_bins]->Draw();
	cont[Num_pT_bins]->Draw("same");
	sec[Num_pT_bins]->Draw("same");
	pri[Num_pT_bins]->Draw("same");
	htot[Num_pT_bins]->Draw("same");
	data[Num_pT_bins]->Draw("same");
	LegendTot->Draw("same");

	cTot->SaveAs("$FEMTO_OUTPUT/Scripts/TemplateFitCPA.png");
	gStyle->SetLineWidth(1);

	pri[Num_pT_bins]->SetLineWidth(pri[Num_pT_bins]->GetLineWidth()/2.5);
	mat[Num_pT_bins]->SetLineWidth(mat[Num_pT_bins]->GetLineWidth()/2.5);
	sec[Num_pT_bins]->SetLineWidth(sec[Num_pT_bins]->GetLineWidth()/2.5);
	cont[Num_pT_bins]->SetLineWidth(cont[Num_pT_bins]->GetLineWidth()/2.5);
	htot[Num_pT_bins]->SetLineWidth(htot[Num_pT_bins]->GetLineWidth()/2.5);
	data[Num_pT_bins]->SetLineWidth(data[Num_pT_bins]->GetLineWidth()/2.5);

	cTot->SaveAs("$FEMTO_OUTPUT/Scripts/TemplateFitCPA.pdf");
	gStyle->SetLineWidth(2.5);
	pri[Num_pT_bins]->SetLineWidth(pri[Num_pT_bins]->GetLineWidth()*2.5);
	mat[Num_pT_bins]->SetLineWidth(mat[Num_pT_bins]->GetLineWidth()*2.5);
	sec[Num_pT_bins]->SetLineWidth(sec[Num_pT_bins]->GetLineWidth()*2.5);
	cont[Num_pT_bins]->SetLineWidth(cont[Num_pT_bins]->GetLineWidth()*2.5);
	htot[Num_pT_bins]->SetLineWidth(htot[Num_pT_bins]->GetLineWidth()*2.5);
	data[Num_pT_bins]->SetLineWidth(data[Num_pT_bins]->GetLineWidth()*2.5);

	TCanvas* cpT = new TCanvas("cpT","cpT",1);

	cpT->cd(0);
	cpT->SetCanvasSize(2560, 1080);

	cpT->Divide(2,1);
	cpT->cd(1);
	cpT->cd(1)->SetMargin(0.15,0.05,0.2,0.05);//lrbt
	cpT->cd(1)->SetLogy();

	pT_Dist->SetStats(false);
	pT_Dist->SetTitle("");
	pT_Dist->GetXaxis()->SetLabelSize(0.055);
	pT_Dist->GetXaxis()->SetTitle("p_{T} (GeV)");
	pT_Dist->GetXaxis()->CenterTitle();
	pT_Dist->GetXaxis()->SetTitleOffset(1.2);
	pT_Dist->GetXaxis()->SetLabelOffset(0.02);
	pT_Dist->GetXaxis()->SetTitleSize(0.055);
	pT_Dist->GetYaxis()->SetLabelSize(0.055);
	pT_Dist->GetYaxis()->SetTitle("dN/dp_{T} (1/GeV)");
	pT_Dist->GetYaxis()->CenterTitle();
	pT_Dist->GetYaxis()->SetTitleOffset(1.25);
	pT_Dist->GetYaxis()->SetTitleSize(0.055);

	//pT_Dist->GetYaxis()->SetRangeUser(1e3,1e7);
	pT_Dist->GetXaxis()->SetNdivisions(506);
	//pT_Dist->GetXaxis()->SetRangeUser(ZoomedMinRange, ZoomedMaxRange);

	pT_Dist->SetLineWidth(0);
	pT_Dist->SetMarkerColor(kBlack);
	pT_Dist->SetMarkerSize(3);
	pT_Dist->SetMarkerStyle(20);

	pT_Dist->Draw("P");

	cpT->cd(2);
	cpT->cd(2)->SetMargin(0.15,0.05,0.2,0.05);//lrbt

	TH1F* hFracAxis = new TH1F("hFracAxis","hFracAxis",pT_Weights->GetNbinsX(),pT_Weights->GetBinLowEdge(1),pT_Weights->GetXaxis()->GetBinUpEdge(pT_Weights->GetNbinsX()));

	hFracAxis->SetStats(false);
	hFracAxis->SetTitle("");
	hFracAxis->GetXaxis()->SetLabelSize(0.055);
	hFracAxis->GetXaxis()->SetTitle("p_{T} (GeV)");
	hFracAxis->GetXaxis()->CenterTitle();
	hFracAxis->GetXaxis()->SetTitleOffset(1.2);
	hFracAxis->GetXaxis()->SetLabelOffset(0.02);
	hFracAxis->GetXaxis()->SetTitleSize(0.055);
	hFracAxis->GetYaxis()->SetLabelSize(0.055);
	hFracAxis->GetYaxis()->SetTitle("Fraction");
	hFracAxis->GetYaxis()->CenterTitle();
	hFracAxis->GetYaxis()->SetTitleOffset(1.25);
	hFracAxis->GetYaxis()->SetTitleSize(0.055);

	//hFracAxis->GetYaxis()->SetRangeUser(1e3,1e7);
	hFracAxis->GetXaxis()->SetNdivisions(506);
	//hFracAxis->GetXaxis()->SetRangeUser(ZoomedMinRange, ZoomedMaxRange);

	graph_pri->SetLineWidth(0);
	graph_pri->SetMarkerColor(kRed+2);
	graph_pri->SetMarkerSize(3);
	graph_pri->SetMarkerStyle(20);

	graph_mat->SetLineWidth(0);
	graph_mat->SetMarkerColor(kGreen+3);
	graph_mat->SetMarkerSize(3);
	graph_mat->SetMarkerStyle(21);

	graph_sec->SetLineWidth(0);
	graph_sec->SetMarkerColor(kBlue+1);
	graph_sec->SetMarkerSize(3);
	graph_sec->SetMarkerStyle(22);

	graph_cont->SetLineWidth(0);
	graph_cont->SetMarkerColor(kOrange+1);
	graph_cont->SetMarkerSize(3);
	graph_cont->SetMarkerStyle(23);

	primAvg->SetLineColor(kRed+2);
	primAvg->SetLineWidth(3);
	primAvg->SetLineStyle(7);

    TLegend* LegendFrac = new TLegend(0.42,0.4,0.93,0.65);//lbrt
    LegendFrac->SetName(TString::Format("LegendFrac"));
    LegendFrac->SetTextSize(0.034);
    LegendFrac->AddEntry(primAvg,"Weighted average of primaries");
    LegendFrac->AddEntry(graph_pri,"Primary (#Lambda or #Sigma^{0}#rightarrow#Lambda)");
    LegendFrac->AddEntry(graph_sec,"Secondary (#Xi#rightarrow#Lambda)");
    LegendFrac->AddEntry(graph_mat,"Material #Lambda");
    LegendFrac->AddEntry(graph_cont,"Fake #Lambda");

	hFracAxis->Draw("axis");
	graph_pri->Draw("same,P");
	graph_mat->Draw("same,P");
	graph_sec->Draw("same,P");
	graph_cont->Draw("same,P");
	primAvg->Draw("same");
	LegendFrac->Draw("same");

	cpT->SaveAs("$FEMTO_OUTPUT/Scripts/pT_Distribution.png");
	gStyle->SetLineWidth(1);
	pT_Dist->SetLineWidth(pT_Dist->GetLineWidth()/2.5);
	primAvg->SetLineWidth(primAvg->GetLineWidth()/2.5);
	cpT->SaveAs("$FEMTO_OUTPUT/Scripts/pT_Distribution.pdf");
	gStyle->SetLineWidth(2.5);
	pT_Dist->SetLineWidth(pT_Dist->GetLineWidth()*2.5);
	primAvg->SetLineWidth(primAvg->GetLineWidth()*2.5);

	//delete [] pri;
	//delete [] mat;
	//delete [] sec;
	//delete [] cont;

	//delete InputFileName;
	//delete fOutPlot;
}
