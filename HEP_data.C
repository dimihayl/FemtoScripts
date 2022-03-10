

//CollEnergy: in GeV !!!
void Convert_TH1F_Yaml(const TString OutputFileName, TGrap* hInput,
  const double& CollEnergy, const TString xaxis, const TString yaxis, , const double& xmin, const double& xmax){
  ofstream myfile (OutputFileName.Data(),ios::out);
  myfile << "dependent_variables:" << endl;
  myfile << TString::Format("- header: {name: %s}",yaxis.Data()).Data() << endl;
  myfile << "  qualifiers:" << endl;
  myfile << TString::Format("  - {name: SQRT(S), units: GeV, value: '%.1f'}",CollEnergy).Data() << endl;
  myfile << "  values:" << endl;
  double SystErr;
  double StatErr;
  double xValue;
  double yValue;

  for(unsigned uBin=1; uBin<hInput->GetNbinsX(); uBin++){
    xValue = hInput->GetBinCenter(uBin);
    yValue = hInput->GetBinContent(uBin);
    if(xValue>xmax) break;
    if(xValue<xmin) break;
    SystErr = Tgraph_syserror_LL_ALAL->GetErrorY(uBin-1);
    StatErr = DataHisto->GetBinError(uBin);
    myfile << "  - errors:" << endl;
    myfile << "    - {label: stat, symerror: "<<StatErr<<"}" << endl;
    myfile << "    - {label: sys, symerror: "<<SystErr<<"}" << endl;
    myfile << "    value: "<<CkValue << endl;
  }
  myfile << "independent_variables:" << endl;
  myfile << "- header: {name: k* (MeV/c)}" << endl;
  myfile << "  values:" << endl;
  for(unsigned uBin=1; uBin<DataHisto->GetNbinsX(); uBin++){
      kValue = DataHisto->GetBinCenter(uBin);
      if(kValue>MAX_PLOT_K) break;
      myfile << "  - {value: "<<kValue<<"}" << endl;
  }
  myfile.close();
}

void Convert_TGraph_Yaml(){





}

void Convert_TGraphErrors_Yaml(){

}

void Convert_TGraphAsymmErrors_Yaml(){

}


void HEP_data(){


}
