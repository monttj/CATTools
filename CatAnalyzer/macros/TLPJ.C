#include "MySample.C"

void TLPJ(){

  const double LUMI = 1263;
 
  vector<MySample> MySamples;
 
  MySamples.push_back( MySample("Data", "vallot_SingleMuon.root", LUMI) ) ; 
  MySamples.push_back( MySample("TTbar", "vallot_TTbar.root", 831.7, 42784971.0, 2, LUMI) ) ; 
  MySamples.push_back( MySample("WJets", "vallot_WJets.root", 61526.7, 23577656.0, 3, LUMI) ) ; 
  
  const int n = MySamples.size();

  TTree * tree[n];

  for(int i = 0; i < n; i++){
    tree[i] = (TTree*) MySamples[i].GetFile()->Get("ttbarSingleLepton/tree");
  }

  TH1 * h_MET[n] ;

  for(int i = 0; i < n ; i++){
    h_MET[i] = new TH1F(Form("h_MET_%d",i),"MET distribution",20,0,200);
    tree[i]->Project(Form("h_MET_%d",i),"MET","lepton_iso < 0.12");
  }

  for(int i = 0; i < n; i++){
    h_MET[i]->Scale( MySamples[i].GetScale() );
  }

  TCanvas * c1 = new TCanvas();
  c1->cd();
  THStack *hs = new THStack("hs","");
 
  for(int i = 1; i < n; i++){
    hs->Add(h_MET[i]);
    h_MET[i]->SetFillColor( MySamples[i].GetColor() );
  }

  hs->Draw();

  cout << "draw data" << endl;

  h_MET[0]->Draw("SameP");
  h_MET[0]->SetMarkerSize(1.2);
  h_MET[0]->SetMarkerStyle(20);

}
