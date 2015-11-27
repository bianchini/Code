
{

  // X axis
  int binX_L = 1;
  int binX_H = -1;

  // Y axis
  int binY_L = 1;
  int binY_H = -1;

  TFile* f = TFile::Open("csv.root", "READ");

  TCanvas *c1  = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(1,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");
  c1->Divide(1,2);

  TLegend* leg = new TLegend(0.34,0.72,0.70,0.89, NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);

  TString title = "TTbar powheg simulation, p_{T}>20 GeV, |#eta|<2.5";
  std::vector<TString> histos_light = {
    "csv_l_pt_eta",
    "csv_g_pt_eta",
    "csv_s_pt_eta",
    "csv_ud_pt_eta",
    "csv_pu_pt_eta",
  };
  std::vector<TString> histos_b = {
    "csv_b_pt_eta",
    "csv_1b_pt_eta",
    "csv_2b_pt_eta"
  };
  std::vector<TString> histos_c = {
    "csv_c_pt_eta",
    "csv_1c_pt_eta",
    "csv_2c_pt_eta"
  };
  std::vector<TString> histos_l_dR = {
    "csv_l_pt_eta",
    "csv_l_dR_0p0_1p0_pt_eta",
    "csv_l_dR_1p0_2p0_pt_eta",
    "csv_l_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_g_dR = {
    "csv_g_pt_eta",
    "csv_g_dR_0p0_1p0_pt_eta",
    "csv_g_dR_1p0_2p0_pt_eta",
    "csv_g_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_ud_dR = {
    "csv_ud_pt_eta",
    "csv_ud_dR_0p0_1p0_pt_eta",
    "csv_ud_dR_1p0_2p0_pt_eta",
    "csv_ud_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_s_dR = {
    "csv_s_pt_eta",
    "csv_s_dR_0p0_1p0_pt_eta",
    "csv_s_dR_1p0_2p0_pt_eta",
    "csv_s_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_b_dR = {
    "csv_b_pt_eta",
    "csv_b_dR_0p0_1p0_pt_eta",
    "csv_b_dR_1p0_2p0_pt_eta",
    "csv_b_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_c_dR = {
    "csv_c_pt_eta",
    "csv_c_dR_0p0_1p0_pt_eta",
    "csv_c_dR_1p0_2p0_pt_eta",
    "csv_c_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_1b_dR = {
    "csv_1b_pt_eta",
    "csv_1b_dR_0p0_1p0_pt_eta",
    "csv_1b_dR_1p0_2p0_pt_eta",
    "csv_1b_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_2b_dR = {
    "csv_2b_pt_eta",
    "csv_2b_dR_0p0_1p0_pt_eta",
    "csv_2b_dR_1p0_2p0_pt_eta",
    "csv_2b_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_1c_dR = {
    "csv_1c_pt_eta",
    "csv_1c_dR_0p0_1p0_pt_eta",
    "csv_1c_dR_1p0_2p0_pt_eta",
    "csv_1c_dR_2p0_Inf_pt_eta"
  };
  std::vector<TString> histos_2c_dR = {
    "csv_2c_pt_eta",
    "csv_2c_dR_0p0_1p0_pt_eta",
    "csv_2c_dR_1p0_2p0_pt_eta",
    "csv_2c_dR_2p0_Inf_pt_eta"
  };

  std::vector<TString> histos = histos_b_dR;



  TH1D* h0 = 0;
  for(int i = 0 ; i < histos.size() ; ++i){
    TH3D* h3 = (TH3D*)f->Get(histos[i]);
    if(h3==0){
      cout << "No histo for " << histos[i] << endl;
      continue;      
    }

  if(binX_H<0) binX_H =  h3->GetXaxis()->GetNbins();
  if(binY_H<0) binY_H =  h3->GetYaxis()->GetNbins();
  TH1D* h1 = h3->ProjectionZ(Form("%d_pz",i), binX_L, binX_H, binY_L, binY_H );
    h1->Rebin(2);
    h1->SetLineWidth(2);
    h1->SetLineColor(i+1);
    leg->AddEntry(h1, histos[i], "L");
    if(i==0){
      h1->SetTitle(title);
      c1->cd(1);
      h1->SetMaximum(h1->GetMaximum()*1.5);
      h1->DrawNormalized("HIST");
      c1->cd(2);
      h0 = (TH1D*)h1->Clone(Form("%s_copy",h1->GetName()));
      h0->Scale(1./h0->Integral());
      TF1* line = new TF1("line", "1.0", 0,1);
      line->SetLineColor(i+1);
      line->SetLineWidth(2);
      line->Draw();
    }
    else{
      c1->cd(1);
      h1->DrawNormalized("HISTESAME");
      c1->cd(2);
      TH1D* h1_copy = (TH1D*)h1->Clone(Form("%s_copy",h1->GetName()));
      h1_copy->Scale(1./h1_copy->Integral());
      h1_copy->Divide(h0);
      if(i==1){
	h1_copy->SetTitle("Ratio to first histogram");
	h1_copy->SetMinimum(0);
	h1_copy->SetMaximum(4);
	h1_copy->Draw("SAME");
      }
      else{
	h1_copy->Draw("SAME");
      }
    }
  }

  c1->cd(1);
  leg->Draw();
  

}
