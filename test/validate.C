

{
  bool DOGLOBAL     = true;
  bool DOUNWEIGHTED = false;

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
  
  TLegend* leg = new TLegend(0.34,0.72,0.70,0.89, NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04);
  
  TString names[] = {
    //
    //"tth_13tev_amcatnlo_pu20bx25_0" 
    //,
    //"tth_13tev_amcatnlo_pu20bx25_1"
    //,
    "ttjets_13tev_madgraph_pu20bx25_phys14_0ByDecay"
    //,
    //"ttjets_13tev_madgraph_pu20bx25_phys14_1"
  };
  TString samples[] = {
    //
    //"t#bar{t}H, reco CSV"
    //,
    //"t#bar{t}H, closure test"
    //,
    "t#bar{t}+jets, reco CSV"
    //,
    //"t#bar{t}+jets, closure test"
  };
  
  for( int it_name = 0 ; it_name < 4 ; ++it_name){

    TString name   = names[it_name];
    TString sample = samples[it_name];  

    if(DOUNWEIGHTED) sample = sample+" (unweighted)";

    TFile* fout = new TFile("test/validate_"+name+".root","RECREATE");
    fout->mkdir("csv");
    fout->mkdir("csv_ratio");
    fout->mkdir("pt");
    fout->mkdir("pt_ratio");
    fout->mkdir("eta");
    fout->mkdir("eta_ratio");
    fout->mkdir("dR");
    fout->mkdir("dR_ratio");
    fout->mkdir("Mjj");
    fout->mkdir("Mjj_ratio");
    fout->mkdir("global");
    fout->mkdir("global_ratio");
    
    //TFile* f = TFile::Open("test/btag_out_"+name+".root");
    //TFile* f = TFile::Open("/shome/bianchi/TTH-72X-heppy/CMSSW/src/TTH/MEIntegratorStandalone/test/btag_out_"+name+".root");
    TFile* f = TFile::Open("/scratch/bianchi/74X/btag_out_"+name+".root");

    TTree* t = (TTree*)f->Get("tree");
    TH1F*  h = new TH1F("h","",10,0,10);
    TH1F*  h_yields_inp   = new TH1F("h_yields_inp","",15,0,15);
    TH1F*  h_yields_rnd   = new TH1F("h_yields_rnd","",15,0,15);
    TH1F*  h_yields_ratio = new TH1F("h_yields_ratio","",15,0,15);
    TH1F*  h_frac_inp     = new TH1F("h_frac_inp","",15,0,15);
    TH1F*  h_frac_rnd     = new TH1F("h_frac_rnd","",15,0,15);
    TH1F*  h_frac_ratio   = new TH1F("h_frac_ratio","",15,0,15);
    

    string globals[3] = {"met","HT","ttCls"};
    //char* globals[2];

    for( int njet = 4 ; njet <=6 ; ++njet ){
      for( int i = 0 ; i <=4 ; ++i ){
	h->Reset();
	h->Sumw2();
	t->Draw("njet>>h",Form("njet==%d && pass[%d]",njet,i));
	double pass_err=0.; 
	int pass_entries = h->GetEntries();
	float pass = h->IntegralAndError(1, h->GetNbinsX(), pass_err); 
	
	h->Reset();
	t->Draw("njet>>h",Form("(njet==%d && pass_rnd[%d])*pcat[%d]", njet,i,i));
	double weight_err=0.; 
	int weight_entries = h->GetEntries();
	float weight = h->IntegralAndError(1, h->GetNbinsX(), weight_err); 

	float ratio = pass>0 ? (weight-pass)/pass  : 0.;
	printf("Jet: %d, Tag: %d\n", njet, i);
	printf("\tPASS:   %.0f +/- %.1f\n", pass, pass_err);
	printf("\tWEIGHT: %.0f +/- %.1f\n", weight, weight_err);
	printf("\tDiff: %.1f \n",   ratio*100);
	
	h_yields_ratio->GetXaxis()->SetBinLabel((njet-4)*5+i+1, Form("(%d,%d)", njet, i));
	h_yields_inp->SetBinContent((njet-4)*5+i+1, pass);
	h_yields_inp->SetBinError((njet-4)*5+i+1, pass_err);
	h_yields_rnd->SetBinContent((njet-4)*5+i+1, weight);
	h_yields_rnd->SetBinError((njet-4)*5+i+1, weight_err);
	h_frac_ratio->GetXaxis()->SetBinLabel((njet-4)*5+i+1, Form("(%d,%d)", njet, i));
	h_frac_inp->SetBinContent((njet-4)*5+i+1, pass_entries);
	h_frac_inp->SetBinError((njet-4)*5+i+1, sqrt(pass_entries));
	h_frac_rnd->SetBinContent((njet-4)*5+i+1, weight_entries);
	h_frac_rnd->SetBinError((njet-4)*5+i+1, sqrt(weight_entries));

	if(DOGLOBAL){	
	for(int gl = 0 ; gl < 3 ; ++gl){
	  string glo = globals[gl];
	  
	  int nbins = 10;
	  if( gl==2 ) nbins = 30;
	  float xmax = 400;
	  if( gl==1 ) xmax = 2000;
	  if( gl==2 ) xmax = 30;

	  TH1F* h_inp   = new TH1F(Form("h_%s_inp_%d_%d",  glo.c_str(),njet,i),Form("N_{jet}=%d, N_{tag}=%d; %s", njet, i, glo.c_str()), nbins,0, xmax);
	  TH1F* h_rnd   = new TH1F(Form("h_%s_rnd_%d_%d",  glo.c_str(),njet,i),Form("N_{jet}=%d, N_{tag}=%d; %s", njet, i, glo.c_str()), nbins,0, xmax);
	  TH1F* h_ratio = new TH1F(Form("h_%s_ratio_%d_%d",glo.c_str(),njet,i),Form("N_{jet}=%d, N_{tag}=%d; %s; (weighted-cut)/cut", njet, i, glo.c_str()), nbins, 0, xmax);
	  h_inp->Sumw2();
	  h_rnd->Sumw2();
	  h_ratio->Sumw2();
	  t->Draw(Form("%s>>h_%s_inp_%d_%d",glo.c_str(),glo.c_str(),njet, i), Form("njet==%d && pass[%d]",njet,i));
	  t->Draw(Form("%s>>h_%s_rnd_%d_%d",glo.c_str(),glo.c_str(),njet, i), Form("(njet==%d && pass_rnd[%d])*pcat[%d]", njet,i,i));
	  h_inp->SetLineColor(kBlue);
	  h_inp->SetLineWidth(3);
	  h_rnd->SetLineWidth(3);
	  h_rnd->SetLineColor(kRed);
	  h_ratio->SetLineColor(kBlack);      
	  fout->cd("global");
	  //h_inp->Write();
	  //h_rnd->Write();
	  c1->Clear();
	  c1->SetName(Form("h_%s_comp_%d_%d", glo.c_str(),njet,i));
	  leg->Clear();
	  c1->cd();
	  h_inp->SetMaximum( TMath::Max( h_inp->GetMaximum(),  h_rnd->GetMaximum() )*1.3 );
	  h_inp->Draw("HISTE");	
	  h_rnd->Draw("HISTESAME");	
	  float KS   = h_inp->KolmogorovTest( h_rnd );
	  float CHI2 = h_inp->Chi2Test( h_rnd, "UW" );
	  leg->AddEntry(h_inp, Form("Full cuts: %.0f +/- %.0f", pass, pass_err), "L");
	  leg->AddEntry(h_rnd, Form("Weighted : %.0f +/- %.0f", weight, weight_err), "L");
	  leg->SetHeader(sample+", "+TString(Form("KS: %.2f, #chi2: %.2f", KS, CHI2)));
	  leg->Draw();
	  c1->Write(); 
	  c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	fout->cd("global_ratio");
	h_rnd->Add(h_inp, -1.);
	h_ratio->Divide(h_rnd,h_inp,1.0,1.0);
	//h_ratio->Write();
	h_ratio->SetMaximum(+2.0);
	h_ratio->SetMinimum(-2.0);
	c1->Clear();
	c1->SetName(Form("h_%s_ratio_%d_%d", glo.c_str(),njet,i));
	c1->cd();
	leg->Clear();
	leg->SetHeader(sample);     
	leg->AddEntry(h_ratio, "Ratio" ,"L");	
	h_ratio->SetLineWidth(3); 
	h_ratio->SetLineColor(kBlue); 
	h_ratio->Draw();
	leg->Draw();
	TF1* line = new TF1("line","0",  h_ratio->GetXaxis()->GetXmin(),  h_ratio->GetXaxis()->GetXmax());
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->Draw("SAME");
	c1->Write();
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");
	delete h_inp; delete h_rnd; delete h_ratio; delete line;
      }


      for(int k = 0 ; k < 4 ; ++k){
	string title = "";
	if(k==0) title = "min for tagged jets";
	if(k==1) title = "min for untagged jets";
	if(k==2) title = "max for tagged jets";
	if(k==3) title = "max for untagged jets";

	TH1F* h_inp   = new TH1F(Form("h_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; #DeltaR %s", njet, i,title.c_str()),10,0,6);
	TH1F* h_rnd   = new TH1F(Form("h_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; #DeltaR %s", njet, i,title.c_str()),10,0,6);
	TH1F* h_ratio = new TH1F(Form("h_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; #DeltaR %s; (weighted-cut)/cut", njet, i,title.c_str()),10,0,6);
	h_inp->Sumw2();
	h_rnd->Sumw2();
	h_ratio->Sumw2();
	t->Draw(Form("corr_inp_%d[%d]>>h_inp_%d_%d_%d",i, k, njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("corr_rnd_%d[%d]>>h_rnd_%d_%d_%d",i, k, njet, i,k), Form("(njet==%d && pass_rnd[%d])*pcat[%d]", njet,i,i));
	if(DOUNWEIGHTED){
	  h_rnd->Reset();
	  t->Draw(Form("corr_rnd_%d[%d]>>h_rnd_%d_%d_%d",i, k, njet, i,k), Form("(njet==%d)", njet));
	  h_rnd->Scale(h_inp->Integral()/h_rnd->Integral());
	}

	h_inp->SetLineColor(kBlue);
	h_rnd->SetLineColor(kRed);
	h_inp->SetLineWidth(3);
	h_rnd->SetLineWidth(3);
	h_ratio->SetLineColor(kBlack);      
	fout->cd("dR");
	//h_inp->Write();
	//h_rnd->Write();
	c1->Clear();
	c1->SetName(Form("h_dR%d_comp_%d_%d", k,njet,i));
	leg->Clear();
	c1->cd();
	h_inp->SetMaximum( TMath::Max( h_inp->GetMaximum(),  h_rnd->GetMaximum() )*1.3 );
	h_inp->Draw("HISTE");	
	h_rnd->Draw("HISTESAME");	
	float KS   = h_inp->KolmogorovTest( h_rnd );
	float CHI2 = h_inp->Chi2Test( h_rnd, "UW" );
	leg->AddEntry(h_inp, Form("Full cuts: %.0f +/- %.0f", pass, pass_err), "L");
	leg->AddEntry(h_rnd, Form("Weighted : %.0f +/- %.0f", weight, weight_err), "L");
	leg->SetHeader(sample+", "+TString(Form("KS: %.2f, #chi2: %.2f", KS, CHI2)));
	leg->Draw();
	c1->Write(); 
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	fout->cd("dR_ratio");
	h_rnd->Add(h_inp, -1.);
	h_ratio->Divide(h_rnd,h_inp,1.0,1.0);
	//h_ratio->Write();
	h_ratio->SetMaximum(+2.0);
	h_ratio->SetMinimum(-2.0);
	c1->Clear();
	c1->SetName(Form("h_dR%d_ratio_%d_%d", k,njet,i));
	c1->cd();
	leg->Clear();
	leg->SetHeader(sample);     
	leg->AddEntry(h_ratio, "Ratio" ,"L");	
	h_ratio->SetLineWidth(3); 
	h_ratio->SetLineColor(kBlue); 
	h_ratio->Draw();
	leg->Draw();
	TF1* line = new TF1("line","0",  h_ratio->GetXaxis()->GetXmin(),  h_ratio->GetXaxis()->GetXmax());
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->Draw("SAME");
	c1->Write();
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	delete h_inp; delete h_rnd; delete h_ratio; delete line;
      }

      for(int k = 4 ; k < 8 ; ++k){
	string title = "";
	if(k==4) title = "min for tagged jets";
	if(k==5) title = "min for untagged jets";
	if(k==6) title = "max for tagged jets";
	if(k==7) title = "max for untagged jets";

	TH1F* h_inp   = new TH1F(Form("h_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; dijet mass %s", njet, i,title.c_str()),20,0, k<6 ? 200 : 600);
	TH1F* h_rnd   = new TH1F(Form("h_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; dijet mass %s", njet, i,title.c_str()),20,0, k<6 ? 200 : 600);
	TH1F* h_ratio = new TH1F(Form("h_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; dijet mass %s; (weighted-cut)/cut", njet, i,title.c_str()),20,0, k<6 ? 200 : 600);
	h_inp->Sumw2();
	h_rnd->Sumw2();
	h_ratio->Sumw2();
	t->Draw(Form("corr_inp_%d[%d]>>h_inp_%d_%d_%d",i, k, njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("corr_rnd_%d[%d]>>h_rnd_%d_%d_%d",i, k, njet, i,k), Form("(njet==%d && pass_rnd[%d])*pcat[%d]", njet,i,i));
	h_inp->SetLineColor(kBlue);
	h_rnd->SetLineColor(kRed);
	h_ratio->SetLineColor(kBlack);      
	h_inp->SetLineWidth(3);
	h_rnd->SetLineWidth(3);
	fout->cd("Mjj");
	//h_inp->Write();
	//h_rnd->Write();
	c1->Clear();
	c1->SetName(Form("h_Mjj%d_comp_%d_%d", k,njet,i));
	leg->Clear();
	c1->cd();
	h_inp->SetMaximum( TMath::Max( h_inp->GetMaximum(),  h_rnd->GetMaximum() )*1.3 );
	h_inp->Draw("HISTE");	
	h_rnd->Draw("HISTESAME");	
	float KS   = h_inp->KolmogorovTest( h_rnd );
        float CHI2 = h_inp->Chi2Test( h_rnd, "UW" );
	leg->AddEntry(h_inp, Form("Full cuts: %.0f +/- %.0f", pass, pass_err), "L");
	leg->AddEntry(h_rnd, Form("Weighted : %.0f +/- %.0f", weight, weight_err), "L");
	leg->SetHeader(sample+", "+TString(Form("KS: %.2f, #chi2: %.2f", KS, CHI2)));
	leg->Draw();
	c1->Write(); 
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	fout->cd("Mjj_ratio");
	h_rnd->Add(h_inp, -1.);
	h_ratio->Divide(h_rnd,h_inp,1.0,1.0);
	//h_ratio->Write();
	h_ratio->SetMaximum(+2.0);
	h_ratio->SetMinimum(-2.0);
	c1->Clear();
	c1->SetName(Form("h_Mjj%d_ratio_%d_%d", k,njet,i));
	c1->cd();
	leg->Clear();
	leg->SetHeader(sample);     
	leg->AddEntry(h_ratio, "Ratio" ,"L");	
	h_ratio->SetLineWidth(3); 
	h_ratio->SetLineColor(kBlue); 
	h_ratio->Draw();
	leg->Draw();
	TF1* line = new TF1("line","0",  h_ratio->GetXaxis()->GetXmin(),  h_ratio->GetXaxis()->GetXmax());
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->Draw("SAME");
	c1->Write();
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	delete h_inp; delete h_rnd; delete h_ratio; delete line;
      }
	}

      for(int k = 0 ; k <6 ; ++k){
	if(k>=njet) continue;
	TH1F* h_csv_inp   = new TH1F(Form("h_csv_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; CSV", njet, i, k),15,0,1);
	TH1F* h_csv_rnd   = new TH1F(Form("h_csv_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; CSV", njet, i, k),15,0,1);
	TH1F* h_csv_ratio = new TH1F(Form("h_csv_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; CSV; (weighted-cut)/cut", njet, i, k),15,0,1);

	TH1F* h_pt_inp   = new TH1F(Form("h_pt_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; p_{T}", njet, i, k),20,0,400);
	TH1F* h_pt_rnd   = new TH1F(Form("h_pt_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; p_{T}", njet, i, k),20,0,400);
	TH1F* h_pt_ratio = new TH1F(Form("h_pt_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; p_{T}; (weighted-cut)/cut", njet, i, k),20,0,400);

	TH1F* h_eta_inp   = new TH1F(Form("h_eta_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; #eta", njet, i, k),27,-2.7,2.7);
	TH1F* h_eta_rnd   = new TH1F(Form("h_eta_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; #eta", njet, i, k),27,-2.7,2.7);
	TH1F* h_eta_ratio = new TH1F(Form("h_eta_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; #eta", njet, i, k),27,-2.7,2.7);

	h_csv_inp->Sumw2();
	h_csv_rnd->Sumw2();
	h_csv_ratio->Sumw2();
	t->Draw(Form("csv_inp_%d[%d]>>h_csv_inp_%d_%d_%d",i,k,njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("csv_rnd_%d[%d]>>h_csv_rnd_%d_%d_%d",i,k,njet, i,k), Form("(njet==%d && pass_rnd[%d])*pcat[%d]", njet,i,i));
	if(DOUNWEIGHTED){
	  h_csv_rnd->Reset();
	  t->Draw(Form("csv_rnd_%d[%d]>>h_csv_rnd_%d_%d_%d",i,k,njet, i,k), Form("(njet==%d)", njet));
	  h_csv_rnd->Scale(h_csv_inp->Integral()/h_csv_rnd->Integral());
	}

	h_pt_inp->Sumw2();
	h_pt_rnd->Sumw2();
	h_pt_ratio->Sumw2();
	t->Draw(Form("pt[%d]>>h_pt_inp_%d_%d_%d",k,njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("pt[%d]>>h_pt_rnd_%d_%d_%d",k,njet, i,k), Form("(njet==%d && pass_rnd[%d])*pcat[%d]", njet,i,i));

	h_eta_inp->Sumw2();
	h_eta_rnd->Sumw2();
	h_eta_ratio->Sumw2();
	t->Draw(Form("eta[%d]>>h_eta_inp_%d_%d_%d",k,njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("eta[%d]>>h_eta_rnd_%d_%d_%d",k,njet, i,k), Form("(njet==%d && pass_rnd[%d])*pcat[%d]", njet,i,i));

	fout->cd();            
	h_csv_inp->SetLineColor(kBlue);
	h_csv_rnd->SetLineColor(kRed);
	h_csv_inp->SetLineWidth(3);
	h_csv_rnd->SetLineWidth(3);

	h_pt_inp->SetLineColor(kBlue);
	h_pt_rnd->SetLineColor(kRed);
	h_pt_inp->SetLineWidth(3);
	h_pt_rnd->SetLineWidth(3);

	h_eta_inp->SetLineColor(kBlue);
	h_eta_rnd->SetLineColor(kRed);
	h_eta_inp->SetLineWidth(3);
	h_eta_rnd->SetLineWidth(3);


	fout->cd("pt");
	//h_pt_inp->Write();
	//h_pt_rnd->Write();
	c1->Clear();
	c1->SetName(Form("h_JetPt%d_comp_%d_%d", k,njet,i));
	leg->Clear();
	c1->cd();
	h_pt_inp->SetMaximum( TMath::Max( h_pt_inp->GetMaximum(),  h_pt_rnd->GetMaximum() )*1.3 );
	h_pt_inp->Draw("HISTE");	
	h_pt_rnd->Draw("HISTESAME");	
	float KS   = h_pt_inp->KolmogorovTest( h_pt_rnd );
	float CHI2 = h_pt_inp->Chi2Test( h_pt_rnd, "UW" );
	leg->AddEntry(h_pt_inp, Form("Full cuts: %.0f +/- %.0f", pass, pass_err), "L");
	leg->AddEntry(h_pt_rnd, Form("Weighted : %.0f +/- %.0f", weight, weight_err), "L");
	leg->SetHeader(sample+", "+TString(Form("KS: %.2f, #chi2: %.2f", KS, CHI2)));
	leg->Draw();
	c1->Write(); 
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	fout->cd("eta");
	//h_eta_inp->Write();
	//h_eta_rnd->Write();
	c1->Clear();
	c1->SetName(Form("h_JetEta%d_comp_%d_%d", k,njet,i));
	leg->Clear();
	c1->cd();
	h_eta_inp->SetMaximum( TMath::Max( h_eta_inp->GetMaximum(),  h_eta_rnd->GetMaximum() )*1.3 );
	h_eta_inp->Draw("HISTE");	
	h_eta_rnd->Draw("HISTESAME");	
	KS   = h_eta_inp->KolmogorovTest( h_eta_rnd );
	CHI2 = h_eta_inp->Chi2Test( h_eta_rnd, "UW" );
	leg->AddEntry(h_eta_inp, Form("Full cuts: %.0f +/- %.0f", pass, pass_err), "L");
	leg->AddEntry(h_eta_rnd, Form("Weighted : %.0f +/- %.0f", weight, weight_err), "L");
	leg->SetHeader(sample+", "+TString(Form("KS: %.2f, #chi2: %.2f", KS, CHI2)));
	leg->Draw();
	c1->Write(); 
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");


	fout->cd("csv");
	//h_csv_inp->Write();
	//h_csv_rnd->Write();
	c1->Clear();
	c1->SetName(Form("h_JetCSV%d_comp_%d_%d", k,njet,i));
	leg->Clear();
	c1->cd();
	h_csv_inp->SetMaximum( TMath::Max( h_csv_inp->GetMaximum(),  h_csv_rnd->GetMaximum() )*1.3 );
	h_csv_inp->Draw("HISTE");	
	h_csv_rnd->Draw("HISTESAME");	
	KS   = h_csv_inp->KolmogorovTest( h_csv_rnd );
	CHI2 = h_csv_inp->Chi2Test( h_csv_rnd, "UW" );
	leg->AddEntry(h_csv_inp, Form("Full cuts: %.0f +/- %.0f", pass, pass_err), "L");
	leg->AddEntry(h_csv_rnd, Form("Weighted : %.0f +/- %.0f", weight, weight_err), "L");
	leg->SetHeader(sample+", "+TString(Form("KS: %.2f, #chi2: %.2f", KS, CHI2)));
	leg->Draw();
	c1->Write(); 
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	fout->cd("csv_ratio");
	h_csv_rnd->Add(h_csv_inp,-1.);
	h_csv_ratio->Divide(h_csv_rnd,h_csv_inp,1.0,1.0);
	//h_csv_ratio->Write();
	h_csv_ratio->SetMaximum(+2.0);
	h_csv_ratio->SetMinimum(-2.0);
	c1->Clear();
	c1->SetName(Form("h_JetCSV%d_ratio_%d_%d", k,njet,i));
	c1->cd();
	leg->Clear();
	leg->SetHeader(sample);     
	leg->AddEntry(h_csv_ratio, "Ratio" ,"L");	
	h_csv_ratio->SetLineWidth(3); 
	h_csv_ratio->SetLineColor(kBlue); 
	h_csv_ratio->Draw();
	leg->Draw();
	TF1* line = new TF1("line","0",  h_csv_ratio->GetXaxis()->GetXmin(),  h_csv_ratio->GetXaxis()->GetXmax());
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->Draw("SAME");
	c1->Write();
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	fout->cd("pt_ratio");
	h_pt_rnd->Add(h_pt_inp,-1.);
	h_pt_ratio->Divide(h_pt_rnd,h_pt_inp,1.0,1.0);
	//h_pt_pt_ratio->Write();
	h_pt_ratio->SetMaximum(+2.0);
	h_pt_ratio->SetMinimum(-2.0);
	c1->Clear();
	c1->SetName(Form("h_JetPt%d_ratio_%d_%d", k,njet,i));
	c1->cd();
	leg->Clear();
	leg->SetHeader(sample);     
	leg->AddEntry(h_pt_ratio, "Ratio" ,"L");	
	h_pt_ratio->SetLineWidth(3); 
	h_pt_ratio->SetLineColor(kBlue); 
	h_pt_ratio->Draw();
	leg->Draw();
	line = new TF1("line","0",  h_pt_ratio->GetXaxis()->GetXmin(),  h_pt_ratio->GetXaxis()->GetXmax());
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->Draw("SAME");
	c1->Write();
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	fout->cd("eta_ratio");
	h_eta_rnd->Add(h_eta_inp,-1.);
	h_eta_ratio->Divide(h_eta_rnd,h_eta_inp,1.0,1.0);
	//h_eta_eta_ratio->Write();
	h_eta_ratio->SetMaximum(+2.0);
	h_eta_ratio->SetMinimum(-2.0);
	c1->Clear();
	c1->SetName(Form("h_JetEta%d_ratio_%d_%d", k,njet,i));
	c1->cd();
	leg->Clear();
	leg->SetHeader(sample);     
	leg->AddEntry(h_eta_ratio, "Ratio" ,"L");	
	h_eta_ratio->SetLineWidth(3); 
	h_eta_ratio->SetLineColor(kBlue); 
	h_eta_ratio->Draw();
	leg->Draw();
	line = new TF1("line","0",  h_eta_ratio->GetXaxis()->GetXmin(),  h_eta_ratio->GetXaxis()->GetXmax());
	line->SetLineWidth(3);
	line->SetLineColor(kBlack);
	line->SetLineStyle(kDashed);
	line->Draw("SAME");
	c1->Write();
	c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

	delete h_pt_inp; delete h_pt_rnd; delete h_pt_ratio;
	delete h_eta_inp; delete h_eta_rnd; delete h_eta_ratio;
	delete h_csv_inp; delete h_csv_rnd; delete h_csv_ratio;

	f->cd();
      }
    }
    }

    fout->cd();    
    h_yields_rnd->Add(h_yields_inp, -1.);
    h_yields_ratio->Divide(h_yields_rnd,h_yields_inp,1.0,1.0);
    h_yields_ratio->SetMaximum(+.3);
    h_yields_ratio->SetMinimum(-.3);
    c1->Clear();
    c1->SetName("yields");
    c1->cd();
    leg->Clear();
    leg->SetHeader(sample);     
    leg->AddEntry(h_yields_ratio, "#frac{N_{weight}-N_{cut}}{N_{cut}}" ,"L");	
    h_yields_ratio->SetLineWidth(3); 
    h_yields_ratio->SetLineColor(kBlue); 
    h_yields_ratio->Draw();
    leg->Draw();
    TF1* line = new TF1("line","0",  h_yields_ratio->GetXaxis()->GetXmin(),  h_yields_ratio->GetXaxis()->GetXmax());
    line->SetLineWidth(3);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw("SAME");
    c1->Write();
    c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");

    fout->cd();    
    h_frac_ratio->Divide(h_frac_inp,h_frac_rnd,1.0,1.0);
    h_frac_ratio->SetMaximum(0.7);
    h_frac_ratio->SetMinimum(0.0);
    c1->Clear();
    c1->SetName("fraction");
    c1->cd();
    leg->Clear();
    leg->SetHeader(sample);     
    leg->AddEntry(h_frac_ratio, "Fraction of passing events" ,"L");	
    h_frac_ratio->SetLineWidth(3); 
    h_frac_ratio->SetLineColor(kBlue); 
    h_frac_ratio->Draw();
    leg->Draw();
    line = new TF1("line","0",  h_frac_ratio->GetXaxis()->GetXmin(),  h_frac_ratio->GetXaxis()->GetXmax());
    line->SetLineWidth(3);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw("SAME");
    c1->Write();
    c1->SaveAs("./test/plots/"+TString(c1->GetName())+"_"+name+".png");


    fout->cd();
    fout->Close();
    f->Close();

    delete f; delete fout;
  }

}
