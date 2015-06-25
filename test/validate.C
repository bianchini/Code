

{

  //TString name = "tth_13tev_amcatnlo_pu20bx25_0";
  //TString name = "tth_13tev_amcatnlo_pu20bx25_1";
  TString name = "ttjets_13tev_madgraph_pu20bx25_phys14_0";
  //TString name = "ttjets_13tev_madgraph_pu20bx25_phys14_1";

  TFile* fout = new TFile("test/validate_"+name+".root","RECREATE");
  fout->mkdir("csv");
  fout->mkdir("csv_ratio");
  fout->mkdir("pt");
  fout->mkdir("pt_ratio");
  fout->mkdir("dR");
  fout->mkdir("dR_ratio");
  fout->mkdir("Mjj");
  fout->mkdir("Mjj_ratio");
  
  TFile* f = TFile::Open("test/btag_out_"+name+".root");
  TTree* t = (TTree*)f->Get("tree");
  TH1F*  h = new TH1F("h","",10,0,10);


  for( int njet = 4 ; njet <=6 ; ++njet ){
    for( int i = 0 ; i <=4 ; ++i ){
      h->Reset();
      h->Sumw2();
      t->Draw("njet>>h",Form("njet==%d && pass[%d]",njet,i));
      double pass_err=0.; 
      float pass = h->IntegralAndError(1, h->GetNbinsX(), pass_err); 
      
      h->Reset();
      t->Draw("njet>>h",Form("(njet==%d)*pcat[%d]", njet,i));
      double weight_err=0.; 
      float weight = h->IntegralAndError(1, h->GetNbinsX(), weight_err); 

      float ratio = pass>0 ? (weight-pass)/pass  : 0.;
      printf("Jet: %d, Tag: %d\n", njet, i);
      printf("\tPASS:   %.0f +/- %.1f\n", pass, pass_err);
      printf("\tWEIGHT: %.0f +/- %.1f\n", weight, weight_err);
      printf("\tDiff: %.1f\%\n",   ratio*100);

      for(int k = 0 ; k <4 ; ++k){
	string title = "";
	if(k==0) title = "min for tagged jets";
	if(k==1) title = "min for untagged jets";
	if(k==2) title = "max for tagged jets";
	if(k==3) title = "max for untagged jets";

	TH1F* h_inp   = new TH1F(Form("h_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; #DeltaR %s", njet, i,title.c_str()),20,0,6);
	TH1F* h_rnd   = new TH1F(Form("h_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; #DeltaR %s", njet, i,title.c_str()),20,0,6);
	TH1F* h_ratio = new TH1F(Form("h_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; #DeltaR %s", njet, i,title.c_str()),20,0,6);
	h_inp->Sumw2();
	h_rnd->Sumw2();
	h_ratio->Sumw2();
	t->Draw(Form("corr_inp_%d[%d]>>h_inp_%d_%d_%d",i, k, njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("corr_rnd_%d[%d]>>h_rnd_%d_%d_%d",i, k, njet, i,k), Form("(njet==%d)*pcat[%d]", njet,i));
	fout->cd();            
	//for( int b = 1 ; b <= h_ratio->GetNbinsX() ; ++b ){
	//  float x =  (h_rnd->GetBinContent(b)>0.) ? h_rnd->GetBinError(b)/h_rnd->GetBinContent(b) : 0.  ;
	//  float y =  (h_inp->GetBinContent(b)>0.) ? h_inp->GetBinError(b)/h_inp->GetBinContent(b) : 0.  ;
	//  h_ratio->SetBinError(b,  h_ratio->GetBinContent(b)*sqrt( x*x + y*y) );
	//}

	h_inp->SetLineColor(kBlue);
	h_rnd->SetLineColor(kRed);
	h_ratio->SetLineColor(kBlack);      
	fout->cd("dR");
	h_inp->Write();
	h_rnd->Write();
	fout->cd("dR_ratio");
	h_rnd->Add(h_inp, -1.);
	h_ratio->Divide(h_rnd,h_inp,1.0,1.0);
	h_ratio->Write();
	delete h_inp; delete h_rnd; delete h_ratio;
      }

      for(int k = 4 ; k <8 ; ++k){
	string title = "";
	if(k==4) title = "min for tagged jets";
	if(k==5) title = "min for untagged jets";
	if(k==6) title = "max for tagged jets";
	if(k==7) title = "max for untagged jets";

	TH1F* h_inp   = new TH1F(Form("h_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; dijet mass %s", njet, i,title.c_str()),20,0, k<6 ? 200 : 600);
	TH1F* h_rnd   = new TH1F(Form("h_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; dijet mass %s", njet, i,title.c_str()),20,0, k<6 ? 200 : 600);
	TH1F* h_ratio = new TH1F(Form("h_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d; dijet mass %s", njet, i,title.c_str()),20,0, k<6 ? 200 : 600);
	h_inp->Sumw2();
	h_rnd->Sumw2();
	h_ratio->Sumw2();
	t->Draw(Form("corr_inp_%d[%d]>>h_inp_%d_%d_%d",i, k, njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("corr_rnd_%d[%d]>>h_rnd_%d_%d_%d",i, k, njet, i,k), Form("(njet==%d)*pcat[%d]", njet,i));
	fout->cd();            
	//for( int b = 1 ; b <= h_ratio->GetNbinsX() ; ++b ){
	//  float x =  (h_rnd->GetBinContent(b)>0.) ? h_rnd->GetBinError(b)/h_rnd->GetBinContent(b) : 0.  ;
	//  float y =  (h_inp->GetBinContent(b)>0.) ? h_inp->GetBinError(b)/h_inp->GetBinContent(b) : 0.  ;
	//  h_ratio->SetBinError(b,  h_ratio->GetBinContent(b)*sqrt( x*x + y*y) );
	//}

	h_inp->SetLineColor(kBlue);
	h_rnd->SetLineColor(kRed);
	h_ratio->SetLineColor(kBlack);      
	fout->cd("Mjj");
	h_inp->Write();
	h_rnd->Write();
	fout->cd("Mjj_ratio");
	h_rnd->Add(h_inp, -1.);
	h_ratio->Divide(h_rnd,h_inp,1.0,1.0);
	h_ratio->Write();
	delete h_inp; delete h_rnd; delete h_ratio;
      }

      for(int k = 0 ; k <6 ; ++k){
	if(k>=njet) continue;
	TH1F* h_csv_inp   = new TH1F(Form("h_csv_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; CSV", njet, i, k),10,0,1);
	TH1F* h_csv_rnd   = new TH1F(Form("h_csv_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; CSV", njet, i, k),10,0,1);
	TH1F* h_csv_ratio = new TH1F(Form("h_csv_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; CSV", njet, i, k),10,0,1);

	TH1F* h_pt_inp   = new TH1F(Form("h_pt_inp_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; PT", njet, i, k),20,0,300);
	TH1F* h_pt_rnd   = new TH1F(Form("h_pt_rnd_%d_%d_%d",  njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; PT", njet, i, k),20,0,300);
	TH1F* h_pt_ratio = new TH1F(Form("h_pt_ratio_%d_%d_%d",njet,i,k),Form("N_{jet}=%d, N_{tag}=%d, Jet %d; PT", njet, i, k),20,0,300);

	h_csv_inp->Sumw2();
	h_csv_rnd->Sumw2();
	h_csv_ratio->Sumw2();
	t->Draw(Form("csv_inp_%d[%d]>>h_csv_inp_%d_%d_%d",i,k,njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("csv_rnd_%d[%d]>>h_csv_rnd_%d_%d_%d",i,k,njet, i,k), Form("(njet==%d)*pcat[%d]", njet,i));

	h_pt_inp->Sumw2();
	h_pt_rnd->Sumw2();
	h_pt_ratio->Sumw2();
	t->Draw(Form("pt[%d]>>h_pt_inp_%d_%d_%d",k,njet, i,k), Form("njet==%d && pass[%d]",njet,i));
	t->Draw(Form("pt[%d]>>h_pt_rnd_%d_%d_%d",k,njet, i,k), Form("(njet==%d)*pcat[%d]", njet,i));

	fout->cd();            
	h_csv_inp->SetLineColor(kBlue);
	h_csv_rnd->SetLineColor(kRed);
	h_csv_ratio->SetLineColor(kBlack);
	h_pt_inp->SetLineColor(kBlue);
	h_pt_rnd->SetLineColor(kRed);
	h_pt_ratio->SetLineColor(kBlack);

	fout->cd("pt");
	h_pt_inp->Write();
	h_pt_rnd->Write();

	fout->cd("csv");
	h_csv_inp->Write();
	h_csv_rnd->Write();

	fout->cd("csv_ratio");
	h_csv_rnd->Add(h_csv_inp,-1.);
	h_csv_ratio->Divide(h_csv_rnd,h_csv_inp,1.0,1.0);
	h_csv_ratio->Write();

	fout->cd("pt_ratio");
	h_pt_rnd->Add(h_pt_inp,-1.);
	h_pt_ratio->Divide(h_pt_rnd,h_pt_inp,1.0,1.0);
	h_pt_ratio->Write();
	
	delete h_pt_inp; delete h_pt_rnd; delete h_pt_ratio;
	delete h_csv_inp; delete h_csv_rnd; delete h_csv_ratio;

	f->cd();
      }
    }
  }

  fout->cd();
  fout->Close();

}
