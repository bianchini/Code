

{

  //TString name = "tth_13tev_amcatnlo_pu20bx25_0";
  //TString name = "tth_13tev_amcatnlo_pu20bx25_1";
  //TString name = "ttjets_13tev_madgraph_pu20bx25_phys14_0";
  TString name = "ttjets_13tev_madgraph_pu20bx25_phys14_1";

  TFile* fout = new TFile("test/validate_"+name+".root","RECREATE");
  fout->mkdir("csv");
  fout->mkdir("csv_ratio");
  fout->mkdir("pt");
  fout->mkdir("pt_ratio");
  
  TFile* f = TFile::Open("test/btag_out_"+name+".root");
  TTree* t = (TTree*)f->Get("tree");
  TH1F*  h = new TH1F("h","",10,0,10);


  for( int njet = 4 ; njet <=6 ; ++njet ){
    for( int i = 0 ; i <=4 ; ++i ){
      h->Reset();
      t->Draw("njet>>h",Form("njet==%d && pass[%d]",njet,i));
      float pass = h->Integral(); 
      
      h->Reset();
      t->Draw("njet>>h",Form("(njet==%d)*pcat[%d]", njet,i));
      float weight = h->Integral(); 

      float ratio = pass>0 ? (weight-pass)/pass  : 0.;
      printf("Jet: %d, Tag: %d\n", njet, i);
      printf("\tPASS:   %.0f +/- %.0f\n", pass, sqrt(pass));
      printf("\tWEIGHT: %.0f\n", weight);
      printf("\tDiff: %.1f\%\n",   ratio*100);
  
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
	h_csv_ratio->Divide(h_csv_rnd,h_csv_inp,1.0,1.0);
	h_csv_inp->SetLineColor(kBlue);
	h_csv_rnd->SetLineColor(kRed);
	h_csv_ratio->SetLineColor(kBlack);
	h_pt_ratio->Divide(h_pt_rnd,h_pt_inp,1.0,1.0);
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
	h_csv_ratio->Write();

	fout->cd("pt_ratio");
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
