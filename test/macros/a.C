

{

  TFile* f = TFile::Open("/scratch/bianchi/74X/btag_out_ttjets_13tev_madgraph_pu20bx25_phys14_0.root","READ");//btag_out_tth_13tev_amcatnlo_pu20bx25_0.root"); //btag_out_ttjets_13tev_madgraph_pu20bx25_phys14_0.root","READ");
  TTree* t = (TTree*)f->Get("tree");

  float pt [6];
  float eta[6];
  float csv[6];
  int pdgid[6];
  int  mcid[6];
  int njet;
  int ttCls_;
  int event_;

  t->SetBranchAddress("event",  &event_);
  t->SetBranchAddress("njet",  &njet);
  t->SetBranchAddress("pt",    pt);
  t->SetBranchAddress("eta",   eta);
  t->SetBranchAddress("csv_inp_4",   csv);
  t->SetBranchAddress("pdgid", pdgid);
  t->SetBranchAddress("mcid",  mcid);
  t->SetBranchAddress("ttCls",  &ttCls_);

  Long64_t nentries = t->GetEntries();
  cout << "Total entries: " << nentries << endl;


  TFile* out  = TFile::Open("a.root","RECREATE");
  out->cd();
  TTree* tout = new TTree("tree","tree");
  float pt1, pt2, eta1, eta2, csv1, csv2;
  int mcid1, mcid2, pdgid1, pdgid2;
  int ttCls;
  int event;
  tout->Branch("event", &event, "event/I");
  tout->Branch("pt1", &pt1, "pt1/F");
  tout->Branch("pt2", &pt2, "pt2/F");
  tout->Branch("csv1", &csv1, "csv1/F");
  tout->Branch("csv2", &csv2, "csv2/F");
  tout->Branch("eta1", &eta1, "eta1/F");
  tout->Branch("eta2", &eta2, "eta2/F");
  tout->Branch("mc1", &mcid1, "mcid1/I");
  tout->Branch("mc2", &mcid2, "mcid2/I");
  tout->Branch("pdg1", &pdgid1, "pdgid1/I");
  tout->Branch("pdg2", &pdgid2, "pdgid2/I");
  tout->Branch("ttCls", &ttCls, "ttCls/I");


  for (Long64_t i = 0; i < nentries; ++i){
    

    if(i%1000==0) cout << i << "/" << nentries << endl;

    t->GetEntry(i);
    
    if(njet<6) continue;
    
    ttCls = ttCls_;

    for( int j = 0; j < njet-1 ; ++j){
      
      event = event_;

      if(event%2==0){
      pt1    = pt[j];
      eta1   = TMath::Abs(eta[j]);
      mcid1  = abs(mcid[j]);
      pdgid1 = abs(pdgid[j]);
      csv1   = TMath::Max(TMath::Min(csv[j], float(1.0)), float(0.0));
      }
      else{
      pt2    = pt[j];
      eta2   = TMath::Abs(eta[j]);
      mcid2  = abs(mcid[j]);
      pdgid2 = abs(pdgid[j]);
      csv2   = TMath::Max(TMath::Min(csv[j], float(1.0)), float(0.0));
      }

      for( int k = j+1; k < njet ; ++k){	
      if(event%2==0){
	pt2    = pt[k];
	eta2   = TMath::Abs(eta[k]);
	mcid2  = abs(mcid[k]);
	pdgid2 = abs(pdgid[k]);       	
	csv2   = TMath::Max(TMath::Min(csv[k], float(1.0)), float(0.0));
      }
      else{
	pt1    = pt[k];
	eta1   = TMath::Abs(eta[k]);
	mcid1  = abs(mcid[k]);
	pdgid1 = abs(pdgid[k]);       	
	csv1   = TMath::Max(TMath::Min(csv[k], float(1.0)), float(0.0));
      }

	tout->Fill();
      }
    }
    
  }


  tout->Write();
  out->Close();

}
