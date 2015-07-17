

{

  TFile* f = TFile::Open("/scratch/bianchi/new/ttjets_13tev_madgraph_pu20bx25_phys14.root");
  TTree* t = (TTree*)f->Get("tree");
  int numJets;
  double jets_pt   [15];
  double jets_mcPt [15];
  double jets_eta  [15];
  double jets_mcEta[15];
  double jets_btagCSV [15];
  int   jets_mcFlavour[15];
  int   jets_mcMatchId[15];

  t->SetBranchAddress("numJets",         &numJets);
  t->SetBranchAddress("jets_pt",         jets_pt);
  t->SetBranchAddress("jets_mcPt",       jets_mcPt);
  t->SetBranchAddress("jets_eta",        jets_eta);
  t->SetBranchAddress("jets_mcEta",      jets_mcEta);
  t->SetBranchAddress("jets_btagCSV",    jets_btagCSV);
  t->SetBranchAddress("jets_mcFlavour",  jets_mcFlavour);
  t->SetBranchAddress("jets_mcMatchId",  jets_mcMatchId);

  Long64_t nentries = t->GetEntries();
  cout << "Total entries: " << nentries << endl;

  TFile* out  = TFile::Open("csv_byFlavour.root","RECREATE");
  out->cd();

  vector<TString> fl = {"b","c","l", "b_t", "b_g", "c_t", "c_g", "s", "u", "g"};
  map<TString, TH3D*> hmap_rec;
  map<TString, TH3D*> hmap_gen;

  map<TString, TH2D*> heff_rec;

  const int nBinsX = 21;
  double binsX[nBinsX+1] = {20., 25., 30., 35., 40., 45., 50., 60., 70., 80., 90., 100., 120., 140., 160., 180., 200., 250., 300., 350., 400., 500.}; 
  
  const int nBinsY = 8;
  double binsY[nBinsY+1] = { 0.0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.50 };

  const int nBinsZ = 100;
  double binsZ[nBinsZ+1];
  for(int k = 0 ; k <= nBinsZ ; ++k) binsZ[k] = k/100.;

  for( int f = 0 ; f < fl.size() ; ++f){
    hmap_rec[fl[f]]  = new TH3D("csv_"+fl[f]+"_pt_eta",     "csv_"+fl[f]+"_pt_eta",     nBinsX, binsX, nBinsY, binsY, nBinsZ, binsZ);
    heff_rec[fl[f]]  = new TH2D("eff_"+fl[f]+"_pt_eta",     "eff_"+fl[f]+"_pt_eta",     nBinsX, binsX, nBinsY, binsY);
    hmap_gen[fl[f]]  = new TH3D("csv_"+fl[f]+"_mcpt_mceta", "csv_"+fl[f]+"_mcpt_mceta", nBinsX, binsX, nBinsY, binsY, nBinsZ, binsZ);
  }

  for (Long64_t i = 0; i < nentries; ++i){
    
    if(i%5000==0) cout << i << "/" << nentries << " [ " << float(i)/nentries*100 << "% ]" << endl;

    t->GetEntry(i);
    
    for( int j = 0; j < numJets ; ++j){
      TString fl = "l";
      int pdg = std::abs( jets_mcFlavour[j] );
      if     (pdg==5) fl = "b";
      else if(pdg==4) fl = "c";
      hmap_rec[fl]->Fill( jets_pt[j],   jets_eta[j],   jets_btagCSV[j]);
      hmap_gen[fl]->Fill( jets_mcPt[j], jets_mcEta[j], jets_btagCSV[j]);
      if( fl=="b" ){
	if( std::abs(jets_mcMatchId[j])==6 ) fl += "_t";
	else fl += "_g";
	hmap_rec[fl]->Fill( jets_pt[j],   jets_eta[j],   jets_btagCSV[j]);
	hmap_gen[fl]->Fill( jets_mcPt[j], jets_mcEta[j], jets_btagCSV[j]);	
      }
      if( fl=="c" ){
	if( std::abs(jets_mcMatchId[j])==0 ) fl += "_g";
	else fl += "_t";
	hmap_rec[fl]->Fill( jets_pt[j],   jets_eta[j],   jets_btagCSV[j]);
	hmap_gen[fl]->Fill( jets_mcPt[j], jets_mcEta[j], jets_btagCSV[j]);	
      }
      if( pdg<4 || pdg==21){
	if (pdg==1 || pdg==2){
	  fl = "u";
	  hmap_rec[fl]->Fill( jets_pt[j],   jets_eta[j],   jets_btagCSV[j]);
	  hmap_gen[fl]->Fill( jets_mcPt[j], jets_mcEta[j], jets_btagCSV[j]);
	}
	if (pdg==3){
	  fl = "s";
	  hmap_rec[fl]->Fill( jets_pt[j],   jets_eta[j],   jets_btagCSV[j]);
	  hmap_gen[fl]->Fill( jets_mcPt[j], jets_mcEta[j], jets_btagCSV[j]);
	}
	if (pdg==21){
	  fl = "g";
	  hmap_rec[fl]->Fill( jets_pt[j],   jets_eta[j],   jets_btagCSV[j]);
	  hmap_gen[fl]->Fill( jets_mcPt[j], jets_mcEta[j], jets_btagCSV[j]);
	}
      }
    }
  }
 
  out->cd();
  for( int f = 0 ; f < fl.size() ; ++f){

    TH3D* h = hmap_rec[fl[f]];
    for( int x = 0 ; x <=h->GetNbinsX()+1; ++x ){
      for( int y = 0 ; y <=h->GetNbinsY()+1; ++y ){
	float tot  = 0.;
	float pass = 0.;
	for( int z = 0 ; z <=h->GetNbinsZ()+1; ++z ){
	  float add = h->GetBinContent(x,y,z);
	  tot += add;
	  if( z >= h->GetZaxis()->FindBin(0.814) ) pass += add;
	}
	heff_rec[fl[f]]->SetBinContent(x,y, tot>0. ? pass/tot : 0.);      
	heff_rec[fl[f]]->SetBinError  (x,y, tot>0. ? sqrt(pass/tot*(1-pass/tot)/tot) : 0.);      
      }
    }
    

    heff_rec[fl[f]]->Write();   
    hmap_rec[fl[f]]->Write();   
    hmap_gen[fl[f]]->Write();   
  }
  out->Close();

}
