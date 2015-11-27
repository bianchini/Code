

{

  Long64_t nmax = 5e+06;

  TString path = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV14/";
  vector<TString> samples;
  
  //samples.push_back("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV14/ZJetsToNuNu_HT-100To200_13TeV-madgraph/VHBB_HEPPY_V14_ZJetsToNuNu_HT-100To200_13TeV-madgraph__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_083120/0000/");
  samples.push_back("TT_TuneCUETP8M1_13TeV-powheg-pythia8/VHBB_HEPPY_V14_TT_TuneCUETP8M1_13TeV-powheg-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_212301/0000/");
  //samples.push_back("WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V14_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_220138/0000/");
  //samples.push_back("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV14/ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8/VHBB_HEPPY_V14_ZH_HToBB_ZToNuNu_M125_13TeV_amcatnloFXFX_madspin_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_161706/0000/");

  TChain* ch = new TChain("tree");
  for(unsigned int ns = 0; ns<samples.size(); ++ns){
    TString sample = samples[ns];
    cout << sample << endl;
    for(int nf=1; nf < 50; ++nf)
      ch->AddFile(path+sample+"/tree_"+TString(Form("%d",nf))+".root");
  }

  int nJet;
  int nGenJet;
  float Jet_pt   [15];
  float Jet_mcPt [15];
  float Jet_eta  [15];
  float Jet_phi  [15];
  float Jet_mcEta[15];
  float Jet_btagCSV [15];
  int Jet_hadronFlavour[15];
  int Jet_mcFlavour[15];
  int Jet_mcIdx[15];
  int GenJet_numBHadrons[15];
  int GenJet_numBHadronsFromTop[15];
  int GenJet_numBHadronsAfterTop[15];
  int GenJet_numCHadrons[15];
  int GenJet_numCHadronsFromTop[15];
  int GenJet_numCHadronsAfterTop[15];

  ch->SetBranchAddress("nJet",           &nJet);
  ch->SetBranchAddress("nGenJet",        &nGenJet);
  ch->SetBranchAddress("Jet_pt",         Jet_pt);
  ch->SetBranchAddress("Jet_mcPt",       Jet_mcPt);
  ch->SetBranchAddress("Jet_eta",        Jet_eta);
  ch->SetBranchAddress("Jet_phi",        Jet_phi);
  ch->SetBranchAddress("Jet_mcEta",      Jet_mcEta);
  ch->SetBranchAddress("Jet_btagCSV",    Jet_btagCSV);
  ch->SetBranchAddress("Jet_hadronFlavour",  Jet_hadronFlavour);
  ch->SetBranchAddress("Jet_mcFlavour",      Jet_mcFlavour);
  ch->SetBranchAddress("Jet_mcIdx",          Jet_mcIdx);
  ch->SetBranchAddress("GenJet_numBHadrons",              GenJet_numBHadrons);
  ch->SetBranchAddress("GenJet_numBHadronsAfterTop",      GenJet_numBHadronsAfterTop);
  ch->SetBranchAddress("GenJet_numBHadronsFromTop",       GenJet_numBHadronsFromTop);
  ch->SetBranchAddress("GenJet_numCHadrons",              GenJet_numCHadrons);
  ch->SetBranchAddress("GenJet_numCHadronsAfterTop",      GenJet_numCHadronsAfterTop);
  ch->SetBranchAddress("GenJet_numCHadronsFromTop",       GenJet_numCHadronsFromTop);


  Long64_t nentries = ch->GetEntries();
  cout << "Total entries: " << nentries << endl;

  TFile* out  = TFile::Open("csv.root","RECREATE");
  out->cd();

  vector<TString> fl = {"b", "b_dR_0p0_1p0",  "b_dR_1p0_2p0",  "b_dR_2p0_Inf",
			"1b","1b_dR_0p0_1p0", "1b_dR_1p0_2p0", "1b_dR_2p0_Inf",
			"2b","2b_dR_0p0_1p0", "2b_dR_1p0_2p0", "2b_dR_2p0_Inf",
			"c", "c_dR_0p0_1p0",  "c_dR_1p0_2p0",  "c_dR_2p0_Inf",
			"1c","1c_dR_0p0_1p0", "1c_dR_1p0_2p0", "1c_dR_2p0_Inf",
			"2c","2c_dR_0p0_1p0", "2c_dR_1p0_2p0", "2c_dR_2p0_Inf",
			"l", "l_dR_0p0_1p0",  "l_dR_1p0_2p0",  "l_dR_2p0_Inf",
			"s", "s_dR_0p0_1p0",  "s_dR_1p0_2p0",  "s_dR_2p0_Inf",
			"ud","ud_dR_0p0_1p0", "ud_dR_1p0_2p0", "ud_dR_2p0_Inf",
			"g", "g_dR_0p0_1p0",  "g_dR_1p0_2p0",  "g_dR_2p0_Inf",
			"pu", "pu_dR_0p0_1p0",  "pu_dR_1p0_2p0",  "pu_dR_2p0_Inf"
  };

  map<TString, TH3D*> hmap_rec;
  map<TString, TH3D*> hmap_gen;

  map<TString, TH2D*> heff_rec;

  const int nBinsX = 12;
  double binsX[nBinsX+1] = {20., 25., 30., 35., 40., 45., 50., 60., 80., 100., 200., 300., 500.}; 
  
  const int nBinsY = 4;
  double binsY[nBinsY+1] = { 0.0, 1.0, 1.5, 2.0, 2.50 };

  const int nBinsZ = 102;
  double binsZ[nBinsZ+1];
  for(int k = 0 ; k <= nBinsZ ; ++k) binsZ[k] = -1./100 + k/100.;

  for( int f = 0 ; f < fl.size() ; ++f){
    hmap_rec[fl[f]]  = new TH3D("csv_"+fl[f]+"_pt_eta",     "csv_"+fl[f]+"_pt_eta",     nBinsX, binsX, nBinsY, binsY, nBinsZ, binsZ);
    heff_rec[fl[f]]  = new TH2D("eff_"+fl[f]+"_pt_eta",     "eff_"+fl[f]+"_pt_eta",     nBinsX, binsX, nBinsY, binsY);
    hmap_gen[fl[f]]  = new TH3D("csv_"+fl[f]+"_mcpt_mceta", "csv_"+fl[f]+"_mcpt_mceta", nBinsX, binsX, nBinsY, binsY, nBinsZ, binsZ);
  }

  for (Long64_t i = 0; i < std::min(nentries,nmax); ++i){
    
    if(i%5000==0) cout << i << "/" << nentries << " [ " << float(i)/nentries*100 << "% ]" << endl;

    ch->GetEntry(i);
    
    for( int j = 0; j < nJet ; ++j){

      float min_dR = 999.;
      for( int jj = 0; jj < nJet ; ++jj){
	if( Jet_pt[jj]<20. || std::abs(Jet_eta[jj])>2.5 ) continue;
	if(jj==j) continue;
	float dR = TMath::Sqrt( TMath::Power(Jet_eta[j]-Jet_eta[jj], 2) +  TMath::Power(Jet_phi[j]-Jet_phi[jj], 2)  );
	if( dR < min_dR ) min_dR = dR;
      }
      TString dR_bin = "dR_0p0_1p0"; 
      if     ( min_dR < 1.0 ) dR_bin = "dR_0p0_1p0";
      else if( min_dR < 2.0 ) dR_bin = "dR_1p0_2p0";
      else                    dR_bin = "dR_2p0_Inf";

      TString fl = "";
      int hadflav = std::abs( Jet_hadronFlavour[j] );
      int mcflav  = std::abs( Jet_mcFlavour[j] );
      
      float csv = Jet_btagCSV[j];
      if(csv<0.) csv = -0.005;
      if(csv>1.) csv = +1.005;

      int mcIdx = Jet_mcIdx[j]>=0 && Jet_mcPt[j]>20. && std::abs(Jet_mcEta[j])<2.4? Jet_mcIdx[j] : -1;

      if( Jet_pt[j]<20. || std::abs(Jet_eta[j])>2.5 ) continue;

      // PU jets
      if( Jet_mcPt[j]<5. ) {
	fl = "pu";
	hmap_rec[fl]->Fill( Jet_pt[j], Jet_eta[j],   csv);
	hmap_gen[fl]->Fill( Jet_pt[j], Jet_mcEta[j], csv);
	hmap_rec[fl+"_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j], csv);
	hmap_gen[fl+"_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j], csv);	
	continue;
      }

      switch(hadflav){
      case 5:
	fl = "b";
	if( mcIdx>=0 && GenJet_numBHadrons[mcIdx] == 1 )
	  fl = "1b";
	else if( mcIdx>=0 && GenJet_numBHadrons[mcIdx] >= 2 )
	  fl = "2b";
	hmap_rec["b"]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	hmap_gen["b"]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);	
	hmap_rec["b_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	hmap_gen["b_"+dR_bin]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);		
	if(fl!="b"){
	  hmap_rec[fl]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	  hmap_gen[fl]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);
	  hmap_rec[fl+"_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	  hmap_gen[fl+"_"+dR_bin]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);
	}
	break;
      case 4:
	fl = "c";
	if( mcIdx>=0 && GenJet_numCHadrons[mcIdx]==1 )
	  fl = "1c";
	if( mcIdx>=0 && GenJet_numCHadrons[mcIdx]>=2 )
	  fl = "2c";
	hmap_rec["c"]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	hmap_gen["c"]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);	
	hmap_rec["c_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	hmap_gen["c_"+dR_bin]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);		
	if(fl!="c"){
	  hmap_rec[fl]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	  hmap_gen[fl]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);
	  hmap_rec[fl+"_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	  hmap_gen[fl+"_"+dR_bin]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);
	}
	break;
      case 0:
	fl = "l";
	if     ( mcflav==3 )              fl = "s";
	else if( mcflav==2 || mcflav==1 ) fl = "ud";
	else if( mcflav==21 )             fl = "g";
	hmap_rec["l"]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	hmap_gen["l"]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);	
	hmap_rec["l_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	hmap_gen["l_"+dR_bin]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);
	if(fl!="l"){
	  hmap_rec[fl]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	  hmap_gen[fl]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);
	  hmap_rec[fl+"_"+dR_bin]->Fill( Jet_pt[j],   Jet_eta[j],   csv);
	  hmap_gen[fl+"_"+dR_bin]->Fill( Jet_mcPt[j], Jet_mcEta[j], csv);
	}
	break;
      default:
	cout << "Not a valid hadronFlavour" << endl;
	break;
      }

      //cout << fl << ", " << Jet_pt[j] << ", " << Jet_eta[j] << ", " << csv << endl;       

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
	  if( z >= h->GetZaxis()->FindBin(0.890) ) pass += add;
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
