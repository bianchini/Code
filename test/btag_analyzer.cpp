#include "interface/BTagRandomizer.h"
#include "TFile.h"
#include "TH3D.h"
#include "TTree.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include<iostream>

using namespace std;
using namespace MEM;

//TFile* fin = TFile::Open("/scratch/bianchi//tth_13tev_amcatnlo_pu20bx25.root");
//TFile* fin = TFile::Open("/scratch/bianchi///ttjets_13tev_madgraph_pu20bx25_phys14.root");

int main(int argc, char *argv[]){

  if(argc<3) return 0;

  double csv_cut{0.814};
  
  TFile* fout  = new TFile( Form("test/btag_out_%s_%d.root", argv[1], atoi(argv[2])),"RECREATE");
  TTree* tout  = new TTree("tree", "");
  int njet, ntag, ncat, event;
  float pcat        [99];
  int   pass        [99];
  int   pass_rnd    [99];
  float csv_rnd_0   [99];
  float csv_inp_0   [99];
  float csv_rnd_1   [99];
  float csv_inp_1   [99];
  float csv_rnd_2   [99];
  float csv_inp_2   [99];
  float csv_rnd_3   [99];
  float csv_inp_3   [99];
  float csv_rnd_4   [99];
  float csv_inp_4   [99];

  float corr_rnd_0   [8];
  float corr_inp_0   [8];
  float corr_rnd_1   [8];
  float corr_inp_1   [8];
  float corr_rnd_2   [8];
  float corr_inp_2   [8];
  float corr_rnd_3   [8];
  float corr_inp_3   [8];
  float corr_rnd_4   [8];
  float corr_inp_4   [8];

  float HT, met;
  float pt    [99];
  float eta   [99];
  int   pdgid [99];

  map<int, float*> csv_rnd_map;
  map<int, float*> csv_inp_map;
  csv_inp_map[0] = csv_inp_0;
  csv_rnd_map[0] = csv_rnd_0;
  csv_inp_map[1] = csv_inp_1;
  csv_rnd_map[1] = csv_rnd_1;
  csv_inp_map[2] = csv_inp_2;
  csv_rnd_map[2] = csv_rnd_2;
  csv_inp_map[3] = csv_inp_3;
  csv_rnd_map[3] = csv_rnd_3;
  csv_inp_map[4] = csv_inp_4;
  csv_rnd_map[4] = csv_rnd_4;

  map<int, float*> corr_rnd_map;
  map<int, float*> corr_inp_map;
  corr_inp_map[0] = corr_inp_0;
  corr_rnd_map[0] = corr_rnd_0;
  corr_inp_map[1] = corr_inp_1;
  corr_rnd_map[1] = corr_rnd_1;
  corr_inp_map[2] = corr_inp_2;
  corr_rnd_map[2] = corr_rnd_2;
  corr_inp_map[3] = corr_inp_3;
  corr_rnd_map[3] = corr_rnd_3;
  corr_inp_map[4] = corr_inp_4;
  corr_rnd_map[4] = corr_rnd_4;


  tout->Branch("event",    &event,    "event/I");
  tout->Branch("njet",     &njet,     "njet/I");
  tout->Branch("ncat",     &ncat,     "ncat/I");
  tout->Branch("ntag",     &ntag,     "ntag/I");
  tout->Branch("pcat",     pcat,      "pcat[ncat]/F");
  tout->Branch("pass",     pass,      "pass[ncat]/I");
  tout->Branch("pass_rnd", pass_rnd,  "pass_rnd[ncat]/I");
  tout->Branch("pt",       pt,        "pt[njet]/F");
  tout->Branch("eta",      eta,       "eta[njet]/F");
  tout->Branch("pdgid",    pdgid,     "pdgid[njet]/I");
  tout->Branch("HT",      &HT,      "HT/F");
  tout->Branch("met",     &met,     "met/F");

  for(auto it = csv_rnd_map.begin() ; it!=csv_rnd_map.end(); ++it)
    tout->Branch(Form("csv_rnd_%d",it->first),  it->second,   Form("csv_rnd_%d[njet]/F",it->first));

  for(auto it = csv_inp_map.begin() ; it!=csv_inp_map.end(); ++it)
    tout->Branch(Form("csv_inp_%d",it->first),  it->second,   Form("csv_inp_%d[njet]/F",it->first));

  for(auto it = corr_rnd_map.begin() ; it!=corr_rnd_map.end(); ++it)
    tout->Branch(Form("corr_rnd_%d",it->first),  it->second,   Form("corr_rnd_%d[8]/F",it->first));

  for(auto it = corr_inp_map.begin() ; it!=corr_inp_map.end(); ++it)
    tout->Branch(Form("corr_inp_%d",it->first),  it->second,   Form("corr_inp_%d[8]/F",it->first));


  TFile* f = TFile::Open("../MEAnalysis/root/ControlPlotsV6_finerPt.root","READ");
  TH3D* h3_b = (TH3D*)f->Get("csv_b_pt_eta");
  TH3D* h3_c = (TH3D*)f->Get("csv_c_pt_eta");
  TH3D* h3_l = (TH3D*)f->Get("csv_l_pt_eta");
  std::map<MEM::DistributionType::DistributionType, TH3D> btag_pdfs;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_b] = *h3_b ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_c] = *h3_c ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_l] = *h3_l ;

  BTagRandomizer* rnd = new BTagRandomizer(DebugVerbosity::init
					   //|DebugVerbosity::init_more
					   //|DebugVerbosity::event
					   , -1            // seed
   					   , btag_pdfs     // pdfs
					   , atoi(argv[2]) // assign random csv
					   , 1             // compress csv [0,1]
					   , 5000          // max number of toys
					   );

  
  TCut cut = "numJets<=6";
  TFile* fin = TFile::Open(Form("/scratch/bianchi/%s.root", argv[1]));
  TTree* tin = (TTree*)fin->Get("tree");

  TTreeFormula* treeformula = new TTreeFormula("selection", cut , tin );

  int njets;
  int is_sl;
  int evt;
  double nBCSVM;
  double jets_pt     [99];
  double jets_eta    [99];
  double jets_mass   [99];
  double jets_phi    [99];
  double jets_btagCSV[99];
  double jets_mcPt   [99];
  int jets_mcFlavour [99];
  double btag_LR_4b_2b;
  double met_pt;

  tin->SetBranchAddress("evt",         &evt);
  tin->SetBranchAddress("njets",       &njets);
  tin->SetBranchAddress("is_sl",       &is_sl);
  tin->SetBranchAddress("nBCSVM",      &nBCSVM);
  tin->SetBranchAddress("jets_pt",     jets_pt );
  tin->SetBranchAddress("jets_eta",    jets_eta  );
  tin->SetBranchAddress("jets_mass",   jets_mass  );
  tin->SetBranchAddress("jets_phi",    jets_phi  );
  tin->SetBranchAddress("jets_btagCSV",jets_btagCSV);
  tin->SetBranchAddress("jets_mcPt",   jets_mcPt   );
  tin->SetBranchAddress("jets_mcFlavour",  jets_mcFlavour   );
  tin->SetBranchAddress("btag_LR_4b_2b",   &btag_LR_4b_2b);
  tin->SetBranchAddress("met_pt",   &met_pt);

  Long64_t nentries = tin->GetEntries();
  cout << "Total entries: " << nentries << endl;
  int count_pass{0};

  for (Long64_t i = 0; i < nentries; i++){

    if( count_pass>100000 ) break;

    tin->GetEntry(i);

    tin->LoadTree(i);
    if( treeformula->EvalInstance() == 0){
      continue;
    }
    if(i%500==0) cout << "Event " << i << " (" << evt << ")" << endl;

    vector<JetCategory> cats;
    vector<Object*> objects_mem;
    vector<JetCategory*> cats_mem;

    HT  = 0.;
    met = met_pt; 
    for( int j = 0 ; j < njets ; ++j){

      pt[j]    = jets_pt[j];
      HT      += jets_pt[j];
      eta[j]   = jets_eta[j];
      pdgid[j] = jets_mcFlavour[j];

      TLorentzVector lv;
      lv.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_mass[j]);
      Object* jet = new Object( lv, ObjectType::Jet );
      jet->addObs( Observable::BTAG,  jets_btagCSV[j]>csv_cut   );  
      jet->addObs( Observable::PDGID, jets_mcFlavour[j]);  
      jet->addObs( Observable::CSV,   jets_btagCSV[j]  );  
      rnd->push_back_object( jet );

      if( j<5 ){
	int ntags_l = j; 
	int ntags_h = ntags_l<4 ? j : -1;
	JetCategory* cat = new JetCategory(ntags_l, ntags_h, csv_cut, j, 
					   string(Form("%dj%s%dt", njets, ntags_h>=0 ? "" : "ge" , ntags_l)));
	cats.push_back( *cat );
	cats_mem.push_back( cat );
      }

      objects_mem.push_back( jet );
    }

    if( njets==4 ){
      int ntags_l = 4; 
      int ntags_h = -1;
      JetCategory* cat = new JetCategory(ntags_l, ntags_h, csv_cut, 4, 
					 string(Form("%dj%s%dt", njets, ntags_h>=0 ? "" : "ge" , ntags_l)));
      cats.push_back( *cat );
      cats_mem.push_back( cat );
    }

    event = evt;
    njet  = njets;
    ncat  = int(cats.size());    
    ntag  = int(nBCSVM);
    

    //cout << "\tRunning BTagRandomizerOutput with " << cats.size() << " confs" << endl;
    vector<BTagRandomizerOutput> out_all = rnd->run_all( cats );
    for(int o = 0; o < int(out_all.size()) ; ++o ){
      pcat[o]     = out_all[o].p;
      pass[o]     = out_all[o].pass;
      pass_rnd[o] = out_all[o].pass_rnd;
      //o.print(cout);    

      for( int j = 0 ; j < njets ; ++j){
	(csv_rnd_map.at(o))[j] = (out_all[o].rnd_btag)[j];
	(csv_inp_map.at(o))[j] = (out_all[o].input_btag)[j];
      }

      //// Input
      float dR_min_bb =  99.;
      float dR_min_jj =  99.;
      float dR_max_bb = -99.;
      float dR_max_jj = -99.;
      float m_min_bb  =  99.;
      float m_min_jj  =  99.;
      float m_max_bb  = -99.;
      float m_max_jj  = -99.;
      for( int j = 0 ; j < njets-1 ; ++j){
	TLorentzVector lv_j;
	lv_j.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_mass[j]);
	int tag_j= (out_all[o].input_btag)[j]>csv_cut;
	for( int k = j+1 ; k < njets ; ++k){
	  TLorentzVector lv_k;
	  lv_k.SetPtEtaPhiM(jets_pt[k], jets_eta[k], jets_phi[k], jets_mass[k]);	
	  int tag_k = (out_all[o].input_btag)[k]>csv_cut;
	  
	  float dR = lv_j.DeltaR( lv_k );
	  float m  = (lv_k+lv_j).M();
	  
	  if( tag_j && tag_k ){
	    if( dR <= dR_min_bb ) dR_min_bb = dR;
	    if( m  <= m_min_bb  ) m_min_bb  = m;
	    if( dR >= dR_max_bb ) dR_max_bb = dR;
	    if( m  >= m_max_bb  ) m_max_bb  = m;
	  }
	  if( !tag_j && !tag_k ){
	    if( dR <= dR_min_jj ) dR_min_jj = dR;
	    if( m  <= m_min_jj  ) m_min_jj  = m;
	    if( dR >= dR_max_jj ) dR_max_jj = dR;
	    if( m  >= m_max_jj  ) m_max_jj  = m;
	  }
	  
	}
      }
      
      (corr_inp_map.at(o))[0] = dR_min_bb;
      (corr_inp_map.at(o))[1] = dR_min_jj;
      (corr_inp_map.at(o))[2] = dR_max_bb;
      (corr_inp_map.at(o))[3] = dR_max_jj;
      (corr_inp_map.at(o))[4] = m_min_bb;
      (corr_inp_map.at(o))[5] = m_min_jj;
      (corr_inp_map.at(o))[6] = m_max_bb;
      (corr_inp_map.at(o))[7] = m_max_jj;

      ///////// Rand
      dR_min_bb =  99.;
      dR_min_jj =  99.;
      dR_max_bb = -99.;
      dR_max_jj = -99.;
      m_min_bb  =  99.;
      m_min_jj  =  99.;
      m_max_bb  = -99.;
      m_max_jj  = -99.;
      for( int j = 0 ; j < njets-1 ; ++j){
	TLorentzVector lv_j;
	lv_j.SetPtEtaPhiM(jets_pt[j], jets_eta[j], jets_phi[j], jets_mass[j]);
	int tag_j = (out_all[o].rnd_btag)[j]>csv_cut;
	for( int k = j+1 ; k < njets ; ++k){
	  TLorentzVector lv_k;
	  lv_k.SetPtEtaPhiM(jets_pt[k], jets_eta[k], jets_phi[k], jets_mass[k]);	
	  int tag_k = (out_all[o].rnd_btag)[k]>csv_cut;
	  
	  float dR = lv_j.DeltaR( lv_k );
	  float m  = (lv_k+lv_j).M();
	  
	  if( tag_j && tag_k ){
	    if( dR <= dR_min_bb ) dR_min_bb = dR;
	    if( m  <= m_min_bb  ) m_min_bb  = m;
	    if( dR >= dR_max_bb ) dR_max_bb = dR;
	    if( m  >= m_max_bb  ) m_max_bb  = m;
	  }
	  if( !tag_j && !tag_k ){
	    if( dR <= dR_min_jj ) dR_min_jj = dR;
	    if( m  <= m_min_jj  ) m_min_jj  = m;
	    if( dR >= dR_max_jj ) dR_max_jj = dR;
	    if( m  >= m_max_jj  ) m_max_jj  = m;
	  }
	  
	}
      }
      
      (corr_rnd_map.at(o))[0] = dR_min_bb;
      (corr_rnd_map.at(o))[1] = dR_min_jj;
      (corr_rnd_map.at(o))[2] = dR_max_bb;
      (corr_rnd_map.at(o))[3] = dR_max_jj;
      (corr_rnd_map.at(o))[4] = m_min_bb;
      (corr_rnd_map.at(o))[5] = m_min_jj;
      (corr_rnd_map.at(o))[6] = m_max_bb;
      (corr_rnd_map.at(o))[7] = m_max_jj;
    }


    rnd->next_event();
    for( auto c : cats_mem    ) delete c;
    for( auto o : objects_mem ) delete o;

    ++count_pass;
    tout->Fill();
  }
    

  cout << "Done!" << endl;
  delete rnd;
  f->Close();
  fin->Close();
 
  fout->cd();
  tout->Write();
  fout->Close();
}
