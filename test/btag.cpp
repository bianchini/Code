#include "interface/BTagRandomizer.h"
#include "TFile.h"
#include "TH3D.h"

#include<iostream>

using namespace std;
using namespace MEM;

int main(int argc, char *argv[]){

  TFile* fout  = new TFile("out.root","RECREATE");
  vector<TH1F*> h_inp;
  vector<TH1F*> h_ran;
  vector<TH1F*> h_ratio;
  for(int j = 0 ; j < 6 ; ++j){
    TH1F* inp = new TH1F(Form("h%d_inp",j), "", 20,0,1); 
    TH1F* ran = new TH1F(Form("h%d_ran",j), "", 20,0,1); 
    TH1F* ratio = new TH1F(Form("h%d_ratio",j), "", 20,0,1); 
    inp->Sumw2();
    ran->Sumw2();
    ratio->Sumw2();
    h_inp.push_back( inp );
    h_ran.push_back( ran );
    h_ratio.push_back( ratio );
  }

  //TFile* f = TFile::Open("../MEAnalysis/root/ControlPlotsV6.root","READ");
  TFile* f = TFile::Open("/shome/bianchi/TTH-74X-heppy//CMSSW/src/TTH/MEAnalysis/root/csv.root");
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
					   , -1
					   , btag_pdfs 
					   , 1
					   , 5000
					   );

  TLorentzVector lv_j1;
  lv_j1.SetPtEtaPhiM(242.816604614, -0.107542805374, 1.25506973267, 24.5408706665);
  Object j1( lv_j1, ObjectType::Jet );
  j1.addObs( Observable::BTAG,  0.   );  
  j1.addObs( Observable::PDGID, 5    );  
  j1.addObs( Observable::CSV,   0.9  );  
  j1.addObs( Observable::IGNORE_FOR_RND,  1 );  
  //j1.addObs( Observable::BTAGPROB, 0.02 );  

  TLorentzVector lv_j2;
  lv_j2.SetPtEtaPhiM(35.6511192322, 0.566395223141, -2.51394343376, 8.94268417358);
  Object j2( lv_j2, ObjectType::Jet );
  j2.addObs( Observable::BTAG, 1. );
  j2.addObs( Observable::PDGID, 3 );  
  j2.addObs( Observable::CSV,   0.2  );  
  j2.addObs( Observable::IGNORE_FOR_RND,  1 );  
  //j2.addObs( Observable::BTAGPROB, 0.70 );  

  TLorentzVector lv_j3;
  lv_j3.SetPtEtaPhiM(77.6708831787, -0.709680855274, -2.53739523888, 10.4904966354);
  Object j3( lv_j3, ObjectType::Jet );
  j3.addObs( Observable::BTAG, 0. );
  j3.addObs( Observable::PDGID, 4 );
  j3.addObs( Observable::CSV,   0.4 );  
  j3.addObs( Observable::IGNORE_FOR_RND,  1 );  
  //j3.addObs( Observable::BTAGPROB, 0.02 );    

  TLorentzVector lv_j4;
  lv_j4.SetPtEtaPhiM(52.0134391785, -0.617823541164, -1.23360788822, 6.45914268494);
  Object j4( lv_j4, ObjectType::Jet );
  j4.addObs( Observable::BTAG, 0. );
  j4.addObs( Observable::PDGID, 21 );  
  j4.addObs( Observable::CSV, 0.12 );  
  j4.addObs( Observable::IGNORE_FOR_RND,  1 );  
  //j4.addObs( Observable::BTAGPROB, 0.70 );    

  TLorentzVector lv_j5;
  lv_j5.SetPtEtaPhiM( 235.892044067, -0.997860729694, -2.10646605492, 27.9887943268 );
  Object j5( lv_j5, ObjectType::Jet );
  j5.addObs( Observable::BTAG, 1. );
  j5.addObs( Observable::PDGID, 5 );
  j5.addObs( Observable::CSV,  0.7 );  
  j5.addObs( Observable::IGNORE_FOR_RND,  1 );  
  //j5.addObs( Observable::BTAGPROB, 0.70 );       

  TLorentzVector lv_j6;
  lv_j6.SetPtEtaPhiM(191.423553467, -0.46368226409, 0.750520706177, 30.5682048798);
  Object j6( lv_j6, ObjectType::Jet );
  j6.addObs( Observable::BTAG, 1. );
  j6.addObs( Observable::PDGID, -5 );
  j6.addObs( Observable::CSV,  0.89 );  
  j6.addObs( Observable::IGNORE_FOR_RND,  1 );  
  //j6.addObs( Observable::BTAGPROB, 0.70 );       


  int    nmax = argc>1 ? std::atoi(argv[1]) : 1;
  int    n_pass{0};
  int    n_pass_rnd{0};
  double n_pass_weight{0.};

  double n_pass_weight2{0.};

  /*
  if(argc>4)
    rnd->set_condition(std::atoi(argv[2]),std::atoi(argv[3]), std::atof(argv[4]));
  else
    rnd->set_condition(3, -1,  0.814);
  */

  for(int i = 0 ; i < nmax ; ++i){


    if(i%1000==0) cout << "Toy " << i << endl;

    rnd->push_back_object( &j1 );
    rnd->push_back_object( &j2 );
    rnd->push_back_object( &j3 );
    rnd->push_back_object( &j4 );
    rnd->push_back_object( &j5 );
    rnd->push_back_object( &j6 );
    
    //JetCategory cat0 = JetCategory(0, 0, 0.879, 0, "6j0t");
    //JetCategory cat1 = JetCategory(1, 1, 0.814, 1, "6j1t");
    //JetCategory cat2 = JetCategory(2, 2, 0.879, 2, "6j2t");
    //JetCategory cat3 = JetCategory(3, 3, 0.879, 3, "6j3t");
    JetCategory cat4 = JetCategory(4, 6, 0.879, 4, "6jge4t");
    vector<BTagRandomizerOutput> out_all = rnd->run_all(vector<JetCategory>{/*cat0, cat1, cat2, cat3, cat4*/ cat4});
    //for(auto o : out_all ) o.print(cout);

    //BTagRandomizerOutput out = rnd->run();
    BTagRandomizerOutput out = out_all[0];
    if(nmax<10) out.print(cout);

    if(out.pass){ 
      ++n_pass;
      for(int j = 0 ; j < 1 ; ++j) h_inp[j]->Fill((out.input_btag)[j]);
    }
    if(out.pass_rnd){
      ++n_pass_rnd;
      n_pass_weight  += out.p;
      n_pass_weight2 += out.p*out.p; 
      for(int j = 0 ; j < 1 ; ++j) h_ran[j]->Fill((out.rnd_btag)[j], out.p);
    }
    
    rnd->next_event();
  }




  for(int j = 0 ; j < 1 ; ++j) h_ratio[j]->Divide(h_inp[j], h_ran[j], 1./h_inp[j]->Integral(), 1./h_ran[j]->Integral());

  double pass_prob = double(n_pass)/nmax;
  cout << "N toys        = " << nmax << endl;
  cout << "Passing       = " << n_pass << " +/- " << TMath::Sqrt(pass_prob*(1-pass_prob)*nmax) << endl;
  cout << "Passing (rnd) = " << n_pass_rnd << endl;
  cout << "Weighted      = " << n_pass_weight << " +/- " << TMath::Sqrt(n_pass_weight2) << endl;

  cout << "Done!" << endl;
  delete rnd;
  f->Close();
  
  fout->cd();
  for(int j = 0 ; j < 1 ; ++j) h_ratio[j]->Write();
  for(int j = 0 ; j < 1 ; ++j) h_inp[j]->Write();
  for(int j = 0 ; j < 1 ; ++j) h_ran[j]->Write();
  fout->Close();
}
