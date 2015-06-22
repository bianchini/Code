#include "interface/BTagRandomizer.h"
#include "TFile.h"
#include "TH3D.h"

#include<iostream>

using namespace std;
using namespace MEM;

int main(int argc, char *argv[]){

  TFile* fout  = new TFile("out.root","RECREATE");
  TH1F* h0_ran   = new TH1F("h0_ran", "", 100,0,1); 
  TH1F* h0_in    = new TH1F("h0_in",  "", 100,0,1); 

  TFile* f = TFile::Open("../MEAnalysis/root/ControlPlotsV6.root","READ");
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
					   );

  TLorentzVector lv_j1;
  lv_j1.SetPtEtaPhiM(242.816604614, -0.107542805374, 1.25506973267, 24.5408706665);
  Object j1( lv_j1, ObjectType::Jet );
  j1.addObs( Observable::BTAG,  0.   );  
  j1.addObs( Observable::PDGID, 1    );  
  j1.addObs( Observable::CSV,   0.1  );  
  //j1.addObs( Observable::BTAGPROB, 0.02 );  

  TLorentzVector lv_j2;
  lv_j2.SetPtEtaPhiM(35.6511192322, 0.566395223141, -2.51394343376, 8.94268417358);
  Object j2( lv_j2, ObjectType::Jet );
  j2.addObs( Observable::BTAG, 1. );
  j2.addObs( Observable::PDGID, 1 );  
  j2.addObs( Observable::CSV,   0.9  );  
  //j2.addObs( Observable::BTAGPROB, 0.70 );  

  TLorentzVector lv_j3;
  lv_j3.SetPtEtaPhiM(77.6708831787, -0.709680855274, -2.53739523888, 10.4904966354);
  Object j3( lv_j3, ObjectType::Jet );
  j3.addObs( Observable::BTAG, 0. );
  j3.addObs( Observable::PDGID, 1 );
  j3.addObs( Observable::CSV,   0.2  );  
  //j3.addObs( Observable::BTAGPROB, 0.02 );    

  TLorentzVector lv_j4;
  lv_j4.SetPtEtaPhiM(52.0134391785, -0.617823541164, -1.23360788822, 6.45914268494);
  Object j4( lv_j4, ObjectType::Jet );
  j4.addObs( Observable::BTAG, 1. );
  j4.addObs( Observable::PDGID, 1 );  
  j4.addObs( Observable::CSV, 0.7 );  
  //j4.addObs( Observable::BTAGPROB, 0.70 );    

  TLorentzVector lv_j5;
  lv_j5.SetPtEtaPhiM( 235.892044067, -0.997860729694, -2.10646605492, 27.9887943268 );
  Object j5( lv_j5, ObjectType::Jet );
  j5.addObs( Observable::BTAG, 1. );
  j5.addObs( Observable::PDGID, 5 );
  j5.addObs( Observable::CSV,  0.95 );  
  //j5.addObs( Observable::BTAGPROB, 0.70 );       

  TLorentzVector lv_j6;
  lv_j6.SetPtEtaPhiM(191.423553467, -0.46368226409, 0.750520706177, 30.5682048798);
  Object j6( lv_j6, ObjectType::Jet );
  j6.addObs( Observable::BTAG, 1. );
  j6.addObs( Observable::PDGID, -5 );
  j6.addObs( Observable::CSV,  0.89 );  
  //j6.addObs( Observable::BTAGPROB, 0.70 );       


  int    nmax = argc>1 ? std::atoi(argv[1]) : 1;
  int    n_pass{0};
  int    n_pass_rnd{0};
  double n_pass_weight{0.};
  double n_pass_weight2{0.};

  if(argc>4)
    rnd->set_condition(std::atoi(argv[2]),std::atoi(argv[3]), std::atof(argv[4]));
  else
    rnd->set_condition(4, -1,  0.879);

  for(int i = 0 ; i < nmax ; ++i){

    rnd->push_back_object( &j1 );
    rnd->push_back_object( &j2 );
    rnd->push_back_object( &j3 );
    rnd->push_back_object( &j4 );
    rnd->push_back_object( &j5 );
    rnd->push_back_object( &j6 );
    
    BTagRandomizerOutput out = rnd->run();
    if(nmax<10) out.print(cout);

    if(out.pass){ 
      ++n_pass;
      h0_in ->Fill((out.input_btag)[0]);
    }
    if(out.pass_rnd || (!out.pass_rnd && out.pass)){
      ++n_pass_rnd;
      n_pass_weight  += out.p;
      n_pass_weight2 += out.p*out.p; 
      h0_ran->Fill((out.rnd_btag)[0], out.p);   
    }

    rnd->next_event();
  }

  double pass_prob = double(n_pass)/nmax;
  cout << "N toys        = " << nmax << endl;
  cout << "Passing       = " << n_pass << " +/- " << TMath::Sqrt(pass_prob*(1-pass_prob)*nmax) << endl;
  cout << "Passing (rnd) = " << n_pass_rnd << endl;
  cout << "Weighted      = " << n_pass_weight << " +/- " << TMath::Sqrt(n_pass_weight2) << endl;

  cout << "Done!" << endl;
  delete rnd;
  f->Close();
  
  fout->cd();
  h0_ran->Write();
  h0_in->Write();
  fout->Close();
}
