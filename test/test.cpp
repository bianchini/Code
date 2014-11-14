#include<iostream>
#include<memory> 

#include "interface/HypoTester.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"


using namespace std;

int main(){

  TTree* t    = new TTree("tree", "tree");
  Algo::HypoTester* tester = new Algo::HypoTester(t) ;
  
  // test
  TLorentzVector j1;
  j1.SetPtEtaPhiM( 50., 0., TMath::Pi()*2/3, 0.);
  tester->push_back_object( j1  , 'j');

  TLorentzVector j2;
  j2.SetPtEtaPhiM( 30., 0.,  0., 0.);
  tester->push_back_object( j2  , 'j');

  TLorentzVector j3;
  j3.SetPtEtaPhiM( 50., 0., -TMath::Pi()/2, 0.);
  tester->push_back_object( j3  , 'j');

  TLorentzVector j4;
  j4.SetPtEtaPhiM(100., 0., -TMath::Pi()/4, 0.);
  //tester->push_back_object( j4  , 'j');

  TLorentzVector j5;
  j5.SetPtEtaPhiM(80., 0., -TMath::Pi()/5, 0.);
  //tester->push_back_object( j5  , 'j');

  TLorentzVector j6;
  j6.SetPtEtaPhiM(120., 0., -TMath::Pi()/6, 0.);
  //tester->push_back_object( j6  , 'j');

  TLorentzVector lep;
  lep.SetPtEtaPhiM( 50., 1., -TMath::Pi()/3, 0.);
  //tester->push_back_object( lep  , 'l');

  TLorentzVector met;
  met.SetPtEtaPhiM( 0., 0., 0., 0.);
  //tester->push_back_object( met  , 'm');

  map<string, vector<Algo::Decay> > hypotheses;
  hypotheses["H0"] = {Algo::Decay::TopHad};
  hypotheses["H1"] = {Algo::Decay::Radiation, Algo::Decay::Radiation, Algo::Decay::Radiation};
  hypotheses["H2"] = {Algo::Decay::WHad, Algo::Decay::Radiation};
  tester->test( hypotheses );

  TFile* fout = new TFile("test/Test.root", "RECREATE");
  fout->cd();
  t->Write("", TObject::kOverwrite);
  fout->Close();
  cout << "ROOT file saved" << endl;

  delete tester;

  /*
  // assumptions
  //tester->assume( Algo::Decay::TopHad  );
  //tester->assume( Algo::Decay::WHad  );
  tester->assume( Algo::Decay::TopLep   );
  //tester->assume( Algo::Decay::Higgs );

  // printout
  tester->print(cout);

  // run the algo
  tester->init();
 
  // run the algo
  tester->run();
  */

}
