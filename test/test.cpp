#include<iostream>
#include<memory> 

#include "interface/HypoTester.h"

#include "TMath.h"
#include "TLorentzVector.h"


using namespace std;

int main(){

  unique_ptr<Algo::HypoTester> tester{ new Algo::HypoTester() };
  
  // test
  TLorentzVector j1;
  j1.SetPtEtaPhiM( 60., 1.6, 0., 0.);
  tester->push_back_object( j1  , 'j');

  TLorentzVector j2;
  j2.SetPtEtaPhiM( 40., 0., -TMath::Pi(), 0.);
  tester->push_back_object( j2  , 'j');

  TLorentzVector j3;
  j3.SetPtEtaPhiM( 50., 0., +TMath::Pi()/2, 0.);
  tester->push_back_object( j3  , 'j');

  TLorentzVector j4;
  j4.SetPtEtaPhiM(100., 0., -TMath::Pi()/4, 0.);
  tester->push_back_object( j4  , 'j');

  TLorentzVector lep;
  lep.SetPtEtaPhiM( 50., 1., +TMath::Pi()/3, 0.);
  tester->push_back_object( lep  , 'l');

  TLorentzVector met;
  met.SetPtEtaPhiM( 30., 0., 0., 0.);
  tester->push_back_object( met  , 'm');

  // assumptions
  tester->assume( Algo::Decay::TopHad  );
  //tester->assume( Algo::Decay::WHad  );
  tester->assume( Algo::Decay::TopLep   );
  //tester->assume( Algo::Decay::HiggsHad );

  // printout
  tester->print(cout);

  // run the algo
  tester->init();
 
  // run the algo
  tester->run();
  

}
