#include<iostream>
#include<memory> 

#include "interface/HypoTester.h"

#include "TLorentzVector.h"


using namespace std;

int main(){

  unique_ptr<Algo::HypoTester> tester{ new Algo::HypoTester() };
  
  // test
  tester->push_back_object( TLorentzVector(70.,10.,0., 100 ) , 'j');
  tester->push_back_object( TLorentzVector(70.,10.,0., 100 ) , 'j');
  tester->push_back_object( TLorentzVector(70.,10.,0., 100 ) , 'j');
  tester->push_back_object( TLorentzVector(70.,10.,0., 100 ) , 'j');
  //tester->push_back_object( TLorentzVector(70.,10.,0., 100 ) , 'j');
  //tester->push_back_object( TLorentzVector(70.,10.,0., 100 ) , 'j');
  //tester->push_back_object( TLorentzVector(20.,0.,0.,  21  ) , 'j');
  //tester->push_back_object( TLorentzVector(50.,0.,0.,  50  ) , 'l');
  tester->push_back_object( TLorentzVector(20.,0.,0.,  20  ) , 'm');

  // assumptions
  tester->assume( Algo::Decay::TopHad   );
  //tester->assume( Algo::Decay::TopLep   );
  //tester->assume( Algo::Decay::HiggsHad );

  // run the algo
  tester->run();

  // printout
  tester->print(cout);

}
