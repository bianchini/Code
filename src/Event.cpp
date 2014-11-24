#include "interface/Event.h"

Algo::Event::Event( TTree* t ){
  outTree = t;
  createBranches() ;
  reset();
}


Algo::Event::~Event(){ };

void Algo::Event::fillTree(){
  treeStruct.print(cout);
  outTree->Fill();
}

void Algo::Event::reset(){
  treeStruct.n_h   = 0;
  treeStruct.n_dim = 0;  
}

void Algo::Event::createBranches(){  
  outTree->Branch( "test__n_hyp",   &(treeStruct.n_h),   "test__n_hyp/i");
  outTree->Branch( "test__nll",     &(treeStruct.nll),   "test__nll[test__n_hyp]/D");
  outTree->Branch( "test__n_dim",   &(treeStruct.n_dim), "test__n_dim/i");
  outTree->Branch( "test__dim",     &(treeStruct.dim),   "test__dim[test__n_dim]/i");
  outTree->Branch( "test__param",   &(treeStruct.param), "test__param[test__n_dim]/D");
  outTree->Branch( "test__obs",     &(treeStruct.obs),   "test__obs[test__n_dim]/D");
}


                                        
