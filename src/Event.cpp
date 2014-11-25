#include "interface/Event.h"

Algo::Event::Event( TTree* t ){
  outTree = t;
  createBranches() ;
  reset();
}


Algo::Event::~Event(){ };

void Algo::Event::fillTree(){
  outTree->Fill();
}

void Algo::Event::printTree(){
  treeStruct.print(cout);
}

void Algo::Event::reset(){
  treeStruct.n_h       = 0;
  treeStruct.n_dim     = 0;  
  treeStruct.all_time  = 0;
}

void Algo::Event::createBranches(){  
  outTree->Branch( "test__n_hyp",     &(treeStruct.n_h),       "test__n_hyp/i");
  outTree->Branch( "test__all_time",  &(treeStruct.all_time),  "test__all_time/I");
  outTree->Branch( "test__nll",       &(treeStruct.nll),       "test__nll[test__n_hyp]/D");
  outTree->Branch( "test__status",    &(treeStruct.status),    "test__status[test__n_hyp]/I");
  outTree->Branch( "test__min_time",  &(treeStruct.min_time),  "test__min_time[test__n_hyp]/I");
  outTree->Branch( "test__n_perm",    &(treeStruct.n_perm),    "test__n_perm[test__n_hyp]/i");
  outTree->Branch( "test__n_dim",     &(treeStruct.n_dim),     "test__n_dim/i");
  outTree->Branch( "test__dim",       &(treeStruct.dim),       "test__dim[test__n_dim]/i");
  outTree->Branch( "test__param",     &(treeStruct.param),     "test__param[test__n_dim]/D");
  outTree->Branch( "test__obs",       &(treeStruct.obs),       "test__obs[test__n_dim]/D");
  outTree->Branch( "test__obs_BTAG",  &(treeStruct.obs_BTAG),  "test__obs_BTAG[test__n_dim]/D");
}


                                        
