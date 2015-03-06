#include "interface/Event.h"

Algo::Event::Event( TTree* t ){
  outTree = t;
  createBranches() ;
  reset();
}


Algo::Event::Event(){
  outTree = nullptr;
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
  treeStruct.n_btag    = 0; 
  treeStruct.n_jet     = 0; 
  treeStruct.n_lep     = 0; 
}

void Algo::Event::createBranches(){  
  outTree->Branch( "test__n_hyp",     &(treeStruct.n_h),       "test__n_hyp/i");
  outTree->Branch( "test__n_dim",     &(treeStruct.n_dim),     "test__n_dim/i");
  outTree->Branch( "test__all_time",  &(treeStruct.all_time),  "test__all_time/I");
  outTree->Branch( "test__n_btag",    &(treeStruct.n_btag),    "test__n_btag/I");
  outTree->Branch( "test__n_jet",     &(treeStruct.n_jet),     "test__n_jet/I");
  outTree->Branch( "test__n_lep",     &(treeStruct.n_lep),     "test__n_lep/I");
  outTree->Branch( "test__nll",       &(treeStruct.nll),       "test__nll[test__n_hyp]/D");
  outTree->Branch( "test__status",    &(treeStruct.status),    "test__status[test__n_hyp]/I");
  outTree->Branch( "test__strategy",  &(treeStruct.strategy),  "test__strategy[test__n_hyp]/I");
  outTree->Branch( "test__min_time",  &(treeStruct.min_time),  "test__min_time[test__n_hyp]/I");
  outTree->Branch( "test__dim",       &(treeStruct.dim),       "test__dim[test__n_hyp]/I");
  outTree->Branch( "test__perm",      &(treeStruct.perm),      "test__perm[test__n_hyp]/I");
  outTree->Branch( "test__param",     &(treeStruct.param),     "test__param[test__n_dim]/D");
  outTree->Branch( "test__obs_e",     &(treeStruct.obs_e),     "test__obs_e[test__n_dim]/D");
  outTree->Branch( "test__obs_pt",    &(treeStruct.obs_pt),    "test__obs_pt[test__n_dim]/D");
  outTree->Branch( "test__obs_eta",   &(treeStruct.obs_eta),   "test__obs_eta[test__n_dim]/D");
  outTree->Branch( "test__obs_phi",   &(treeStruct.obs_phi),   "test__obs_phi[test__n_dim]/D");
  outTree->Branch( "test__obs_btag",  &(treeStruct.obs_btag),  "test__obs_btag[test__n_dim]/D");
}


                                        
