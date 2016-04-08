#include "TTH/MEIntegratorStandalone/interface/FitUtils.h"
#include<iostream>
#include "TH1D.h"
#include "TTree.h"
#include "TRandom3.h"

using namespace std;

int main(int argc, char *argv[]){

  TRandom3* ran = new TRandom3();  

  FitUtils* fit = new FitUtils();
  fit->save_plots(100);

  TTree* tree_sgn = new TTree("tree_sgn", "signal tree");
  TTree* tree_bkg = new TTree("tree_bkg", "background tree");
  TTree* tree_tra = new TTree("tree_tra", "transform tree");
  double x;
  tree_sgn->Branch("x", &x, "x/D");
  tree_bkg->Branch("x", &x, "x/D");
  tree_tra->Branch("x", &x, "x/D");
  for(int i = 0 ; i < 10000; ++i){
    x = ran->Gaus(+2.,2.);
    tree_sgn->Fill();
    x = ran->Gaus(-2.,2.);
    tree_bkg->Fill();
    x = ran->Gaus(-1.5,2.);
    tree_tra->Fill();
  }
  fit->addPdf(tree_tra, Process::Transform);
  fit->addPdf(tree_bkg, Process::Background, 1000.);
  fit->addPdf(tree_sgn, Process::Signal, 200.);
  
  fit->run_test(1000, Option::Binned, 4);

  delete fit;
  delete ran;

}
