#include<iostream>
#include<memory>
#include<cstdlib>

#include "interface/ToyGenerator.h"
//#include "interface/HypoTester.h"

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"


using namespace std;

int main(int argc, char *argv[]){

  TTree* t    = new TTree("tree", "tree");
  Algo::HypoTester* tester = new Algo::HypoTester(t) ;

  Algo::ToyGenerator* toyGenerator = new Algo::ToyGenerator(0);

  vector<Algo::Decay> decays;
  decays.push_back( Algo::Decay::TopHad );
  decays.push_back( Algo::Decay::WHad  );

  const int ntoy  = argc>1 ? atoi(argv[1]) : 1 ;
  const int smear = argc>2 ? atoi(argv[2]) : 0 ;

  int itoy = 0;
  while( itoy < ntoy ){
    vector<pair<char,TLorentzVector>> out = 
      toyGenerator->generate( decays , smear );

    int count_j = 0;
    int count_l = 0;
    int count_m = 0;
    TLorentzVector invisible(0.,0.,0.,0.);

    for( auto fs : out ){
      if( fs.first=='j' && fs.second.Pt()>30 && TMath::Abs(fs.second.Eta())<2.5 ) ++count_j;
      if( fs.first=='l' && fs.second.Pt()>20 && TMath::Abs(fs.second.Eta())<2.5 ) ++count_l;
      if( fs.first=='m' && fs.second.Pt()>0. ){
	invisible += fs.second;
	++count_m;
      }
    }
    cout << "Generate event " << itoy << "/" << ntoy << endl;
    cout << "\tNumber of jets: " << count_j << endl;
    cout << "\tNumber of lept: " << count_l << endl;
    cout << "\tNumber of nus : " << count_m << endl;

    if( count_j==5 ) {
      for( auto fs : out ){
	if( fs.first=='j' && fs.second.Pt()>30 && TMath::Abs(fs.second.Eta())<2.5 ){
	  tester->push_back_object( fs.second  , 'j');
	}
      }

      map<string, vector<Algo::Decay> > hypotheses;
      hypotheses["H0"] = {Algo::Decay::TopHad, Algo::Decay::Higgs};
      hypotheses["H1"] = {Algo::Decay::TopHad, Algo::Decay::Radiation, Algo::Decay::Radiation};
      tester->test( hypotheses );

      ++itoy;
    }

  }

  TFile* fout = new TFile("test/Test_3.root", "RECREATE");
  fout->cd();
  t->Write("", TObject::kOverwrite);
  fout->Close();
  cout << "ROOT file saved" << endl;

  delete toyGenerator;
  delete tester;
}
