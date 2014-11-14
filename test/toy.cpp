#include<iostream>
#include<memory>
#include<cstdlib>
#include <getopt.h>

#include "interface/ToyGenerator.h"
#include "interface/HypoTester.h"

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"


using namespace std;

int main(int argc, char *argv[]){

  string fname = "DUMMY.root";
  int ntoys    = 0;  
  int smear    = 0;
  int verbose  = 0;

  static struct option long_options[] = {
    {"toys",     required_argument,   0, 't' },
    {"output",   required_argument,   0, 'o' },
    {"smear",    required_argument,   0, 's' },
    {"verbose",  optional_argument,   0, 'v' },
    {0, 0, 0, 0 }
  };

  int long_index  = 0;
  int opt;
  while ((opt = getopt_long(argc, argv,"t:o:s:v:",
			    long_options, &long_index )) != -1) {
    switch (opt) {
    case 't' : ntoys   = atoi(optarg);
      break;
    case 'o' : fname   = string(optarg);
      break;
    case 's' : smear   = atoi(optarg);
      break;
    case 'v' : verbose = atoi(optarg);                                 
      break;  
    case '?':
      cout << "Unknown option " << opt << endl;
      break;
    default:
      cout << "Invalid usage of toy" << endl;
      return 0;
    }
  }


  TTree* t = new TTree("tree", "tree");
  Algo::HypoTester* tester = new Algo::HypoTester(t) ;

  Algo::ToyGenerator* toyGenerator = new Algo::ToyGenerator(verbose);

  vector<Algo::Decay> decays;
  decays.push_back( Algo::Decay::TopHad );
  decays.push_back( Algo::Decay::WHad  );

  //const int ntoy  = argc>1 ? atoi(argv[1]) : 1 ;
  //const int smear = argc>2 ? atoi(argv[2]) : 0 ;

  int itoy = 0;
  while( itoy < ntoys ){
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
    cout << "Generate event " << itoy << "/" << ntoys << endl;
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

  TFile* fout = new TFile(("test/"+fname).c_str(), "RECREATE");
  fout->cd();
  t->Write("", TObject::kOverwrite);
  fout->Close();
  cout << "ROOT file saved" << endl;

  delete toyGenerator;
  delete tester;
}
