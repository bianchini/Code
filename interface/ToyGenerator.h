#include "TRandom3.h"
#include "TLorentzRotation.h"
#include "interface/HypoTester.h"
#include "TMath.h"
#include<iostream>

using namespace std;

typedef TLorentzVector LV;

namespace Algo {

  class ToyGenerator{

  public:
    ToyGenerator(const int);
    ~ToyGenerator();
    vector<pair<char,TLorentzVector>> generate( const vector<Decay>& , const int&);
  private:

    int verbose;
    void generate_hypo( vector<pair<char,TLorentzVector>>& , Decay, const int&);
    TRandom3* ran;

  };

}
