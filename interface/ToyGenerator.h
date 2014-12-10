#ifndef TOYGENERATOR_H
#define TOYGENERATOR_H

#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TLorentzRotation.h"
#include "interface/Utils.h"
#include "interface/StandardIncludes.h"
#include<sstream>

#include<iostream>

using namespace std;

#define  SMEAR_EREL_MIN     0.2
#define  SMEAR_EREL_MAX      5.
#define  SMEAR_MET_PX_MIN -200.
#define  SMEAR_MET_PX_MAX +200.
#define  SMEAR_MET_PY_MIN -200.
#define  SMEAR_MET_PY_MAX +200.



typedef TLorentzVector LV;

namespace Algo {

  class ToyGenerator{

  public:
    ToyGenerator(const int);
    ToyGenerator(const int, const unsigned int);
    ~ToyGenerator();
    vector<Algo::Object> generate( const vector<Algo::Decay>& , const int&, const int&);

  private:
    int verbose;
    void generate_hypo( vector<Algo::Object>& , Algo::Decay, const int&, const int&);
    void assign_rnd_btag( const Algo::QuarkType, Algo::Object&, const int&);
    void smear_by_TF( TLorentzVector&, const char);
    TRandom3* ran;

  };

}

#endif
