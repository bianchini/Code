#ifndef Utils_h
#define Utils_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TFormula.h"
#include "TMath.h"

#include<cmath>
#include<cstddef>
#include<iostream>
#include<cstddef>
#include<cstdio>
#include<map> 
#include<string>
#include<algorithm>
#include<assert.h>
#include<memory> 
#include<limits> 

typedef TLorentzVector LV;

using namespace std;

namespace Algo {

  constexpr double MTOP = 174.3;
  constexpr double MB   = 4.8;
  constexpr double MW   = 80.19;
  constexpr double DM2  = (MTOP*MTOP-MB*MB-MW*MW)*0.5;
  constexpr double MH   = 125.;
  constexpr double DMH2 = (MH*MH-2*MB*MB)*0.5;
  
  const string TF_Q = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]*[2]+[3]*[3]/y+[4]*[4]/y/y), 1)";
  const string TF_B = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]*[2]+[3]*[3]/y+[4]*[4]/y/y), 1)";
  const double TF_Q_param[2][5] =
    { { 0.0e+00, 1.0e+00, 0.0e+00, 1.5e+00, 0.0e+00 },
      { 0.0e+00, 1.0e+00, 1.3e-01, 1.5e+00, 0.0e+00 }
    };
  const double TF_B_param[2][5] =
    { { 0.0e+00, 1.0e+00, 0.0e+00, 1.5e+00, 0.0e+00 },
      { 0.0e+00, 1.0e+00, 1.3e-01, 1.5e+00, 0.0e+00 }
    };
  const string TF_MET = "TMath::Gaus(x,0.,20, 1)*TMath::Gaus(y,0.,20, 1)";
  
  const string BTAG_U = 
    "([0]>30. && [0]<9999.)&&([1]>=0.0 && [1]< 1.0)*( (x<0.5)*(0.98*0.5 + 0.80*0.5) + (x>0.5)*(0.02*0.5 + 0.20*0.5) )"
    "+"
    "([0]>30. && [0]<9999.)&&([1]>=1.0 && [1]< 2.5)*( (x<0.5)*(0.98*0.5 + 0.80*0.5) + (x>0.5)*(0.02*0.5 + 0.20*0.5) )";
  const string BTAG_B = 
    "([0]>30. && [0]<9999.)&&([1]>=0.0 && [1]< 1.0)*( (x<0.5)*0.30 + (x>0.5)*0.70 )"
    "+"
    "([0]>30. && [0]<9999.)&&([1]>=1.0 && [1]<=2.5)*( (x<0.5)*0.30 + (x>0.5)*0.70 )";
  const string BTAG_D =
    "([0]>30. && [0]<9999.)&&([1]>=0.0 && [1]< 1.0)*( (x<0.5)*(0.98*0.5 + 0.98*0.5) + (x>0.5)*(0.02*0.5 + 0.02*0.5) )"
    "+"
    "([0]>30. && [0]<9999.)&&([1]>=1.0 && [1]< 2.5)*( (x<0.5)*(0.98*0.5 + 0.98*0.5) + (x>0.5)*(0.02*0.5 + 0.02*0.5) )";


  size_t eta_to_bin( const LV& );
 
  enum class Decay { TopLep, TopHad, WHad, Higgs, Radiation_q, Radiation_b, MET, UNKNOWN };
  string translateDecay(Decay&);


  enum class FinalState { TopLep_l=0, TopLep_b,
      TopHad_q,   TopHad_qbar,  TopHad_b,
      WHad_q,     WHad_qbar,
      Higgs_b, Higgs_bbar,
      Radiation_q, Radiation_b };

  // compare hypos                                                                                                                                   
  struct CompFinalState {
    bool operator()(pair<FinalState,size_t> a, pair<FinalState,size_t> b){
      return (a.first>b.first) || (a.first==b.first && a.second>b.second);
    }
  };
  
  struct Object {
    
    LV p4;                   // the four momentum                                                                                                    
    map<string,double> obs;  // observables                                                                                                          
    char type;

    void init(const LV& p, const char typ){
      p4   = p;
      type = typ;
    }    
    void addObs(const string& name, double val){
      obs.insert( make_pair(name, val) );
    }
    
  };

  bool isSame( const std::vector<std::pair<FinalState,size_t>>&, const std::vector<std::pair<FinalState,size_t>>&);

}

#endif
