#ifndef Utils_h
#define Utils_h

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TFormula.h"
#include "TF1.h"
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
#include<chrono>

typedef TLorentzVector LV;

using namespace std;
using namespace std::chrono;

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

  const double BTAG_Q_param[2][2] = 
    { {0.98, 0.02},
      {0.98, 0.02}
    };
  const double BTAG_C_param[2][2] = 
    { {0.80, 0.20},     
      {0.80, 0.20} 
    };    
  const double BTAG_B_param[2][2] =
     { {0.30, 0.70},
       {0.30, 0.70}
     };

  enum QuarkType : int { QuarkTypeUp=0, QuarkTypeDown=1, QuarkTypeBottom=2 };

  double pdf_btag(double*, double*);
  size_t eta_to_bin( const double& );
  size_t eta_to_bin( const LV& );
 
  enum class Decay { TopLep, TopHad, WHad, Higgs, Radiation_u, Radiation_d, Radiation_b, Radiation_g, MET, UNKNOWN };
  string translateDecay(Decay&);

  enum class FinalState { TopLep_l=0, TopLep_b,
      TopHad_q,   TopHad_qbar,  TopHad_b,
      WHad_q,     WHad_qbar,
      Higgs_b, Higgs_bbar,
      Radiation_u, Radiation_d, Radiation_b, Radiation_g };

  struct CompFinalState {
    bool operator()(pair<FinalState,size_t> a, pair<FinalState,size_t> b){
      return (a.first>b.first) || (a.first==b.first && a.second>b.second);
    }
  };
  
  struct Object {   
    LV p4;
    map<string,double> obs; 
    char type;
    void init(const LV& p, const char typ)     { p4 = p; type = typ; }    
    void addObs(const string& name, double val){ obs.insert( make_pair(name, val) ); }    
  };

  bool isSame( const std::vector<std::pair<FinalState,size_t>>&, 
	       const std::vector<std::pair<FinalState,size_t>>&);

  enum Strategy : int { FirstTrial=0, Coarser=1 };

}

#endif
