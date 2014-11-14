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
  
  const string TF_Q = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]+[3]/y+[4]/y/y), 1)";
  const string TF_B = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]+[3]/y+[4]/y/y), 1)";
  const double TF_Q_param[2][5] =
    { { 0.0e+00, 1.0e+00, 0.0e+00, 1.5e+00, 0.0e+00 },
      { 0.0e+00, 1.0e+00, 1.3e+01, 1.5e+00, 0.0e+00 }
    };
  const double TF_B_param[2][5] =
    { { 0.0e+00, 1.0e+00, 0.0e+00, 1.5e+00, 0.0e+00 },
      { 0.0e+00, 1.0e+00, 1.3e+01, 1.5e+00, 0.0e+00 }
    };

  const string TF_MET = "TMath::Gaus(x,0.,100, 1)*TMath::Gaus(y,0.,100, 1)";
  
  size_t eta_to_bin( const LV& );
 
  enum Decay { TopLep, TopHad, WHad, Higgs, Radiation, MET, UNKNOWN };
  string translateDecay(Decay&);


  enum FinalState { TopLep_l=0, TopLep_b,
                    TopHad_q,   TopHad_qbar,  TopHad_b,
                    WHad_q,     WHad_qbar,
                    Higgs_b, Higgs_bbar,
                    Radiation_q };


  bool isSame( const std::vector<std::pair<FinalState,size_t>>&, const std::vector<std::pair<FinalState,size_t>>&);

}
