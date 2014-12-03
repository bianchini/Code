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
  constexpr double PTTHRESHOLD = 30.;
  
  const string TF_Q     = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]*[2]+[3]*[3]/y+[4]*[4]/y/y), 1)";
  const string TF_Q_CUM = "(1 - 0.5*(TMath::Erf(  (x-([0] + [1]*y)) / (y*TMath::Sqrt([2]*[2]+[3]*[3]/y+[4]*[4]/y/y)) ) + 1 ))";

  const string TF_B = 
    "[10]*TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]*[2]+[3]*[3]/y+[4]*[4]/y/y), 1) + "
    "(1-[10])*TMath::Gaus(x, [5] + [6]*y, y*TMath::Sqrt([7]*[7]+[8]*[8]/y+[9]*[9]/y/y), 1)";
  const double TF_Q_param[2][5] =
    { { 0.00e+00, 1.00e+00, 0.00e+00, 1.56e+00, 0.00e+00 },
      { 0.00e+00, 1.00e+00, 1.30e-01, 1.52e+00, 0.00e+00 }
    };
  const double TF_B_param[2][11] = 
    { { -3.60e+00, 1.00e+00, 0.00e+00, 0.99e+00, 5.70e+00,-3.30e+00, 0.94e+00, 0.16e+00, 1.70e+00, 6.60e+00, 0.65e+00 },
      { -4.30e+00, 0.98e+00, 0.00e+00, 1.90e+00, 6.00e+00, 0.91e+01, 0.87e+00, 0.23e+00, 1.10e+00, 0.00e+00, 0.65e+00 },
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
 
  enum class Decay { TopLep, TopHad, TopHadLost, WHad, Higgs, Radiation_u, Radiation_d, Radiation_b, Radiation_g, Lepton, MET, UNKNOWN };
  string translateDecay(Decay&);

  enum class FinalState { 
    TopLep_l,      // lep   from top decay   () 
      TopLep_b,    // b     from top decay   (TF and btag[B])
      TopHad_q,    // q     from top decay   (TF and btag[0.5*LF+0.5*C])
      TopHad_qbar, // qbar  from top decay   (TF and btag[LF]) 
  TopHadLost_qbar, // lost qbar  from top decay   (Acceptance) 
      TopHad_b,    // b     from top decay   (TF and btag[B])
      WHad_q,      // q     from top decay   (TF and btag[0.5*LF+0.5*C])     
      WHad_qbar,   // qbar  from top decay   (TF and btag[LF])
      Higgs_b,     // b     from Higgs decay (TF and btag[B]) 
      Higgs_bbar,  // bbar  from Higgs decay (TF and btag[B]) 
      Radiation_u, // extra quark            (TF and btag[0.5*LF+0.5*C])
      Radiation_d, // extra quark            (TF and btag[LF])
      Radiation_b, // extra quark            (TF and btag[B])  
      Radiation_g  // extra quark            (TF, no flavour assignment)  
      };

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

  bool filter_by_btag( const std::vector<std::pair<FinalState,size_t>>&,
		       const vector<Algo::Object>&  );

  enum Strategy : int { FirstTrial=0, StartFromLastMinimum=1 };

}

#endif
