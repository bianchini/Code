#ifndef PARAMETERS_H
#define PARAMETERS_H

#ifndef DEBUG_MODE 
#define DEBUG_MODE
#endif

// std library
#include<iostream>
#include<map>
#include<unordered_map>
#include<assert.h>
#include<vector>
#include<bitset> 
#include<algorithm>  
#include<initializer_list>
#include<limits> 
#include<chrono>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF1.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"

typedef TLorentzVector LV;
using namespace std;
using namespace std::chrono;

namespace MEM {
  
  size_t eta_to_bin( const double& );

  constexpr double PI   = 3.14159265359;
  constexpr double MTOP = 174.3;
  constexpr double MTOP2= MTOP*MTOP;
  constexpr double MTOP4= MTOP2*MTOP2;
  constexpr double GTOP = 2.0;
  constexpr double MB   = 4.8;
  constexpr double MB2  = MB*MB;
  constexpr double MUB  = MB2/MTOP2;
  constexpr double MQ   = 0.;
  constexpr double ML   = 0.;
  constexpr double MW   = 80.19;
  constexpr double MW2  = MW*MW;
  constexpr double GW   = 2.08;
  constexpr double DMT2 = (MTOP*MTOP-MB*MB-MW*MW)*0.5;
  constexpr double MH   = 125.;
  constexpr double MH2  = MH*MH;
  constexpr double PSHBB= 1-4*MB2/MH2;
  constexpr double GH   = 0.00407;
  constexpr double DMH2 = (MH*MH-2*MB*MB)*0.5;
  constexpr double DMW2 = (MW*MW)*0.5;
  constexpr double GFERMI = 1.166e-05;
  constexpr double GEWK2  = GFERMI/sqrt(2)*8*MW2;
  constexpr double GEWK4  = GEWK2*GEWK2;
  constexpr double YB     = MB*sqrt(sqrt(2)*GFERMI);
  constexpr double YB2    = YB*YB;
  constexpr double BWTOP  = PI/MTOP/GTOP;
  constexpr double BWH    = PI/MH/GH;
    
  const double TF_Q_param[2][5] =
    { { 0.00e+00, 1.00e+00, 0.00e+00, 1.56e+00, 0.00e+00 },
      { 0.00e+00, 1.00e+00, 1.30e-01, 1.52e+00, 0.00e+00 }
    };
  const double TF_B_param[2][11] = 
    { { -3.60e+00, 1.00e+00, 0.00e+00, 0.99e+00, 5.70e+00,-3.30e+00, 0.94e+00, 0.16e+00, 1.70e+00, 6.60e+00, 0.65e+00 },
      { -4.30e+00, 0.98e+00, 0.00e+00, 1.90e+00, 6.00e+00, 0.91e+01, 0.87e+00, 0.23e+00, 1.10e+00, 0.00e+00, 0.65e+00 },
    };
  const double TF_MET_param[2] = {20., 20.};
  const double TF_ACC_param[3] = {2.5, 30., 1.0};

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

  enum DebugVerbosity { silent=1, input=2, init=4, init_more=8, event=16, integration=32};
 

  enum class TFType {bReco=0, qReco=1, bLost=2, qLost=3, muReco=4, elReco=5, MET=6, Recoil=7, Unknown=8};

  bool isQuark   (const TFType&);
  bool isNeutrino(const TFType&);
  bool isLepton  (const TFType&);
  double transfer_function( double*,  double*, const TFType&, const int&);

  pair<double, double> get_support( double*, const TFType&, const double&, const int&);

  enum class ObjectType { Jet=0, Lepton=1, MET=2, Recoil=3};

  enum class Observable { E_LOW_Q=0, E_HIGH_Q=1, E_LOW_B=2, E_HIGH_B=3, BTAG=4, CHARGE=5};

  class ObsHash{
  public:
    std::size_t operator()(const Observable& s) const {
      std::size_t h1 = std::hash<size_t>()(static_cast<size_t>(s));
      return h1;
    }
  };

  class ObsEqual{
  public:
    bool operator()( const Observable& a, const Observable& b ) const {
      return static_cast<size_t>(a)==static_cast<size_t>(b);
    }
  };

  class Object {       
  public:
    Object(const LV&, const ObjectType&); 
    ~Object(); 
    LV p4() const;
    double getObs(const Observable&) const; 
    bool isSet(const Observable&) const;
    void addObs(const Observable&, const double&);
    void print(ostream& os) const;
  private:
    LV p;
    ObjectType type;
    std::unordered_map<Observable,double, ObsHash, ObsEqual > obs; 
  };  
 
  
  enum class PSVar { E_q1=0,     cos_q1=1,     phi_q1=2,  
      E_qbar1=3,  cos_qbar1=4,  phi_qbar1=5,  
      E_b1=6,     cos_b1=7,     phi_b1=8,
      E_q2=9,     cos_q2=10,    phi_q2=11,  
      E_qbar2=12, cos_qbar2=13, phi_qbar2=14,  
      E_b2=15,    cos_b2=16,    phi_b2=17,
      E_b=18,     cos_b=19,     phi_b=20,  
      E_bbar=21,  cos_bbar=22,  phi_bbar=23,
      P_t=24,     cos_t=25,     phi_t=26,
      P_tbar=27,  cos_tbar=28,  phi_tbar=29,
      P_h=30,     cos_h=31,     phi_h=32,
      Px_h=33,    Py_h=34,      Pz_h=35};
  
  enum class PSPart {
      q1=0, qbar1=1, b1=2,
      q2=3, qbar2=4, b2=5,
      b =6, bbar =7, 
      t=8,  tbar=9, h=10 };

  class PSVarHash{
  public:
    std::size_t operator()(const PSVar& s) const {
      std::size_t h1 = std::hash<size_t>()(static_cast<size_t>(s));
      return h1;
    }
  };

  class PSVarEqual{
  public:
    bool operator()( const PSVar& a, const PSVar& b ) const {
      return static_cast<size_t>(a)==static_cast<size_t>(b);
    }
  };

  class PSPartHash{
  public:
    std::size_t operator()(PSPart const& s) const {
      std::size_t h1 = std::hash<size_t>()(static_cast<size_t>(s));
      return h1;
    }
  };

  class PSPartEqual{
  public:
    bool operator()( const PSPart& a, const PSPart& b ) const {
      return static_cast<size_t>(a)==static_cast<size_t>(b);
    }
  };


  struct GenPart {
    GenPart(){ lv = LV(); type = TFType::Unknown;}
    GenPart(const LV& a, const TFType& b, const int c){ lv = a; type = b; charge=c;}    
    LV lv;
    TFType type;
    int charge;
  };

  typedef std::unordered_map<MEM::PSPart, MEM::GenPart, MEM::PSPartHash, MEM::PSPartEqual> PSMap;
  //typedef std::map<MEM::PSPart, MEM::GenPart> PSMap; 

  class PS {
  public: 
    PS( size_t=0);
    ~PS();
    PSMap::const_iterator begin() const;
    PSMap::const_iterator end() const;
    LV lv(const PSPart&) const;  
    int charge(const PSPart& ) const;
    TFType type(const PSPart&) const;
    void set(const PSPart&, const GenPart&);    
    void print(ostream&) const;
  private:
    size_t dim;
    PSMap val;
  };
  
  enum class Hypothesis { TTH=0, TTBB=1, Undefined=3 };

  enum class FinalState { HH=0, LH=1, LL=2, TTH=3, Undefined=4};

  enum class Assumption { ZeroQuarkLost=0, OneQuarkLost, TwoQuarkLost};

  enum class Permutations { BTagged=0, QUntagged, QQbarSymmetry, BBbarSymmetry};

  enum IntegrandType { Constant=1, Jacobian=2, Transfer=4, ScattAmpl=8, DecayAmpl=16, PDF=32 };

  struct CompPerm {
    bool operator()(int a, int b){
      return a>b;
    }
  };
  
}

#endif
