#ifndef PARAMETERS_H
#define PARAMETERS_H

#ifndef DEBUG_MODE 
#define DEBUG_MODE
#endif

// std library
#include<iostream>
#include<map>
#include<string>


#include<assert.h>

//This is not C++11 safe, we need to run gccxml which does not support it
//#include <unordered_map>

//This is C++11 safe
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
 
#include<assert.h>
#include<vector>
#include<bitset> 
#include<algorithm>
#include<limits>

//Not available in old C++
#ifndef __GCCXML__
#include<chrono>
#endif


// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF1.h"
#include "TH3D.h"
#include "TF2.h"
#include "TRandom3.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"

typedef TLorentzVector LV;
using namespace std;

#ifndef __GCCXML__
using namespace std::chrono;
#endif

namespace MEM {
  
  int eta_to_bin( const double&, bool =false );
  double deltaR( const LV&, const LV&);

  bool descending (double, double);

  std::vector<std::size_t> get_sorted_indexes(const std::vector<double>&, const double&);
  bool is_in( const vector<std::size_t>&, const std::size_t&);

  const double PI   = 3.14159265359;
  const double MTOP = 174.3;
  const double MTOP2= MTOP*MTOP;
  const double MTOP4= MTOP2*MTOP2;
  const double GTOP = 2.0;
  const double MB   = 4.8;
  const double MB2  = MB*MB;
  const double MUB  = MB2/MTOP2;
  const double MQ   = 0.;
  const double ML   = 0.;
  const double MW   = 80.19;
  const double MW2  = MW*MW;
  const double GW   = 2.08;
  const double DMT2 = (MTOP*MTOP-MB*MB-MW*MW)*0.5;
  const double MH   = 125.;
  const double MH2  = MH*MH;
  const double PSHBB= 1-4*MB2/MH2;
  const double GH   = 0.00407;
  const double DMH2 = (MH*MH-2*MB*MB)*0.5;
  const double DMW2 = (MW*MW)*0.5;
  const double GFERMI = 1.166e-05;
  const double GEWK2  = GFERMI/sqrt(2)*8*MW2;
  const double GEWK4  = GEWK2*GEWK2;
  const double YB     = MB*sqrt(sqrt(2)*GFERMI);
  const double YB2    = YB*YB;
  const double BWTOP  = PI/MTOP/GTOP;
  const double BWH    = PI/MH/GH;
    
  const double TF_Q_param[2][5] =
    { { 0.00e+00, 1.00e+00, 0.00e+00, 1.56e+00, 0.00e+00 },
      { 0.00e+00, 1.00e+00, 1.30e-01, 1.52e+00, 0.00e+00 }
    };
  const double TF_B_param[2][11] = 
    { { -3.60e+00, 1.00e+00, 0.00e+00, 0.99e+00, 5.70e+00,-3.30e+00, 0.94e+00, 0.16e+00, 1.70e+00, 6.60e+00, 0.65e+00 },
      { -4.30e+00, 0.98e+00, 0.00e+00, 1.90e+00, 6.00e+00, 0.91e+01, 0.87e+00, 0.23e+00, 1.10e+00, 0.00e+00, 0.65e+00 },
    };
  const double TF_MET_param   [3] = {30.,  30., 0.};
  const double TF_RECOIL_param[3] = {4.1,  1.35,  9999.0};
  const double TF_ACC_param   [3] = {2.5, 30.};

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
     
     
  //int getEtaBin(double);

  enum DebugVerbosity { output=1, input=2, init=4, init_more=8, event=16, integration=32};
 
  namespace TFType {
    enum TFType {bReco=0, qReco=1, bLost=2, qLost=3, muReco=4, elReco=5, MET=6, Recoil=7, Unknown=8};
  }
  
  namespace DistributionType {
    enum DistributionType {csv_b=0, csv_c=1, csv_l=2};
  }
  
  namespace TFMethod {
    enum TFMethod {Builtin=0, External=1, Unknown=2};
  }
  
  bool isQuark   (const TFType::TFType&);
  bool isNeutrino(const TFType::TFType&);
  bool isLepton  (const TFType::TFType&);
  double Chi2Corr(const double&, const double&, const double&, const double&, const double&);
  double Chi2(const double&, const double&, const double&);

  double transfer_function( double*,  double*, const TFType::TFType&, int&, const double&, const int&);

  double transfer_function_smear(double*, double* );
  
  namespace ObjectType {
    enum ObjectType { Jet=0, Lepton=1, MET=2, Recoil=3, Unknown=4};
  } 
  
  namespace Observable {
    enum Observable { E_LOW_Q=0, E_HIGH_Q=1, E_LOW_B=2, E_HIGH_B=3, BTAG=4, CHARGE=5, PDGID=6, CSV=7};
  }
  
  class ObsHash{
  public:
    std::size_t operator()(const Observable::Observable& s) const {
      std::size_t h1 = boost::hash<std::size_t>()(static_cast<std::size_t>(s));
      return h1;
    }
  };

  class ObsEqual{
  public:
    bool operator()( const Observable::Observable& a, const Observable::Observable& b ) const {
      return static_cast<std::size_t>(a)==static_cast<std::size_t>(b);
    }
  };
  
  class TFTypeHash{
  public:
    std::size_t operator()(const TFType::TFType& s) const {
      std::size_t h1 = boost::hash<std::size_t>()(static_cast<std::size_t>(s));
      return h1;
    }
  };

  class TFTypeEqual{
  public:
    bool operator()( const TFType::TFType& a, const TFType::TFType& b ) const {
      return static_cast<std::size_t>(a)==static_cast<std::size_t>(b);
    }
  };

  class Object {       
  public:
    Object(const LV&, const ObjectType::ObjectType&); 
    Object(); 
    ~Object(); 
    LV p4() const;
    void setp4(const LV&);
    ObjectType::ObjectType type() const;
    double getObs(const Observable::Observable&) const; 
    TF1* getTransferFunction(const TFType::TFType&); 
    std::size_t getNumTransferFunctions() const; 
    bool isSet(const Observable::Observable&) const;
    void addObs(const Observable::Observable&, const double&);
    void addTransferFunction(const TFType::TFType&, TF1*);
    void print(ostream& os) const;
  private:
    LV p;
    ObjectType::ObjectType t;
    boost::unordered_map<const Observable::Observable, double, ObsHash, ObsEqual> obs; 
    boost::unordered_map<const TFType::TFType, TF1*, TFTypeHash, TFTypeEqual> transfer_funcs; 
  };  
  pair<double, double> get_support( double*, const TFType::TFType&, const double&, const int&, Object* = nullptr);

  namespace PSVar {
  enum PSVar { E_q1=0,     cos_q1=1,     phi_q1=2,  
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
  }
  
  namespace PSPart {
  enum PSPart {
      q1=0, qbar1=1, b1=2,
      q2=3, qbar2=4, b2=5,
      b =6, bbar =7, 
      t=8,  tbar=9, h=10 };
  }
  
  class PSVarHash{
  public:
    std::size_t operator()(const PSVar::PSVar& s) const {
      std::size_t h1 = boost::hash<std::size_t>()(static_cast<std::size_t>(s));
      return h1;
    }
  };

  class PSVarEqual{
  public:
    typedef PSVar::PSVar value_type;
    bool operator()( const PSVar::PSVar& a, const PSVar::PSVar& b ) const {
      return static_cast<std::size_t>(a)==static_cast<std::size_t>(b);
    }
  };

  class PSPartHash{
  public:
    std::size_t operator()(PSPart::PSPart const& s) const {
      std::size_t h1 = boost::hash<std::size_t>()(static_cast<std::size_t>(s));
      return h1;
    }
  };

  class PSPartEqual{
  public:
    bool operator()( const PSPart::PSPart& a, const PSPart::PSPart& b ) const {
      return static_cast<std::size_t>(a)==static_cast<std::size_t>(b);
    }
  };

  struct GenPart {
    GenPart(){ lv = LV(); type = TFType::Unknown;}
    GenPart(const LV& a, const TFType::TFType& b, const int c){ lv = a; type = b; charge=c;}    
    LV lv;
    TFType::TFType type;
    int charge;
  };

  typedef boost::unordered_map<MEM::PSPart::PSPart, MEM::GenPart, MEM::PSPartHash, MEM::PSPartEqual> PSMap;

  class PS {
  public: 
    PS( std::size_t=0);
    ~PS();
    PSMap::const_iterator begin() const;
    PSMap::const_iterator end() const;
    LV lv(const PSPart::PSPart&) const;  
    int charge(const PSPart::PSPart& ) const;
    TFType::TFType type(const PSPart::PSPart&) const;
    void set(const PSPart::PSPart&, const GenPart&);    
    void print(ostream&) const;
  private:
    std::size_t dim;
    PSMap val;
  };
  
  namespace Hypothesis {
    enum Hypothesis { TTH=0, TTBB=1, Undefined=2 };
  }

  namespace PermConstants {
    enum PermConstants { btag_TTBB=0, btag_TTCC=1, btag_TTJJ=2, VarTransf=4 };
  }
  
  namespace FinalState {
    enum FinalState { HH=0, LH=1, LL=2, TTH=3, Undefined=4};
  }
  
  namespace Assumption {
    enum Assumption { ZeroQuarkLost=0, OneQuarkLost=1, TwoQuarkLost=2};
  }
  
  namespace Permutations {
    enum Permutations { BTagged=0, QUntagged, QQbarSymmetry, BBbarSymmetry, QQbarBBbarSymmetry, HEPTopTagged, HEPTopTaggedNoPrefix, HiggsTagged, FirstRankedByBTAG, FirstTwoRankedByBTAG, FirstThreeRankedByBTAG};
  }
  
  namespace IntegrandType {
    enum IntegrandType { Constant=1, Jacobian=2, Transfer=4, ScattAmpl=8, DecayAmpl=16, PDF=32, Sudakov=64, Recoil=128, SmearJets=256, SmearMET=512 };
  }

  struct CompPerm {
    CompPerm(int order=0){
      highpt_first = order;
    }
    CompPerm& operator=(const CompPerm& cmp){
      highpt_first = cmp.highpt_first;
      return *this;
    }
    bool operator()(int a, int b){
      if(highpt_first)
	return a<b;
      else
	return a>b;
    }
    int highpt_first;
  };

  
  struct MEMConfig{
    MEMConfig( int    =4000,             // num of int points
	       double =1.e-12,           // absolute tol.
	       double =1.e-5,            // relative tol.
	       int    =0,                // int_code
	       int    =0,                // =0 <=> Int{ Perm }; =1 <=> Perm{ Int }
	       double =13000.,           // c.o.m. energy
	       double =8000.,            // max energy for integration over momenta
	       string ="cteq65.LHgrid",  // PDF set
	       double =0.98,             // light quark energy CL
	       double =0.98,             // heavy quark energy CL
	       double =0.98,             // nu phi CL
	       int    =0,                // skip matrix evaluation if some TF are evaluated art chi2>...
	       double =6.6,              // ... ( <=> TMath::ChisquareQuantile(0.99, 1)=6.6 )
	       bool   =false,            // restrict tf to same range used for quark energy integration
	       int    =1,                // use highest pT jets for E_q/E_b,
               TFMethod::TFMethod =TFMethod::Builtin,
	       int    =0,                // do minimisation instead of integration
	       int    =0,                // do runtime pruning of permutations
	       double =1e-03,            // pruning accuracy  
	       int    =0                 // prefit
	       );

    void defaultCfg(float nCallsMultiplier=1.0);
    void setNCalls(FinalState::FinalState, Hypothesis::Hypothesis, Assumption::Assumption, int);
    int getNCalls(FinalState::FinalState, Hypothesis::Hypothesis, Assumption::Assumption);
    void set_tf_global(TFType::TFType type, int etabin, TF1 tf);
    void add_distribution_global(DistributionType::DistributionType type, TH3D tf);

    // optionally this can be called instead of the built-in array
    int n_max_calls;

    // "map" between an integration type and the number of function calls
    // FinalState vs Hypothesis vs Assumption
    int calls[4][2][3];

    // the VEGAS options
    double rel;
    double abs;

    // what to include into the integrand
    int int_code;
    
    // sqrt of S and maximum quark energy
    double sqrts;
    double emax;

    // pdf set (to be initiated once)
    string pdfset;

    // strategy to prune permutations
    std::vector<Permutations::Permutations> perm_pruning;

    // if true, use n_max_calls instead of the built-in array
    bool is_default;    
    
    // CL of jet and met range
    double j_range_CL;
    double b_range_CL;
    double m_range_CL;
   
    // do sum over permutations inside or outside integral
    int perm_int;

    // the number of jets for which the tf,
    // being evaluated at a too unlikely value,
    // triggers the return of a 0.
    int    tf_suppress;

    // the maximum value of a chi2 in a tf evaluation
    // such that the tf is truncated
    double tf_offscale;

    // restrict the tf to be in the same range
    // used for integration over the quark energy
    bool tf_in_range;    

    // use high pT jets first
    int highpt_first;
    
    TFMethod::TFMethod transfer_function_method;

    // do minimzation, not integration
    int do_minimize;
    
    // do runtime pruning of permutations
    int do_perm_filtering;

    // pruning accuracy
    double perm_filtering_rel;

    // do a pre-fit to filter permutations
    int do_prefit;

    std::map<std::pair<TFType::TFType, int>, TF1> tf_map;
    
    std::map<DistributionType::DistributionType, TH3D> btag_pdfs;
  };

  struct MEMOutput{
    double p;
    double p_err;
    double chi2;
    int    error_code;
    int    time;	
    int    num_max_calls;	
    int    num_calls;
    float  efficiency;
    int    prefit_code;
    double btag_weights[3];
    std::size_t num_perm;    
    std::size_t assumption;
    FinalState::FinalState final_state;
    Hypothesis::Hypothesis hypothesis;
    void print(std::ostream& os){
      os.precision(3);
      os << "\t**************** MEM output ****************" << endl;
      os << "\tProbability             = (" << p << " +/- " << p_err << ")" << endl;
      os.precision(2);
      os << "\tRelative precision      = " << (p_err/p)*100 << "%" << endl;
      os.precision(3);
      os << "\tChi2                    = " << chi2 << endl;
      os << "\tFinal state             = " << static_cast<std::size_t>(final_state) << endl;
      os << "\tHypothesis              = " << static_cast<std::size_t>(hypothesis) << endl;
      os << "\tNumber of unreco's jets = " << assumption << endl;
      os << "\tNumber of permutations  = " << num_perm << endl;
      os << "\tTotal number of calls   = " << num_calls << endl;
      os << "\tMaximum number of calls = " << num_max_calls << endl;
      os << "\tPhase-space efficiency  = " << efficiency*100 << "%" << endl;
      os << "\tError code              = " << error_code << endl;
      os << "\tPre-fit code            = " << prefit_code << endl;
      os << "\tB-tag weights           = " 
	 << "(" <<  btag_weights[0] << ", " <<  btag_weights[1] << ", " <<  btag_weights[2] << ")" << endl;
      os << "\tJob done in..............." << time*0.001 << " seconds" << endl;
      os << "\t********************************************" << endl;
      os.precision(8);
    }
  };
  
  double transfer_function2(void*, const double*, const TFType::TFType&, int&, const double&, const bool&, const int&);
  
}

#endif
