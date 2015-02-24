#ifndef INTEGRAND_H
#define INTEGRAND_H

// std library
#include<iostream>
#include<map>
#include<assert.h>
#include<vector>
#include<bitset> 
#include<algorithm>  
#include<initializer_list>

// ROOT
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"

// user headers
#include "interface/Utils.h"


typedef TLorentzVector LV;
using namespace std;

namespace MEM {

  struct Object {   
    enum Type { jet=0, lepton, met, recoil};

    LV p4;
    Type type;
    map<string,double> obs; 
    void init(const LV& p, const Type t) { p4 = p; type = t; }    
    void addObs(const string& name, double val){ obs.insert( make_pair(name, val) ); }    
    bool yesNo(const string& var) const{
      if( obs.find(var)==obs.end() )        return false;
      if( obs.find(var+"_RND")!=obs.end() ) return false;
      return true;
    }
    void print(ostream& os){
      os << "\tType: " << static_cast<int>(type) << ", p4=(Pt, Eta, Phi, M)=(" 
	 << p4.Pt() << ", " << p4.Eta() << ", " << p4.Phi() << ", " << p4.M() 
	 << ")" << endl; 
    }
  };  
  
  enum DebugVerbosity { silent=0, input=2, init=4, event=8, integration=16};
  
  enum class PSVar { E_q1=0,     cos_q1=1,     phi_q1=2,  
	       E_qbar1=3,  cos_qbar1=4,  phi_qbar1=5,  
	       E_b1=6,     cos_b1=7,     phi_b1=8,
	       E_q2=9,     cos_q2=10,    phi_q2=11,  
	       E_qbar2=12, cos_qbar2=13, phi_qbar2=14,  
	       E_b2=15,    cos_b2=16,    phi_b2=17,
	       E_b=18,     cos_b=19,     phi_b=20,  
	       E_bbar=21,  cos_bbar=22,  phi_bbar=23};
  
  struct PS {
    double E_q1, cos_q1, phi_q1,  
      E_qbar1,  cos_qbar1,  phi_qbar1,  
      E_b1,     cos_b1,     phi_b1,
      E_q2,     cos_q2,     phi_q2,  
      E_qbar2,  cos_qbar2,  phi_qbar2,  
      E_b2,     cos_b2,     phi_b2,
      E_b,      cos_b,      phi_b,  
      E_bbar,   cos_bbar,   phi_bbar;
    void init(){
      E_q1=0.;     cos_q1=0.;     phi_q1=0.;  
      E_qbar1=0.;  cos_qbar1=0.;  phi_qbar1=0.;  
      E_b1=0.;     cos_b1=0.;     phi_b1=0.;
      E_q2=0.;     cos_q2=0.;     phi_q2=0.;  
      E_qbar2=0.;  cos_qbar2=0.;  phi_qbar2=0.;  
      E_b2=0.;     cos_b2=0.;     phi_b2=0.;
      E_b=0.;      cos_b=0.;      phi_b=0.;  
      E_bbar=0.;   cos_bbar=0.;   phi_bbar=0.;
    }
    PS(){ init(); }
    PS(PS&& a){
      E_q1=a.E_q1;       cos_q1=a.cos_q1;        phi_q1=a.phi_q1;  
      E_qbar1=a.E_qbar1; cos_qbar1=a.cos_qbar1;  phi_qbar1=a.phi_qbar1;  
      E_b1=a.E_b1;       cos_b1=a.cos_b1;        phi_b1=a.phi_b1;
      E_q2=a.E_q2;       cos_q2=a.cos_q2;        phi_q2=a.phi_q2;  
      E_qbar2=a.E_qbar2; cos_qbar2=a.cos_qbar2;  phi_qbar2=a.phi_qbar2;  
      E_b2=a.E_b2;       cos_b2=a.cos_b2;        phi_b2=a.phi_b2;
      E_b=a.E_b;         cos_b=a.cos_b;          phi_b=a.phi_b;  
      E_bbar=a.E_bbar;   cos_bbar=a.cos_bbar;    phi_bbar=a.phi_bbar;     
    }
    void print(ostream& os){
      os << "\tPS[ " << endl;
      os << "\t\t[" << E_q1    << ", " << cos_q1    << ", " << phi_q1     << "]," << endl; 
      os << "\t\t[" << E_qbar1 << ", " << cos_qbar1 << ", " << phi_qbar1  << "]," << endl; 
      os << "\t\t[" << E_b1    << ", " << cos_b1    << ", " << phi_b1     << "]," << endl; 
      os << "\t\t[" << E_q2    << ", " << cos_q2    << ", " << phi_q2     << "]," << endl; 
      os << "\t\t[" << E_qbar2 << ", " << cos_qbar2 << ", " << phi_qbar2  << "]," << endl; 
      os << "\t\t[" << E_b2    << ", " << cos_b2    << ", " << phi_b2     << "]," << endl; 
      os << "\t\t[" << E_b     << ", " << cos_b     << ", " << phi_b      << "]," << endl; 
      os << "\t\t[" << E_bbar  << ", " << cos_bbar  << ", " << phi_bbar   << "]"  << endl; 
      os << "\t]" << endl;
    }
  };

  enum class Hypothesis { TTH=0, TTBB, Undefined };

  enum class FinalState { LH=0, LL, HH, Undefined};

  enum class Assumption { ZeroQuarkLost=0, OneQuarkLost, TwoQuarkLost};

  struct CompPerm {
    bool operator()(int a, int b){
      return a>b;
    }
  };

  
  class Integrand {

  public:
    
    // constructor
    Integrand(Hypothesis =Hypothesis::TTH, int =0);
    
    // detstructor
    ~Integrand();

    // add objects (jets, leptons, MET)
    void push_back_object( const LV& , const Object::Type);

    // initialise (once for event)
    void init();

    // create a map between variable names and positions
    void fill_map( const initializer_list<PSVar>& );

    // main method
    void run();

    // make assumption
    double make_assumption( initializer_list<PSVar>&& );

    // test if given assunption is viable
    bool test_assumption( const size_t&, size_t& );

    // main method. Needed by GSLMCIntegrator
    double Eval(const double*) const;
    
    // create PS point
    PS evaluate_PS(const double*, const vector<int>&) const;

    // solve for Wqbar energy
    double solve_for_W_E_qbar(const PS&, const size_t) const;

    // solve for Tb energy
    double solve_for_T_E_b(const PS&, const size_t) const;

    // solve for Hbbar energy
    double solve_for_H_E_bbar(const PS&, const size_t) const;

    // get lower integration edges
    void get_xL(double*, const initializer_list<PSVar>&, const size_t);

    // get upper integration edges
    void get_xU(double*, const initializer_list<PSVar>&, const size_t);

    // get widths
    double get_width(const double*, const double*, const size_t);

  private:

    // report an error
    int error_code;

    // debug
    int debug_code;

    // integration type
    Hypothesis hypo;

    // number of unknowns
    size_t num_of_vars;

    // measured objects
    vector<MEM::Object> obs_jets;
    vector<MEM::Object> obs_leptons;
    vector<MEM::Object> obs_mets;

    // final state
    FinalState fs;

    // contain indexes of obs_jets that need permutations
    vector< vector<int> > perm_indexes;
    vector< vector<int> > perm_indexes_assumption;

    // map between parameter names (physical) and positions in
    // VEGAS space
    map< PSVar, size_t > map_to_var;

    // VEGAS integrator
    ROOT::Math::GSLMCIntegrator* ig2;

  };

}



#endif
