#ifndef HypoTester_h
#define HypoTester_h

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "interface/StandardIncludes.h"
#include "interface/Event.h"

namespace Algo {


  class HypoTester {
    
  public:
   
    // constructor
    HypoTester();
    
    // constructor with output file
    HypoTester(TTree*);
    
    // destructor
    ~HypoTester();

    // add objects (jets, leptons, MET)
    void push_back_object       ( const LV& ,    char);

    // add objects observables
    // they will be accessible in the form of a map
    void add_object_observables ( const string&, const double& , const char);
     
    // add multiple hypotheses
    void test( const map< string,vector<Decay>>& );

    // clean content
    void reset();

    // prepare for next event
    void next_event();

    // add hypotheses
    void assume(Decay);  

    // read
    void init();
      
    // run
    void run();

    // print objects
    void print(ostream&);

    // set verbosity
    void set_verbosity(const int&);

    // set minimizer
    void setup_minimizer( const Algo::Strategy );

  private:

    // unpack hypotheses
    void unpack_assumptions();

    // group particles
    void group_particles(vector<Algo::DecayBuilder*>& );

    // eval
    double eval(const double* );

    // check whether a given variable is in at least one permutation
    bool is_variable_used( const size_t );

    ROOT::Math::Minimizer* minimizer;     

    vector<Algo::Object> p4_Jet;  
    vector<Algo::Object> p4_Lepton;  
    vector<Algo::Object> p4_MET;  

    vector<Decay> decays;
    vector<pair<FinalState,size_t>> particles;

    vector<Algo::CombBuilder*> permutations;

    size_t nParam_j;
    size_t nParam_n;
    size_t nParam_m;

    int    count_hypo;
    int    count_perm;
    size_t count_TopHad;
    size_t count_TopHadLost;
    size_t count_WHad;
    size_t count_TopLep;
    size_t count_Higgs;
    size_t count_Radiation_u;
    size_t count_Radiation_d;
    size_t count_Radiation_b;
    size_t count_Radiation_g;
    size_t invisible;
    int verbose;

    Event* event;

    CompFinalState MyComp;

  };


}

#endif
