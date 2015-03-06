#ifndef HYPOTESTER_H
#define HYPOTESTER_H

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "TTH/MEIntegratorStandalone/interface/StandardIncludes.h"
#include "TTH/MEIntegratorStandalone/interface/Event.h"

namespace Algo {


  class HypoTester {
    
  public:
   
    // constructor with output file
    HypoTester(TTree* =nullptr);
    
    // destructor
    ~HypoTester();

    // add objects (jets, leptons, MET)
    void push_back_object       ( const LV& ,    char);

    // add objects observables
    // they will be accessible in the form of a map
    void add_object_observables ( const string&, const double& , const char);
     
    // method called by user
    void test( const map< string,vector<Decay::Decay> >& );

    // print objects
    void print(ostream&);

    // set verbosity
    void set_verbosity(const int&);

    // validity
    int get_status();
    
  private:

    // Add hypothesis. Called by test()
    void assume(Decay::Decay); 

    // Read the parameters and do permutations. Uses:
    //   - unpack_assumptions()
    //   - group_particles()
    void init(); 

    // unpack decays providing list of particles
    void unpack_assumptions();

    // assign jet to particles
    vector<Algo::DecayBuilder*> group_particles();

    // clean content, but not input parameters 
    void reset();

    // clean content, including input parameters
    // prepare for following event
    void next_event();

    // call to Minui2 minimize method. Uses:
    //  - setup_minimizer() 
    void run();

    // eval method (called by Minuit)
    double eval(const double* );

    // set minimizer options. Uses:
    //  - is_variable_used()
    void setup_minimizer( const Algo::Strategy );    

    // check whether a given variable is 
    // in at least one permutation
    bool is_variable_used( const std::size_t );

    // print permutation  
    void print_permutation(const vector<std::pair<FinalState::FinalState,std::size_t> >&) const;

    // filter out permutations. Uses:
    //  - print_permutation()
    bool go_to_next(vector<vector<std::pair<FinalState::FinalState,std::size_t> > > &);

    // internal method to save global variables
    void save_global_variables();

    // the minimizer
    ROOT::Math::Minimizer* minimizer;     

    // vector of jets, leptons, and MET
    vector<Algo::Object> p4_Jet;  
    vector<Algo::Object> p4_Lepton;  
    vector<Algo::Object> p4_MET;  

    // filled by test()
    vector<Decay::Decay> decays;

    // vector or particles. Filled by unpack_assumptions()
    // it gets permutated by init()
    vector<pair<FinalState::FinalState, std::size_t> > particles;

    // the actual permutations
    vector<Algo::CombBuilder*> permutations;

    std::size_t nParam_j;
    std::size_t nParam_n;
    std::size_t nParam_m;
    int    count_hypo;
    int    count_perm;
    std::size_t count_TopHad;
    std::size_t count_TopHadLost;
    std::size_t count_WHad;
    std::size_t count_TopLep;
    std::size_t count_Higgs;
    std::size_t count_Radiation_u;
    std::size_t count_Radiation_d;
    std::size_t count_Radiation_c;
    std::size_t count_Radiation_b;
    std::size_t count_Radiation_g;
    std::size_t invisible;
    int verbose;
    int error_code;

    Event* event;

    CompFinalState MyComp;

  };


}

#endif
