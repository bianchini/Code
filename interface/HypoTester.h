#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "interface/StandardIncludes.h"


namespace Algo {


  class HypoTester {
    
  public:
   
    struct Object {
      
      LV p4;                   // the four momentum
      map<string,double> obs;  // observables
      
      void init(const LV& p){
	p4 = p;      
      }
      
      void addObs(const string& name, double val){
	obs.insert( make_pair(name, val) );
      }
      
    };

    // compare hypos
    struct CompFinalState {
      bool operator()(pair<FinalState,size_t> a, pair<FinalState,size_t> b){
	return (a.first>b.first) || (a.first==b.first && a.second>b.second);
      }
    } MyComp;


    // constructor
    HypoTester();

    // destructor
    ~HypoTester();

    // add objects (jets, leptons, MET)
    void push_back_object       ( const LV& ,    char);

    // add objects observables
    // they will be accessible in the form of a map
    void add_object_observables ( const string&, const double , char);
     
    // add hypotheses
    void assume(Decay);  

    // read
    void init();
      
    // run
    void run();

    // print objects
    void print(ostream&);

  private:

    // unpack hypotheses
    void unpack_assumptions();

    // group particles
    void group_particles(vector<Algo::DecayBuilder*>& );

    // eval
    double eval(const double* );

    ROOT::Math::Minimizer* minimizer;     

    vector<Object> p4_Jet;  
    vector<Object> p4_Lepton;  
    vector<Object> p4_MET;  

    vector<Decay> decays;
    vector<pair<FinalState,size_t>> particles;

    vector<Algo::CombBuilder*> permutations;

    size_t nParam_j;
    size_t nParam_n;

    int    count_perm;
    size_t count_TopHad;
    size_t count_WHad;
    size_t count_TopLep;
    size_t count_HiggsHad;
    size_t count_Radiation;
    size_t invisible;

    int verbose;
    
  };


}
