#include<iostream>
#include<cstddef>
#include<map> 
#include<string>
#include<algorithm>
#include<assert.h>

#include "TLorentzVector.h"

using namespace std;

typedef TLorentzVector LV;


namespace Algo {

  // define it within Algo namespace, but outside class
  enum Decay { TopLep, TopHad, HiggsHad, Radiation };
  string translateDecay(Decay&);


  class HypoTester {
    
  public:

    enum FinalState { TopLep_b=0, 
		      TopHad_q=1,    TopHad_qbar=2, TopHad_b=3, 
		      HiggsHad_b=4, HiggsHad_bbar=5, 
		      Radiation_q=6 };

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

    // constructor
    HypoTester();

    // destructor
    ~HypoTester();

    // add objects (jets, leptons, MET)
    void push_back_object       ( const LV& ,    char);

    // add objects observables
    // they will be accessible in the form of a map
    void add_object_observables ( const string&, const double , char);
  
    // print objects
    void print(ostream&);

    // add hypotheses
    void assume(Decay);

    // unpack hypotheses
    void unpack_assumptions();

    // compare hypos
    struct CompFinalState {
      bool operator()(pair<FinalState,int> a, pair<FinalState,int> b){
	return (a.first>b.first) || (a.first==b.first && a.second>b.second);
      }
    } MyComp;


    struct TopHad {

      LV p4_q;
      LV p4_qbar;
      LV p4_b;

      void init( const LV&, const LV&, const LV&);
      double get_e_qbar(const double&);
      double get_e_b   (const double&);

    };



    // run
    void run();

    // group particles
    void group_particles();

  private:
      
    vector<Object> p4_Jet;  
    vector<Object> p4_Lepton;  
    vector<Object> p4_MET;  

    vector<Decay> decays;
    vector<pair<FinalState,int>> particles;

    int    count_perm;
    size_t count_TopHad;
    size_t count_TopLep;
    size_t count_HiggsHad;
    size_t count_Radiation;
    
  };


}
