#include<iostream>
#include<cstddef>
#include<map> 
#include<string>
#include<algorithm>
#include<assert.h>
#include<memory> 
#include<limits> 


#include "TLorentzVector.h"
#include "TFormula.h"

using namespace std;

typedef TLorentzVector LV;


namespace Algo {

  constexpr double MTOP = 174.3;
  constexpr double MB   = 4.8;
  constexpr double MW   = 80.19;
  constexpr double DM2  = (MTOP*MTOP-MB*MB-MW*MW)*0.5;
  constexpr double MH   = 125.;
  constexpr double DMH2 = (MH*MH-2*MB*MB)*0.5;

  const string TF_Q = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]+[3]/y+[4]/y/y))";
  const string TF_B = "TMath::Gaus(x, [0] + [1]*y, y*TMath::Sqrt([2]+[3]/y+[4]/y/y))";
  const double TF_Q_param[2][5] = 
    {  { 0.0e+00, 1.0e+00, 0.0e+00, 1.5e+00, 0.0e+00 },
       { 0.0e+00, 1.0e+00, 1.3e+01, 1.5e+00, 0.0e+00 } 
    };

  const string TF_MET = "TMath::Gaus(x,0.,20)*TMath::Gaus(y,0.,20)";


  // define it within Algo namespace, but outside class
  enum Decay { TopLep, TopHad, WHad, HiggsHad, Radiation };
  string translateDecay(Decay&);


  class HypoTester {
    
  public:
    
    enum FinalState { TopLep_b=0, 
		      TopHad_q=1,   TopHad_qbar=2,  TopHad_b=3, 
		      WHad_q=4,     WHad_qbar=5,
		      HiggsHad_b=6, HiggsHad_bbar=7, 
		      Radiation_q=8 };

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

    // unpack hypotheses
    void create_tf_met();

    // compare hypos
    struct CompFinalState {
      bool operator()(pair<FinalState,size_t> a, pair<FinalState,size_t> b){
	return (a.first>b.first) || (a.first==b.first && a.second>b.second);
      }
    } MyComp;


  
    class TransferFunction  {
      
    public:
      TransferFunction(const string&, const string&);
      ~TransferFunction();
      void init(const double*);
      const string getFormula() const; 
      double eval (const double& , const double&) const ;
    private:
      string formula;
      TFormula* f;
    };


    class DecayBuilder {
    public:
      virtual ~DecayBuilder() {};
      virtual double eval( const double* , LV&) = 0; 
      virtual void print(ostream&) = 0;     
    };


    class TopHadBuilder: public DecayBuilder {
      
    public:
      TopHadBuilder();
      ~TopHadBuilder();
      void init(const FinalState& , const LV&, const size_t&);
      double eval ( const double* , LV&) ;
      void print(ostream&);     

    private:
      
      LV p4_q;
      LV p4_qbar;
      LV p4_b;
      size_t index_q;
      size_t index_qbar;
      size_t index_b;
      TransferFunction* tf_q;
      TransferFunction* tf_qbar;
      TransferFunction* tf_b;
      Decay decay;
      size_t errFlag;
    };


    class WHadBuilder: public DecayBuilder {
      
    public:
      WHadBuilder();
      ~WHadBuilder();
      void init(const FinalState& , const LV&, const size_t&);
      double eval ( const double* , LV&) ;
      void print(ostream&);     

    private:
      
      LV p4_q;
      LV p4_qbar;
      size_t index_q;
      size_t index_qbar;
      TransferFunction* tf_q;
      TransferFunction* tf_qbar;
      Decay decay;
      size_t errFlag;
    };







    // read
    void read();

    // run
    double run(const double* );

    // group particles
    void group_particles(vector<DecayBuilder*>& );

  private:
     

    vector<Object> p4_Jet;  
    vector<Object> p4_Lepton;  
    vector<Object> p4_MET;  

    vector<Decay> decays;
    vector<pair<FinalState,size_t>> particles;

    int    count_perm;
    size_t count_TopHad;
    size_t count_WHad;
    size_t count_TopLep;
    size_t count_HiggsHad;
    size_t count_Radiation;

    TransferFunction* tf_met;

    int verbose;
    
  };


}
