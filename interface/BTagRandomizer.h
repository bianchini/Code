#ifndef BTAGRANDOMIZER_H
#define BTAGRANDOMIZER_H

#include "TTH/MEIntegratorStandalone/interface/Parameters.h"
#include "TTH/MEIntegratorStandalone/interface/Utils.h"
#include <random>

namespace MEM {

  const std::map<DistributionType::DistributionType, TH3D> CONSTMAP = std::map<DistributionType::DistributionType, TH3D>();

  struct Lock{
    Lock(const std::size_t& s){
      size = s;
      lock       = new vector<bool>(s);
      haslock    = new vector<bool>(s);
      alwayslock = new vector<bool>(s);
      for(size_t l = 0 ; l < s; ++l) (*lock)[l]       = false;
      for(size_t l = 0 ; l < s; ++l) (*haslock)[l]    = false;
      for(size_t l = 0 ; l < s; ++l) (*alwayslock)[l] = false;
    }
    ~Lock(){
      lock->clear();
      haslock->clear();
      alwayslock->clear();
      delete lock;
      delete haslock;
      delete alwayslock;
    }
    void set_lock(const size_t& j, const bool& val){
      (*lock)[j] = val;
    }
    void set_haslock(const size_t& j, const bool& val){
      (*haslock)[j] = val;
    }
    void set_alwayslock(const size_t& j, const bool& val){
      (*alwayslock)[j] = val;
    }
    bool allset(){
      bool out = true;
      for(size_t l = 0 ; l < size; ++l){
	if((*haslock)[l]) out &= (*lock)[l];
      }
      return out;
    }
    bool get_lock(const size_t& j){
      return (*lock)[j];
    }
    bool get_haslock(const size_t& j){
      return (*haslock)[j];
    }
    bool get_alwayslock(const size_t& j){
      return (*alwayslock)[j];
    }
    void reset(){
      for(size_t l = 0 ; l < size; ++l){
	(*lock)[l]       = false;
	(*haslock)[l]    = false;
	(*alwayslock)[l] = false;
      }
    }
    void reset_locked(){
      for(size_t l = 0 ; l < size; ++l){
	if( (*haslock)[l] && !(*alwayslock)[l] ) (*lock)[l] = false;
      }
    }
    vector<bool>* lock;
    vector<bool>* haslock;
    vector<bool>* alwayslock;
    std::size_t size;
  };

  struct BTagRandomizerOutput {
    BTagRandomizerOutput(){
      p        = 0.;
      ntoys    = 0;
      err      = 0;
      n_b      = 0;
      n_c      = 0;
      n_l      = 0;
      pass     = 0;
      pass_rnd = 0;
      n_jets   = 0;
      n_tags_l = 0;
      n_tags_h = 0;
      tag_id   = 0;
      tag_name = "";
      seed     = 0;
    }
    ~BTagRandomizerOutput(){
      input_btag.clear();
      rnd_btag.clear();
    }
    double p;
    int n_jets;
    int n_tags_l;
    int n_tags_h;
    int ntoys;
    int pass_rnd;
    int pass;
    int err;
    int n_b;
    int n_c;
    int n_l;
    int tag_id;
    unsigned long seed;
    string tag_name;
    vector<double> input_btag; 
    vector<double> rnd_btag; 
    void print(ostream& os){
      os << "\t************ BTag output ************" << endl;
      os << "\tPassing probability: " << p << endl;
      os << "\tNumber of jets,tags: " << n_jets << ", [" << n_tags_l << ", " << n_tags_h << "]" << endl;
      os << "\tPass:                " << pass << endl;
      os << "\tPass (random):       " << pass_rnd << endl;
      os << "\tError code:          " << err << endl;
      os << "\tB,C,L:               " << n_b << "," << n_c << "," << n_l << endl;
      os << "\tNumber of toys:      " << ntoys << endl;
      os << "\tTag ID:              " << tag_id << endl;
      os << "\tTag name:            " << tag_name << endl;
      os << "\tSeed:                " << seed << endl;
      os << "\tInput:    --------- Random:" << endl;
      for(size_t i = 0 ; i < rnd_btag.size() ; ++i)
	printf( "\t[ %.3f ] --------- [ %.3f ]\n", input_btag[i], rnd_btag[i]);
      os << "\t*************************************" << endl;
    }
  };

  struct JetCategory {
    JetCategory( const int& ntl =0, const int& nth =0, const double& ct =0., const int& t =0, const string& name ="", const unsigned long& sd=4357){
      ntags_l  = ntl;
      ntags_h  = nth;
      cut      = ct;
      tag      = t;
      name_tag = name;
      seed     = sd;
    }
    int tag;
    string name_tag;
    int ntags_l;
    int ntags_h;
    double cut;    
    unsigned long seed;
  };


  class BTagRandomizer{
    
  public:
    BTagRandomizer(int =0, int =0, const std::map<DistributionType::DistributionType, TH3D>& =CONSTMAP, 
		   int =0, 
		   int =1,
		   int =5000);
    ~BTagRandomizer();
    BTagRandomizerOutput run();
    std::vector<BTagRandomizerOutput> run_all( const vector<JetCategory>& );
    void push_back_object( Object* );
    void next_event();
    void next_category();
    void set_condition(const int&, const int&, const double&, const unsigned long&);

  private:

    void fill_perm();
    bool consider(const std::vector<int>&, const std::size_t&) const;
    bool check_ignore_jets ();
    bool check_restore_jets(int&);
    vector<int> get_permutation(const std::size_t&,const int&);
    void init_pdfs();

    TRandom3* ran;
    unsigned long seed;
    bool use_random_seed;
    int init;
    int debug_code;
    int n_tags_l;
    int n_tags_h;
    double cut_val;
    int n_jets;
    int error_code;
    int n_max_toys;
    int count_b;
    int count_c;
    int count_l;
    int assign_rnd;  
    int tag_id;
    int compress_csv;
    string tag_name;
    std::vector<MEM::Object*> jets;
    std::map<std::size_t, TH1D* >  pdfs;
    std::map<std::size_t, double > effs;
    std::map<std::size_t, double > vals;
    std::vector<int> perm_index;
    std::size_t n_perm_max;
    std::map<DistributionType::DistributionType, TH3D> btag_pdfs;
  };

}



#endif
