#ifndef BTAGRANDOMIZER_H
#define BTAGRANDOMIZER_H

#include "interface/Parameters.h"
#include "interface/Utils.h"
#include <random>

namespace MEM {

  const std::map<DistributionType::DistributionType, TH3D> CONSTMAP = std::map<DistributionType::DistributionType, TH3D>();

  struct Lock{
    Lock(const std::size_t& s){
      size = s;
      lock    = new vector<bool>(s);
      haslock = new vector<bool>(s);
      for(size_t l = 0 ; l < s; ++l) (*lock)[l]    = false;
      for(size_t l = 0 ; l < s; ++l) (*haslock)[l] = false;
    }
    ~Lock(){
      lock->clear();
      haslock->clear();
      delete lock;
      delete haslock;
    }
    void set_lock(const size_t& j, const bool& val){
      (*lock)[j] = val;
    }
    void set_haslock(const size_t& j, const bool& val){
      (*haslock)[j] = val;
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
    void reset(){
      for(size_t l = 0 ; l < size; ++l) (*lock)[l] = false;
    }
    vector<bool>* lock;
    vector<bool>* haslock;
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
      n_tags_l = 0;
      n_tags_h = 0;
    }
    ~BTagRandomizerOutput(){
      input_btag.clear();
      rnd_btag.clear();
    }
    double p;
    int n_tags_l;
    int n_tags_h;
    int ntoys;
    int pass_rnd;
    int pass;
    int err;
    int n_b;
    int n_c;
    int n_l;
    vector<double> input_btag; 
    vector<double> rnd_btag; 
    void print(ostream& os){
      os << "\t************ BTag output ************" << endl;
      os << "\tPassing probability: " << p << endl;
      os << "\tNumber of tags:      " << n_tags_l << ", " << n_tags_h << endl;
      os << "\tPass:                " << pass << endl;
      os << "\tPass (random):       " << pass_rnd << endl;
      os << "\tError code:          " << err << endl;
      os << "\tB,C,L:               " << n_b << "," << n_c << "," << n_l << endl;
      os << "\tNumber of toys:      " << ntoys << endl;
      os << "\tInput:    --------- Random:" << endl;
      for(size_t i = 0 ; i < rnd_btag.size() ; ++i)
	printf( "\t[ %.3f ] --------- [ %.3f ]\n", input_btag[i], rnd_btag[i]);
      os << "\t*************************************" << endl;
    }
  };

  class BTagRandomizer{
    
  public:
    BTagRandomizer(int =0, int =0, const std::map<DistributionType::DistributionType, TH3D>& =CONSTMAP, int =0, int =5000);
    ~BTagRandomizer();
    BTagRandomizerOutput run();
    void push_back_object( Object* );
    void next_event();
    void set_condition(const int&, const int&, const double&);

  private:

    vector<int> get_permutation(const std::size_t&,const int&);
    void init_pdfs();

    TRandom3* ran;
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
