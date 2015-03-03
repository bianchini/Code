#ifndef INTEGRAND_H
#define INTEGRAND_H

// user headers
#include "interface/Parameters.h"


extern "C" {
  void pphttxcallme2born_( double*, double[20], double*, double* );
}
extern "C" {
  void ppttxbbxcallme2born_( double*, double[24], double*, double* );
}


namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF (int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

namespace MEM {
  
  class Integrand {

  public:
    
    // constructor
    Integrand(Hypothesis =Hypothesis::TTH, int =0);
    
    // detstructor
    ~Integrand();

    // add objects (jets, leptons, MET)
    void push_back_object( const LV& , const ObjectType&);

    // initialise (once for event)
    void init();

    // create a map between variable names and positions
    void fill_map( const initializer_list<PSVar>& );

    // main method
    void run();

    // make assumption
    double make_assumption( initializer_list<PSVar>&& );

    // clear after each event
    void clear();

    // test if given assunption is viable
    bool test_assumption( const size_t&, size_t& );

    // main method. Needed by GSLMCIntegrator
    double Eval(const double*) const;
    
    // create PS point
    void create_PS   (MEM::PS&, const double*, const vector<int>&) const;
    void create_PS_LH(MEM::PS&, const double*, const vector<int>&) const;
    void create_PS_LL(MEM::PS&, const double*, const vector<int>&) const;
    void create_PS_HH(MEM::PS&, const double*, const vector<int>&) const;

    void extend_PS(PS&, const PSPart&, const double& , const double& , const TVector3&, const int&, const PSVar&,const PSVar&,const PSVar&,const TFType&) const;

    double probability(const double*, const vector<int>&) const;

    double t_decay_amplitude(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&) const;

    double H_decay_amplitude(const TLorentzVector&, const TLorentzVector&) const;

    double pdf(const double&, const double&, const double&) const;

    double scattering(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, 
		      double&, double&) const;

    // solve for energy given the masses and angles
    double solve( const LV&,  const double& ,  const double& , const TVector3&, const double&) const;

    // get integration edges
    void get_edges(double*, const initializer_list<PSVar>&, const size_t&, const size_t&);

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

    size_t naive_jet_counting;

    // measured objects
    vector<MEM::Object*> obs_jets;
    vector<MEM::Object*> obs_leptons;
    vector<MEM::Object*> obs_mets;

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
