#ifndef INTEGRAND_H
#define INTEGRAND_H

// user headers
#include "TTH/MEIntegratorStandalone/interface/Parameters.h"

//extern "C" {
//  void pphttxcallme2born_( double*, double[20], double*, double* );
//}
//extern "C" {
//  void ppttxbbxcallme2born_( double*, double[24], double*, double* );
//}
//extern "C" {
//  void ppttxjcallme2born_( double*, double[20], double*, double* );
//}

// Interface
extern "C" {
    void ol_setparameter_int(const char* param, int val);
    void ol_setparameter_double(const char* param, double val);
    void ol_setparameter_string(const char* param, const char* val);
    int ol_register_process(const char* process, int amptype);
    int ol_n_external(int id);
    void ol_phase_space_point(int id, double sqrt_s, double* pp);
    void ol_start();
    void ol_finish();
    void ol_evaluate_tree(int id, double* pp, double* m2_tree);
    void ol_evaluate_loop(int id, double* pp, double* m2_tree, double* m2_loop, double* acc);
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

    /* Public interface */
  public:    

    // constructor (initialise with verbosity)
    Integrand( int =0, const MEMConfig& =MEMConfig() );
    
    // detstructor
    ~Integrand();

    // add objects (jets, leptons, MET)
    // (ObjectType defined in Parameters.h)
    void push_back_object( const LV& , const ObjectType::ObjectType&);    
    void push_back_object( Object* );

    // filter out permutations
    //void set_permutation_strategy(const std::vector<MEM::Permutations>&);
    void set_permutation_strategy(const std::vector<MEM::Permutations::Permutations>&);

    // choose what to include into the integrand
    void set_integrand(const int =0);
    
    // choose c.o.m. energy
    void set_sqrts(const double&);

    // choose n calls
    void set_ncalls(const std::size_t&);

    // set main cfg
    void set_cfg(const MEMConfig&);

    // add object information
    // WARNING: to be called just after push_back of related object!
    // (using .back() method of std::vector<>)
    void add_object_observable(const pair<Observable::Observable, double>&, const ObjectType::ObjectType& );

    // main method: call it to have the ME calculated
    //void run( const FinalState =FinalState::LH, const Hypothesis =Hypothesis::TTH, const std::vector<PSVar> ={});
    MEMOutput run( const FinalState::FinalState=FinalState::LH,
		   const Hypothesis::Hypothesis =Hypothesis::TTH,
		   const std::vector<PSVar::PSVar> = std::vector<PSVar::PSVar>(),
		   const std::vector<PSVar::PSVar> = std::vector<PSVar::PSVar>(),
                   int ncalls=-1
		   );

    // clear containers and counters after each event
    void next_event();
    
    const TF1* get_tf_global(TFType::TFType type, int etabin) const;

    /* Used internally */
  private:

    // initialise (once for event)
    void init( const FinalState::FinalState =FinalState::LH,
        const Hypothesis::Hypothesis=Hypothesis::TTH);

    // create a map between variable names and positions
    //void fill_map( const std::vector<PSVar>& );
    void fill_map( const std::vector<PSVar::PSVar>& );

    // make assumption
    void make_assumption(const std::vector<PSVar::PSVar>&, const std::vector<PSVar::PSVar>&, MEMOutput& );

    void do_integration (const std::size_t&, double*, double*, double&, double&, double&);
    void do_minimization(const std::size_t&, double*, double*, double&, double&, double&);

    // clear containers before new hypothesis
    void next_hypo();

    // test if given assunption is viable
    bool test_assumption( const std::size_t&);

    // filter out permutations
    bool accept_perm( const std::vector<int>&, const std::vector<Permutations::Permutations>& ) const;

    // filter out permutations
    std::vector<std::size_t> get_permutations_byBTAG( const std::size_t& ) const;

    // remove permutations not passing get_permutations_byBTAG()
    void filter_permutations_byBTAG();

    // a constanta value for each permutation
    std::map<MEM::PermConstants::PermConstants, double> get_permutation_constants(  const std::vector<int>& ) const;

    // main method. Needed by GSLMCIntegrator
    double Eval(const double*) const;
    
    // create PS point
    int create_PS    (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_LH (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_LL (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_HH (MEM::PS&, const double*, const std::vector<int>&) const;
    int create_PS_TTH(MEM::PS&, const double*, const std::vector<int>&) const;

    void extend_PS(PS&, const PSPart::PSPart&, const double& , const double& , const TVector3&, const int&, const PSVar::PSVar&,const PSVar::PSVar&,const PSVar::PSVar&,const TFType::TFType&, const int =0) const;
    void extend_PS_nodebug(PS&, const PSPart::PSPart&, const double& , const double& , const TVector3&) const;

    double probability(const double*, const std::size_t&) const;

    double constants() const;

    double t_decay_amplitude(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, const int&) const;

    double H_decay_amplitude(const TLorentzVector&, const TLorentzVector&) const;

    double pdf(const double&, const double&, const double&) const;

    double scattering(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, 
		      double&, double&) const;

    // evaluate TF
    double transfer(const PS&, const vector<int>&, int&) const;

    // evaluate ME
    double matrix  (const PS&) const;

    // evaluate ME
    double matrix_nodecay (const PS&) const;

    // solve for energy given the masses and angles
    double solve( const LV&,  const double& ,  const double& , const TVector3&, const double&, int&) const;

    // get integration edges
    void get_edges(double*, const std::vector<PSVar::PSVar>&, const std::size_t&, const std::size_t&);

    // get widths
    double get_width(const double*, const double*, const std::size_t);

    // get permutation number n
    vector<int> get_permutation(const std::size_t&);

    // setup the minimzer
    void setup_minimizer();

    // improve minimization
    void refine_minimization(std::size_t&, const ROOT::Math::Functor&, const std::size_t&, double*, double*);

    // smear MET
    void smear_met();

    // smear jets
    void smear_jet(MEM::Object*, const bool&);

    // report an error
    int error_code;

    // count function calls
    int n_calls;

    // count maximum function calls
    int n_max_calls;

    // count number of invalid phase space points
    int n_skip;

    // debug
    int debug_code;
    
    // OpenLoops process ID-s
    int proc_id0, proc_id1;

    // main configurator
    MEMConfig cfg;

    // integration type
    Hypothesis::Hypothesis hypo;

    // number of unknowns
    std::size_t num_of_vars;

    // number of original dimensions
    std::size_t ps_dim;

    // keep track of howm many jets one would expect
    std::size_t naive_jet_counting;
    std::size_t extra_jets;

    // measured objects
    std::vector<MEM::Object*> obs_jets;
    std::vector<MEM::Object*> obs_leptons;
    std::vector<MEM::Object*> obs_mets;

    // final state
    FinalState::FinalState fs;

    // the comparateor between permutations
    CompPerm comparator;

    // contain indexes of obs_jets that need permutations
    std::vector<int> perm_index;
    std::size_t n_perm_max;
    std::vector<std::vector<int> > perm_indexes_assumption;
    std::size_t this_perm;
    std::vector< double >          perm_const_assumption;
    std::vector< double >          perm_btag_assumption;
    std::vector< double >          perm_tmpval_assumption;
    std::vector< std::size_t >     perm_pruned;

    std::vector< double >          perm_btag_bb_assumption;
    std::vector< double >          perm_btag_cc_assumption;
    std::vector< double >          perm_btag_jj_assumption;
    
    // status of the prefit (0= not run, 1= run and succesfull, -1= run and unsuccesfull)
    int prefit_code;

    // the btag likelihood (summed over permutations and for the maximal permutation only)
    double btag_weights[3];

    // map between parameter names (physical) and positions in
    // VEGAS space
    boost::unordered_map< PSVar::PSVar, std::size_t, PSVarHash, PSVarEqual> map_to_var;

    // map between a particle and the jet position in obs_jets/obs_leptons
    // to which the particle is matched
    boost::unordered_map< PSPart::PSPart, std::size_t, PSPartHash, PSPartEqual> map_to_part;

    // VEGAS integrator
    ROOT::Math::GSLMCIntegrator* ig2;

    // MINUIT2 minimizer
    ROOT::Math::Minimizer* minimizer;
        
    // an intenal step to perform the pre-fit
    std::size_t prefit_step;
    
    //Stores the global cumulative transfer functions for jet reconstruction efficiency
    const std::map<std::pair<TFType::TFType, int>, TF1*> tf_map;
    
    //the PDF-s for the b-tagger discriminants
    std::map<DistributionType::DistributionType, TH3D*> btag_pdfs;


  };

}

#endif
