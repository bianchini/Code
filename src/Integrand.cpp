#include "interface/Integrand.h"

MEM::Integrand::Integrand(int debug, const MEMConfig& config) :
tf_map(config.tf_map),
btag_pdfs(config.btag_pdfs) {

  // establish invariants
  debug_code         = debug;
  error_code         = 0;
  num_of_vars        = 0;
  ps_dim             = 0;
  naive_jet_counting = 0;
  extra_jets         = 0;
  fs                 = FinalState::Undefined;
  hypo               = Hypothesis::Undefined;
  ig2                = nullptr;
  minimizer          = nullptr;
  prefit_step        = 0;
  n_calls            = 0;
  n_max_calls        = 0;
  n_skip             = 0;
  this_perm          = 0;
  n_perm_max         = 0;
  prefit_code        = 0;
  for(size_t b = 0 ; b < 3; ++b) btag_weights[b] = 0.;
  cfg                = config;
  comparator         = CompPerm();

  // init PDF set
  LHAPDF::initPDFSet(1, cfg.pdfset);

  if( debug_code&DebugVerbosity::init )  
    cout << "Integrand::Integrand(): START" << endl;
}

MEM::Integrand::~Integrand(){
  if( debug_code&DebugVerbosity::init )  
    cout << "Integrand::~Integrand()" << endl;    
  obs_jets.clear();
  obs_leptons.clear();
  obs_mets.clear();
  perm_index.clear();
  for(auto p : perm_indexes_assumption) p.clear();
  perm_indexes_assumption.clear();
  perm_const_assumption.clear();
  perm_btag_assumption.clear();
  perm_btag_bb_assumption.clear();
  perm_btag_jj_assumption.clear();
  perm_btag_cc_assumption.clear();
  perm_tmpval_assumption.clear();
  perm_pruned.clear();
  map_to_var.clear();
  map_to_part.clear();
  if(ig2!=nullptr){
    delete ig2;
    ig2 = nullptr;
  }
  if(minimizer!=nullptr){
    delete minimizer;
    minimizer = nullptr;
  }
}


/* 
   Initialise parameters (***once per event***)
   - determine final state 
   - save jet informations 
   - create list of permutations
   - determine number of variables
*/
void MEM::Integrand::init( const MEM::FinalState::FinalState f, const MEM::Hypothesis::Hypothesis h){

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::init(): START" << endl;
  }    

  // set hypothesis to be tested
  hypo = h;
  fs   = f;

  assert( int(obs_leptons.size())==( 0*(fs==FinalState::HH || fs==FinalState::TTH) 
				     + 1*(fs==FinalState::LH) 
				     + 2*(fs==FinalState::LL) ) );

  // a class member: used to keep track of how many quarks are expected
  naive_jet_counting = 
    8*(fs==FinalState::HH) + 
    6*(fs==FinalState::LH) + 
    4*(fs==FinalState::LL) + 
    0*(fs==FinalState::TTH);

  // deal with jets:
  // if less jets are recorded than the naive_jet_counting,
  // fill in perm_index with -1;  
  size_t n_jets = obs_jets.size();

  // fill the permutation word
  perm_index.clear();
  for(size_t id = 0; id < n_jets ; ++id) perm_index.push_back( id );
  while( perm_index.size() < naive_jet_counting) perm_index.push_back( -1 );

  // smear Sum Nu's by TF
  if(cfg.int_code&IntegrandType::SmearMET)
    smear_met();
  
  // calculate upper / lower edges
  for( auto j : obs_jets ){

    const bool external_tf = (
      cfg.transfer_function_method == TFMethod::External &&
      j->getNumTransferFunctions()>0
    );

    // smear gen jets by TF
    if(cfg.int_code&IntegrandType::SmearJets)
      smear_jet(j, external_tf);    

    Object* obj = nullptr;
    double y[2] = { j->p4().E(), j->p4().Eta() };
    
    //In case using external TF, use Pt (FIXME: unify with internal TF).
    if (external_tf) {
      y[0] = j->p4().Pt();
      obj  = j;
    }

    pair<double, double> edges;    
    if( !j->isSet(Observable::E_LOW_Q) || !j->isSet(Observable::E_HIGH_Q) ){
      edges = get_support( y, TFType::qReco, cfg.j_range_CL,  debug_code, obj) ;
      j->addObs( Observable::E_LOW_Q,  edges.first  );
      j->addObs( Observable::E_HIGH_Q, edges.second );
    }
    if( !j->isSet(Observable::E_LOW_B) || !j->isSet(Observable::E_HIGH_B) ){
      edges = get_support( y, TFType::bReco, cfg.b_range_CL , debug_code, obj) ;
      j->addObs( Observable::E_LOW_B,  edges.first  );
      j->addObs( Observable::E_HIGH_B, edges.second );    
    }
  }

  // sort permutations and compute maximum number of permutations
  comparator = CompPerm(cfg.highpt_first);
  sort( perm_index.begin(), perm_index.end(), comparator );  
  vector<int> perm_index_copy = perm_index;
  n_perm_max = 0;
  do{ ++n_perm_max; } 
  while( next_permutation( perm_index_copy.begin(), perm_index_copy.end(), comparator) );

  if( debug_code&DebugVerbosity::init ){
    cout << "\tMaximum of " << n_perm_max << " permutation(s) considered" << endl;
  }

  // Formula to get the number of unknowns
  // The number of variables is equal to npar - 2*extra_jets
  int unstable = (fs==FinalState::LH || fs==FinalState::LL || fs==FinalState::HH);
  ps_dim = 8*unstable + (3*(hypo==Hypothesis::TTH) + 4*(hypo==Hypothesis::TTBB))*(!unstable);
  
  num_of_vars = 
    // dimension of the phase-space
    3*ps_dim
    // leptons
    -3*obs_leptons.size()               
    // jet directions
    //-2*( TMath::Min(obs_jets.size(), naive_jet_counting) ) 
    -2*naive_jet_counting 
    // top/W mass
    -4*(unstable)               
    // H mass
    -1*(hypo==Hypothesis::TTH && unstable)
    // Px/Py
    -2*(obs_mets.size()==0);

  if( debug_code&DebugVerbosity::init ){
    cout << "\tTotal of " << num_of_vars << " unknowns (does not take into account lost jets)" << endl;
    cout << "\tIntegration code: " << cfg.int_code << endl;
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::init(): END" << endl;
  }
  return;
}

vector<int> MEM::Integrand::get_permutation(const std::size_t& n){
  vector<int> perm_index_copy = perm_index;
  std::size_t n_perm{0};
  do{
    if( n==n_perm ){
      if( debug_code&DebugVerbosity::init_more ) {
	cout << "\tperm. " << n_perm << ": [ ";
	for( auto ind : perm_index_copy ) cout << ind << " ";
	cout << "]" << endl;
      }
      return perm_index_copy;    
    }
    ++n_perm;
  } while( next_permutation( perm_index_copy.begin(), perm_index_copy.end(), comparator) );

  return vector<int>{};
} 


void MEM::Integrand::get_edges(double* lim, const std::vector<PSVar::PSVar>& lost, const size_t& nvar, const size_t& edge){

  // convention is: 
  //   even <=> cosTheta [-1,   +1]
  //   odd  <=> phi      [-PI, +PI]
  size_t count_extra{0};
  pair<double, double> phi_edges = {-TMath::Pi(),+TMath::Pi()};
  double y[2] = { obs_mets[0]->p4().Px(), obs_mets[0]->p4().Py()} ;

  switch( fs ){
  case FinalState::LH:
    if( cfg.m_range_CL<1. )
      phi_edges = get_support( y, TFType::MET, (edge? +cfg.m_range_CL : -cfg.m_range_CL),  debug_code ) ;
    lim[map_to_var[PSVar::E_q1]]      =  edge ?  1. :  0.;
    lim[map_to_var[PSVar::cos_qbar2]] =  edge ? +1  : -1.;
    lim[map_to_var[PSVar::phi_qbar2]] =  edge ? phi_edges.second : phi_edges.first;
    lim[map_to_var[PSVar::E_b]]       =  edge ?  1. :  0.;
    if( hypo==Hypothesis::TTBB ) 
      lim[map_to_var[PSVar::E_bbar]]  =  edge ?  1. :  0.;
    for( auto l = lost.begin() ; l !=lost.end(); ++l ){
      if(edge)
	lim[map_to_var[*l]] = count_extra%2==0 ? +1 : +TMath::Pi(); 
      else
	lim[map_to_var[*l]] = count_extra%2==0 ? -1 : -TMath::Pi(); 
      ++count_extra;
    }          
    break;
  case FinalState::LL:
    lim[map_to_var[PSVar::cos_qbar1]] =  edge ? +1  : -1.;
    lim[map_to_var[PSVar::phi_qbar1]] =  edge ? +TMath::Pi() : -TMath::Pi();
    lim[map_to_var[PSVar::cos_qbar2]] =  edge ? +1  : -1.;
    lim[map_to_var[PSVar::phi_qbar2]] =  edge ? +TMath::Pi() : -TMath::Pi();
    lim[map_to_var[PSVar::E_b]]       =  edge ?  1. :  0.;
    if( hypo==Hypothesis::TTBB ) 
      lim[map_to_var[PSVar::E_bbar]]  =  edge ?  1. :  0.;
    for( auto l = lost.begin() ; l !=lost.end(); ++l ){
      if(edge)
	lim[map_to_var[*l]] = count_extra%2==0 ? +1 : +TMath::Pi(); 
      else
	lim[map_to_var[*l]] = count_extra%2==0 ? -1 : -TMath::Pi(); 
      ++count_extra;
    }          
    break;
  case FinalState::HH:
    lim[map_to_var[PSVar::E_q1]]      =  edge ?  1. :  0.;
    lim[map_to_var[PSVar::E_q2]]      =  edge ?  1. :  0.;
    lim[map_to_var[PSVar::E_b]]       =  edge ?  1. :  0.;
    if( hypo==Hypothesis::TTBB ) 
      lim[map_to_var[PSVar::E_bbar]]  =  edge ?  1. :  0.;
    for( auto l = lost.begin() ; l !=lost.end(); ++l ){
      if(edge)
	lim[map_to_var[*l]] = count_extra%2==0 ? +1 : +TMath::Pi(); 
      else
	lim[map_to_var[*l]] = count_extra%2==0 ? -1 : -TMath::Pi(); 
      ++count_extra;
    }          
    break;
  case FinalState::TTH:
    lim[map_to_var[PSVar::P_t]]      =  edge ?  cfg.emax :  0.;
    lim[map_to_var[PSVar::cos_t]]    =  edge ? +0.99  : -0.99; //0.9 corresponds to |eta|<4.5
    lim[map_to_var[PSVar::phi_t]]    =  edge ? +TMath::Pi() : -TMath::Pi();
    lim[map_to_var[PSVar::P_tbar]]   =  edge ?  cfg.emax :  0.;
    lim[map_to_var[PSVar::cos_tbar]] =  edge ? +0.99  : -0.99;
    lim[map_to_var[PSVar::phi_tbar]] =  edge ? +TMath::Pi() : -TMath::Pi();
    lim[map_to_var[PSVar::Pz_h]]     =  edge ?  cfg.emax/2 : -cfg.emax/2;
    break;
  default:
    break;
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "\tIntegrand::get_edges(): SUMMARY" << endl;
    cout << (edge ? "\t\tH" : "\t\tL") << " edges: [ ";
    for(size_t i = 0 ; i < nvar ; ++i) cout << lim[i] << " " ;
    cout << "]" << endl;
  }
 
}

double MEM::Integrand::get_width(const double* xL, const double* xU, const size_t nvar){
  double out{1.};
  for(size_t i = 0; i < nvar ; ++i) out *= TMath::Abs( xU[i]-xL[i] );
  return out;
}

void MEM::Integrand::fill_map(const std::vector<PSVar::PSVar>& lost){
  
  size_t count_extra{0};
  switch( fs ){
  case FinalState::LH:
    // jet <-> quark
    map_to_part[PSPart::q2]         = 0; // lepton
    map_to_part[PSPart::q1]         = 0; // jet
    map_to_part[PSPart::qbar1]      = 1; // jet
    map_to_part[PSPart::b1]         = 2; // jet
    map_to_part[PSPart::b2]         = 3; // jet
    map_to_part[PSPart::b]          = 4; // jet
    map_to_part[PSPart::bbar]       = 5; // jet

    // PS variable <-> VEGAS variable
    map_to_var[PSVar::E_q1]         = 0; //Eq 
    map_to_var[PSVar::cos_qbar2]    = 1; //cosNu
    map_to_var[PSVar::phi_qbar2]    = 2; //phiNu
    map_to_var[PSVar::E_b]          = 3; //Eb
    if( hypo==Hypothesis::TTBB ) map_to_var[PSVar::E_bbar] = 4;
    count_extra = 4+(hypo==Hypothesis::TTBB);
    for( auto l = lost.begin() ; l!=lost.end() ; ++l ) map_to_var[*l] = count_extra++;
    break;

  case FinalState::LL:
    // jet <-> quark
    map_to_part[PSPart::q1]         = 0; // lepton
    map_to_part[PSPart::q2]         = 1; // lepton                                                                                           
    map_to_part[PSPart::b1]         = 0; // jet
    map_to_part[PSPart::b2]         = 1; // jet
    map_to_part[PSPart::b]          = 2; // jet
    map_to_part[PSPart::bbar]       = 3; // jet

    // PS variable <-> VEGAS variable
    map_to_var[PSVar::cos_qbar1]    = 0; //cosNu
    map_to_var[PSVar::phi_qbar1]    = 1; //phiNu
    map_to_var[PSVar::cos_qbar2]    = 2; //cosNu
    map_to_var[PSVar::phi_qbar2]    = 3; //phiNu
    map_to_var[PSVar::E_b]          = 4; //Eb
    if( hypo==Hypothesis::TTBB ) map_to_var[PSVar::E_bbar] = 5;
    count_extra = 5+(hypo==Hypothesis::TTBB);
    for( auto l = lost.begin() ; l!=lost.end() ; ++l ) map_to_var[*l] = count_extra++;
    break;

  case FinalState::HH:
    // jet <-> quark
    map_to_part[PSPart::q1]         = 0; // jet
    map_to_part[PSPart::qbar1]      = 1; // jet
    map_to_part[PSPart::b1]         = 2; // jet
    map_to_part[PSPart::q2]         = 3; // jet
    map_to_part[PSPart::qbar2]      = 4; // jet
    map_to_part[PSPart::b2]         = 5; // jet
    map_to_part[PSPart::b]          = 6; // jet
    map_to_part[PSPart::bbar]       = 7; // jet

    // PS variable <-> VEGAS variable
    map_to_var[PSVar::E_q1]         = 0; //Eq 
    map_to_var[PSVar::E_q2]         = 1; //Eq
    map_to_var[PSVar::E_b]          = 2; //Eb
    if( hypo==Hypothesis::TTBB ) map_to_var[PSVar::E_bbar] = 3;
    count_extra = 3+(hypo==Hypothesis::TTBB);
    for( auto l = lost.begin() ; l!=lost.end() ; ++l ) map_to_var[*l] = count_extra++;
    break;

  case FinalState::TTH:
    // PS variable <-> VEGAS variable
    map_to_var[PSVar::P_t]          = 0;
    map_to_var[PSVar::cos_t]        = 1;
    map_to_var[PSVar::phi_t]        = 2;
    map_to_var[PSVar::P_tbar]       = 3;
    map_to_var[PSVar::cos_tbar]     = 4;
    map_to_var[PSVar::phi_tbar]     = 5;
    map_to_var[PSVar::Pz_h]         = 6;
  default:
    break;
  }

  if( debug_code&DebugVerbosity::init_more ){
    cout << "\tIntegrand::fill_map(): SUMMARY" << endl;
    cout << "\tMapping between phase-space variables and VEGAS variables:" << endl;
    for( auto iter = map_to_var.begin() ; iter!=map_to_var.end() ; ++iter ){
      cout << "\t\tPS[" << static_cast<int>(iter->first) << "] maps to x[" << iter->second << "]" << endl;
    }
    cout << "\tMapping between quarks and objects in the permuted obs_ collections:" << endl;
    for( auto iter = map_to_part.begin() ; iter!=map_to_part.end() ; ++iter ){
      cout << "\t\tParticle (" << static_cast<int>(iter->first) << ") maps to position " << iter->second << endl;
    }
  }

  return;
}

void MEM::Integrand::push_back_object(const LV& p4,  const MEM::ObjectType::ObjectType& type){

  Object* obj = new Object(p4, type);

  switch( type ){
  case ObjectType::Jet:
    obs_jets.push_back( obj ); 
    break;
  case ObjectType::Lepton:
    obs_leptons.push_back( obj ); 
    break;
  case ObjectType::MET:
    obs_mets.push_back( obj ); 
    break;
  default:
    cout << "*** MEM::Intgrator::push_back_object(): Unknown type of object added" << endl;
    break;
  }

  if( debug_code&DebugVerbosity::init_more ){
    cout << "Integrand::push_back_object(): SUMMARY" << endl;
    obj->print(cout);
  }
  
  return;
}

void MEM::Integrand::push_back_object(MEM::Object* obj){

  switch( obj->type() ){
  case ObjectType::Jet:
    obs_jets.push_back( obj );
    break;
  case ObjectType::Lepton:
    obs_leptons.push_back( obj ); 
    break;
  case ObjectType::MET:
    obs_mets.push_back( obj ); 
    break;
  default:
    cout << "*** MEM::Intgrator::push_back_object(): Unknown type of object added" << endl;
    break;
  }

  if( debug_code&DebugVerbosity::init_more ){
    cout << "Integrand::fill_map(): SUMMARY" << endl;
    obj->print(cout);
  }
  
  return;
}

void MEM::Integrand::add_object_observable( const std::pair<MEM::Observable::Observable, double>& obs, const ObjectType::ObjectType& type ){
  
  switch( type ){
  case ObjectType::Jet :
    if(obs_jets.size()>0)    
      (obs_jets.back())->addObs( obs.first, obs.second );
    break;
  case ObjectType::Lepton :
    if(obs_leptons.size()>0) 
      (obs_leptons.back())->addObs( obs.first, obs.second );
    break;
  case ObjectType::MET:
    if(obs_mets.size()>0) 
      (obs_mets.back())->addObs( obs.first, obs.second );
    break;
  default:
    cout << "Integrand::add_object_observables(): Unknown type of object added" << endl;
    break;
  }

  return;
}

void MEM::Integrand::set_integrand(const int code){
  cfg.int_code = code;
}

void MEM::Integrand::set_ncalls(const size_t& n){
  cfg.n_max_calls = n;
  cfg.is_default  = false;
}

void MEM::Integrand::set_sqrts(const double& s){
  cfg.sqrts = s;
}

void MEM::Integrand::set_cfg(const MEMConfig& config){
  cfg = config;
}

//void MEM::Integrand::set_permutation_strategy(const std::vector<MEM::Permutations>& str){
void MEM::Integrand::set_permutation_strategy(const std::vector<MEM::Permutations::Permutations>& str){
  cfg.perm_pruning = str; 
}


MEM::MEMOutput MEM::Integrand::run( const MEM::FinalState::FinalState f, const MEM::Hypothesis::Hypothesis h, const std::vector<MEM::PSVar::PSVar> missed, const std::vector<MEM::PSVar::PSVar> any){

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::run(): START" << endl;
  }

  // the output for this evaluation
  MEMOutput out;

  // start the clock....
  auto t0 = high_resolution_clock::now();
  
  // prepare permutation, count variables
  init(f,h);

  std::vector<PSVar::PSVar> list;
  for(auto it : missed) list.push_back(it); // quarks out-of-acceptance
  for(auto it : any)    list.push_back(it); // quark directions to be marginalised

  // number of calls
  n_max_calls = cfg.is_default ? 
    cfg.calls[static_cast<std::size_t>(fs)][static_cast<std::size_t>(h)][list.size()/2] : 
    cfg.n_max_calls;

  if( debug_code&DebugVerbosity::init ){
    cout << "\tcfg.calls[" << static_cast<std::size_t>(fs) 
	 << "][" << static_cast<std::size_t>(h) << "][" << list.size()/2 << "] = " 
	 << cfg.calls[static_cast<std::size_t>(fs)][static_cast<std::size_t>(h)][list.size()/2]
	 << " (n_max_calls = " << n_max_calls << ")" << endl;
  }

  // create integrator
  ig2 = new ROOT::Math::GSLMCIntegrator(ROOT::Math::IntegrationMultiDim::kVEGAS, cfg.abs, cfg.rel, n_max_calls);

  if( debug_code&DebugVerbosity::init_more ){
    ig2->Options().Print(std::cout);
    ig2->ExtraOptions()->Print(std::cout);
  }

  // start the clock....
  auto t1 = high_resolution_clock::now();

  // do the calculation
  make_assumption( missed, any, out ); 

  // stop the clock!
  auto t2 = high_resolution_clock::now();

  out.time          = static_cast<int>(duration_cast<milliseconds>(t2-t1).count());
  out.num_perm      = perm_pruned.size()>0 ? perm_pruned.size() : perm_indexes_assumption.size();
  out.final_state   = fs;
  out.hypothesis    = h;
  out.assumption    = list.size()/2;
  out.num_max_calls = n_max_calls;
  out.num_calls     = n_calls;
  out.efficiency    = float(n_calls)/(n_calls+n_skip);
  out.error_code    = error_code;
  out.prefit_code   = prefit_code;
  for(size_t b = 0 ; b < 3; ++b) out.btag_weights[b] = btag_weights[b];

  if( debug_code&DebugVerbosity::output ){
    out.print(cout);
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::run(): DONE in " << static_cast<int>(duration_cast<milliseconds>(t2-t0).count())*0.001 << " sec" << endl;
  }

  // delete stuff and prepare for new hypothesis
  next_hypo();

  return out;
}
void MEM::Integrand::next_event(){
  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::next_event(): START" << endl;
  }
  obs_jets.clear();
  obs_leptons.clear();
  obs_mets.clear();  
  error_code         = 0;
  num_of_vars        = 0;
  ps_dim             = 0;
  naive_jet_counting = 0;
  extra_jets         = 0;
  prefit_step        = 0;
  n_calls            = 0;
  n_skip             = 0;
  cfg.is_default     = true;
  n_perm_max         = 0;
  prefit_code        = 0;
  for(size_t b = 0 ; b < 3; ++b) btag_weights[b] = 0.;
  perm_index.clear();
  for(auto p : perm_indexes_assumption) p.clear();
  perm_indexes_assumption.clear();
  perm_const_assumption.clear();
  perm_btag_assumption.clear();
  perm_btag_bb_assumption.clear();
  perm_btag_jj_assumption.clear();
  perm_btag_cc_assumption.clear();
  perm_tmpval_assumption.clear();
  perm_pruned.clear();
  map_to_var.clear();
  map_to_part.clear();
  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::next_event(): END" << endl;
  }
}

void MEM::Integrand::next_hypo(){
  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::next_hypo(): START" << endl;
  }
  if(ig2!=nullptr) {
    delete ig2;
    ig2 = nullptr;
  }
  if(minimizer!=nullptr) {
    delete minimizer;
    minimizer = nullptr;
  }
  perm_index.clear();
  for(auto p : perm_indexes_assumption) p.clear();
  perm_indexes_assumption.clear();
  perm_const_assumption.clear();
  perm_btag_assumption.clear();
  perm_btag_bb_assumption.clear();
  perm_btag_jj_assumption.clear();
  perm_btag_cc_assumption.clear();
  perm_tmpval_assumption.clear();
  perm_pruned.clear();
  map_to_var.clear();
  map_to_part.clear();
  prefit_step= 0;
  n_calls    = 0;
  n_skip     = 0;
  n_perm_max = 0;
  prefit_code= 0;
  for(size_t b = 0 ; b < 3; ++b) btag_weights[b] = 0.;
  cfg.is_default = true;
  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::next_hypo(): END" << endl;
  }
}

bool MEM::Integrand::test_assumption( const size_t& lost){

  if( (obs_jets.size()+lost) < naive_jet_counting){
    if( debug_code&DebugVerbosity::init ) cout << "\t This assumption cannot be made: too few jets" << endl;
    return false;
  }
  extra_jets = (obs_jets.size()+lost-naive_jet_counting);
  return true;
}

void MEM::Integrand::make_assumption( const std::vector<MEM::PSVar::PSVar>& missed, 
				      const std::vector<MEM::PSVar::PSVar>& any, 
				      MEMOutput& out){

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::make_assumption(): START" << endl;
  }

  double prob{0.};
  double err2{0.};
  double chi2{0.};

  std::vector<PSVar::PSVar> lost;
  for(auto it : missed) lost.push_back( it );
  for(auto it : any)    lost.push_back( it );

  // an assumption may not be consistent with the number of observed jets
  // E.g.: assume 1 lost quark but 2 jets missing wrt to expectation
  // N.B. extra_jets filled here!!!
  if(!test_assumption(lost.size()/2)) return;

  perm_indexes_assumption.clear();

  // extra variables to integrate over  
  fill_map( lost );  

  // Remove unwanted permutations:
  //    CASE (1) ==> perm contains already -1: then -1 must be aligned with the lost quark
  //    CASE (2) ==> perm does not contain -1: then set the correct index to -1
  for( std::size_t n_perm = 0 ; n_perm < n_perm_max ; ++n_perm ){    
    
    auto perm = get_permutation(n_perm);
    if(perm.size()==0) continue;

    // - *it gives the integ. var. position in PSVar
    // - provide first cosTheta: then *it-1 gives the position of E
    // - (*it-1) / 3 gives particle position (0=q1,1=qbar1,2=b1,...)
    size_t count{0};
    for(auto it = missed.begin() ; it!=missed.end() ; ++count, ++it ){
      size_t lost_particle = (static_cast<size_t>(*it)-1)/3;
      if(count%2==0) perm[ map_to_part[ static_cast<PSPart::PSPart>(lost_particle)] ] = -1;      
    }
    for(auto it = any.begin() ; it!=any.end() ; ++count, ++it ){
      size_t lost_particle = (static_cast<size_t>(*it)-1)/3;
      if(count%2==0) perm[ map_to_part[ static_cast<PSPart::PSPart>(lost_particle)] ] = -2;      
    }

    // count the number of lost quarks as assumed in perm
    // if it turns out to be equal to the assumed (lost.size()/2), 
    // then push back the permutation
    count = 0;
    for(auto ind : perm){ if(ind<0) ++count; }

    if(count!=(lost.size()/2))  continue;
    if( !accept_perm( perm, cfg.perm_pruning )) continue;
    
    if( debug_code&DebugVerbosity::init_more ) {
      cout << "\tAdd permutation [ ";
      for( auto ind : perm ) cout << ind << " ";
      cout << "]" << endl;
    }
    perm_indexes_assumption.push_back( perm ); 
    std::map<PermConstants::PermConstants, double> cperm = get_permutation_constants(perm);
    perm_const_assumption.push_back   ( cperm.at(PermConstants::PermConstants::VarTransf) );

    if( hypo==Hypothesis::Hypothesis::TTH || hypo==Hypothesis::Hypothesis::TTBB)
      perm_btag_assumption.push_back    ( cperm.at(PermConstants::PermConstants::btag_TTBB) );
    else
      perm_btag_assumption.push_back(1.0);

    perm_btag_bb_assumption.push_back ( cperm.at(PermConstants::PermConstants::btag_TTBB) );
    perm_btag_jj_assumption.push_back ( cperm.at(PermConstants::PermConstants::btag_TTCC) );
    perm_btag_cc_assumption.push_back ( cperm.at(PermConstants::PermConstants::btag_TTJJ) );
    perm_tmpval_assumption.push_back( 0. );
  }  

  // filter permutations by btag:
  //   - loop over permutations and rank them by decreasing value of btag likelihood
  //   - if Permutations::First*RankedByBTAG is among the perm_pruning strategies, remover underranked permutations
  //   - clean perm_*_assumptions
  //   - if Permutations::First*RankedByBTAG is NOT among the perm_pruning strategies, perm_btag_assumption is filled with 1.0's
  filter_permutations_byBTAG();

  // save the btag weights
  for( auto b_perm : perm_btag_bb_assumption ) btag_weights[static_cast<std::size_t>(PermConstants::PermConstants::btag_TTBB)] += b_perm;
  for( auto b_perm : perm_btag_cc_assumption ) btag_weights[static_cast<std::size_t>(PermConstants::PermConstants::btag_TTCC)] += b_perm;
  for( auto b_perm : perm_btag_jj_assumption ) btag_weights[static_cast<std::size_t>(PermConstants::PermConstants::btag_TTJJ)] += b_perm;

  if( debug_code&DebugVerbosity::init ) {
    size_t n_perm{0};
    for( auto perm : perm_indexes_assumption ){
      cout << "\tperm. " << n_perm << ", k-factor=" << perm_const_assumption[n_perm] 
	   << ", b-tag= " << perm_btag_assumption[n_perm] << ": [ ";
      ++n_perm;
      for( auto ind : perm )
	cout << ind << " ";
      cout << "]" << endl;
    }
  }
  if( debug_code&DebugVerbosity::init ) 
    cout << "\tA total of " << perm_indexes_assumption.size() << " permutations have been considered for this assumption" << endl;

  // create integration ranges
  size_t npar = num_of_vars + lost.size();

  double xL[npar], xU[npar];
  get_edges(xL, lost, npar, 0);
  get_edges(xU, lost, npar, 1);      

  double volume = get_width(xL,xU,npar);

  if(cfg.do_minimize){
    if( debug_code&DebugVerbosity::init ) cout << "\t>>> Minimising..." << endl;
    do_minimization( npar, xL, xU, prob, err2, chi2);
  }

  else if(!cfg.do_minimize && !cfg.do_prefit){
    if( debug_code&DebugVerbosity::init ) cout << "\t>>> Marginalising..." << endl;
    do_integration ( npar, xL, xU, prob, err2, chi2);
  }

  else if(!cfg.do_minimize && cfg.do_prefit){
    if( debug_code&DebugVerbosity::init ) cout << "\t>>> Marginalising with prefit..." << endl;
    prefit_step = 0;
    if( debug_code&DebugVerbosity::init ) cout << "\tSTEP...." << prefit_step << ": global minimisation" << endl;
    do_minimization( npar, xL, xU, prob, err2, chi2);
    prefit_step = 2;
    if( debug_code&DebugVerbosity::init ) cout << "\tSTEP...." << prefit_step << ": marginalisation" << endl;
    n_calls     = 0;
    n_skip      = 0;
    prob        = 0.;
    err2        = 0.;
    error_code  = 0;
    chi2        = 0.;
    //debug_code |= DebugVerbosity::integration;
    do_integration ( npar, xL, xU, prob, err2, chi2);
  }

  else{ /*...*/ }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::make_assumption(): END" << endl;
  }

  if(!cfg.int_code) prob /= (volume*perm_indexes_assumption.size());  

  if( TMath::IsNaN(prob) ){
    cout << "\tdo_integration() returned a NaN" << endl;
    prob = 0.;
    err2 = 0.;
    chi2 = 99.;
    error_code = 1;
    return;
  }

  out.p     = prob;
  out.p_err = sqrt(err2);
  out.chi2  = chi2;

  return;
}

void  MEM::Integrand::do_integration(const std::size_t& npar, double* xL, double* xU, double& prob, double& err2, double& chi2){

  // function
  ROOT::Math::Functor toIntegrate(this, &MEM::Integrand::Eval, npar);  
  ig2->SetFunction(toIntegrate);
  
  // do the integral permutation by permutation
  if( cfg.perm_int ){
    for( std::size_t n_perm = 0; n_perm < perm_indexes_assumption.size() ; ++n_perm ){
      this_perm = n_perm;
      // create integrator
      delete ig2;
      ig2 = new ROOT::Math::GSLMCIntegrator(ROOT::Math::IntegrationMultiDim::kVEGAS, cfg.abs, cfg.rel, n_max_calls);
      ig2->SetFunction(toIntegrate);
      double n_prob = ig2->Integral(xL,xU);
      if( debug_code&DebugVerbosity::output ){
	cout << "\tPermutation num. " << this_perm << " returned p=" << n_prob << endl;
      }
      prob += n_prob;
      err2 += TMath::Power(ig2->Error(),2.);
      chi2 += ig2->ChiSqr()/perm_indexes_assumption.size();
    }
  }

  // do the integral over the sum of permutations
  else{
    double p     = ig2->Integral(xL,xU); 
    double p_err = TMath::Power(ig2->Error(),2.);
    double c2    = ig2->ChiSqr(); 
    if( TMath::IsNaN(p) ){
      cout << "\tIntegral() returned a NaN..." << endl;
      p     = 0.;
      p_err = 0.;
      c2    = 99.;
      error_code = 1;
    }
    prob += p;
    err2 += p_err;
    chi2 += c2;
  }

  return;
}

void  MEM::Integrand::do_minimization(const std::size_t& npar, double* xL, double* xU, double& prob, double& err2, double& chi2){

  ROOT::Math::Functor toIntegrate(this, &MEM::Integrand::Eval, npar);  
  //ig2->SetFunction(toIntegrate);

  // do a global minimization of the integrand  
  if( !cfg.perm_int || (cfg.do_prefit && prefit_step==0)){
    
    size_t num_trials{0};

    // setup the minimzer
    if(minimizer!=nullptr) {
      delete minimizer;
      minimizer = nullptr;
    }
    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    setup_minimizer();
    minimizer->SetFunction(toIntegrate);

    // init variables
    for(size_t np = 0; np < npar ; ++np){
      string var_name = "par_"+std::to_string(np);
      minimizer->SetLimitedVariable(np, var_name.c_str() , (xU[np]+xL[np])/2 , 5e-02 , xL[np] , xU[np] );
      if(debug_code&DebugVerbosity::init ){
	printf("\tParam[%lu] = %s set to %.2f. Range: [%.2f,%.2f]\n", np, var_name.c_str(),  (xU[np]+xL[np])/2 , xL[np], xU[np] ); 
      }
    }

    // run!
    ++num_trials;
    minimizer->Minimize();
    if(cfg.do_prefit && minimizer->Status()!=0) 
      refine_minimization(num_trials, toIntegrate, npar, xL, xU);

    double nll       = minimizer->MinValue();
    const double *xs = minimizer->X(); 

    if(debug_code&DebugVerbosity::init ){
      cout << "\tStatus = " << minimizer->Status()  << " after " << num_trials << " trials"
	   << ", Minimum nll = " <<  nll 
	   << " (p = " << TMath::Exp(-nll)  << ")" << endl;
      for( size_t var = 0 ; var < npar ; ++var)
	cout << "\tVar[" << var << "] = " << xs[var] << endl;
    }
    
    if(cfg.do_prefit && minimizer->Status()==0){
      prefit_code = 1;
      prefit_step = 1;
      if( debug_code&DebugVerbosity::init ) cout << "\tSTEP...." << prefit_step << ": evaluate permutations at minimum" << endl;

      for( std::size_t n_perm = 0; n_perm < perm_indexes_assumption.size() ; ++n_perm ){
	this_perm   = n_perm;
	perm_tmpval_assumption[n_perm] = TMath::Exp(-Eval( xs ));
      }
      perm_pruned = get_sorted_indexes(perm_tmpval_assumption, cfg.perm_filtering_rel);
      
      if( debug_code&DebugVerbosity::init ){
	cout << "\tPruning the " <<  perm_indexes_assumption.size() << " permutations using the pre-fit" << endl;
	for( size_t it = 0; it <  perm_tmpval_assumption.size() ; ++it ){
	  if(is_in(perm_pruned,it)) 
	    cout << "\t\t" << perm_tmpval_assumption[it] << " <==" << endl;    
	  else
	    cout << "\t\t" << perm_tmpval_assumption[it] << endl;    
	}
	cout << "\tTotal of " << perm_pruned.size() << " permutations filtered." << endl;
      }

    }
    else{
      prefit_code = -1;
    }

    // fill variables
    prob      += nll;
    error_code = minimizer->Status();
  }

  // do a permutation-by-permutation minimization of the integrand
  else{

    double nll_min{ numeric_limits<double>::max() };
    for( std::size_t n_perm = 0; n_perm < perm_indexes_assumption.size() ; ++n_perm ){
      this_perm = n_perm;

      // setup the minimzer
      if(n_perm!=0) delete minimizer;
      minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
      setup_minimizer();

      // init variables
      for(size_t np = 0; np < npar ; ++np){
	string var_name = "par_"+std::to_string(np);
	minimizer->SetLimitedVariable(np, var_name.c_str() , (xU[np]+xL[np])/2 , 5e-02 , xL[np] , xU[np] );
	if(debug_code&DebugVerbosity::init ){
	  printf("\tParam[%lu] = %s set to %.2f. Range: [%.2f,%.2f]\n", np, var_name.c_str(),  (xU[np]+xL[np])/2 , xL[np], xU[np] ); 
	}
      }

      // run!
      minimizer->SetFunction(toIntegrate);
      minimizer->Minimize();
      double nll = minimizer->MinValue();
      if(debug_code&DebugVerbosity::init ){
	const double *xs = minimizer->X();
	cout << "\tPermutation num. " << this_perm << " returned "
	     << "Status = " << minimizer->Status() 
	     << ", Minimum nll = " <<  nll 
	     << " (p = " << TMath::Exp(-nll)  << ")" << endl;
	if(debug_code&DebugVerbosity::init_more ){
	  for( size_t var = 0 ; var < npar ; ++var)
	    cout << "\tVar[" << var << "] = " << xs[var] << endl;
	}
      }
      
      error_code += minimizer->Status();
      if(nll < nll_min) prob = nll;
    }
  }
  
  return;
}


// N.B. like this, the code will loop over all jets and choose the combination(s) with largest likelihood.
// ==> always 12*rank_max combinations
// if uncommented, one gets 12*rank_max*(# of combinations of jets.size() taken [...] at the time), where [...] depends on the hypothesis
// ==> this may be VERY CPU intensive
std::vector<std::size_t> MEM::Integrand::get_permutations_byBTAG( const std::size_t& rank_max ) const {

  if( debug_code&DebugVerbosity::init_more ){
    cout << "Integrand::get_permutations_byBTAG(): START" << endl;
  }


  /*
  // these are the indexes that identify different choice of jets 
  vector<std::size_t> indexes_to_check;
  if( fs==FinalState::LH ) indexes_to_check = vector<size_t> {map_to_part.find(PSPart::q1)->second, 
							      map_to_part.find(PSPart::qbar1)->second,
							      map_to_part.find(PSPart::b1)->second, 
							      map_to_part.find(PSPart::b2)->second, 
							      map_to_part.find(PSPart::b)->second, 
							      map_to_part.find(PSPart::bbar)->second};
  if( fs==FinalState::LL ) indexes_to_check = vector<size_t> {map_to_part.find(PSPart::b1)->second, 
							      map_to_part.find(PSPart::b2)->second, 
							      map_to_part.find(PSPart::b)->second, 
							      map_to_part.find(PSPart::bbar)->second};
  if( fs==FinalState::HH ) indexes_to_check  = vector<size_t>{map_to_part.find(PSPart::q1)->second,
							      map_to_part.find(PSPart::qbar1)->second,
							      map_to_part.find(PSPart::b1)->second, 
							      map_to_part.find(PSPart::q2)->second,
							      map_to_part.find(PSPart::qbar2)->second,
							      map_to_part.find(PSPart::b2)->second,
							      map_to_part.find(PSPart::b)->second, 
							      map_to_part.find(PSPart::bbar)->second};
  */

  // store here the indices to perm_indexes_assumption that pass the selection
  vector<std::size_t> pass;

  for( size_t p = 0 ; p < perm_indexes_assumption.size() ; ++p){
    const std::vector<int> perm_p = perm_indexes_assumption[p];
    
    double btag_p = perm_btag_assumption[p];
    double tmpMax = btag_p;

    size_t rank{0};

    if( debug_code&DebugVerbosity::init_more ) 
      cout << "\tPermutation " << p << " has btag likelihood = " << btag_p << ", tmp maximum at " << tmpMax << endl; 
  
    for( size_t q = 0 ; q < perm_indexes_assumption.size() && rank < rank_max; ++q){
      if( q==p ) continue;
      
      /*
      const std::vector<int> perm_q = perm_indexes_assumption[q];
      bool hasSameJets{true};
      for( auto ind_p : indexes_to_check){
	int indx_p = perm_p[ind_p];
	if(indx_p<0) continue;
	bool matched{false};
	for( auto ind_q : indexes_to_check){
	  int indx_q = perm_q[ind_q];
	  if(indx_q==indx_p) matched = true;       
	}
	hasSameJets &= matched;
      }
      if(!hasSameJets) continue;
      */

      double btag_q = perm_btag_assumption[q];
      
      if( (btag_q-tmpMax)/tmpMax > 1e-03 ){
	tmpMax          = btag_q;
	++rank;
	if( debug_code&DebugVerbosity::init_more ) 
	  cout << "\t\t\tnew maximum at " << tmpMax << ", tmp ranking of p is " << rank << endl;
      }      
    }
    
    if(rank<rank_max){
      pass.push_back(p);
      if( debug_code&DebugVerbosity::init_more ) cout << "\tPermutation " << p << " has been accepted!" << endl; 
    }
  }  

  if( debug_code&DebugVerbosity::init_more ){
    cout << "Integrand::get_permutations_byBTAG(): END" << endl;
  }

  return pass;
} 

void MEM::Integrand::filter_permutations_byBTAG(){

  if( debug_code&DebugVerbosity::init){
    cout << "Integrand::filter_permutations_byBTAG(): START" << endl;
  }

  // the indices to perm_indexes_assumption that correspond to the selected permutations
  std::vector<std::size_t> pruned_by_btag;

  // find out which stragey was requested (3 options for the moment)
  bool take_first_1_by_btag{false};
  bool take_first_2_by_btag{false};
  bool take_first_3_by_btag{false};

  for( auto strat : cfg.perm_pruning ){
    if(strat == Permutations::Permutations::FirstRankedByBTAG)      take_first_1_by_btag = true;
    if(strat == Permutations::Permutations::FirstTwoRankedByBTAG)   take_first_2_by_btag = true;
    if(strat == Permutations::Permutations::FirstThreeRankedByBTAG) take_first_3_by_btag = true;
  }

  if( take_first_1_by_btag || take_first_2_by_btag || take_first_3_by_btag ){

    if( take_first_1_by_btag && 
	!take_first_2_by_btag) pruned_by_btag = get_permutations_byBTAG(1); 
    if( take_first_2_by_btag &&
	!take_first_3_by_btag) pruned_by_btag = get_permutations_byBTAG(2); 
    if( take_first_3_by_btag ) pruned_by_btag = get_permutations_byBTAG(3); 

    if( debug_code&DebugVerbosity::init) 
      cout << "\tPruning by btag selected " << pruned_by_btag.size() << " permutations." << endl;
    
    // temporary containers for the permutations
    std::vector<std::vector<int> > perm_indexes_assumption_tmp;
    std::vector< double >          perm_const_assumption_tmp;
    std::vector< double >          perm_btag_assumption_tmp;
    std::vector< double >          perm_btag_bb_assumption_tmp;
    std::vector< double >          perm_btag_cc_assumption_tmp;
    std::vector< double >          perm_btag_jj_assumption_tmp;
    std::vector< double >          perm_tmpval_assumption_tmp;

    // store filtered permutations
    for(auto p : pruned_by_btag ){
      if( debug_code&DebugVerbosity::init ) {
	cout << "\tAccept permutation [ ";
	for( auto ind : perm_indexes_assumption[p] ) cout << ind << " ";
	cout << "]" << endl;
      }
      perm_indexes_assumption_tmp.push_back( perm_indexes_assumption[p] );
      perm_const_assumption_tmp.push_back  ( perm_const_assumption[p]   );
      perm_btag_assumption_tmp.push_back   ( perm_btag_assumption[p]    );
      perm_btag_bb_assumption_tmp.push_back( perm_btag_bb_assumption[p] );
      perm_btag_cc_assumption_tmp.push_back( perm_btag_cc_assumption[p] );
      perm_btag_jj_assumption_tmp.push_back( perm_btag_jj_assumption[p] );
      perm_tmpval_assumption_tmp.push_back ( perm_tmpval_assumption[p]  );
    }    

    // clear old vectors to store the new ones
    for(auto p : perm_indexes_assumption) p.clear();
    perm_indexes_assumption.clear();
    perm_const_assumption.clear();
    perm_btag_assumption.clear();
    perm_btag_bb_assumption.clear();
    perm_btag_jj_assumption.clear();
    perm_btag_cc_assumption.clear();    
    perm_tmpval_assumption.clear();

    // better to use a std::move here to speed up...
    perm_indexes_assumption   = std::move(perm_indexes_assumption_tmp);
    perm_const_assumption     = std::move(perm_const_assumption_tmp);
    perm_btag_assumption      = std::move(perm_btag_assumption_tmp);
    perm_btag_bb_assumption   = std::move(perm_btag_bb_assumption_tmp);
    perm_btag_cc_assumption   = std::move(perm_btag_cc_assumption_tmp);
    perm_btag_jj_assumption   = std::move(perm_btag_jj_assumption_tmp);
    perm_tmpval_assumption    = std::move(perm_tmpval_assumption_tmp);
  }

  // if no pruning based on the b likelihood is requested, then do not use b tagging for probability()
  else{
    for( auto& p : perm_btag_assumption ) p = 1.0; 
  }

  if( debug_code&DebugVerbosity::init){
    cout << "Integrand::filter_permutations_byBTAG(): END" << endl;
  }

  return;
}


bool MEM::Integrand::accept_perm( const vector<int>& perm, const std::vector<MEM::Permutations::Permutations>& strategies ) const {

  if( debug_code&DebugVerbosity::init_more ){
    cout << "Integrand::accept_perm(): START" << endl;
  }

  // helper containers
  vector<size_t> indexes1;
  vector<size_t> indexes2;

  // loop over strategies bto filter out permutations
  for( auto s = strategies.begin() ; s != strategies.end() ; ++s){

    if( debug_code&DebugVerbosity::init_more ){
      cout << "\tUsing strategy Permutations::" << static_cast<size_t>(*s) << "..." << endl;
    }

    switch( *s ){

      // Require all b quarks to be matched to tagged jets
    case Permutations::BTagged:
      indexes1 =  vector<size_t>{map_to_part.find(PSPart::b1)->second, 
				 map_to_part.find(PSPart::b2)->second, 
				 map_to_part.find(PSPart::b)->second, 
				 map_to_part.find(PSPart::bbar)->second};
      for( auto ind : indexes1 ){ 
	if( perm[ind]>=0 && 
	    obs_jets[perm[ind]]->isSet(Observable::BTAG) && 
	    obs_jets[perm[ind]]->getObs(Observable::BTAG)<0.5){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation: obs_jets[ perm[" << ind << "] ] has BTAG=" << obs_jets[perm[ind]]->getObs(Observable::BTAG) << endl;
	  }
	  return false;
	}
      }
      break;

      // Require all non-b quarks to be matched to untagged jets
    case Permutations::QUntagged:
      indexes1.clear();
      if( fs==FinalState::LH ) indexes1 = vector<size_t>{map_to_part.find(PSPart::q1)->second, 
							 map_to_part.find(PSPart::qbar1)->second};
      if( fs==FinalState::LL ) indexes1 = vector<size_t>{};
      if( fs==FinalState::HH ) indexes1  = vector<size_t>{map_to_part.find(PSPart::q1)->second,
							  map_to_part.find(PSPart::qbar1)->second, 
							  map_to_part.find(PSPart::q2)->second,
                                                          map_to_part.find(PSPart::qbar2)->second};
      for( auto ind : indexes1 ){
	if( perm[ind]>=0 &&
            obs_jets[perm[ind]]->isSet(Observable::BTAG) &&
            obs_jets[perm[ind]]->getObs(Observable::BTAG)>0.5){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation: obs_jets[ perm[" << ind << "] ] has BTAG=" << obs_jets[perm[ind]]->getObs(Observable::BTAG) << endl;
	  }
	  return false;  
	}
      }
      break;
      
      // require that no other permutations has been already considered
      // differing from perm by swapping the W quarks
    case Permutations::QQbarSymmetry:
      indexes1.clear();
      indexes2.clear();
      if( fs==FinalState::LH ){
	// the  symmetric part
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::q1)->second, 
				   map_to_part.find(PSPart::qbar1)->second};
	// the asymmetric part
	indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second, 
				   map_to_part.find(PSPart::b)->second,
				   map_to_part.find(PSPart::bbar)->second};  
      }
      if( fs==FinalState::LL ){
	// the  symmetric part                                                                                                                      
        // the asymmetric part                                                                                                                      
        indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second, 
				   map_to_part.find(PSPart::b)->second,
				   map_to_part.find(PSPart::bbar)->second};
      }
      if( fs==FinalState::HH ){
	// the  symmetric part
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::q1)->second, 
				   map_to_part.find(PSPart::qbar1)->second, 
				   map_to_part.find(PSPart::q2)->second,
				   map_to_part.find(PSPart::qbar2)->second};
	// the asymmetric part
	indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second,
				   map_to_part.find(PSPart::b)->second,
				   map_to_part.find(PSPart::bbar)->second};  
      }

      for( auto visited : perm_indexes_assumption ){ 	
	bool asymmetric_part{true};
	bool symmetric_part {true};
	
	// loop over quark positions that should be matched to same jet
	for(auto same : indexes2 ) 
	  if(visited[same]!=perm[same]) asymmetric_part = false;

	// loop over quark position pairs searching for a swap
	for(size_t i = 0 ; i < (indexes1.size()/2) ; ++i){
	  bool same = (visited[ indexes1[2*i] ]==perm[indexes1[2*i]])   && (visited[ indexes1[2*i+1] ]==perm[indexes1[2*i+1]]);
	  bool swap = (visited[ indexes1[2*i] ]==perm[indexes1[2*i+1]]) && (visited[ indexes1[2*i+1] ]==perm[indexes1[2*i]]);
	  if( !(same || swap) ) symmetric_part = false;
	}

	if( asymmetric_part && symmetric_part ){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation: a (Q,QBAR) swap has been found" << endl;
	  }
	  return false;
	}
      }
      break;

    case Permutations::BBbarSymmetry:
      indexes1.clear();
      indexes2.clear();
      if( fs==FinalState::LH ){
	// the  symmetric part
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::b)->second, 
				   map_to_part.find(PSPart::bbar)->second};
	// the asymmetric part
	indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second,
				   map_to_part.find(PSPart::q1)->second, 
				   map_to_part.find(PSPart::qbar1)->second};  
      }
      if( fs==FinalState::LL ){
	// the  symmetric part                                                                                                                      
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::b)->second, 
				   map_to_part.find(PSPart::bbar)->second};
        // the asymmetric part                                                                                                                      
        indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second}; 
      }
      if( fs==FinalState::HH ){
	// the  symmetric part
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::b)->second, 
				   map_to_part.find(PSPart::bbar)->second};
	// the asymmetric part
	indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second, 
				   map_to_part.find(PSPart::q1)->second, 
				   map_to_part.find(PSPart::qbar1)->second};  
      }

      for( auto visited : perm_indexes_assumption ){ 	
	bool asymmetric_part{true};
	bool symmetric_part {true};
	
	// loop over quark positions that should be matched to same jet
	for(auto same : indexes2 ) 
	  if(visited[same]!=perm[same]) asymmetric_part = false;

	// loop over quark position pairs searching for a swap
	for(size_t i = 0 ; i < (indexes1.size()/2) ; ++i){
	  bool same = (visited[ indexes1[2*i] ]==perm[indexes1[2*i]])   && (visited[ indexes1[2*i+1] ]==perm[indexes1[2*i+1]]);
	  bool swap = (visited[ indexes1[2*i] ]==perm[indexes1[2*i+1]]) && (visited[ indexes1[2*i+1] ]==perm[indexes1[2*i]]);
	  if( !(same || swap) ) symmetric_part = false;
	}

	if( asymmetric_part && symmetric_part ){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation: a (B,BBAR) swap has been found" << endl;
	  }
	  return false;
	}
      }
      break;

    case Permutations::QQbarBBbarSymmetry:
      indexes1.clear();
      indexes2.clear();
      if( fs==FinalState::LH ){
	// the  symmetric part
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::q1)->second, 
				   map_to_part.find(PSPart::qbar1)->second,
				   map_to_part.find(PSPart::b)->second,
				   map_to_part.find(PSPart::bbar)->second};
	// the asymmetric part
	indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second};  
      }
      if( fs==FinalState::LL ){
	// the  symmetric part                                                                                                                      
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::b)->second,
				   map_to_part.find(PSPart::bbar)->second};
        // the asymmetric part                                                                                                                      
        indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second};
      }
      if( fs==FinalState::HH ){
	// the  symmetric part
	indexes1 = vector<size_t>{ map_to_part.find(PSPart::q1)->second, 
				   map_to_part.find(PSPart::qbar1)->second, 
				   map_to_part.find(PSPart::q2)->second,
				   map_to_part.find(PSPart::qbar2)->second,
				   map_to_part.find(PSPart::b)->second,
				   map_to_part.find(PSPart::bbar)->second};
	// the asymmetric part
	indexes2 = vector<size_t>{ map_to_part.find(PSPart::b1)->second,
				   map_to_part.find(PSPart::b2)->second};
      }

      for( auto visited : perm_indexes_assumption ){ 	
	bool asymmetric_part{true};
	bool symmetric_part {true};
	
	// loop over quark positions that should be matched to same jet
	for(auto same : indexes2 ) 
	  if(visited[same]!=perm[same]) asymmetric_part = false;

	// loop over quark position pairs searching for a swap
	for(size_t i = 0 ; i < (indexes1.size()/2) ; ++i){
	  bool same = (visited[ indexes1[2*i] ]==perm[indexes1[2*i]])   && (visited[ indexes1[2*i+1] ]==perm[indexes1[2*i+1]]);
	  bool swap = (visited[ indexes1[2*i] ]==perm[indexes1[2*i+1]]) && (visited[ indexes1[2*i+1] ]==perm[indexes1[2*i]]);
	  if( !(same || swap) ) symmetric_part = false;
	}

	if( asymmetric_part && symmetric_part ){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation: a (Q,QBAR) swap has been found" << endl;
	  }
	  return false;
	}
      }
      break;

      // assume one can match quarks to HEPTopTag subjets 
    case Permutations::HEPTopTagged:
      // deal with the first top
      if(!(fs==FinalState::HH || fs==FinalState::LH)) continue;
      indexes1 =  vector<size_t>{map_to_part.find(PSPart::q1)->second, 
				 map_to_part.find(PSPart::qbar1)->second, 
				 map_to_part.find(PSPart::b1)->second};
      for(size_t idx_hep = 0 ; idx_hep<3; ++idx_hep ){
	int idx = indexes1[idx_hep];
	if( perm[idx]>=0 && obs_jets[perm[idx]]->isSet(Observable::PDGID) &&
	    obs_jets[perm[idx]]->getObs(Observable::PDGID)!=((idx_hep<2 ? +1 : +5) ) ){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation because jet[" << perm[idx] 
		 << "] should be associated with a pdgid=" <<(idx_hep<2 ? +1 : +5)  
		 << " quark from top, but its pdgid from HEPTopTagging is " 
		 <<  obs_jets[perm[idx]]->getObs(Observable::PDGID) << endl;
	  }	
	  return false;
	}
      }

      // deal with the second top (if any)
      if(!(fs==FinalState::HH)) continue;
      indexes2 =  vector<size_t>{map_to_part.find(PSPart::q2)->second, 
				 map_to_part.find(PSPart::qbar2)->second, 
				 map_to_part.find(PSPart::b2)->second};
      for(size_t idx_hep = 0 ; idx_hep<3; ++idx_hep ){
	int idx = indexes2[idx_hep];
	if( perm[idx]>=0 && obs_jets[perm[idx]]->isSet(Observable::PDGID) &&
	    obs_jets[perm[idx]]->getObs(Observable::PDGID)!=(idx_hep<2 ? -1 : -5) ){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation because jet[" << perm[idx] 
		 << "] should be associated with a pdgid=" <<  (idx_hep<2 ? -1 : -5) 
		 << " quark from anti-top, but its pdgid from HEPTopTagging is " 
		 <<  obs_jets[perm[idx]]->getObs(Observable::PDGID) << endl;
	  }	
	  return false;
	}
      }
      break;

      // assume one can match quarks to HiggsTagger subjets 
    case Permutations::HiggsTagged:
      indexes1 =  vector<size_t>{map_to_part.find(PSPart::b)->second, 
				 map_to_part.find(PSPart::bbar)->second};
      for(size_t idx_hep = 0 ; idx_hep<2; ++idx_hep ){
	int idx = indexes1[idx_hep];
	if( perm[idx]>=0 && obs_jets[perm[idx]]->isSet(Observable::PDGID) &&
	    obs_jets[perm[idx]]->getObs(Observable::PDGID)!=22 ){
	  if( debug_code&DebugVerbosity::init_more ){
	    cout << "\t\tDiscard permutation because jet[" << perm[idx] 
		 << "] should be associated with a pdgid=22" 
		 << " quark from Higgs, but its pdgid from HiggsTagging is " 
		 <<  obs_jets[perm[idx]]->getObs(Observable::PDGID) << endl;
	  }	
	  return false;
	}
      }
      break;      

    default:
      break;
    }
    
  }
  
  return true;
}

std::map<MEM::PermConstants::PermConstants, double> MEM::Integrand::get_permutation_constants( const vector<int>& perm ) const {

  if( debug_code&DebugVerbosity::init_more ){
    cout << "Integrand::get_permutation_constants(): START" << endl;
  }

  std::map<PermConstants::PermConstants, double> out;

  double p     {1.};
  double DeltaE{1.};
  double bb    {1.};
  double bp    {1.};

  // aditional likelihoods (for TT_JJ/BJ/CC)
  double jj    {1.};
  double cc    {1.};

  size_t pos;
  MEM::Object* obj = nullptr;

  switch( fs ){
  case FinalState::LH:
    // PSVar::E_q1
    pos = map_to_part.find(PSPart::q1)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      DeltaE = (obj->getObs(Observable::E_HIGH_Q) - obj->getObs(Observable::E_LOW_Q));
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp =  0.5*h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) )) 
	  + 0.5*h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", l/c) = " << bp << endl;
      }      
    }
    else{
      DeltaE = cfg.emax - MQ;
      bp     = 1.0;
    }
    if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_q = " << DeltaE << " GeV" << endl;
    p *= DeltaE;
    bb*= bp;

    // PSVar::E_qbar1    
    pos = map_to_part.find(PSPart::qbar1)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", l) = " << bp << endl;
      }      
    }
    else bp = 1.;    
    bb*= bp;

    // PSVar::E_b1    
    pos = map_to_part.find(PSPart::b1)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_b2    
    pos = map_to_part.find(PSPart::b2)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	bp *=  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_b    
    jj = bb;
    cc = bb;
    pos = map_to_part.find(PSPart::b)->second;
    if( perm[ pos ]>=0 ){
      obj = obs_jets[ perm[pos] ];
      bp  = 1.;
      DeltaE = (obj->getObs(Observable::E_HIGH_B) - obj->getObs(Observable::E_LOW_B)); 
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	TH3D h3 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp  = h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	jj *= h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	cc *= h3.GetBinContent( h3.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else{
      DeltaE = cfg.emax - MB;
      bp     = 1.;
    }
    if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_b = " << DeltaE << " GeV" << endl;
    p *= DeltaE;
    bb*= bp;


    // PSVar::E_bbar
    if( hypo==Hypothesis::TTBB ){
      pos = map_to_part.find(PSPart::bbar)->second;
      if( perm[ pos ]>=0 ){
	obj = obs_jets[ perm[pos] ];
	DeltaE = (obj->getObs(Observable::E_HIGH_B) - obj->getObs(Observable::E_LOW_B));
      }
      else{
	DeltaE = cfg.emax - MB;
      }
      if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_bbar = " << DeltaE << " GeV" << endl;
      p *= DeltaE;
    }

    // PSVar::E_bbar 
    pos = map_to_part.find(PSPart::bbar)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	TH3D h3 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp  = h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	jj *= h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	cc *= h3.GetBinContent( h3.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;
    
    break;
  case FinalState::LL:
    // PSVar::E_b1    
    pos = map_to_part.find(PSPart::b1)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_b2    
    pos = map_to_part.find(PSPart::b2)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_b
    jj = bb;
    cc = bb;
    pos = map_to_part.find(PSPart::b)->second;
    if( perm[ pos ]>=0 ){      
      obj = obs_jets[ perm[pos] ];
      bp  = 1.;
      DeltaE = (obj->getObs(Observable::E_HIGH_B) - obj->getObs(Observable::E_LOW_B));
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	TH3D h3 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp  = h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	jj *= h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	cc *= h3.GetBinContent( h3.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else{
      DeltaE = cfg.emax - MB;
      bp     = 1.;
    }
    if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_b = " << DeltaE << " GeV" << endl;
    p *= DeltaE;
    bb*= bp;

    // PSVar::E_bbar 
    if( hypo==Hypothesis::TTBB ){
      pos = map_to_part.find(PSPart::bbar)->second;
      if( perm[ pos ]>=0 ){
	obj = obs_jets[ perm[pos] ];
	DeltaE = (obj->getObs(Observable::E_HIGH_B) - obj->getObs(Observable::E_LOW_B));
      }
      else{
	DeltaE = cfg.emax - MB;
      }
      if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_bbar = " << DeltaE << " GeV" << endl;
      p *= DeltaE;
    }
    // PSVar::E_bbar 
    pos = map_to_part.find(PSPart::bbar)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	TH3D h3 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp  = h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	jj *= h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	cc *= h3.GetBinContent( h3.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    break;

  case FinalState::HH:
    // PSVar::E_q1
    pos = map_to_part.find(PSPart::q1)->second;
    if( perm[ pos ]>=0 ){
      obj = obs_jets[ perm[pos] ];
      bp  = 1.;
      DeltaE = (obj->getObs(Observable::E_HIGH_Q) - obj->getObs(Observable::E_LOW_Q));
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
 	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
 	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp =  h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) )) +
	  h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));;
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", l/c) = " << bp << endl;
      }      
    }
    else{
      DeltaE = cfg.emax - MQ;
      bp     = 1.;
    }
    if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_q = " << DeltaE << " GeV" << endl;
    p *= DeltaE;
    bb*= bp;

    // PSVar::E_qbar1    
    pos = map_to_part.find(PSPart::qbar1)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", l) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_b1    
    pos = map_to_part.find(PSPart::b1)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_q2
    pos = map_to_part.find(PSPart::q2)->second;
    if( perm[ pos ]>=0 ){
      obj = obs_jets[ perm[pos] ];
      bp  = 1.;
      DeltaE = (obj->getObs(Observable::E_HIGH_Q) - obj->getObs(Observable::E_LOW_Q));
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", l) = " << bp << endl;
      }      
    }
    else{
      DeltaE = cfg.emax - MQ;
      bp     = 1.;
    }
    if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_q = " << DeltaE << " GeV" << endl;
    p *= DeltaE;
    bb*= bp;

    // PSVar::E_qbar2    
    pos = map_to_part.find(PSPart::qbar2)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", l) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_b2    
    pos = map_to_part.find(PSPart::b1)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	bp =  h.GetBinContent( h.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    // PSVar::E_b
    jj = bb;
    cc = bb;
    pos = map_to_part.find(PSPart::b)->second;
    if( perm[ pos ]>=0 ){
      obj = obs_jets[ perm[pos] ];
      bp  = 1.;
      DeltaE = (obj->getObs(Observable::E_HIGH_B) - obj->getObs(Observable::E_LOW_B)); 
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	TH3D h3 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp  = h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	jj *= h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	cc *= h3.GetBinContent( h3.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else{
      DeltaE = cfg.emax - MB;
      bp     = 1.;
    }
    if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_b = " << DeltaE << " GeV" << endl;
    p *= DeltaE;
    bb*= bp;

    // PSVar::E_bbar
    if( hypo==Hypothesis::TTBB ){
      pos = map_to_part.find(PSPart::bbar)->second;
      if( perm[ pos ]>=0 ){
	obj = obs_jets[ perm[pos] ];
	DeltaE = (obj->getObs(Observable::E_HIGH_B) - obj->getObs(Observable::E_LOW_B));
      }
      else{
	DeltaE = cfg.emax - MB;
      }
      if( debug_code&DebugVerbosity::init_more ) cout << "\tdE_bbar = " << DeltaE << " GeV" << endl;
      p *= DeltaE;
    }

    // PSVar::E_bbar
    pos = map_to_part.find(PSPart::bbar)->second;
    if( perm[ pos ]>=0 ){
      obj    = obs_jets[ perm[pos] ];
      bp     = 1.;
      if( btag_pdfs.size()>0 && obj->isSet(Observable::CSV) ){
	TH3D h1 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_b);
	TH3D h2 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_l);
	TH3D h3 = btag_pdfs.at(MEM::DistributionType::DistributionType::csv_c);
	bp  = h1.GetBinContent( h1.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	jj *= h2.GetBinContent( h2.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	cc *= h3.GetBinContent( h3.FindBin(obj->p4().Pt(), std::abs(obj->p4().Eta()), obj->getObs(Observable::CSV) ));
	if( debug_code&DebugVerbosity::init_more ) cout << "\tbtag(" <<  obj->getObs(Observable::CSV) << ", " << obj->p4().Pt() << ", " <<  obj->p4().Eta()  << ", b) = " << bp << endl;
      }      
    }
    else bp = 1.;
    bb*= bp;

    break;

  default:
    break;
  }

  if( debug_code&DebugVerbosity::init_more ){
    cout << "\tVariable transformation permutation-dependent returned p*b = " << p << " * " << bb << " = " << p*bb << endl;
  }
  
  if( debug_code&DebugVerbosity::init_more ){
    cout << "Integrand::get_permutation_constants(): END" << endl;
  }

  out[PermConstants::PermConstants::VarTransf] = p;
  out[PermConstants::PermConstants::btag_TTBB] = bb;
  out[PermConstants::PermConstants::btag_TTCC] = cc;
  out[PermConstants::PermConstants::btag_TTJJ] = jj;

  return out;
}


double MEM::Integrand::Eval(const double* x) const{

#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration ){
      cout << "\tIntegrand::Eval(): START" << endl;
      cout << "\t\tFunction call num. " << n_calls << endl;
    }
#endif

  double p{0.};

  // strategy to speed-up 
  if( !cfg.do_minimize && cfg.do_perm_filtering && n_calls==n_max_calls ){
    (const_cast<Integrand*>(this))->perm_pruned = 
      get_sorted_indexes(perm_tmpval_assumption, cfg.perm_filtering_rel);
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::init ){
      cout << "\tPruning the " <<  perm_indexes_assumption.size() << " permutations at the " 
	   << n_calls << "th function call. By decreasing value:" << endl;
      //sort( ((const_cast<Integrand*>(this))->perm_tmpval_assumption).begin(),
      //    ((const_cast<Integrand*>(this))->perm_tmpval_assumption).end(), MEM::descending );
      for( size_t it = 0; it <  perm_tmpval_assumption.size() ; ++it ){
	if(is_in(perm_pruned,it)) 
	  cout << "\t\t" << perm_tmpval_assumption[it] << " <==" << endl;    
	else
	  cout << "\t\t" << perm_tmpval_assumption[it] << endl;    
      }
      cout << "\tTotal of " << perm_pruned.size() << " permutations filtered." << endl;
    }
#endif
  }
  
  for( std::size_t n_perm = 0; n_perm < perm_indexes_assumption.size() ; ++n_perm ){

    if( cfg.perm_int || (cfg.do_prefit && prefit_step==1) ){
      if(n_perm != this_perm) continue;
    }
    if( cfg.do_perm_filtering || (cfg.do_prefit && prefit_step==2) ){
      if( !is_in(perm_pruned,n_perm) ) continue;
    }
    

    double p0 = probability(x, n_perm );
    double p1 = cfg.int_code>0 ? perm_const_assumption[n_perm]*perm_btag_assumption[n_perm] : 1.0;
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration ){
      cout << "\t\tPermutation #" << n_perm 
	   << " => p = (" << p0 << "*" << p1 << ") = " << (p0*p1) << endl;
      cout << "\t\t\tP --> " << p << " + " << (p0*p1) << endl;
    }
#endif
    
    // strategy to speed-up 
    if( cfg.do_perm_filtering && n_calls<n_max_calls)
      ((const_cast<Integrand*>(this))->perm_tmpval_assumption)[ n_perm ] += (p0*p1);
    
    p += (p0*p1);
  }
  
  if( cfg.do_minimize || (cfg.do_prefit && prefit_step<2) ){
    if( p>0. ) p = -TMath::Log(p);
    else p = numeric_limits<double>::max();
  } 

  ++(const_cast<Integrand*>(this)->n_calls);

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::Eval(): END" << endl;
  }
#endif

  if( TMath::IsNaN(p) ){
    cout << "\tEval() returned a NaN..." << endl;
    return 0.;
  }


  return p;
}

int MEM::Integrand::create_PS(MEM::PS& ps, const double* x, const vector<int>& perm ) const {

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS(): START" << endl;
  }
#endif

  switch( fs ){
  case FinalState::LH:
    return create_PS_LH(ps, x, perm);
    break;
  case FinalState::LL:
    return create_PS_LL(ps, x, perm);
    break;
  case FinalState::HH:
    return create_PS_HH(ps, x, perm);
    break;
  case FinalState::TTH:
    return create_PS_TTH(ps, x, perm);
    break;
  default:
    break;
  }

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS(): END" << endl;
  }
#endif

  return 0;
}


int MEM::Integrand::create_PS_LH(MEM::PS& ps, const double* x, const vector<int>& perm ) const {

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_LH(): START" << endl;
  }
#endif

  // corrupted phase space
  int accept{0};

  // store temporary values to build four-vectors
  double E     {0.};
  double E_LOW {0.};
  double E_HIGH{0.};
  double E_REC {numeric_limits<double>::max()};
  TVector3 dir (1.,0.,0.);
  TFType::TFType tftype = TFType::Unknown;

  // map a quark to an index inside the obs_ collections
  size_t nj_q1    =  map_to_part.find(PSPart::q1)->second; 
  size_t nj_qbar1 =  map_to_part.find(PSPart::qbar1)->second; 
  size_t nj_b1    =  map_to_part.find(PSPart::b1)->second; 
  size_t nl_q2    =  map_to_part.find(PSPart::q2)->second; 
  size_t nj_b2    =  map_to_part.find(PSPart::b2)->second; 
  size_t nj_b     =  map_to_part.find(PSPart::b)->second;
  size_t nj_bbar  =  map_to_part.find(PSPart::bbar)->second;

  /////  PSPart::q1
  if( perm[ nj_q1 ]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_q1] ]; 
    dir     = obj->p4().Vect().Unit();
    E_LOW   = obj->getObs(Observable::E_LOW_Q) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_Q);
    tftype  = TFType::qReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_q1)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_q1)->second ] );
    E_LOW   = MQ;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) ); // Restrict to 2*E threshol
    tftype  = (perm[nj_q1]==-1 ? TFType::qLost : TFType::Unknown);
  }
  E       = E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_q1)->second ]);
  extend_PS( ps, PSPart::q1, E, MQ, dir, perm[nj_q1], PSVar::cos_q1, PSVar::phi_q1, PSVar::E_q1, tftype ); 

  /////  PSPart::qbar1
  if( perm[nj_qbar1]>=0 ){
    dir     = obs_jets[ perm[nj_qbar1] ]->p4().Vect().Unit();
    tftype  = TFType::qReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_qbar1)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_qbar1)->second ] );
    tftype  = (perm[nj_qbar1]==-1 ? TFType::qLost : TFType::Unknown);
  }
  E    = solve( ps.lv(PSPart::q1), DMW2 , MQ, dir, E_REC, accept );
  extend_PS( ps, PSPart::qbar1, E, MQ, dir, perm[nj_qbar1], PSVar::cos_qbar1, PSVar::phi_qbar1, PSVar::E_qbar1, tftype); 

  /////  PSPart::b1
  if( perm[nj_b1]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b1] ]; 
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b1)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b1)->second ] );
    tftype  = (perm[nj_b1]==-1 ? TFType::qLost : TFType::Unknown);
  }
  E       = solve(ps.lv(PSPart::q1) + ps.lv(PSPart::qbar1), DMT2, MB, dir, E_REC, accept);
  extend_PS( ps, PSPart::b1, E, MB , dir, perm[nj_b1], PSVar::cos_b1, PSVar::phi_b1, PSVar::E_b1, tftype ); 
  
  /////  PSPart::q2
  MEM::Object* lep = obs_leptons[ nl_q2 ]; 
  dir     = lep->p4().Vect().Unit(); 
  E       = lep->p4().E();
  extend_PS( ps, PSPart::q2, E, ML, dir, nl_q2, PSVar::cos_q2, PSVar::phi_q2, PSVar::E_q2, TFType::muReco, int(lep->getObs(Observable::CHARGE)) ); 
    
  /////  PSPart::qbar2
  dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_qbar2)->second ]) );
  double phi_n = x[ map_to_var.find(PSVar::phi_qbar2)->second ] + obs_mets[0]->p4().Phi();
  if     ( phi_n>+PI ) phi_n -= 2*PI;
  else if( phi_n<-PI ) phi_n += 2*PI;
  dir.SetPhi  ( phi_n );
  E    = solve( ps.lv(PSPart::q2), DMW2, ML, dir, E_REC, accept );
  extend_PS( ps, PSPart::qbar2, E, 0., dir, -1, PSVar::cos_qbar2, PSVar::phi_qbar2, PSVar::E_qbar2, TFType::MET  ); 

  /////  PSPart::b2
  if( perm[nj_b2]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b2] ];
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b2)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b2)->second ] );
    tftype  = (perm[nj_b2]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = solve( ps.lv(PSPart::q2) + ps.lv(PSPart::qbar2), DMT2, MB, dir, E_REC, accept );
  extend_PS( ps, PSPart::b2, E, MB, dir, perm[nj_b2], PSVar::cos_b2, PSVar::phi_b2, PSVar::E_b2,  tftype); 
  
  /////  PSPart::b
  if( perm[nj_b]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b] ];
    dir     = obj->p4().Vect().Unit();
    E_LOW   = obj->getObs(Observable::E_LOW_B) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_B);    
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b)->second ] );
    E_LOW   = MB;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_b]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_b)->second ]);
  extend_PS( ps, PSPart::b, E, MB, dir, perm[nj_b], PSVar::cos_b, PSVar::phi_b, PSVar::E_b,  tftype ); 

  /////  PSPart::bbar   
  if( perm[nj_bbar]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_bbar] ];
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    E_LOW   = obj->getObs(Observable::E_LOW_B) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_B);        
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_bbar)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_bbar)->second ] );
    E_LOW   = MB;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_bbar]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E    = hypo==Hypothesis::TTBB ?  
    E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_bbar)->second ]) : 
    solve( ps.lv(PSPart::b), DMH2, MB, dir, E_REC, accept);
  extend_PS( ps, PSPart::bbar, E, MB, dir, perm[nj_bbar], PSVar::cos_bbar, PSVar::phi_bbar, PSVar::E_bbar,  tftype ); 

  // protect against collinear radiation for TTBB
  if( hypo==Hypothesis::TTBB ){
    LV lv_b     = ps.lv(PSPart::b);
    LV lv_bbar  = ps.lv(PSPart::bbar);
    if( lv_b.Pt()<0. || lv_bbar.Pt()<0. || deltaR(lv_b, lv_bbar)<0.3 ){
      accept = -1;
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\tSkip this PS because of collinearity: [" 
	     << lv_b.Pt() << ", " << lv_bbar.Pt() << ", " << deltaR(lv_b, lv_bbar) << "]" 
	     << endl;
      }
#endif
    }
  }

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_LH(): END" << endl;
  }
#endif

  return accept;
}

void MEM::Integrand::extend_PS(MEM::PS& ps, const MEM::PSPart::PSPart& part, 
			       const double& E,  const double& M, const TVector3& dir,
			       const int& pos,  const PSVar::PSVar& var_cos, const PSVar::PSVar& var_phi, const PSVar::PSVar& var_E, 
			       const TFType::TFType& type, const int charge) const {

  double E_phys  = TMath::Max(E,M); 
  double P       = sqrt(E_phys*E_phys - M*M);
  ps.set(part,  MEM::GenPart(TLorentzVector( dir*P, E_phys ), type, charge ) );

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\t\tExtend phase-space point: adding variable " << static_cast<size_t>(part) << endl;
    if(  map_to_var.find(var_E)!=map_to_var.end() )
      cout << "\t\tE   = x[" <<  map_to_var.find(var_E)->second << "] = " << E << " GeV" << endl;
    else
      cout << "\t\tE   = SOLVE() = " << E << " GeV" << endl;
    if(   map_to_var.find(var_cos)!= map_to_var.end() && map_to_var.find(var_phi)!= map_to_var.end()  ){
      cout << "\t\tcos = x[" << map_to_var.find(var_cos)->second << "] = " << TMath::Cos(dir.Theta()) << endl;    
      cout << "\t\tphi = x[" << map_to_var.find(var_phi)->second << "] = " << dir.Phi() << endl;
    }
    else{
      cout << "\t\tUsing obs[" << pos << "]" << endl;
    }
  }
#endif

}

void MEM::Integrand::extend_PS_nodebug(MEM::PS& ps, const MEM::PSPart::PSPart& part, 
				       const double& E,  const double& M, const TVector3& dir) const {
  double E_phys  = TMath::Max(E,M); 
  double P       = sqrt(E_phys*E_phys - M*M);
  ps.set(part,  MEM::GenPart(TLorentzVector( dir*P, E_phys ), TFType::Unknown, 0 ) );
}


int MEM::Integrand::create_PS_LL(MEM::PS& ps, const double* x, const vector<int>& perm ) const {

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_LL(): START" << endl;
  }
#endif

  // corrupted phase space
  int accept{0};

  // store temporary values to build four-vectors
  double E     {0.};
  double E_LOW {0.};
  double E_HIGH{0.};
  double E_REC {numeric_limits<double>::max()};
  TVector3 dir (1.,0.,0.);
  TFType::TFType tftype = TFType::Unknown;

  // map a quark to an index inside the obs_ collections
  size_t nl_q1    =  map_to_part.find(PSPart::q1)->second; 
  size_t nj_b1    =  map_to_part.find(PSPart::b1)->second; 
  size_t nl_q2    =  map_to_part.find(PSPart::q2)->second; 
  size_t nj_b2    =  map_to_part.find(PSPart::b2)->second; 
  size_t nj_b     =  map_to_part.find(PSPart::b)->second;
  size_t nj_bbar  =  map_to_part.find(PSPart::bbar)->second;

  /////  PSPart::q1
  MEM::Object* lep1 = obs_leptons[ nl_q1 ];
  dir     = lep1->p4().Vect().Unit();
  E       = lep1->p4().E();
  extend_PS( ps, PSPart::q1, E, ML, dir, nl_q1, PSVar::cos_q1, PSVar::phi_q1, PSVar::E_q1, TFType::muReco, int(lep1->getObs(Observable::CHARGE)) ); 

  /////  PSPart::qbar1
  dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_qbar1)->second ]) );
  dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_qbar1)->second ] );
  E    = solve( ps.lv(PSPart::q1), DMW2, ML, dir, E_REC, accept );
  extend_PS( ps, PSPart::qbar1, E, 0., dir, -1, PSVar::cos_qbar1, PSVar::phi_qbar1, PSVar::E_qbar1, TFType::MET  ); 

  /////  PSPart::b1
  if( perm[nj_b1]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b1] ];
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b1)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b1)->second ] );
    tftype  = (perm[nj_b1]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = solve(ps.lv(PSPart::q1) + ps.lv(PSPart::qbar1), DMT2, MB, dir, E_REC, accept);
  extend_PS( ps, PSPart::b1, E, MB , dir, perm[nj_b1], PSVar::cos_b1, PSVar::phi_b1, PSVar::E_b1, tftype ); 
  
  /////  PSPart::q2
  MEM::Object* lep2 = obs_leptons[ nl_q2 ];
  dir     = lep2->p4().Vect().Unit();
  E       = lep2->p4().E();
  extend_PS( ps, PSPart::q2, E, ML, dir, nl_q2, PSVar::cos_q2, PSVar::phi_q2, PSVar::E_q2, TFType::muReco, int(lep2->getObs(Observable::CHARGE)) ); 

  /////  PSPart::qbar2
  dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_qbar2)->second ]) );
  dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_qbar2)->second ] );
  E    = solve( ps.lv(PSPart::q2), DMW2, ML, dir, E_REC, accept);
  extend_PS( ps, PSPart::qbar2, E, 0., dir, -1, PSVar::cos_qbar2, PSVar::phi_qbar2, PSVar::E_qbar2, TFType::MET  ); 

  //  PSVar::cos_b2, PSVar::phi_b2, PSVar::E_b2   
  if( perm[nj_b2]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b2] ];
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b2)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b2)->second ] );
    tftype  = (perm[nj_b2]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = solve( ps.lv(PSPart::q2) + ps.lv(PSPart::qbar2), DMT2, MB, dir, E_REC, accept );
  extend_PS( ps, PSPart::b2, E, MB, dir, perm[nj_b2], PSVar::cos_b2, PSVar::phi_b2, PSVar::E_b2, tftype ); 
  
  /////  PSPart::b 
  if( perm[nj_b]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b] ];
    dir     = obj->p4().Vect().Unit();
    E_LOW   = obj->getObs(Observable::E_LOW_B) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_B);
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b)->second ] );
    E_LOW   = MB;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_b]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_b)->second ]);
  extend_PS( ps, PSPart::b, E, MB, dir, perm[nj_b], PSVar::cos_b, PSVar::phi_b, PSVar::E_b,  tftype ); 

  /////  PSPart::bbar  
  if( perm[nj_bbar]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_bbar] ];
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    E_LOW   = obj->getObs(Observable::E_LOW_B) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_B);
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_bbar)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_bbar)->second ] );
    E_LOW   = MB;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_bbar]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E    = hypo==Hypothesis::TTBB ? 
    E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_bbar)->second ]) : 
    solve( ps.lv(PSPart::b), DMH2, MB, dir, E_REC, accept);
  extend_PS( ps, PSPart::bbar, E, MB, dir, perm[nj_bbar], PSVar::cos_bbar, PSVar::phi_bbar, PSVar::E_bbar,  tftype ); 

  // protect against collinear radiation for TTBB
  if( hypo==Hypothesis::TTBB ){
    LV lv_b     = ps.lv(PSPart::b);
    LV lv_bbar  = ps.lv(PSPart::bbar);
    if( lv_b.Pt()<0. || lv_bbar.Pt()<0. || deltaR(lv_b, lv_bbar)<0.3 ){
      accept = -1;
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\tSkip this PS because of collinearity: [" 
	     << lv_b.Pt() << ", " << lv_bbar.Pt() << ", " << deltaR(lv_b, lv_bbar) << "]" 
	     << endl;
      }
#endif
    }
  }


#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_LHL(): END" << endl;
  }
#endif

  return accept;
}

int MEM::Integrand::create_PS_HH(MEM::PS& ps, const double* x, const vector<int>& perm ) const { 
  
#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_HH(): START" << endl;
  }
#endif

  // corrupted phase space
  int accept{0};

  // store temporary values to build four-vectors
  double E     {0.};
  double E_LOW {0.};
  double E_HIGH{0.};
  double E_REC {numeric_limits<double>::max()};
  TVector3 dir (1.,0.,0.);
  TFType::TFType tftype;

  // map a quark to an index inside the obs_ collections
  size_t nj_q1    =  map_to_part.find(PSPart::q1)->second; 
  size_t nj_qbar1 =  map_to_part.find(PSPart::qbar1)->second; 
  size_t nj_b1    =  map_to_part.find(PSPart::b1)->second; 
  size_t nj_q2    =  map_to_part.find(PSPart::q2)->second; 
  size_t nj_qbar2 =  map_to_part.find(PSPart::qbar2)->second; 
  size_t nj_b2    =  map_to_part.find(PSPart::b2)->second; 
  size_t nj_b     =  map_to_part.find(PSPart::b)->second;
  size_t nj_bbar  =  map_to_part.find(PSPart::bbar)->second;

  /////  PSPart::q1
  if( perm[ nj_q1 ]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_q1] ]; 
    dir     = obj->p4().Vect().Unit();
    E_LOW   = obj->getObs(Observable::E_LOW_Q) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_Q);
    tftype  = TFType::qReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_q1)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_q1)->second ] );
    E_LOW   = MQ;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_q1]==-1 ? TFType::qLost : TFType::Unknown);
  }
  E       = E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_q1)->second ]);
  extend_PS( ps, PSPart::q1, E, MQ, dir, perm[nj_q1], PSVar::cos_q1, PSVar::phi_q1, PSVar::E_q1, tftype); 

  /////  PSPart::qbar1
  if( perm[nj_qbar1]>=0 ){
    dir     = obs_jets[ perm[nj_qbar1] ]->p4().Vect().Unit();
    tftype  = TFType::qReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_qbar1)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_qbar1)->second ] );
    tftype  = (perm[nj_qbar1]==-1 ? TFType::qLost : TFType::Unknown);
  }
  E    = solve( ps.lv(PSPart::q1), DMW2 , MQ, dir, E_REC, accept );
  extend_PS( ps, PSPart::qbar1, E, MQ, dir, perm[nj_qbar1], PSVar::cos_qbar1, PSVar::phi_qbar1, PSVar::E_qbar1,  tftype); 

  /////  PSPart::b1
  if( perm[nj_b1]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b1] ]; 
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b1)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b1)->second ] );
    tftype  = (perm[nj_b1]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = solve(ps.lv(PSPart::q1) + ps.lv(PSPart::qbar1), DMT2, MB, dir, E_REC, accept);
  extend_PS( ps, PSPart::b1, E, MB , dir, perm[nj_b1], PSVar::cos_b1, PSVar::phi_b1, PSVar::E_b1, tftype); 
  
  /////  PSPart::q2
  if( perm[ nj_q2 ]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_q2] ]; 
    dir     = obj->p4().Vect().Unit();
    E_LOW   = obj->getObs(Observable::E_LOW_Q) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_Q);
    tftype  = TFType::qReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_q2)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_q2)->second ] );
    E_LOW   = MQ;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_q2]==-1 ? TFType::qLost : TFType::Unknown);
  }
  E       = E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_q2)->second ]);
  extend_PS( ps, PSPart::q2, E, MQ, dir, perm[nj_q2], PSVar::cos_q2, PSVar::phi_q2, PSVar::E_q2, tftype ); 

  /////  PSPart::qbar2
  if( perm[nj_qbar2]>=0 ){
    dir = obs_jets[ perm[nj_qbar2] ]->p4().Vect().Unit();
    tftype  = TFType::qReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_qbar2)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_qbar2)->second ] );
    tftype  = (perm[nj_qbar2]==-1 ? TFType::qLost : TFType::Unknown);
  }
  E    = solve( ps.lv(PSPart::q2), DMW2 , MQ, dir, E_REC, accept );
  extend_PS( ps, PSPart::qbar2, E, MQ, dir, perm[nj_qbar2], PSVar::cos_qbar2, PSVar::phi_qbar2, PSVar::E_qbar2,  tftype ); 
  
  /////  PSPart::b2
  if( perm[nj_b2]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b2] ];
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b2)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b2)->second ] );
    tftype  = (perm[nj_b2]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = solve( ps.lv(PSPart::q2) + ps.lv(PSPart::qbar2), DMT2, MB, dir, E_REC, accept );
  extend_PS( ps, PSPart::b2, E, MB, dir, perm[nj_b2], PSVar::cos_b2, PSVar::phi_b2, PSVar::E_b2,  tftype); 
  
  /////  PSPart::b
  if( perm[nj_b]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_b] ];
    dir     = obj->p4().Vect().Unit();
    E_LOW   = obj->getObs(Observable::E_LOW_B) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_B);    
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_b)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_b)->second ] );
    E_LOW   = MB;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_b]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E       = E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_b)->second ]);
  extend_PS( ps, PSPart::b, E, MB, dir, perm[nj_b], PSVar::cos_b, PSVar::phi_b, PSVar::E_b,  tftype ); 

  /////  PSPart::bbar   
  if( perm[nj_bbar]>=0 ){
    MEM::Object* obj = obs_jets[ perm[nj_bbar] ];
    dir     = obj->p4().Vect().Unit();
    E_REC   = obj->p4().E();
    E_LOW   = obj->getObs(Observable::E_LOW_B) ;
    E_HIGH  = obj->getObs(Observable::E_HIGH_B);        
    tftype  = TFType::bReco;
  }
  else{
    dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_bbar)->second ]) );
    dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_bbar)->second ] );
    E_LOW   = MB;
    E_HIGH  = TMath::Min(cfg.emax, 2*MEM::TF_ACC_param[1]/TMath::Sin(dir.Theta()) );
    tftype  = (perm[nj_bbar]==-1 ? TFType::bLost : TFType::Unknown);
  }
  E    = hypo==Hypothesis::TTBB ?  
    E_LOW + (E_HIGH-E_LOW)*(x[ map_to_var.find(PSVar::E_bbar)->second ]) : 
    solve( ps.lv(PSPart::b), DMH2, MB, dir, E_REC, accept);
  extend_PS( ps, PSPart::bbar, E, MB, dir, perm[nj_bbar], PSVar::cos_bbar, PSVar::phi_bbar, PSVar::E_bbar,  tftype); 

  // protect against collinear radiation for TTBB
  if( hypo==Hypothesis::TTBB ){
    LV lv_b     = ps.lv(PSPart::b);
    LV lv_bbar  = ps.lv(PSPart::bbar);
    if( lv_b.Pt()<0. || lv_bbar.Pt()<0. || deltaR(lv_b, lv_bbar)<0.3 ){
      accept = -1;
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\tSkip this PS because of collinearity: [" 
	     << lv_b.Pt() << ", " << lv_bbar.Pt() << ", " << deltaR(lv_b, lv_bbar) << "]" 
	     << endl;
      }
#endif
    }
  }
  
#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_HH(): END" << endl;
  }
#endif

  return accept;
}


int MEM::Integrand::create_PS_TTH(MEM::PS& ps, const double* x, const vector<int>& perm ) const {

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_TTH(): START" << endl;
  }
#endif

  // corrupted phase space
  int accept{0};

  // store temporary values to build four-vectors
  double P{0.};
  TVector3 dir(1.,0.,0.);

  dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_t)->second ]) );
  dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_t)->second ] );
  P = x[ map_to_var.find(PSVar::P_t)->second ];
  extend_PS_nodebug( ps, PSPart::t, sqrt(P*P+MTOP2), MTOP, dir);

  dir.SetTheta( TMath::ACos( x[ map_to_var.find(PSVar::cos_tbar)->second ]) );
  dir.SetPhi  ( x[ map_to_var.find(PSVar::phi_tbar)->second ] );
  P = x[ map_to_var.find(PSVar::P_tbar)->second ];
  extend_PS_nodebug( ps, PSPart::tbar, sqrt(P*P+MTOP2), MTOP, dir);

  double Px = -ps.lv(PSPart::t).Px()-ps.lv(PSPart::tbar).Px();
  double Py = -ps.lv(PSPart::t).Py()-ps.lv(PSPart::tbar).Py();
  dir = TVector3(Px,Py, x[ map_to_var.find(PSVar::Pz_h)->second] );
  extend_PS_nodebug( ps, PSPart::h, sqrt(dir.Mag2()+MH2), MH, dir.Unit());


#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::create_PS_LH(): END" << endl;
  }
#endif

  return accept;
}



double MEM::Integrand::probability(const double* x, const std::size_t& n_perm ) const {

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::probability(): START" << endl;
  }
#endif

  // the total probability
  double p{1.0};
  if( !cfg.int_code ) return p;

  // create phas-space point and test if it is physical
  PS ps(ps_dim);
  int accept = create_PS(ps, x, perm_indexes_assumption[n_perm]);

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ) ps.print(cout);
#endif
  
  if( accept<0 ){
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration){
      cout << "\tCORRUPTED PS (no solution): return 0." << endl;      
    }
#endif
    ++(const_cast<Integrand*>(this)->n_skip);
    return 0.;
  }

  p *= constants();

  p *= transfer(ps, perm_indexes_assumption[n_perm], accept);
  if(cfg.do_prefit>1 && prefit_step==0) return p;

  if(cfg.tf_suppress && accept>=cfg.tf_suppress){
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration ){
      cout << "\tTransfer functions out-of-range " << accept << " times: return from this PS before calculatin matrix()" << endl;
    }
#endif
    return 0.;
  }

  p *= matrix(ps);

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::probability(): END" << endl;
  }
#endif

  return p;
}

double MEM::Integrand::matrix(const PS& ps) const {

  if( fs==FinalState::TTH ){ return matrix_nodecay(ps); }

  double m{1.};

  double x1  {0.}; // x1 fraction
  double x2  {0.}; // x2 fraction

  LV lv_q1    = ps.lv(PSPart::q1);
  LV lv_qbar1 = ps.lv(PSPart::qbar1);
  LV lv_b1    = ps.lv(PSPart::b1);
  LV lv_q2    = ps.lv(PSPart::q2);
  LV lv_qbar2 = ps.lv(PSPart::qbar2);
  LV lv_b2    = ps.lv(PSPart::b2);
  LV lv_b     = ps.lv(PSPart::b);
  LV lv_bbar  = ps.lv(PSPart::bbar);

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\t\tFilling m..." << endl;
    cout << "\t\tCheck masses: m(W1)=" << (lv_q1+lv_qbar1).M() << ", m(t1)=" << (lv_q1+lv_qbar1+lv_b1).M()
	 << ", m(W2)=" << (lv_q2+lv_qbar2).M() << ", m(t2)=" << (lv_q2+lv_qbar2+lv_b2).M() 
	 << ", m(H)=" << (lv_b+lv_bbar).M() << endl;
  }
#endif

  m *= t_decay_amplitude(lv_q1, lv_qbar1, lv_b1, ps.charge(PSPart::q1) );
  m *= t_decay_amplitude(lv_q2, lv_qbar2, lv_b2, ps.charge(PSPart::q2) );
  m *= H_decay_amplitude(lv_b, lv_bbar);
  m *= scattering( lv_q1+lv_qbar1+lv_b1, lv_q2+lv_qbar2+lv_b2, lv_b, lv_bbar, x1, x2);
  m *= pdf( x1, x2 , lv_b.Pt()+lv_bbar.Pt() );

  if( TMath::IsNaN(m) ){
    cout << "\tmatrix() returned a NaN..." << endl;
    return 0.;
  }

  return m;
}

double MEM::Integrand::matrix_nodecay(const PS& ps) const {
  double m{1.};

  double x1  {0.}; // x1 fraction
  double x2  {0.}; // x2 fraction

  LV lv_t     = ps.lv(PSPart::t);
  LV lv_tbar  = ps.lv(PSPart::tbar);
  LV lv_b     = fs==FinalState::TTH ? ps.lv(PSPart::h)      : ps.lv(PSPart::b);
  LV lv_bbar  = fs==FinalState::TTH ? LV(1e-06,0.,0.,1e-06) : ps.lv(PSPart::bbar);
  LV lv_h     = lv_b+lv_bbar;

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\t\tFilling m..." << endl;
    cout << "\t\tCheck masses: m(t1)=" << lv_t.M() << ", m(t2)=" << lv_tbar.M() << ", m(h)= " << (lv_b+lv_bbar).M() << endl;
  }
#endif

  m *= scattering( lv_t, lv_tbar, lv_b, lv_bbar, x1, x2);
  m *= pdf( x1, x2 , lv_h.Pt() );

  // check this: it won't work for tt+bb
  double Jac = lv_t.Beta()*lv_t.Vect().Mag()/2 * lv_tbar.Beta()*lv_tbar.Vect().Mag()/2 * (1./lv_h.E()/2);
  m *= Jac;

  if( TMath::IsNaN(m) ){
    cout << "\tA NaN occurred while evaluation m..." << endl;
    return 0.;
  }

  return m;
}

double MEM::Integrand::transfer(const PS& ps, const vector<int>& perm, int& accept) const {

  double w{1.};
  if( !(cfg.int_code&IntegrandType::Transfer) ) return w;
  if( fs==FinalState::TTH ) return w;

  double nu_x {0.}; // total nu's px
  double nu_y {0.}; // total nu's py
  double corr_nu_x {0.}; // sum of dPx
  double corr_nu_y {0.}; // sum of dPy
  double rho_x{0.}; // recoil px
  double rho_y{0.}; // recoil py
  double pT_x {0.}; // pT px
  double pT_y {0.}; // pT py

  // subtract MET from the recoil
  rho_x -= obs_mets[0]->p4().Px();
  rho_x -= obs_mets[0]->p4().Py();

  // Dealing with jets and leptons
  PSMap::const_iterator p;
  for( p = ps.begin() ; p != ps.end() ; ++p ){    

    MEM::Object* obj = nullptr;
    if( isLepton  ( p->second.type ) ){
      obj = obs_leptons[ map_to_part.find(p->first)->second ];
      // subtract from recoil
      rho_x -= obj->p4().Px();
      rho_y -= obj->p4().Py();

      // subtract from total pT
      pT_x  -= p->second.lv.Px();
      pT_y  -= p->second.lv.Py();
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\tDealing with a lepton..." << endl;
	cout << "\t\trho_x -= " << obj->p4().Px() << ", pT_x -= " << p->second.lv.Px() << endl;
	cout << "\t\trho_y -= " << obj->p4().Py() << ", pT_y -= " << p->second.lv.Py() << endl;
      }
#endif
      continue;
    }

    if( isNeutrino( p->second.type ) ){
      // add up to neutrino
      nu_x += p->second.lv.Px();
      nu_y += p->second.lv.Py();

      // subtract from total pT
      pT_x -= p->second.lv.Px();
      pT_y -= p->second.lv.Py();
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\tDealing with a neutrino..." << endl;
	cout << "\t\tnu_x -= " <<  p->second.lv.Px() << ", pT_x -= " << p->second.lv.Px() << endl;
        cout << "\t\tnu_y -= " <<  p->second.lv.Py() << ", pT_y -= " << p->second.lv.Py() << endl;
      }
#endif
      continue;
    }      

    // observables and generated-level quantities
    // if the parton is matched, test value of jet energy

    // subtract from total pT
    pT_x  -= p->second.lv.Px();
    pT_y  -= p->second.lv.Py();    

    const double e_gen  {p->second.lv.E()}; 
    const double pt_gen {p->second.lv.Pt()}; 
    const double eta_gen{p->second.lv.Eta()};     
    int jet_indx  = perm[ map_to_part.find(p->first)->second ];
    double e_rec{0.};

    if( jet_indx>=0 ) {
      obj    = obs_jets[ jet_indx ];
      e_rec  = obj->p4().E();
      rho_x -= obj->p4().Px();
      rho_y -= obj->p4().Py();
      
      //Enux = Etot * sin(theta) * cos(phi)
      //Px = |P| * sin(theta) * cos(phi) -> Px/P = sin(theta) * cos(phi)
      corr_nu_x += (e_rec-e_gen)*obj->p4().Px()/obj->p4().P();
      corr_nu_y += (e_rec-e_gen)*obj->p4().Py()/obj->p4().P();

      // if this flag is true, the TF are multiplied by the range step-functions
      // N.B false by default
      if(cfg.tf_in_range){
	double e_gen_min, e_gen_max;
	if( obj->isSet(Observable::BTAG) && 
	    obj->getObs(Observable::BTAG)<0.5 ){
	  e_gen_min = obj->getObs(Observable::E_LOW_Q); 
	  e_gen_max = obj->getObs(Observable::E_HIGH_Q); 
	}
	else if(obj->isSet(Observable::BTAG) && 
		obj->getObs(Observable::BTAG)>0.5){
	  e_gen_min = obj->getObs(Observable::E_LOW_B); 
	  e_gen_max = obj->getObs(Observable::E_HIGH_B); 	  
	}
	else{
	  e_gen_min = 0.;
	  e_gen_max = numeric_limits<double>::max();
	}
	if( !(e_gen>=e_gen_min && e_gen<=e_gen_max) ){
	  w *= 0.;
	  return w;
	}
      }

#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\tDealing with a jet..." << endl;
	cout << "\t\trho_x -= " << obj->p4().Px() << ", pT_x -= " << p->second.lv.Px() << endl;
	cout << "\t\trho_y -= " << obj->p4().Py() << ", pT_y -= " << p->second.lv.Py() << endl;
	cout << "\t\tdE_x   = " << (e_rec-e_gen)*obj->p4().Px()/obj->p4().P() << endl;
	cout << "\t\tdE_y   = " << (e_rec-e_gen)*obj->p4().Py()/obj->p4().P() << endl;
      }
#endif
    } //jet index >= 0
    else{
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\tDealing with a missed jet..." << endl;
	cout << "\t\trho_x -= " << 0 << ", pT_x -= " << p->second.lv.Px() << endl;
	cout << "\t\trho_y -= " << 0 << ", pT_y -= " << p->second.lv.Py() << endl;
      }
#endif
    }

    // build x,y vectors 
    double y[1] = { e_rec };
    double x[2] = { e_gen, eta_gen };
    
    //Try to calculate using externally supplied transfer functions
    //N.B.: new transfer functions are functions of jet pt
    if (cfg.transfer_function_method == TFMethod::External && 
	obj!=nullptr && obj->getNumTransferFunctions()>0) {
      x[0] = pt_gen;
      double _w = transfer_function2( obj, x, p->second.type, accept, cfg.tf_offscale, false, debug_code );
      w *= _w;
    } 

    //Calculate using reco efficiency
    //N.B.: new transfer functions are functions of jet pt
    else if (cfg.transfer_function_method == TFMethod::External && 
	     obj==nullptr) {

      int eta_bin = eta_to_bin(eta_gen, true);
      // if outside acceptance, return 1.0
      if( eta_bin<0 ){
	w *= 1.0;
	continue;
      }

      //pass the efficienty function as a pointer, 
      //need to remove const modifier as TF1::Eval does not specify const
      x[0] = pt_gen;      
      const TF1* tf = get_tf_global(p->second.type, eta_bin);
      assert(tf != nullptr);
      double _w = transfer_function2( const_cast<TF1*>(tf), x, p->second.type, accept, cfg.tf_offscale, false, debug_code );
      w *= _w;
    } 

    //Calculate using internal transfer functions
    else {
      double _w = transfer_function( y, x, p->second.type, accept, cfg.tf_offscale, debug_code );
      w *= _w;
    }
  }

  // Dealing with the MET
  double y_MET[2] = { obs_mets[0]->p4().Px(), obs_mets[0]->p4().Py() };
  double x_Nu [2] = { nu_x-corr_nu_x, nu_y-corr_nu_y };
  if( !(cfg.int_code&IntegrandType::Recoil) ){
    x_Nu[0] += corr_nu_x;
    x_Nu[1] += corr_nu_y;
  }
  w *= transfer_function( y_MET, x_Nu, TFType::MET, accept, cfg.tf_offscale, debug_code );

  // Dealing with the recoil
  double y_rho[1] = { (extra_jets>0 ? TF_RECOIL_param[2]+1. : sqrt(rho_x*rho_x + rho_y*rho_y)) };
  double x_pT [1] = { sqrt(pT_x*pT_x + pT_y*pT_y)  };
  if( cfg.int_code&IntegrandType::Sudakov ) 
    w *= transfer_function( y_rho, x_pT, TFType::Recoil, accept, cfg.tf_offscale, debug_code );

  if( TMath::IsNaN(w) ){
    cout << "\ttransfer() returned a NaN..." << endl;
    return 0.;
  }

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration){
    cout << "\tTotal transfer function: " << w 
	 << " (" << accept << " functions are out-of-range by more than " << cfg.tf_offscale << " sigmas" << endl;
  }
#endif

  return w;
}


double MEM::Integrand::scattering(const TLorentzVector& top, const TLorentzVector& atop, const TLorentzVector& b1, const TLorentzVector& b2,
				  double& x1, double& x2) const{

  // return value (address passed to OpenLoops)
  double M2{1.};  

  // temporary objects
  TLorentzVector t, tx, b, bx, h, sum;
  t.SetPtEtaPhiM ( top.Pt(),     top.Eta(),     top.Phi(),    MTOP);
  tx.SetPtEtaPhiM( atop.Pt(),    atop.Eta(),    atop.Phi(),   MTOP);
  b.SetPtEtaPhiM ( b1.Pt(),      b1.Eta(),      b1.Phi(),       0.);
  bx.SetPtEtaPhiM( b2.Pt(),      b2.Eta(),      b2.Phi(),       0.);
  h.SetPtEtaPhiM ( (b1+b2).Pt(), (b1+b2).Eta(), (b1+b2).Phi(),  MH);

  // the total sum (needed to get the boost factor);
  TLorentzVector vSum = hypo==Hypothesis::TTH ? t+tx+h : t+tx+b+bx;

  if(vSum.E()>cfg.sqrts){
    x1 = .99;
    x2 = .99;
    return 0.;
  }

  // boost such that SumPx = SumPy = 0
  TVector3 boostPt( vSum.Px()/vSum.E(), vSum.Py()/vSum.E(), 0.0 );
  bool apply_boost{vSum.Px()>1. || vSum.Py()>1.};

  if( apply_boost ){ // some tolerance!
    t.Boost  ( -boostPt );
    tx.Boost ( -boostPt );
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration && apply_boost){
      cout << "\t\tBoost system along the (x,y) plane: beta = (" 
	   << boostPt.Px() << ", " << boostPt.Py() << ", " << boostPt.Pz() << ")" << endl;
    }
#endif
  }

  if(hypo==Hypothesis::TTH){
    h.Boost  ( -boostPt );
    // fix for rounding
    double hPx = -(t.Px() + tx.Px());
    double hPy = -(t.Py() + tx.Py());
    double hPz = h.Pz();
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration && apply_boost){
      cout << "\t\tRounding Higgs px by " << (hPx-h.Px())/h.Px()*100 << "%" << endl;
      cout << "\t\tRounding Higgs py by " << (hPy-h.Py())/h.Py()*100 << "%" << endl;
    }
#endif
    h.SetPxPyPzE( hPx, hPy, hPz, sqrt(hPx*hPx + hPy*hPy + hPz*hPz + MH*MH) );    
    sum = t+tx+h;
  }
  else{
    b.Boost  ( -boostPt );
    bx.Boost ( -boostPt );
    // fix for rounding
    double bPx = -(t.Px() + tx.Px() + bx.Px());
    double bPy = -(t.Py() + tx.Py() + bx.Py());
    double bPz = b.Pz();
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration && apply_boost){
      cout << "\t\tRounding b px by " << (bPx-b.Px())/b.Px()*100 << "%" << endl;
      cout << "\t\tRounding b py by " << (bPy-b.Py())/b.Py()*100 << "%" << endl;
    }
#endif
    b.SetPxPyPzE( bPx, bPy, bPz, sqrt(bPx*bPx + bPy*bPy + bPz*bPz ) );
    sum =  t+tx+b+bx;
  }
  
  // update x1 and x2
  double E  = sum.E();
  double Pz = sum.Pz();
  x1 = ( Pz + E)/cfg.sqrts;
  x2 = (-Pz + E)/cfg.sqrts;
  if( !(cfg.int_code&IntegrandType::ScattAmpl) ) return M2;

  // create gluon p4s
  TLorentzVector g1 = TLorentzVector(0.,0.,  (E+Pz)/2., (E+Pz)/2.);
  TLorentzVector g2 = TLorentzVector(0.,0., -(E-Pz)/2., (E-Pz)/2.);

  // needed to interface with OpenLoops
  double ccP_0[20] = {
    g1.E(), g1.Px(), g1.Py(), g1.Pz(),
    g2.E(), g2.Px(), g2.Py(), g2.Pz(),
    h.E(),  h.Px(),  h.Py(),  h.Pz(),
    t.E(),  t.Px(),  t.Py(),  t.Pz(),
    tx.E(), tx.Px(), tx.Py(), tx.Pz()
  };

  double ccP_1[24] = {
    g1.E(), g1.Px(), g1.Py(), g1.Pz(),
    g2.E(), g2.Px(), g2.Py(), g2.Pz(),
    t.E(),  t.Px(),  t.Py(),  t.Pz(),
    tx.E(), tx.Px(), tx.Py(), tx.Pz(),
    b.E(),  b.Px(),  b.Py(),  b.Pz(),
    bx.E(), bx.Px(), bx.Py(), bx.Pz(),
  };

  // call OpenLoops functions
  switch(hypo){
  case Hypothesis::TTH:
    pphttxcallme2born_  (const_cast<double*>(&M2), ccP_0, const_cast<double*>(&MTOP), const_cast<double*>(&MH));
    break; 
  case Hypothesis::TTBB:
    ppttxbbxcallme2born_(const_cast<double*>(&M2), ccP_1, const_cast<double*>(&MTOP), const_cast<double*>(&MH));
    break;
  default:
    break;
  }

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\t\tTotal (px,py,pz,E) = (" <<  sum.Px() << "," <<  sum.Py()<< "," <<  sum.Pz()<< "," <<  sum.E() << ")" << endl;
    cout << "\t\tGluons (x1,x2)     = (" << x1 << "," << x2 << ")" << endl;
    cout << "\t\tM2 (OpenLoops)     = " << M2 << endl;
  }
#endif

  return M2;
}

double MEM::Integrand::pdf(const double& x1, const double& x2, const double& dynamical ) const{
  double p{1.};
  if( !(cfg.int_code&IntegrandType::PDF) ) return p;

  if(x1>0.99 || x2>0.99) return 0.;

  double Q{2*MTOP};
  switch( hypo ){
  case Hypothesis::TTH:
    Q = (2*MTOP + MH)/2;
    break;
  case Hypothesis::TTBB:
    Q = TMath::Sqrt( 4*MTOP*MTOP + TMath::Power(dynamical, 2) );
    break;
  default:
    break;
  }

  double f1 =  LHAPDF::xfx(1, x1, Q, 0)/x1;
  double f2 =  LHAPDF::xfx(1, x2, Q, 0)/x2;

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\t\tPDF(x1,Q)*PDF(x2,Q) = " << f1 << "*" << f2 << endl;
  }
#endif

  p = (f1*f2)/(x1*x2);

  return p;
}

double MEM::Integrand::constants() const {
  double p{1.};
  if( !(cfg.int_code&IntegrandType::Constant) ) return p;
  return TMath::Power((2.*PI), int(4-3*ps_dim))/(cfg.sqrts*cfg.sqrts*cfg.sqrts*cfg.sqrts);
}

double MEM::Integrand::t_decay_amplitude(const TLorentzVector& q, const TLorentzVector& qbar, const TLorentzVector& b, const int& charge_q) const{
  double p{1.};
  if( !(cfg.int_code&IntegrandType::DecayAmpl) ) return p;

  p *= BWTOP;

  TLorentzVector w = q+qbar;
  TLorentzVector t = w+b;
  double InvJac = TMath::Abs(2*MW2/qbar.E()*( w.E() - w.Vect().Dot( b.Vect().Unit() )/b.Beta() ));
  double Jac    = (1./InvJac) * q.Vect().Mag() * qbar.Vect().Mag() * b.Vect().Mag() / (2*2*2);
  if( cfg.int_code&IntegrandType::Jacobian ) 
    p *= Jac;

  double x_e1 = 2*(q*t)/MTOP2;
  double x_e2 = 2*(qbar*t)/MTOP2;

  // if the flavour is determined, use formula. Otherwise, take average.
  double m2   = charge_q!=0 ? x_e1*(1-MUB-x_e1) : 0.5*(x_e1*(1-MUB-x_e1) + x_e2*(1-MUB-x_e2));
  m2 *= (32*PI*MTOP4*GEWK4/(MW*GW));
  if(m2<0){
    cout << "\tIntegrand::t_decay_amplitude() returned negative |M2|..." << endl;
    return 0.;
  }
  p *= m2;

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::t_decay_amplitude():" << endl;
    cout << "\t\tBreit-Wigner top = " << BWTOP << " GeV^-2" << endl;
    cout << "\t\tJacobian (Eqbar,Eb) -> (m2_qq, m2_qqb) = " << Jac << " GeV" << endl;
    cout << "\t\t|M2|(t->bqq') = " << m2 << (charge_q==0?" (charge symmetrised)" : "")  << endl;
    cout << "\t\tTotal = " << p << " GeV^-1" << endl;
  }
#endif

  return p;
}

double MEM::Integrand::H_decay_amplitude(const TLorentzVector& b, const TLorentzVector& bbar) const{
  double p{1.};
  if( !(cfg.int_code&IntegrandType::DecayAmpl) ) return p;

  double InvJac{1.0};
  double m2{1.0};
  if( hypo==Hypothesis::TTH ){
    p *= BWH;
    InvJac = TMath::Abs(2*(b.E() - b.Vect().Dot( bbar.Vect().Unit() )/bbar.Beta()));
    m2 = 2*YB2*MH2*PSHBB;
  }
  double Jac    = (1./InvJac) * b.Vect().Mag() * bbar.Vect().Mag() / (2*2);
  if( cfg.int_code&IntegrandType::Jacobian ) 
    p *= Jac; 

  p *= m2;

#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::H_decay_amplitude():" << endl;
    if( hypo==Hypothesis::TTH ){
      cout << "\t\tBreit-Wigner Higgs = " << BWH << " GeV^-2" << endl;
      cout << "\t\tJacobian (Eb,Ebbar) -> (E_b, m2_bb) = " << Jac << " GeV" << endl;
      cout << "\t\t|M2|(t->bqq') = " << m2 << " GeV^2" << endl;
      cout << "\t\tTotal = " << p << " GeV" << endl;
    }
    else{
      cout << "\t\tJacobian = " << Jac << " GeV^2" << endl;
    }
  }
#endif

  return p;
}

double MEM::Integrand::solve(const LV& p4_w, const double& DM2, const double& M, const TVector3& e_b, const double& target, int& accept ) const {

  double a     = DM2/p4_w.E();
  double b     = TMath::Cos(p4_w.Angle(e_b));
  if( M<1e-03 ){
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration ){
      cout << "\t\tUse masless formula: " << a << "/(1-" << b << ")=" << a/(1-b) << endl;  
    }
#endif
    if( b<1. ) return a/(1-b);
    else{
      accept = -1;
      return numeric_limits<double>::max();
    } 
  }
  
  // use adimensional 'a', account for velocity<1
  a /= M;
  b *= p4_w.Beta();
  double a2    = a*a;
  double b2    = b*b;

  // this is needed to test the solutions
  double discr = a2 + b2 - a2*b2 - 1; 

  // make sure there is >0 solutions
  if( (a2 + b2 - 1)<0. ){
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration ){
      cout << "\t\t(a2 + b2 - 1)<0. return max()" << endl;
    }
#endif
    accept = -1;
    return numeric_limits<double>::max();
  }

  // the roots
  double g_p = (a + TMath::Abs(b)*sqrt(a2 + b2 - 1))/(1-b2);
  double g_m = (a - TMath::Abs(b)*sqrt(a2 + b2 - 1))/(1-b2);

  // make sure this is >1 ( g_m < g_p )
  if( g_p < 1.0 ){
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration ){
      cout << "\t\tg_p=" << g_p << ": return max()" << endl;
    }
#endif
    accept = -1;
    return numeric_limits<double>::max();
  }

  // remove unphysical root
  if( g_m < 1.0) g_m = g_p;

  // test for the roots
  switch( b>0 ){
  case true :
    if( discr<0 ){
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\t\tb>0 AND discr<0: return root closest to target" << endl;
	LV p4_b(e_b*(sqrt(g_p*g_p-1)*M), g_p*M);
	cout << "\t\tTwo solutions: " << endl;
	cout << "\t\tEb+ = " <<  g_p*M << endl;
	cout << "\t\t\tp4_w = (" <<  p4_w.Px() << ", " << p4_w.Py() << ", " << p4_w.Pz() << ", " << p4_w.E() << "), M=" << p4_w.M() << endl;
	cout << "\t\t\tp4_b = (" <<  p4_b.Px() << ", " << p4_b.Py() << ", " << p4_b.Pz() << ", " << p4_b.E() << "), M=" << p4_b.M() << endl;
	cout << "\t\t\tp4_t = (" <<  (p4_w+p4_b).Px() << ", " << (p4_w+p4_b).Py() << ", " << (p4_w+p4_b).Pz() << ", " << (p4_w+p4_b).E() << "), M=" << (p4_w+p4_b).M() << endl;
	cout << "\t\t\tAngle=" << p4_w.Vect().Angle(p4_b.Vect()) << endl;
	p4_b = LV(e_b*(sqrt(g_m*g_m-1)*M), g_m*M);
	cout << "\t\tEb- = " <<  g_m*M << endl;
	cout << "\t\t\tp4_w = (" <<  p4_w.Px() << ", " << p4_w.Py() << ", " << p4_w.Pz() << ", " << p4_w.E() << "), M=" << p4_w.M() << endl;
	cout << "\t\t\tp4_b = (" <<  p4_b.Px() << ", " << p4_b.Py() << ", " << p4_b.Pz() << ", " << p4_b.E() << "), M=" << p4_b.M() << endl;
	cout << "\t\t\tp4_t = (" <<  (p4_w+p4_b).Px() << ", " << (p4_w+p4_b).Py() << ", " << (p4_w+p4_b).Pz() << ", " << (p4_w+p4_b).E() << "), M=" << (p4_w+p4_b).M() << endl;
	cout << "\t\t\tAngle=" << p4_w.Vect().Angle(p4_b.Vect()) << endl;
      }
#endif
      return ( TMath::Abs(target-g_p)<TMath::Abs(target-g_m) ? g_p*M : g_m*M );
    }
#ifdef DEBUG_MODE
    if( debug_code&DebugVerbosity::integration ){
      cout << "\t\tb>0 AND discr>0: return g_p*M" << endl;
    }
#endif
    return g_p*M;
    break;
  case false:
    if( discr>0 ){
#ifdef DEBUG_MODE
      if( debug_code&DebugVerbosity::integration ){
	cout << "\t\tb<0 AND discr>0: return g_m*M" << endl;
      }
#endif
      return g_m*M;
    }
    break;
  }
#ifdef DEBUG_MODE
  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::solve(): END" << endl;
  }
#endif

  accept = -1;
  return numeric_limits<double>::max();
}

void MEM::Integrand::setup_minimizer(){
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(0);
  return;
}

void MEM::Integrand::refine_minimization(std::size_t& num_trials, const ROOT::Math::Functor& toIntegrate, const std::size_t& npar, double* xL, double* xU){

  if(debug_code&DebugVerbosity::init ){
    cout << "\trefine_minimization(): try to set quark energies to constant" << endl;
  }

  while(num_trials<2){
    
    if(minimizer!=nullptr)
      delete minimizer ;

    minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    setup_minimizer();
    minimizer->SetFunction(toIntegrate);
	
    for(size_t np = 0; np < npar ; ++np){
      string var_name = "par_"+std::to_string(np);
      if( xU[np]==1. && xL[np]==0. ){
	minimizer->SetFixedVariable(np, var_name, (xU[np]+xL[np])/2. );
	if(debug_code&DebugVerbosity::init )
	  printf("\tParam[%lu] = %s fixed to %.2f\n", np, var_name.c_str(),  (xU[np]+xL[np])/2); 
      }
      else{
	minimizer->SetLimitedVariable(np, var_name.c_str() , (xU[np]+xL[np])/2 , 5e-02 , xL[np] , xU[np] );
	if(debug_code&DebugVerbosity::init )
	  printf("\tParam[%lu] = %s set to %.2f. Range: [%.2f,%.2f]\n", np, var_name.c_str(),  (xU[np]+xL[np])/2 , xL[np], xU[np] ); 	    
      }
    }
    ++num_trials;
    minimizer->Minimize(); 
  }
  
  return;
}


void MEM::Integrand::smear_met(){
  
  if( debug_code&DebugVerbosity::init ){
      cout << "\tInput MET vector will be smeared using MEM::transfer_function_smear:" << endl;
  }

  // needed to get a fresh seed every time
  gRandom->SetSeed(0);
  
  TF2 tf2("met_tf", MEM::transfer_function_smear,  
	  obs_mets[0]->p4().Px()-100,  obs_mets[0]->p4().Px()+100, 
	  obs_mets[0]->p4().Py()-100,  obs_mets[0]->p4().Py()+100, 
	  3);  
  tf2.SetNpx(400); tf2.SetNpy(400);
  tf2.SetParameter(0, obs_mets[0]->p4().Px() );
  tf2.SetParameter(1, obs_mets[0]->p4().Py() );
  tf2.SetParameter(2, static_cast<int>(TFType::MET) );

  double met_x{0.};
  double met_y{0.};
  tf2.GetRandom2(met_x, met_y);    
  
  if( debug_code&DebugVerbosity::init ){
    cout << "\t\t MET (Px,Py,Pz,E): (" << obs_mets[0]->p4().Px() << ", "
	 << obs_mets[0]->p4().Py() << ", " << obs_mets[0]->p4().Pz() << ", "
	 << obs_mets[0]->p4().E() << ") --> " ;
  }

  obs_mets[0]->setp4( LV(met_x, met_y, 0., sqrt(met_x*met_x + met_y*met_y)) );

  if( debug_code&DebugVerbosity::init ){
    cout << "(" << obs_mets[0]->p4().Px() << ", " << obs_mets[0]->p4().Py() << ", " 
	 << obs_mets[0]->p4().Pz() << ", " << obs_mets[0]->p4().E() << ")" << endl;
  }
  return;
}

void MEM::Integrand::smear_jet(MEM::Object* j, const bool& external_tf){

  // needed to get a fresh seed every time
  gRandom->SetSeed(0);

  TFType::TFType jet_type = (j->isSet(Observable::BTAG) && j->getObs(Observable::BTAG)>0.5) 
    ? TFType::bReco : TFType::qReco;
  
  int jet_type_int = static_cast<int>(jet_type);
  
  if( debug_code&DebugVerbosity::init ){
    cout << "\tInput jet of type (" << jet_type_int << ") will be smeared using " <<
      (external_tf? "MEM::transfer_function2" : "MEM::transfer_function_smear") << endl;
    cout << "\t\t Jet (Px,Py,Pz,E): (" << j->p4().Px() << ", " << j->p4().Py() 
	 << ", " << j->p4().Pz() << ", " << j->p4().E() << ") --> " ;	
  }
  
  // if using builtin method, wrap the transfer_function inside a TF1/TF2 
  if( !external_tf || cfg.transfer_function_method==TFMethod::Builtin){

    TF1 tf1("jet_tf", MEM::transfer_function_smear, j->p4().E()/5., j->p4().E()*5., 3);
    tf1.SetNpx(1000);
    tf1.SetParameter(0, j->p4().E()   );
    tf1.SetParameter(1, j->p4().Eta() );
    tf1.SetParameter(2, jet_type_int );

    double e_ran = tf1.GetRandom();
    int count{0};
    while( count<100 && e_ran*TMath::Sin(j->p4().Theta()) < MEM::TF_ACC_param[1] ){
      e_ran = tf1.GetRandom();
      ++count;
    }

    if( debug_code&DebugVerbosity::init )
      cout << " (x " << (e_ran/j->p4().E()) << " after " << count << " trials) --> ";

    // smear p4...
    j->setp4( (e_ran/j->p4().E())*j->p4() );
  }
  
  // if using external TFs, call directly transfer_function2 asking for the RND value (not the density)
  else if( external_tf ){

    double x[2] = {j->p4().Pt(), j->p4().Eta()};
    int accept{0};
    // calling with 'true' returns the randomised pT 
    double pt_ran = transfer_function2( j, x, jet_type , accept, cfg.tf_offscale, true, debug_code );

    int count{0};
    while( count<100 && pt_ran<MEM::TF_ACC_param[1] ){
      pt_ran = transfer_function2( j, x, jet_type , accept, cfg.tf_offscale, true, debug_code );
      ++count;
    }

    if( debug_code&DebugVerbosity::init )
      cout << " (x " << (pt_ran/j->p4().Pt()) << ", after " << count << " trials) --> ";

    // smear p4...
    j->setp4( (pt_ran/j->p4().Pt())*j->p4() ); 
  }
  else{ /*...*/ }
  
  if( debug_code&DebugVerbosity::init ){
    cout << "(" << j->p4().Px() << ", " << j->p4().Py() << ", " 
	 << j->p4().Pz() << ", " << j->p4().E() << ")" << endl;
  }
  
  return;
}


const TF1* MEM::Integrand::get_tf_global(TFType::TFType type, int etabin) const {
    std::pair<TFType::TFType, int> p = std::make_pair(type, etabin);
    if (tf_map.find(p) != tf_map.end()) {
        return &(tf_map.at(p));
    } else {
        std::cerr << "could not find tf for " << type << " " << etabin << std::endl;
        return nullptr;
    }
}
