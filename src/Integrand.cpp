#include "interface/Integrand.h"

MEM::Integrand::Integrand(MEM::Hypothesis hyp, int debug){
  hypo        = hyp;
  debug_code  = debug;
  num_of_vars = 1;
  fs          = FinalState::Undefined;
  ig2         = nullptr;

  if( debug_code&DebugVerbosity::init )  
    cout << "Integrand::Integrand()" << endl;
}

MEM::Integrand::~Integrand(){
  if( debug_code&DebugVerbosity::init )  
    cout << "Integrand::~Integrand()" << endl;
}

void MEM::Integrand::init(){

  CompPerm comparator;

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::init(): Start" << endl;
  }  

  // deal with jets
  vector<int> perm_index;
  size_t n_jets = obs_jets.size();
  for(size_t id = 0; id < n_jets ; ++id) perm_index.push_back( id );
  while( perm_index.size()<6 ) perm_index.push_back( -1 );

  if( debug_code&DebugVerbosity::init ){
    cout << "\tIndexes to be permuted: [ " ;
    for( auto ind : perm_index )
      cout << ind << " ";
    cout << "]" << endl;
  }

  sort( perm_index.begin(), perm_index.end(), comparator);

  size_t n_perm{0};
  do{
    perm_indexes.push_back( perm_index );
    if( debug_code&DebugVerbosity::init ) {
      cout << "\tperm. " << n_perm << ": [ ";
      for( auto ind : perm_index )
	cout << ind << " ";
      cout << "]" << endl;
    }
    ++n_perm;
  } while( next_permutation( perm_index.begin(), perm_index.end(), comparator ) );
  if( debug_code&DebugVerbosity::init ){
    cout << "\tTotal of " << n_perm << " permutation(s) created" << endl;
  }
  
  // deal with leptons
  switch( obs_leptons.size() ){
  case 0:
    fs = FinalState::HH;
    break;
  case 1:
    fs = FinalState::LH;
    break;
  case 2:
    fs = FinalState::LL;
    break;
  default:
    cout << "*** MEM::Intgrator::init(): Unsupported final state" << endl;
    break;
  }
 
  // formula to get the nuber of unknowns
  num_of_vars = 
    24                                   // dimension of the phas-space
    -3*obs_leptons.size()                // leptons
    -2*( TMath::Max(obs_jets.size(), size_t(6)) ) // jet directions
    -4                                   // top/W mass
    -(hypo==Hypothesis::TTH);            // H mass
  if( debug_code&DebugVerbosity::init ){
    cout << "\tTotal of " << num_of_vars << " unknowns" << endl;
  }

  fill_map();

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::init(): End" << endl;
  }
 
  return;
}

void MEM::Integrand::get_xL(double* xL, const size_t nvar){

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::get_xL(): Start" << endl;
  }

  for(size_t i = 0 ; i < nvar ; ++i) xL[i]=0.;

  switch( fs ){
  case FinalState::LH:
    try{
      xL[map_to_var[PSVar::E_q1]]      =  0.;
      xL[map_to_var[PSVar::cos_qbar2]] = -1.;
      xL[map_to_var[PSVar::phi_qbar2]] = -TMath::Pi();
      xL[map_to_var[PSVar::E_b]]       =  0.;
    }
    catch(...){
      cout << "*** Integrand::get_xL(): Maybe out-of-range...check" << endl; 
    }
    break;
  default:
    break;
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::get_xL()" << endl;
    cout << "\tLower edges: [ ";
    for(size_t i = 0 ; i < nvar ; ++i)
      cout << xL[i] << " " ;
    cout << "]" << endl;
  }
 
}

void MEM::Integrand::get_xU(double* xU, const size_t nvar){

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::get_xU(): Start" << endl;
  }

  for(size_t i = 0 ; i < nvar ; ++i) xU[i]=0.;

  switch( fs ){
  case FinalState::LH:
    try{
      xU[map_to_var[PSVar::E_q1]]      =  1.;
      xU[map_to_var[PSVar::cos_qbar2]] = +1.;
      xU[map_to_var[PSVar::phi_qbar2]] = +TMath::Pi();
      xU[map_to_var[PSVar::E_b]]       =  1.;
    }
    catch(...){
      cout << "*** Integrand::get_xU(): Maybe out-of-range...check" << endl; 
    }
    break;
  default:
    break;
  }
 
  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::get_xU()" << endl;
    cout << "\tUpper edges: [ ";
    for(size_t i = 0 ; i < nvar ; ++i)
      cout << xU[i] << " " ;
    cout << "]" << endl;
  }
  
}

double MEM::Integrand::get_width(const double* xL, const double* xU, const size_t nvar){
  double out{1.};
  for(size_t i = 0; i < nvar ; ++i) out *= TMath::Abs( xU[i]-xL[i] );
  return out;
}

void MEM::Integrand::fill_map(){
  
  switch( fs ){
  case FinalState::LH:
    map_to_var[PSVar::E_q1]         = 0; //Eq 
    map_to_var[PSVar::cos_qbar2]    = 1; //cosNu
    map_to_var[PSVar::phi_qbar2]    = 2; //phiNu
    map_to_var[PSVar::E_b]          = 3; //Eb
    map_to_var[PSVar::E_bbar]       = hypo==Hypothesis::TTBB ? 4 : 23 ;   
    map_to_var[PSVar::cos_qbar1]    = 4  + (hypo==Hypothesis::TTBB); 
    map_to_var[PSVar::phi_qbar1]    = 5  + (hypo==Hypothesis::TTBB); 
    map_to_var[PSVar::cos_b1]       = 6  + (hypo==Hypothesis::TTBB); 
    map_to_var[PSVar::phi_b1]       = 7  + (hypo==Hypothesis::TTBB); 
    map_to_var[PSVar::cos_b2]       = 8  + (hypo==Hypothesis::TTBB); 
    map_to_var[PSVar::phi_b2]       = 9  + (hypo==Hypothesis::TTBB); 
    map_to_var[PSVar::cos_bbar]     = 10 + (hypo==Hypothesis::TTBB); 
    map_to_var[PSVar::phi_bbar]     = 11 + (hypo==Hypothesis::TTBB); 
    break;
  default:
    break;
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::fill_map()" << endl;
    for( auto iter = map_to_var.begin() ; iter!=map_to_var.end() ; ++iter ){
      cout << static_cast<int>(iter->first) << " <=> " << iter->second << endl;
    }
  }

  return;
}

void MEM::Integrand::push_back_object(const LV& p4,  const Object::Type type){

  MEM::Object obj;
  obj.init( p4, type );
  if( debug_code&DebugVerbosity::input ) obj.print(cout);

  switch( type ){
  case Object::Type::jet:
    obs_jets.push_back( obj ); 
    break;
  case Object::Type::lepton:
    obs_leptons.push_back( obj ); 
    break;
  case Object::Type::met:
  case Object::Type::recoil:
    obs_mets.push_back( obj ); 
    break;
  default:
    cout << "*** MEM::Intgrator::push_back_object(): Unknown type of object added" << endl;
    break;
  }

  return;
}

void MEM::Integrand::run(){

  ROOT::Math::Functor toIntegrate(this, &MEM::Integrand::Eval, num_of_vars);
  ig2 = new ROOT::Math::GSLMCIntegrator(ROOT::Math::IntegrationMultiDim::kVEGAS, 1.e-12, 1.e-5, 4);
  ig2->SetFunction(toIntegrate);

  double xL[num_of_vars];
  double xU[num_of_vars];
  get_xL(xL, num_of_vars);
  get_xU(xU, num_of_vars);
  cout << ig2->Integral(xL,xU)/get_width(xL,xU,num_of_vars) << endl;

  delete ig2;
  return;
}

double MEM::Integrand::Eval(const double* x) const{
  double prob{1.};

  for( auto perm : perm_indexes ){
    PS ps = evaluate_PS( x, perm);
    ps.print(cout);
  }

  return prob;
}

MEM::PS MEM::Integrand::evaluate_PS(const double* x, const vector<int>& perm ) const {

  PS ps;

  // jet,lepton counters
  size_t nj{0};
  size_t nl{0};

  switch( fs ){
  case FinalState::LH:
    // t . udb
    ps.E_q1       = x[ map_to_var.find(PSVar::E_q1)->second ];
    ps.cos_q1     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_q1)->second ];
    ps.phi_q1     = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_q1)->second ];
    ++nj;
    ps.cos_qbar1  = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_qbar1)->second ];
    ps.phi_qbar1  = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_qbar1)->second ];
    ps.E_qbar1    = solve_for_W_E_qbar(ps, 1);
    ++nj;
    ps.cos_b1     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_b1)->second ];
    ps.phi_b1     = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_b1)->second ];
    ps.E_b1       = solve_for_T_E_b(ps, 1);
    ++nj;

    // t . lnb
    ps.E_q2       = obs_jets[ perm[nl] ].p4.E();
    ps.cos_q2     = TMath::Cos(obs_jets[ perm[nl] ].p4.Theta());
    ps.phi_q2     = obs_jets[ perm[nl] ].p4.Phi();
    ++nl;
    ps.cos_qbar2  = x[ map_to_var.find(PSVar::cos_qbar2)->second ];
    ps.phi_qbar2  = x[ map_to_var.find(PSVar::phi_qbar2)->second ];
    ps.E_qbar2    = solve_for_W_E_qbar(ps, 2);
    ps.cos_b1     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_b2)->second ];
    ps.phi_b1     = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_b2)->second ];
    ps.E_b1       = solve_for_T_E_b(ps, 2);
    ++nj;

    // H . bb
    ps.E_b       = x[ map_to_var.find(PSVar::E_b)->second ];
    ps.cos_b     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_b)->second ];
    ps.phi_b     = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_b)->second ];
    ++nj;
    ps.cos_bbar  = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_bbar)->second ];
    ps.phi_bbar  = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_bbar)->second ];
    ps.E_bbar    = hypo==Hypothesis::TTBB ? x[ map_to_var.find(PSVar::E_bbar)->second ]   : solve_for_H_E_bbar(ps, 1);
    ++nj;
    break;
  default:
    break;
    }

  return ps;
}

double MEM::Integrand::solve_for_W_E_qbar(const MEM::PS& ps, const size_t n) const {
  double out{1.};
  return out;
}

double MEM::Integrand::solve_for_T_E_b(const MEM::PS& ps, const size_t n) const {
  double out{1.};
  return out;
}

double MEM::Integrand::solve_for_H_E_bbar(const MEM::PS& ps, const size_t n) const {
  double out{1.};
  return out;
}
