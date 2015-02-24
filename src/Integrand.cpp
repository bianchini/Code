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
    if( debug_code&DebugVerbosity::init&0 ) {
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
    -2*( TMath::Min(obs_jets.size(), size_t(6)) ) // jet directions
    -4                                   // top/W mass
    -1*(hypo==Hypothesis::TTH);          // H mass
  if( debug_code&DebugVerbosity::init ){
    cout << "\tTotal of " << num_of_vars << " unknowns" << endl;
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::init(): End" << endl;
  }
 
  return;
}

void MEM::Integrand::get_xL(double* xL, const initializer_list<PSVar>& lost, const size_t nvar){

  for(size_t i = 0 ; i < nvar ; ++i) xL[i]=-99.;

  switch( fs ){
  case FinalState::LH:
    try{
      xL[map_to_var[PSVar::E_q1]]      =  0.;
      xL[map_to_var[PSVar::cos_qbar2]] = -1.;
      xL[map_to_var[PSVar::phi_qbar2]] = -TMath::Pi();
      xL[map_to_var[PSVar::E_b]]       =  0.;
      size_t count_extra{0};
      for( auto l = lost.begin() ; l !=lost.end(); ++l ){
	xL[map_to_var[*l]] = count_extra%2==0 ? -1 : -TMath::Pi(); 
	++count_extra;
      }      
    }
    catch(...){
      cout << "*** Integrand::get_xL(): Maybe out-of-range...check" << endl; 
    }
    break;
  default:
    break;
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::get_xL(): Lower edges: [ ";
    for(size_t i = 0 ; i < nvar ; ++i)
      cout << xL[i] << " " ;
    cout << "]" << endl;
  }
 
}

void MEM::Integrand::get_xU(double* xU, const initializer_list<PSVar>& lost, const size_t nvar){

  for(size_t i = 0 ; i < nvar ; ++i) xU[i]=-99.;

  switch( fs ){
  case FinalState::LH:
    try{
      xU[map_to_var[PSVar::E_q1]]      =  1.;
      xU[map_to_var[PSVar::cos_qbar2]] = +1.;
      xU[map_to_var[PSVar::phi_qbar2]] = +TMath::Pi();
      xU[map_to_var[PSVar::E_b]]       =  1.;
      size_t count_extra{0};
      for( auto l = lost.begin() ; l!=lost.end() ; ++l ){
	xU[map_to_var[*l]] = count_extra%2==0 ? +1 : +TMath::Pi(); 
	++count_extra;
      }      
    }
    catch(...){
      cout << "*** Integrand::get_xU(): Maybe out-of-range...check" << endl; 
    }
    break;
  default:
    break;
  }
 
  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::get_xU(): Upper edges: [ ";
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

void MEM::Integrand::fill_map(const initializer_list<PSVar>& lost){
  
  size_t count_extra{0};
  switch( fs ){
  case FinalState::LH:
    map_to_var[PSVar::E_q1]         = 0; //Eq 
    map_to_var[PSVar::cos_qbar2]    = 1; //cosNu
    map_to_var[PSVar::phi_qbar2]    = 2; //phiNu
    map_to_var[PSVar::E_b]          = 3; //Eb
    if( hypo==Hypothesis::TTBB ) map_to_var[PSVar::E_bbar] = 4;
    count_extra = 4+(hypo==Hypothesis::TTBB);
    for( auto l = lost.begin() ; l!=lost.end() ; ++l ) map_to_var[*l] = count_extra++;
    break;
  case FinalState::LL:
  case FinalState::HH:
  default:
    break;
  }

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::fill_map()" << endl;
    for( auto iter = map_to_var.begin() ; iter!=map_to_var.end() ; ++iter ){
      cout << "PS[" << static_cast<int>(iter->first) << "] maps to x[" << iter->second << "]" << endl;
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
  
  double prob{0.};

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::run(): Start" << endl;
  }

  ig2 = new ROOT::Math::GSLMCIntegrator(ROOT::Math::IntegrationMultiDim::kVEGAS, 1.e-12, 1.e-5, 10);

  //prob += make_assumption( initializer_list<PSVar>{},  );
  //prob += make_assumption( initializer_list<PSVar>{PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //prob += make_assumption( initializer_list<PSVar>{PSVar::cos_b1, PSVar::phi_b1} );
  prob += make_assumption( initializer_list<PSVar>{PSVar::cos_qbar1, PSVar::phi_qbar1,PSVar::cos_b1, PSVar::phi_b1} );
  
  cout << "Probability = " << prob << endl;


  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::run(): End" << endl;
  }
  delete ig2;
  return;
}

bool MEM::Integrand::test_assumption( const size_t& lost, size_t& extra_jets){

  if( (fs==FinalState::LH && obs_jets.size()+lost<6) ||
      (fs==FinalState::LL && obs_jets.size()+lost<4) ||
      (fs==FinalState::HH && obs_jets.size()+lost<8) ){
    if( debug_code&DebugVerbosity::init )
      cout << "\t This assumption cannot be made: too few jets" << endl;       
    return false;
  }
  if( fs==FinalState::LH ) extra_jets = (obs_jets.size()+lost-6);
  if( fs==FinalState::LL ) extra_jets = (obs_jets.size()+lost-4);
  if( fs==FinalState::HH ) extra_jets = (obs_jets.size()+lost-8);

  return true;
}

double MEM::Integrand::make_assumption( initializer_list<MEM::PSVar>&& lost){

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::make_assumption(): Start" << endl;
  }

  double prob{0.};
  size_t extra_jets{0};
  if(!test_assumption(lost.size()/2, extra_jets)) return prob;

  perm_indexes_assumption.clear();

  // extra variables to integrate over  
  fill_map( lost );  

  // remove unwanted permutations
  for( auto perm : perm_indexes ){    
    auto good_perm = perm;
    size_t count{0};
    for(auto it = lost.begin() ; it!=lost.end() ; ++count, ++it ){
      if(count%2==0) perm[ (static_cast<size_t>(*it)-1)/3 ] = -1;      
    }
    count = 0;
    for(auto ind : perm){if(ind<0) ++count;}
    if(count==(lost.size()/2))
      perm_indexes_assumption.push_back( perm );    
  }  

  if( debug_code&DebugVerbosity::init ) {
    cout << "\tAssumption: Total of " << perm_indexes_assumption.size() << " considered" << endl;
    size_t n_perm{0};
    for( auto perm : perm_indexes_assumption ){
      cout << "\tperm. " << n_perm++ << ": [ ";
      for( auto ind : perm )
	cout << ind << " ";
      cout << "]" << endl;
    }
  }


  // create integration ranges
  size_t npar = num_of_vars+2*extra_jets;

  double xL[npar];
  double xU[npar];
  get_xL(xL, lost, npar);
  get_xU(xU, lost, npar);      
  
  // function
  ROOT::Math::Functor toIntegrate(this, &MEM::Integrand::Eval, npar);  
  ig2->SetFunction(toIntegrate);
  
  // do the integral
  prob = ig2->Integral(xL,xU)/get_width(xL,xU,npar) ;
  
  return prob;

}

double MEM::Integrand::Eval(const double* x) const{
  double prob{1.};

  size_t n_perm{0};
  for( auto perm : perm_indexes_assumption ){
    if(n_perm>0) break;
    PS ps = evaluate_PS( x, perm );
    //ps.print(cout);
    ++n_perm;
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
    if( debug_code&DebugVerbosity::integration ){       
      cout << "\tFilling q1..." << endl;
      cout << "\t\tE_q1  =x[" <<  map_to_var.find(PSVar::E_q1)->second << "]" << endl;
      if( perm[nj]<0 ){
	cout << "\t\tcos_q1=x[" << map_to_var.find(PSVar::cos_q1)->second << "]" << endl;
	cout << "\t\tphi_q1=x[" << map_to_var.find(PSVar::phi_q1)->second << "]" << endl;
      }
      else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
    }
    ++nj;

    ps.cos_qbar1  = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_qbar1)->second ];
    ps.phi_qbar1  = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_qbar1)->second ];
    ps.E_qbar1    = solve_for_W_E_qbar(ps, 1);
    if( debug_code&DebugVerbosity::integration ){       
      cout << "\tFilling qbar1..." << endl;
      cout << "\t\tE_q1  = solve" << endl;
      if( perm[nj]<0 ){
	cout << "\t\tcos_qbar1=x[" << map_to_var.find(PSVar::cos_qbar1)->second << "]" << endl;
	cout << "\t\tphi_qbar1=x[" << map_to_var.find(PSVar::phi_qbar1)->second << "]" << endl;
      }
      else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
    }
    ++nj;

    ps.cos_b1     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_b1)->second ];
    ps.phi_b1     = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_b1)->second ];
    ps.E_b1       = solve_for_T_E_b(ps, 1);
    if( debug_code&DebugVerbosity::integration ){       
      cout << "\tFilling b1..." << endl;
      cout << "\t\tE_b1  = solve" << endl;
      if( perm[nj]<0 ){
	cout << "\t\tcos_b1=x[" << map_to_var.find(PSVar::cos_b1)->second << "]" << endl;
	cout << "\t\tphi_b1=x[" << map_to_var.find(PSVar::phi_b1)->second << "]" << endl;
      }
      else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
    }
    ++nj;

    // t . lnb
    ps.E_q2       = obs_leptons[ nl].p4.E();
    ps.cos_q2     = TMath::Cos(obs_leptons[ nl ].p4.Theta());
    ps.phi_q2     = obs_leptons[ nl ].p4.Phi();
    if( debug_code&DebugVerbosity::integration ){       
      cout << "\tFilling q2..." << endl;
      cout << "\t\tusing obs_leptons[" << nl << "]" << endl;
    }    
    ++nl;

    ps.cos_qbar2  = x[ map_to_var.find(PSVar::cos_qbar2)->second ];
    ps.phi_qbar2  = x[ map_to_var.find(PSVar::phi_qbar2)->second ];
    ps.E_qbar2    = solve_for_W_E_qbar(ps, 2);
    if( debug_code&DebugVerbosity::integration ){ 
      cout << "\tFilling qbar2..." << endl;
	cout << "\t\tcos_qbar2=x[" << map_to_var.find(PSVar::cos_qbar2)->second << "]" << endl;
	cout << "\t\tphi_qbar2=x[" << map_to_var.find(PSVar::phi_qbar2)->second << "]" << endl;
    }

    ps.cos_b2     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_b2)->second ];
    ps.phi_b2     = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_b2)->second ];
    ps.E_b2       = solve_for_T_E_b(ps, 2);
    if( debug_code&DebugVerbosity::integration ){       
      cout << "\tFilling b2..." << endl;
      cout << "\t\tE_b2  = solve" << endl;
      if( perm[nj]<0 ){
	cout << "\t\tcos_b2=x[" << map_to_var.find(PSVar::cos_b2)->second << "]" << endl;
	cout << "\t\tphi_b2=x[" << map_to_var.find(PSVar::phi_b2)->second << "]" << endl;
      }
      else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
    }
    ++nj;

    // H . bb
    ps.E_b       = x[ map_to_var.find(PSVar::E_b)->second ];
    ps.cos_b     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_b)->second ];
    ps.phi_b     = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_b)->second ];
    if( debug_code&DebugVerbosity::integration ){       
      cout << "\tFilling b..." << endl;
      cout << "\t\tE_b  = x[" << map_to_var.find(PSVar::E_b)->second << "]" << endl;
      if( perm[nj]<0 ){
	cout << "\t\tcos_b=x[" << map_to_var.find(PSVar::cos_b)->second << "]" << endl;
	cout << "\t\tphi_b=x[" << map_to_var.find(PSVar::phi_b)->second << "]" << endl;
      }
      else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
    }
    ++nj;

    ps.cos_bbar  = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ].p4.Theta()) : x[ map_to_var.find(PSVar::cos_bbar)->second ];
    ps.phi_bbar  = perm[nj]>=0 ? obs_jets[ perm[nj] ].p4.Phi()               : x[ map_to_var.find(PSVar::phi_bbar)->second ];
    ps.E_bbar    = hypo==Hypothesis::TTBB ? x[ map_to_var.find(PSVar::E_bbar)->second ]   : solve_for_H_E_bbar(ps, 1);
    if( debug_code&DebugVerbosity::integration ){       
      cout << "\tFilling bbar..." << endl;
      if(hypo==Hypothesis::TTBB) cout << "\t\tE_bbar  = x[" <<  map_to_var.find(PSVar::E_bbar)->second << "]" << endl;
      else cout << "\t\tE_bbar  = solve" << endl;
      if( perm[nj]<0 ){
	cout << "\t\tcos_bbar=x[" << map_to_var.find(PSVar::cos_bbar)->second << "]" << endl;
	cout << "\t\tphi_bbar=x[" << map_to_var.find(PSVar::phi_bbar)->second << "]" << endl;
      }
      else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
    }
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
