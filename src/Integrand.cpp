#include "interface/Integrand.h"



MEM::Integrand::Integrand(MEM::Hypothesis hyp, int debug){
  hypo               = hyp;
  debug_code         = debug;
  num_of_vars        = 1;
  naive_jet_counting = 0;
  fs                 = FinalState::Undefined;
  ig2                = nullptr;

  if( debug_code&DebugVerbosity::init )  
    cout << "Integrand::Integrand()" << endl;
}

MEM::Integrand::~Integrand(){
  if( debug_code&DebugVerbosity::init )  
    cout << "Integrand::~Integrand()" << endl;    
  clear();
}

void MEM::Integrand::init(){

  CompPerm comparator;

  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::init(): Start" << endl;
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

  naive_jet_counting = 8*(fs==FinalState::HH) + 6*(fs==FinalState::LH) + 4*(fs==FinalState::LL);

  // deal with jets
  vector<int> perm_index;
  size_t n_jets = obs_jets.size();
  for(size_t id = 0; id < n_jets ; ++id) perm_index.push_back( id );
  while( perm_index.size() < naive_jet_counting)
  perm_index.push_back( -1 );
  // calculate upper / lower edges
  for( auto j : obs_jets ){
    double y[2] = { j->p4().E(), j->p4().Eta() };
    double* x;
    x = get_support( y, TFType::qReco ) ;
    j->addObs( Observable::E_LOW_Q,  x[0] );
    j->addObs( Observable::E_HIGH_Q, x[1] );
    x = get_support( y, TFType::bReco ) ;
    j->addObs( Observable::E_LOW_B,  x[0] );
    j->addObs( Observable::E_HIGH_B, x[1] );    
  }


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
  
  // formula to get the nuber of unknowns
  num_of_vars = 
    24                                   // dimension of the phas-space
    -3*obs_leptons.size()                // leptons
    -2*( TMath::Min(obs_jets.size(), naive_jet_counting) ) // jet directions
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

void MEM::Integrand::push_back_object(const LV& p4,  const MEM::ObjectType& type){

  Object* obj = new Object(p4, type);
  if( debug_code&DebugVerbosity::input ){
    cout << "Integrand::fill_map()" << endl;
    obj->print(cout);
  }

  switch( type ){
  case ObjectType::Jet:
    obs_jets.push_back( obj ); 
    break;
  case ObjectType::Lepton:
    obs_leptons.push_back( obj ); 
    break;
  case ObjectType::MET:
  case ObjectType::Recoil:
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

  prob += make_assumption( initializer_list<PSVar>{}  );
  //prob += make_assumption( initializer_list<PSVar>{PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //prob += make_assumption( initializer_list<PSVar>{PSVar::cos_b1, PSVar::phi_b1} );
  //prob += make_assumption( initializer_list<PSVar>{PSVar::cos_qbar1, PSVar::phi_qbar1,PSVar::cos_b1, PSVar::phi_b1} );
  
  cout << "Probability = " << prob << endl;


  if( debug_code&DebugVerbosity::init ){
    cout << "Integrand::run(): End" << endl;
  }
  delete ig2;

  clear();

  return;
}

void MEM::Integrand::clear(){
  for( auto j : obs_jets )     delete j;
  for( auto l : obs_leptons )  delete l;
  for( auto m : obs_mets )     delete m;
  obs_jets.clear();
  obs_leptons.clear();
  obs_mets.clear();
}

bool MEM::Integrand::test_assumption( const size_t& lost, size_t& extra_jets){

  if( (obs_jets.size()+lost) < naive_jet_counting){
    if( debug_code&DebugVerbosity::init ) cout << "\t This assumption cannot be made: too few jets" << endl;
    return false;
  }
    extra_jets = (obs_jets.size()+lost-naive_jet_counting);
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
  double prob{0.};

  size_t n_perm{0};
  for( auto perm : perm_indexes_assumption ){
    if(n_perm>0) break;
    double p{0.};
    p += probability(x, perm );
    prob += p;
    ++n_perm;
  }

  return prob;
}

void MEM::Integrand::evaluate_PS(MEM::PS& ps, const double* x, const vector<int>& perm ) const {

  switch( fs ){
  case FinalState::LH:
    return evaluate_PS_LH(ps, x, perm);
    break;
  case FinalState::LL:
    return evaluate_PS_LL(ps, x, perm);
    break;
  case FinalState::HH:
    return evaluate_PS_HH(ps, x, perm);
    break;
  default:
    break;
  }
}


void MEM::Integrand::evaluate_PS_LH(MEM::PS& ps, const double* x, const vector<int>& perm ) const {

  // jet,lepton counters
  size_t nj{0};
  size_t nl{0};

  // store temporary values to build four-vectors
  double E, cos, phi, sin;

  // t -> udb
  E       = x[ map_to_var.find(PSVar::E_q1)->second ];
  cos     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ]->p4().Theta()) : x[ map_to_var.find(PSVar::cos_q1)->second ];
  phi     = perm[nj]>=0 ? obs_jets[ perm[nj] ]->p4().Phi()               : x[ map_to_var.find(PSVar::phi_q1)->second ];
  sin     = sqrt(1-cos*cos);
  ps.set(PSPart::q1,  GenPart(TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*E, E ), 
			      perm[nj]>=0 ? TFType::qReco : TFType::qLost ) );    
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
  
  cos  = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ]->p4().Theta()) : x[ map_to_var.find(PSVar::cos_qbar1)->second ];
  phi  = perm[nj]>=0 ? obs_jets[ perm[nj] ]->p4().Phi()               : x[ map_to_var.find(PSVar::phi_qbar1)->second ];
  E    = solve_for_W_E_qbar(ps, 1);
  sin     = sqrt(1-cos*cos);
  ps.set(PSPart::qbar1, GenPart(  TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*E, E ),
				  perm[nj]>=0 ? TFType::qReco : TFType::qLost ));
  
  if( debug_code&DebugVerbosity::integration ){       
    cout << "\tFilling qbar1..." << endl;
    cout << "\t\tE_q1  = SOLVE" << endl;
    if( perm[nj]<0 ){
      cout << "\t\tcos_qbar1=x[" << map_to_var.find(PSVar::cos_qbar1)->second << "]" << endl;
      cout << "\t\tphi_qbar1=x[" << map_to_var.find(PSVar::phi_qbar1)->second << "]" << endl;
    }
    else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
  }
  ++nj;
  
  cos     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ]->p4().Theta()) : x[ map_to_var.find(PSVar::cos_b1)->second ];
  phi     = perm[nj]>=0 ? obs_jets[ perm[nj] ]->p4().Phi()               : x[ map_to_var.find(PSVar::phi_b1)->second ];
  E       = solve_for_T_E_b(ps, 1);
  sin     = sqrt(1-cos*cos);
  ps.set(PSPart::b1,  GenPart(  TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*sqrt(E*E - MB*MB), E ),
				perm[nj]>=0 ? TFType::bReco : TFType::bLost ));
  
  if( debug_code&DebugVerbosity::integration ){       
    cout << "\tFilling b1..." << endl;
    cout << "\t\tE_b1  = SOLVE" << endl;
    if( perm[nj]<0 ){
      cout << "\t\tcos_b1=x[" << map_to_var.find(PSVar::cos_b1)->second << "]" << endl;
      cout << "\t\tphi_b1=x[" << map_to_var.find(PSVar::phi_b1)->second << "]" << endl;
    }
    else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
  }
  ++nj;
  
  // t . lnb
  E       = obs_leptons[ nl]->p4().E();
  cos     = TMath::Cos(obs_leptons[ nl ]->p4().Theta());
  phi     = obs_leptons[ nl ]->p4().Phi();
  sin     = sqrt(1-cos*cos);
  ps.set(PSPart::q2,  GenPart(  TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*E, E ),
				TFType::muReco  ));
  
  if( debug_code&DebugVerbosity::integration ){       
    cout << "\tFilling q2..." << endl;
    cout << "\t\tusing obs_leptons[" << nl << "]" << endl;
  }    
  ++nl;
  
  cos  = x[ map_to_var.find(PSVar::cos_qbar2)->second ];
  phi  = x[ map_to_var.find(PSVar::phi_qbar2)->second ];
  E    = solve_for_W_E_qbar(ps, 2);
  sin     = sqrt(1-cos*cos);
  ps.set(PSPart::qbar2, GenPart(  TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*E, E ),
				  TFType::MET  ));
  
  if( debug_code&DebugVerbosity::integration ){ 
    cout << "\tFilling qbar2..." << endl;
    cout << "\t\tcos_qbar2=x[" << map_to_var.find(PSVar::cos_qbar2)->second << "]" << endl;
    cout << "\t\tphi_qbar2=x[" << map_to_var.find(PSVar::phi_qbar2)->second << "]" << endl;
  }
  
  cos     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ]->p4().Theta()) : x[ map_to_var.find(PSVar::cos_b2)->second ];
  phi     = perm[nj]>=0 ? obs_jets[ perm[nj] ]->p4().Phi()               : x[ map_to_var.find(PSVar::phi_b2)->second ];
  E       = solve_for_T_E_b(ps, 2);
  sin     = sqrt(1-cos*cos);
  ps.set(PSPart::b2,  GenPart(  TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*sqrt(E*E - MB*MB), E ),
				perm[nj]>=0 ? TFType::bReco : TFType::bLost  ));
  
  if( debug_code&DebugVerbosity::integration ){       
    cout << "\tFilling b2..." << endl;
    cout << "\t\tE_b2  = SOLVE" << endl;
    if( perm[nj]<0 ){
      cout << "\t\tcos_b2=x[" << map_to_var.find(PSVar::cos_b2)->second << "]" << endl;
      cout << "\t\tphi_b2=x[" << map_to_var.find(PSVar::phi_b2)->second << "]" << endl;
    }
    else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
  }
  ++nj;
  
  // H . bb
  E       = x[ map_to_var.find(PSVar::E_b)->second ];
  cos     = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ]->p4().Theta()) : x[ map_to_var.find(PSVar::cos_b)->second ];
  phi     = perm[nj]>=0 ? obs_jets[ perm[nj] ]->p4().Phi()               : x[ map_to_var.find(PSVar::phi_b)->second ];
  ps.set(PSPart::b, GenPart(  TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*sqrt(E*E - MB*MB), E ),
			      perm[nj]>=0 ? TFType::bReco : TFType::bLost  ));
  
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
  
  cos  = perm[nj]>=0 ? TMath::Cos(obs_jets[ perm[nj] ]->p4().Theta()) : x[ map_to_var.find(PSVar::cos_bbar)->second ];
  phi  = perm[nj]>=0 ? obs_jets[ perm[nj] ]->p4().Phi()               : x[ map_to_var.find(PSVar::phi_bbar)->second ];
  E    = hypo==Hypothesis::TTBB ? x[ map_to_var.find(PSVar::E_bbar)->second ]   : solve_for_H_E_bbar(ps, 1);
  ps.set(PSPart::bbar, GenPart(  TLorentzVector( TVector3( sin*TMath::Cos(phi),sin*TMath::Sin(phi),cos )*sqrt(E*E - MB*MB), E ),
				 perm[nj]>=0 ? TFType::bReco : TFType::bLost  ));  
  
  if( debug_code&DebugVerbosity::integration ){       
    cout << "\tFilling bbar..." << endl;
    if(hypo==Hypothesis::TTBB) cout << "\t\tE_bbar  = x[" <<  map_to_var.find(PSVar::E_bbar)->second << "]" << endl;
    else cout << "\t\tE_bbar  = SOLVE" << endl;
    if( perm[nj]<0 ){
      cout << "\t\tcos_bbar=x[" << map_to_var.find(PSVar::cos_bbar)->second << "]" << endl;
      cout << "\t\tphi_bbar=x[" << map_to_var.find(PSVar::phi_bbar)->second << "]" << endl;
    }
    else{ cout << "\t\tusing obs_jets[" << perm[nj] << "]" << endl;}
  }
  ++nj;  
}

void MEM::Integrand::evaluate_PS_LL(MEM::PS& ps, const double* x, const vector<int>& perm ) const {}

void MEM::Integrand::evaluate_PS_HH(MEM::PS& ps, const double* x, const vector<int>& perm ) const {}

double MEM::Integrand::probability(const double* x, const vector<int>& perm ) const {

  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::probability(): Start" << endl;
  }

  PS ps;
  evaluate_PS(ps, x, perm);
  if( debug_code&DebugVerbosity::integration ) ps.print(cout);

  double m   {1.}; // matrix element
  double w   {1.}; // transfer function
  double nu_x{0.}; // total nu's px
  double nu_y{0.}; // total nu's py
  size_t indx{0};  // quark counter 

  LV lv_q1    = ps.lv(PSPart::q1);
  LV lv_qbar1 = ps.lv(PSPart::qbar1);
  LV lv_b1    = ps.lv(PSPart::b1);
  LV lv_q2    = ps.lv(PSPart::q2);
  LV lv_qbar2 = ps.lv(PSPart::qbar2);
  LV lv_b2    = ps.lv(PSPart::b2);
  LV lv_b     = ps.lv(PSPart::b);
  LV lv_bbar  = ps.lv(PSPart::bbar);

  if( debug_code&DebugVerbosity::integration ){
    cout << "\tFilling m..." << endl;
  }

  m *= t_decay_amplitude(lv_q1, lv_qbar1, lv_b1);
  m *= t_decay_amplitude(lv_q2, lv_qbar2, lv_b2);
  m *= (hypo==Hypothesis::TTH ? H_decay_amplitude(lv_b, lv_bbar) : 1.0 );
  m *= scattering( lv_q1+lv_qbar1+lv_b1, lv_q2+lv_qbar2+lv_b2, lv_b, lv_bbar);
  m *= pdf( lv_q1+lv_qbar1+lv_b1+lv_q2+lv_qbar2+lv_b2+lv_b+lv_bbar );

  if( debug_code&DebugVerbosity::integration ){
    cout << "\tFilling p..." << endl;
  }

  map<MEM::PSPart, MEM::GenPart>::const_iterator p;
  for( p = ps.begin() ; p != ps.end() ; ++p ){
    if( isLepton  ( p->second.type ) ) continue;
    if( isNeutrino( p->second.type ) ){
      nu_x += p->second.lv.Px();
      nu_y += p->second.lv.Py();
      continue;
    }      
    // observables and generated-level quantities
    // if the parton is matched, test value of jet energy
    double e_gen  {p->second.lv.E()}; 
    double eta_gen{p->second.lv.Eta()};     
    double e_rec  = perm[indx]>=0 ? obs_jets[perm[indx]]->p4().E() : 0.;

    // build x,y vectors 
    double y[1] = { e_rec };
    double x[2] = { e_gen, eta_gen };
    w *= transfer_function( y, x, p->second.type ); 

    // increment quark counter
    ++indx;
  }

  if( fs!=FinalState::HH ){
    double y[2] = { obs_mets[0]->p4().Px(), obs_mets[0]->p4().Py() };
    double x[2] = { nu_x, nu_y };
    w *= transfer_function( y, x, TFType::MET );
  }

  if( debug_code&DebugVerbosity::integration ){
    cout << "\tIntegrand::probability(): End" << endl;
  }

  return m*w;
}


double MEM::Integrand::scattering(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&, const TLorentzVector&) const{
  double p{1.};
  return p;
}

double MEM::Integrand::pdf(const TLorentzVector& lv) const{
  double p{1.};
  return p;
}

double MEM::Integrand::t_decay_amplitude(const TLorentzVector&, const TLorentzVector&, const TLorentzVector&) const{
  double p{1.};
  return p;
}

double MEM::Integrand::H_decay_amplitude(const TLorentzVector&, const TLorentzVector&) const{
  double p{1.};
  return p;
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
