#include "interface/HypoTester.h"


Algo::HypoTester::HypoTester(TTree* t){

  // reset variables
  nParam_j          = 0;
  nParam_n          = 0;
  nParam_m          = 0;
  count_hypo        = 0;
  count_perm        = 0;
  count_TopHad      = 0;
  count_TopHadLost  = 0;
  count_TopLep      = 0;
  count_Higgs       = 0;
  count_Radiation_u = 0;
  count_Radiation_d = 0;
  count_Radiation_c = 0;
  count_Radiation_b = 0;
  count_Radiation_g = 0;
  invisible         = 0;
  verbose           = 0;
  error_code        = 0;
  minimizer         = nullptr;

  if(t!=nullptr){
    cout << "Algo::HypoTester::HypoTester(): Creating branches to output file" << endl;
    event  = new Event(t);
    event->createBranches();
    event->reset();
  }
  else{
    event = nullptr;
  }
}

Algo::HypoTester::~HypoTester(){
  cout << "Algo::HypoTester::~HypoTester(): Removing HypoTester" << endl;
  for( auto perm : permutations ) delete perm;
  delete event;
}


void Algo::HypoTester::push_back_object( const LV& p4, char type ){

  Algo::Object obj;
  obj.init( p4, type );
  
  switch( type ){
  case 'j':
    p4_Jet.push_back( obj ); 
    break;
  case 'l':
    p4_Lepton.push_back( obj ); 
    break;
  case 'm':
    p4_MET.push_back( obj ); 
    break;
  default:
    cout << "Algo::HypoTester::push_back_object(): Unknown type of object added" << endl;
    break;
  }

  if(event!=nullptr){
    for(auto jet : p4_Jet){
      double btag = (jet.obs).find("BTAG")!=(jet.obs).end() ? (jet.obs).find("BTAG")->second : 0.; 
      if(btag>0.5) ++(event->treeStruct.n_btag);
      ++(event->treeStruct.n_jet);
    }
    for(auto lep : p4_Lepton){ 
      ++(event->treeStruct.n_lep); 
    }
  }
	
}

void Algo::HypoTester::add_object_observables( const string& name, const double& val, const char type){

  switch( type ){
  case 'j':
    if(p4_Jet.size()>0) (p4_Jet.back()).addObs( name, val );
    break;
  case 'l':
    if(p4_Lepton.size()>0) (p4_Lepton.back()).addObs( name, val );
    break;
  case 'm':
    if(p4_MET.size()>0) (p4_MET.back()).addObs( name, val );
    break;
  default:
    cout << "Algo::HypoTester::add_object_observables(): Unknown type of object added" << endl;
    break;
  }

}


void Algo::HypoTester::test( const map<string, vector<Decay::Decay>>& all ){

  // benchmark time: start clock for global timing
  auto t0 = high_resolution_clock::now();

  if(event!=nullptr) event->reset();

  // save event-dependent variables 
  save_global_variables();

  // create minimizer once for all hypotheses
  // the function is however set inside the run() block
  minimizer =  ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

  // make sure we start from zero
  count_hypo = 0;
  for( auto hypo : all ){
    if(verbose>0)
      cout << "Algo::HypoTester::test(): Testing hypothesis " << count_hypo 
	   << " with name \"" << hypo.first << "\"" << endl;
    reset();
    for( auto decay : hypo.second ){
      assume( decay );
    }
    try{    
      init();
    }
    catch(...){
      cout << "Algo::HypoTester::init() has thrown an exception:"
	" cannot compute all permutations. Return" << endl;      
      error_code = 1;
      next_event();
      return;
    }
    run();
    ++count_hypo;
  }
  
  // stop clock
  auto t1 = high_resolution_clock::now();
  if(event!=nullptr){
    int time_0 = static_cast<int>(duration_cast<milliseconds>(t1-t0).count());
    event->treeStruct.all_time = time_0; 
    event->fillTree();
    if(verbose>0) event->printTree();
  }
  next_event();
}

void Algo::HypoTester::reset(){

  // remove old permutations
  // perm calls destructor of block
  for( auto perm : permutations ) delete perm;
  
  decays.clear();
  particles.clear();
  permutations.clear();
  nParam_j          = 0;
  nParam_n          = 0;
  nParam_m          = 0;
  count_perm        = 0;
  count_TopHad      = 0;
  count_TopHadLost  = 0;
  count_WHad        = 0;
  count_TopLep      = 0;
  count_Higgs       = 0;
  count_Radiation_u = 0;
  count_Radiation_d = 0;
  count_Radiation_c = 0;
  count_Radiation_b = 0;
  count_Radiation_g = 0;
  invisible         = 0;

  // reset minimzer for recursive minimizations
  if(minimizer!=nullptr) minimizer->Clear();

}

void Algo::HypoTester::next_event(){

  reset();
  p4_Jet.clear();
  p4_Lepton.clear();
  p4_MET.clear();

  // delete the minimizer: no longer needed
  delete minimizer;
}

void Algo::HypoTester::assume( Decay::Decay decay ){
  decays.push_back( decay );
}


void Algo::HypoTester::unpack_assumptions(){

  if(verbose>0){ cout << "Algo::HypoTester::unpack_assumptions()" << endl; }

  // reset
  nParam_j          = 0;
  nParam_n          = 0;
  nParam_m          = 0;
  count_TopHad      = 0;
  count_TopHadLost  = 0;
  count_WHad        = 0;
  count_TopLep      = 0;
  count_Higgs       = 0;
  count_Radiation_u = 0;
  count_Radiation_d = 0;
  count_Radiation_c = 0;
  count_Radiation_b = 0;
  count_Radiation_g = 0;

  for( auto decay : decays ){

   switch( decay ){

   case Algo::Decay::Decay::TopHad:
     particles.push_back( make_pair( FinalState::FinalState::TopHad_q,    count_TopHad) );
     particles.push_back( make_pair( FinalState::FinalState::TopHad_qbar, count_TopHad) );
     particles.push_back( make_pair( FinalState::FinalState::TopHad_b,    count_TopHad) );
     ++count_TopHad;
     nParam_j += 3;
     if(verbose>0){ cout << "\tAdded TopHad" << endl; }
     break;

   case Algo::Decay::Decay::TopHadLost:
     particles.push_back( make_pair( FinalState::FinalState::TopHad_q,        count_TopHadLost) );
     particles.push_back( make_pair( FinalState::FinalState::TopHad_b,        count_TopHadLost) );
     ++count_TopHadLost;
     nParam_j += 2;
     nParam_m += 2;
     if(verbose>0){ cout << "\tAdded TopHadLost" << endl; }
     break;

   case Algo::Decay::Decay::WHad:
     particles.push_back( make_pair( FinalState::FinalState::WHad_q,    count_WHad) );
     particles.push_back( make_pair( FinalState::FinalState::WHad_qbar, count_WHad) );
     ++count_WHad;
     nParam_j += 2;
     if(verbose>0){ cout << "\tAdded WHad" << endl; }
     break;

   case Algo::Decay::Decay::TopLep:
     particles.push_back( make_pair( FinalState::FinalState::TopLep_b,    count_TopLep) );
     ++count_TopLep;
     ++invisible;
     nParam_j += 1;
     nParam_n += 2;
     if(verbose>0){ cout << "\tAdded TopLep" << endl; }
     break;

   case Algo::Decay::Decay::Higgs:
     particles.push_back( make_pair( FinalState::FinalState::Higgs_b,    count_Higgs) );
     particles.push_back( make_pair( FinalState::FinalState::Higgs_bbar, count_Higgs) );
     ++count_Higgs;
     nParam_j += 2;
     if(verbose>0){ cout << "\tAdded Higgs" << endl; }
     break;

   case Algo::Decay::Decay::Radiation_u:
     particles.push_back( make_pair( FinalState::FinalState::Radiation_u,   count_Radiation_u) );
     if(verbose>0){ cout << "\tAdded Radiation (u)" << endl; }
     ++count_Radiation_u;
     nParam_j += 1;
     break;

   case Algo::Decay::Decay::Radiation_d:
     particles.push_back( make_pair( FinalState::FinalState::Radiation_d,   count_Radiation_d) );
     if(verbose>0){ cout << "\tAdded Radiation (d)" << endl; }
     ++count_Radiation_d;
     nParam_j += 1;
     break;

   case Algo::Decay::Decay::Radiation_c:
     particles.push_back( make_pair( FinalState::FinalState::Radiation_c,   count_Radiation_c) );
     if(verbose>0){ cout << "\tAdded Radiation (c)" << endl; }
     ++count_Radiation_c;
     nParam_j += 1;
     break;

   case Algo::Decay::Decay::Radiation_b:
     particles.push_back( make_pair( FinalState::FinalState::Radiation_b,   count_Radiation_b) );
     if(verbose>0){ cout << "\tAdded Radiation (b)" << endl; }
     ++count_Radiation_b;
     nParam_j += 1;
     break;

   case Algo::Decay::Decay::Radiation_g:
     particles.push_back( make_pair( FinalState::FinalState::Radiation_g,   count_Radiation_g) );
     if(verbose>0){ cout << "\tAdded Radiation (g)" << endl; }
     ++count_Radiation_g;
     nParam_j += 1;
     break;

   default:     
     break;

   }    

  }

  // complement with Radiation
  if( particles.size() < p4_Jet.size() ){

    size_t addradiation = p4_Jet.size()-particles.size();
    while( addradiation > 0){
      particles.push_back( make_pair( FinalState::FinalState::Radiation_g,   count_Radiation_g) );
      if(verbose>0){ cout << "\tAdded Radiation (g)" << endl; }
      ++count_Radiation_g;
      nParam_j += 1;
      --addradiation;
    }

  }
    
  if(verbose>0){ cout << "\tTotal number of parameters: nParam_j=" << nParam_j << ", nParam_n=" << nParam_n << ", nParam_m=" << nParam_m << endl; }

}


void Algo::HypoTester::init(){

  if(verbose>0){ cout << "Algo::HypoTester::init()" << endl; }

  // first, unpack assumptions
  unpack_assumptions();

  // number of jets at least as large as num of quarks
  assert( particles.size()<=p4_Jet.size() );

  // number of leptons at least as large as num of leptonic tops
  assert( count_TopLep<=p4_Lepton.size()  );
  if( count_TopLep>0 && count_TopLep<p4_Lepton.size() ) 
    cout << "Algo::HypoTester::init(): Warning: only the first " 
	 << count_TopLep << " leptons will be used (no permutations among leptons)" << endl;

  // if >0 topLep, need MET
  assert( count_TopLep==0 || (count_TopLep>0 && p4_MET.size()>0)  );
  
  // order particles
  sort( particles.begin(), particles.end(), MyComp );

  // permutations visited
  count_perm = 0;

  // bookkeep to remove unnecessary permutations
  vector<vector<std::pair<FinalState::FinalState,size_t>>> logbook;

  // permutations using std library
  do {            
    if(verbose>2){ cout << "\tPermutation number " << count_perm << ":" << endl; }    
    if( go_to_next(logbook) ) continue;
    permutations.push_back( new Algo::CombBuilder( group_particles() , verbose) );    
    ++count_perm;    
  } while ( next_permutation(particles.begin(), particles.end(), MyComp  ) );

  if(verbose>0) {
    cout << "\tTotal of  " << count_perm << " permutations created:" << endl;
    size_t count {0};
    for( auto perm : permutations){
      cout << "(" << count << ")";
      perm->print(cout);
      ++count;
    }
  }
  
}

vector<Algo::DecayBuilder*> Algo::HypoTester::group_particles(){
  
  if(verbose>2){ cout << "Algo::HypoTester::group_particles()" << endl; }
  
  vector<Algo::DecayBuilder*> decayed;
  
  //  get all hadronically decaying tops
  for( size_t t_had = 0; t_had < count_TopHad; ++t_had ){
    if(verbose>1) cout << "\tProcessing " << t_had << "th TopHad" << endl;
    Algo::TopHadBuilder* topHad = new Algo::TopHadBuilder(verbose);    
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == t_had )  
	topHad->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( topHad );
  }

  //  get all leptonically decaying tops
  for( size_t t_had_lost = 0; t_had_lost < count_TopHadLost; ++t_had_lost ){
    if(verbose>1) cout << "\tProcessing " << t_had_lost << "th TopHadLost" << endl;
    Algo::TopHadBuilder* topHadLost = new Algo::TopHadBuilder(verbose);
    Object dummy;
    dummy.init( LV(0.,0.,0.,0.), 'j');
    topHadLost->init( FinalState::FinalState::TopHadLost_qbar , dummy , 2*t_had_lost + nParam_j + nParam_n ); // offset
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == t_had_lost ){
	  topHadLost->init( part.first, p4_Jet[pos] , pos );
      }
      ++pos;
    }
    decayed.push_back( topHadLost );
  }

  // get all hadronically decaying tops
  for( size_t w_had = 0; w_had < count_WHad; ++w_had ){
    if(verbose>1) cout << "\tProcessing " << w_had << "th WHad" << endl;
    Algo::WHadBuilder* wHad = new Algo::WHadBuilder(verbose);    
    size_t pos = 0;
    for( auto part : particles ){      
      if( part.second == w_had ) 
	wHad->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( wHad );
  }

  //  get all leptonically decaying tops
  for( size_t t_lep = 0; t_lep < count_TopLep; ++t_lep ){
    if(verbose>1) cout << "\tProcessing " << t_lep << "th TopLep" << endl;
    Algo::TopLepBuilder* topLep = new Algo::TopLepBuilder(verbose);   
    topLep->init( FinalState::FinalState::TopLep_l , p4_Lepton[t_lep] , 2*t_lep + nParam_j ); // offset
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == t_lep ) 
	topLep->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( topLep );
  }

  //  get all radiation
  for( size_t r_had_u = 0; r_had_u < count_Radiation_u; ++r_had_u ){
    if(verbose>1) cout << "\tProcessing " << r_had_u << "th Radiaton (u)" << endl;
    Algo::RadiationBuilder* rad = new Algo::RadiationBuilder(verbose, Decay::Decay::Radiation_u);    
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == r_had_u && part.first==FinalState::FinalState::Radiation_u )  
	rad->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( rad );
  }

  for( size_t r_had_d = 0; r_had_d < count_Radiation_d; ++r_had_d ){
    if(verbose>1) cout << "\tProcessing " << r_had_d << "th Radiaton (d)" << endl;
    Algo::RadiationBuilder* rad = new Algo::RadiationBuilder(verbose, Decay::Decay::Radiation_d);
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == r_had_d && part.first==FinalState::FinalState::Radiation_d )
        rad->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( rad );
  }

  for( size_t r_had_c = 0; r_had_c < count_Radiation_c; ++r_had_c ){
    if(verbose>1) cout << "\tProcessing " << r_had_c << "th Radiaton (c)" << endl;
    Algo::RadiationBuilder* rad = new Algo::RadiationBuilder(verbose, Decay::Decay::Radiation_c);
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == r_had_c && part.first==FinalState::FinalState::Radiation_c )
        rad->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( rad );
  }
  
  //  get all radiation                                                                                                                  
  for( size_t r_had_b = 0; r_had_b < count_Radiation_b; ++r_had_b ){
    if(verbose>1) cout << "\tProcessing " << r_had_b << "th Radiaton (b)" << endl;
    Algo::RadiationBuilder* rad = new Algo::RadiationBuilder(verbose, Decay::Decay::Radiation_b);
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == r_had_b && part.first==FinalState::FinalState::Radiation_b )
        rad->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( rad );
  }

  for( size_t r_had_g = 0; r_had_g < count_Radiation_g; ++r_had_g ){
    if(verbose>1) cout << "\tProcessing " << r_had_g << "th Radiaton (g)" << endl;
    Algo::RadiationBuilder* rad = new Algo::RadiationBuilder(verbose, Decay::Decay::Radiation_g);
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second == r_had_g && part.first==FinalState::FinalState::Radiation_g )
        rad->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( rad );
  }

  // get all hadronically decaying higgs
  for( size_t higgs = 0; higgs < count_Higgs; ++higgs ){
    if(verbose>1) cout << "\tProcessing " << higgs << "th Higgs" << endl;
    Algo::HiggsBuilder* hig = new Algo::HiggsBuilder(verbose);   
    size_t pos = 0;
    for( auto part : particles ){      
      if( part.second == higgs ) 
	hig->init( part.first, p4_Jet[pos] , pos );
      ++pos;
    }
    decayed.push_back( hig );
  }

  /* Add other blocks here */

  // if there are invisible particles, add MET as last element                                                                            
  if(invisible>0){
    Algo::METBuilder* met = new Algo::METBuilder(verbose);
    assert( p4_MET.size()>0 );
    met->init( p4_MET[0] );
    decayed.push_back( met );
  }
  // if there are not invisible particles, but MET is an input, eval MET tf at (0.,0.)                                                      
  else if( p4_MET.size()>0 ){
    Algo::METBuilder* met = new Algo::METBuilder(verbose);
    met->init( p4_MET[0] );
    met->fix_vars();
    decayed.push_back( met );
  }
  else{ /* do nothing */ }

  return decayed;

}


void Algo::HypoTester::setup_minimizer( const Algo::Strategy str){

  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(0);

  // initial values
  double step_E  {1.0};
  double step_Phi{0.2};
  double step_Cos{0.2};

  // count parameters
  size_t count_param{0};

  switch( str ){
  case Strategy::StartFromLastMinimum:    
    cout << "\t>>> Start from last minimum";    
  case Strategy::FirstTrial:    

    for( size_t p = 0 ; p < nParam_j ; p++){
      string name = "jet"+to_string(p);
      double inVal    {p4_Jet[p].p4.E()};
      double inVal_lo {inVal/3.0};
      double inVal_hi {inVal*3.0};      

      if( !is_variable_used(p) ){
	minimizer->SetFixedVariable(p, name.c_str(), inVal);
	if(verbose>0)
	  printf("\tParam[%lu] = %s set to FIXED value %.0f\n", p,name.c_str(), inVal);      
      }
      else{
	minimizer->SetLimitedVariable(p, name.c_str(), inVal  , step_E, inVal_lo , inVal_hi);
	if(verbose>0)
	  printf("\tParam[%lu] = %s set to %.0f. Range: [%.0f,%.0f]\n", p,name.c_str(), inVal, inVal_lo, inVal_hi );      
      }

      if(event!=nullptr){
	double btag = (p4_Jet[p].obs).find("BTAG")!=(p4_Jet[p].obs).end() ?  (p4_Jet[p].obs).find("BTAG")->second : 0.; 
	event->treeStruct.obs_e   [ event->treeStruct.n_dim + count_param ] = inVal;     
	event->treeStruct.obs_pt  [ event->treeStruct.n_dim + count_param ] = p4_Jet[p].p4.Pt();     
	event->treeStruct.obs_eta [ event->treeStruct.n_dim + count_param ] = p4_Jet[p].p4.Eta();     
	event->treeStruct.obs_phi [ event->treeStruct.n_dim + count_param ] = p4_Jet[p].p4.Phi();     
	event->treeStruct.obs_btag[ event->treeStruct.n_dim + count_param ] = btag;
      }
      ++count_param;
    }

    for( size_t p = 0 ; p < nParam_n ; p++){

      string name;
      double inVal, inVal_lo, inVal_hi, step;
      if( p%2==0 ){
	name = "phi"+to_string(nParam_j+p);
	inVal    = p4_MET[0].p4.Phi();
	inVal_lo = -TMath::Pi() ;
	inVal_hi = +TMath::Pi() ;
	step     = step_Phi;     
      }
      else{
	name = "cos"+to_string(nParam_j+p/2);
	inVal    = 0.;
	inVal_lo = -1.;
	inVal_hi = +1.;
	step     = step_Cos;   
      }
      
      if( !is_variable_used(nParam_j+p) ){
        minimizer->SetFixedVariable(nParam_j+p, name.c_str(), inVal );
	if(verbose>0)
	  printf("\tParam[%lu] = %s set to FIXED value %.0f\n", nParam_j+p, name.c_str(), inVal);
      }
      else{
	minimizer->SetLimitedVariable(nParam_j+p, name.c_str(),  inVal , step, inVal_lo , inVal_hi);   
	if(verbose>0)
	  printf("\tParam[%lu] = %s set to %.0f. Range: [%.2f,%.2f]\n", nParam_j+p,name.c_str(), inVal, inVal_lo, inVal_hi );	
      }
      
      if(event!=nullptr){
	event->treeStruct.obs_e   [ event->treeStruct.n_dim + count_param ] = p4_MET[0].p4.E(); 
	event->treeStruct.obs_pt  [ event->treeStruct.n_dim + count_param ] = p4_MET[0].p4.Pt();
	event->treeStruct.obs_eta [ event->treeStruct.n_dim + count_param ] = 0. ;
	event->treeStruct.obs_phi [ event->treeStruct.n_dim + count_param ] = 
	  (p%2==0 ? p4_MET[0].p4.Phi() : TMath::Cos(p4_MET[0].p4.Theta()) );
	event->treeStruct.obs_btag[ event->treeStruct.n_dim + count_param ] = 0.;
      }
      ++count_param;
    }


    for( size_t p = 0 ; p < nParam_m ; p++){
      
      string name;
      double inVal, inVal_lo, inVal_hi, step;
      if( p%2==0 ){
	name     = "phi"+to_string(nParam_j+nParam_n+p);
        inVal    = 0.;
        inVal_lo = -TMath::Pi() ;
        inVal_hi = +TMath::Pi() ;
        step     = step_Phi;
      }
      else{
	name     = "cos"+to_string(nParam_j+nParam_n+p/2);
        inVal    = 0.;
        inVal_lo = -1.;
        inVal_hi = +1.;
        step     = step_Cos;
      }

      if( !is_variable_used(nParam_j+nParam_n+p) ){
        minimizer->SetFixedVariable(nParam_j+nParam_n+p, name.c_str(), inVal );
        if(verbose>0)
          printf("\tParam[%lu] = %s set to FIXED value %.0f\n", nParam_j+nParam_n+p, name.c_str(), inVal);
      }
      else{
        minimizer->SetLimitedVariable(nParam_j+nParam_n+p, name.c_str(),  inVal , step, inVal_lo , inVal_hi);
        if(verbose>0)
          printf("\tParam[%lu] = %s set to %.0f. Range: [%.2f,%.2f]\n", nParam_j+nParam_n+p,name.c_str(), inVal, inVal_lo, inVal_hi );
      }

      if(event!=nullptr){
        event->treeStruct.obs_e   [ event->treeStruct.n_dim + count_param ] = 0.;
        event->treeStruct.obs_pt  [ event->treeStruct.n_dim + count_param ] = 0.;
        event->treeStruct.obs_eta [ event->treeStruct.n_dim + count_param ] = 0.;
        event->treeStruct.obs_phi [ event->treeStruct.n_dim + count_param ] = 0.;
        event->treeStruct.obs_btag[ event->treeStruct.n_dim + count_param ] = 0.;
      }
      ++count_param;   
    }
        
    // sanity check
    assert( minimizer->NDim()==count_param );
    break;
  }
  
}


void Algo::HypoTester::run(){

  if(verbose>0){ cout << "Algo::HypoTester::run()" << endl; }

  if(minimizer==nullptr){
    cout << "\tNullptr for minimizer: run() returns" << endl;
    return;
  }
  if(permutations.size()==0){
    cout << "\tNo permutations: run() returns" << endl; 
    return;
  }

  ROOT::Math::Functor f0( this , &Algo::HypoTester::eval, nParam_j + nParam_n + nParam_m); 
  minimizer->SetFunction(f0);

  Strategy strategy = Strategy::FirstTrial;

  // minimize nll
  setup_minimizer( strategy );
  auto t0 = high_resolution_clock::now();
  minimizer->Minimize(); 
  if( minimizer->Status()!=0 ){
    cout << "Minimizer failed with strategy " 
	 <<  static_cast<int>(Strategy::FirstTrial) 
	 << ": try with strategy " 
	 << static_cast<int>(Strategy::StartFromLastMinimum) << endl;
    strategy =  Strategy::StartFromLastMinimum;
    setup_minimizer( strategy );
    minimizer->Minimize();
    if( minimizer->Status()==0 )  cout << "...success" << endl;
    else cout << "...failure" << endl; 
  }
  auto t1 = high_resolution_clock::now();  

  int time_0 = static_cast<int>(duration_cast<milliseconds>(t1-t0).count());
  if(verbose>0){
    cout << "\tMinimization done in " << time_0 << " msec" << endl;
  }
  
  // get the result
  double nll        = minimizer->MinValue();
  int status        = minimizer->Status();
  const size_t ndim = minimizer->NDim();

  // get the results
  const double *xs = minimizer->X();

  if(verbose>0){
    cout << "\tStatus " << status << ", Minimum = " << nll  << endl;
    for( size_t var = 0 ; var < ndim ; ++var)
      cout << "\tVar[" << var << "] = " << xs[var] << endl;
  }
  
  if(event!=nullptr){

    // hypothesis counter
    ++(event->treeStruct.n_h);

    // add nll and number of dimension per each hypo
    if(count_hypo<HMAX){
      event->treeStruct.nll     [(size_t)count_hypo] = nll;
      event->treeStruct.status  [(size_t)count_hypo] = status;
      event->treeStruct.strategy[(size_t)count_hypo] = static_cast<int>(strategy);
      event->treeStruct.min_time[(size_t)count_hypo] = time_0;
      event->treeStruct.dim     [(size_t)count_hypo] = (int)ndim;    
      event->treeStruct.perm    [(size_t)count_hypo] = (int)count_perm;    
    }

    // add value of parameters per each hypo:
    // [ [p0_0 p0_1 ... p0_N ], ... , [ [pH_0 pH_1 ... pH_N] ]
    for( size_t p = 0 ; p < ndim && p < PMAX; ++p){
      event->treeStruct.param[ event->treeStruct.n_dim + p ] = xs[p];
    }

    // increment total dim counter
    (event->treeStruct.n_dim) += ndim;
  }

}


double Algo::HypoTester::eval(const double* xx){

  double val {0.};

  if(verbose>2){
    cout << "#iter: " ;
    for( size_t dim = 0 ; dim <  minimizer->NDim() ; ++dim)
      cout << "xx[" << dim << "] = " << xx[dim] << ", ";    
    cout << endl;
  }

  int count {0};
  for( auto perm : permutations){
    if(verbose>2) cout << "Algo::HypoTester::eval(): Eval perm " << count << endl;
    val += perm->eval( xx );
    ++count;
  }
  
  if( TMath::IsNaN(val) )
    return numeric_limits<double>::max();
  else if( val<=0 )
    return numeric_limits<double>::max();
  else{
    double ret_val = -TMath::Log(count_perm>0 ? val/count_perm : val);
    if(verbose>2) cout << "-Log(f) = " << ret_val << endl;
    return ret_val;
  }
  
  
  return val;

}

bool Algo::HypoTester::is_variable_used( const size_t pos ){

  for( auto perm : permutations ){
    for(size_t i = 0 ; i < perm->size() ; ++i){

      Algo::DecayBuilder* decay = perm->at(i);
      vector<size_t> vars ;

      switch( decay->get_decay() ){
      case Decay::Decay::TopHad:
	vars = (static_cast<TopHadBuilder*>(decay))->get_variables();
	if( pos==vars[0] ||  pos==vars[1] ||  pos==vars[2] ) return true;
	break;
      case Decay::Decay::TopHadLost:
	vars = (static_cast<TopHadBuilder*>(decay))->get_variables();
	if( pos==vars[0] ) return true;
	break;
      case Decay::Decay::WHad:
	vars = (static_cast<WHadBuilder*>(decay))->get_variables();
	if( pos==vars[0] ) return true; 
	break;
      case Decay::Decay::Higgs:
	vars = (static_cast<HiggsBuilder*>(decay))->get_variables();  
	if( pos==vars[0] ) return true; 
	break;
      case Decay::Decay::TopLep:
	vars = (static_cast<TopLepBuilder*>(decay))->get_variables(); 
	if( pos==vars[0] || pos==vars[1] ) return true;
	break;
      case Decay::Decay::Radiation_u:
      case Decay::Decay::Radiation_d:
      case Decay::Decay::Radiation_c:
      case Decay::Decay::Radiation_g:
      case Decay::Decay::Radiation_b:
	vars = (static_cast<RadiationBuilder*>(decay))->get_variables();
        if( pos==vars[0] ) return true;
        break;
      default:
	break;
      }
    }
  }

  return false;
}


void Algo::HypoTester::print(ostream& os){

  cout << "Algo::HypoTester::print()" << endl;

  os << " -Content:" << endl;
  os << "\t -jets:" << endl;
  int count_j = 0;
  for( auto jet : p4_Jet ){
    os << "\t\tjet[" << count_j << "]: (" << jet.p4.Pt() << "," << jet.p4.Eta() << "," << jet.p4.Phi() << "," << jet.p4.M() << ")" << endl;
    for( auto it = (jet.obs).begin() ; it !=  (jet.obs).end(); ++it)
      os << "\t\t" << it->first << " = " << it->second << endl;
    ++count_j;
  }

  os << "\t -leptons:" << endl;
  int count_l = 0;
  for( auto lepton : p4_Lepton ){
    os << "\t\tlepton[" << count_l << "]: (" << lepton.p4.Pt() << "," << lepton.p4.Eta() << "," << lepton.p4.Phi() << "," << lepton.p4.M() << ")" << endl;
    for( auto it = (lepton.obs).begin() ; it !=  (lepton.obs).end(); ++it)
      os << "\t\t" << it->first << " = " << it->second << endl;
    ++count_l;
  }

  os << "\t -MET:" << endl;
  int count_m = 0;
  for( auto MET : p4_MET ){
    os << "\t\tMET[" << count_l << "]: (" << MET.p4.Pt() << "," << MET.p4.Eta() << "," << MET.p4.Phi() << "," << MET.p4.M() << ")" << endl;
    os << "\t\t\t cos(t)=" << TMath::Cos( MET.p4.Theta() ) << endl;
    for( auto it = (MET.obs).begin() ; it !=  (MET.obs).end(); ++it)
      os << "\t\t" << it->first << " = " << it->second << endl;
    ++count_m;
  }

  os << " -Decay::Decays:" << endl;
  for( auto decay : decays )
    os << "\t" << int(decay) << " (=" << Algo::Decay::translateDecay(decay)  << ")" << endl;

  os << " -Partons:" << endl;
  int count_all = 0;  
  for( auto particle : particles ){
    os << "\t" << static_cast<int>(particle.first) << " (x" << particle.second << ")" << endl;
    ++count_all;
  }
  os << "\tTotal = " << count_all << " partons" << endl;

  os << "****************************************************" << endl;
}

void Algo::HypoTester::set_verbosity(const int& verb){
  verbose = verb;
}

int Algo::HypoTester::get_status(){
  return error_code;
}

void Algo::HypoTester::save_global_variables() {
  if(verbose>0){ cout << "Algo::HypoTester::save_global_variables()" << endl; }
  if(event==nullptr) return;

  for(auto jet : p4_Jet){
    double btag = (jet.obs).find("BTAG")!=(jet.obs).end() ? (jet.obs).find("BTAG")->second : 0.;
    if(btag>0.5) ++(event->treeStruct.n_btag);
    ++(event->treeStruct.n_jet);
  }
  for(auto lep : p4_Lepton){
    ++(event->treeStruct.n_lep);
  } 
}

void Algo::HypoTester::print_permutation(const vector<std::pair<Algo::FinalState::FinalState,size_t>>& perm) const {
  size_t count{0};
  cout << "\t[";
  for(auto p : perm){
    cout << "(" << p.second << ", " <<  static_cast<int>(p.first) << "), ";
    ++count;
  }
  cout << "]" << endl;
}

bool Algo::HypoTester::go_to_next(vector<vector<std::pair<Algo::FinalState::FinalState,size_t>>> & logbook){

  // filter out permutations with invalid quark <-> jet assignment                                                                                     
  // (only if BTAG_RND is filled and its value less than 0.5)                                                                                          
  if( !filter_by_btag( particles , p4_Jet ) ) return true;
  
  // skip unnecessary permutations: this can take some time                                                                                            
  int isnew {1};
  if(logbook.size()==0)
    logbook.push_back(particles);
  else{
    size_t count{0};
    for( auto log : logbook ){
      if( !(isnew = diff( log, particles , p4_Jet)) ) break;
      if(verbose>1){
	cout << "\t" << count << ": This permutation is new (code:" << isnew << ")" << endl;
	cout << "\t...Current: " << endl;
	print_permutation(particles);
	cout << "\t...Logbook: " << endl;
	print_permutation(logbook[count]);
      }
      ++count;
    }
    if(isnew){
      logbook.push_back(particles);
    }
    else{
      if(verbose>1){
	cout << "\tDo not consider this permutation, match between " << endl;
	cout << "\t...Current: " << endl;
	print_permutation(particles);
	cout << "\t...Logbook: " << endl;
	print_permutation(logbook[count]);
      }
      return true;
    }
  }
  return false;
}
