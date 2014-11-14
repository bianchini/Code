#include "interface/HypoTester.h"


Algo::HypoTester::HypoTester(){

  // reset variables
  nParam_j          = 0;
  nParam_n          = 0;
  count_hypo        = 0;
  count_perm        = 0;
  count_TopHad      = 0;
  count_TopLep      = 0;
  count_Higgs       = 0;
  count_Radiation   = 0;
  invisible         = 0;
  verbose           = 0;
  event             = nullptr;
}

Algo::HypoTester::HypoTester(TTree* t){

  // reset variables                                                                                                                                
  nParam_j          = 0;
  nParam_n          = 0;
  count_hypo        = 0;
  count_perm        = 0;
  count_TopHad      = 0;
  count_TopLep      = 0;
  count_Higgs       = 0;
  count_Radiation   = 0;
  invisible         = 0;
  verbose           = 0;

  cout << "Creating branches to output file" << endl;
  event  = new Event(t);
  event->createBranches();
  event->reset();
}

Algo::HypoTester::~HypoTester(){
  cout << "Removing HypoTester" << endl;
  for( auto perm : permutations ) delete perm;
  //delete minimizer;
  delete event;
}


void Algo::HypoTester::push_back_object( const LV& p4, char type ){

  Object obj;
  obj.init( p4 );

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
    cout << "Algo::HypoTester::push_back_object: Unknown type of object added" << endl;
    break;
  }

}

void Algo::HypoTester::add_object_observables( const string& name, const double val, char type){

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
    cout << "Algo::HypoTester::add_object_observables: Unknown type of object added" << endl;
    break;
  }

}


void Algo::HypoTester::test( const map<string, vector<Decay>>& all ){

  if(event!=nullptr) event->reset();

  // make sure we start from zero
  count_hypo = 0;
  for( auto hypo : all ){
    cout << "Testing hypothesis " << count_hypo 
	 << " with name \"" << hypo.first << "\"" << endl;
    reset();
    for( auto decay : hypo.second ){
      assume( decay );
    }    
    init();
    run();
    ++count_hypo;
  }

  if(event!=nullptr) event->fillTree();
  next_event();
}

void Algo::HypoTester::reset(){

  // remove old permutations
  // perm calls destructor of block
  for( auto perm : permutations ) delete perm;
  
  decays.clear();
  particles.clear();
  permutations.clear();
  nParam_j        = 0;
  nParam_n        = 0;
  count_perm      = 0;
  count_TopHad    = 0;
  count_WHad      = 0;
  count_TopLep    = 0;
  count_Higgs     = 0;
  count_Radiation = 0;
  invisible       = 0;
}

void Algo::HypoTester::next_event(){
  reset();
  p4_Jet.clear();
  p4_Lepton.clear();
  p4_MET.clear();
}

void Algo::HypoTester::assume( Decay decay ){
  decays.push_back( decay );
}


void Algo::HypoTester::unpack_assumptions(){

  if(verbose>0){ cout << "HypoTester::unpack_assumptions()" << endl; }

  // reset
  nParam_j          = 0;
  nParam_n          = 0;
  count_TopHad      = 0;
  count_WHad        = 0;
  count_TopLep      = 0;
  count_Higgs       = 0;
  count_Radiation   = 0;

  for( auto decay : decays ){

   switch( decay ){

   case Algo::Decay::TopHad:
     particles.push_back( make_pair( FinalState::TopHad_q,    count_TopHad) );
     particles.push_back( make_pair( FinalState::TopHad_qbar, count_TopHad) );
     particles.push_back( make_pair( FinalState::TopHad_b,    count_TopHad) );
     ++count_TopHad;
     nParam_j += 3;
     if(verbose>0){ cout << "\tAdded TopHad" << endl; }
     break;

   case Algo::Decay::WHad:
     particles.push_back( make_pair( FinalState::WHad_q,    count_WHad) );
     particles.push_back( make_pair( FinalState::WHad_qbar, count_WHad) );
     ++count_WHad;
     nParam_j += 2;
     if(verbose>0){ cout << "\tAdded WHad" << endl; }
     break;

   case Algo::Decay::TopLep:
     particles.push_back( make_pair( FinalState::TopLep_b,    count_TopLep) );
     ++count_TopLep;
     ++invisible;
     nParam_j += 1;
     nParam_n += 2;
     if(verbose>0){ cout << "\tAdded TopLep" << endl; }
     break;

   case Algo::Decay::Higgs:
     particles.push_back( make_pair( FinalState::Higgs_b,    count_Higgs) );
     particles.push_back( make_pair( FinalState::Higgs_bbar, count_Higgs) );
     ++count_Higgs;
     nParam_j += 2;
     if(verbose>0){ cout << "\tAdded Higgs" << endl; }
     break;

   case Algo::Decay::Radiation:
     particles.push_back( make_pair( FinalState::Radiation_q,   count_Radiation) );
     if(verbose>0){ cout << "\tAdded Radiation" << endl; }
     ++count_Radiation;
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
      particles.push_back( make_pair( FinalState::Radiation_q,   count_Radiation) );
      if(verbose>0){ cout << "\tAdded Radiation" << endl; }
      ++count_Radiation;
      nParam_j += 1;
      --addradiation;
    }

  }
    
  if(verbose>0){ cout << "\tTotal number of parameters: nParam_j=" << nParam_j << ", nParam_n=" << nParam_n << endl; }

}


void Algo::HypoTester::init(){

  if(verbose>2){ cout << "HypoTester::init()" << endl; }

  // first, unpack assumptions
  this->unpack_assumptions();

  // number of jets at least as large as num of quarks
  assert( particles.size()<=p4_Jet.size() );

  // if >0 topLep, need leptons
  assert( count_TopLep==p4_Lepton.size()  );

  // if >0 topLep, need MET
  assert( count_TopLep==0 || (count_TopLep>0 && p4_MET.size()>0)  );
  
  // order particles
  sort( particles.begin(), particles.end(), MyComp );

  // permutations visited
  count_perm = 0;

  // bookkeep to remove unnecessary permutations
  vector<vector<std::pair<FinalState,size_t>>> logbook;

  // permutatios using std library
  do {
    
    // skip unnecessary permutations
    bool skip = false;
    if(logbook.size()==0) 
      logbook.push_back(particles);
    else{
      for( auto log : logbook ){
	if( isSame( log, particles ) ){
	  skip = true;
	}
      }
      if(!skip) 
	logbook.push_back(particles);             
      else{
	if(verbose>2) cout << "\tDo not consider this permutation  " << endl;
	continue;
      }
    }
    

    // this vector will contain all blocks 
    // E.g. block = TopHad * TopLep * Rad * Rad * ... * MET
    vector<Algo::DecayBuilder*> decays; 

    // pass by reference: here the group is formed
    group_particles( decays );
    
    // add block list of permutations
    // a block is assembled by CombBuilder
    permutations.push_back( new Algo::CombBuilder(decays, verbose) );
    
    // increment counter
    ++count_perm;
    
  } while ( next_permutation(particles.begin(), particles.end(), MyComp  ) );


   cout << "\tTotal of  " << count_perm << " permutations created:" << endl;
   if(verbose>0) {
     size_t count = 0;
     for( auto perm : permutations){
       cout << "(" << count << ")";
       perm->print(cout);
       ++count;
     }
   }

}



void Algo::HypoTester::group_particles(vector<DecayBuilder*>& decayed){

  if(verbose>2){ cout << "HypoTester::group_particles()" << endl; }

  //  get all hadronically decaying tops
  for( size_t t_had = 0; t_had < count_TopHad; ++t_had ){

    if(verbose>1) cout << "\tProcessing " << t_had << "th TopHad" << endl;

    Algo::TopHadBuilder* topHad = new Algo::TopHadBuilder(verbose);
    
    size_t pos = 0;
    for( auto part : particles ){

      if( part.second == t_had )  
	topHad->init( part.first, p4_Jet[pos].p4 , pos );
      ++pos;
    }

    decayed.push_back( topHad );

  }

  // get all hadronically decaying tops
  for( size_t w_had = 0; w_had < count_WHad; ++w_had ){

    if(verbose>1) cout << "\tProcessing " << w_had << "th WHad" << endl;

    Algo::WHadBuilder* wHad = new Algo::WHadBuilder(verbose);
    
    size_t pos = 0;
    for( auto part : particles ){
      
      if( part.second == w_had ) 
	wHad->init( part.first, p4_Jet[pos].p4 , pos );
      ++pos;
    }

    decayed.push_back( wHad );

  }

  //  get all leptonically decaying tops
  for( size_t t_lep = 0; t_lep < count_TopLep; ++t_lep ){

    if(verbose>1) cout << "\tProcessing " << t_lep << "th TopLep" << endl;

    Algo::TopLepBuilder* topLep = new Algo::TopLepBuilder(verbose);
    
    topLep->init( FinalState::TopLep_l , p4_Lepton[t_lep].p4 , t_lep + nParam_j );

    size_t pos = 0;
    for( auto part : particles ){

      if( part.second == t_lep ) 
	topLep->init( part.first, p4_Jet[pos].p4 , pos );
      ++pos;
    }


    decayed.push_back( topLep );

  }


  //  get all radiation
  for( size_t r_had = 0; r_had < count_Radiation; ++r_had ){
    
    if(verbose>1) cout << "\tProcessing " << r_had << "th Radiaton" << endl;

    Algo::RadiationBuilder* rad = new Algo::RadiationBuilder(verbose);
    
    size_t pos = 0;
    for( auto part : particles ){

      if( part.second == r_had )  
	rad->init( part.first, p4_Jet[pos].p4 , pos );
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
	hig->init( part.first, p4_Jet[pos].p4 , pos );
      ++pos;
    }

    decayed.push_back( hig );

  }


  /* Add other blocks here */


  // if there are invisible particles, add MET as last element                                                                            
  if(invisible>0){

    Algo::METBuilder* met = new Algo::METBuilder(verbose);
    assert( p4_MET.size()>0 );
    met->init( p4_MET[0].p4 );
    decayed.push_back( met );
  }
  // if there are not insible particles, but MET is an input, eval MET tf at (0.,0.)                                                      
  else if( p4_MET.size()>0 ){

    Algo::METBuilder* met = new Algo::METBuilder(verbose);
    met->init( p4_MET[0].p4 );
    met->fix_vars();
    decayed.push_back( met );
  }
  else{ /* do nothing */ }


}

void Algo::HypoTester::run(){

  if(verbose>2){ cout << "HypoTester::run()" << endl; }

  minimizer =  ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetMaxIterations(10000);
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(0);

  ROOT::Math::Functor f0( this , &Algo::HypoTester::eval, nParam_j + nParam_n); 
  minimizer->SetFunction(f0);

  for( size_t p = 0 ; p < nParam_j ; p++){
    char name[6];
    sprintf(name, "j%lu", p);

    double inVal    {p4_Jet[p].p4.E()};
    double inVal_lo {0.};
    double inVal_hi {inVal*3.0};
    double step     {0.5};

    minimizer->SetLimitedVariable(p, name, inVal  , step, inVal_lo , inVal_hi);

    if(verbose>2)
      printf("\tParam[%lu] = %s set to %.0f. Range: [%.0f,%.0f]\n", p,name, inVal, inVal_lo, inVal_hi );
    
  }

  for( size_t p = 0 ; p < nParam_n ; p++){

    char name[6];
    double inVal, inVal_lo, inVal_hi, step;

    if( p%2==0 ){
      sprintf(name, "phi%lu", p);
      inVal    = p4_MET[0].p4.Phi();
      inVal_lo = -TMath::Pi() ;
      inVal_hi = +TMath::Pi() ;
      step     = 0.1;     
    }
    else{
      sprintf(name, "cos%lu", p/2);
      inVal    = 0.;
      inVal_lo = -1. ;
      inVal_hi = +1. ;
      step     = 0.1;   
    }

    minimizer->SetLimitedVariable(nParam_j+p, name,  inVal , step, inVal_lo , inVal_hi);
   
    if(verbose>2)
      printf("\tParam[%lu] = %s set to %.0f. Range: [%.2f,%.2f]\n", nParam_j+p,name, inVal, inVal_lo, inVal_hi );
    
  }
  
  // minimie nll
  minimizer->Minimize(); 

  // get the result
  double nll  = minimizer->MinValue();
  const size_t ndim = minimizer->NDim();
 
  const double *xs = minimizer->X();

  if(verbose>-1){
    cout << "\tMinimum = " << nll  << endl;
    for( size_t var = 0 ; var < ndim ; ++var)
      cout << "\tVar[" << var << "] = " << xs[var] << endl;
  }
  
  if(event!=nullptr){
    // hypothesis counter
    ++(event->treeStruct.n_h);

    // add nll and number of dimension per each hypo
    if(count_hypo<HMAX){
      event->treeStruct.nll[(size_t)count_hypo] = nll;
      event->treeStruct.dim[(size_t)count_hypo] = ndim;    
    }

    // add value of parameters per each hypo:
    // [ [p0_0 p0_1 ... p0_N ], ... , [ [pH_0 pH_1 ... pH_N] ]
    for( size_t p = 0 ; p < ndim && p < PMAX; ++p){
      event->treeStruct.param[ event->treeStruct.n_dim + p ] = xs[p];
    }

    // increment total dim counter
    (event->treeStruct.n_dim) += ndim;
  }

  delete minimizer;

}


double Algo::HypoTester::eval(const double* xx){

  double val {0.};

  int count {0};
  for( auto perm : permutations){
    if(verbose>2) cout << "Eval perm " << count << endl;
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


void Algo::HypoTester::print(ostream& os){

  cout << "HypoTester::print()" << endl;

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
    for( auto it = (MET.obs).begin() ; it !=  (MET.obs).end(); ++it)
      os << "\t\t" << it->first << " = " << it->second << endl;
    ++count_m;
  }

  os << " -Decays:" << endl;
  for( auto decay : decays )
    os << "\t" << int(decay) << " (=" << Algo::translateDecay(decay)  << ")" << endl;

  os << " -Partons:" << endl;
  int count_all = 0;  
  for( auto particle : particles ){
    os << "\t" << particle.first << " (x" << particle.second << ")" << endl;
    ++count_all;
  }
  os << "\tTotal = " << count_all << " partons" << endl;

  os << "****************************************************" << endl;
}

