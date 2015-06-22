#include "interface/BTagRandomizer.h"


MEM::BTagRandomizer::BTagRandomizer(int debug, int seed, const std::map<DistributionType::DistributionType, TH3D>& pdf, int assignrnd, int nmax){
  btag_pdfs  = pdf;
  debug_code = debug;
  n_tags_l   = 0;
  n_tags_h   = 0;
  n_jets     = 0;
  cut_val    = 0.;
  error_code = 0;
  n_max_toys = nmax;
  count_b    = 0;
  count_c    = 0;
  count_l    = 0;
  assign_rnd = assignrnd;
  ran = new TRandom3();
  if(seed<0) gRandom->SetSeed(0);
  if(seed>0) gRandom->SetSeed(seed);
}

MEM::BTagRandomizer::~BTagRandomizer(){
  jets.clear();
  pdfs.clear();
  effs.clear();
  vals.clear();
  perm_index.clear();
  delete ran;
}

void MEM::BTagRandomizer::next_event(){
  jets.clear();
  pdfs.clear();
  effs.clear();
  vals.clear();
  perm_index.clear();
  n_perm_max = 0;
  n_jets     = 0;
  count_b    = 0;
  count_c    = 0;
  count_l    = 0;  
  error_code = 0;
  return;
}

void MEM::BTagRandomizer::init_pdfs(){

  if( btag_pdfs.size()<1 ){
    cout << "BTagRandomizer::init_pdfs(): no btag pdf's provided" << endl;   
    error_code = 1;
    return;
  }

  for(size_t j = 0; j < size_t(n_jets) ; ++j){

    if( !jets[j]->isSet(Observable::PDGID) ) {
      cout << "BTagRandomizer::init_pdfs(): jet[" << j << "] has no no PDG ID info" << endl;   
      error_code = 1; 
      continue;
    }

    if( !jets[j]->isSet(Observable::CSV) ) {
      if(!assign_rnd){
	cout << "BTagRandomizer::init_pdfs(): jet[" << j << "] has no no CSV discriminant info" << endl;   
	error_code = 1; 
	continue;
      }
      else{
	if( debug_code&DebugVerbosity::init_more)
	  cout << "BTagRandomizer::init_pdfs(): get random CSV" << endl;
      }
    }

    double discr = jets[j]->isSet(Observable::CSV) ? jets[j]->getObs(Observable::CSV) : -99.; 
    int pdg      = int(jets[j]->getObs(Observable::PDGID));

    if( debug_code&DebugVerbosity::init_more)
      cout << "BTagRandomizer::init_pdfs(): init jet[" << j << "] with (Pt,Eta,PDGID,csv) = (" 
	   << jets[j]->p4().Pt() << ", " << jets[j]->p4().Eta() << ", " << pdg
	   << ", " << discr << ")" << endl;   
    
    DistributionType::DistributionType type = DistributionType::DistributionType::csv_l;

    if( std::abs(pdg)<4 || std::abs(pdg)==21 ){
      type = DistributionType::DistributionType::csv_l;
      ++count_l;
    }
    else if( std::abs(pdg)==4 ){
      type = DistributionType::DistributionType::csv_c;
      ++count_c;
    }
    else if( std::abs(pdg)==5 ){
      type = DistributionType::DistributionType::csv_b;
      ++count_b;
    }
    else{
      cout << "BTagRandomizer::init_pdfs(): no pdf available for jet[" << j 
	   << "] (PDGID=" << pdg << ")" << endl;   
      error_code = 1; 
      continue;      
    }
  
    TH3D  h3 = btag_pdfs.at(type);
    h3.SetDefaultSumw2();
    int binX = h3.GetXaxis()->FindBin(jets[j]->p4().Pt());
    int binY = h3.GetYaxis()->FindBin(std::abs(jets[j]->p4().Eta()));
    TH1D* h1 = h3.ProjectionZ(Form("%d_pz",int(j)), binX,binX,binY,binY);
    if( assign_rnd ){
      for(int k = 0 ; k < 2 ; ++k) discr = h1->GetRandom(); 
    }

    if(!h1){
      cout << "BTagRandomizer::init_pdfs(): no histogram for jet[" << j << "] " << endl;   
      error_code = 1; 
      continue;
    }

    int bin = h1->FindBin(cut_val);
    double pass = h1->Integral( bin , h1->GetNbinsX() );
    pass -= (cut_val-h1->GetBinLowEdge( bin ))/h1->GetBinWidth( bin )*h1->GetBinContent(bin);
    pass /= h1->Integral(0,  h1->GetNbinsX() );

    pdfs.insert( make_pair(j, h1)    );
    effs.insert( make_pair(j, pass)  );
    vals.insert( make_pair(j, discr) );
    if( debug_code&DebugVerbosity::init_more)
      cout << "BTagRandomizer::init_pdfs(): adding <TH1D*, " << pass << "," << discr << "> for jet[" << j << "], bin(" << binX << "," << binY << ")" << endl;   
  }

  
}


MEM::BTagRandomizerOutput MEM::BTagRandomizer::run(){

  if( debug_code&DebugVerbosity::init_more )
    cout << "BTagRandomizer::run(): START" << endl;   

  BTagRandomizerOutput out = BTagRandomizerOutput();

  n_jets = int(jets.size());
  assert( n_jets >= n_tags_l );
  
  init_pdfs();
  if(error_code>0) return out;

  perm_index.clear();
  for(size_t id = 0; id < size_t(n_jets) ; ++id) perm_index.push_back( id );
  
  comparator = CompPerm(1);
  sort( perm_index.begin(), perm_index.end(), comparator );  
  vector<int> perm_index_copy = perm_index;
  n_perm_max = 0;
  do{ ++n_perm_max; } 
  while( next_permutation( perm_index_copy.begin(), perm_index_copy.end(), comparator) );

  if( debug_code&DebugVerbosity::init_more ){
    cout << "\tMaximum of " << n_perm_max << " permutation(s) considered" << endl;
  }

  double pass{0.};
  for(int n_tags = n_tags_l ; n_tags <= (n_tags_h>0 ? TMath::Min(n_tags_h,n_jets) : n_jets) ; ++n_tags){
    for( std::size_t n_perm = 0 ; n_perm < n_perm_max ; ++n_perm ){        

      double p{1.};
      if( debug_code&DebugVerbosity::event ) cout << "\tPermutation " << n_perm  << endl;
      auto perm = get_permutation(n_perm);
      if(perm.size()==0) continue;
      
      for( int j = 0 ; j < n_jets ; ++j){

	double p_t = 1.;
	if( jets.at(perm[j])->isSet( Observable::Observable::BTAGPROB ) )
	  p_t = jets.at(perm[j])->getObs( Observable::Observable::BTAGPROB );
	else{
	  p_t = effs.at(perm[j]);
	}

	if(j<n_tags){
	  if( debug_code&DebugVerbosity::event ) cout << "\tPass: " << p_t  << endl;
	  p *= p_t;
	}
	else{
	  if( debug_code&DebugVerbosity::event ) cout << "\tFail: " << (1-p_t)  << endl;
	  p *= (1.-p_t);
	}
      }        

      if( debug_code&DebugVerbosity::event ) cout << "\tp = : " << p  << endl;

      p /= (TMath::Factorial( n_tags )*TMath::Factorial( n_jets-n_tags ));
      pass += p;
    }
  }

  
  const size_t n = size_t(n_jets);
  MEM::Lock lk(n);
  vector<double> rnd_btag   = vector<double>(n);
  vector<double> input_btag = vector<double>(n);

  int ntoys = 0;
  while(ntoys<n_max_toys){

    ++ntoys;
    if( debug_code&DebugVerbosity::event ) 
      cout << "\tToy: " << ntoys << endl;

    int count_pass{0};
    for(size_t j = 0; j < std::size_t(n_jets) ; ++j){

      if( lk.get(j) ){
	if( debug_code&DebugVerbosity::event ) 
	cout << "\tRandom value jet[" << j << "]: " << rnd_btag[j] << " (locked)"  << endl;
	continue;
      }

      double rnd{0.};
      if(ntoys<2){ 
	rnd           = vals.at(j);
	input_btag[j] = rnd;
      }
      else{
	rnd         = pdfs.at(j)->GetRandom();
	rnd_btag[j] = rnd;
      }

      if( debug_code&DebugVerbosity::event ){
	if(ntoys<2) cout << "\tInput value  jet[" << j << "]: " << input_btag[j] << endl;
	else        cout << "\tRandom value jet[" << j << "]: " << rnd_btag[j] << endl;
      }

      if( rnd >= cut_val ){
	lk.set(j,true);
	++count_pass;
      }
    }
    if( debug_code&DebugVerbosity::event ) 
      cout << "\tPassing: " << count_pass << endl;

    // burn first toy
    if( ntoys < 2){
      if( (n_tags_h<0         && count_pass>=n_tags_l) ||
	  (n_tags_h==n_tags_l && count_pass==n_tags_l) ||
	  (n_tags_h>n_tags_l  && count_pass>=n_tags_l && count_pass<=n_tags_h) )
	out.pass = 1;
      if( debug_code&DebugVerbosity::event ) 
	cout << "\tBurn this toy....continue" << endl; 
      lk.reset(); 
      continue;
    }
      
    if( (n_tags_h<0         && count_pass>=n_tags_l) ||
	(n_tags_h==n_tags_l && count_pass==n_tags_l) ||
	(n_tags_h>n_tags_l  && count_pass>=n_tags_l && count_pass<=n_tags_h) ||
	lk.allset() ){
      out.pass_rnd = 1;
      break;
    }
    if(n_tags_h==n_tags_l && count_pass>n_tags_l)  lk.reset();
    if(n_tags_h>n_tags_l  && count_pass>n_tags_h)  lk.reset();
  }

  out.ntoys      = ntoys;  
  if( out.pass_rnd ){
    out.p          = pass;  
    out.rnd_btag   = rnd_btag;
  }
  else{
    out.p          = 1.;  
    out.rnd_btag   = input_btag;
  }
  out.input_btag = input_btag;
  out.err        = error_code;  
  out.n_b        = count_b;  
  out.n_c        = count_c;  
  out.n_l        = count_l;  
  out.n_tags_l   = n_tags_l;
  out.n_tags_h   = n_tags_h;

  return out;
}

void MEM::BTagRandomizer::set_condition(const int& ntags_l, const int& ntags_h, const double& v){
  n_tags_l = ntags_l;
  n_tags_h = ntags_h;
  cut_val  = v;
}

void MEM::BTagRandomizer::push_back_object( Object* obj){

  jets.push_back( obj );
  
  if( debug_code&DebugVerbosity::init_more ){
    cout << "BTagRandomizer::push_back_object()" << endl;
    obj->print(cout);
  }

  return;
}

std::vector<int> MEM::BTagRandomizer::get_permutation(const std::size_t& n){
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
