#include "interface/BTagRandomizer.h"


MEM::BTagRandomizer::BTagRandomizer(int debug, int seed, 
				    const std::map<DistributionType::DistributionType, TH3D>& pdf, 
				    int assignrnd,
				    int compress,
				    int nmax){
  btag_pdfs    = pdf;
  init         = 0;
  debug_code   = debug;
  n_tags_l     = 0;
  n_tags_h     = 0;
  n_jets       = 0;
  cut_val      = 0.;
  error_code   = 0;
  n_max_toys   = nmax;
  count_b      = 0;
  count_c      = 0;
  count_l      = 0;
  assign_rnd   = assignrnd;
  compress_csv = compress;
  tag_id       = 0;
  tag_name     = "";

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
  init       = 0;
  n_perm_max = 0;
  n_jets     = 0;
  count_b    = 0;
  count_c    = 0;
  count_l    = 0;  
  error_code = 0;
  tag_id     = 0;
  tag_name   = "";
  return;
}

void MEM::BTagRandomizer::next_category(){
  n_perm_max = 0;
  error_code = 0;
  tag_id     = 0;
  tag_name   = "";  
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
    if(compress_csv) discr = TMath::Max(discr,  0.);
    if(compress_csv) discr = TMath::Min(discr,  0.999999);

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

    // protection against underflow
    TH1D* h1 = h3.ProjectionZ(Form("%d_pz",int(j)), 
			      TMath::Max(binX,1),TMath::Max(binX,1),
			      TMath::Max(binY,1),TMath::Max(binY,1));

    if(!h1){
      cout << "BTagRandomizer::init_pdfs(): no histogram for jet[" << j << "] " << endl;   
      error_code = 1; 
      continue;
    }

    //deal with underflow
    if(compress_csv){
      if(h1->GetBinContent(0)>0. )
	h1->SetBinContent(1, h1->GetBinContent(0) + h1->GetBinContent(1));
      if(h1->GetBinContent(h1->GetNbinsX()+1)>0. )
	h1->SetBinContent(h1->GetNbinsX(), h1->GetBinContent(h1->GetNbinsX()) + h1->GetBinContent(h1->GetNbinsX()+1));
    }

    if( assign_rnd ){
      for(int k = 0 ; k < 2 ; ++k) discr = h1->GetRandom(); 
      if( debug_code&DebugVerbosity::init_more)
	cout << "BTagRandomizer::init_pdfs(): assigning random output " << discr << endl;
    }

    int bin     = h1->FindBin(cut_val);
    double pass = h1->Integral( bin , compress_csv ? h1->GetNbinsX() :  h1->GetNbinsX()+1 );
    pass -= (cut_val-h1->GetBinLowEdge( bin ))/h1->GetBinWidth( bin )*h1->GetBinContent(bin);
    pass /= h1->Integral(compress_csv ? 1 : 0,  compress_csv ? h1->GetNbinsX() : h1->GetNbinsX()+1 );

    pdfs.insert( make_pair(j, h1)    );
    effs.insert( make_pair(j, pass)  );
    vals.insert( make_pair(j, discr) );
    if( debug_code&DebugVerbosity::init_more)
      cout << "BTagRandomizer::init_pdfs(): adding <TH1D*, " << pass << "," << discr << "> for jet[" << j << "], bin(" << binX << "," << binY << ")" << endl;   
  }

  init = 1;
  return;
}

vector<MEM::BTagRandomizerOutput> MEM::BTagRandomizer::run_all(const vector<JetCategory>& cats){

  vector<BTagRandomizerOutput> out;

  for(size_t c = 0 ; c < cats.size() ; ++c){
    JetCategory cat = cats[c];
    tag_id   = cat.tag;
    tag_name = cat.name_tag;
    set_condition( cat.ntags_l, cat.ntags_h, cat.cut );
    BTagRandomizerOutput out_c = run();
    out.push_back(out_c);
    next_category();
  }

  return out;
}

MEM::BTagRandomizerOutput MEM::BTagRandomizer::run(){

  if( debug_code&DebugVerbosity::init_more )
    cout << "BTagRandomizer::run(): START" << endl;   

  BTagRandomizerOutput out = BTagRandomizerOutput();

  n_jets = int(jets.size());
  assert( n_jets >= n_tags_l );
  assert( n_tags_h<0 || n_tags_h>=n_tags_l );
  
  if(!init) init_pdfs();
  if(error_code>0 || !init) return out;

  perm_index.clear();
  for(size_t id = 0; id < size_t(n_jets) ; ++id) perm_index.push_back( id );
  
  sort( perm_index.begin(), perm_index.end(),  std::less<int>() );  
  vector<int> perm_index_copy = perm_index;

  n_perm_max = 0;
  std::map<int,size_t> n_perms_max;
  for(int n_tags = n_tags_l ; n_tags <= (n_tags_h>=0 ? TMath::Min(n_tags_h,n_jets) : n_jets) ; ++n_tags){
    n_perms_max[n_tags] = 0;
    do{ 
      if( debug_code&DebugVerbosity::init_more ) {
	cout << "\tperm. " << n_perm_max << ": [ ";
	for( auto ind : perm_index_copy ) cout << ind << " ";
	cout << "]" << endl;
      }
      ++n_perm_max; 
      ++(n_perms_max.at(n_tags));
    } 
    while( next_combination( perm_index_copy.begin(), perm_index_copy.begin()+n_tags, perm_index_copy.end(), std::less<int>() ) );
  }

  if( debug_code&DebugVerbosity::init_more ){
    cout << "\tMaximum of " << n_perm_max << " permutation(s) considered" << endl;
  }

  TH1D h_combinations("h_combinations","", n_perm_max+1, 0, n_perm_max+1);
  double pass{0.};
  size_t count_perm{0};
  for(int n_tags = n_tags_l ; n_tags <= (n_tags_h>=0 ? TMath::Min(n_tags_h,n_jets) : n_jets) ; ++n_tags){
    for( std::size_t n_perm = 0 ; n_perm < n_perms_max.at(n_tags) ; ++n_perm ){        

      double p{1.};
      if( debug_code&DebugVerbosity::event ) cout << "\tPermutation " << n_perm  << endl;
      auto perm = get_permutation(n_perm, n_tags);
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

      //p /= (TMath::Factorial( n_tags )*TMath::Factorial( n_jets-n_tags ));
      if( debug_code&DebugVerbosity::event ) cout << "\tFilling bin " << count_perm << endl;
      h_combinations.Fill( count_perm, p);
      pass += p;
      ++count_perm;
    }
  }

  size_t rnd_perm = h_combinations.GetBinLowEdge(h_combinations.FindBin(h_combinations.GetRandom()));
  int n_tags_rnd {0};
  int n_perm_rnd {0};
  count_perm = 0;
  for(int n_tags = n_tags_l ; n_tags <= (n_tags_h>=0 ? TMath::Min(n_tags_h,n_jets) : n_jets) ; ++n_tags){
    for( std::size_t n_perm = 0 ; n_perm < n_perms_max.at(n_tags) ; ++n_perm ){        
      if( count_perm==rnd_perm){
	n_tags_rnd = n_tags;
	n_perm_rnd = n_perm;
      }
      ++count_perm;
    }
  }


  if( debug_code&DebugVerbosity::init_more ) 
    cout << "Random permutation: " 
	 << rnd_perm << ", n tags: " << n_tags_rnd << ", n_perm: " << n_perm_rnd  << endl;

  if( debug_code&DebugVerbosity::init_more ) 
    h_combinations.Print("all");
  

  const size_t n = size_t(n_jets);
  MEM::Lock lk(n);
  auto perm = get_permutation(n_perm_rnd, n_tags_rnd);
  for(int t = 0 ; t < n_tags_rnd ; ++t){
    if( debug_code&DebugVerbosity::init_more ) cout << "\tLocking jet index [" << perm[t] << "]" << endl;
    lk.set_haslock( size_t(perm[t]), true );
  }

  vector<double> rnd_btag   = vector<double>(n);
  vector<double> input_btag = vector<double>(n);

  int ntoys = 0;
  while(ntoys<n_max_toys){

    ++ntoys;
    if( debug_code&DebugVerbosity::event ) 
      cout << "\tToy: " << ntoys << endl;

    int count_pass  {0};
    int count_pass_t{0};
    int count_pass_u{0};

    for(size_t j = 0; j < std::size_t(n_jets) ; ++j){

      if( ntoys>=2 && lk.get_haslock(j) && lk.get_lock(j) ){
	if( debug_code&DebugVerbosity::event ) 
	cout << "\tRandom value jet[" << j << "]: " << rnd_btag[j] << " (locked)"  << endl;
	++count_pass;
	++count_pass_t;	
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

      if(  rnd >= cut_val ){
	++count_pass;
	if( ntoys>=2 && lk.get_haslock(j) ){
	  lk.set_lock(j, true); 
	  ++count_pass_t;
	}
	if( ntoys>=2 && !lk.get_haslock(j) ) ++count_pass_u;
      }
    }

    if( debug_code&DebugVerbosity::event ) 
      cout << "\tPassing: " << count_pass << " (" << count_pass_t << "," << count_pass_u << ")" << endl;

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
      
    if( count_pass_t == n_tags_rnd && count_pass_u==0 ){
      out.pass_rnd = 1;
      break;
    }
    if(count_pass_t>n_tags_rnd)  lk.reset();
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
  out.n_jets     = n_jets;
  out.n_tags_l   = n_tags_l;
  out.n_tags_h   = n_tags_h;
  out.tag_id     = tag_id;
  out.tag_name   = tag_name;
  
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

std::vector<int> MEM::BTagRandomizer::get_permutation(const std::size_t& n, const int& ntag){
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
  } while( next_combination( perm_index_copy.begin(), perm_index_copy.begin()+ntag, perm_index_copy.end(), std::less<int>()) );

  return vector<int>{};
} 
