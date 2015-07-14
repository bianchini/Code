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

  if( debug_code&DebugVerbosity::init){
    cout << "BTagRandomizer::BTagRandomizer(): " << endl;
    cout << "\tRandom seed:        " << gRandom->GetSeed() << endl;
    cout << "\tPdfs initilialized: " << (btag_pdfs.size()>0) << endl;
    cout << "\tMax number of toys: " << n_max_toys << endl;
    if(assign_rnd)   cout << "\tWARNING: Sampling csv output from pdfs" << endl;
    if(compress_csv) cout << "\tWARNING: Compressing csv output in [0,1]" << endl;
  }

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
  for(auto pdf : pdfs)
    delete pdf.second;
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

    int pdg      = std::abs(int(jets[j]->getObs(Observable::PDGID)));
    int mcmatch  = jets[j]->isSet(Observable::MCMATCH) ? std::abs(int(jets[j]->getObs(Observable::MCMATCH))) : -1;

    if( debug_code&DebugVerbosity::init_more)
      cout << "BTagRandomizer::init_pdfs(): init jet[" << j << "] with (Pt,Eta,PDGID,csv) = (" 
	   << jets[j]->p4().Pt() << ", " << jets[j]->p4().Eta() << ", " << pdg
	   << ", " << discr << ")" << endl;   
    
    DistributionType::DistributionType type = DistributionType::DistributionType::csv_l;

    if( pdg<4 || pdg==21 ){
      type = DistributionType::DistributionType::csv_l;
      ++count_l;
      if     ( pdg==3  && btag_pdfs.find(DistributionType::DistributionType::csv_s) != btag_pdfs.end() )
	type = DistributionType::DistributionType::csv_s;
      else if( pdg<2   && btag_pdfs.find(DistributionType::DistributionType::csv_u) != btag_pdfs.end() )
	type = DistributionType::DistributionType::csv_u;
      else if( pdg==21 && btag_pdfs.find(DistributionType::DistributionType::csv_g) != btag_pdfs.end() )
	type = DistributionType::DistributionType::csv_g;
      else{ /*...*/ }
    }
    else if( pdg==4 ){
      type = DistributionType::DistributionType::csv_c;
      ++count_c;
      if     ( mcmatch==0 && btag_pdfs.find(DistributionType::DistributionType::csv_c_g) != btag_pdfs.end() )
	type = DistributionType::DistributionType::csv_c_g; 
      else if( mcmatch>0  && btag_pdfs.find(DistributionType::DistributionType::csv_c_t) != btag_pdfs.end() )
	type = DistributionType::DistributionType::csv_c_t;
      else{ /*...*/ }  
    }
    else if( pdg==5 ){
      type = DistributionType::DistributionType::csv_b;
      ++count_b;
      if     ( mcmatch==0 && btag_pdfs.find(DistributionType::DistributionType::csv_b_g) != btag_pdfs.end() )
	type = DistributionType::DistributionType::csv_b_g; 
      else if( mcmatch>0  && btag_pdfs.find(DistributionType::DistributionType::csv_b_t) != btag_pdfs.end() )
	type = DistributionType::DistributionType::csv_b_t; 
      else{ /*...*/ }  
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
      for(int k = 0 ; k < int(ran->Uniform(1,10)) ; ++k) discr = h1->GetRandom(); 
      if( debug_code&DebugVerbosity::init_more)
	cout << "BTagRandomizer::init_pdfs(): assigning random output " << discr << endl;
    }

    int bin     = h1->FindBin(cut_val);
    double pass = h1->Integral( bin , compress_csv ? h1->GetNbinsX() :  h1->GetNbinsX()+1 );
    pass -= (cut_val-h1->GetBinLowEdge( bin ))/h1->GetBinWidth( bin )*h1->GetBinContent(bin);
    double norm = h1->Integral(compress_csv ? 1 : 0,  compress_csv ? h1->GetNbinsX() : h1->GetNbinsX()+1 );
    pass /= (norm>0. ? norm : 1.0);

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

bool MEM::BTagRandomizer::check_ignore_jets(){

  int count_n{0};
  int count_t{0};
  
  for(size_t j = 0; j < size_t(n_jets) ; ++j){  
    if( jets.at(j)->isSet(Observable::IGNORE_FOR_RND) && jets.at(j)->getObs(Observable::IGNORE_FOR_RND)>0){
      if( vals.at(j) >=cut_val)
	++count_t;	
      else
	++count_n;   
    }
  }

  int residual = n_jets-(count_n+count_t);
  if( (count_t+residual)<n_tags_l){
    if( debug_code&DebugVerbosity::init_more )
      cout << "BTagRandomizer::check_ignore_jets(): only " << count_t << " tagged + " << residual << " residual jets, while a minimum of " 
	   << n_tags_l << " required" << endl;
    return false;
  }
  if( n_tags_h>=0 && count_t>n_tags_h){
    if( debug_code&DebugVerbosity::init_more )
      cout << "BTagRandomizer::check_ignore_jets(): already " << count_t << " tagged jets, while a maximum of " << n_tags_h << " required" << endl;
    return false;
  }

  // that many tagged jets are needed now...
  if( count_t <= n_tags_l){
    n_tags_l -= count_t;  
    n_tags_h -= (n_tags_h>=0 ? count_t : 0 ); 
  }
  else if( count_t > n_tags_l ){
    n_tags_l  = 0;
    n_tags_h -= (n_tags_h>=0 ? count_t : 0 ); 
  }

  return true;
}
 

bool MEM::BTagRandomizer::check_restore_jets(int& nr){

  int count_t{0};
  for(size_t j = 0; j < size_t(n_jets) ; ++j){  
    if( jets.at(j)->isSet(Observable::IGNORE_FOR_RND) && jets.at(j)->getObs(Observable::IGNORE_FOR_RND)>0 &&
	vals.at(j)>=cut_val)
      ++count_t;	
  }

  n_tags_l += count_t;
  nr       += count_t;
  n_tags_h += (n_tags_h>=0 ? count_t : 0 );

  return true;
}


void MEM::BTagRandomizer::fill_perm(){
  perm_index.clear();
  for(size_t id = 0; id < size_t(n_jets) ; ++id){
    if( jets.at(id)->isSet(Observable::IGNORE_FOR_RND) && jets.at(id)->getObs(Observable::IGNORE_FOR_RND)>0 ){
      perm_index.push_back( -1 );
    }
    else
      perm_index.push_back( id );
  }
  return;
}


MEM::BTagRandomizerOutput MEM::BTagRandomizer::run(){

  if( debug_code&DebugVerbosity::init_more )
    cout << "BTagRandomizer::run(): START" << endl;   

  const std::size_t n = jets.size();
  n_jets = int(n);
  assert( n_jets >= n_tags_l );
  assert( n_tags_h<0 || n_tags_h>=n_tags_l );

  BTagRandomizerOutput out = BTagRandomizerOutput();
  if(!init) init_pdfs();
  if(error_code>0 || !init){
    out.err = error_code;
    return out;
  }
  
  out.n_b        = count_b;  
  out.n_c        = count_c;  
  out.n_l        = count_l;  
  out.n_jets     = n_jets;
  out.n_tags_l   = n_tags_l;
  out.n_tags_h   = n_tags_h;
  out.tag_id     = tag_id;
  out.tag_name   = tag_name;

  vector<double> input_btag = vector<double>(n);
  vector<double> rnd_btag   = vector<double>(n);

  int count{0};
  for(size_t j = 0; j < n ; ++j){
    double val = vals.at(j); 
    if( val >= cut_val ) ++count;
    input_btag[j] = val;
    rnd_btag[j]   = val;
    if( debug_code&DebugVerbosity::event )
      printf("\tInput value  jet[%d]: %.3f %s\n", 
	     int(j), input_btag[j], 
	     ((jets.at(j)->isSet(Observable::IGNORE_FOR_RND) &&	jets.at(j)->getObs(Observable::IGNORE_FOR_RND)>0 ) ? " ( excluded )" : ""));
  }
  
  out.input_btag = input_btag;
  out.rnd_btag   = rnd_btag;
  out.pass       = int( count>=n_tags_l && (n_tags_h<0 || count<=n_tags_h) );
  

  // check if there are jets that I need to ignore  
  // It can be that given the ignored jets, the event cannot fall into the category
  if( !check_ignore_jets() ){
    if( debug_code&DebugVerbosity::init_more )
      cout << "This event should not be processed further: it " 
	   << "cannot pass the category cuts given the input. Return p=0" << endl;
    return out;
  }


  // fill the permutation word
  fill_perm();

  bool usable{false};
  for( auto p : perm_index ) usable |= (p>=0);
  if(!usable){
    out.pass_rnd = out.pass;
    out.p        = 1.0;
    return out;
  }

  sort( perm_index.begin(), perm_index.end(),  std::less<int>() );  
  vector<int> perm_index_copy = perm_index;

  n_perm_max = 0;
  std::map<int,size_t> n_perms_max;
  for(int n_tags = n_tags_l ; n_tags <= (n_tags_h>=0 ? TMath::Min(n_tags_h,n_jets) : n_jets) ; ++n_tags){
    n_perms_max[n_tags] = 0;
    if( debug_code&DebugVerbosity::init_more ) 
      cout << "Tags: " << n_tags << endl;
    do{ 
      if( consider(perm_index_copy, n_tags) ){
	if( debug_code&DebugVerbosity::init_more ) {
	  int count{0};
	  cout << "\tperm. " << n_perm_max << ": [ ";
	  for( auto ind : perm_index_copy ){
	    if(count==n_tags) cout << " | ";
	    cout << ind << " ";
	    ++count;
	  }
	  cout << "]" << endl;
	}
	++n_perm_max; 
	++(n_perms_max.at(n_tags));
      }
    } 
    while( next_combination( perm_index_copy.begin(), perm_index_copy.begin()+n_tags, perm_index_copy.end(), std::less<int>() ) );
  }

  if( debug_code&DebugVerbosity::init_more ){
    cout << "Maximum of " << n_perm_max << " permutation(s) considered" << endl;
  }

  TH1D h_combinations("h_combinations","", n_perm_max, 0, n_perm_max);
  double pass{0.};
  int count_perm{0};
  for(int n_tags = n_tags_l ; n_tags <= (n_tags_h>=0 ? TMath::Min(n_tags_h,n_jets) : n_jets) ; ++n_tags){
    for( std::size_t n_perm = 0 ; n_perm < n_perms_max.at(n_tags) ; ++n_perm ){        

      double p{1.};
      if( debug_code&DebugVerbosity::event ) cout << "\tPermutation " << n_perm  << endl;
      auto perm = get_permutation(n_perm, n_tags);
      if(perm.size()==0) continue;
      
      for( int j = 0 ; j < n_jets ; ++j){

	if( perm[j]<0 ) continue;

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

      if( debug_code&DebugVerbosity::event ) cout << "\tFilling perm " << count_perm << endl;
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
      if( count_perm==int(rnd_perm)){
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
  

  MEM::Lock lk(n);
  auto perm = get_permutation(n_perm_rnd, n_tags_rnd);
  for(int t = 0 ; t < n_tags_rnd ; ++t){
    if( debug_code&DebugVerbosity::init_more ) cout << "\tEnable lock for jet index (" << perm[t] << ")" << endl;
    lk.set_haslock( size_t(perm[t]), true );
  }

  for(size_t j = 0; j < n ; ++j){  
    if(jets.at(j)->isSet(Observable::IGNORE_FOR_RND) && jets.at(j)->getObs(Observable::IGNORE_FOR_RND)>0 ){
      if( debug_code&DebugVerbosity::init_more ) cout << "\tIgnoring jet index (" << j << ")" << endl;
      //lk.set_haslock( j, true );
      //lk.set_lock   ( j, true );
      lk.set_alwayslock( j, true );
    }
  }

  // restore n_tags_l/n_tags_h and n_tags_rnd
  check_restore_jets(n_tags_rnd);

  int ntoys = 0;
  while(ntoys<n_max_toys){

    ++ntoys;
    if( debug_code&DebugVerbosity::event ) 
      cout << "\tToy: " << ntoys << endl;

    // number of passing jets 
    //int count_pass  {0};

    // number of passing jets among those that have to pass ( specified in perm )
    int count_pass_t{0};

    // number of passing jets among those that DO NOT have to pass ( specified in perm ) 
    int count_pass_u{0};

    for(size_t j = 0; j < n ; ++j){


      if( lk.get_alwayslock(j) ){
	rnd_btag[j] = vals.at(j);
	if( rnd_btag[j]>=cut_val ) ++count_pass_t;
	if( debug_code&DebugVerbosity::event ) 
	  printf("\t\tRandom value jet[%d]: %.3f (excluded)\n", int(j), rnd_btag[j]);
	continue;
      }

      if( lk.get_haslock(j) ){

	if( !lk.get_lock(j) ){ 
	  rnd_btag[j] = pdfs.at(j)->GetRandom();
	  if( rnd_btag[j] >= cut_val ){
	    ++count_pass_t;  
	    lk.set_lock( j, true ); 
	    if( debug_code&DebugVerbosity::event ) 
	      printf("\t\tRandom value jet[%d]: %.3f (has_lock /  locked)\n", int(j), rnd_btag[j]);
	    continue;
	  }
	  if( debug_code&DebugVerbosity::event ) 
	    printf("\t\tRandom value jet[%d]: %.3f (has_lock /        )\n", int(j), rnd_btag[j]);
	  continue;
	}

	++count_pass_t;  
	if( debug_code&DebugVerbosity::event ) 
	  printf("\t\tRandom value jet[%d]: %.3f (has_lock /  locked)\n", int(j), rnd_btag[j]);
	continue;
      }

      if( !lk.get_haslock(j) ){ 
	
	if( !lk.get_lock(j) ){ 
	  rnd_btag[j] = pdfs.at(j)->GetRandom();
	  if( rnd_btag[j] >= cut_val ){
	    ++count_pass_u;  
	    if( debug_code&DebugVerbosity::event ) 
	      printf("\t\tRandom value jet[%d]: %.3f (no_lock  /       )\n", int(j), rnd_btag[j]);
	    continue;
	  }
	  lk.set_lock( j, true ); 
	  if( debug_code&DebugVerbosity::event ) 
	    printf("\t\tRandom value jet[%d]: %.3f (no_lock  /  locked)\n", int(j), rnd_btag[j]);	  
	  continue;
	}

	if( debug_code&DebugVerbosity::event ) 
	  printf("\t\tRandom value jet[%d]: %.3f (no_lock  /  locked)\n", int(j), rnd_btag[j]);
	continue;
      }


      /*
      // needed when I have to unlock all jets
      if(jets.at(j)->isSet(Observable::IGNORE_FOR_RND) && 
	 jets.at(j)->getObs(Observable::IGNORE_FOR_RND)>0 ) lk.set_lock( j, true );

      // check among the locked jets
      if( lk.get_haslock(j) && lk.get_lock(j) ){
	if( debug_code&DebugVerbosity::event ) 
	  printf("\t\tRandom value jet[%d]: %.3f (locked)\n", int(j), rnd_btag[j]);
	if( rnd_btag[j]>cut_val ){
	  ++count_pass;
	  ++count_pass_t;	
	}
	continue;
      }
      
      double rnd  = lk.get_lock(j) ? rnd_btag[j] : pdfs.at(j)->GetRandom();
      rnd_btag[j] = rnd;

      if(  rnd >= cut_val ){
	++count_pass;
	if( lk.get_haslock(j) ){
	  lk.set_lock(j, true); 
	  ++count_pass_t;
	}
	if( !lk.get_haslock(j) ) ++count_pass_u;
      }
      else{
	if( !lk.get_haslock(j) ) lk.set_lock(j, true);
      }
            
      if( debug_code&DebugVerbosity::event )
	printf("\t\tRandom value jet[%d]: %.3f %s\n", int(j), rnd_btag[j], 
	       lk.get_lock(j) ? "(locked)" : "" );      
      */
    }

    if( debug_code&DebugVerbosity::event ) 
      cout << "\tPassing:  in perm...." << count_pass_t << ", out perm...." << count_pass_u << ". Requested...." << n_tags_rnd << endl;

    // condition to pass rnd selection
    if( count_pass_t == n_tags_rnd && count_pass_u==0 ){
      out.pass_rnd = 1;
      break;
    }

    // if the number of passing jets for this toy exceeds the maximum number, reset
    if(count_pass_t>n_tags_rnd){
      if( debug_code&DebugVerbosity::event ) cout << "\tUnlock!" << endl; 
      lk.reset_locked();
    }

  }

  if( out.pass_rnd ){
    out.p          = pass;  
    out.rnd_btag   = rnd_btag;
  }
  out.ntoys      = ntoys;  
  out.err        = error_code;  
  
  return out;
}


void MEM::BTagRandomizer::set_condition(const int& ntags_l, const int& ntags_h, const double& v){
  n_tags_l = ntags_l;
  n_tags_h = ntags_h;
  cut_val  = v;
  if( debug_code&DebugVerbosity::init_more ){
    cout << "BTagRandomizer::set_condition(): add category (" << n_tags_l << "," << n_tags_h << ")" << endl;
  }
}

void MEM::BTagRandomizer::push_back_object( Object* obj){

  jets.push_back( obj );
  
  if( debug_code&DebugVerbosity::init_more ){
    cout << "BTagRandomizer::push_back_object()" << endl;
    obj->print(cout);
  }

  return;
}


bool MEM::BTagRandomizer::consider(const std::vector<int>& perm, const std::size_t& ntag) const {
  for( std::size_t p = 0 ; p < ntag ; ++p ){
    if( perm[p]<0 ) return false;
  }
  return true;
}


std::vector<int> MEM::BTagRandomizer::get_permutation(const std::size_t& n, const int& ntag){
  vector<int> perm_index_copy = perm_index;
  std::size_t n_perm{0};
  do{
    bool is_ok = consider(perm_index_copy, ntag);
    if( is_ok && n==n_perm ){
      if( debug_code&DebugVerbosity::init_more ) {
	int count{0};
	cout << "\tperm. " << n_perm << ": [ ";
	for( auto ind : perm_index_copy ){ 
	  if( count==ntag ) cout << " | ";
	  cout << ind << " ";
	  ++count;
	}
	cout << "]" << endl;
      }
      return perm_index_copy;    
    }
    if( is_ok ) ++n_perm;
  } while( next_combination( perm_index_copy.begin(), perm_index_copy.begin()+ntag, perm_index_copy.end(), std::less<int>()) );

  return vector<int>{};
} 
