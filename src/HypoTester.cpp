#include "interface/HypoTester.h"
#include "interface/StandardIncludes.h"



string Algo::translateDecay(Algo::Decay& decay){

    string name = "Unknown";
    switch( decay ){
    case Algo::Decay::TopHad:
      name = "TopHad";
      break;
    case Algo::Decay::WHad:
      name = "WHad";
      break;
    case Algo::Decay::TopLep:
      name = "TopLep";
      break;
    case Algo::Decay::HiggsHad:
      name = "HiggsHad";
      break;
    default:
      name = "Radiation";
      break;
    }
    
    return name;
  }


Algo::HypoTester::HypoTester(){

  // reset variables
  count_perm        = 0;
  count_TopHad      = 0;
  count_TopLep      = 0;
  count_HiggsHad    = 0;
  count_Radiation   = 0;

  verbose = 2;

  minimizer =  ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
  minimizer->SetMaxFunctionCalls(1000000); 
  minimizer->SetMaxIterations(10000);  
  minimizer->SetTolerance(0.001);
  minimizer->SetPrintLevel(0);

  tf_met    = nullptr;
}


Algo::HypoTester::~HypoTester(){
  cout << "Removing HypoTester" << endl;
  if(tf_met != nullptr) delete tf_met;
  delete minimizer;
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

void Algo::HypoTester::print(ostream& os){

  os << "****************************************************" << endl;
  os << "Content of this HypoTester:" << endl;
  os << " -- Jets:" << endl;
  int count_j = 0;
  for( auto jet : p4_Jet ){
    os << "\tjet[" << count_j << "]: (" << jet.p4.Pt() << "," << jet.p4.Eta() << "," << jet.p4.Phi() << "," << jet.p4.M() << ")" << endl;
    for( auto it = (jet.obs).begin() ; it !=  (jet.obs).end(); ++it)
      os << "\t\t" << it->first << " = " << it->second << endl;
    ++count_j;
  }

  os << " -- Leptons:" << endl;
  int count_l = 0;
  for( auto lepton : p4_Lepton ){
    os << "\tlepton[" << count_l << "]: (" << lepton.p4.Pt() << "," << lepton.p4.Eta() << "," << lepton.p4.Phi() << "," << lepton.p4.M() << ")" << endl;
    for( auto it = (lepton.obs).begin() ; it !=  (lepton.obs).end(); ++it)
      os << "\t\t" << it->first << " = " << it->second << endl;
    ++count_l;
  }

  os << " -- MET:" << endl;
  int count_m = 0;
  for( auto MET : p4_MET ){
    os << "\tMET[" << count_l << "]: (" << MET.p4.Pt() << "," << MET.p4.Eta() << "," << MET.p4.Phi() << "," << MET.p4.M() << ")" << endl;
    for( auto it = (MET.obs).begin() ; it !=  (MET.obs).end(); ++it)
      os << "\t\t" << it->first << " = " << it->second << endl;
    ++count_m;
  }

  os << " -- Assuming the following decays:" << endl;
  for( auto decay : decays )
    os << "\t" << int(decay) << " (=" << Algo::translateDecay(decay)  << ")" << endl;
  os << " -- Assumed partons:" << endl;
  int count_all = 0;  
  for( auto particle : particles ){
    os << "\t" << particle.first << " (x" << particle.second << ")" << endl;
    ++count_all;
  }
  os << "\tTotal = " << count_all << " partons" << endl;

  os << "****************************************************" << endl;
}

void Algo::HypoTester::assume( Decay decay ){
  decays.push_back( decay );
}


void Algo::HypoTester::unpack_assumptions(){

  // reset
  count_TopHad      = 0;
  count_WHad      = 0;
  count_TopLep      = 0;
  count_HiggsHad    = 0;
  count_Radiation   = 0;

  for( auto decay : decays ){

   switch( decay ){

   case Algo::Decay::TopHad:
     particles.push_back( make_pair( FinalState::TopHad_q,    count_TopHad) );
     particles.push_back( make_pair( FinalState::TopHad_qbar, count_TopHad) );
     particles.push_back( make_pair( FinalState::TopHad_b,    count_TopHad) );
     ++count_TopHad;
     if(verbose>0){ cout << "Added TopHad" << endl; }
     break;

   case Algo::Decay::WHad:
     particles.push_back( make_pair( FinalState::WHad_q,    count_WHad) );
     particles.push_back( make_pair( FinalState::WHad_qbar, count_WHad) );
     ++count_WHad;
     if(verbose>0){ cout << "Added WHad" << endl; }
     break;

   case Algo::Decay::TopLep:
     particles.push_back( make_pair( FinalState::TopLep_b,    count_TopLep) );
     ++count_TopLep;
     if(verbose>0){ cout << "Added TopLep" << endl; }
     break;

   case Algo::Decay::HiggsHad:
     particles.push_back( make_pair( FinalState::HiggsHad_b,    count_HiggsHad) );
     particles.push_back( make_pair( FinalState::HiggsHad_bbar, count_HiggsHad) );
     ++count_HiggsHad;
     if(verbose>0){ cout << "Added HiggsHad" << endl; }
     break;

   case Algo::Decay::Radiation:
     particles.push_back( make_pair( FinalState::Radiation_q,   count_Radiation) );
     if(verbose>0){ cout << "Added Radiation" << endl; }
     ++count_Radiation;
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
      ++count_Radiation;
      --addradiation;
    }

  }
    

}


void Algo::HypoTester::read(){

  // first, unpack assumptions
  this->unpack_assumptions();

  // create TF for invisible (once for good)
  this->create_tf_met();

  assert( particles.size()<=p4_Jet.size() );
  assert( count_TopLep==p4_Lepton.size()  );
  
  sort( particles.begin(), particles.end(), MyComp );
}


double Algo::HypoTester::eval(const double* xx){

  double val {0.};

  count_perm = 0;

  do {

    double val_perm{1.};

    if(verbose>1){
      cout << count_perm << "th perm: [ " ; 
      for( auto p : particles )
	cout << "(" << p.first << "," << p.second << ") " ;
      cout << "]" << endl;
    }
    
    vector<DecayBuilder*> decayed; 
    group_particles( decayed );
    
    LV invisible;
    LV MET = p4_MET.size() ? p4_MET[0].p4 : LV();
    for( auto dec : decayed ){  
      if( verbose > 0) dec->print(cout) ;
      double val_perm_tmp = dec->eval( xx ,invisible );
      val_perm *= val_perm_tmp;
    }
    if(invisible.E()>1e-02){
      if( verbose>1 ) cout << "Evaluating tf_met..." << endl;
      val_perm *= tf_met->eval( invisible.Px()-MET.Px(),invisible.Py()-MET.Py() );
    }
    
    val += val_perm;

    ++count_perm;

    // clean
    for( auto dec : decayed ){      
      delete dec;
    }

  } while ( next_permutation(particles.begin(), particles.end(), MyComp  ) );
  
  if( TMath::IsNaN(val) )
    return numeric_limits<double>::min();
  else if( val<=0 )
    return numeric_limits<double>::min();
  else 
    return -TMath::Log(count_perm>0 ? val/count_perm : val);
  
  return val;

  cout << "Run finished with " << count_perm << " permutations" << endl;
}


void Algo::HypoTester::group_particles(vector<DecayBuilder*>& decayed){

  if(verbose>2){ cout << "Start grouping" << endl; }

  // first get all hadronically decaying tops
  for( size_t t_had = 0; t_had < count_TopHad; ++t_had ){

    if(verbose>1) cout << "Processing " << t_had << "th TopHad" << endl;

    TopHadBuilder* topHad = new TopHadBuilder();
    
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second != t_had ) continue;
      
      // no effect if the type does not match
      topHad->init( part.first, p4_Jet[pos].p4 , pos );
      ++pos;
    }

    decayed.push_back( topHad );

  }

  // first get all hadronically decaying tops
  for( size_t w_had = 0; w_had < count_WHad; ++w_had ){

    if(verbose>1) cout << "Processing " << w_had << "th WHad" << endl;

    WHadBuilder* wHad = new WHadBuilder();
    
    size_t pos = 0;
    for( auto part : particles ){
      if( part.second != w_had ) continue;
      
      // no effect if the type does not match
      wHad->init( part.first, p4_Jet[pos].p4 , pos );
      ++pos;
    }

    decayed.push_back( wHad );

  }


  /* ... */  

}


void Algo::HypoTester::run(){

  ROOT::Math::Functor f0( this , &Algo::HypoTester::eval, 2); 
  minimizer->SetFunction(f0);

  minimizer->SetVariable(0,"j1", p4_Jet[0].p4.E() , 1.);
  minimizer->SetVariable(1,"j2", p4_Jet[1].p4.E() , 1.);

  verbose = 0;
  minimizer->Minimize(); 
  ++verbose;

  double min0 = minimizer->MinValue();

  const double *xs = minimizer->X();
  if(verbose>0)  
    cout << "Minimum --> f(" << xs[0] << "," << xs[1] << ") => " << min0  << endl;

}


void Algo::HypoTester::create_tf_met(){  
  tf_met = new TransferFunction("tf_met", TF_MET );
}


//////////////////////////////////////////////////////

Algo::HypoTester::TransferFunction::TransferFunction(const string& name, const string& form){
  formula = form;
  f       = new TFormula(name.c_str(), formula.c_str());
}

Algo::HypoTester::TransferFunction::~TransferFunction(){
  //cout << "Destroy tf " << string(f->GetName()) << endl;
  delete f;
}

const string Algo::HypoTester::TransferFunction::getFormula() const {
  return string(f->GetExpFormula("p"));
}

double Algo::HypoTester::TransferFunction::eval(const double& rec, const double& gen) const {
  return f->Eval(rec,gen);
}

void Algo::HypoTester::TransferFunction::init(const double* param){
  assert(f!=nullptr);
  f->SetParameters( param );
}






//////////////////////////////////////////////////////

Algo::HypoTester::TopHadBuilder::TopHadBuilder () {
  decay   = Decay::TopHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
  tf_b    = nullptr;
};

Algo::HypoTester::TopHadBuilder::~TopHadBuilder() {
  //cout << "Destroy TopHadBuilder" << endl;
  if( tf_q   != nullptr ) delete tf_q;
  if( tf_qbar!= nullptr ) delete tf_qbar;
  if( tf_b   != nullptr ) delete tf_b;
};

void Algo::HypoTester::TopHadBuilder::print(ostream& os){
  os << "\tDecay: " <<  Algo::translateDecay(decay) << endl; 
  os << "\tq    [" << index_q    << "]: p4 = (" << p4_q.Pt() << ", " << p4_q.Eta()  << ", " << p4_q.Phi() << ", " << p4_q.M() << ")" << endl; 
  os << "\tqbar [" << index_qbar << "]: p4 = (" << p4_qbar.Pt() << ", " << p4_qbar.Eta()  << ", " << p4_qbar.Phi() << ", " << p4_qbar.M() << ")" <<endl; 
  os << "\tb    [" << index_b    << "]: p4 = (" << p4_b.Pt() << ", " << p4_b.Eta()  << ", " << p4_b.Phi() << ", " << p4_b.M() << ")" <<endl; 
  if(tf_q!=nullptr)    os << "\tTF q   : " << tf_q->getFormula() << endl;
  if(tf_qbar!=nullptr) os << "\tTF qbar: " << tf_qbar->getFormula() << endl;
  if(tf_b!=nullptr)    os << "\tTF b   : " << tf_b->getFormula() << endl;
}



void Algo::HypoTester::TopHadBuilder::init( const FinalState& fs, const LV& lv, const size_t& sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::TopHad_q:
    p4_q    = lv;
    index_q = sz;
    tf = new TransferFunction("tf_q", TF_Q);
    tf->init( TF_Q_param[eta_to_bin(lv)] );
    tf_q    = tf;
    break;
  case FinalState::TopHad_qbar:
    p4_qbar    = lv;
    index_qbar = sz;
    tf = new TransferFunction("tf_qbar", TF_Q);
    tf->init( TF_Q_param[eta_to_bin(lv)] );
    tf_qbar    = tf;
    break;
  case FinalState::TopHad_b:
    p4_b    = lv;
    index_b = sz;
    tf = new TransferFunction("tf_b", TF_B);
    tf->init( TF_Q_param[eta_to_bin(lv)] );
    tf_b    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::HypoTester::TopHadBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val = 0.;

  double E1 = xx[ index_q ];
  double E2, E3;

  int nSol = 0;

  double a12 = p4_q.Angle(p4_qbar.Vect());
  double a13 = p4_q.Angle(p4_b.Vect());
  double a23 = p4_qbar.Angle(p4_b.Vect());
  
  E2 = MW*MW/E1/(4*TMath::Sin(a12/2.)*TMath::Sin(a12/2.));

  double a = E1+E2;
  double b = E1*TMath::Cos(a13)+E2*TMath::Cos(a23);
  if( (DM2*DM2 - (a*a - b*b)*MB*MB) < 0){
    errFlag = 1;
    return val;
  }

  double E3_1 = (a*DM2 + b*TMath::Sqrt(DM2*DM2 - (a*a - b*b)*MB*MB))/(a*a - b*b) ;
  double E3_2 = (a*DM2 - b*TMath::Sqrt(DM2*DM2 - (a*a - b*b)*MB*MB))/(a*a - b*b) ;
  double E3tmp1 = -999.;
  double E3tmp2 = -999.;

  if( b>0 ){
    if(E3_1>DM2/a){
      E3tmp1 = E3_1;
      nSol++;
    }
    if(E3_2>DM2/a){
      E3tmp2 = E3_2;
      nSol++;
    }
  }
  else{
    if(E3_1<DM2/a){
      E3tmp1 = E3_1;
      nSol++;
    }
    if(E3_2<DM2/a){
      E3tmp2 = E3_2;
      nSol++;
    }
  }
  if( E3tmp1>0 && E3tmp2>0 )
    E3 = TMath::Max( E3tmp1,E3tmp2 );
  else if( E3tmp1>0 && E3tmp2<0)
    E3 = E3tmp1;
  else if( E3tmp1<0 && E3tmp2>0)
    E3 = E3tmp2;
  else{
    errFlag = 1;
    return val;
  }
  if(E3<MB){
    errFlag = 1;
    return val;
  }

  invisible += LV(0.,0.,0.,0.);

  val += tf_q->eval(p4_q.E(), E1) * tf_qbar->eval(p4_qbar.E(), E2) * tf_b->eval(p4_b.E(), E3);

  return val;

  /*
  TLorentzVector w1 ( p4_q.Vect().Unit()*E1, E1);
  TLorentzVector w2 ( p4_qbar.Vect().Unit()*E2, E2);
  TLorentzVector blv( p4_b.Vect().Unit()*(TMath::Sqrt(E3*E3 - MB*MB)), E3);
  */
  
}


//////////////////////////////////////////////////////

Algo::HypoTester::WHadBuilder::WHadBuilder () {
  decay   = Decay::WHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
};

Algo::HypoTester::WHadBuilder::~WHadBuilder() {
  //cout << "Destroy WHadBuilder" << endl;
  if( tf_q   != nullptr ) delete tf_q;
  if( tf_qbar!= nullptr ) delete tf_qbar;
};

void Algo::HypoTester::WHadBuilder::print(ostream& os){
  os << "\tDecay: " << Algo::translateDecay(decay) << endl; 
  os << "\tq    [" << index_q    << "]: p4 = (" << p4_q.Pt() << ", " << p4_q.Eta()  << ", " << p4_q.Phi() << ", " << p4_q.M() << ")" << endl; 
  os << "\tqbar [" << index_qbar << "]: p4 = (" << p4_qbar.Pt() << ", " << p4_qbar.Eta()  << ", " << p4_qbar.Phi() << ", " << p4_qbar.M() << ")" <<endl; 
  if(tf_q!=nullptr)    os << "\tTF q   : " << tf_q->getFormula() << endl;
  if(tf_qbar!=nullptr) os << "\tTF qbar: " << tf_qbar->getFormula() << endl;
}



void Algo::HypoTester::WHadBuilder::init( const FinalState& fs, const LV& lv, const size_t& sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::WHad_q:
    p4_q    = lv;
    index_q = sz;
    tf = new TransferFunction("tf_q", TF_Q);
    tf->init( TF_Q_param[eta_to_bin(lv)] );
    tf_q    = tf;
    break;
  case FinalState::WHad_qbar:
    p4_qbar    = lv;
    index_qbar = sz;
    tf = new TransferFunction("tf_qbar", TF_Q);
    tf->init( TF_Q_param[eta_to_bin(lv)] );
    tf_qbar    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::HypoTester::WHadBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val = 0.;

  double E1 = xx[ index_q ];
  double E2;

  double a12 = p4_q.Angle(p4_qbar.Vect());
  
  E2 = MW*MW/E1/(4*TMath::Sin(a12/2.)*TMath::Sin(a12/2.));

  invisible += LV(0.,0.,0.,0.);

  val += tf_q->eval(p4_q.E(), E1) * tf_qbar->eval(p4_qbar.E(), E2) ;

  return val;

  /*
  TLorentzVector w1 ( p4_q.Vect().Unit()*E1, E1);
  TLorentzVector w2 ( p4_qbar.Vect().Unit()*E2, E2);
  TLorentzVector blv( p4_b.Vect().Unit()*(TMath::Sqrt(E3*E3 - MB*MB)), E3);
  */
  
}

