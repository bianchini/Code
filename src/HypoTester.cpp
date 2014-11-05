#include "interface/HypoTester.h"


Algo::HypoTester::HypoTester(){

  // reset variables
  count_perm        = 0;
  count_TopHad      = 0;
  count_TopLep      = 0;
  count_HiggsHad    = 0;
  count_Radiation   = 0;
}


Algo::HypoTester::~HypoTester(){
  cout << "Removing HypoTester" << endl;
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
     break;

   case Algo::Decay::TopLep:
     particles.push_back( make_pair( FinalState::TopLep_b,    count_TopLep) );
     ++count_TopLep;
     break;

   case Algo::Decay::HiggsHad:
     particles.push_back( make_pair( FinalState::HiggsHad_b,    count_HiggsHad) );
     particles.push_back( make_pair( FinalState::HiggsHad_bbar, count_HiggsHad) );
     ++count_HiggsHad;
     break;

   case Algo::Decay::Radiation:
     particles.push_back( make_pair( FinalState::Radiation_q,   count_Radiation) );
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


void Algo::HypoTester::run(){

  // first, unpack assumptions
  this->unpack_assumptions();

  assert( particles.size()<=p4_Jet.size() );
  assert( count_TopLep==p4_Lepton.size()  );
  
  sort( particles.begin(), particles.end(), MyComp );
  count_perm = 0;
  do {
    cout << count_perm << "th perm: [ " ; 
    for( auto p : particles )
      cout << "(" << p.first << "," << p.second << ") " ;
    cout << "]" << endl;

    group_particles();

    ++count_perm;
  } while ( next_permutation(particles.begin(), particles.end(), MyComp  ) );
  
  cout << "Run finished with " << count_perm << " permutations" << endl;
}


void Algo::HypoTester::group_particles(){

  vector<pair<FinalState,int>> this_permut;

  // get a copy
  this_permut.assign( particles.begin(), particles.end() );

  // first get all hadronically decaying tops
  for( size_t t_had = 0; t_had < count_TopHad; ++t_had ){
    /*...*/
  }

}


void Algo::HypoTester::TopHad::init( const LV& a, const LV& b , const LV& c){
  p4_q    = a; 
  p4_qbar = b; 
  p4_b    = c;
}

double Algo::HypoTester::TopHad::get_e_qbar(const double& e_q){
  double val{0.};
  /* do somewthing */
  return val;
}

double Algo::HypoTester::TopHad::get_e_b(const double& e_q){
  double val{0.};
  /* do somewthing */
  return val;
}


string Algo::translateDecay(Algo::Decay& decay){

    string name = "Unknown";
    switch( decay ){
    case Algo::Decay::TopHad:
      name = "TopHad";
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

