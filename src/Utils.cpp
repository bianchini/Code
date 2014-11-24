#include "interface/Utils.h"

size_t Algo::eta_to_bin( const LV& lv ){
  if( fabs(lv.Eta())<1.0 ) return 0;
  if( fabs(lv.Eta())>1.0 ) return 1;
  return -99;
}

string Algo::translateDecay(Algo::Decay& decay){
  
  string name = "";
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
  case Algo::Decay::Higgs:
    name = "Higgs";
    break;
  case Algo::Decay::Radiation_q:
    name = "Radiation_q";
    break;
  case Algo::Decay::Radiation_b:
    name = "Radiation_b";
    break;
  case Algo::Decay::MET:
    name = "MET";
    break;
  default:
    name = "UNKNOWN";
    break;
  }
  
  return name;
}


bool Algo::isSame( const std::vector<std::pair<FinalState,size_t>>& a, const std::vector<std::pair<FinalState,size_t>>& b){

  if(a.size()!=b.size()) return false;

  for(size_t i = 0 ; i < a.size() ; ++i ){

    //printf("A: [%d,%d]\t", int(a[i].first), int(a[i].second)) ;
    //printf("B: [%d,%d]\n", int(b[i].first), int(b[i].second)) ;

    if( !( a[i].first == FinalState::Radiation_q || a[i].first == FinalState::Radiation_b )){
      if( a[i].first != b[i].first || a[i].second != b[i].second ) return false;
    }
    else{
      if( a[i].first != b[i].first ) return false;
    }

  }
  
  return true;

}


//////////////////////////////////////////////////////
