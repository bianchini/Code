#include "interface/Utils.h"


double Algo::pdf_btag(double* x, double* par){

  // par[0] = flavour: 0=u, 1=d, 2=b
  // par[1] = pt
  // par[2] = eta

  double val{1.0};
  if(x==nullptr || par==nullptr) return val;

  size_t bin = Algo::eta_to_bin(par[2]);
  if(bin<0) return val;

  Algo::QuarkType type = static_cast<Algo::QuarkType>(static_cast<int>(par[0])); 
  double btag = x[0];
  double pt   = par[1];

  switch(type){
  case Algo::QuarkTypeUp: // 50% mixture of C and D quarks
    if( pt>20. && pt<9999.) 
      val = 
	(btag<0.5)*(Algo::BTAG_Q_param[bin][0]*0.5 + Algo::BTAG_C_param[bin][0]*0.5)  +
	(btag>0.5)*(Algo::BTAG_Q_param[bin][1]*0.5 + Algo::BTAG_C_param[bin][1]*0.5);
    break;
  case Algo::QuarkTypeDown: // pure D quarks
    if( pt>20. && pt<9999.) 
      val =
	(btag<0.5)*(Algo::BTAG_Q_param[bin][0]) + (btag>0.5)*(Algo::BTAG_Q_param[bin][1]);
    break;
  case Algo::QuarkTypeBottom: // pure B quarks
    if( pt>20. && pt<9999.) 
      val =
	(btag<0.5)*(Algo::BTAG_B_param[bin][0]) + (btag>0.5)*(Algo::BTAG_B_param[bin][1]);
    break;
  default:
    break;
  }

  return val;
}

size_t Algo::eta_to_bin( const double& eta ){
  if( fabs(eta)<1.0 ) return 0;
  if( fabs(eta)>1.0 ) return 1;
  return -99;
}

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
  case Algo::Decay::Radiation_u:
    name = "Radiation_u";
    break;
  case Algo::Decay::Radiation_d:
    name = "Radiation_d";  
    break;
  case Algo::Decay::Radiation_b:
    name = "Radiation_b";
    break;
  case Algo::Decay::Radiation_g:
    name = "Radiation_g";
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

    if( !( a[i].first == FinalState::Radiation_u || a[i].first == FinalState::Radiation_d || a[i].first == FinalState::Radiation_b || a[i].first == FinalState::Radiation_g )){
      if( a[i].first != b[i].first || a[i].second != b[i].second ) return false;
    }
    else{
      if( a[i].first != b[i].first ) return false;
    }

  }
  
  return true;

}


//////////////////////////////////////////////////////
