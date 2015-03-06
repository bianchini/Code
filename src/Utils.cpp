#include "interface/Utils.h"


double Algo::pdf_btag(double* x, double* par){

  // par[0] = flavour: 0=u, 1=d, 2=c, b=3
  // par[1] = pt
  // par[2] = eta

  double val{1.0};
  if(x==nullptr || par==nullptr) return val;

  size_t bin = Algo::eta_to_bin(par[2]);

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
  case Algo::QuarkTypeCharm: // pure C quarks
    if( pt>20. && pt<9999.) 
      val =
	(btag<0.5)*(Algo::BTAG_C_param[bin][0]) + (btag>0.5)*(Algo::BTAG_C_param[bin][1]);
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

string Algo::Decay::translateDecay(Algo::Decay::Decay& decay){
  
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
  case Algo::Decay::Radiation_c:
    name = "Radiation_c";
    break;
  case Algo::Decay::Radiation_b:
    name = "Radiation_b";
    break;
  case Algo::Decay::Radiation_g:
    name = "Radiation_g";
    break;
  case Algo::Decay::Lepton:
    name = "Lepton";
    break;
  case Algo::Decay::MET:
    name = "MET";
    break;    
  case Algo::Decay::UNKNOWN:
  default:
    name = "UNKNOWN";
    break;
  }
  
  return name;
}


int Algo::diff( const std::vector<std::pair<FinalState::FinalState,size_t>>& a, 
		const std::vector<std::pair<FinalState::FinalState,size_t>>& b, 
		const vector<Algo::Object>& jets){
  
  if(a.size()!=b.size())    return -1;
  if(jets.size()!=a.size()) return -1;

  for(size_t i = 0 ; i < a.size() ; ++i ){

    auto fs_a = a[i].first;
    auto fs_b = b[i].first;
    auto id_a = a[i].second;
    auto id_b = b[i].second;

    FinalState::FinalState fs_abar = Algo::partner(fs_a);
    bool btagd  = (jets[i]).isDiscriminating("BTAG");
    bool same   = (id_a==id_b);
    bool matchL = (fs_abar==fs_b);
    bool matchT = (fs_a==fs_b);
    bool israd  = isRadiation(fs_b);

    //matchL = matchT;

    switch(fs_a){
    case FinalState::FinalState::Radiation_g:
    case FinalState::FinalState::Radiation_u:
    case FinalState::FinalState::Radiation_d:
    case FinalState::FinalState::Radiation_c:
    case FinalState::FinalState::Radiation_b:      
      if( btagd && !(matchT || matchL)) return 1;
      if(!btagd && !israd)  return 2;
      continue;
      break;
    case FinalState::FinalState::TopHad_q: 
    case FinalState::FinalState::TopHad_qbar:  
    case FinalState::FinalState::WHad_q:   
    case FinalState::FinalState::WHad_qbar:   
      if(!same)             return 3;
      if( btagd && !matchT) return 4;
      if(!btagd && !(matchL || matchT)) return 5;
      continue;
      break;       
    case FinalState::FinalState::Higgs_b: 
    case FinalState::FinalState::Higgs_bbar: 
      if(!same)                return 6;
      if(!(matchL || matchT) ) return 7;        
      continue;
      break;
    default:
      if( !same || !matchT )   return 8;
      continue;
      break;
    }

    /*
    if( !( a[i].first == FinalState::FinalState::Radiation_u || a[i].first == FinalState::FinalState::Radiation_d || 
	   a[i].first == FinalState::FinalState::Radiation_b || a[i].first == FinalState::FinalState::Radiation_g )){
      if( a[i].first != b[i].first || a[i].second != b[i].second ) return false;
    }
    else{
      if( a[i].first != b[i].first ) return false;
    }
    */
  }
  
  return 0;

}


/*
  Decide whether a permutation assigns quark flavours to correct jets
  (b <-> tagged, u/d/g <-> untagged)
  The decision is taken ONLY if "BTAG" and "BTAG_RND" are filled, and "BTAG_RND"=0
   Input:  a = test, b = target                                                                                                                        
   Output: true if it is a good permutation, false otherwise 
*/
bool Algo::filter_by_btag( const std::vector<std::pair<FinalState::FinalState,size_t>>& particles, const vector<Algo::Object>& jets ){
  
  bool passes {true};

  if( jets.size()==0 ) return passes;

  size_t count {0};
  for( auto p : particles ){
    if( (jets[count].obs).find("BTAG")    ==(jets[count].obs).end() ) continue;
    if( (jets[count].obs).find("BTAG_RND")==(jets[count].obs).end() ) continue;
    if( (jets[count].obs).find("BTAG_RND")!=(jets[count].obs).end() &&
	(jets[count].obs).find("BTAG_RND")->second>0.5 )              continue;

    double btag = (jets[count].obs).find("BTAG")->second;
    switch( p.first ){
    case FinalState::FinalState::TopHad_q:
    case FinalState::FinalState::TopHad_qbar:
    case FinalState::FinalState::WHad_q:
    case FinalState::FinalState::WHad_qbar:
    case FinalState::FinalState::Radiation_u:
    case FinalState::FinalState::Radiation_d:
    case FinalState::FinalState::Radiation_c:
    case FinalState::FinalState::Radiation_g:
      if( btag>0.5 ){
	passes = false;
	return passes;
      }
      break;
    case FinalState::FinalState::TopLep_b:
    case FinalState::FinalState::TopHad_b:
    case FinalState::FinalState::Higgs_b:
    case FinalState::FinalState::Higgs_bbar:
    case FinalState::FinalState::Radiation_b:
      if( btag<0.5 ){
	passes = false;
	return passes;
      }
      break;
    default:
      break;      
    }
    ++count;
  }

  return passes;
}

Algo::FinalState::FinalState Algo::partner( const Algo::FinalState::FinalState fs){
  switch( fs ){
  case FinalState::FinalState::TopHad_q:
    return FinalState::FinalState::TopHad_qbar;
    break;
  case FinalState::FinalState::TopHad_qbar:
    return FinalState::FinalState::TopHad_q;
    break;
  case FinalState::FinalState::WHad_q:
    return FinalState::FinalState::WHad_qbar;
    break;
  case FinalState::FinalState::WHad_qbar:
    return FinalState::FinalState::WHad_q;
    break;
  case FinalState::FinalState::Higgs_b:
    return FinalState::FinalState::Higgs_bbar;
    break;
  case FinalState::FinalState::Higgs_bbar:
    return FinalState::FinalState::Higgs_b;
    break;
  default:
    return fs;
    break;
  }
  return fs;
}

bool Algo::isRadiation( const Algo::FinalState::FinalState fs ){
  return (fs==FinalState::FinalState::Radiation_u || fs==FinalState::FinalState::Radiation_d || fs==FinalState::FinalState::Radiation_c || fs==FinalState::FinalState::Radiation_b || fs==FinalState::FinalState::Radiation_g);
}


//////////////////////////////////////////////////////
