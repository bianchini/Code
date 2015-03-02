#include "interface/Parameters.h"

size_t MEM::eta_to_bin( const double& eta ){
  if( fabs(eta)<1.0 ) return 0;
  if( fabs(eta)>1.0 ) return 1;
  return -99;
}
bool MEM::isQuark(const MEM::TFType& t) {
  return (t==TFType::bReco || t==TFType::qReco || t==TFType::bLost || t==TFType::qLost);
}
bool MEM::isNeutrino(const MEM::TFType& t) {
    return (t==TFType::MET);
  }
bool MEM::isLepton(const MEM::TFType& t)  {
  return (t==TFType::elReco || t==TFType::muReco);
}


/////////////////////////////////
//   y    := observables
//   x    := gen level quantities
//   type := decides the TF
/////////////////////////////////
double MEM::transfer_function(double* y, double* x, const TFType& type){

  // return value
  double w{1.};

  // parameters
  const double* par;

  switch( type ){
    
  case TFType::bReco:
    // x[0] = parton energy ; x[1] = parton eta
    // y[0] = jet energy
    cout << "\t\ttransfer_function: Evaluate W(" << y[0] << " | E=" << x[0] << ", y=" << x[1] << ", TFType::bReco) = ";
    par = TF_B_param[ eta_to_bin(x[1]) ];    
    w *=  par[10]*TMath::Gaus(y[0], par[0] + par[1]*x[0], x[0]*TMath::Sqrt(par[2]*par[2]+par[3]*par[3]/x[0]+par[4]*par[4]/x[0]/x[0]), 1) + 
      (1-par[10])*TMath::Gaus(y[0], par[5] + par[6]*x[0], x[0]*TMath::Sqrt(par[7]*par[7]+par[8]*par[8]/x[0]+par[9]*par[9]/x[0]/x[0]), 1);
    cout << w << endl;
    break;
  case TFType::qReco:
    // x[0] = parton energy ; x[1] = parton eta
    // y[0] = jet energy
    cout << "\t\ttransfer_function: Evaluate W(" << y[0] << " | E=" << x[0] << ", y=" << x[1] << ", TFType::qReco) = ";
    par = TF_Q_param[ eta_to_bin(x[1]) ];
    w *= TMath::Gaus(y[0], par[0] + par[1]*x[0], x[0]*TMath::Sqrt(par[2]*par[2]+par[3]*par[3]/x[0]+par[4]*par[4]/x[0]/x[0]), 1);
    cout << w << endl;
    break;
  case TFType::MET:
    cout << "\t\ttransfer_function: Evaluate W(" << y[0]-x[0] << " , " << y[1]-x[1] << ", TFType::MET) = ";
    // x[0] = sum nu_x ; x[1] = sum nu_y
    // y[0] = MET_x    ; y[1] = MET_y
    par = TF_MET_param;
    w *= TMath::Gaus(y[0], x[0], par[0], 1)*TMath::Gaus(y[1], x[1], par[1], 1);
    cout << w << endl;
    break;
  case TFType::bLost:
  case TFType::qLost:
    // x[0]     = parton energy ; x[1] = parton eta
    // y[0]     = jet energy
    // param[0] = max eta ; param[1] = min pT; param[2] = acceptance    
    par = TF_ACC_param;
    if( TMath::Abs(x[1])>par[0] ){
      cout << "\t\ttransfer_function: Evaluate W(" << x[0] << ", " << x[1] << ", TFType::qLost) = ";
      w = par[2];
      cout << w << endl;
    } 
    else{
      // x[0]     = parton energy ; x[1] = parton eta
      // y[0]     = 0.
      par = TF_Q_param[ eta_to_bin(x[1]) ];
      cout << "\t\ttransfer_function: Evaluate W(" <<  TF_ACC_param[1] << " | " << x[0] << ", " << x[1] << ", TFType::qLost) = ";
      w *= (1 - 0.5*(TMath::Erf(  (TF_ACC_param[1]*TMath::CosH(x[1]) /* need sin */  - (par[0] + par[1]*x[0])) / (x[0]*TMath::Sqrt(par[2]*par[2]+par[3]*par[3]/x[0]+par[4]*par[4]/x[0]/x[0])) ) + 1 ));
      cout << w << endl;
    }
    break;
  default:
    break;
  }

  return w;
}

/////////////////////////////////                                                                                                                      //   y    := observables                                                                                                                               
//   type := decides the TF                                                                                                                          
///////////////////////////////// 
pair<double, double> MEM::get_support(double* y, const TFType& type, const double& alpha){

  double e_rec   = y[0];
  double eta_rec = y[1];
  double e_L     = e_rec*0.25;
  double e_H     = e_rec*4.0;

  size_t steps     = size_t(e_H-e_L)/2; 
  double step_size = (e_H-e_L)/steps;

  /*
  double tot{0.};
  for(size_t step = 0; step < steps ; ++step){
    double gen[2] = {e_rec-step_size*step, eta_rec};    
    for(size_t i = 0; i < 100 ; ++i){
      double rec[1] = {e_rec+i*step_size};
      tot += transfer_function( rec, gen, type)*step_size;
      if(tot>=(1-alpha)/2){  
	e_L = gen[0];
	tot = 0.;
      break;
      }
    }
  }

  for(size_t step = 0; step < steps ; ++step){
    double gen[2] = {e_H-step_size*step, eta_rec};
    tot += transfer_function( rec, gen, type)*step_size;
    if(tot>=(1-alpha)/2){  
      e_H = gen[0];
      break;
    }
  }

  */

  cout << "MEM::get_support: E(reco) = " << e_rec << " ==> range at " << alpha 
       << " CL is [" << e_L << ", " << e_H << "] (stepping every " << step_size << " GeV)" << endl;

  return make_pair(e_L, e_H);
}


MEM::PS::PS(size_t d){
  dim = d;
}

MEM::PS::~PS(){}

map<MEM::PSPart, MEM::GenPart>::const_iterator MEM::PS::begin() const {
  return val.begin();
}

map<MEM::PSPart, MEM::GenPart>::const_iterator MEM::PS::end() const {
  return val.end();
}

LV MEM::PS::lv(const MEM::PSPart& p) const { 
  return val.find(p)!=val.end() ? (val.find(p)->second).lv : LV() ; 
}

MEM::TFType MEM::PS::type(const MEM::PSPart& p) const {
  return val.find(p)!=val.end() ? (val.find(p)->second).type : TFType::Unknown ;
}

void MEM::PS::set(const MEM::PSPart& a, const MEM::GenPart& b){
  val[a] = b;
}

void MEM::PS::print(ostream& os) const{
  for( auto p = val.begin() ; p != val.end() ; ++p ){
    cout << "\tPS[" << static_cast<size_t>(p->first) << "] : type("
	 << static_cast<size_t>(p->second.type) << ") "
	 << p->second.lv.Pt() << ", "
	 << p->second.lv.Eta() << ", "
	 << p->second.lv.Phi() << ", "
	 << p->second.lv.M()
	 << endl; 
  }
}

MEM::Object::Object(const LV& lv, const MEM::ObjectType& t){ 
  p   = lv; 
  type = t; 
} 

MEM::Object::~Object(){}    

LV MEM::Object::p4() const { return p; }

double MEM::Object::getObs(const Observable& name) const { 
  return (obs.find(name)!=obs.end() ? obs.find(name)->second : -99.);
}

void MEM::Object::addObs(const Observable& name, const double& val){ 
  obs.insert( make_pair(name, val) ); 
}

void MEM::Object::print(ostream& os) const {
  os << "\tType: " << static_cast<int>(type) << ", p=(Pt, Eta, Phi, M)=("
     << p.Pt() << ", " << p.Eta() << ", " << p.Phi() << ", " << p.M()
     << ")" << endl;
}

