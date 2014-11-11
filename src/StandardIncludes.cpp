#include "interface/StandardIncludes.h"


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
  case Algo::Decay::HiggsHad:
    name = "HiggsHad";
    break;
  case Algo::Decay::Radiation:
    name = "Radiation";
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


//////////////////////////////////////////////////////


Algo::TransferFunction::TransferFunction(const string& name, const string& form){
  formula = form;
  f       = new TFormula(name.c_str(), formula.c_str());
}

Algo::TransferFunction::~TransferFunction(){
  if(VERBOSE) cout << "Destroy tf " << string(f->GetName()) << endl;
  delete f;
}

const string Algo::TransferFunction::getFormula() const {
  return string(f->GetExpFormula("p"));
}

double Algo::TransferFunction::eval(const double& rec, const double& gen) const {
  return f->Eval(rec,gen);
}

void Algo::TransferFunction::init(const double* param){
  assert(f!=nullptr);
  f->SetParameters( param );
}


//////////////////////////////////////////////////////

Algo::CombBuilder::CombBuilder() {

}

Algo::CombBuilder::CombBuilder(vector<DecayBuilder*>& comb) {
  combined = comb;
}


Algo::CombBuilder::~CombBuilder() {
  if(VERBOSE) cout << "Destroy CombBuilder" << endl;
  for(auto dec : combined) 
    delete dec;
}


void Algo::CombBuilder::add(DecayBuilder* dec) {
  combined.push_back(dec);
}


double Algo::CombBuilder::eval(const double* xx) {

  double val{1.};

  LV invisible(0.,0.,0.,0.);
  for(auto dec : combined) 
    val *= dec->eval(xx, invisible);

  return val;
}

double Algo::CombBuilder::eval(const double* xx,  LV& lv) {
  return eval( xx ); 
}


void Algo::CombBuilder::print(ostream& os){
  os << "\t\tCombBuilder contains " << combined.size() << " block(s):" << endl; 
  for(auto dec : combined) dec->print(os);
}

//////////////////////////////////////////////////////


Algo::METBuilder::METBuilder() {
  decay = Decay::MET;
  tf_met = nullptr;
}


Algo::METBuilder::~METBuilder() {
   if(VERBOSE) cout << "Destroy METBuilder" << endl;
  if( tf_met!=nullptr) delete tf_met;
}


void Algo::METBuilder::init(const LV& lv) {
  p4_invisible = lv;
  tf_met = new TransferFunction("tf_met", TF_MET );
}


double Algo::METBuilder::eval(const double* xx,  LV& lv) {

  double val{1.};

  val *= tf_met->eval( p4_invisible.Px()-lv.Px(),p4_invisible.Py()-lv.Py() );
  return val;
}


void Algo::METBuilder::print(ostream& os){
  os << "\t\t\tDecay: " <<  Algo::translateDecay(decay) <<  endl; 
  os << "\t\t\tMET:   p4 = (" << p4_invisible.Pt() << ", " << p4_invisible.Phi() << ")" << endl;
  if(tf_met!=nullptr)  os << "\t\t\tTF MET   : " << tf_met->getFormula() << endl;
}


//////////////////////////////////////////////////////

Algo::TopHadBuilder::TopHadBuilder () {
  decay   = Decay::TopHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
  tf_b    = nullptr;
}

Algo::TopHadBuilder::~TopHadBuilder() {
   if(VERBOSE) cout << "Destroy TopHadBuilder" << endl;
  if( tf_q   != nullptr ) delete tf_q;
  if( tf_qbar!= nullptr ) delete tf_qbar;
  if( tf_b   != nullptr ) delete tf_b;
}

void Algo::TopHadBuilder::print(ostream& os){
  os << "\t\t\tDecay: " <<  Algo::translateDecay(decay) << endl; 
  os << "\t\t\tq    [" << index_q    << "]: p4 = (" << p4_q.Pt() << ", " << p4_q.Eta()  << ", " << p4_q.Phi() << ", " << p4_q.M() << ")" << endl; 
  os << "\t\t\tqbar [" << index_qbar << "]: p4 = (" << p4_qbar.Pt() << ", " << p4_qbar.Eta()  << ", " << p4_qbar.Phi() << ", " << p4_qbar.M() << ")" <<endl; 
  os << "\t\t\tb    [" << index_b    << "]: p4 = (" << p4_b.Pt() << ", " << p4_b.Eta()  << ", " << p4_b.Phi() << ", " << p4_b.M() << ")" <<endl; 
  if(tf_q!=nullptr)    os << "\t\t\tTF q   : " << tf_q->getFormula() << endl;
  if(tf_qbar!=nullptr) os << "\t\t\tTF qbar: " << tf_qbar->getFormula() << endl;
  if(tf_b!=nullptr)    os << "\t\t\tTF b   : " << tf_b->getFormula() << endl;
}



void Algo::TopHadBuilder::init( const FinalState& fs, const LV& lv, const size_t& sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::TopHad_q:
    p4_q    = lv;
    index_q = sz;
    tf = new TransferFunction("tf_q", TF_Q);
    tf->init( TF_Q_param[Algo::eta_to_bin(lv)] );
    tf_q    = tf;
    break;
  case FinalState::TopHad_qbar:
    p4_qbar    = lv;
    index_qbar = sz;
    tf = new TransferFunction("tf_qbar", TF_Q);
    tf->init( TF_Q_param[Algo::eta_to_bin(lv)] );
    tf_qbar    = tf;
    break;
  case FinalState::TopHad_b:
    p4_b    = lv;
    index_b = sz;
    tf = new TransferFunction("tf_b", TF_B);
    tf->init( TF_Q_param[Algo::eta_to_bin(lv)] );
    tf_b    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::TopHadBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val {0.};

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

Algo::WHadBuilder::WHadBuilder () {
  decay   = Decay::WHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
};

Algo::WHadBuilder::~WHadBuilder() {
   if(VERBOSE) cout << "Destroy WHadBuilder" << endl;
  if( tf_q   != nullptr ) delete tf_q;
  if( tf_qbar!= nullptr ) delete tf_qbar;
};

void Algo::WHadBuilder::print(ostream& os){
  os << "\t\t\tDecay: " << Algo::translateDecay(decay) << endl; 
  os << "\t\t\tq    [" << index_q    << "]: p4 = (" << p4_q.Pt() << ", " << p4_q.Eta()  << ", " << p4_q.Phi() << ", " << p4_q.M() << ")" << endl; 
  os << "\t\t\tqbar [" << index_qbar << "]: p4 = (" << p4_qbar.Pt() << ", " << p4_qbar.Eta()  << ", " << p4_qbar.Phi() << ", " << p4_qbar.M() << ")" <<endl; 
  if(tf_q!=nullptr)    os << "\t\t\tTF q   : " << tf_q->getFormula() << endl;
  if(tf_qbar!=nullptr) os << "\t\t\tTF qbar: " << tf_qbar->getFormula() << endl;
}



void Algo::WHadBuilder::init( const FinalState& fs, const LV& lv, const size_t& sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::WHad_q:
    p4_q    = lv;
    index_q = sz;
    tf = new TransferFunction("tf_q", TF_Q);
    tf->init( TF_Q_param[Algo::eta_to_bin(lv)] );
    tf_q    = tf;
    break;
  case FinalState::WHad_qbar:
    p4_qbar    = lv;
    index_qbar = sz;
    tf = new TransferFunction("tf_qbar", TF_Q);
    tf->init( TF_Q_param[Algo::eta_to_bin(lv)] );
    tf_qbar    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::WHadBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val {0.};

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


//////////////////////////////////////////////////////

Algo::TopLepBuilder::TopLepBuilder () {
  decay   = Decay::TopLep;
  errFlag = 0;
  tf_b    = nullptr;
}

Algo::TopLepBuilder::~TopLepBuilder() {
   if(VERBOSE) cout << "Destroy TopLepBuilder" << endl;
  if( tf_b   != nullptr ) delete tf_b;
}

void Algo::TopLepBuilder::print(ostream& os){
  os << "\t\t\tDecay: " <<  Algo::translateDecay(decay) << endl; 
  os << "\t\t\tlep  [" << index_l    << "]: p4 = (" << p4_l.Pt() << ", " << p4_l.Eta()  << ", " << p4_l.Phi() << ", " << p4_l.M() << ")" <<endl; 
  os << "\t\t\tb    [" << index_b    << "]: p4 = (" << p4_b.Pt() << ", " << p4_b.Eta()  << ", " << p4_b.Phi() << ", " << p4_b.M() << ")" <<endl; 
  if(tf_b!=nullptr)    os << "\t\t\tTF b   : " << tf_b->getFormula() << endl;
}



void Algo::TopLepBuilder::init( const FinalState& fs, const LV& lv, const size_t& sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::TopLep_l:
    p4_l    = lv;
    index_l = sz;
    break;
  case FinalState::TopLep_b:
    p4_b    = lv;
    index_b = sz;
    tf = new TransferFunction("tf_b", TF_B);
    tf->init( TF_Q_param[Algo::eta_to_bin(lv)] );
    tf_b    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::TopLepBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val {0.};

  int nSol {0};

  double Elep     = p4_l.E();
  double nuPhi    = xx[ index_l   ] ; // fill here
  double nuTheta  = xx[ index_l +1] ; // fill here
  
  double Enu, Eb;

  TVector3 e3(0.,0.,1.); // neutrino
  e3.SetTheta( TMath::ACos( nuTheta ) );
  e3.SetPhi ( nuPhi);
  e3.SetMag ( 1.);
  
  double a12 = p4_l.Angle(p4_b.Vect()); // lep - b
  double a13 = p4_l.Angle(e3);          // lep - nu
  double a23 = p4_b.Angle(e3);          // b - nu
  
  Enu = MW*MW/ Elep / (4*TMath::Sin(a13/2.)*TMath::Sin(a13/2.));
  double a = Elep + Enu;
  double b = Elep*TMath::Cos(a12)+Enu*TMath::Cos(a23);
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
    Eb = TMath::Max( E3tmp1,E3tmp2 );
  else if( E3tmp1>0 && E3tmp2<0)
    Eb = E3tmp1;
  else if( E3tmp1<0 && E3tmp2>0)
    Eb = E3tmp2;
  else{
    errFlag = 1;
    return val;
  }
  if(Eb<MB){
    errFlag = 1;
    return val;
  }
  
  LV wNu ( (e3.Unit())*Enu, Enu);
  invisible += wNu;
  
  val += tf_b->eval(p4_b.E(), Eb);
  
  return val;
  
  /*
    TLorentzVector wLep( p4_l.Vect().Unit()*Elep, Elep);
    TLorentzVector blv ( p4_b.Vect().Unit() *(TMath::Sqrt(Eb*Eb - Mb_*Mb_)), Eb);
  */
  
}


//////////////////////////////////////////////////////
