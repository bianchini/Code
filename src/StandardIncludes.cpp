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
  case Algo::Decay::Higgs:
    name = "Higgs";
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


bool Algo::isSame( const std::vector<std::pair<FinalState,size_t>>& a, const std::vector<std::pair<FinalState,size_t>>& b){

  if(a.size()!=b.size()) return false;

  for(size_t i = 0 ; i < a.size() ; ++i ){

    //printf("A: [%d,%d]\t", int(a[i].first), int(a[i].second)) ;
    //printf("B: [%d,%d]\n", int(b[i].first), int(b[i].second)) ;

    if( a[i].first != FinalState::Radiation_q ){
      if( a[i].first != b[i].first || a[i].second != b[i].second ) return false;
    }
    if( a[i].first == FinalState::Radiation_q ){
      if( a[i].first != b[i].first ) return false;
    }

  }
  
  return true;

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
  verbose = 0;
}

Algo::CombBuilder::CombBuilder(vector<DecayBuilder*>& comb) {
  combined = comb;
  verbose  = 0;
}

Algo::CombBuilder::CombBuilder(vector<DecayBuilder*>& comb, const int& verb) {
  combined = comb;
  verbose  = verb;
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

  int count {0};  
  for(auto dec : combined){
    double val_tmp = dec->eval(xx, invisible);
    val *= val_tmp;
    if(verbose){
      cout << "\tEval block " << count << " += " << (val_tmp>0. ? -TMath::Log( val_tmp ) : 0.) ;
      if( invisible.Px()>0. || invisible.Py()>0.) 
	cout << ", Invisible = " << "(eng=" << invisible.E() << ", eta=" << invisible.Eta() << ", phi=" << invisible.Phi() << ")" << endl;
      else
	cout << endl;
    }
    ++count;
  }

  return val;
}

double Algo::CombBuilder::eval(const double* xx,  LV& lv) {
  return eval( xx ); 
}


void Algo::CombBuilder::print(ostream& os){
  os << "\t\tCombBuilder contains " << combined.size() << " block(s):" << endl; 
  for(auto dec : combined) 
    dec->print(os);
}

//////////////////////////////////////////////////////


Algo::METBuilder::METBuilder() {
  decay    = Decay::MET;
  tf_met   = nullptr;
  verbose  = 0;
  saturate = 0;
}

Algo::METBuilder::METBuilder(const int& verb) {
  decay    = Decay::MET;
  tf_met   = nullptr;
  verbose  = verb;
  saturate = 0;
}


Algo::METBuilder::~METBuilder() {
  if(VERBOSE) cout << "Destroy METBuilder" << endl;
  if( tf_met!=nullptr) delete tf_met;
}


void Algo::METBuilder::init(const LV& lv) {
  p4_invisible = lv;
  tf_met = new TransferFunction("tf_met", TF_MET );
}

void Algo::METBuilder::fix_vars(){
  ++saturate;
}

double Algo::METBuilder::eval(const double* xx,  LV& lv) {

  double val{1.};
  
  if( saturate==0 )
    val *= tf_met->eval( p4_invisible.Px()-lv.Px(), p4_invisible.Py()-lv.Py() );
  else
    val *= tf_met->eval( 0., 0. );

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
  verbose = 0;
}

Algo::TopHadBuilder::TopHadBuilder (const int& verb) {
  decay   = Decay::TopHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
  tf_b    = nullptr;
  verbose = verb;
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
    tf->init( TF_B_param[Algo::eta_to_bin(lv)] );
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

  if(verbose) cout << "\tTopHadBuilder::eval xx[" << index_q << "]=" << E1 << endl;

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

  /*
  TLorentzVector w1 ( p4_q.Vect().Unit()   *E1, E1);
  TLorentzVector w2 ( p4_qbar.Vect().Unit()*E2, E2);
  TLorentzVector blv( p4_b.Vect().Unit()   *(TMath::Sqrt(E3*E3 - MB*MB)), E3);  
  cout << (w1+w2).M() << ", " << (w1+w2+blv).M() << endl;
  */

  return val;
}


//////////////////////////////////////////////////////

Algo::WHadBuilder::WHadBuilder () {
  decay   = Decay::WHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
  verbose = 0;
};

Algo::WHadBuilder::WHadBuilder (const int& verb) {
  decay   = Decay::WHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
  verbose = verb;
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

  if(verbose) cout << "\tWHadBuilder::eval xx[" << index_q << "]=" << E1 << endl;
  

  double a12 = (p4_q.Vect()).Angle(p4_qbar.Vect());
  
  E2 = MW*MW/E1/(4*TMath::Sin(a12/2.)*TMath::Sin(a12/2.));

  invisible += LV(0.,0.,0.,0.);

  val += tf_q->eval(p4_q.E(), E1) * tf_qbar->eval(p4_qbar.E(), E2) ;

  /*
  TLorentzVector w1 ( p4_q.Vect().Unit()   *E1, E1);
  TLorentzVector w2 ( p4_qbar.Vect().Unit()*E2, E2);
  */

  return val;

}

//////////////////////////////////////////////////////

Algo::HiggsBuilder::HiggsBuilder () {
  decay   = Decay::Higgs;
  errFlag = 0;
  tf_b    = nullptr;
  tf_bbar = nullptr;
  verbose = 0;
};

Algo::HiggsBuilder::HiggsBuilder (const int& verb) {
  decay   = Decay::Higgs;
  errFlag = 0;
  tf_b    = nullptr;
  tf_bbar = nullptr;
  verbose = verb;
};

Algo::HiggsBuilder::~HiggsBuilder() {
   if(VERBOSE) cout << "Destroy HiggsBuilder" << endl;
  if( tf_b   != nullptr ) delete tf_b;
  if( tf_bbar!= nullptr ) delete tf_bbar;
};

void Algo::HiggsBuilder::print(ostream& os){
  os << "\t\t\tDecay: " << Algo::translateDecay(decay) << endl; 
  os << "\t\t\tb    [" << index_b    << "]: p4 = (" << p4_b.Pt() << ", " << p4_b.Eta()  << ", " << p4_b.Phi() << ", " << p4_b.M() << ")" << endl; 
  os << "\t\t\tbbar [" << index_bbar << "]: p4 = (" << p4_bbar.Pt() << ", " << p4_bbar.Eta()  << ", " << p4_bbar.Phi() << ", " << p4_bbar.M() << ")" <<endl; 
  if(tf_b!=nullptr)    os << "\t\t\tTF b   : " << tf_b->getFormula() << endl;
  if(tf_bbar!=nullptr) os << "\t\t\tTF bbar: " << tf_bbar->getFormula() << endl;
}



void Algo::HiggsBuilder::init( const FinalState& fs, const LV& lv, const size_t& sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::Higgs_b:
    p4_b    = lv;
    index_b = sz;
    tf = new TransferFunction("tf_b", TF_B);
    tf->init( TF_B_param[Algo::eta_to_bin(lv)] );
    tf_b    = tf;
    break;
  case FinalState::Higgs_bbar:
    p4_bbar    = lv;
    index_bbar = sz;
    tf = new TransferFunction("tf_bbar", TF_B);
    tf->init( TF_B_param[Algo::eta_to_bin(lv)] );
    tf_bbar    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::HiggsBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val {0.};

  double E1 = xx[ index_b ];
  double E2;

  if(verbose) cout << "\tHiggsBuilder::eval xx[" << index_b << "]=" << E1 << endl;


  int nSol = 0;

  if(E1<MB){
    errFlag = 1;
    return val;
  }

  double a12 = (p4_b.Vect()).Angle(p4_bbar.Vect());
  double a = E1;
  double b = TMath::Sqrt(E1*E1 - MB*MB)*TMath::Cos(a12);

  if( (DMH2*DMH2 - (a*a - b*b)*MB*MB) < 0){
    errFlag = 1;
    return val;
  }
  double E2_1 = (a*DMH2 + b*TMath::Sqrt(DMH2*DMH2 - (a*a - b*b)*MB*MB))/(a*a - b*b);
  double E2_2 = (a*DMH2 - b*TMath::Sqrt(DMH2*DMH2 - (a*a - b*b)*MB*MB))/(a*a - b*b);
  double E2tmp1 = -999.;
  double E2tmp2 = -999.;

  if( b>0 ){
    if(E2_1>DMH2/a){
      E2tmp1 = E2_1;
      nSol++;
    }
    if(E2_2>DMH2/a){
      E2tmp2 = E2_2;
      nSol++;
    }
  }
  else{
    if(E2_1<DMH2/a){
      E2tmp1 = E2_1;
      nSol++;
    }
    if(E2_2<DMH2/a){
      E2tmp2 = E2_2;
      nSol++;
    }
  }
  if( E2tmp1>0 && E2tmp2>0 )
    E2 = TMath::Max( E2tmp1,E2tmp2 );
  else if( E2tmp1>0 && E2tmp2<0)
    E2 = E2tmp1;
  else if( E2tmp1<0 && E2tmp2>0)
    E2 = E2tmp2;
  else{
    errFlag = 1;
    return val;
  }
  if(E2<MB){
    errFlag = 1;
    return val;
  }

  /*
    TLorentzVector b1 ( p4_b.Vect().Unit()   *(TMath::Sqrt(E1*E1 - MB*MB)), E1);
    TLorentzVector b2 ( p4_bbar.Vect().Unit()*(TMath::Sqrt(E2*E2 - MB*MB)), E2);
    cout << (b1+b2).M() << endl;
  */

  invisible += LV(0.,0.,0.,0.);

  val += tf_b->eval(p4_b.E(), E1) * tf_bbar->eval(p4_bbar.E(), E2) ;

  return val;

}


//////////////////////////////////////////////////////

Algo::TopLepBuilder::TopLepBuilder () {
  decay   = Decay::TopLep;
  errFlag = 0;
  tf_b    = nullptr;
  verbose = 0;
}


Algo::TopLepBuilder::TopLepBuilder (const int& verb) {
  decay   = Decay::TopLep;
  errFlag = 0;
  tf_b    = nullptr;
  verbose = verb;
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
    tf->init( TF_B_param[Algo::eta_to_bin(lv)] );
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
  
  if(verbose) cout << "\tTopLepBuilder::eval xx[" << index_l << "]=" << nuPhi << ", xx[" << index_l +1 << "]=" << nuTheta << endl;

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
  
  /*
  TLorentzVector wLep( p4_l.Vect().Unit()*Elep, Elep);
  TLorentzVector blv ( p4_b.Vect().Unit() *(TMath::Sqrt(Eb*Eb - MB*MB)), Eb);  
  cout << (wNu+wLep).M() << ", " << (wNu+wLep+blv).M() << endl;
  */

  return val;
  
  
}


//////////////////////////////////////////////////////



Algo::RadiationBuilder::RadiationBuilder () {
  decay   = Decay::Radiation;
  tf_g    = nullptr;
  verbose = 0;
};

Algo::RadiationBuilder::RadiationBuilder (const int& verb) {
  decay   = Decay::Radiation;
  tf_g    = nullptr;
  verbose = verb;
};

Algo::RadiationBuilder::~RadiationBuilder() {
  if(VERBOSE) cout << "Destroy RadiationBuilder" << endl;
  if( tf_g   != nullptr ) delete tf_g;
};

void Algo::RadiationBuilder::print(ostream& os){
  os << "\t\t\tDecay: " << Algo::translateDecay(decay) << endl; 
  os << "\t\t\tq    [" << index_g    << "]: p4 = (" << p4_g.Pt() << ", " << p4_g.Eta()  << ", " << p4_g.Phi() << ", " << p4_g.M() << ")" << endl; 
  if(tf_g!=nullptr)    os << "\t\t\tTF g   : " << tf_g->getFormula() << endl;
}


void Algo::RadiationBuilder::init( const FinalState& fs, const LV& lv, const size_t& sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::Radiation_q:
    p4_g    = lv;
    index_g = sz;
    tf = new TransferFunction("tf_g", TF_Q);
    tf->init( TF_Q_param[Algo::eta_to_bin(lv)] );
    tf_g    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::RadiationBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val {0.};

  double E = xx[ index_g ];

  if(verbose) cout << "\tRadiationBuilder::eval xx[" << index_g << "]=" << E << endl;


  val += tf_g->eval(p4_g.E(), E);

  return val;
  
}
