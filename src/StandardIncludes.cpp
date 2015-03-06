#include "interface/StandardIncludes.h"


Algo::TransferFunction::TransferFunction(const string& name, const string& form, const int verb){
  formula   = form;
  f         = new TFormula(name.c_str(), formula.c_str());
  verbose   = verb;
  threshold = 0.;
}

Algo::TransferFunction::~TransferFunction(){
  if(verbose>2) cout << "Destroy tf " << string(f->GetName()) << endl;
  delete f;
}

const string Algo::TransferFunction::getFormula() const {
  return string(f->GetExpFormula("p"));
}

double Algo::TransferFunction::eval(const double& rec, const double& gen) const {
  double val{1.0};
  val *= get_pdfs();
  val *= f->Eval(rec,gen);
  return val;
}

double Algo::TransferFunction::get_threshold() const{
  return threshold;
}

void  Algo::TransferFunction::set_threshold(const double& val){
  threshold = val;
}

void Algo::TransferFunction::add_pdf_obs( const string& name, const Object& obj, const int& type){

  if( obj.obs.find(name)==obj.obs.end() ){
    cout << "Algo::TransferFunction::add_pdf_obs(): Observable " << name << " is not available" << endl;
    return;
  }

  if( name == "BTAG" && pdfs.find(name) == pdfs.end() ){
      TF1 pdf(("pdf_"+name).c_str(), Algo::pdf_btag , 0., 1. , 3);
      pdf.SetParameter(0, type);
      pdf.SetParameter(1, obj.p4.Pt()  );
      pdf.SetParameter(2, obj.p4.Eta() );
      // RND=0 means that btagging is used as YES/NO variable, then do not use btag likelihood  
      if( (obj.obs).find(name+"_RND")!=(obj.obs).end() && (obj.obs).find(name+"_RND")->second<0.5 )
	pdfs[name] = 1.0;
      else	
	pdfs[name] = pdf.Eval( (obj.obs).find(name)->second );

      if( verbose>1 ){
	cout << "\t\tAlgo::TransferFunction::add_pdf_obs(): init pdf name " << name << " with function: " << endl;
	cout << "\t\tpdf(x=" << (obj.obs).find(name)->second << ") = " <<  pdfs[name] << endl;
      }    
  }

  /* implement others here */
  
}

double Algo::TransferFunction::get_pdfs() const {
  double val{1.0};
  for( auto pdf : pdfs) val *= pdf.second;
  return val;    
}

void Algo::TransferFunction::init(const double* param){
  assert(f!=nullptr);
  f->SetParameters( param );
}


//////////////////////////////////////////////////////

Algo::CombBuilder::CombBuilder(vector<DecayBuilder*> comb, int verb) {
  combined = comb;
  verbose  = verb;
}

Algo::CombBuilder::~CombBuilder() {
  if(verbose>2) cout << "Destroy CombBuilder" << endl;
  for(auto comb : combined) 
    delete comb;
}

void Algo::CombBuilder::add(DecayBuilder* comb) {
  combined.push_back(comb);
}

double Algo::CombBuilder::eval(const double* xx) {

  if(verbose>1) cout << "\tCombBuilder::eval(): BEGIN" << endl;

  double val{1.};

  // this vector contains the total 4-momentum from invisible particles
  LV invisible(0.,0.,0.,0.);

  int count {0};  
  for(auto comb : combined){
    double val_tmp = comb->eval(xx, invisible);
    val *= val_tmp;
    if(verbose>1){
      cout << "\tEval block " << count << " += " << (val_tmp>0. ? -TMath::Log( val_tmp ) : numeric_limits<double>::max() ) ;
      cout << "; Invisible = " << "(E=" << invisible.E() << ", cos(theta)=" 
	   << TMath::Cos(invisible.Theta()) << ", phi=" << invisible.Phi() << ")" << endl;
    }
    ++count;
  }

  if(verbose>1) cout << "\tCombBuilder::eval(): END" << endl;

  return val;
}

double Algo::CombBuilder::eval(const double* xx,  LV& lv) {
  return eval( xx ); 
}


void Algo::CombBuilder::print(ostream& os){
  os << "\t\tCombBuilder contains " << combined.size() << " block(s):" << endl; 
  for(auto comb : combined) 
    comb->print(os);
}

Algo::DecayBuilder* Algo::CombBuilder::at( const std::size_t pos ){
  return ( (pos<combined.size()) ? combined[pos] : nullptr);
}

std::size_t Algo::CombBuilder::size(){
  return combined.size();
}

Algo::Decay::Decay Algo::CombBuilder::get_decay(){
  return Decay::Decay::UNKNOWN;
}

//////////////////////////////////////////////////////


Algo::METBuilder::METBuilder(int verb) {
  decay    = Decay::Decay::MET;
  tf_met   = nullptr;
  verbose  = verb;
  saturate = 0;
}


Algo::METBuilder::~METBuilder() {
  if(verbose>2) cout << "Destroy METBuilder" << endl;
  if( tf_met!=nullptr) delete tf_met;
}


void Algo::METBuilder::init(const Object& obj) {
  p4_invisible = obj.p4;
  tf_met = new TransferFunction("tf_met", TF_MET , verbose);
}

void Algo::METBuilder::fix_vars(){
  ++saturate;
}

double Algo::METBuilder::eval(const double* xx,  LV& lv) {

  if(verbose>1) cout << "\tMETBuilder::eval()" << endl;

  double val{1.};
  
  if( saturate==0 )
    val *= tf_met->eval( p4_invisible.Px()-lv.Px(), p4_invisible.Py()-lv.Py() );
  else
    val *= tf_met->eval( 0., 0. );

  return val;
}


void Algo::METBuilder::print(ostream& os){
  os << "\t\t\tDecay: " <<  Algo::Decay::translateDecay(decay) <<  endl; 
  os << "\t\t\tMET:   p4 = (" << p4_invisible.Pt() << ", " << p4_invisible.Phi() << ")" << endl;
  if(tf_met!=nullptr)  os << "\t\t\tTF MET   : " << tf_met->getFormula() << endl;
}

Algo::Decay::Decay Algo::METBuilder::get_decay(){
  return decay;
}

//////////////////////////////////////////////////////

Algo::TopHadBuilder::TopHadBuilder (int verb) {
  decay     = Decay::Decay::TopHad;
  errFlag   = 0;
  tf_q      = nullptr;
  tf_qbar   = nullptr;
  tf_b      = nullptr;
  verbose   = verb;
  qbar_lost = 0;
}

Algo::TopHadBuilder::~TopHadBuilder() {
  if(verbose>2) cout << "Destroy TopHadBuilder" << endl;
  if( tf_q   != nullptr ) delete tf_q;
  if( tf_qbar!= nullptr ) delete tf_qbar;
  if( tf_b   != nullptr ) delete tf_b;
}

void Algo::TopHadBuilder::print(ostream& os){
  os << "\t\t\tDecay: " <<  Algo::Decay::translateDecay(decay) << endl; 
  os << "\t\t\tq    [" << index_q    << "]: p4 = (" << p4_q.Pt() << ", " << p4_q.Eta()  << ", " << p4_q.Phi() << ", " << p4_q.M() << ")" << endl; 
  if(qbar_lost==0)
    os << "\t\t\tqbar [" << index_qbar << "]: p4 = (" << p4_qbar.Pt() << ", " << p4_qbar.Eta()  << ", " << p4_qbar.Phi() << ", " << p4_qbar.M() << ")" <<endl; 
  os << "\t\t\tb    [" << index_b    << "]: p4 = (" << p4_b.Pt() << ", " << p4_b.Eta()  << ", " << p4_b.Phi() << ", " << p4_b.M() << ")" <<endl; 
  if(tf_q!=nullptr)    os << "\t\t\tTF q   : " << tf_q->getFormula() << endl;
  if(tf_qbar!=nullptr) os << "\t\t\tTF qbar: " << tf_qbar->getFormula() << endl;
  if(tf_b!=nullptr)    os << "\t\t\tTF b   : " << tf_b->getFormula() << endl;
}



void Algo::TopHadBuilder::init( const FinalState::FinalState& fs, const Object& obj, const std::size_t sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::FinalState::TopHad_q:
    p4_q    = obj.p4;
    index_q = sz;
    tf = new TransferFunction("tf_q", TF_Q, verbose);    
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeUp );      
    tf_q    = tf;
    break;
  case FinalState::FinalState::TopHad_qbar:
    p4_qbar    = obj.p4;
    index_qbar = sz;
    tf = new TransferFunction("tf_qbar", TF_Q , verbose);
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeDown );      
    tf_qbar    = tf;
    break;
  case FinalState::FinalState::TopHadLost_qbar:
    index_qbar = sz;
    tf = new TransferFunction("tfLost_qbar", TF_Q_CUM , verbose);
    tf->set_threshold(PTTHRESHOLD);
    tf_qbar    = tf;
    ++qbar_lost;
    break;
  case FinalState::FinalState::TopHad_b:
    p4_b    = obj.p4;
    index_b = sz;
    tf = new TransferFunction("tf_b", TF_B , verbose);
    tf->init( TF_B_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeBottom );      
    tf_b    = tf;
    break;
  default:
    break;
  }

}
 
double Algo::TopHadBuilder::eval( const double *xx, LV& invisible ) {

  // return value
  double val {0.};
  int nSol   {0};

  double E1 = xx[ index_q ];

  if(verbose>1) cout << "\tTopHadBuilder::eval() xx[" << index_q << "]=" << E1;

  if( TMath::IsNaN(E1) ){
    if(verbose>1) cout << "xx[" << index_q << "] is nan" << endl;
    return val;
  }

  // is the second quark reconstructed ?
  if( qbar_lost ){
    double qbarPhi   =  xx[ index_qbar ] ;
    double qbarTheta =  xx[ index_qbar + 1] ;
    if(verbose>1) cout << ", xx[" << index_qbar << "]=" << qbarPhi << ", xx[" << index_qbar+1 << "]=" << qbarTheta <<  endl;
    TVector3 e2(0.,0.,1.);
    e2.SetTheta( TMath::ACos( qbarTheta ) );
    e2.SetPhi ( qbarPhi );
    e2.SetMag ( 1.);
    p4_qbar = LV(e2, 1.0 );
    tf_qbar->init( TF_Q_param[Algo::eta_to_bin(p4_qbar)] );
  }
  else{
    if(verbose>1) cout << endl;
  }

  double E2, E3;
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
  
  if(qbar_lost){
    double acc = 1.0;
    if( p4_qbar.Eta()<2.5 )
      acc = tf_qbar->eval( tf_qbar->get_threshold()/TMath::Sin(p4_qbar.Theta()) , E2) ;
    val += tf_q->eval(p4_q.E(), E1) * acc * tf_b->eval(p4_b.E(), E3);
    if(verbose>2) cout << "\tAcceptance: " << acc << " for E=" <<  E2 << " (Pt=" <<  E2*TMath::Sin(p4_qbar.Theta()) << ")" << endl;
    return val;
  }

  val += tf_q->eval(p4_q.E(), E1) * tf_qbar->eval(p4_qbar.E(), E2) * tf_b->eval(p4_b.E(), E3);

  /*
  TLorentzVector w1 ( p4_q.Vect().Unit()   *E1, E1);
  TLorentzVector w2 ( p4_qbar.Vect().Unit()*E2, E2);
  TLorentzVector blv( p4_b.Vect().Unit()   *(TMath::Sqrt(E3*E3 - MB*MB)), E3);  
  cout << (w1+w2).M() << ", " << (w1+w2+blv).M() << endl;
  */

  return val;
}

Algo::Decay::Decay Algo::TopHadBuilder::get_decay(){
  return decay;
}

vector<std::size_t> Algo::TopHadBuilder::get_variables(){
  vector<std::size_t> out;
  out.push_back( index_q );
  if(qbar_lost){
    out.push_back( index_qbar );
    out.push_back( index_qbar + 1 );
  }
  return out;
}
//////////////////////////////////////////////////////

Algo::WHadBuilder::WHadBuilder (int verb) {
  decay   = Decay::Decay::WHad;
  errFlag = 0;
  tf_q    = nullptr;
  tf_qbar = nullptr;
  verbose = verb;
};

Algo::WHadBuilder::~WHadBuilder() {
  if(verbose>2) cout << "Destroy WHadBuilder" << endl;
  if( tf_q   != nullptr ) delete tf_q;
  if( tf_qbar!= nullptr ) delete tf_qbar;
};

void Algo::WHadBuilder::print(ostream& os){
  os << "\t\t\tDecay: " << Algo::Decay::translateDecay(decay) << endl; 
  os << "\t\t\tq    [" << index_q    << "]: p4 = (" << p4_q.Pt() << ", " << p4_q.Eta()  << ", " << p4_q.Phi() << ", " << p4_q.M() << ")" << endl; 
  os << "\t\t\tqbar [" << index_qbar << "]: p4 = (" << p4_qbar.Pt() << ", " << p4_qbar.Eta()  << ", " << p4_qbar.Phi() << ", " << p4_qbar.M() << ")" <<endl; 
  if(tf_q!=nullptr)    os << "\t\t\tTF q   : " << tf_q->getFormula() << endl;
  if(tf_qbar!=nullptr) os << "\t\t\tTF qbar: " << tf_qbar->getFormula() << endl;
}


void Algo::WHadBuilder::init( const FinalState::FinalState& fs, const Object& obj, const std::size_t sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::FinalState::WHad_q:
    p4_q    = obj.p4;
    index_q = sz;
    tf = new TransferFunction("tf_q", TF_Q , verbose);
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeUp );      
    tf_q    = tf;
    break;
  case FinalState::FinalState::WHad_qbar:
    p4_qbar    = obj.p4;
    index_qbar = sz;
    tf = new TransferFunction("tf_qbar", TF_Q , verbose);
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeDown );      
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

  if(verbose>1) cout << "\tWHadBuilder::eval() xx[" << index_q << "]=" << E1 << endl;
  
  if( TMath::IsNaN(E1) ){
    if(verbose>1) cout << "xx[" << index_q << "] is nan" << endl;
    return val;
  }

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

Algo::Decay::Decay Algo::WHadBuilder::get_decay(){
  return decay;
}

vector<std::size_t> Algo::WHadBuilder::get_variables(){
  vector<std::size_t> out;
  out.push_back( index_q );
  return out;
}

//////////////////////////////////////////////////////

Algo::HiggsBuilder::HiggsBuilder (int verb) {
  decay   = Decay::Decay::Higgs;
  errFlag = 0;
  tf_b    = nullptr;
  tf_bbar = nullptr;
  verbose = verb;
};

Algo::HiggsBuilder::~HiggsBuilder() {
  if(verbose>2) cout << "Destroy HiggsBuilder" << endl;
  if( tf_b   != nullptr ) delete tf_b;
  if( tf_bbar!= nullptr ) delete tf_bbar;
};

void Algo::HiggsBuilder::print(ostream& os){
  os << "\t\t\tDecay: " << Algo::Decay::translateDecay(decay) << endl; 
  os << "\t\t\tb    [" << index_b    << "]: p4 = (" << p4_b.Pt() << ", " << p4_b.Eta()  << ", " << p4_b.Phi() << ", " << p4_b.M() << ")" << endl; 
  os << "\t\t\tbbar [" << index_bbar << "]: p4 = (" << p4_bbar.Pt() << ", " << p4_bbar.Eta()  << ", " << p4_bbar.Phi() << ", " << p4_bbar.M() << ")" <<endl; 
  if(tf_b!=nullptr)    os << "\t\t\tTF b   : " << tf_b->getFormula() << endl;
  if(tf_bbar!=nullptr) os << "\t\t\tTF bbar: " << tf_bbar->getFormula() << endl;
}



void Algo::HiggsBuilder::init( const FinalState::FinalState& fs, const Object& obj, const std::size_t sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::FinalState::Higgs_b:
    p4_b    = obj.p4;
    index_b = sz;
    tf = new TransferFunction("tf_b", TF_B , verbose);
    tf->init( TF_B_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeBottom );      
    tf_b    = tf;
    break;
  case FinalState::FinalState::Higgs_bbar:
    p4_bbar    = obj.p4;
    index_bbar = sz;
    tf = new TransferFunction("tf_bbar", TF_B , verbose);
    tf->init( TF_B_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeBottom );      
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

  if(verbose>1) cout << "\tHiggsBuilder::eval() xx[" << index_b << "]=" << E1 << endl;

  if( TMath::IsNaN(E1) ){
    if(verbose>1) cout << "xx[" << index_b << "] is nan" << endl;
    return val;
  }
  
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

Algo::Decay::Decay Algo::HiggsBuilder::get_decay(){
  return decay;
}

vector<std::size_t> Algo::HiggsBuilder::get_variables(){
  vector<std::size_t> out;
  out.push_back( index_b );
  return out;
}

//////////////////////////////////////////////////////

Algo::TopLepBuilder::TopLepBuilder (int verb) {
  decay   = Decay::Decay::TopLep;
  errFlag = 0;
  tf_b    = nullptr;
  verbose = verb;
}

Algo::TopLepBuilder::~TopLepBuilder() {
  if(verbose>2) cout << "Destroy TopLepBuilder" << endl;
  if( tf_b   != nullptr ) delete tf_b;
}

void Algo::TopLepBuilder::print(ostream& os){
  os << "\t\t\tDecay: " <<  Algo::Decay::translateDecay(decay) << endl; 
  os << "\t\t\tlep  [" << index_l    << "]: p4 = (" << p4_l.Pt() << ", " << p4_l.Eta()  << ", " << p4_l.Phi() << ", " << p4_l.M() << ")" <<endl; 
  os << "\t\t\tb    [" << index_b    << "]: p4 = (" << p4_b.Pt() << ", " << p4_b.Eta()  << ", " << p4_b.Phi() << ", " << p4_b.M() << ")" <<endl; 
  if(tf_b!=nullptr)    os << "\t\t\tTF b   : " << tf_b->getFormula() << endl;
}



void Algo::TopLepBuilder::init( const FinalState::FinalState& fs, const Object& obj, const std::size_t sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::FinalState::TopLep_l:
    p4_l    = obj.p4;
    index_l = sz;
    break;
  case FinalState::FinalState::TopLep_b:
    p4_b    = obj.p4;
    index_b = sz;
    tf = new TransferFunction("tf_b", TF_B , verbose);
    tf->init( TF_B_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeBottom );      
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
  double nuPhi    = xx[ index_l    ] ;
  double nuTheta  = xx[ index_l + 1] ;
  
  if(verbose>1) cout << "\tTopLepBuilder::eval() xx[" << index_l << "]=" << nuPhi << ", xx[" << (index_l+1) << "]=" << nuTheta << endl;

  if( TMath::IsNaN(nuPhi) ){
    if(verbose>1) cout << "xx[" << index_l << "] is nan" << endl;
    return val;
  }
  if( TMath::IsNaN(nuTheta) ){
    if(verbose>1) cout << "xx[" << index_l +1 << "] is nan" << endl;
    return val;
  }

  double Enu, Eb;

  TVector3 e3(0.,0.,1.); // neutrino
  e3.SetTheta( TMath::ACos( nuTheta ) );
  e3.SetPhi ( nuPhi );
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

Algo::Decay::Decay Algo::TopLepBuilder::get_decay(){
  return decay;
}

vector<std::size_t> Algo::TopLepBuilder::get_variables(){
  vector<std::size_t> out;
  out.push_back( index_l );
  out.push_back( index_l + 1);
  return out;
}
//////////////////////////////////////////////////////

Algo::RadiationBuilder::RadiationBuilder (int verb, Decay::Decay type) { 
  decay   = type;
  verbose = verb;
  tf_g    = nullptr;   
}; 

Algo::RadiationBuilder::~RadiationBuilder() {
  if(verbose>2) cout << "Destroy RadiationBuilder" << endl;
  if( tf_g   != nullptr ) delete tf_g;
};

Algo::Decay::Decay Algo::RadiationBuilder::get_decay() {
  return decay;
}

vector<std::size_t> Algo::RadiationBuilder::get_variables(){
  vector<std::size_t> out;
  out.push_back( index_g );
  return out;
}

void Algo::RadiationBuilder::print(ostream& os){
  os << "\t\t\tDecay: " << Algo::Decay::translateDecay(decay) << endl; 
  os << "\t\t\trad" << "    [" << index_g    << "]: p4 = (" << p4_g.Pt() << ", " << p4_g.Eta()  << ", " << p4_g.Phi() << ", " << p4_g.M() << ")" << endl; 
  if(tf_g!=nullptr)    os << "\t\t\tTF " << tf_g->getFormula() << endl;
}


void Algo::RadiationBuilder::init( const FinalState::FinalState& fs, const Object& obj, const std::size_t sz ){

  TransferFunction* tf = nullptr;

  switch( fs ){
  case FinalState::FinalState::Radiation_u:
    p4_g    = obj.p4;
    index_g = sz;
    tf = new TransferFunction("tf_q", TF_Q , verbose);
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeUp );      
    tf_g    = tf;
    break;
  case FinalState::FinalState::Radiation_d:
    p4_g    = obj.p4;
    index_g = sz;
    tf = new TransferFunction("tf_q", TF_Q , verbose);
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeDown );
    tf_g    = tf;
    break;
  case FinalState::FinalState::Radiation_c:
    p4_g    = obj.p4;
    index_g = sz;
    tf = new TransferFunction("tf_q", TF_Q , verbose);
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeCharm );
    tf_g    = tf;
    break;
  case FinalState::FinalState::Radiation_b:
    p4_g    = obj.p4;
    index_g = sz;
    tf = new TransferFunction("tf_b", TF_B , verbose);
    tf->init( TF_B_param[Algo::eta_to_bin(obj.p4)] );
    for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeBottom );      
    tf_g    = tf;
    break;
  case FinalState::FinalState::Radiation_g: //do not assume the parton flavour
    p4_g    = obj.p4;
    index_g = sz;
    tf = new TransferFunction("tf_q", TF_Q , verbose);
    tf->init( TF_Q_param[Algo::eta_to_bin(obj.p4)] );    
    //for(auto iobs : obj.obs ) tf->add_pdf_obs( iobs.first, obj, Algo::QuarkTypeDown );
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

  if(verbose>1) cout << "\tRadiationBuilder::eval() xx[" << index_g << "]=" << E << endl;

  if( TMath::IsNaN(E) ){
    if(verbose>1) cout << "xx[" << index_g << "] is nan" << endl;
    return val;
  }

  val += tf_g->eval(p4_g.E(), E);

  return val;
  
}
