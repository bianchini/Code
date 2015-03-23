#include "interface/Parameters.h"

size_t MEM::eta_to_bin( const double& eta ){
  if( fabs(eta)<1.0 ) return 0;
  if( fabs(eta)>1.0 ) return 1;
  return -99;
}
bool MEM::isQuark(const MEM::TFType::TFType& t) {
  return (t==TFType::bReco || t==TFType::qReco || t==TFType::bLost || t==TFType::qLost);
}
bool MEM::isNeutrino(const MEM::TFType::TFType& t) {
    return (t==TFType::MET);
  }
bool MEM::isLepton(const MEM::TFType::TFType& t)  {
  return (t==TFType::elReco || t==TFType::muReco);
}

double MEM::Chi2(const double& x, const double& m, const double& s){
  return s>0. ? (x-m)*(x-m)/s/s : 99.;
}

double MEM::Chi2Corr(const double& x, const double& y, const double& sx, const double& sy, const double& rho){
  return  1./(1-rho*rho)*( Chi2(x,0.,sx) +  Chi2(y,0.,sy) - 2*rho*x*y/sx/sy ) ;
}

/////////////////////////////////
//   y    := observables
//   x    := gen level quantities
//   type := decides the TF
/////////////////////////////////
double MEM::transfer_function(double* y, double* x, const TFType::TFType& type, int& out_of_range, const double& cutoff, const int& debug){

  // return value
  double w{1.};

  // temporary values;
  double E, H;
  double m1,s1, m2, s2, f, rho, c1, c2;

  // parameters
  const double* par;

  switch( type ){
    
  case TFType::bReco:
    // x[0] = parton energy ; 
    // x[1] = parton eta;
    // y[0] = jet energy;
    E  = x[0];
    H  = x[1];
    par = TF_B_param[ eta_to_bin(H) ];    

    f  = par[10];
    m1 = par[0] + par[1]*E;
    m2 = par[5] + par[6]*E;
    s1 = E*TMath::Sqrt(par[2]*par[2] + par[3]*par[3]/E + par[4]*par[4]/E/E);
    s2 = E*TMath::Sqrt(par[7]*par[7] + par[8]*par[8]/E + par[9]*par[9]/E/E);
    c1 = Chi2(y[0], m1, s1);
    c2 = Chi2(y[0], m2, s2);
    if( c1>cutoff && c2>cutoff ) ++out_of_range;
    w *= (1./sqrt(2*PI) * (f/s1*TMath::Exp(-0.5*c1) + (1-f)/s2*TMath::Exp(-0.5*c2) ));
#ifdef DEBUG_MODE
    if( debug&DebugVerbosity::integration) 
      cout << "\t\ttransfer_function: Evaluate W(" << y[0] << " | E=" << E << ", y=" << H << ", TFType::bReco) = " << w << endl;
#endif
    break;
    
  case TFType::qReco:
    // x[0] = parton energy ; 
    // x[1] = parton eta;
    // y[0] = jet energy;
    E   = x[0];
    H   = x[1];
    par = TF_Q_param[ eta_to_bin(H) ];
    m1  = par[0] + par[1]*E;
    s1  = E*TMath::Sqrt(par[2]*par[2] + par[3]*par[3]/E + par[4]*par[4]/E/E);
    c1  = Chi2(y[0], m1, s1);
    if( c1>cutoff ) ++out_of_range;

    w *= (1./sqrt(2*PI)/s1*TMath::Exp(-0.5*c1));
#ifdef DEBUG_MODE
    if( debug&DebugVerbosity::integration) 
      cout << "\t\ttransfer_function: Evaluate W(" << y[0] << " | E=" << E << ", y=" << H << ", TFType::qReco) = " << w << endl;
#endif
    break;
    
  case TFType::MET:
    // x[0] = sum nu_x ; x[1] = sum nu_y
    // y[0] = MET_x    ; y[1] = MET_y

    par = TF_MET_param;
    s1  = par[0];
    s2  = par[1];
    rho = par[2];
    c1  = Chi2Corr(y[0]-x[0], y[1]-x[1], s1, s2, rho);

    if( c1/2>cutoff ) ++out_of_range;

    w *= 1./(2*PI)/s1/s2/sqrt(1.-rho*rho)*TMath::Exp( -0.5*c1 );
#ifdef DEBUG_MODE
    if( debug&DebugVerbosity::integration) 
      cout << "\t\ttransfer_function: Evaluate W(" << y[0] << "-" << x[0] << " , " << y[1] << "-" << x[1] << ", TFType::MET) = " << w << endl;
#endif
    break;

  case TFType::Recoil:
    // Sudakov factor
    // x[0] = pT
    // y[0] = rhoT if extra_jets==0, else  par[2]+1GeV 

    par = TF_RECOIL_param;    
    m1 = par[0];
    s1 = par[1];
    if( y[0] < par[2] )
      w *= TMath::Gaus( log(x[0]), m1, s1, 1 );
    else 
      w *= 1.;
#ifdef DEBUG_MODE
    if( debug&DebugVerbosity::integration) 
      cout << "\t\ttransfer_function: Evaluate W( log(" << x[0] << "); TFType::Recoil) = " << w << endl;
#endif
    break;

  case TFType::bLost:
  case TFType::qLost:
    // x[0]     = parton energy ;
    // x[1]     = parton eta;
    // y[0]     = jet energy
    // par: [0]-> eta acceptance, [1]-> pT cut, [2]-> E max, [3]->acceptance (cos*phi)

    if( TMath::Abs(x[1])>TF_ACC_param[0] ){
      w *= 1.;      
#ifdef DEBUG_MODE
      if( debug&DebugVerbosity::integration) 
	cout << "\t\ttransfer_function: Evaluate W(" << x[0] << ", " << x[1] << ", TFType::qLost) = " << w << endl;
#endif
    }
    else{
      // x[0]     = parton energy ; 
      // x[1]     = parton eta
      // y[0]     = 0.
      E   = x[0];
      H   = x[1];
      par = TF_Q_param[ eta_to_bin(H) ];
      double mean_pt  = (par[0] + par[1]*E)/TMath::CosH(H);
      double sigma_pt = E*TMath::Sqrt(par[2]*par[2] + par[3]*par[3]/E + par[4]*par[4]/E/E)/TMath::CosH(H);
      c1 = Chi2( TF_ACC_param[1], mean_pt, sigma_pt);
      if( c1>cutoff ) ++out_of_range;

      w *= 0.5*(TMath::Erf( sqrt(c1) ) + 1 ) ;    
#ifdef DEBUG_MODE
      if( debug&DebugVerbosity::integration) 
	cout << "\t\ttransfer_function: Evaluate W(" <<  TF_ACC_param[1] << " | " << E << ", " << H << ", TFType::qLost) = " << w << endl;
#endif
    }
    break;
  default:
    break;
  }

  return w;
}

/////////////////////////////////                                                                                                                      
//   y     := observables                                                                                                                               
//   type  := decides the TF  
//   alpha := CL (e.g. 0.95, 0.98, ...)
///////////////////////////////// 
pair<double, double> MEM::get_support(double* y, const TFType::TFType& type, const double& alpha, const int& debug){

  if( type==TFType::TFType::MET ){

    double alpha_n = TMath::Abs(alpha);
    int sign       = alpha>0?1:0;

    // the MET px and py
    double Px        = y[0];
    double Py        = y[1];

    // return values
    double xLowPhi   = -TMath::Pi();
    double xHighPhi  = +TMath::Pi();

    double phiStep   = 0.04;

    // phi of MET vector
    double Phi       = Py>0?TMath::ACos(Px/sqrt(Px*Px+Py*Py)):2*TMath::Pi()-TMath::ACos(Px/sqrt(Px*Px+Py*Py));
    if( debug&DebugVerbosity::init_more) cout << "MET phi at " << Phi << endl;

    // elements of the MET cov matrix
    double Vx        = TF_MET_param[0]*TF_MET_param[0];
    double Vy        = TF_MET_param[1]*TF_MET_param[1];
    double rho       = TF_MET_param[2];

    // chi2 cut to find the CL
    double chi2Cut   = TMath::ChisquareQuantile(alpha_n ,2);    

    // MET TF at zero
    double tfAtZero = Chi2Corr(Px,Py,sqrt(Vx),sqrt(Vy),rho);

    // nothing to do...
    if( tfAtZero <= chi2Cut ){
      if( debug&DebugVerbosity::init_more)
	cout << "(0,0) is inside the 2-sigma CL => integrate over -TMath::Pi()/+TMath::Pi()" << endl;
    }

    // search for boundaries
    else{

      if( debug&DebugVerbosity::init_more)
	cout << "(0,0) is outside the 2-sigma CL => find phi-window with interpolation..." << endl;	
      
      for( int dir = 0; dir < 2 ; ++dir){	

	if(dir!=sign) continue;
	if( debug&DebugVerbosity::init_more) cout << "Doing scan along " << (dir?"+":"-") << " direction" << endl;

	bool stopPhiScan = false;
	for(std::size_t step = 0; step <= (std::size_t)(TMath::Pi()/phiStep) && !stopPhiScan; ++step){
	  double phi = Phi + double(2.*dir-1)*phiStep*step;
	  //if(phi<0.) phi += 2*TMath::Pi();
	  //else if(phi>2*TMath::Pi())  phi -= 2*TMath::Pi();
	  if( debug&DebugVerbosity::init_more) cout << "\tScan phi=" << phi << endl;
	  double sin    = TMath::Sin(phi);
	  double cos    = TMath::Cos(phi);
	  bool crossing = false;
	  //bool exceeded = false;
	  //bool alreadyInTheBox = PxMax*PxMin<=0. && PyMax*PyMin<=0.;

	  double p_step = 2.;
	  for(std::size_t stepP = 0; stepP<200 && !crossing && /*!exceeded && !crossing &&*/ !stopPhiScan; ++stepP){
	    double Px_P = stepP*p_step*cos;
	    double Py_P = stepP*p_step*sin;	    
	    //if( alreadyInTheBox && (Px_P>PxMax || Px_P<PxMin || Py_P>PyMax || Py_P<PyMin)){
	    //exceeded = true;
	    //if( debug&DebugVerbosity::init_more) cout << "\tWas in box, and got out at P=" << stepP*p_step << endl;
	    //continue;
	    //}
	    //else if( !alreadyInTheBox && (Px_P>PxMax || Px_P<PxMin || Py_P>PyMax || Py_P<PyMin)){
	      //if( debug&DebugVerbosity::init_more) cout << "\tWas not in the box, and I am still out at P=" << stepP*5 << endl;
	      //continue;
	    //}
	    
	    if( Chi2Corr(Px_P-Px, Py_P-Py, sqrt(Vx), sqrt(Vy), rho) < chi2Cut ) {
	      crossing = true;	 
	      if( debug&DebugVerbosity::init_more) 
		cout << "\tWas not in the box, and found crossing at (" << Px_P << "," << Py_P << ")" << endl;
	    }
	  } // end loop over |P|
	  
	  if(!crossing){
	    if( debug&DebugVerbosity::init_more) cout << "\tNo crossing at " << phi << " => stop phi scan" << endl;
	    if(dir==0) xLowPhi  = phi+0.5*phiStep;
	    if(dir==1) xHighPhi = phi-0.5*phiStep;
	    stopPhiScan = true;
	  }
	}
      }
      
      //xLowPhi  = -TMath::ACos(TMath::Cos( Phi - xLowPhi ));
      //xHighPhi = +TMath::ACos(TMath::Cos( Phi - xHighPhi));
      xLowPhi  -= Phi;
      xHighPhi -= Phi;
    }
    
    return make_pair(xLowPhi,xHighPhi);
  } 

  // the reconstructed values
  double e_rec   = y[0];
  double eta_rec = y[1];

  // start with reconstructed value
  double e_L{e_rec};
  double e_H{e_rec};

  // granularity
  double step_size{2.5};

  double tot{1.};
  int accept{0};
  double cutoff{99.};
  while( tot>(1-alpha)/2 && e_L>0. ){
    tot = 0.;
    for(size_t i = 0; i < 500.; ++i){
      double gen[2] = {e_L, eta_rec};
      double rec[1] = {e_rec+i*step_size};
      tot += transfer_function(rec,gen,type,accept,cutoff,debug)*step_size;
      if(  tot>(1-alpha)/2 ) break;
    }
    e_L -= step_size;
  }
  if(e_L<0.) e_L=0.;

  tot = 1.;
  while( tot>(1-alpha)/2 ){
    tot = 0.;
    for(size_t i = 0; i < 500.; ++i){
      double gen[2] = {e_H, eta_rec};
      double rec[1] = {e_rec-i*step_size};
      if(rec[0]<0.) continue;
      tot += transfer_function(rec,gen,type,accept,cutoff,debug)*step_size;
      if(  tot>(1-alpha)/2 ) break;
    }
    e_H += step_size;
  }
#ifdef DEBUG_MODE
  if( debug&DebugVerbosity::init_more) 
    cout << "MEM::get_support: E(reco) = " << e_rec << " ==> range at " << alpha 
	 << " CL is [" << e_L << ", " << e_H << "] (stepping every " << step_size << " GeV)" << endl;
#endif
  return make_pair(e_L, e_H);
}


MEM::PS::PS(size_t d){
  dim = d;
}

MEM::PS::~PS(){}

MEM::PSMap::const_iterator MEM::PS::begin() const {
  return val.begin();
}

MEM::PSMap::const_iterator MEM::PS::end() const {
  return val.end();
}

LV MEM::PS::lv(const MEM::PSPart::PSPart& p) const { 
  return val.find(p)!=val.end() ? (val.find(p)->second).lv : LV() ; 
}

int MEM::PS::charge(const MEM::PSPart::PSPart& p) const {
  return val.find(p)!=val.end() ? (val.find(p)->second).charge : 0 ;
}

MEM::TFType::TFType MEM::PS::type(const MEM::PSPart::PSPart& p) const {
  return val.find(p)!=val.end() ? (val.find(p)->second).type : TFType::Unknown ;
}

void MEM::PS::set(const MEM::PSPart::PSPart& a, const MEM::GenPart& b){
  val[a] = b;
}

void MEM::PS::print(ostream& os) const{
  cout << "\tContent of this PS: dim(PS)=" << dim << "..." << endl;
  for( auto p = val.begin() ; p != val.end() ; ++p ){
    cout << "\t\tPS[" << static_cast<size_t>(p->first) << "] : type("
	 << static_cast<size_t>(p->second.type) << "), (pT,h,phi,M)=("
	 << p->second.lv.Pt() << ", "
	 << p->second.lv.Eta() << ", "
	 << p->second.lv.Phi() << ", "
	 << p->second.lv.M() << "), (px,py,pz,E)=("
	 << p->second.lv.Px() << ", "
	 << p->second.lv.Py() << ", "
	 << p->second.lv.Pz() << ", "
	 << p->second.lv.E() << ")"
	 << endl; 
  }
}

MEM::Object::Object(const LV& lv, const MEM::ObjectType::ObjectType& ty){ 
  p  = lv; 
  t  = ty; 
} 

MEM::Object::Object(){ 
  p = LV(1e-06,0.,0.,1e-06); 
  t = ObjectType::Unknown; 
} 

MEM::Object::~Object(){}    

LV MEM::Object::p4() const { return p; }

MEM::ObjectType::ObjectType MEM::Object::type() const { return t; }

double MEM::Object::getObs(const MEM::Observable::Observable& name) const { 
  return (obs.find(name)!=obs.end() ? obs.find(name)->second : 0.);
}

bool MEM::Object::isSet(const MEM::Observable::Observable& name) const { 
  return obs.find(name)!=obs.end();
}

void MEM::Object::addObs(const MEM::Observable::Observable& name, const double& val){ 
  obs.insert( make_pair(name, val) ); 
}

void MEM::Object::print(ostream& os) const {
  os << "\tType: " << static_cast<int>(t) << ", p=(Pt, Eta, Phi, M)=("
     << p.Pt() << ", " << p.Eta() << ", " << p.Phi() << ", " << p.M()
     << ")" << endl;
}

MEM::MEMConfig::MEMConfig(int nmc, 
			  double ab, double re, 
			  int ic, int pi, 
			  double s, double e, 
			  std::string pdf, 
			  double jCL, double bCL, double mCL,
			  int tfsupp, double tfoff,
			  int hpf){
  n_max_calls  = nmc;
  abs          = ab;
  rel          = re;
  int_code     = ic;
  perm_int     = pi;
  sqrts        = s;
  emax         = e;
  pdfset       = pdf;
  is_default   = true;
  j_range_CL   = jCL;
  b_range_CL   = bCL;
  m_range_CL   = mCL;
  tf_suppress  = tfsupp;
  tf_offscale  = tfoff;
  highpt_first = hpf;
  for( int i = 0; i < 4 ; ++i){
    for( int j = 0; j < 2 ; ++j){
      for( int k = 0; k < 2 ; ++k){
	calls[i][j][k] = 2000;
      }
    }
  }
  perm_pruning = {};
}

void MEM::MEMConfig::defaultCfg(float nCallsMultiplier){
  
  // FinalState::LH
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTH)]
    [ static_cast<std::size_t>(Assumption::Assumption::ZeroQuarkLost)] = 2000;     
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTBB)]
    [ static_cast<std::size_t>(Assumption::Assumption::ZeroQuarkLost)] = 2000;
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTH)]
    [ static_cast<std::size_t>(Assumption::Assumption::OneQuarkLost)] = 4000;     
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTBB)]
    [ static_cast<std::size_t>(Assumption::Assumption::OneQuarkLost)] = 4000;

  // FinalState::LL
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LL) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTH)]
    [ static_cast<std::size_t>(Assumption::Assumption::ZeroQuarkLost)] = 10000;     
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LL) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTBB)]
    [ static_cast<std::size_t>(Assumption::Assumption::ZeroQuarkLost)] = 10000;      
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LL) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTH)]
    [ static_cast<std::size_t>(Assumption::Assumption::OneQuarkLost)]  = 20000;     
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::LL) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTBB)]
    [ static_cast<std::size_t>(Assumption::Assumption::OneQuarkLost)]  = 20000;      

  // FinalState::HH
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::HH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTH)]
    [ static_cast<std::size_t>(Assumption::Assumption::ZeroQuarkLost)] = 1000;     
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::HH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTBB)]
    [ static_cast<std::size_t>(Assumption::Assumption::ZeroQuarkLost)] = 1000;
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::HH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTH)]
    [ static_cast<std::size_t>(Assumption::Assumption::OneQuarkLost)] = 4000;     
  calls
    [ static_cast<std::size_t>(FinalState::FinalState::HH) ]
    [ static_cast<std::size_t>(Hypothesis::Hypothesis::TTBB)]
    [ static_cast<std::size_t>(Assumption::Assumption::OneQuarkLost)] = 4000;
  
  if (nCallsMultiplier != 1.0) {
    for (int i=0; i<4; i++) {
      for (int j=0; j<2; j++) {
        for (int k=0; k<2; k++) {
          calls[i][j][k] = nCallsMultiplier * calls[i][j][k];
        }
      }
    }
  }
  
  int_code = 
    IntegrandType::IntegrandType::Constant
    |IntegrandType::IntegrandType::ScattAmpl
    |IntegrandType::IntegrandType::DecayAmpl
    |IntegrandType::IntegrandType::Jacobian
    |IntegrandType::IntegrandType::PDF
    |IntegrandType::IntegrandType::Transfer;
    //|IntegrandType::IntegrandType::Sudakov
    //|IntegrandType::IntegrandType::Recoil;
  
  perm_pruning = {Permutations::BTagged, Permutations::QUntagged,
		  Permutations::QQbarSymmetry, Permutations::BBbarSymmetry};
}


void MEM::MEMConfig::setNCalls(FinalState::FinalState f, Hypothesis::Hypothesis h, Assumption::Assumption a, int n) {
  calls
    [ static_cast<std::size_t>(f) ]
    [ static_cast<std::size_t>(h)]
    [ static_cast<std::size_t>(a)] = n;
}

int MEM::MEMConfig::getNCalls(FinalState::FinalState f, Hypothesis::Hypothesis h, Assumption::Assumption a) {
  return calls
    [ static_cast<std::size_t>(f) ]
    [ static_cast<std::size_t>(h)]
    [ static_cast<std::size_t>(a)];
}
