#include "interface/ToyGenerator.h"


Algo::ToyGenerator::ToyGenerator(const int verb) {
  verbose = verb;
  ran     = new TRandom3();
}

Algo::ToyGenerator::ToyGenerator(const int verb, const unsigned int seed) {
  verbose = verb;
  ran     = new TRandom3();
  ran->SetSeed(seed);
}

Algo::ToyGenerator::~ToyGenerator(){
  delete ran;
}

vector<Algo::Object> Algo::ToyGenerator::generate( const vector<Decay>& decays, const int& smear, const int& btag){
  
  vector<Algo::Object> output;

  for( auto decay : decays ){
    generate_hypo( output, decay, smear, btag);
  }

  return output;

}


void Algo::ToyGenerator::generate_hypo( vector<Algo::Object>& out, Decay type, const int& smear, const int& btag ){
  
  double pt_top  = 20.;
  double eta_top = ran->Uniform(-2.5,2.5);
  double phi_top = ran->Uniform(-TMath::Pi(),TMath::Pi());
  LV p4_top;
  p4_top.SetPtEtaPhiM( pt_top, eta_top, phi_top, MTOP );

  double pt_h  = 20.;
  double eta_h = ran->Uniform(-2.5,2.5);
  double phi_h = ran->Uniform(-TMath::Pi(),TMath::Pi());
  LV p4_h;
  p4_h.SetPtEtaPhiM( pt_h, eta_h, phi_h, MH );


  double phiS_w    = ran->Uniform(-TMath::Pi(),TMath::Pi());
  double cthetaS_w = ran->Uniform(-1.,1.);
  double sthetaS_w = TMath::Sqrt(1-cthetaS_w*cthetaS_w);
  double ES_w      = (MTOP*MTOP + MW*MW - MB*MB)/2./MTOP;
  double ES_b      = (MTOP*MTOP - MW*MW + MB*MB)/2./MTOP;

  double P_w  = TMath::Sqrt(ES_w*ES_w - MW*MW); 
  double Px_w = sthetaS_w*TMath::Cos(phiS_w)*P_w;
  double Py_w = sthetaS_w*TMath::Sin(phiS_w)*P_w;
  double Pz_w = cthetaS_w*P_w;
  LV p4_w;
  p4_w.SetPxPyPzE( Px_w, Py_w, Pz_w, ES_w);
  LV p4_b;
  p4_b.SetPxPyPzE( -Px_w, -Py_w, -Pz_w, ES_b);

  p4_w.Boost( p4_top.BoostVector() );
  p4_b.Boost( p4_top.BoostVector() );

  double phiS_q    = ran->Uniform(-TMath::Pi(),TMath::Pi());
  double cthetaS_q = ran->Uniform(-1.,1.);
  double sthetaS_q = TMath::Sqrt(1-cthetaS_q*cthetaS_q);
  double ES_q      = MW/2.;

  double P_q       = ES_q;
  double Px_q      = sthetaS_q*TMath::Cos(phiS_q)*P_q;
  double Py_q      = sthetaS_q*TMath::Sin(phiS_q)*P_q;
  double Pz_q      = cthetaS_q*P_q;

  LV p4_q, p4_qbar;
  p4_q.SetPxPyPzE   (  Px_q,  Py_q,  Pz_q, P_q );
  p4_qbar.SetPxPyPzE( -Px_q, -Py_q, -Pz_q, P_q );

  p4_q.Boost   (p4_w.BoostVector());
  p4_qbar.Boost(p4_w.BoostVector());

  double phiS_bh    = ran->Uniform(-TMath::Pi(),TMath::Pi());
  double cthetaS_bh = ran->Uniform(-1.,1.);
  double sthetaS_bh = TMath::Sqrt(1-cthetaS_w*cthetaS_w);
  double ES_bh      = MH/2.;

  double P_bh  = TMath::Sqrt(ES_bh*ES_bh - MB*MB);
  double Px_bh = sthetaS_bh*TMath::Cos(phiS_bh)*P_bh;
  double Py_bh = sthetaS_bh*TMath::Sin(phiS_bh)*P_bh;
  double Pz_bh = cthetaS_bh*P_bh;
  LV p4_bh;
  p4_bh.SetPxPyPzE   ( Px_bh, Py_bh, Pz_bh, ES_bh);
  LV p4_bbarh;
  p4_bbarh.SetPxPyPzE( -Px_bh, -Py_bh, -Pz_bh, ES_bh);

  p4_bh.Boost   ( p4_h.BoostVector() );
  p4_bbarh.Boost( p4_h.BoostVector() );

  if(verbose>2){
    cout << "Top E    = " << (p4_q + p4_qbar + p4_b).E() << " (" << p4_top.E() << ")" << endl;
    cout << "Top Px   = " << (p4_q + p4_qbar + p4_b).Px() << " (" << p4_top.Px() << ")" << endl;
    cout << "Top Py   = " << (p4_q + p4_qbar + p4_b).Py() << " (" << p4_top.Py() << ")" << endl;
    cout << "Top Pz   = " << (p4_q + p4_qbar + p4_b).Pz() << " (" << p4_top.Pz() << ")" << endl;
    cout << "Higgs E  = " << (p4_bh + p4_bbarh ).E() << " (" << p4_h.E() << ")" << endl;
    cout << "Higgs Px = " << (p4_bh + p4_bbarh ).Px() << " (" << p4_h.Px() << ")" << endl;
    cout << "Higgs Py = " << (p4_bh + p4_bbarh ).Py() << " (" << p4_h.Py() << ")" << endl;
    cout << "Higgs Pz = " << (p4_bh + p4_bbarh ).Pz() << " (" << p4_h.Pz() << ")" << endl;
  }

  // at most three outputs
  Algo::Object out1, out2, out3; 

  switch(type){

  case Decay::TopHad:
    if(smear){
      smear_by_TF(p4_q,   'q');
      smear_by_TF(p4_qbar,'q');
      smear_by_TF(p4_b,   'b');
    }
    out1.init(p4_q,    'q');
    out2.init(p4_qbar, 'q');
    out3.init(p4_b,    'b');
    if(btag){
      assign_rnd_btag( Algo::QuarkTypeUp,     out1 );
      assign_rnd_btag( Algo::QuarkTypeDown,   out2 );
      assign_rnd_btag( Algo::QuarkTypeBottom, out3 );
    }
    out.push_back( out1 );
    out.push_back( out2 );
    out.push_back( out3 );
    break;
  case Decay::TopLep:
    if(smear){
      smear_by_TF(p4_b,   'b');
      smear_by_TF(p4_qbar,'m');
    }
    out1.init( p4_b,    'b');
    out2.init( p4_q,    'l');
    out3.init( p4_qbar, 'm');   
    if(btag){
      assign_rnd_btag( Algo::QuarkTypeBottom, out1 ); 
    }
    out.push_back( out1 );
    out.push_back( out2 );
    out.push_back( out3 );
    break;
  case Decay::WHad:
    if(smear){
      smear_by_TF(p4_q,   'q');
      smear_by_TF(p4_qbar,'q');
    }
    out1.init( p4_q,    'q');
    out2.init( p4_qbar, 'q');
    if(btag){
      assign_rnd_btag( Algo::QuarkTypeUp,   out1 );
      assign_rnd_btag( Algo::QuarkTypeDown, out2 );
    }
    out.push_back( out1 );
    out.push_back( out2 );
    break;
  case Decay::Higgs:
    if(smear){
      smear_by_TF(p4_bh,   'b');
      smear_by_TF(p4_bbarh,'b');
    }
    out1.init( p4_bh,    'b');
    out2.init( p4_bbarh, 'b');
    if(btag){
      assign_rnd_btag( Algo::QuarkTypeBottom,   out1 ); 
      assign_rnd_btag( Algo::QuarkTypeBottom,   out2 ); 
    }
    out.push_back( out1 );
    out.push_back( out2 );
    break;
  case Decay::Radiation_q:
    if(smear) 
      smear_by_TF(p4_q, 'q');
    out1.init( p4_q,    'q');
    if(btag){
      assign_rnd_btag( Algo::QuarkTypeDown,   out1 ); 
    }
    out.push_back( out1 );
    break;
  case Decay::Radiation_b:
    if(smear)
      smear_by_TF(p4_b, 'b');
    out1.init( p4_b,    'b');
    if(btag){
      assign_rnd_btag( Algo::QuarkTypeBottom,   out1 ); 
    }
    out.push_back( out1 );
    break;

  default:
    break;
  }

  return;
}

void Algo::ToyGenerator::assign_rnd_btag( const Algo::QuarkType type, Object& obj){
  TF1 pdf("pdf_btag", Algo::pdf_btag , 0., 1., 3);
  pdf.SetParameter(0, static_cast<int>(type) );
  pdf.SetParameter(1, obj.p4.Pt()  ); 
  pdf.SetParameter(2, obj.p4.Eta() ); 
  obj.addObs( "BTAG", pdf.GetRandom() );
}


void Algo::ToyGenerator::smear_by_TF(TLorentzVector& lv, char type){
  
  string func;

  switch(type){
  case 'q':
    func = TF_Q;
    break;
  case 'b':
    func = TF_B;
    break;
  case 'm':
    func = TF_MET;
    break;
  default:
    cout << "Unknown smear type" << endl;
    return;
    break;
  }

  TransferFunction tf("tf", func, verbose); 

  if(type=='q' || type=='b'){

    if(type=='q') 
      tf.init( TF_Q_param[Algo::eta_to_bin(lv)] );
    else
      tf.init( TF_B_param[Algo::eta_to_bin(lv)] );
    
    string func_eval = tf.getFormula();
    std::ostringstream os;
    os << lv.E();
    while( func_eval.find("y")!=string::npos ){
      func_eval.replace(func_eval.find("y"), 1 , os.str());
    }
    TF1 f("f", func_eval.c_str(), lv.E() * SMEAR_EREL_MIN, lv.E() * SMEAR_EREL_MAX );
    double E_smear = TMath::Max( f.GetRandom(), 0.);

    if(verbose>0)
      cout << func_eval << ": E := " << lv.E() ;
    lv *= E_smear/lv.E() ;
    if(verbose>0)
      cout << " --> " << lv.E() << endl;
  }


  if(type=='m'){
    string func_eval = tf.getFormula(); 
    TF2 f("f", func_eval.c_str(), SMEAR_MET_PX_MIN, SMEAR_MET_PX_MAX, SMEAR_MET_PY_MIN, SMEAR_MET_PY_MAX);
    double Px{0.};
    double Py{0.};
    f.GetRandom2(Px,Py);
    if(verbose>0)   
      cout << func_eval << ": (Px,Py) := (" << lv.Px() << "," << lv.Py() ;
    lv += TLorentzVector(Px,Py,0,Px+Py);
    if(verbose>0)   
      cout << ") --> (" << lv.Px() << "," << lv.Py() << ")" << endl;
  }

  return;
}




