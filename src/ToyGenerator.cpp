#include "interface/ToyGenerator.h"


Algo::ToyGenerator::ToyGenerator(const int verb) {
  verbose = verb;
  ran     = new TRandom3();
}

Algo::ToyGenerator::~ToyGenerator(){
  delete ran;
}

vector<pair<char,TLorentzVector>> Algo::ToyGenerator::generate( const vector<Decay>& decays, const int& smear){

  vector<pair<char,TLorentzVector>> output;

  for( auto decay : decays ){
    generate_hypo( output, decay, smear);
  }

  return output;

}


void Algo::ToyGenerator::generate_hypo( vector<pair<char,TLorentzVector>>& out, Decay type, const int& smear ){
  
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

  if(verbose){
    cout << "Top E    = " << (p4_q + p4_qbar + p4_b).E() << " (" << p4_top.E() << ")" << endl;
    cout << "Top Px   = " << (p4_q + p4_qbar + p4_b).Px() << " (" << p4_top.Px() << ")" << endl;
    cout << "Top Py   = " << (p4_q + p4_qbar + p4_b).Py() << " (" << p4_top.Py() << ")" << endl;
    cout << "Top Pz   = " << (p4_q + p4_qbar + p4_b).Pz() << " (" << p4_top.Pz() << ")" << endl;
    cout << "Higgs E  = " << (p4_bh + p4_bbarh ).E() << " (" << p4_h.E() << ")" << endl;
    cout << "Higgs Px = " << (p4_bh + p4_bbarh ).Px() << " (" << p4_h.Px() << ")" << endl;
    cout << "Higgs Py = " << (p4_bh + p4_bbarh ).Py() << " (" << p4_h.Py() << ")" << endl;
    cout << "Higgs Pz = " << (p4_bh + p4_bbarh ).Pz() << " (" << p4_h.Pz() << ")" << endl;
  }


  switch(type){
  case Decay::TopHad:
    out.push_back( make_pair('j', p4_b) );
    out.push_back( make_pair('j', p4_q) );
    out.push_back( make_pair('j', p4_qbar) );
    break;
  case Decay::TopLep:
    out.push_back( make_pair('j', p4_b) );
    out.push_back( make_pair('l', p4_q) );
    p4_qbar.SetPz(0.);
    out.push_back( make_pair('m', p4_qbar) );
    break;
  case Decay::WHad:
    out.push_back( make_pair('j', p4_q) );
    out.push_back( make_pair('j', p4_qbar) );
    break;
  case Decay::Higgs:
    out.push_back( make_pair('j', p4_bh) );
    out.push_back( make_pair('j', p4_bbarh) );
    break;
  case Decay::Radiation:
    out.push_back( make_pair('j', p4_q) );
    break;
  default:
    break;
  }

  return;
}




