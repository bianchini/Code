#include "interface/Integrand.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"

#include<iostream>

using namespace std;
using namespace MEM;

int main(){

  MEMConfig cfg;
  cfg.defaultCfg();
  //cfg.do_perm_filtering  = 1;
  cfg.do_prefit = 1;
  //cfg.perm_filtering_rel = 1e-03;
  //cfg.transfer_function_method = TFMethod::External;
  //cfg.do_minimize = 1;
  //cfg.perm_int    = 0;
  //cfg.defaultCfg(2.0);
  //cfg.tf_suppress = 0;
  //cfg.tf_offscale = 5.;
  //cfg.tf_in_range = true;
  //cfg.j_range_CL = 0.98;
  //cfg.b_range_CL = 0.99;
  //cfg.highpt_first = 0;


  Integrand* integrand = new Integrand(  DebugVerbosity::output
					 |DebugVerbosity::init
					 //|DebugVerbosity::input
					 //|DebugVerbosity::init_more
					 //|DebugVerbosity::input
					 //|DebugVerbosity::integration				        
					 ,cfg);
  TLorentzVector lv_j1;
  lv_j1.SetPtEtaPhiM(242.816604614, -0.107542805374, 1.25506973267, 24.5408706665);
  Object j1( lv_j1, ObjectType::Jet );
  j1.addObs( Observable::BTAG, 0. );  
  j1.addObs( Observable::PDGID, 1 );  

  TLorentzVector lv_j2;
  lv_j2.SetPtEtaPhiM(35.6511192322, 0.566395223141, -2.51394343376, 8.94268417358);
  Object j2( lv_j2, ObjectType::Jet );
  j2.addObs( Observable::BTAG, 1. );
  j2.addObs( Observable::PDGID, 5 );  

  TLorentzVector lv_j3;
  lv_j3.SetPtEtaPhiM(77.6708831787, -0.709680855274, -2.53739523888, 10.4904966354);
  Object j3( lv_j3, ObjectType::Jet );
  j3.addObs( Observable::BTAG, 0. );
  j3.addObs( Observable::PDGID, 1 );  

  TLorentzVector lv_j4;
  lv_j4.SetPtEtaPhiM(52.0134391785, -0.617823541164, -1.23360788822, 6.45914268494);
  Object j4( lv_j4, ObjectType::Jet );
  j4.addObs( Observable::BTAG, 1. );
  j4.addObs( Observable::PDGID, 22 );  

  TLorentzVector lv_j5;
  lv_j5.SetPtEtaPhiM( 235.892044067, -0.997860729694, -2.10646605492, 27.9887943268 );
  Object j5( lv_j5, ObjectType::Jet );
  j5.addObs( Observable::BTAG, 1. );
  j5.addObs( Observable::PDGID, 22 );   

  TLorentzVector lv_j6;
  lv_j6.SetPtEtaPhiM(191.423553467, -0.46368226409, 0.750520706177, 30.5682048798);
  Object j6( lv_j6, ObjectType::Jet );
  j6.addObs( Observable::BTAG, 1. );
  j6.addObs( Observable::PDGID, -5 );

  TLorentzVector lv_j7;
  lv_j6.SetPtEtaPhiM(30, 0, 0, 0.);
  Object j7( lv_j7, ObjectType::Jet );
  j7.addObs( Observable::BTAG, 0. );
  j7.addObs( Observable::PDGID, -1 );

  TLorentzVector lv_j8;
  lv_j8.SetPtEtaPhiM(50, 1, -1, 0.);
  Object j8( lv_j8, ObjectType::Jet );
  j8.addObs( Observable::BTAG, 0. );
  j8.addObs( Observable::PDGID, -1 );

  TLorentzVector lv_l1;
  lv_l1.SetPtEtaPhiM(52.8751449585, -0.260020583868, -2.55171084404, 0.139569997787);
  Object l1( lv_l1, ObjectType::Lepton );
  l1.addObs( Observable::CHARGE, +1. );
  Object l2( TLorentzVector(70,-10, -20, sqrt(70*70+10*10+20*20)), ObjectType::Lepton );
  l2.addObs( Observable::CHARGE, -1. );

  TLorentzVector lv_met;
  lv_met.SetPtEtaPhiM( 92.1731872559,0., -1.08158898354, 0.);
  Object met( lv_met, ObjectType::MET );
  /*
  integrand->push_back_object( &j1 );
  integrand->push_back_object( &j2 );
  integrand->push_back_object( &j3 );
  integrand->push_back_object( &j4 );
  integrand->push_back_object( &j5 );
  integrand->push_back_object( &j6 );
  integrand->push_back_object( &j7 );
  //integrand->push_back_object( &j8 );
  integrand->push_back_object( &l1 );
  //integrand->push_back_object( &l2 );
  integrand->push_back_object( &met );
  */
  
  integrand->set_permutation_strategy
    (  { Permutations::BTagged
	,Permutations::QUntagged 
	,Permutations::QQbarBBbarSymmetry
	//,Permutations::HEPTopTagged
	//,Permutations::HiggsTagged
	} 
      );

  integrand->set_integrand(IntegrandType::Constant
			   |IntegrandType::ScattAmpl
			   |IntegrandType::DecayAmpl
			   |IntegrandType::Jacobian
			   |IntegrandType::PDF
			   |IntegrandType::Transfer
			   //|IntegrandType::SmearJets
			   //|IntegrandType::SmearMET
			   //|IntegrandType::Sudakov
			   //|IntegrandType::Recoil
			   );
  //integrand->set_integrand(0);

  //integrand->set_ncalls(16000);
  //integrand->set_sqrts (13000.);  

  MEMOutput res;			   
  //res = integrand->run( FinalState::LH, Hypothesis::TTH,  {} );
  //res = integrand->run( FinalState::LH, Hypothesis::TTBB, {} );
  //res = integrand->run( FinalState::HH, Hypothesis::TTH,  {} );
  //res = integrand->run( FinalState::LL, Hypothesis::TTH,  {} );

  //res = integrand->run( FinalState::LH, Hypothesis::TTH,  {}, {PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //res = integrand->run( FinalState::LH, Hypothesis::TTBB,  {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );
  for(int i = 0; i < 100; ++i){

    integrand->push_back_object( &j1 );
    integrand->push_back_object( &j2 );
    integrand->push_back_object( &j3 );
    integrand->push_back_object( &j4 );
    integrand->push_back_object( &j5 );
    integrand->push_back_object( &j6 );
    integrand->push_back_object( &j7 );
    //integrand->push_back_object( &j8 );
    integrand->push_back_object( &l1 );
    //integrand->push_back_object( &l2 );
    integrand->push_back_object( &met );
    
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {} );
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {} ); 
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_qbar1, PSVar::phi_qbar1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_qbar1, PSVar::phi_qbar1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_b1,    PSVar::phi_b1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_b1,    PSVar::phi_b1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_b2,    PSVar::phi_b2} );
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_b2,    PSVar::phi_b2} );
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_qbar1, PSVar::phi_qbar1, PSVar::cos_b1, PSVar::phi_b1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_qbar1, PSVar::phi_qbar1, PSVar::cos_b1, PSVar::phi_b1} );
    res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_qbar1, PSVar::phi_qbar1, PSVar::cos_b2, PSVar::phi_b2} );
    res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_qbar1, PSVar::phi_qbar1, PSVar::cos_b2, PSVar::phi_b2} );
    integrand->next_event();
  }
  //res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_b1, PSVar::phi_b1} );

  //integrand->run( FinalState::TTH, Hypothesis::TTH,  {} );
  integrand->next_event();

  delete integrand;
}
