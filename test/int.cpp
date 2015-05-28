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
  cfg.perm_filtering_rel = 1e-03;
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

  Object j1( TLorentzVector(30,30, 10, sqrt(30*30+10*10)), ObjectType::Jet );
  j1.addObs( Observable::BTAG, 0. );  
  j1.addObs( Observable::PDGID, 1 );  
  Object j2( TLorentzVector(20,50, 20, sqrt(50*50+20*20)), ObjectType::Jet );
  j2.addObs( Observable::BTAG, 1. );
  j2.addObs( Observable::PDGID, 5 );  
  Object j3( TLorentzVector(20,30, 40, sqrt(30*30+30*30+40*40)), ObjectType::Jet );
  j3.addObs( Observable::BTAG, 0. );
  j3.addObs( Observable::PDGID, 1 );  
  Object j4( TLorentzVector(10,20, 10, sqrt(70*70+20*20+10*10)), ObjectType::Jet );
  j4.addObs( Observable::BTAG, 1. );
  j4.addObs( Observable::PDGID, 22 );  
  Object j5( TLorentzVector(20,50, 10, sqrt(20*20+50*50+10*10)), ObjectType::Jet );
  j5.addObs( Observable::BTAG, 1. );
  j5.addObs( Observable::PDGID, 22 );   
  Object j6( TLorentzVector(100,10, 20, sqrt(100*100+10*10+20*20)), ObjectType::Jet );
  j6.addObs( Observable::BTAG, 1. );
  j6.addObs( Observable::PDGID, -5 );
  Object j7( TLorentzVector(100,30, 50, sqrt(100*100+30*30+50*50)), ObjectType::Jet );
  j7.addObs( Observable::BTAG, 0. );
  j7.addObs( Observable::PDGID, -1 );
  Object j8( TLorentzVector(100,-30, -50, sqrt(100*100+30*30+50*50)), ObjectType::Jet );
  j8.addObs( Observable::BTAG, 0. );
  j8.addObs( Observable::PDGID, -1 );

  Object l1( TLorentzVector(70,10, 20, sqrt(70*70+10*10+20*20)), ObjectType::Lepton );
  l1.addObs( Observable::CHARGE, +1. );
  Object l2( TLorentzVector(70,-10, -20, sqrt(70*70+10*10+20*20)), ObjectType::Lepton );
  l2.addObs( Observable::CHARGE, -1. );

  Object met( TLorentzVector(40,20,0,sqrt(40*40+20*20)), ObjectType::MET );

  integrand->push_back_object( &j1 );
  integrand->push_back_object( &j2 );
  integrand->push_back_object( &j3 );
  integrand->push_back_object( &j4 );
  integrand->push_back_object( &j5 );
  integrand->push_back_object( &j6 );
  //integrand->push_back_object( &j7 );
  //integrand->push_back_object( &j8 );
  integrand->push_back_object( &l1 );
  //integrand->push_back_object( &l2 );
  integrand->push_back_object( &met );

  
  integrand->set_permutation_strategy
    (  { Permutations::BTagged
	//,Permutations::QUntagged 
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

  //res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //res = integrand->run( FinalState::LH, Hypothesis::TTBB,  {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );

  res = integrand->run( FinalState::LH, Hypothesis::TTBB,  {PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );

  //integrand->run( FinalState::TTH, Hypothesis::TTH,  {} );
  integrand->next_event();

  delete integrand;
}
