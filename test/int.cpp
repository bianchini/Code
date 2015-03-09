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

double PI = TMath::Pi();

int main(){

  Integrand* integrand = new Integrand(DebugVerbosity::init
				       //|DebugVerbosity::init_more
				       //|DebugVerbosity::input
				       //|DebugVerbosity::integration
				       );

  integrand->push_back_object     ( TLorentzVector(50,0, 10, sqrt(50*50+10*10)),    ObjectType::Jet );
  integrand->add_object_observable( make_pair(Observable::BTAG, 0.), ObjectType::Jet);

  integrand->push_back_object     ( TLorentzVector(0,50, 20,  sqrt(50*50+20*20)),   ObjectType::Jet );
  integrand->add_object_observable( make_pair(Observable::BTAG, 1.), ObjectType::Jet);

  integrand->push_back_object     ( TLorentzVector(30,30,40,  sqrt(30*30+30*30+40*40)),   ObjectType::Jet );
  integrand->add_object_observable( make_pair(Observable::BTAG, 0.), ObjectType::Jet);

  integrand->push_back_object     ( TLorentzVector(70,20,10, sqrt(70*70+20*20+10*10)),   ObjectType::Jet );
  integrand->add_object_observable( make_pair(Observable::BTAG, 1.), ObjectType::Jet);

  integrand->push_back_object     ( TLorentzVector(20,50,10,  sqrt(20*20+50*50+10*10)),   ObjectType::Jet );
  integrand->add_object_observable( make_pair(Observable::BTAG, 1.), ObjectType::Jet);

  integrand->push_back_object     ( TLorentzVector(100,10,20, sqrt(100*100+10*10+20*20)),  ObjectType::Jet );
  integrand->add_object_observable( make_pair(Observable::BTAG, 1.), ObjectType::Jet);

  integrand->push_back_object     ( TLorentzVector(70,10.,20.,  sqrt(70*70+20*20+10*10)),      ObjectType::Lepton );
  integrand->add_object_observable( make_pair(Observable::CHARGE, +1.), ObjectType::Lepton);

  integrand->push_back_object     ( TLorentzVector(30,0,0,30),        ObjectType::MET );

  integrand->set_permutation_strategy( {Permutations::BTagged, Permutations::QUntagged, Permutations::QQbarSymmetry, Permutations::BBbarSymmetry});
  integrand->set_integrand(IntegrandType::Constant
			   |IntegrandType::Jacobian
			   |IntegrandType::ScattAmpl
			   |IntegrandType::DecayAmpl
			   |IntegrandType::PDF
			   |IntegrandType::Transfer
			   );

  //integrand->run( Hypothesis::TTH,  {} );
  integrand->run( Hypothesis::TTH,  {PSVar::cos_qbar1, PSVar::phi_qbar1} );
  //integrand->run( Hypothesis::TTBB, {} );
  integrand->next_event();

  delete integrand;
}
