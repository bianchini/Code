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

  Integrand* integrand = new Integrand(Hypothesis::TTH, 
				       //DebugVerbosity::init|DebugVerbosity::input|DebugVerbosity::integration
				       DebugVerbosity::init
				       );

  integrand->push_back_object( TLorentzVector(50,0, 10,  50),   ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(0,50, 20,  50),   ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(30,30,40,  40),   ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(70,10,20,  70),   ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(20,50,10,  50),   ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(100,10,20, 100),  ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(70,20,10,  70),   ObjectType::Lepton );
  integrand->push_back_object( TLorentzVector(30,30,0,30),      ObjectType::MET );
  integrand->push_back_object( TLorentzVector(-1,-1,0,2),       ObjectType::Recoil );

  integrand->run( Hypothesis::TTH,  {} );
  integrand->run( Hypothesis::TTH,  {PSVar::cos_qbar1, PSVar::phi_qbar1} );
  integrand->run( Hypothesis::TTBB, {} );
  integrand->next_event();

  delete integrand;
}
