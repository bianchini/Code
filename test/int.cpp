#include "interface/Integrand.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"

#include<iostream>

using namespace std;

double PI = TMath::Pi();

int main(){

  MEM::Integrand* integrand = new MEM::Integrand(MEM::Hypothesis::TTH, MEM::DebugVerbosity::init|MEM::DebugVerbosity::input|MEM::DebugVerbosity::integration);
  integrand->push_back_object( TLorentzVector(50,0, 10,  50),   MEM::ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(0,50, 20,  50),   MEM::ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(30,30,40,  40),   MEM::ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(70,10,20,  70),   MEM::ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(20,50,10,  50),   MEM::ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(100,10,20, 100),  MEM::ObjectType::Jet );
  integrand->push_back_object( TLorentzVector(70,20,10,  70),   MEM::ObjectType::Lepton );
  integrand->push_back_object( TLorentzVector(30,30,0,30),      MEM::ObjectType::MET );
  integrand->push_back_object( TLorentzVector(-1,-1,0,2),       MEM::ObjectType::Recoil );

  integrand->init();
  integrand->run();

  delete integrand;

}
