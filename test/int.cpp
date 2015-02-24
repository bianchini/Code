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
  integrand->push_back_object( TLorentzVector(1,0,0,1),   MEM::Object::Type::jet );
  integrand->push_back_object( TLorentzVector(1,1,0,2),   MEM::Object::Type::jet );
  integrand->push_back_object( TLorentzVector(1,0,1,2),   MEM::Object::Type::jet );
  integrand->push_back_object( TLorentzVector(1,1,0,2),   MEM::Object::Type::jet );
  //integrand->push_back_object( TLorentzVector(1,1,0,2),   MEM::Object::Type::jet );
  //integrand->push_back_object( TLorentzVector(1,1,0,2),   MEM::Object::Type::jet );
  integrand->push_back_object( TLorentzVector(1,0,1,2),   MEM::Object::Type::lepton );
  integrand->push_back_object( TLorentzVector(1,1,0,2),   MEM::Object::Type::met );
  integrand->push_back_object( TLorentzVector(-1,-1,0,2), MEM::Object::Type::recoil );

  integrand->init();

  /*
  ROOT::Math::Functor toIntegrate(integrand, &MEM::Integrand::Eval, 1);
    
  ROOT::Math::GSLMCIntegrator* ig2 = 
    new ROOT::Math::GSLMCIntegrator(
				    ROOT::Math::IntegrationMultiDim::kVEGAS,
				    1.e-12, //absolute tolerance
				    1.e-5,  //relative tolerance
				    4000    //maximum number of calls
				    );


  ig2->SetFunction(toIntegrate);
  const double xL[24] = {0.};
  const double xU[24] = {1.};

  cout << ig2->Integral(xL,xU) << endl;
  */
  integrand->run();

  delete integrand;

}
