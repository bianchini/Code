#include "interface/Integrand.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/IOptions.h"
#include "Math/IntegratorOptions.h"
#include "Math/AllIntegrationTypes.h"
#include "TFile.h"
#include <iostream>

using namespace std;
using namespace MEM;

// Returns the transfer function corresponding to a jet flavour and eta
TF1* getTransferFunction(TFile* tffile, const char* flavour, double eta) {
    int etabin = 0;
    if (std::abs(eta) > 1.0) {
        etabin = 1;
    }
    stringstream ss;
    ss << "tf_" << flavour << "_etabin" << etabin;
    const char* fname = ss.str().c_str();
    TF1* tf = (TF1*)(tffile->Get(fname));
    if (tf == 0) {
        cerr << "could not get transfer function " << fname << endl;
        cerr << flush;
        throw exception();
    }
    return tf;
}

int main(){

  //Load the transfer functions
  TFile* tffile = new TFile("transfer.root");


  //create a MEM configuration.
  //this needs to be done once per job, not for every event
  MEMConfig cfg;
  cfg.defaultCfg();
  cfg.transfer_function_method = TFMethod::Builtin;

  //Transfer functions for jet reconstruction efficiency
  cfg.set_tf_global(TFType::bLost, 0, *getTransferFunction(tffile, "beff", 0.0));
  cfg.set_tf_global(TFType::bLost, 1, *getTransferFunction(tffile, "beff", 2.0));
  cfg.set_tf_global(TFType::qLost, 0, *getTransferFunction(tffile, "leff", 0.0));
  cfg.set_tf_global(TFType::qLost, 1, *getTransferFunction(tffile, "leff", 2.0));

  //Create the mem integrator, once per job
  Integrand* integrand = new Integrand( 
    DebugVerbosity::output
    //|DebugVerbosity::init
    //|DebugVerbosity::input
    //|DebugVerbosity::init_more
    //|DebugVerbosity::integration				        
    ,cfg
  );

  //Add some objects to the MEM
  TLorentzVector lv_j1;
  lv_j1.SetPtEtaPhiM(242.816604614, -0.107542805374, 1.25506973267, 24.5408706665);
  Object j1( lv_j1, ObjectType::Jet );
  j1.addObs( Observable::BTAG, 0. ); // 0 - jet is assumed to be from a light quark, 1 - a b quark
  j1.addObs( Observable::PDGID, 1 );  // currently not used
  // attach the transfer functions corresponding to the jet
  j1.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j1.Eta()));
  j1.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j1.Eta()));

  TLorentzVector lv_j2;
  lv_j2.SetPtEtaPhiM(191.423553467, -0.46368226409, 0.750520706177, 30.5682048798);
  Object j2( lv_j2, ObjectType::Jet );
  j2.addObs( Observable::BTAG, 1. );
  j2.addObs( Observable::PDGID, 5 );  
  j2.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j2.Eta()));
  j2.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j2.Eta()));

  TLorentzVector lv_j3;
  lv_j3.SetPtEtaPhiM(77.6708831787, -0.709680855274, -2.53739523888, 10.4904966354);
  Object j3( lv_j3, ObjectType::Jet );
  j3.addObs( Observable::BTAG, 0. );
  j3.addObs( Observable::PDGID, 1 );  
  j3.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j3.Eta()));
  j3.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j3.Eta()));

  TLorentzVector lv_j4;
  lv_j4.SetPtEtaPhiM(235.892044067, -0.997860729694, -2.10646605492, 27.9887943268);
  Object j4( lv_j4, ObjectType::Jet );
  j4.addObs( Observable::BTAG, 1. );
  j4.addObs( Observable::PDGID, 22 );  
  j4.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j4.Eta()));
  j4.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j4.Eta()));

  TLorentzVector lv_j5;
  lv_j5.SetPtEtaPhiM( 52.0134391785, -0.617823541164, -1.23360788822, 6.45914268494 );
  Object j5( lv_j5, ObjectType::Jet );
  j5.addObs( Observable::BTAG, 1. );
  j5.addObs( Observable::PDGID, 22 );   
  j5.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j5.Eta()));
  j5.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j5.Eta()));

  TLorentzVector lv_j6;
  lv_j6.SetPtEtaPhiM(35.6511192322, 0.566395223141, -2.51394343376, 8.94268417358);
  Object j6( lv_j6, ObjectType::Jet );
  j6.addObs( Observable::BTAG, 1. );
  j6.addObs( Observable::PDGID, -5 );
  j6.addTransferFunction(TFType::bReco, getTransferFunction(tffile, "b", lv_j6.Eta()));
  j6.addTransferFunction(TFType::qReco, getTransferFunction(tffile, "l", lv_j6.Eta()));

  //create a lepton
  TLorentzVector lv_l1;
  lv_l1.SetPtEtaPhiM(52.8751449585, -0.260020583868, -2.55171084404, 0.139569997787);
  Object l1( lv_l1, ObjectType::Lepton );
  l1.addObs( Observable::CHARGE, +1. );
  // Object l2( TLorentzVector(70,-10, -20, sqrt(70*70+10*10+20*20)), ObjectType::Lepton );
  // l2.addObs( Observable::CHARGE, -1. );

  //create a MET
  TLorentzVector lv_met;
  lv_met.SetPtEtaPhiM( 92.1731872559,0., -1.08158898354, 0.);
  Object met( lv_met, ObjectType::MET );
  
  //add all objects to the MEM integrator
  integrand->push_back_object( &j1 );
  integrand->push_back_object( &j2 );
  integrand->push_back_object( &j3 );
  integrand->push_back_object( &j4 );
  integrand->push_back_object( &j5 );
  integrand->push_back_object( &j6 );
  integrand->push_back_object( &l1 );
  integrand->push_back_object( &met );
  
  //permute light quarks only in untagged jets, b-quarks in tagged
  //remove symmetric permutations
  integrand->set_permutation_strategy
    (  {Permutations::BTagged
        ,Permutations::QUntagged 
        ,Permutations::QQbarBBbarSymmetry
        //,Permutations::HEPTopTagged
        //,Permutations::HiggsTagged
        } 
    );

  MEMOutput res;			   
  
  //Evaluate fully reconstructed hypothesis
  //LH - single-leptonic decay channel
  //LL - dileptonic decay channel
  cout << "Fully reconstructed interpretation" << endl;
  cout << "evaluating tth hypo" << endl;
  //third variable is variables to integrate over
  //if nothing is specified, assume that all jets (4b + 2 light) were reconstructed
  res = integrand->run( FinalState::LH, Hypothesis::TTH,  {} );
  cout << "p = " << res.p << " +- " << res.p_err << endl;
  double p0 = res.p;

  cout << "evaluating ttbb hypo" << endl;
  res = integrand->run( FinalState::LH, Hypothesis::TTBB, {} );
  cout << "p = " << res.p << " +- " << res.p_err << endl;
  double p1 = res.p;

  double mem_w = p0 / (p0 + 0.02*p1);
  cout << "mem 222 discriminator " << mem_w << endl;

  //Evaluate 022 hypothesis. We do not use the information provided by the light quarks.
  cout << "Integrating over light quarks" << endl;
  cout << "evaluating tth hypo" << endl;
  //integrate over the light quark angles
  res = integrand->run( FinalState::LH, Hypothesis::TTH,  {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );
  cout << "p = " << res.p << " +- " << res.p_err << endl;
  p0 = res.p;
  
  cout << "evaluating ttbb hypo" << endl;
  res = integrand->run( FinalState::LH, Hypothesis::TTBB, {PSVar::cos_q1, PSVar::phi_q1, PSVar::cos_qbar1, PSVar::phi_qbar1} );
  cout << "p = " << res.p << " +- " << res.p_err << endl;
  p1 = res.p;
  
  mem_w = p0 / (p0 + 0.02*p1);
  cout << "mem 022 discriminator " << mem_w << endl;

  integrand->next_event();

  delete integrand;
}
