#include<iostream>
#include<memory>
#include<cstdlib>
#include <getopt.h>

#include "interface/ToyGenerator.h"
#include "interface/HypoTester.h"

#include "TMath.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"


enum GEN_HYPO  : int { GEN_tH_wH=0 ,  GEN_tH_tL=1,   GEN_tH_tL_hH=2,  GEN_tL_tL_hH=3,  GEN_tH_tH_hH=4,
		     GEN_tH_bb=5 ,  GEN_udb_tL=6,  GEN_tH_tL_bb=7,  GEN_tL_tL_bb=8,  GEN_tH_tH_bb=9,
		     GEN_tL_ggggg=10, GEN_lm_gggggg=11, 
		     };

enum TEST_HYPO : int { TEST_tH_wH=0 , TEST_tH_tL=1, TEST_tH_tL_hH=2, TEST_tL_tL_hH=3, TEST_tH_tH_hH=4,
		     TEST_tL=5};


using namespace std;

int main(int argc, char *argv[]){

  string fname = "tmp.root";
  int ntoys    = 1;
  int pass     = 1;  
  int smear    = 0;
  int btag     = 0;
  int verbose  = 0;
  int debug    = -1;
  int seed     = 4357; // default seed
  GEN_HYPO  gen_hypo  = GEN_HYPO::GEN_tH_tL;
  TEST_HYPO test_hypo = TEST_HYPO::TEST_tH_tL;

  static struct option long_options[] = {
    {"toys",     optional_argument,   0, 't' },
    {"output",   optional_argument,   0, 'o' },
    {"smear",    optional_argument,   0, 'g' },
    {"btag",     optional_argument,   0, 'b' },
    {"verbose",  optional_argument,   0, 'v' },
    {"seed",     optional_argument,   0, 's' },
    {"pass",     optional_argument,   0, 'p' },
    {"help",     optional_argument,   0, 'h' },
    {"debug",    optional_argument,   0, 'd' },
    {"gen",      optional_argument,   0, 'G' },
    {"test",     optional_argument,   0, 'T' },
    {0, 0, 0, 0 }
  };

  int long_index  = 0;
  int opt;
  while ((opt = getopt_long(argc, argv,"t:o:gb:v:s:phd:G:T:",
			    long_options, &long_index )) != -1) {
    switch (opt) {
    case 't' : ntoys   = atoi(optarg);
      break;
    case 'o' : fname   = string(optarg);
      break;
    case 'g' : smear   = 1;
      break;
    case 'b' : btag    = atoi(optarg);
      break;
    case 'v' : verbose = atoi(optarg);                                 
      break;  
    case 's' : seed    = atoi(optarg);
      break;
    case 'p' : pass    = 1;
      break;
    case 'd' : debug   = atoi(optarg);
      break;
    case 'G' : gen_hypo  = static_cast<GEN_HYPO>(atoi(optarg));
      break;
    case 'T' : test_hypo = static_cast<TEST_HYPO>(atoi(optarg));
      break;
    case 'h' :
      cout << "Usage: toy [-t TOYS] [-o OUTFILENAME.root] [-g] [-b BTAG] [-v VERBOSE] [-s SEED] [-p] [-d DEBUGEVTNUM] [-G GENHYPO] [-T TESTHYPO]" << endl; 
      return 0;
    case '?':
      cout << "Unknown option " << opt << endl;
      break;
    default:
      cout << "Invalid usage of toy" << endl;
      return 0;
    }
  }


  TTree* t = new TTree("tree", "tree");
  Algo::HypoTester* tester = new Algo::HypoTester(t) ;
  tester->set_verbosity(verbose);

  Algo::ToyGenerator* toyGenerator = new Algo::ToyGenerator(verbose, seed);

  vector<Algo::Decay> decays;

  switch( gen_hypo ){
  case GEN_tH_wH:
    decays.push_back( Algo::Decay::TopHad );
    decays.push_back( Algo::Decay::Higgs ); 
    break;
  case GEN_tH_tL:
    decays.push_back( Algo::Decay::TopHad );  
    decays.push_back( Algo::Decay::TopLep );   
    break;
  case GEN_tH_tL_hH:
    decays.push_back( Algo::Decay::TopHad );
    decays.push_back( Algo::Decay::TopLep ); 
    decays.push_back( Algo::Decay::Higgs );
    break;
  case GEN_tL_tL_hH:
    decays.push_back( Algo::Decay::TopLep );  
    decays.push_back( Algo::Decay::TopLep );  
    decays.push_back( Algo::Decay::Higgs );   
    break;
  case GEN_tH_tH_hH:
    decays.push_back( Algo::Decay::TopHad );  
    decays.push_back( Algo::Decay::TopHad );  
    decays.push_back( Algo::Decay::Higgs );  
    break;
  case GEN_tH_bb:
    decays.push_back( Algo::Decay::TopHad ); 
    decays.push_back( Algo::Decay::Radiation_b );  
    decays.push_back( Algo::Decay::Radiation_b );  
    break;
  case GEN_udb_tL:
    decays.push_back( Algo::Decay::Radiation_u );       
    decays.push_back( Algo::Decay::Radiation_d );        
    decays.push_back( Algo::Decay::Radiation_b ); 
    decays.push_back( Algo::Decay::TopLep );  
    break;
  case GEN_tH_tL_bb:
    decays.push_back( Algo::Decay::TopHad );  
    decays.push_back( Algo::Decay::TopLep );        
    decays.push_back( Algo::Decay::Radiation_b );     
    decays.push_back( Algo::Decay::Radiation_b );  
    break;
  case GEN_tL_tL_bb:
    decays.push_back( Algo::Decay::TopLep ); 
    decays.push_back( Algo::Decay::TopLep );  
    decays.push_back( Algo::Decay::Radiation_b );     
    decays.push_back( Algo::Decay::Radiation_b );  
    break;
  case GEN_tH_tH_bb:
    decays.push_back( Algo::Decay::TopHad ); 
    decays.push_back( Algo::Decay::TopHad ); 
    decays.push_back( Algo::Decay::Radiation_b );  
    decays.push_back( Algo::Decay::Radiation_b );  
    break;
  case GEN_tL_ggggg:
    decays.push_back( Algo::Decay::TopLep );
    decays.push_back( Algo::Decay::Radiation_g );
    decays.push_back( Algo::Decay::Radiation_g );
    decays.push_back( Algo::Decay::Radiation_g );
    //decays.push_back( Algo::Decay::Radiation_g );
    //decays.push_back( Algo::Decay::Radiation_g );
    break;
  case GEN_lm_gggggg:
    decays.push_back( Algo::Decay::Lepton );
    decays.push_back( Algo::Decay::MET );
    decays.push_back( Algo::Decay::Radiation_g );
    decays.push_back( Algo::Decay::Radiation_g );
    decays.push_back( Algo::Decay::Radiation_g );
    decays.push_back( Algo::Decay::Radiation_g );
    //decays.push_back( Algo::Decay::Radiation_g );
    //decays.push_back( Algo::Decay::Radiation_g );
    break;
  default:
    break;
  }

  int itoy = 0;
  while( itoy < ntoys ){
    vector<Algo::Object> out = 
      toyGenerator->generate( decays , smear, btag );

    int count_j = 0;
    int count_l = 0;
    int count_m = 0;
    TLorentzVector invisible(0.,0.,0.,0.);

    for( auto fs : out ){
      if( (fs.type=='q' || fs.type=='b') && fs.p4.Pt()>30 && TMath::Abs(fs.p4.Eta())<2.5 ) ++count_j;
      if( fs.type=='l' && fs.p4.Pt()>20 && TMath::Abs(fs.p4.Eta())<2.5 ) ++count_l;
      if( fs.type=='m' && fs.p4.Pt()>0. ){
	invisible += fs.p4;
	++count_m;
      }
    }

    if(verbose>0){
      cout << "*************" << endl;
      for( auto fs : out ){
	cout <<  fs.type << ": " << fs.p4.Pt() << "," << fs.p4.Eta() << "," << fs.p4.Phi() << endl;
      }
    }
  
    int overlap {0};
    for( size_t i = 0 ; i < out.size()-1 ; i++ ){
      if( !( ((out[i].type=='q' || out[i].type=='b') && out[i].p4.Pt()>30 && TMath::Abs(out[i].p4.Eta())<2.5) ||
	     (out[i].type=='l' && out[i].p4.Pt()>20 && TMath::Abs(out[i].p4.Eta())<2.5)
	     ) ) continue;
      for( size_t j = i+1 ; j < out.size() ; j++ ){
	if( !( ((out[j].type=='q' || out[j].type=='b') && out[j].p4.Pt()>30 && TMath::Abs(out[j].p4.Eta())<2.5) || 
	       (out[j].type=='l' && out[j].p4.Pt()>20 && TMath::Abs(out[j].p4.Eta())<2.5)
	       ) ) continue;
	if( TMath::Sqrt( TMath::Power(out[i].p4.Eta()-out[j].p4.Eta(), 2) + TMath::Power(out[i].p4.Phi()-out[j].p4.Phi(), 2) )<0.4 ) {
	  ++overlap;
	  if(verbose>0){
	    cout << "Overlap: " << i << " and " << j << endl;
	    cout << out[i].type << " and " << out[j].type << ": (" << out[i].p4.Eta() << ", " << out[j].p4.Eta() << ")" 
		 << ", (" << out[i].p4.Phi() << ", " << out[j].p4.Phi() << ")" << endl;	    
	  }
	}
      }
    }
    


    if(pass==0) cout << "Generate event " << itoy+1 << "/" << ntoys << endl;
    if(verbose>0){
      cout << "\tNumber of jets: " << count_j << endl;
      cout << "\tNumber of lept: " << count_l << endl;
      cout << "\tNumber of nus : " << count_m << endl;
    }

    // TopHad + TopLep
    if(test_hypo==TEST_HYPO::TEST_tH_tL && count_j==4 && count_l==1 && count_m==1 && overlap==0) {
      if(pass){
	++itoy;
	cout << "Generate event " << itoy << "/" << ntoys << endl;
      }

      if( debug>=0 && itoy!=debug ) continue;
      
      // fill jets                                                                                                             
      for( auto fs : out ){
	if( (fs.type=='q' || fs.type=='b') && fs.p4.Pt()>30 && TMath::Abs(fs.p4.Eta())<2.5 ){
          tester->push_back_object( fs.p4  , 'j');
          if(btag)   tester->add_object_observables( "BTAG",     fs.obs["BTAG"]     , 'j');	 
	  if(btag>1) tester->add_object_observables( "BTAG_RND", fs.obs["BTAG_RND"] , 'j');
        }
	if( fs.type=='l' )
	  tester->push_back_object( fs.p4  , 'l');
      }
      tester->push_back_object( invisible  , 'm');

      map<string, vector<Algo::Decay> > hypotheses;
      hypotheses["H0"] = {Algo::Decay::TopHad, Algo::Decay::TopLep};
      hypotheses["H1"] = {Algo::Decay::TopLep, Algo::Decay::Radiation_u, Algo::Decay::Radiation_d, Algo::Decay::Radiation_b};
      if(verbose>0) tester->print(cout);
      tester->test( hypotheses );
    }


    // TopHad + TopLep + Higgs
    if( test_hypo==TEST_HYPO::TEST_tH_tL_hH && count_j==6 && count_l==1 && count_m==1 && overlap==0) {
      if(pass){
	++itoy;
	cout << "Generate event " << itoy << "/" << ntoys << endl;
      }
      
      if( debug>=0 && itoy!=debug ) continue;

      // fill jets                                                                                                             
      for( auto fs : out ){
	if( (fs.type=='q' || fs.type=='b') && fs.p4.Pt()>30 && TMath::Abs(fs.p4.Eta())<2.5 ){
          tester->push_back_object( fs.p4  , 'j');
          if(btag)   tester->add_object_observables( "BTAG",     fs.obs["BTAG"] ,     'j');
	  if(btag>1) tester->add_object_observables( "BTAG_RND", fs.obs["BTAG_RND"] , 'j');
        }
	if( fs.type=='l' )
	  tester->push_back_object( fs.p4  , 'l');
      }
      tester->push_back_object( invisible  , 'm');

      map<string, vector<Algo::Decay> > hypotheses;
      hypotheses["H0"] = {Algo::Decay::TopHad, Algo::Decay::TopLep, Algo::Decay::Higgs};
      hypotheses["H1"] = {Algo::Decay::TopHad, Algo::Decay::TopLep, Algo::Decay::Radiation_b, Algo::Decay::Radiation_b};
      if(verbose>0) tester->print(cout);
      tester->test( hypotheses );
    }

    // TopLep + TopLep + Higgs                                                                                                                       
    if( test_hypo==TEST_HYPO::TEST_tL_tL_hH && count_j==4 && count_l==2 && count_m==2 && overlap==0) {
      if(pass){
        ++itoy;
        cout << "Generate event " << itoy << "/" << ntoys << endl;
      }

      if( debug>=0 && itoy!=debug ) continue;
      
      // fill jets                                                                                                                                    
      for( auto fs : out ){
        if( (fs.type=='q' || fs.type=='b') && fs.p4.Pt()>30 && TMath::Abs(fs.p4.Eta())<2.5 ){
          tester->push_back_object( fs.p4  , 'j');
          if(btag)   tester->add_object_observables( "BTAG",     fs.obs["BTAG"] ,     'j');
	  if(btag>1) tester->add_object_observables( "BTAG_RND", fs.obs["BTAG_RND"] , 'j');
        }
        if( fs.type=='l' )
          tester->push_back_object( fs.p4  , 'l');
      }
      tester->push_back_object( invisible  , 'm');

      map<string, vector<Algo::Decay> > hypotheses;
      hypotheses["H0"] = {Algo::Decay::TopLep, Algo::Decay::TopLep, Algo::Decay::Higgs};
      hypotheses["H1"] = {Algo::Decay::TopLep, Algo::Decay::TopLep, Algo::Decay::Radiation_b, Algo::Decay::Radiation_b};
      if(verbose>0) tester->print(cout);
      tester->test( hypotheses );
    }

    // TopHad + TopHad + Higgs
    if( test_hypo==TEST_HYPO::TEST_tH_tH_hH && count_j==8 && count_l==0 && count_m==0 && overlap==0 ) {
      if(pass){
        ++itoy;
        cout << "Generate event " << itoy << "/" << ntoys << endl;
      }
      
      if( debug>=0 && itoy!=debug ) continue;

      // fill jets                                                                                                                            
      for( auto fs : out ){
        if( (fs.type=='q' || fs.type=='b') && fs.p4.Pt()>30 && TMath::Abs(fs.p4.Eta())<2.5 ){
          tester->push_back_object( fs.p4  , 'j');
          if(btag)   tester->add_object_observables( "BTAG",     fs.obs["BTAG"] ,     'j');
	  if(btag>1) tester->add_object_observables( "BTAG_RND", fs.obs["BTAG_RND"] , 'j');
        }
      }

      map<string, vector<Algo::Decay> > hypotheses;
      hypotheses["H0"] = {Algo::Decay::TopHad, Algo::Decay::TopHad, Algo::Decay::Higgs};
      hypotheses["H1"] = {Algo::Decay::TopHad, Algo::Decay::TopHad, Algo::Decay::Radiation_b, Algo::Decay::Radiation_b};
      if(verbose>0) tester->print(cout);
      tester->test( hypotheses );
    }


    // TopHad + WHad
    else if( test_hypo==TEST_HYPO::TEST_tH_wH && count_j==5 && count_m==0 && count_l==0 && overlap==0) {
      if(pass){
	++itoy;
	cout << "Generate event " << itoy << "/" << ntoys << endl;
      }

      if( debug>=0 && itoy!=debug ) continue;

      // fill jets
      for( auto fs : out ){
	if( (fs.type=='q' || fs.type=='b') && fs.p4.Pt()>30 && TMath::Abs(fs.p4.Eta())<2.5 ){
	  tester->push_back_object( fs.p4  , 'j');
          if(btag)   tester->add_object_observables( "BTAG",     fs.obs["BTAG"] ,     'j');
	  if(btag>1) tester->add_object_observables( "BTAG_RND", fs.obs["BTAG_RND"] , 'j');
	}
      }

      map<string, vector<Algo::Decay> > hypotheses;
      hypotheses["H0"] = {Algo::Decay::TopHad, Algo::Decay::Higgs};
      hypotheses["H1"] = {Algo::Decay::TopHad, Algo::Decay::Radiation_b, Algo::Decay::Radiation_b};
      if(verbose>0) tester->print(cout);
      tester->test( hypotheses );
    }

    //TopLep + other jets

    else if(test_hypo==TEST_HYPO::TEST_tL && count_j==4 && count_m>=0 && count_l>=0 && overlap==0){
      if(pass){
        ++itoy;
        cout << "Generate event " << itoy << "/" << ntoys << endl;
      }

      if( debug>=0 && itoy!=debug ) continue;
      
      // fill jets                                                                                                                             
      for( auto fs : out ){
        if( (fs.type=='q' || fs.type=='b') && fs.p4.Pt()>30 && TMath::Abs(fs.p4.Eta())<2.5 ){
          tester->push_back_object( fs.p4  , 'j');
          if(btag)   tester->add_object_observables( "BTAG",     fs.obs["BTAG"] ,     'j');
          if(btag>1) tester->add_object_observables( "BTAG_RND", fs.obs["BTAG_RND"] , 'j');
        }
	if( fs.type=='l' )
          tester->push_back_object( fs.p4  , 'l');
      }
      tester->push_back_object( invisible  , 'm');

      map<string, vector<Algo::Decay> > hypotheses;
      hypotheses["H0"] = { Algo::Decay::TopLep };
      hypotheses["H1"] = { Algo::Decay::Radiation_g, Algo::Decay::Radiation_g,Algo::Decay::Radiation_g,
			   Algo::Decay::Radiation_g /*, Algo::Decay::Radiation_g,Algo::Decay::Radiation_g*/};
      if(verbose>0) tester->print(cout);
      tester->test( hypotheses );
    }


    // other final state
    else{ /* ... */ }


    if(pass==0) ++itoy;

  }

  TFile* fout = new TFile(("test/"+fname).c_str(), "RECREATE");
  fout->cd();
  t->Write("", TObject::kOverwrite);
  fout->Close();
  if(verbose>0) cout << "ROOT file saved" << endl;

  delete toyGenerator;
  delete tester;
}
