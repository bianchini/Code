#ifndef EVENT_H
#define EVENT_H

#include <TTree.h>
#include <iostream>

#define HMAX 10
#define PMAX 100


using namespace std;

namespace Algo {

  struct TreeStruct {    
    std::size_t n_h;             // num of hypo
    std::size_t n_dim;           // total number of dimensions
    int    all_time;        // total CPU time (msec)
    int    n_btag;          // n btag
    int    n_jet;           // n jet
    int    n_lep;           // n lep
    double nll      [HMAX]; // nll per hypo
    int    status   [HMAX]; // status per hypo
    int    strategy [HMAX]; // strategy per hypo
    int    min_time [HMAX]; // minimization time per hypo
    int    dim      [HMAX]; // num of parameters per hypo 
    int    perm     [HMAX]; // num of permutations per hypo
    double param    [PMAX]; // minimzed parameters [ [p_0,...,p_dim[0]], ... ]
    double obs_e    [PMAX]; // input observables
    double obs_pt   [PMAX]; // input observables
    double obs_eta  [PMAX]; // input observables
    double obs_phi  [PMAX]; // input observables
    double obs_btag [PMAX]; // input observables
    void print(ostream& os){
      os << "TreeStruct contains: "   << endl;
      os << "\tall_time = " << all_time << " msec" << endl;
      os << "\tn_h      = " << n_h    << endl;
      os << "\tn_dim    = " << n_dim  << endl;
      os << "\tn_btag   = " << n_btag << endl;
      os << "\tn_jet    = " << n_jet  << endl;
      os << "\tn_lep    = " << n_lep  << endl;
      std::size_t it_p = 0;
      for( std::size_t i = 0 ; (i < n_h && i < (std::size_t)HMAX) ; ++i){
	os << "\tHypo "     << i << ": "   << endl;
	os << "\tstatus["     << i << "] = " << status[i];
	os << ", strategy[" << i << "] = " << strategy[i];
	os << ", perm["     << i << "] = " << perm[i];
	os << ", nll["      << i << "] = " << nll[i];
	os << ", min_time[" << i << "] = " << min_time[i];
	os << ", dim["      << i << "] = " << dim[i] << endl;
	for(int j = 0 ; j < dim[i] ; ++j){	  
	  os << "\t\tparam["  << j << "] = " << param[it_p] 
	     << ", obs_e["    << j << "] = " << obs_e[it_p]   
	     << ", obs_pt["   << j << "] = " << obs_pt[it_p]   
	     << ", obs_eta["  << j << "] = " << obs_eta[it_p]   
	     << ", obs_phi["  << j << "] = " << obs_phi[it_p]   
	     << ", obs_btag[" << j << "] = " << obs_btag[it_p] << endl; 
	  ++it_p;
	}
      }
    }
  };

  struct Event {
    Event (TTree*) ;
    Event();
    ~Event();
    void fillTree () ;
    void printTree() ;
    void createBranches () ;
    void reset();
    TreeStruct treeStruct ;
    private :
    TTree* outTree ;
  };

}

#endif
