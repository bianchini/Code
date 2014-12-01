#include <TTree.h>
#include <iostream>

#define HMAX 10
#define PMAX 100


using namespace std;

namespace Algo {

  struct TreeStruct {    
    size_t n_h;             // num of hypo
    size_t n_dim;           // total number of dimensions
    int    all_time;        // total CPU time (msec)
    double nll      [HMAX]; // nll per hypo
    int    status   [HMAX]; // status per hypo
    int    strategy [HMAX]; // strategy per hypo
    int    min_time [HMAX]; // minimization time per hypo
    int    dim      [HMAX]; // num of parameters per hypo 
    int    perm     [HMAX]; // num of permutations per hypo
    double param    [PMAX]; // minimzed parameters [ [p_0,...,p_dim[0]], ... ]
    double obs      [PMAX]; // input observables
    double obs_BTAG [PMAX]; // input observables
    void print(ostream& os){
      os << "TreeStruct contains: " << endl;
      os << "\tall_time = " << all_time << " msec" << endl;
      os << "\tn_h      = " << n_h << endl;
      os << "\tn_dim    = " << n_dim << endl;
      size_t it_p {0};
      for( size_t i = 0 ; (i < n_h && i < (size_t)HMAX) ; ++i){
	os << "\tHypo "     << i << ": "   << endl;
	os << "status["     << i << "] = " << status[i];
	os << ", strategy[" << i << "] = " << strategy[i];
	os << ", perm["     << i << "] = " << perm[i];
	os << ", nll["      << i << "] = " << nll[i];
	os << ", min_time[" << i << "] = " << min_time[i];
	os << ", dim["      << i << "] = " << dim[i] << endl;
	for(int j = 0 ; j < dim[i] ; ++j){	  
	  os << "\t\tparam["  << j << "] = " << param[it_p] 
	     << ", obs["      << j << "] = " << obs[it_p]   
	     << ", obs_BTAG[" << j << "] = " << obs_BTAG[it_p] << endl; 
	  ++it_p;
	}
      }
      //for( size_t i = 0 ; (i < n_dim && i < (size_t)PMAX) ; ++i)
      //os << "\tparam[" << i << "] = " << param[i] << endl;
    }
  };

  struct Event {

    Event (TTree*) ;

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
