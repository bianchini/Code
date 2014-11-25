#include <TTree.h>
#include <iostream>

#define HMAX 10
#define PMAX 100


using namespace std;

namespace Algo {

  struct TreeStruct {
    double nll      [HMAX];
    int    status   [HMAX];
    int    min_time [HMAX];
    size_t dim      [HMAX];
    double param    [PMAX];
    double obs      [PMAX];
    double obs_BTAG [PMAX];
    size_t n_h;
    size_t n_dim;
    int    all_time;
    void print(ostream& os){
      os << "TreeStruct contains: " << endl;
      os << "\tall_time = " << all_time << " msec" << endl;
      os << "\tn_h      = " << n_h << endl;
      os << "\tn_dim    = " << n_dim << endl;
      size_t it_p {0};
      for( size_t i = 0 ; (i < n_h && i < (size_t)HMAX) ; ++i){
	os << "\tHypo " << i << ":" << endl;
	os << "\t\tstatus[" << i << "] = " << status[i];
	os << "\t\tnll[" << i << "] = " << nll[i];
	os << "\t\tmin_time[" << i << "] = " << min_time[i];
	os << ", dim[" << i << "] = " << dim[i] << endl;
	for(size_t j = 0 ; j < dim[i] ; ++j){	  
	  os << "\t\tparam[" << j << "] = " 
	     << param[it_p] << "  obs[" << j << "] = " 
	     << obs[it_p] << ", obs_BTAG[" << j << "] = " 
	     << obs_BTAG[it_p] << endl; 
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
