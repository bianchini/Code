#include <TTree.h>
#include <iostream>

#define HMAX 10
#define PMAX 100


using namespace std;

namespace Algo {

  struct TreeStruct {
    double nll[HMAX];
    size_t n_h;
    size_t n_dim;
    size_t dim  [HMAX];
    double param[PMAX];
    void print(ostream& os){
      os << "TreeStruct contains: " << endl;
      os << "\tn_h   = " << n_h << endl;
      os << "\tn_dim = " << n_dim << endl;
      size_t it_p {0};
      for( size_t i = 0 ; (i < n_h && i < (size_t)HMAX) ; ++i){
	os << "\tHypo " << i << ":" << endl;
	os << "\t\tnll[" << i << "] = " << nll[i];
	os << ", dim[" << i << "] = " << dim[i] << endl;
	for(size_t j = 0 ; j < dim[i] ; ++j){	  
	  os << "\t\tparam[" << j << "] = " << param[it_p] << endl; 
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
    void createBranches () ;
    void reset();

    TreeStruct treeStruct ;

    private :
    TTree* outTree ;

  };

}
