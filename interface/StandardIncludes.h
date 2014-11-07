#include "TLorentzVector.h"
#include<cmath>
#include<cstddef>

using namespace std;

size_t eta_to_bin( const TLorentzVector& lv ){
  if( fabs(lv.Eta())<1.0 ) return 0;
  if( fabs(lv.Eta())>1.0 ) return 1;
  return -99;
}
