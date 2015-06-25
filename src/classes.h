#include "TTH/MEIntegratorStandalone/interface/Integrand.h"
#include "TTH/MEIntegratorStandalone/interface/BTagRandomizer.h"

//template class map<string, vector<Algo::Decay::Decay> >;
//template class vector<Algo::Decay::Decay>;
namespace {
    namespace {
        std::vector<MEM::Permutations::Permutations> _i1;
        std::vector<MEM::PSVar::PSVar> _i2;
        std::pair<MEM::TFType::TFType,int> _o1;
        std::pair<MEM::DistributionType::DistributionType,TH3D> _o2;
        std::map<std::pair<MEM::TFType::TFType,int>,TF1> _i3;
        std::map<MEM::DistributionType::DistributionType,TH3D> _i4;
        // MEM::Integrand _i1;
        // MEM::Object _i2;
    }
}
