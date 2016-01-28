#include "TTH/MEIntegratorStandalone/interface/JetLikelihood.h"
#include <iostream>

using namespace std;
using namespace MEM;

double JetLikelihood::calcProbability(
    MEM::JetInterpretation::JetInterpretation hypo1,
    MEM::JetInterpretation::JetInterpretation hypo2,
    unsigned int nhypo1,
    vector<unsigned int>& outBestPerm
    ) {
    vector<unsigned int> perm_index_copy;

    //cout << "calcProbability njets=" << jets.size() << " nhypo1=" << nhypo1 << endl;
    for (unsigned int i=0; i<jets.size(); i++) {
        perm_index_copy.push_back(i);
    }

    int nperms = 0;
    double P = 0.0;

    //use negative value, as any probability will certainly be positive
    double maxCombP = -1.0;
    vector<unsigned int> bestPerm;

    //In case less jets specified than in hypo, truncate
    if (jets.size() < nhypo1) {
        //cout << "truncating " << nhypo1 << " to " << jets.size() << endl;
        nhypo1 = jets.size();
    }

    do {
        unsigned int ip = 0;
        double combP = 1.0;
        for (auto& p : perm_index_copy) {
            const auto& jet = jets.at(p);
            double _p = jet.probas.at(ip < nhypo1 ? hypo1 : hypo2);
            //we don't want a combination ruined just because a probability happened to be 0 
            if (_p > 0) {
                //std::cout << "_p=" << _p << " ip=" << ip << " p=" << p << std::endl;
                combP *= _p;
            //multiply by some small value. FIXME: verify that this does not change results
            } else {
                combP *= 0.00001;
            }
            ip += 1;
        }
        if (combP > maxCombP) {
            maxCombP = combP;
            bestPerm = perm_index_copy;
        }
        P += combP;

        nperms += 1;
    } while( next_combination(
        perm_index_copy.begin(),
        perm_index_copy.begin() + nhypo1,
        perm_index_copy.end(),
        std::less<int>()
        )
    );
    assert(bestPerm.size() == jets.size());
    outBestPerm = bestPerm;
    return P;
}
