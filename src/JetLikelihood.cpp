#include "interface/JetLikelihood.h"
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
    double maxCombP = 0;
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
            combP *= _p;
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
    outBestPerm = bestPerm;
    return P;
}
