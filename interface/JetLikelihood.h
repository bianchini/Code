#ifndef JETLIKELIHOOD_H
#define JETLIKELIHOOD_H

#include "TTH/MEIntegratorStandalone/interface/Utils.h"
#include <map>

namespace MEM {

    //possible ways a jet can be interepreted
    namespace JetInterpretation {
        enum JetInterpretation {
            l,
            c,
            g,
            b
        };
    }

    //collects probabilities of a jet arising from different interpretations
    class JetProbability {
        public:
            std::map<JetInterpretation::JetInterpretation, double> probas;
            void setProbability(JetInterpretation::JetInterpretation h, double p) {
                probas[h] = p;
            }
    };


    class JetLikelihood {
    public:

        // calculates a probability for a sample of jets being divided between Hypo1 and Hypo2
        // with the assumption that nhypo1 of the jets are from Hypo1
        // e.g. Hypo1 = b, Hypo2 = l, nhypo1=4 calculates the probability that 4 of the jets were b quarks.
        double calcProbability(
            MEM::JetInterpretation::JetInterpretation hypo1,
            MEM::JetInterpretation::JetInterpretation hypo2,
            unsigned int nhypo1,
            std::vector<unsigned int>& outBestPerm
        );
        void push_back_object( const JetProbability& jp ) {
            jets.push_back(jp);
        }
        void next_event() {
            jets.clear();
        }

        std::vector<JetProbability> jets;
    };
}

#endif
