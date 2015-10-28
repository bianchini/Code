import ROOT
ROOT.gSystem.Load("libTTHMEIntegratorStandalone")
from ROOT import MEM
jlh = MEM.JetLikelihood()

Cvectoruint = getattr(ROOT, "std::vector<unsigned int>")

for i in range(2):
    jp1 = ROOT.MEM.JetProbability()
    jp1.setProbability(MEM.JetInterpretation.l, 0.2)
    jp1.setProbability(MEM.JetInterpretation.b, 0.8)

    jp2 = ROOT.MEM.JetProbability()
    jp2.setProbability(MEM.JetInterpretation.l, 0.2)
    jp2.setProbability(MEM.JetInterpretation.b, 0.8)

    jlh.push_back_object(jp1)
    jlh.push_back_object(jp2)

    bperm = Cvectoruint()

    print jlh.calcProbability(MEM.JetInterpretation.b, MEM.JetInterpretation.l, 0, bperm)
    print jlh.calcProbability(MEM.JetInterpretation.b, MEM.JetInterpretation.l, 1, bperm)
    print jlh.calcProbability(MEM.JetInterpretation.b, MEM.JetInterpretation.l, 2, bperm)
    print bperm.size(), [bperm.at(i) for i in range(bperm.size())]
    jlh.next_event()
