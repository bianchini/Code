import ROOT
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gROOT.ProcessLine('AutoLibraryLoader::enable();')
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gSystem.Load("libCintex")
ROOT.gROOT.ProcessLine('ROOT::Cintex::Cintex::Enable();')
ROOT.gSystem.Load("libTTHMEIntegratorStandalone")
from ROOT import Algo
from ROOT import TLorentzVector

from ROOT.std import map
tf = ROOT.TFile("out.root", "RECREATE")
tf.cd()
tt = ROOT.TTree("a", "a")

j1 = TLorentzVector()
j1.SetPtEtaPhiM(50., 0., ROOT.TMath.Pi()*2/3, 0.);

j2 = TLorentzVector()
j2.SetPtEtaPhiM( 30., 0.,  0., 0.);

j3 = TLorentzVector()
j3.SetPtEtaPhiM( 50., 0., -ROOT.TMath.Pi()/2, 0.);

j4 = TLorentzVector()
j4.SetPtEtaPhiM(100., 0., -ROOT.TMath.Pi()/4, 0.);

j5 = TLorentzVector()
j5.SetPtEtaPhiM(80., 0., -ROOT.TMath.Pi()/5, 0.);

j6 = TLorentzVector()
j6.SetPtEtaPhiM(120., 0., -ROOT.TMath.Pi()/6, 0.);

lep = TLorentzVector()
lep.SetPtEtaPhiM( 50., 1., ROOT.TMath.Pi()/3, 0.);

met = TLorentzVector()
met.SetPtEtaPhiM( 0., 0., 0., 0.);

ht = Algo.HypoTester(tt)
#ht = Algo.HypoTester()
for j in [j1, j2, j3, j4, j5, j6]:
    ht.push_back_object(j, 'j')
ht.push_back_object(lep, 'l')
ht.push_back_object(met, 'm')

hmap = getattr(ROOT, "map<string,vector<Algo::Decay::Decay> >")()
#ROOT.gROOT.ProcessLine("typedef std::map<std::string, std::vector<Algo::Decay> > HypoMap;")
#hmap = ROOT.HypoMap()
vec0 = getattr(ROOT, "vector<Algo::Decay::Decay>")()
vec0.push_back(ROOT.Algo.Decay.TopHad)

hmap[ROOT.std.string("h0")] = vec0

print "Testing hypos"
ht.test(hmap)

tf.Write()
print "All done"
