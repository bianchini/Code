import sys, pickle

#One cool hack to make TFClasses visible for pickle
import TTH.MEAnalysis.TFClasses as TFClasses
sys.modules["TFClasses"] = TFClasses

import ROOT
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gROOT.ProcessLine('AutoLibraryLoader::enable();')
ROOT.gSystem.Load("libFWCoreFWLite")
ROOT.gSystem.Load("libCintex")
ROOT.gROOT.ProcessLine('ROOT::Cintex::Cintex::Enable();')
ROOT.gSystem.Load("libTTHMEIntegratorStandalone")
from ROOT import MEM
from ROOT import TLorentzVector
import math

CvectorPermutations = getattr(ROOT, "std::vector<MEM::Permutations::Permutations>")
CvectorPSVar = getattr(ROOT, "std::vector<MEM::PSVar::PSVar>")

cfg = MEM.MEMConfig()
cfg.defaultCfg()
cfg.transfer_function_method = MEM.TFMethod.Builtin

infiles = sys.argv[1:]
#enum DebugVerbosity { output=1, input=2, init=4, init_more=8, event=16, integration=32};
verbosity = 3+4

try:
    pi_file = open("testdata/transfers.pickle", 'rb')
    tf_matrix = pickle.load(pi_file)
except Exception as e:
    print e
    tf_matrix = None

for nb in [0, 1]:
    for fl1, fl2 in [('b', MEM.TFType.bLost), ('l', MEM.TFType.qLost)]:
        tf = tf_matrix[fl1][nb].Make_CDF()
        #set pt cut for efficiency function
        tf.SetParameter(0, 30)
        tf.SetNpx(10000)
        tf.SetRange(0,400)
        cfg.set_tf_global(fl2, nb, tf)

def add_tf(jet_eta, obj):
    jet_eta_bin = 0
    if abs(jet_eta)>1.0:
        jet_eta_bin = 1
    tf_b = tf_matrix['b'][jet_eta_bin].Make_Formula(False)
    tf_l = tf_matrix['l'][jet_eta_bin].Make_Formula(False)

    tf_b.SetNpx(10000)
    tf_b.SetRange(0,400)

    tf_l.SetNpx(10000)
    tf_l.SetRange(0,400)

    obj.addTransferFunction(MEM.TFType.bReco, tf_b)
    obj.addTransferFunction(MEM.TFType.qReco, tf_l)

#Note: must create integrand __after__ setting global transfer functions
mem = MEM.Integrand(verbosity, cfg)

for inf in infiles:
    print inf
    inf = open(inf, "r").readlines()

    for line in inf:
        line = line.strip()
        if line.startswith("#"):
            continue
        elif line.startswith("fstate"):
            fstate = int(line.split()[1])
        elif line.startswith("integ"):
            integ = int(line.split()[1])
        elif line.startswith("hypo"):
            hypo = int(line.split()[1])
        elif line.startswith("tf_method"):
            tf_method = int(line.split()[1])
            cfg.transfer_function_method = tf_method
        elif line.startswith("int_code"):
            print "cfg.int_code", cfg.int_code
            int_code = int(line.split()[1])
            cfg.int_code = int_code
        elif line.startswith("bq") or line.startswith("lq"):
            jet_type, pt, eta, phi, m, btagflag, match = line.split()
            pt, eta, phi, m, btagflag = tuple(map(float, [pt, eta, phi, m, btagflag]))
            v = TLorentzVector()
            v.SetPtEtaPhiM(pt, eta, phi, m)
            o = MEM.Object(v, MEM.ObjectType.Jet)
            o.addObs(MEM.Observable.BTAG, btagflag)
            add_tf(eta, o)
            mem.push_back_object(o)
            #print "jet", pt, eta, phi, m, btagflag
        elif line.startswith("lp"):
            lep_type, pt, eta, phi, m, charge = line.split()
            pt, eta, phi, m, charge = tuple(map(float, [pt, eta, phi, m, charge]))
            v = TLorentzVector()
            v.SetPtEtaPhiM(pt, eta, phi, m)
            o = MEM.Object(v, MEM.ObjectType.Lepton)
            o.addObs(MEM.Observable.CHARGE, charge)
            add_tf(eta, o)
            mem.push_back_object(o)
            #print "lep", pt, eta, phi
        elif line.startswith("mt"):
            met_type, pt, phi = line.split()
            pt, phi = tuple(map(float, [pt, phi]))
            v = TLorentzVector()
            v.SetPtEtaPhiM(pt, 0.0, phi, 0.0)
            o = MEM.Object(v, MEM.ObjectType.MET)
            mem.push_back_object(o)
            #print "met", pt, eta, phi

    pvec = CvectorPermutations()
    pvec.push_back(MEM.Permutations.BTagged)
    pvec.push_back(MEM.Permutations.QUntagged)
    mem.set_permutation_strategy(pvec)

    psvar_vec = CvectorPSVar()
    if integ == 1:
        psvar_vec.push_back(MEM.PSVar.cos_qbar1)
        psvar_vec.push_back(MEM.PSVar.phi_qbar1)
    mem.set_cfg(cfg)

    print "Calling MEM::run", fstate, hypo
    r = mem.run(fstate, hypo, psvar_vec)
    mem.next_event()
    print r.p
