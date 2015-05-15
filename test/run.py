import sys
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
verbosity = 1
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
        elif line.startswith("bq") or line.startswith("lq"):
            jet_type, pt, eta, phi, m, btagflag = line.split()
            pt, eta, phi, m, btagflag = tuple(map(float, [pt, eta, phi, m, btagflag]))
            v = TLorentzVector()
            v.SetPtEtaPhiM(pt, eta, phi, m)
            o = MEM.Object(v, MEM.ObjectType.Jet)
            o.addObs(MEM.Observable.BTAG, btagflag)
            mem.push_back_object(o)
            #print "jet", pt, eta, phi, m, btagflag
        elif line.startswith("lp"):
            lep_type, pt, eta, phi, m, charge = line.split()
            pt, eta, phi, m, charge = tuple(map(float, [pt, eta, phi, m, charge]))
            v = TLorentzVector()
            v.SetPtEtaPhiM(pt, eta, phi, m)
            o = MEM.Object(v, MEM.ObjectType.Lepton)
            o.addObs(MEM.Observable.CHARGE, charge)
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

    #print "Calling MEM::run", fstate, hypo
    r = mem.run(fstate, hypo, psvar_vec)
    mem.next_event()
    print r.p
