import ROOT, sys

of = ROOT.TFile("a.root", "RECREATE")

tf2 = ROOT.TFile("/scratch/bianchi/74X/btag_out_ttjets_13tev_madgraph_pu20bx25_phys14_0.root")
tt = tf2.Get("tree")
tt.SetBranchStatus("*", False)
tt.SetBranchStatus("njet", True)
tt.SetBranchStatus("pt", True)
tt.SetBranchStatus("csv_inp_*", True)

def hist(x):
    h = (ROOT.TH1D(x, x, 10, 0, 1))
    h.GetXaxis().SetTitle("CSV discriminator")
    return h
    
def hist2d(x):
    h = ROOT.TH2D(x, x, 10,0,1,10,0,1)
    h.GetXaxis().SetTitle("Jet 1")
    h.GetYaxis().SetTitle("Jet 2")
    return h

of.cd()
h_tt = hist2d("corr_tt")
h_bb = hist2d("corr_bb")
h_b  = hist("b")
h_t  = hist("t")
h_t_ttbb  = hist("t_ttbb")
to        = ROOT.TTree("tree")

for iev in range(tt.GetEntries()):
    tt.GetEntry(iev)
    ev = tt

    #if not ( getattr(ev,"pass")[4] and ev.njet==6 ):
    if not ( ev.njet==6 ):
        continue

    for ji in range(ev.njet):
        pt_i, csv_i, pdgid_i, mcid_i = ev.pt[ji], ev.csv_inp_4[ji], ev.pdgid[ji], ev.mcid[ji]

        if abs(pdgid_i)==5 and abs(mcid_i)==0 and pt_i>30 and pt_i<50:
            h_b.Fill(csv_i)
        if abs(pdgid_i)==5 and abs(mcid_i)==6 and pt_i>30 and pt_i<50:
            h_t.Fill(csv_i)
        if abs(pdgid_i)==5 and abs(mcid_i)==6 and pt_i>30 and pt_i<50 and ev.ttCls==53:
            h_t_ttbb.Fill(csv_i)
            

        for jj in range(ev.njet):
            if ji==jj:
                continue
            pt_j, csv_j, pdgid_j, mcid_j = ev.pt[jj], ev.csv_inp_4[jj], ev.pdgid[jj], ev.mcid[jj]

            if abs(pdgid_i)==5 and abs(pdgid_j)==5 and abs(mcid_i)==0 and abs(mcid_j)==0  and ev.ttCls==53:
                h_bb.Fill( csv_i, csv_j)

            if abs(pdgid_i)==5 and abs(pdgid_j)==5 and abs(mcid_i)==6 and abs(mcid_j)==6  and ev.ttCls==53:
                h_tt.Fill( csv_i, csv_j)

of.Write("", ROOT.TObject.kOverwrite)
of.Close()

tf2.Close()
