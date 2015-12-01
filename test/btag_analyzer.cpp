#include "interface/BTagRandomizer.h"
#include "TFile.h"
#include "TH3D.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include<iostream>

using namespace std;
using namespace MEM;


//
// btag_analyzer SAMPLE ASSIGN_RND LABEL EVT_MIN EVT_MAX OPTION
//  > OPTION: 0 = b,c,l,pu
//  > OPTION: 1 = b/1b/2b,c,l,pu
//  > OPTION: 2 = b/1b/2b,c/1c/2c,l,pu
//  > OPTION: 3 = b/1b/2b,c/1c/2c,ud,g,s,pu
//  > OPTION: 4 = b/1b/2b,c/1c/2c, l/dR, pu
//  > OPTION: 5 = b/dR, 1b/2b, c/dR, 1c/2c, l/dR, pu
//  > OPTION: 6 = b[1b/2b], c[1c/2c], l, pu 
//

int main(int argc, char *argv[]){

  if(argc<5) return 0;

  int option = atoi(argv[6]);

  double csv_cut{0.890};
  
  //TFile* fout  = new TFile( Form("/scratch/bianchi/vhbb74X/btag_out_%s_%d%s_%d_%d.root", 
  TFile* fout  = new TFile( Form("./btag_out_%s_rnd%d_csv%s_files%d_%d_option%d.root", 
				 argv[1], atoi(argv[2]), argc>3 ? argv[3] : "", atoi(argv[4]), atoi(argv[5]), option),"RECREATE");
  TTree* tout  = new TTree("tree", "");
  int njet, ntag, ncat, event, nB, nC, nL, ttCls;
  float pcat        [5];
  int   pass        [5];
  int   pass_rnd    [5];
  float csv_rnd_0   [8];
  float csv_inp_0   [8];
  float csv_rnd_1   [8];
  float csv_inp_1   [8];
  float csv_rnd_2   [8];
  float csv_inp_2   [8];
  float csv_rnd_3   [8];
  float csv_inp_3   [8];
  float csv_rnd_4   [8];
  float csv_inp_4   [8];

  float corr_rnd_0   [8];
  float corr_inp_0   [8];
  float corr_rnd_1   [8];
  float corr_inp_1   [8];
  float corr_rnd_2   [8];
  float corr_inp_2   [8];
  float corr_rnd_3   [8];
  float corr_inp_3   [8];
  float corr_rnd_4   [8];
  float corr_inp_4   [8];

  float HT, met;
  float pt    [8];
  float mcpt  [8];
  float eta   [8];
  int   pdgid [8];

  map<int, float*> csv_rnd_map;
  map<int, float*> csv_inp_map;
  csv_inp_map[0] = csv_inp_0;
  csv_rnd_map[0] = csv_rnd_0;
  csv_inp_map[1] = csv_inp_1;
  csv_rnd_map[1] = csv_rnd_1;
  csv_inp_map[2] = csv_inp_2;
  csv_rnd_map[2] = csv_rnd_2;
  csv_inp_map[3] = csv_inp_3;
  csv_rnd_map[3] = csv_rnd_3;
  csv_inp_map[4] = csv_inp_4;
  csv_rnd_map[4] = csv_rnd_4;

  map<int, float*> corr_rnd_map;
  map<int, float*> corr_inp_map;
  corr_inp_map[0] = corr_inp_0;
  corr_rnd_map[0] = corr_rnd_0;
  corr_inp_map[1] = corr_inp_1;
  corr_rnd_map[1] = corr_rnd_1;
  corr_inp_map[2] = corr_inp_2;
  corr_rnd_map[2] = corr_rnd_2;
  corr_inp_map[3] = corr_inp_3;
  corr_rnd_map[3] = corr_rnd_3;
  corr_inp_map[4] = corr_inp_4;
  corr_rnd_map[4] = corr_rnd_4;


  //tout->Branch("event",    &event,    "event/I");
  tout->Branch("njet",     &njet,     "njet/I");
  //tout->Branch("nB",       &nB,       "nB/I");
  //tout->Branch("nC",       &nC,       "nC/I");
  //tout->Branch("nL",       &nL,       "nL/I");
  tout->Branch("ttCls",    &ttCls,    "ttCls/I");
  tout->Branch("ncat",     &ncat,     "ncat/I");
  tout->Branch("ntag",     &ntag,     "ntag/I");
  tout->Branch("pcat",     pcat,      "pcat[ncat]/F");
  tout->Branch("pass",     pass,      "pass[ncat]/I");
  tout->Branch("pass_rnd", pass_rnd,  "pass_rnd[ncat]/I");
  //tout->Branch("pt",       pt,        "pt[njet]/F");
  //tout->Branch("mcpt",     mcpt,       "mcpt[njet]/F");
  //tout->Branch("eta",      eta,       "eta[njet]/F");
  //tout->Branch("pdgid",    pdgid,     "pdgid[njet]/I");
  //tout->Branch("HT",       &HT,       "HT/F");
  //tout->Branch("met",      &met,      "met/F");

  //for(auto it = csv_rnd_map.begin() ; it!=csv_rnd_map.end(); ++it)
  //  tout->Branch(Form("csv_rnd_%d",it->first),  it->second,   Form("csv_rnd_%d[njet]/F",it->first));

  for(auto it = csv_inp_map.begin() ; it!=csv_inp_map.end(); ++it)
    tout->Branch(Form("csv_inp_%d",it->first),  it->second,   Form("csv_inp_%d[njet]/F",it->first));

  //for(auto it = corr_rnd_map.begin() ; it!=corr_rnd_map.end(); ++it)
  //  tout->Branch(Form("corr_rnd_%d",it->first),  it->second,   Form("corr_rnd_%d[8]/F",it->first));

  //for(auto it = corr_inp_map.begin() ; it!=corr_inp_map.end(); ++it)
  //  tout->Branch(Form("corr_inp_%d",it->first),  it->second,   Form("corr_inp_%d[8]/F",it->first));


  TFile* f = TFile::Open("./csv"+TString(argv[3])+".root");
  TH3D* h3_b  = (TH3D*)f->Get("csv_b_pt_eta");
  TH3D* h3_c  = (TH3D*)f->Get("csv_c_pt_eta");
  TH3D* h3_l  = (TH3D*)f->Get("csv_l_pt_eta");
  TH3D* h3_g  = (TH3D*)f->Get("csv_g_pt_eta");
  TH3D* h3_s  = (TH3D*)f->Get("csv_s_pt_eta");
  TH3D* h3_ud = (TH3D*)f->Get("csv_ud_pt_eta");
  TH3D* h3_pu = (TH3D*)f->Get("csv_pu_pt_eta");
  TH3D* h3_1b = (TH3D*)f->Get("csv_1b_pt_eta");
  TH3D* h3_1c = (TH3D*)f->Get("csv_1c_pt_eta");
  TH3D* h3_2b = (TH3D*)f->Get("csv_2b_pt_eta");
  TH3D* h3_2c = (TH3D*)f->Get("csv_2c_pt_eta");
  TH3D* h3_b_dR_0p0_1p0 = (TH3D*)f->Get("csv_b_dR_0p0_1p0_pt_eta");
  TH3D* h3_b_dR_1p0_2p0 = (TH3D*)f->Get("csv_b_dR_1p0_2p0_pt_eta");
  TH3D* h3_b_dR_2p0_Inf = (TH3D*)f->Get("csv_b_dR_2p0_Inf_pt_eta");
  TH3D* h3_1b_dR_0p0_1p0 = (TH3D*)f->Get("csv_1b_dR_0p0_1p0_pt_eta");
  TH3D* h3_1b_dR_1p0_2p0 = (TH3D*)f->Get("csv_1b_dR_1p0_2p0_pt_eta");
  TH3D* h3_1b_dR_2p0_Inf = (TH3D*)f->Get("csv_1b_dR_2p0_Inf_pt_eta");
  TH3D* h3_2b_dR_0p0_1p0 = (TH3D*)f->Get("csv_2b_dR_0p0_1p0_pt_eta");
  TH3D* h3_2b_dR_1p0_2p0 = (TH3D*)f->Get("csv_2b_dR_1p0_2p0_pt_eta");
  TH3D* h3_2b_dR_2p0_Inf = (TH3D*)f->Get("csv_2b_dR_2p0_Inf_pt_eta");
  TH3D* h3_1c_dR_0p0_1p0 = (TH3D*)f->Get("csv_1c_dR_0p0_1p0_pt_eta");
  TH3D* h3_1c_dR_1p0_2p0 = (TH3D*)f->Get("csv_1c_dR_1p0_2p0_pt_eta");
  TH3D* h3_1c_dR_2p0_Inf = (TH3D*)f->Get("csv_1c_dR_2p0_Inf_pt_eta");
  TH3D* h3_2c_dR_0p0_1p0 = (TH3D*)f->Get("csv_2c_dR_0p0_1p0_pt_eta");
  TH3D* h3_2c_dR_1p0_2p0 = (TH3D*)f->Get("csv_2c_dR_1p0_2p0_pt_eta");
  TH3D* h3_2c_dR_2p0_Inf = (TH3D*)f->Get("csv_2c_dR_2p0_Inf_pt_eta");
  TH3D* h3_c_dR_0p0_1p0 = (TH3D*)f->Get("csv_c_dR_0p0_1p0_pt_eta");
  TH3D* h3_c_dR_1p0_2p0 = (TH3D*)f->Get("csv_c_dR_1p0_2p0_pt_eta");
  TH3D* h3_c_dR_2p0_Inf = (TH3D*)f->Get("csv_c_dR_2p0_Inf_pt_eta");
  TH3D* h3_pu_dR_0p0_1p0 = (TH3D*)f->Get("csv_pu_dR_0p0_1p0_pt_eta");
  TH3D* h3_pu_dR_1p0_2p0 = (TH3D*)f->Get("csv_pu_dR_1p0_2p0_pt_eta");
  TH3D* h3_pu_dR_2p0_Inf = (TH3D*)f->Get("csv_pu_dR_2p0_Inf_pt_eta");
  TH3D* h3_s_dR_0p0_1p0 = (TH3D*)f->Get("csv_s_dR_0p0_1p0_pt_eta");
  TH3D* h3_s_dR_1p0_2p0 = (TH3D*)f->Get("csv_s_dR_1p0_2p0_pt_eta");
  TH3D* h3_s_dR_2p0_Inf = (TH3D*)f->Get("csv_s_dR_2p0_Inf_pt_eta");
  TH3D* h3_ud_dR_0p0_1p0 = (TH3D*)f->Get("csv_ud_dR_0p0_1p0_pt_eta");
  TH3D* h3_ud_dR_1p0_2p0 = (TH3D*)f->Get("csv_ud_dR_1p0_2p0_pt_eta");
  TH3D* h3_ud_dR_2p0_Inf = (TH3D*)f->Get("csv_ud_dR_2p0_Inf_pt_eta");
  TH3D* h3_g_dR_0p0_1p0 = (TH3D*)f->Get("csv_g_dR_0p0_1p0_pt_eta");
  TH3D* h3_g_dR_1p0_2p0 = (TH3D*)f->Get("csv_g_dR_1p0_2p0_pt_eta");
  TH3D* h3_g_dR_2p0_Inf = (TH3D*)f->Get("csv_g_dR_2p0_Inf_pt_eta");
  TH3D* h3_l_dR_0p0_1p0 = (TH3D*)f->Get("csv_l_dR_0p0_1p0_pt_eta");
  TH3D* h3_l_dR_1p0_2p0 = (TH3D*)f->Get("csv_l_dR_1p0_2p0_pt_eta");
  TH3D* h3_l_dR_2p0_Inf = (TH3D*)f->Get("csv_l_dR_2p0_Inf_pt_eta");


  std::map<MEM::DistributionType::DistributionType, TH3D> btag_pdfs;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_b]  = *h3_b ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_c]  = *h3_c ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_l]  = *h3_l ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_ud] = *h3_ud ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_g]  = *h3_g ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_s]  = *h3_s ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_pu] = *h3_pu ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1b] = *h3_1b ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2b] = *h3_2b ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1c] = *h3_1c ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2c] = *h3_2c ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_b_dR_0p0_1p0]  = *h3_b_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_b_dR_1p0_2p0]  = *h3_b_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_b_dR_2p0_Inf]  = *h3_b_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1b_dR_0p0_1p0] = *h3_1b_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1b_dR_1p0_2p0] = *h3_1b_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1b_dR_2p0_Inf] = *h3_1b_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2b_dR_0p0_1p0] = *h3_2b_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2b_dR_1p0_2p0] = *h3_2b_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2b_dR_2p0_Inf] = *h3_2b_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_c_dR_0p0_1p0]  = *h3_c_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_c_dR_1p0_2p0]  = *h3_c_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_c_dR_2p0_Inf]  = *h3_c_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1c_dR_0p0_1p0] = *h3_1c_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1c_dR_1p0_2p0] = *h3_1c_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_1c_dR_2p0_Inf] = *h3_1c_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2c_dR_0p0_1p0] = *h3_2c_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2c_dR_1p0_2p0] = *h3_2c_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_2c_dR_2p0_Inf] = *h3_2c_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_pu_dR_0p0_1p0]  = *h3_pu_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_pu_dR_1p0_2p0]  = *h3_pu_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_pu_dR_2p0_Inf]  = *h3_pu_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_s_dR_0p0_1p0]  = *h3_s_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_s_dR_1p0_2p0]  = *h3_s_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_s_dR_2p0_Inf]  = *h3_s_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_ud_dR_0p0_1p0]  = *h3_ud_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_ud_dR_1p0_2p0]  = *h3_ud_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_ud_dR_2p0_Inf]  = *h3_ud_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_g_dR_0p0_1p0]  = *h3_g_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_g_dR_1p0_2p0]  = *h3_g_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_g_dR_2p0_Inf]  = *h3_g_dR_2p0_Inf ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_l_dR_0p0_1p0]  = *h3_l_dR_0p0_1p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_l_dR_1p0_2p0]  = *h3_l_dR_1p0_2p0 ;
  btag_pdfs[MEM::DistributionType::DistributionType::csv_l_dR_2p0_Inf]  = *h3_l_dR_2p0_Inf ;

  BTagRandomizer* rnd = new BTagRandomizer(DebugVerbosity::init
					   //|DebugVerbosity::init_more
					   //|DebugVerbosity::event
					   , -1            // seed
   					   , btag_pdfs     // pdfs
					   , atoi(argv[2]) // assign random csv
					   , 0             // compress csv [0,1]
					   , 5000          // max number of toys
					   );

  
  //TCut cut = "Vtype>=0";

  TString path = "dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/t3groups/ethz-higgs/run2/VHBBHeppyV14/";
  vector<TString> samples;
  //samples.push_back("TT_TuneCUETP8M1_13TeV-powheg-pythia8/VHBB_HEPPY_V14_TT_TuneCUETP8M1_13TeV-powheg-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_212301/0000/");
  //samples.push_back("WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V14_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_220559/0000/");
  samples.push_back("TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/VHBB_HEPPY_V14_TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_212453/0000/");
  //samples.push_back("ttHTobb_M125_13TeV_powheg_pythia8/VHBB_HEPPY_V14_ttHTobb_M125_13TeV_powheg_pythia8__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151025_084144/0000/");
  TChain* ch = new TChain("tree");
  for(unsigned int ns = 0; ns<samples.size(); ++ns){
    TString sample = samples[ns];
    cout << sample << endl;
    for(int nf=atoi(argv[4]); nf <=atoi(argv[5]); ++nf)
      ch->AddFile(path+sample+"/tree_"+TString(Form("%d",nf))+".root");
  }

  //TTreeFormula* treeformula = new TTreeFormula("selection", cut , ch );

  int nJet;
  ULong64_t evt;
  float ttCls_;
  float Jet_pt     [15];
  float Jet_eta    [15];
  float Jet_mass   [15];
  float Jet_phi    [15];
  float Jet_btagCSV[15];
  float Jet_mcPt   [15];
  float Jet_mcEta  [15];
  int Jet_mcFlavour [15];
  int Jet_hadronFlavour [15];
  int Jet_mcIdx[15];
  int GenJet_numBHadrons[15];
  int GenJet_numBHadronsFromTop[15];
  int GenJet_numBHadronsAfterTop[15];
  int GenJet_numCHadrons[15];
  int GenJet_numCHadronsFromTop[15];
  int GenJet_numCHadronsAfterTop[15];

  float met_pt;

  ch->SetBranchAddress("evt",         &evt);
  ch->SetBranchAddress("ttCls",       &ttCls_);
  ch->SetBranchAddress("nJet",        &nJet);
  ch->SetBranchAddress("Jet_pt",      Jet_pt );
  ch->SetBranchAddress("Jet_eta",     Jet_eta  );
  ch->SetBranchAddress("Jet_mass",    Jet_mass  );
  ch->SetBranchAddress("Jet_phi",     Jet_phi  );
  ch->SetBranchAddress("Jet_btagCSV", Jet_btagCSV);
  ch->SetBranchAddress("Jet_mcPt",    Jet_mcPt   );
  ch->SetBranchAddress("Jet_mcEta",   Jet_mcEta   );
  ch->SetBranchAddress("Jet_mcFlavour",      Jet_mcFlavour);
  ch->SetBranchAddress("Jet_hadronFlavour",  Jet_hadronFlavour);
  ch->SetBranchAddress("Jet_mcIdx",          Jet_mcIdx);
  ch->SetBranchAddress("GenJet_numBHadrons",              GenJet_numBHadrons);
  ch->SetBranchAddress("GenJet_numBHadronsAfterTop",      GenJet_numBHadronsAfterTop);
  ch->SetBranchAddress("GenJet_numBHadronsFromTop",       GenJet_numBHadronsFromTop);
  ch->SetBranchAddress("GenJet_numCHadrons",              GenJet_numCHadrons);
  ch->SetBranchAddress("GenJet_numCHadronsAfterTop",      GenJet_numCHadronsAfterTop);
  ch->SetBranchAddress("GenJet_numCHadronsFromTop",       GenJet_numCHadronsFromTop);
  ch->SetBranchAddress("met_pt",   &met_pt);

  Long64_t nentries = ch->GetEntries();
  cout << "Total entries: " << nentries << endl;
  int count_pass{0};

  for (Long64_t i = 0; i < nentries; i++){

    ch->GetEntry(i);

    //ch->LoadTree(i);
    //if( treeformula->EvalInstance() == 0){
    //  continue;
    //}
    ++count_pass;

    //if( count_pass < atoi(argv[4]))  continue;
    //if( count_pass > 1000 )  break;


    if(i%500==0) cout << "Event " << i << " (" << evt << ")" << endl;

    vector<JetCategory> cats;
    vector<Object*> objects_mem;
    vector<JetCategory*> cats_mem;

    HT  = 0.;
    met = met_pt; 

    nB = 0;
    nC = 0;
    nL = 0;
    int nBCSVM = 0;

    vector<int> good_jets;
    for( int j = 0 ; j < nJet ; ++j){
      if( Jet_pt[j]<30. || std::abs(Jet_eta[j])>2.5 ) continue;
      good_jets.push_back(j);
    }

    int goodjets = good_jets.size();
    if(goodjets>8 || goodjets<4) continue;

    int count_j = 0;
    for( int gj = 0 ; gj < goodjets ; ++gj){

      int j = good_jets[gj];

      bool marg = false;

      pt[j]    = Jet_pt[j];
      mcpt[j]  = Jet_mcPt[j];
      HT      += Jet_pt[j];
      eta[j]   = Jet_eta[j];
      pdgid[j] = Jet_hadronFlavour[j];

      float min_dR = 999.;
      for( int gk = 0; gk < goodjets ; ++gk){
	int k = good_jets[gk];
	if(k==j) continue;
	float dR = TMath::Sqrt( TMath::Power(Jet_eta[j]-Jet_eta[k], 2) +  TMath::Power(Jet_phi[j]-Jet_phi[k], 2)  );
	if( dR < min_dR ) min_dR = dR;
      }
      //cout << "******************" << endl;
      //cout << "min dR: " << min_dR << endl;

      float csv = Jet_btagCSV[j];
      if(csv<0.) csv = -0.005;
      if(csv>1.) csv = +1.005;


      if ( csv>csv_cut) ++nBCSVM;

      if     ( std::abs(Jet_hadronFlavour[j])==5 ) ++nB;
      else if( std::abs(Jet_hadronFlavour[j])==4 ) ++nC;
      else ++nL;
	
      TLorentzVector lv;
      lv.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);

      int mcIdx = Jet_mcPt[j]>20. && Jet_mcIdx[j]>=0 && std::abs(Jet_mcEta[j])<2.4 ? Jet_mcIdx[j] : -1;

      int numBHadrons =  mcIdx>=0 ? GenJet_numBHadrons[mcIdx] : 0;// + GenJet_numBHadronsAfterTop[mcIdx] + GenJet_numBHadronsFromTop[mcIdx] : 0;
      int numCHadrons =  mcIdx>=0 ? GenJet_numCHadrons[mcIdx] : 0;// + GenJet_numCHadronsAfterTop[mcIdx] + GenJet_numCHadronsFromTop[mcIdx] : 0;

      MEM::DistributionType::DistributionType dt = MEM::DistributionType::DistributionType::csv_l;
      MEM::DistributionType::DistributionType dt_main = dt;
      if(option==4 || option==5){      
	if( min_dR < 1.0 )
	  dt = MEM::DistributionType::DistributionType::csv_l_dR_0p0_1p0;
	else if( min_dR < 2.0 )
	  dt = MEM::DistributionType::DistributionType::csv_l_dR_1p0_2p0;
	else
	  dt = MEM::DistributionType::DistributionType::csv_l_dR_2p0_Inf;
      }

      if( Jet_mcPt[j]<5. ){
	dt =  MEM::DistributionType::DistributionType::csv_pu;
	dt_main = dt;
	/*
	if( min_dR < 1.0 )
	  dt = MEM::DistributionType::DistributionType::csv_pu_dR_0p0_1p0;
	else if( min_dR < 2.0 )
	  dt = MEM::DistributionType::DistributionType::csv_pu_dR_1p0_2p0;
	else
	  dt = MEM::DistributionType::DistributionType::csv_pu_dR_2p0_Inf;
	*/
      }
      else{
	int hadflav = Jet_hadronFlavour[j];
	switch(hadflav){
	case 5:
	  dt =  MEM::DistributionType::DistributionType::csv_b;
	  dt_main = dt;	
	  if(option==5){
	    if( min_dR < 1.0 )
	      dt = MEM::DistributionType::DistributionType::csv_b_dR_0p0_1p0;
	    else if( min_dR < 2.0 )
	      dt = MEM::DistributionType::DistributionType::csv_b_dR_1p0_2p0;
	    else
	      dt = MEM::DistributionType::DistributionType::csv_b_dR_2p0_Inf;
	  }
	  if( mcIdx>=0 && numBHadrons == 1 ){
	    if(option>=1)
	      dt =  MEM::DistributionType::DistributionType::csv_1b;
	    if(option==6) marg = true;
	    /*
	    if( min_dR < 1.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_1b_dR_0p0_1p0;
	    else if( min_dR < 2.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_1b_dR_1p0_2p0;
	    else
	      dt = MEM::DistributionType::DistributionType::csv_1b_dR_2p0_Inf;
	    */
	  }
	  else if( mcIdx>=0 && numBHadrons >= 2 ){
	    if(option>=1)
	      dt =  MEM::DistributionType::DistributionType::csv_2b;
	    if(option==6) marg = true;
	    /*
	    if( min_dR < 1.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_2b_dR_0p0_1p0;
	    else if( min_dR < 2.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_2b_dR_1p0_2p0;
	    else
	      dt = MEM::DistributionType::DistributionType::csv_2b_dR_2p0_Inf;	    
	    */
	  }
	  break;
	case 4:
	  dt =  MEM::DistributionType::DistributionType::csv_c;
	  dt_main = dt;
	  if(option==5){  
	    if( min_dR < 1.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_c_dR_0p0_1p0;
	    else if( min_dR < 2.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_c_dR_1p0_2p0;
	    else
	      dt = MEM::DistributionType::DistributionType::csv_c_dR_2p0_Inf;	    
	  }

	  if( mcIdx>=0 && numCHadrons==1 ){
	    if(option>=2)
	      dt =  MEM::DistributionType::DistributionType::csv_1c;
	    if(option==6) marg = true;
	    /*
	    if( min_dR < 1.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_1c_dR_0p0_1p0;
	    else if( min_dR < 2.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_1c_dR_1p0_2p0;
	    else
	      dt = MEM::DistributionType::DistributionType::csv_1c_dR_2p0_Inf;
	    */
	  }
	  if( mcIdx>=0 && numCHadrons>=2 ){
	    if(option>=2)
	      dt =  MEM::DistributionType::DistributionType::csv_2c;
	    if(option==6) marg = true;
	    /*
	    if( min_dR < 1.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_2c_dR_0p0_1p0;
	    else if( min_dR < 2.0 ) 
	      dt = MEM::DistributionType::DistributionType::csv_2c_dR_1p0_2p0;
	    else
	      dt = MEM::DistributionType::DistributionType::csv_2c_dR_2p0_Inf;	    
	    */
	  }
	  break;
	default:
	  if(option==3){
	    if     (  std::abs( Jet_mcFlavour[j] )==3 ){             
	      dt =  MEM::DistributionType::DistributionType::csv_s;
	      /*
		if( min_dR < 1.0 )
		dt = MEM::DistributionType::DistributionType::csv_s_dR_0p0_1p0;
		else if( min_dR < 2.0 )
		dt = MEM::DistributionType::DistributionType::csv_s_dR_1p0_2p0;
		else
		dt = MEM::DistributionType::DistributionType::csv_s_dR_2p0_Inf;
	      */
	    }
	    else if(  std::abs( Jet_mcFlavour[j] )==2 ||  std::abs( Jet_mcFlavour[j] )==1 ){
	      dt =  MEM::DistributionType::DistributionType::csv_ud;
	      /*
		if( min_dR < 1.0 )
		dt = MEM::DistributionType::DistributionType::csv_ud_dR_0p0_1p0;
		else if( min_dR < 2.0 )
		dt = MEM::DistributionType::DistributionType::csv_ud_dR_1p0_2p0;
		else
		dt = MEM::DistributionType::DistributionType::csv_ud_dR_2p0_Inf;
	      */
	    }
	    else if(  std::abs( Jet_mcFlavour[j] )==21 ){            
	      dt =  MEM::DistributionType::DistributionType::csv_g;
	      /*
		if( min_dR < 1.0 )
		dt = MEM::DistributionType::DistributionType::csv_g_dR_0p0_1p0;
		else if( min_dR < 2.0 )
		dt = MEM::DistributionType::DistributionType::csv_g_dR_1p0_2p0;
		else
		dt = MEM::DistributionType::DistributionType::csv_g_dR_2p0_Inf;
	      */
	    }
	  }
	  break;
	}
      }
      
      //cout << static_cast<int>(dt) << ": mcFlavour: " << Jet_mcFlavour[j]  << ", hadronFlavour: " << Jet_hadronFlavour[j] << endl;
      Object* jet = new Object( lv, ObjectType::Jet, dt, dt_main );
      jet->addObs( Observable::BTAG, csv>csv_cut   );  
      jet->addObs( Observable::CSV, csv );  

      if( marg ){
	jet->addObs( Observable::IGNORE_FOR_RND,  1.0);
      }
      rnd->push_back_object( jet );

      if( count_j<=4 ){
	int ntags_l = count_j; 
	int ntags_h = ntags_l<4 ? count_j : -1;
	JetCategory* cat = new JetCategory(ntags_l, ntags_h, csv_cut, goodjets, 
					   string(Form("%dj%s%dt", goodjets, ntags_h>=0 ? "" : "ge" , ntags_l)));
	cats.push_back( *cat );
	cats_mem.push_back( cat );
      }
      
      objects_mem.push_back( jet );
      ++count_j;
    }

    
    if( goodjets==4 ){
      int ntags_l = 4; 
      int ntags_h = -1;
      JetCategory* cat = new JetCategory(ntags_l, ntags_h, csv_cut, 4, 
					 string(Form("%dj%s%dt", goodjets, ntags_h>=0 ? "" : "ge" , ntags_l)));
      cats.push_back( *cat );
      cats_mem.push_back( cat );
    }
    

    event = evt;
    njet  = goodjets;
    ncat  = int(cats.size());    
    ntag  = int(nBCSVM);
    ttCls = ttCls_;    

    //cout << "\tRunning BTagRandomizerOutput with " << cats.size() << " confs" << endl;
    vector<BTagRandomizerOutput> out_all = rnd->run_all( cats );
    for(int o = 0; o < int(out_all.size()) ; ++o ){
      pcat[o]     = out_all[o].p;
      pass[o]     = out_all[o].pass;

      pass_rnd[o] = out_all[o].pass_rnd;
      //out_all[o].print(cout);    

      for( int j = 0 ; j < nJet ; ++j){
	(csv_rnd_map.at(o))[j] = (out_all[o].rnd_btag)[j];
	(csv_inp_map.at(o))[j] = (out_all[o].input_btag)[j];
      }

      //// Input
      float dR_min_bb =  9999.;
      float dR_min_jj =  9999.;
      float dR_max_bb = -99.;
      float dR_max_jj = -99.;
      float m_min_bb  =  9999.;
      float m_min_jj  =  9999.;
      float m_max_bb  = -99.;
      float m_max_jj  = -99.;
      for( int gj = 0 ; gj < goodjets-1 ; ++gj){
	int j = good_jets[gj];
	TLorentzVector lv_j;
	lv_j.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	int tag_j= (out_all[o].input_btag)[j]>csv_cut;
	for( int gk = gj+1 ; gk < goodjets ; ++gk){
	  int k = good_jets[gk];
	  TLorentzVector lv_k;
	  lv_k.SetPtEtaPhiM(Jet_pt[k], Jet_eta[k], Jet_phi[k], Jet_mass[k]);	
	  int tag_k = (out_all[o].input_btag)[k]>csv_cut;
	  
	  float dR = lv_j.DeltaR( lv_k );
	  float m  = (lv_k+lv_j).M();
	  
	  if( tag_j && tag_k ){
	    if( dR <= dR_min_bb ) dR_min_bb = dR;
	    if( m  <= m_min_bb  ) m_min_bb  = m;
	    if( dR >= dR_max_bb ) dR_max_bb = dR;
	    if( m  >= m_max_bb  ) m_max_bb  = m;
	  }
	  if( !tag_j && !tag_k ){
	    if( dR <= dR_min_jj ) dR_min_jj = dR;
	    if( m  <= m_min_jj  ) m_min_jj  = m;
	    if( dR >= dR_max_jj ) dR_max_jj = dR;
	    if( m  >= m_max_jj  ) m_max_jj  = m;
	  }
	  
	}
      }
      
      (corr_inp_map.at(o))[0] = dR_min_bb;
      (corr_inp_map.at(o))[1] = dR_min_jj;
      (corr_inp_map.at(o))[2] = dR_max_bb;
      (corr_inp_map.at(o))[3] = dR_max_jj;
      (corr_inp_map.at(o))[4] = m_min_bb;
      (corr_inp_map.at(o))[5] = m_min_jj;
      (corr_inp_map.at(o))[6] = m_max_bb;
      (corr_inp_map.at(o))[7] = m_max_jj;

      ///////// Rand
      dR_min_bb =  9999.;
      dR_min_jj =  9999.;
      dR_max_bb = -99.;
      dR_max_jj = -99.;
      m_min_bb  =  9999.;
      m_min_jj  =  9999.;
      m_max_bb  = -99.;
      m_max_jj  = -99.;

      for( int gj = 0 ; gj < goodjets-1 ; ++gj){
	int j = good_jets[gj];
	TLorentzVector lv_j;
	lv_j.SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	int tag_j= (out_all[o].input_btag)[j]>csv_cut;
	for( int gk = gj+1 ; gk < goodjets ; ++gk){
	  int k = good_jets[gk];
	  TLorentzVector lv_k;
	  lv_k.SetPtEtaPhiM(Jet_pt[k], Jet_eta[k], Jet_phi[k], Jet_mass[k]);	
	  
	  int tag_k = (out_all[o].rnd_btag)[k]>csv_cut;	  
	  float dR = lv_j.DeltaR( lv_k );
	  float m  = (lv_k+lv_j).M();
	  
	  if( tag_j && tag_k ){
	    if( dR <= dR_min_bb ) dR_min_bb = dR;
	    if( m  <= m_min_bb  ) m_min_bb  = m;
	    if( dR >= dR_max_bb ) dR_max_bb = dR;
	    if( m  >= m_max_bb  ) m_max_bb  = m;
	  }
	  if( !tag_j && !tag_k ){
	    if( dR <= dR_min_jj ) dR_min_jj = dR;
	    if( m  <= m_min_jj  ) m_min_jj  = m;
	    if( dR >= dR_max_jj ) dR_max_jj = dR;
	    if( m  >= m_max_jj  ) m_max_jj  = m;
	  }
	  
	}
      }
      
      (corr_rnd_map.at(o))[0] = dR_min_bb;
      (corr_rnd_map.at(o))[1] = dR_min_jj;
      (corr_rnd_map.at(o))[2] = dR_max_bb;
      (corr_rnd_map.at(o))[3] = dR_max_jj;
      (corr_rnd_map.at(o))[4] = m_min_bb;
      (corr_rnd_map.at(o))[5] = m_min_jj;
      (corr_rnd_map.at(o))[6] = m_max_bb;
      (corr_rnd_map.at(o))[7] = m_max_jj;
    }


    rnd->next_event();
    for( auto c : cats_mem    ) delete c;
    for( auto o : objects_mem ) delete o;
    
    
    tout->Fill();
  }
  
  
  cout << "Done!" << endl;
  delete rnd;
  f->Close();
 
  fout->cd();
  tout->Write();
  fout->Close();
}
