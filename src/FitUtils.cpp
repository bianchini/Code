#include "TTH/MEIntegratorStandalone/interface/FitUtils.h"

FitUtils::FitUtils(){

  file = TFile::Open("fit.root","RECREATE");
  tree = new TTree("tree", "tree");

  // save toy
  tree->Branch("nsgn_gen", &nsgn_gen, "nsgn_gen/F");
  tree->Branch("nsgn_fit", &nsgn_fit, "nsgn_fit/F");
  tree->Branch("nsgn_err", &nsgn_err, "nsgn_err/F");
  tree->Branch("nbkg_gen", &nbkg_gen, "nbkg_gen/F");
  tree->Branch("nbkg_fit", &nbkg_fit, "nbkg_fit/F");
  tree->Branch("nbkg_err", &nbkg_err, "nbkg_err/F");
  tree->Branch("minNll", &minNll, "minNll/F");
  tree->Branch("edm", &edm, "edm/F");
  tree->Branch("status", &status, "status/I");

  x = new RooRealVar("x","x",-10.,10.);
  x->setRange(-4.,4.);
  x->setBins(10, "plot");
  x->setBins(20, "pdf");
  x->setBins(20, "fit");

  x_transf = new RooRealVar("x_transf","x_transf",0.,1.);
  x_transf->setRange(0.,1.);
  x_transf->setBins(10, "plot");
  x_transf->setBins(20, "pdf");
  x_transf->setBins(20, "fit");

  pdf_sgn_transf = nullptr;
  pdf_sgn = nullptr;
  pdf_bkg_transf = nullptr;
  pdf_bkg = nullptr;
  transformation = nullptr;

  norm_sgn = 1.0;
  norm_bkg = 1.0;

  savePlots = 100;
}

FitUtils::~FitUtils(){
  delete x;
  delete x_transf;
  delete pdf_sgn;
  delete pdf_bkg;
  delete transformation;
  delete pdf_sgn_transf;
  delete pdf_bkg_transf;

  file->Close();
  delete file;
  //delete tree;
  cout << "Done!" << endl;
}

void FitUtils::save_plots(const int& save){
  savePlots = save;
}

void FitUtils::addPdf(TTree* tr, const Process::Process& proc, const double& norm, const Option::Option& doBinned){

  string proc_name = "";
  switch(proc){
  case Process::Signal:
    proc_name = "sgn";
    break;
  case Process::Background:
    proc_name = "bkg";
    break;
  case Process::Transform:
    proc_name = "transf";
    break;
  default:
    cout << "*** Not a valid process name. Return! ***" << endl;
    return;
    break;
  }

  cout << "FitUtils::addPdf(): " << "adding process name " << proc_name << endl;

  // the p.d.f. of variable x
  RooAbsPdf* pdf = nullptr;
  // the p.d.f. of variable x_transf
  RooAbsPdf* pdf_transf = nullptr;

  RooDataSet data("data","data", tr, *x);
  data.Print("v");
  x->setBins(x->getBins("pdf"));
  RooDataHist data_hist("data_hist", "data_hist", *x, data); 

  if(doBinned==Option::Binned)
    pdf = new RooHistPdf(("pdf_"+proc_name).c_str(),("pdf_"+proc_name).c_str(), *x, data_hist);
  else
    pdf = new RooKeysPdf(("pdf_"+proc_name).c_str(),("pdf_"+proc_name).c_str(), *x, data, RooKeysPdf::MirrorBoth, 1.0);

  if(proc==Process::Transform){
    transformation = pdf->createCdf(*x);
    file->cd();
    if(savePlots){
      RooPlot* frame = x->frame();
      transformation->plotOn(frame);
      frame->Write("Cumulative_Pdf", TObject::kOverwrite);   
    }
    return;
  }
  
  RooDataSet data_transf("data_transf", "data_transf", RooArgSet(*x_transf));
  for(int i = 0 ; i < data.numEntries(); ++i){
    //if(i%500==0) cout << "\tProcessing event... " << i << endl;
    double input = ((RooAbsReal*)(data.get(i)->find("x")))->getVal();
    x->setVal(input);
    double output = transformation->getVal(*x);
    x_transf->setVal(output);
    data_transf.add(*x_transf);
  }

  data_transf.Print("v");
  x_transf->setBins(x_transf->getBins("pdf"));
  RooDataHist* data_transf_hist = new RooDataHist(("data_transf_hist_"+proc_name).c_str(), "data_transf_hist", *x_transf, data_transf); 
  pdf_transf = new RooHistPdf(("pdf_"+proc_name+"_transf").c_str(),("pdf_"+proc_name+"_transf").c_str(), *x_transf, *data_transf_hist);

  switch(proc){
  case Process::Signal:
    norm_sgn = norm;
    pdf_sgn = pdf;
    pdf_sgn_transf = pdf_transf;
    break;
  case Process::Background:
    norm_bkg = norm;
    pdf_bkg = pdf;
    pdf_bkg_transf = pdf_transf;
    break;
  default:
    break;
  }

  if(savePlots){
    file->cd();
    RooPlot* frame = x->frame( Bins(x->getBins("pdf")) );
    RooPlot* frame_transf = x_transf->frame( Bins(x->getBins("pdf")) );
    data.plotOn(frame);
    pdf->plotOn(frame, LineStyle(kSolid)); 
    frame->Write(("Pdf_"+proc_name).c_str(), TObject::kOverwrite);   
    data_transf.plotOn(frame_transf); 
    pdf_transf->plotOn(frame_transf, LineStyle(kSolid)); 
    frame_transf->Write(("Pdf_"+proc_name+"_transf").c_str(), TObject::kOverwrite);   
  }

  return;
}

void FitUtils::run_test(const int& ntoys, const Option::Option& doBinned, const size_t& degree){

  RooRealVar nsig("nsig","true fraction of sgn", norm_sgn);
  RooRealVar nbkg("nbkg","true fraction of bkg", norm_bkg);

  RooExtendPdf epdf_sgn("epdf_sgn", "epdf_sgn", *pdf_sgn, nsig);
  RooExtendPdf epdf_bkg("epdf_bkg", "epdf_bkg", *pdf_bkg, nbkg);

  RooAddPdf model("model","model", RooArgList(epdf_sgn,epdf_bkg) ) ;

  // Analytical Bkg Fit function 
  RooRealVar a0("a0","a0",1,-30.,30.) ;
  RooRealVar a1("a1","a1",1,-30.,30.) ;
  RooRealVar a2("a2","a2",1,-30.,30.) ;
  RooRealVar a3("a3","a3",1,-30.,30.) ;
  RooRealVar a4("a4","a4",1,-30.,30.) ;
  RooArgList coeff = RooArgList();
  switch(degree){
  case 1: coeff.add(RooArgList(a0,a1)); break;
  case 2: coeff.add(RooArgList(a0,a1,a2)); break;
  case 3: coeff.add(RooArgList(a0,a1,a2,a3)); break;
  case 4: coeff.add(RooArgList(a0,a1,a2,a3,a4)); break;
  default:
    cout << "Too large degree!" << endl;
    return;
  }

  RooBernstein fit_pdf_bkg_transf("fit_pdf_bkg_transf","fit_pdf_bkg_transf", *x_transf, coeff);
  //RooChebychev fit_pdf_bkg_transf("fit_pdf_bkg_transf","fit_pdf_bkg_transf", *x_transf, coeff);

  // Sum the signal components into a composite signal p.d.f. (transformed) 
  RooRealVar ns_fit("ns_fit","fraction of signal", norm_sgn, -100, 10000.) ;
  RooRealVar nb_fit("nb_fit","fraction of background", norm_bkg, -100, 10000.) ;

  // The S+B model (transformed) 
  RooAddPdf model_transf("model_transf","model_transf", 
			 RooArgList(*pdf_sgn_transf, fit_pdf_bkg_transf), 
			 RooArgList(ns_fit,nb_fit)) ;


  int toy = 0;
  while(toy<ntoys){
    cout << "Generate toy n. " << toy << endl;

    RooDataSet *data = model.generate(*x, Extended()) ;
    cout << "\tNumber of entries: " << data->numEntries() << endl;

    RooDataSet data_transf("data_tr", "data_tr", RooArgSet(*x_transf));
    for(int i = 0 ; i < data->numEntries(); ++i){
      double input = ((RooAbsReal*)(data->get(i)->find("x")))->getVal();
      x->setVal(input);
      double output = transformation->getVal(*x);
      x_transf->setVal(output);
      data_transf.add(*x_transf);
    }

    RooFitResult* r = nullptr;
    if(doBinned==Option::Binned){
      x_transf->setBins(x_transf->getBins("fit"));
      RooDataHist data_transf_hist("data_transf_hist", "data_transf_hist", *x_transf, data_transf); 
      r = model_transf.fitTo(data_transf_hist, Extended(kTRUE), Save());
    }
    else{
      r = model_transf.fitTo(data_transf, Extended(kTRUE), Save());
    }

    if( r->status()==0 ){
      RooRealVar* ns_fit = (RooRealVar*) model_transf.getVariables()->find("ns_fit");
      RooRealVar* nb_fit = (RooRealVar*) model_transf.getVariables()->find("nb_fit");
      nsgn_gen = nsig.getVal();
      nsgn_fit = ns_fit->getVal();
      nsgn_err = ns_fit->getError();
      nbkg_gen = nbkg.getVal();
      nbkg_fit = nb_fit->getVal();
      nbkg_err = nb_fit->getError();
      minNll = r->minNll();
      edm = r->edm();
      status = r->status();
      tree->Fill();
    }

    if(savePlots && toy%savePlots==0){
      RooPlot* frame_transf = x_transf->frame( Bins(x_transf->getBins("fit")) );
      data_transf.plotOn(frame_transf);
      model_transf.plotOn(frame_transf, LineStyle(kSolid), LineColor(kRed));
      model_transf.plotOn(frame_transf, Components(fit_pdf_bkg_transf), LineStyle(kDashed), LineColor(kGreen));
      model_transf.plotOn(frame_transf, Components(*pdf_sgn_transf), LineStyle(kSolid), LineColor(kBlue));
      file->cd();
      frame_transf->Write(Form("x_transf_frame_toy%d",toy), TObject::kOverwrite);
      RooPlot* frame = x->frame( Bins(x_transf->getBins("fit")) );
      data->plotOn(frame);
      model.plotOn(frame, LineStyle(kSolid), LineColor(kRed));
      model.plotOn(frame, Components(epdf_bkg), LineStyle(kDashed), LineColor(kGreen));
      model.plotOn(frame, Components(epdf_sgn), LineStyle(kSolid), LineColor(kBlue));
      frame->Write(Form("x_frame_toy%d",toy), TObject::kOverwrite);  
    }

    ++toy;
  }

  file->cd();
  tree->Write("", TObject::kOverwrite);

  return;
}
