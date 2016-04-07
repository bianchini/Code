#include "TTH/MEIntegratorStandalone/interface/FitUtils.h"

FitUtils::FitUtils(){

  f = TFile::Open("fit.root","RECREATE");

  x = new RooRealVar("x","x",-10.,10.);
  x->setRange(-4.,4.);
  x->setBins(20, "plot");
  x->setBins(100, "pdf");

  x_transf = new RooRealVar("x_transf","x_transf",0.,1.);
  x_transf->setRange(0.,1.);
  x_transf->setBins(10, "plot");
  x_transf->setBins(20, "pdf");

  h_pdf_sgn_transf = new TH1D("h_pdf_sgn_transf", "h_pdf_sgn_transf", 20, 0., 1.);
  h_pdf_sgn_transf->Sumw2();
  h_pdf_bkg_transf = new TH1D("h_pdf_bkg_transf", "h_pdf_bkg_transf", 50, 0., 1.);
  h_pdf_bkg_transf->Sumw2();
  
  pdf_sgn_transf = nullptr;
  pdf_sgn = nullptr;
  pdf_bkg_transf = nullptr;
  pdf_bkg = nullptr;

  savePlots = true;
}

FitUtils::~FitUtils(){
  delete x;
  delete x_transf;
  delete h_pdf_sgn;
  delete h_pdf_bkg;
  delete h_pdf_sgn_transf;
  delete h_pdf_bkg_transf;
  delete cum_pdf;
  delete pdf_sgn;
  delete pdf_bkg;
  delete transformation;
  delete pdf_sgn_transf;
  delete pdf_bkg_transf;

  f->Close();
  delete f;
  cout << "Done!" << endl;
}

void FitUtils::addPdfSgn(TH1D* h){
  
  h_pdf_sgn = new TH1D(*h);
  for(int i = 0 ; i < 20000 ; ++i){
    double val = h->GetRandom();
    double val_transf = cum_pdf->Interpolate(val);
    h_pdf_sgn_transf->Fill(val_transf);
  }
  h_pdf_sgn_transf->Scale(h_pdf_sgn->Integral()/h_pdf_sgn_transf->Integral());
  cout << "Sgn template filled with...." << h_pdf_sgn->Integral() << " entries" << endl;
  f->cd();
  h_pdf_sgn->Write("Sgn_Pdf", TObject::kOverwrite);
  h_pdf_sgn_transf->Write("Sgn_Pdf_Transformed", TObject::kOverwrite);
}

void FitUtils::addPdfBkg(TH1D* h){
  h_pdf_bkg = new TH1D(*h);
  for(int i = 0 ; i < h_pdf_bkg->GetEntries() ; ++i){
    double val = h->GetRandom();
    double val_transf = cum_pdf->Interpolate(val);
    h_pdf_bkg_transf->Fill(val_transf);
  }
  cout << "Bkg template filled with...." << h_pdf_bkg->Integral() << " entries" << endl;
  f->cd();
  h_pdf_bkg->Write("Bkg_Pdf", TObject::kOverwrite);
  h_pdf_bkg_transf->Write("Bkg_Pdf_Transformed", TObject::kOverwrite);
}

void FitUtils::addTransformPdf(TH1D* h){
  cum_pdf = (TH1D*)h->GetCumulative();
  cum_pdf->Scale(1./h->Integral());
  f->cd();
  h->Write("Transform_Pdf", TObject::kOverwrite);
  cum_pdf->Write("Cumulative_Pdf", TObject::kOverwrite);  
}

void FitUtils::addPdf(TTree* tree, const string& proc, const bool& doBinned){

  cout << "FitUtils::addPdf(): " << "adding process name " << proc << endl;

  // the p.d.f. of variable x
  RooAbsPdf* pdf = nullptr;
  // the p.d.f. of variable x_transf
  RooAbsPdf* pdf_transf = nullptr;

  RooDataSet data("data","data", tree, *x);
  //data.Print("v");
  x->setBins(x->getBins("pdf"));
  RooDataHist data_hist("data_hist", "data_hist", *x, data); 

  if(doBinned)
    pdf = new RooHistPdf(("pdf_"+proc).c_str(),("pdf_"+proc).c_str(), *x, data_hist);
  else
    pdf = new RooKeysPdf(("pdf_"+proc).c_str(),("pdf_"+proc).c_str(), *x, data);

  if(proc=="transf"){
    transformation = pdf->createCdf(*x);
    f->cd();
    RooPlot* frame = x->frame();
    transformation->plotOn(frame);
    frame->Write("Cumulative_Pdf", TObject::kOverwrite);   
    return;
  }
  
  RooDataSet data_transf("data_transf", "data_transf", RooArgSet(*x_transf));
  for(int i = 0 ; i < data.numEntries(); ++i){
    if(i%500==0) cout << "\tProcessing event... " << i << endl;
    double input = ((RooAbsReal*)(data.get(i)->find("x")))->getVal();
    x->setVal(input);
    double output = transformation->getVal(*x);
    x_transf->setVal(output);
    data_transf.add(*x_transf);
  }

  //data_transf.Print("v");
  x_transf->setBins(x_transf->getBins("plot"));
  RooDataHist* data_transf_hist = new RooDataHist(("data_transf_hist_"+proc).c_str(), "data_transf_hist", *x_transf, data_transf); 
  pdf_transf = new RooHistPdf(("pdf_"+proc+"_transf").c_str(),("pdf_"+proc+"_transf").c_str(), *x_transf, *data_transf_hist);

  if(proc=="sgn"){
    pdf_sgn = pdf;
    pdf_sgn_transf = pdf_transf;
    cout << "Pdf_sgn loaded" << endl;
  }
  else if(proc=="bkg"){
    pdf_bkg = pdf;
    pdf_bkg_transf = pdf_transf;
    cout << "Pdf_bkg loaded" << endl;
  }
  else{
    cout << "*** Not a valid process name. Return! ***" << endl;
    return;
  }

  f->cd();
  RooPlot* frame = x->frame();
  RooPlot* frame_transf = x_transf->frame();
  data.plotOn(frame);
  pdf->plotOn(frame, LineStyle(kSolid) ); 
  frame->Write(("Pdf_"+proc).c_str(), TObject::kOverwrite);   
  data_transf.plotOn(frame_transf); 
  pdf_transf->plotOn(frame_transf, LineStyle(kSolid) ); 
  frame_transf->Write(("Pdf_"+proc+"_transf").c_str(), TObject::kOverwrite);   

  return;
}

void FitUtils::run_test(const int& ntoys){

  // pull
  TH1F pull_nsig("pull_nsig", "", 100, -5,5);

  RooRealVar nsig("nsig","true fraction of sgn", 10.);
  RooRealVar nbkg("nbkg","true fraction of bkg", 100.);

  RooExtendPdf epdf_sgn("epdf_sgn", "epdf_sgn", *pdf_sgn, nsig);
  RooExtendPdf epdf_bkg("epdf_bkg", "epdf_bkg", *pdf_bkg, nbkg);

  RooAddPdf model("model","model", RooArgList(epdf_sgn,epdf_bkg) ) ;

  // Bkg pdf (transformed) 
  RooRealVar a0("a0","a0",1,-20.,20.) ;
  RooRealVar a1("a1","a1",1,-20.,20.) ;
  RooRealVar a2("a2","a2",1,-20.,20.) ;
  RooBernstein fit_pdf_bkg_transf("fit_pdf_bkg_transf","fit_pdf_bkg_transf", *x_transf, RooArgList(a0,a1,a2));

  // Sum the signal components into a composite signal p.d.f. (transformed) 
  RooRealVar nsig_fit("nsig_fit","fraction of signal",      10,0., 10000.) ;
  RooRealVar nbkg_fit("nbkg_fit","fraction of background", 100,0., 10000.) ;

  // The S+B model (transformed) 
  RooAddPdf model_transf("model_transf","model_transf", 
			 RooArgList(*pdf_sgn_transf, fit_pdf_bkg_transf), 
			 RooArgList(nsig_fit,nbkg_fit)) ;


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

    if(savePlots){
      RooPlot* frame_transf = x_transf->frame();
      data_transf.plotOn(frame_transf);
      model_transf.plotOn(frame_transf, Components(fit_pdf_bkg_transf), LineStyle(kDashed) );
      model_transf.plotOn(frame_transf, Components(*pdf_sgn_transf), LineStyle(kSolid) );
      f->cd();
      frame_transf->Write(Form("x_transf_frame_toy%d",toy), TObject::kOverwrite);
      RooPlot* frame = x->frame();
      data->plotOn(frame);
      model.plotOn(frame);
      frame->Write(Form("x_frame_toy%d",toy), TObject::kOverwrite);  
    }

    return;

    RooFitResult* r = model_transf.fitTo(data_transf, Extended(kTRUE), Save());
    RooRealVar* ns_fit = (RooRealVar*) model_transf.getVariables()->find("nsig_fit");
    cout << ns_fit->getVal() << endl;    

    if( r->status()==0 )
      pull_nsig.Fill( (ns_fit->getVal() - nsig.getVal() )/ns_fit->getError() );

    ++toy;
  }

  f->cd();
  pull_nsig.Write("pull_nsig",  TObject::kOverwrite);

  return;
}
