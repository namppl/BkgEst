#include "../interface/tdrstyle.C"
#include "../interface/CMS_lumi.C"
using namespace RooFit;

void fitTemplates(const TString& category)
{
	//Get ROOT Files
	TFile *f0 = new TFile("histograms/hist0.root");
	TFile *f1 = new TFile("histograms/hist1.root");
	TFile *f2 = new TFile("histograms/hist2.root");
	TFile *f5 = new TFile("histograms/hist5.root");
	TFile *f6 = new TFile("histograms/hist6.root");
	TFile *f11 = new TFile("histograms/hist11.root");
	TFile *f12 = new TFile("histograms/hist12.root");
	TFile *f13 = new TFile("histograms/hist13.root");

	//Get Histograms
	TH1D *h_ttbar = (TH1D*)f1->Get( "denominator_" + category );
	TH1D *h_WJets = (TH1D*)f2->Get( "denominator_" + category );
	TH1D *h_DYJets = (TH1D*)f0->Get( "denominator_" + category );
	TH1D *h_QCD = (TH1D*)f5->Get( "denominator_" + category );
	TH1D *h_WW = (TH1D*)f11->Get( "denominator_" + category );
	TH1D *h_WZ = (TH1D*)f12->Get( "denominator_" + category );
	TH1D *h_ZZ = (TH1D*)f13->Get( "denominator_" + category );

	TH1D *h_data = (TH1D*)f6->Get( "denominator_" + category );

	//Convert TH1D to RooDataHist
	RooRealVar obs("obs", "TrkIso/p_{T}", 0, 5);

	RooDataHist *RooHist_ttbar = new RooDataHist("RooHist_ttbar", "RooHistogram_ttbar", obs, h_ttbar);
	RooDataHist *RooHist_WJets = new RooDataHist("RooHist_WJets", "RooHistogram_WJets", obs, h_WJets);
	RooDataHist *RooHist_DYJets = new RooDataHist("RooHist_DYJets", "RooHistogram_DYJets", obs, h_DYJets);
	RooDataHist *RooHist_QCD = new RooDataHist("RooHist_QCD", "RooHistogram_QCD", obs, h_QCD);
	RooDataHist *RooHist_data = new RooDataHist("RooHist_data", "RooHistogram_data", obs, h_data);
	RooDataHist *RooHist_WW = new RooDataHist("RooHist_WW", "RooHistogram_WW", obs, h_WW);
	RooDataHist *RooHist_WZ = new RooDataHist("RooHist_WZ", "RooHistogram_WZ", obs, h_WZ);
	RooDataHist *RooHist_ZZ = new RooDataHist("RooHist_ZZ", "RooHistogram_ZZ", obs, h_ZZ);

	//Convert RooDataHist to RooHistPdf
	RooHistPdf *pdf_ttbar = new RooHistPdf("pdf_ttbar", "Template from ttbar MC", obs, *RooHist_ttbar, 0);
	RooHistPdf *pdf_WJets = new RooHistPdf("pdf_WJets", "Template from WJets MC", obs, *RooHist_WJets, 0);
	RooHistPdf *pdf_DYJets = new RooHistPdf("pdf_DYJets", "Template from DYJets MC", obs, *RooHist_DYJets, 0);
	RooHistPdf *pdf_QCD = new RooHistPdf("pdf_QCD", "Template from QCD MC", obs, *RooHist_QCD, 0);
	RooHistPdf *pdf_WW = new RooHistPdf("pdf_WW", "Template from WW MC", obs, *RooHist_WW, 0);
	RooHistPdf *pdf_WZ = new RooHistPdf("pdf_WZ", "Template from WZ MC", obs, *RooHist_WZ, 0);
	RooHistPdf *pdf_ZZ = new RooHistPdf("pdf_ZZ", "Template from ZZ MC", obs, *RooHist_ZZ, 0);

	// Construct model = n_ttbar * ttbar + n_WJets * WJets

	Double_t Lumi = 2765;
	Double_t xsec_ttbar = 831.76; Double_t Nprocessed_ttbar = 187625980; Double_t Npass_ttbar = h_ttbar->Integral();
	Double_t NN_ttbar = ((xsec_ttbar * Lumi) / Nprocessed_ttbar) * Npass_ttbar;
	cout << "N_ttbar: "<< Npass_ttbar << endl;
	cout << "N_ttbar: "<< NN_ttbar << endl;

	Double_t xsec_DYJets = 3*2008.4; Double_t Nprocessed_DYJets = 8.1236e+07; Double_t Npass_DYJets = h_DYJets->Integral();
	Double_t NN_DYJets = ((xsec_DYJets * Lumi) / Nprocessed_DYJets) * Npass_DYJets;
	cout << "N_DYJets: "<< NN_DYJets << endl;

	Double_t xsec_WJets = 60290; Double_t Nprocessed_WJets = 1.65208e+07; Double_t Npass_WJets = h_WJets->Integral();
	Double_t NN_WJets = ((xsec_WJets * Lumi) / Nprocessed_WJets) * Npass_WJets;
	cout << "N_WJets: "<< NN_WJets << endl;

	Double_t xsec_WW = 63.21; Double_t Nprocessed_WW = 988416; Double_t Npass_WW = h_WW->Integral();
	Double_t NN_WW = ((xsec_WW * Lumi) / Nprocessed_WW) * Npass_WW;
	cout << "N_WW: "<< Npass_WW << endl;
	cout << "N_WW: "<< NN_WW << endl;

	Double_t xsec_WZ = 22.82; Double_t Nprocessed_WZ = 999996; Double_t Npass_WZ = h_WZ->Integral();
	Double_t NN_WZ = ((xsec_WZ * Lumi) / Nprocessed_WZ) * Npass_WZ;
	cout << "N_WZ: "<< NN_WZ << endl;

	Double_t xsec_ZZ = 10.32; Double_t Nprocessed_ZZ = 985598; Double_t Npass_ZZ = h_ZZ->Integral();
	Double_t NN_ZZ = ((xsec_ZZ * Lumi) / Nprocessed_ZZ) * Npass_ZZ;
	cout << "N_ZZ: "<< NN_ZZ << endl;

	Double_t NN_QCD = Lumi * h_QCD->Integral();
	cout << "N_QCD: "<< NN_QCD << endl;

	double N_total = NN_ttbar + NN_DYJets + NN_WJets + NN_QCD + NN_WW + NN_WZ + NN_ZZ;
	double N_ttbar = h_data->Integral()*NN_ttbar/N_total;
	double N_DYJets = h_data->Integral()*NN_DYJets/N_total;
	double N_WJets = h_data->Integral()*NN_WJets/N_total;
	double N_QCD = h_data->Integral()*NN_QCD/N_total;
	double N_WW = h_data->Integral()*NN_WW/N_total;
	double N_WZ = h_data->Integral()*NN_WZ/N_total;
	double N_ZZ = h_data->Integral()*NN_ZZ/N_total;
	// Double_t xsec_DYJets = 6025.2; Double_t Nprocessed_DYJets = 2000000; Double_t Npass_DYJets = h_DYJets->Integral();
	// Double_t N_DYJets = ((xsec_DYJets * Lumi) / Nprocessed_DYJets) * Npass_DYJets;
	// cout << "N_DYJets: "<< N_DYJets << " Range: " << N_DYJets*0.5 << " " << N_DYJets*1.5 << endl;

	RooRealVar n_ttbar("n_ttbar", "n_ttbar", N_ttbar, N_ttbar*0.95, N_ttbar*1.05);
	RooRealVar n_WJets("n_WJets", "n_WJets", N_WJets, N_WJets*0.5, N_WJets*1.5);
	RooRealVar n_DYJets("n_DYJets", "n_DYJets", N_DYJets, N_DYJets*0.95, N_DYJets*1.05);
	RooRealVar n_QCD("n_QCD", "n_QCD", N_QCD, N_QCD*0.5, N_QCD*1.5);
	RooRealVar n_WW("n_WW", "n_WW", N_WW, 0.95*N_WW, N_WW*1.05);
	RooRealVar n_WZ("n_WZ", "n_WZ", N_WZ, 0.95*N_WZ, N_WZ*1.05);
	RooRealVar n_ZZ("n_ZZ", "n_ZZ", N_ZZ, 0.95*N_ZZ, N_ZZ*1.05);
  	RooAddPdf model( "model","model",RooArgList(*pdf_QCD, *pdf_WJets, *pdf_DYJets, *pdf_ttbar, *pdf_WW, *pdf_WZ, *pdf_ZZ), RooArgList(n_QCD, n_WJets, n_DYJets, n_ttbar, n_WW, n_WZ, n_ZZ) );
  	//RooAddPdf model( "model","model", RooArgList(*pdf_ttbar, *pdf_WJets, *pdf_QCD), RooArgList(n_ttbar, n_WJets, n_QCD) );

  	//RooFitResult* r = pdf_WJets->fitTo( *RooHist_data, Save() );
  	RooFitResult* r = model.fitTo( *RooHist_data, Save() );

  	TCanvas *c_fit = new TCanvas("c_fit", "", 800, 800);
  	c_fit->cd();

  	//Top Pad
  	TPad *c1_1 = new TPad("padc1_1","padc1_1",0.01,0.01,0.99,0.99);
  	c1_1->Draw();
  	c1_1->cd();
  	c1_1->SetTopMargin(0.01);
  	c1_1->SetBottomMargin(0.25);
  	c1_1->SetRightMargin(0.03);
  	c1_1->SetLeftMargin(0.09);
  	c1_1->SetFillStyle(0.01);
	c1_1->SetLogy();

	RooPlot* frame1 = obs.frame( Title(" ") ) ;
	// pdf_ttbar->plotOn( frame1, LineColor(kOrange) );
	// pdf_WJets->plotOn(frame1, LineColor(kGreen) );
	RooHist_data->plotOn(frame1, DataError(RooAbsData::SumW2));
	model.plotOn(frame1, Components("pdf_ZZ,pdf_WZ,pdf_WW,pdf_ttbar,pdf_DYJets,pdf_WJets,pdf_QCD"), LineColor(0), FillColor(7), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_ZZ,pdf_WZ,pdf_WW,pdf_ttbar,pdf_DYJets,pdf_WJets"), LineColor(0), FillColor(4), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_ZZ,pdf_WZ,pdf_WW,pdf_ttbar,pdf_DYJets"), LineColor(0), FillColor(2), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_ZZ,pdf_WZ,pdf_WW,pdf_ttbar"), LineColor(0), FillColor(3), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_ZZ,pdf_WZ,pdf_WW"), LineColor(0), FillColor(13), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_ZZ,pdf_WZ"), LineColor(0), FillColor(14), DrawOption("F") );
	model.plotOn(frame1, Components("pdf_ZZ"), LineColor(0), FillColor(15), DrawOption("F") );
	//model.plotOn(frame1, Components("pdf_WJets"), LineColor(0), FillColor(kGreen), DrawOption("F") );
	//model.plotOn(frame1, Components("pdf_QCD,pdf_WJets,pdf_ttbar"), LineColor(kRed), FillColor(kRed), DrawOption("F") );
	//model.plotOn(frame1, Components("pdf_QCD,pdf_WJets"), LineColor(kBlue), FillColor(kBlue), DrawOption("F") );
	//model.plotOn(frame1, Components("pdf_QCD"), LineColor(kGreen), FillColor(kGreen), DrawOption("F") );
	//model.plotOn(frame1, LineColor(1) );
	RooHist_data->plotOn(frame1, DataError(RooAbsData::SumW2));
	//model.paramOn(frame1, Layout(0.65,0.9,0.9) );
	frame1->Draw();
	r->Print();

	TLegend *leg1 = new TLegend(0.65,0.7,.95,.97);
	leg1->SetFillColor(kWhite);
	leg1->SetLineColor(kWhite);
	leg1->AddEntry(frame1->nameOf(0),"Data", "EP");
	leg1->AddEntry(frame1->nameOf(1),"QCD","F");
	leg1->AddEntry(frame1->nameOf(2),"WJets","F");
	leg1->AddEntry(frame1->nameOf(3),"DYJets","F");
	leg1->AddEntry(frame1->nameOf(4),"ttbar","F");
	leg1->AddEntry(frame1->nameOf(5),"WW","F");
	leg1->AddEntry(frame1->nameOf(6),"WZ","F");
	leg1->AddEntry(frame1->nameOf(7),"ZZ","F");
	leg1->Draw();


	frame1->GetYaxis()->SetTitle("Entry");
	frame1->GetXaxis()->SetLabelSize(0);

	TH1D *h_MC = (TH1D*)model.createHistogram("h_MC", obs);
	// h_MC->Sumw2();
	// TCanvas *c_MC = new TCanvas("c_MC", "", 700, 700);
	
	Double_t Ndata = h_data->Integral();
	Double_t NMC = h_MC->Integral();
	h_MC->Scale(Ndata / NMC );
	// cout << "# data: " << Ndata << endl;
	// h_MC->Draw(); 
	// h_data->Draw("SAMEEP");


	//Bottom Pad
	TPad *c1_2 = new TPad("padc1_2","padc1_2",0.01,0.01,0.99,0.25);
	c1_2->Draw();
	c1_2->cd();
	c1_2->SetTopMargin(0.1);
	c1_2->SetBottomMargin(0.30);
	c1_2->SetRightMargin(0.02);
	c1_2->SetLeftMargin(0.08);
	c1_2->SetFillStyle(0);
	c1_2->SetGrid();

	//Make ratio plot
	TH1D *h_ratio = (TH1D*)h_data->Clone();
	h_data->Sumw2(); h_MC->Sumw2();
	h_ratio->Divide(h_data, h_MC);
	h_ratio->SetTitle("");
	h_ratio->GetXaxis()->SetMoreLogLabels();
	h_ratio->GetXaxis()->SetNoExponent();
	h_ratio->GetXaxis()->SetTitle( "TrkIso/p_{T}" );
	h_ratio->GetYaxis()->SetTitle("data/MC");
	h_ratio->GetXaxis()->SetTitleSize(0.13);
	h_ratio->GetYaxis()->SetTitleSize(0.09);
	h_ratio->GetYaxis()->SetTitleOffset(0.4);
	h_ratio->GetXaxis()->SetLabelSize(0.11);
	h_ratio->GetYaxis()->SetLabelSize(0.07);
	h_ratio->GetYaxis()->SetTickLength(0.015);
	h_ratio->SetMaximum( 1.3 );
	h_ratio->SetMinimum( 0.7 );
	h_ratio->SetMarkerSize(0.3);
	h_ratio->SetStats(kFALSE);

	h_ratio->Draw("e1p");
	
	TH1D *h_line = (TH1D*)h_data->Clone();
	h_line->Reset("ICES");
	Int_t Nbins = h_line->GetNbinsX();
	for(Int_t i_bin=0; i_bin< Nbins; i_bin++)
		h_line->SetBinContent(i_bin+1, 1);

	h_line->SetLineColor(kRed);
	h_line->Draw("LSAME");

	//leg1->Draw();

	RooAbsReal *chi2 = model.createChi2(*RooHist_data);
	cout << "chi2: " << chi2->getVal() << endl;
	cout << "Normalized chi2: " << chi2->getVal() / ((Double_t)h_data->GetNbinsX()) << endl;

	c_fit->Print("print/fit_"+category+".pdf");
}

