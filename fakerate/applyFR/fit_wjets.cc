using namespace RooFit;

void fit_wjets()
{
	Double_t Lumi = 2759.017;

	//Get ROOT Files
	TFile *ttbar = new TFile("histograms/fake1.root");
	TFile *DYM50 = new TFile("histograms/fake0.root");
	TFile *DYM10to50 = new TFile("histograms/fake100.root");
	TFile *data = new TFile("histograms/fake6.root");

	//Get Histograms
	TH1D *h_ttbar = (TH1D*)ttbar->Get( "fitWJets1" );
	h_ttbar->Scale(Lumi*831.76/187625980);

	TH1D *h_DYJets = (TH1D*)DYM50->Get( "fitWJets1" );
	TH1D *h_DYJets_M10to50 = (TH1D*)DYM10to50->Get( "fitWJets1" );
	h_DYJets->Scale(Lumi*2008.4*3/8.1236e+07);
	h_DYJets_M10to50->Scale(Lumi*18610/4.52115e+07);
	h_DYJets->Add(h_DYJets_M10to50);

	TH1D *h_WJets = (TH1D*)data->Get( "fitSameWJets1" );
	TH1D *h_QCD_SS = (TH1D*)data->Get( "fitSameDijet1" );
	TH1D *h_ttbar_SS = (TH1D*)ttbar->Get( "fitSameWJets1" );
	h_ttbar_SS->Scale(Lumi*831.76/187625980);
	h_WJets->Add(h_QCD_SS,-2.0);
	h_WJets->Add(h_ttbar_SS,-1.0);

	TH1D *h_QCD = (TH1D*)data->Get( "fitDijet1" );
	TH1D *h_DYJets_Dijet = (TH1D*)DYM50->Get( "fitDijet1" );
	TH1D *h_DYJets_M10to50_Dijet = (TH1D*)DYM10to50->Get( "fitDijet1" );
	TH1D *h_ttbar_Dijet = (TH1D*)ttbar->Get( "fitDijet1" );

	h_DYJets_Dijet->Scale(Lumi*2008.4*3/8.1236e+07);
	h_DYJets_M10to50_Dijet->Scale(Lumi*18610/4.52115e+07);
	h_ttbar_Dijet->Scale(Lumi*831.76/187625980);
	h_QCD->Add(h_DYJets_Dijet,-1.0);
	h_QCD->Add(h_DYJets_M10to50_Dijet,-1.0);
	h_QCD->Add(h_ttbar_Dijet,-1.0);

	TH1D *h_data = (TH1D*)data->Get( "fitWJets1" );

	//Convert TH1D to RooDataHist
	RooRealVar mass("mass", "Dimuon mass [GeV]", 15,200);

	RooDataHist *RooHist_ttbar = new RooDataHist("RooHist_ttbar", "RooHistogram_ttbar", mass, h_ttbar);
	RooDataHist *RooHist_DYJets = new RooDataHist("RooHist_DYJets", "RooHistogram_DYJets", mass, h_DYJets);
	RooDataHist *RooHist_WJets = new RooDataHist("RooHist_WJets", "RooHistogram_WJets", mass, h_WJets);
	RooDataHist *RooHist_QCD = new RooDataHist("RooHist_QCD", "RooHistogram_QCD", mass, h_QCD);
	RooDataHist *RooHist_data = new RooDataHist("RooHist_data", "RooHistogram_data", mass, h_data);

	//Convert RooDataHist to RooHistPdf
	RooHistPdf *pdf_ttbar = new RooHistPdf("pdf_ttbar", "Template from ttbar MC", mass, *RooHist_ttbar, 0);
	RooHistPdf *pdf_DYJets = new RooHistPdf("pdf_DYJets", "Template from DYJets MC", mass, *RooHist_DYJets, 0);
	RooHistPdf *pdf_WJets = new RooHistPdf("pdf_WJets", "Template from same-sign WJets", mass, *RooHist_WJets, 0);
	RooHistPdf *pdf_QCD = new RooHistPdf("pdf_QCD", "Template from data-driven dijet", mass, *RooHist_QCD, 0);

	// Construct model = n_ttbar * ttbar + n_WJets * WJets
	double N_ttbar = h_ttbar->Integral();
	double N_DYJets = h_DYJets->Integral();
	double N_WJets = h_WJets->Integral()*3;
	double N_QCD = h_QCD->Integral()*2;
	double ratio = h_data->Integral()/(N_ttbar + N_DYJets + N_WJets + N_QCD);
	N_ttbar *= ratio;
	N_DYJets *= ratio;
	N_WJets *= ratio;
	N_QCD *= ratio;

	RooRealVar n_ttbar("n_ttbar", "n_ttbar", N_ttbar, N_ttbar*0.9, N_ttbar*1.1);
	RooRealVar n_DYJets("n_DYJets", "n_DYJets", N_DYJets, N_DYJets*0.9, N_DYJets*1.1);
	RooRealVar n_WJets("n_WJets", "n_WJets", N_WJets, N_WJets*0.7, N_WJets*1.3);
	RooRealVar n_QCD("n_QCD", "n_QCD", N_QCD, N_QCD*0.7, N_QCD*1.3);
	RooAddPdf model( "model","model",RooArgList(*pdf_WJets, *pdf_QCD, *pdf_DYJets, *pdf_ttbar), RooArgList(n_WJets, n_QCD, n_DYJets, n_ttbar) );

	// Fit to data
	RooFitResult* r = model.fitTo( *RooHist_data, Save() );
	r->Print();
	RooAbsReal *chi2 = model.createChi2(*RooHist_data);
	cout << "chi2ndof: " << chi2->getVal() / ((Double_t)h_data->GetNbinsX()) << endl;
}
