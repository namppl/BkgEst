#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <THStack.h>
#include <TMath.h>
#include <TText.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLorentzVector.h>
#include <TStopwatch.h>
#include <TColor.h>
#include <TLatex.h>
#include <TEfficiency.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

#include "../../interface/tdrstyle.C"
#include "../../interface/CMS_lumi.C"
using namespace std;

void estimateDijet() {

    int W = 1200;
    int H = 1200;

    int H_ref = 1200;
    int W_ref = 1200;

    // references for T, B, L, R
    float T = 0.08*H_ref;
    float B = 0.12*H_ref;
    float L = 0.12*W_ref;
    float R = 0.04*W_ref;

    // UPDATED IN 2017
    lumi_13TeV = "2759 pb^{-1}";
    // lumi_13TeV = "2833 pb^{-1}";
    lumiTextSize = 0.5;
    writeExtraText = true;
    extraText = "Preliminary";
    drawLogo = false;

    int binnum = 43;
    double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1500, 3000};

    TH1D* massFrame = new TH1D("massFrame","",38,15,3000);
    massFrame->SetMinimum(0.001);
    massFrame->SetMaximum(1000000);
    massFrame->SetStats(kFALSE);
    massFrame->GetXaxis()->SetTitle("Mass[GeV]");
    //massFrame->GetXaxis()->CenterTitle(kTRUE);
    //massFrame->GetYaxis()->CenterTitle(kTRUE);
    massFrame->GetYaxis()->SetTitleOffset(1);
    massFrame->GetYaxis()->SetTitle("Number of events");
    massFrame->GetXaxis()->SetTitleSize(0);
    massFrame->GetYaxis()->SetTitleSize(0.05);
    massFrame->GetXaxis()->SetLabelSize(0);
    //massFrame->GetYaxis()->SetLabelSize(0.025); 
    massFrame->GetXaxis()->SetMoreLogLabels();

    TFile* f[14];
    f[0] = new TFile("histograms/fake0.root","READ");
    f[3] = new TFile("histograms/fake100.root","READ");
    f[1] = new TFile("histograms/fake1.root","READ");
    f[6] = new TFile("histograms/fake6.root","READ");

    TH1D* wjets_template[14]; // just for draw MC histograms
    TH1D* dijet_template[14];
    TH1D* dijetSS_template[14];
    TH1D* dijet_ratio[14];
    TH1D* dijetSS_ratio[14];

    dijet_template[5] = (TH1D*)f[6]->Get("histDijet1");
    dijet_template[5]->SetFillColor(7);
    dijet_template[5]->SetStats(kFALSE);
    dijet_template[5]->Sumw2();

    dijetSS_template[5] = (TH1D*)f[6]->Get("histSameDijet1");
    dijetSS_template[1] = (TH1D*)f[1]->Get("histSameDijet1");
    dijetSS_template[5]->SetLineColor(1);
    dijetSS_template[5]->SetMarkerColor(1);
    dijetSS_template[5]->SetMarkerStyle(22);
    dijetSS_template[5]->SetMarkerSize(3);
    dijetSS_template[5]->SetStats(kFALSE);
    dijetSS_template[5]->Sumw2();

    /////////
    wjets_template[0] = (TH1D*)f[0]->Get("histWJets1");
    wjets_template[3] = (TH1D*)f[3]->Get("histWJets1");
    wjets_template[0]->SetFillColor(2);
    wjets_template[0]->SetStats(kFALSE);
    wjets_template[0]->Sumw2();
    wjets_template[3]->Sumw2();

    wjets_template[1] = (TH1D*)f[1]->Get("histWJets1");
    wjets_template[1]->SetFillColor(3);
    wjets_template[1]->SetStats(kFALSE);
    wjets_template[1]->Sumw2();
    /////////

    dijet_template[0] = (TH1D*)f[0]->Get("histDijet1");
    dijet_template[3] = (TH1D*)f[3]->Get("histDijet1");
    dijet_template[0]->SetFillColor(2);
    dijet_template[0]->SetStats(kFALSE);
    dijet_template[0]->Sumw2();
    dijet_template[3]->Sumw2();

    dijet_template[1] = (TH1D*)f[1]->Get("histDijet1");
    dijet_template[1]->SetFillColor(3);
    dijet_template[1]->SetStats(kFALSE);
    dijet_template[1]->Sumw2();
    cout<<"ttbar OS before = "<<dijet_template[1]->Integral()<<endl;
    cout<<"ttbar SS before = "<<dijetSS_template[1]->Integral()<<endl;

    dijet_ratio[5] = (TH1D*)f[6]->Get("histDijet2");
    dijet_ratio[5]->SetFillColor(7);
    dijet_ratio[5]->SetStats(kFALSE);
    dijet_ratio[5]->Sumw2();

    dijetSS_ratio[5] = (TH1D*)f[6]->Get("histSameDijet2");
    dijetSS_ratio[1] = (TH1D*)f[1]->Get("histSameDijet2");
    dijetSS_ratio[5]->SetLineColor(1);
    dijetSS_ratio[5]->SetMarkerColor(1);
    dijetSS_ratio[5]->SetMarkerStyle(22);
    dijetSS_ratio[5]->SetMarkerSize(3);
    dijetSS_ratio[5]->SetStats(kFALSE);
    dijetSS_ratio[5]->Sumw2();

    dijet_ratio[0] = (TH1D*)f[0]->Get("histDijet2");
    dijet_ratio[3] = (TH1D*)f[3]->Get("histDijet2");
    dijet_ratio[0]->SetFillColor(2);
    dijet_ratio[0]->SetStats(kFALSE);
    dijet_ratio[0]->Sumw2();
    dijet_ratio[3]->Sumw2();

    dijet_ratio[1] = (TH1D*)f[1]->Get("histDijet2");
    dijet_ratio[1]->SetFillColor(3);
    dijet_ratio[1]->SetStats(kFALSE);
    dijet_ratio[1]->Sumw2();

    cout<<"2"<<endl;

    //DYMuMu, ttbar, WJets, WW, tautau, QCD
    double nEvts[14] = {8.1236e+07, 187625980, 1.65208e+07, 4.52115e+07, 3268361,2.23076e+06,0,0,0,0,0,988416, 999996,985598};
    double xsec[14] = {2008.4*3, 831.76, 60290, 18610, 1915, 2.23076e+06,0,0,0,0,0,63.21,22.82,10.32};
    double norm[14];
    // UPDATED IN 2017
    double lumi = 2759.017;
    // double lumi = 2832.673;

    for(int i=0;i<14;i++) {
        if(i>5&&i<11) continue;
        norm[i] = (xsec[i]*lumi)/nEvts[i];
        cout<<norm[i]<<endl;
    }
    ////
    wjets_template[0]->Scale(norm[0]);
    wjets_template[3]->Scale(norm[3]);
    wjets_template[0]->Add(wjets_template[3]);
    wjets_template[1]->Scale(norm[1]);
    ////

    dijet_template[0]->Scale(norm[0]);
    dijet_template[3]->Scale(norm[3]);
    dijet_template[0]->Add(dijet_template[3]);
    dijet_template[1]->Scale(norm[1]);
    dijet_ratio[0]->Scale(norm[0]);
    dijet_ratio[3]->Scale(norm[3]);
    dijet_ratio[0]->Add(dijet_ratio[3]);
    dijet_ratio[1]->Scale(norm[1]);
    cout<<"DY(template): "<<dijet_template[0]->Integral()<<endl;
    cout<<"DY(ratio): "<<dijet_template[0]->Integral()<<endl;

    setTDRStyle();
    tdrGrid(true);

    lumiTextSize = 0.6;
    cmsTextSize = 1.0;

    dijet_template[5]->GetXaxis()->SetTitle("Mass[GeV]");
    dijet_template[5]->GetYaxis()->SetTitleOffset(1.5);
    dijet_template[5]->GetYaxis()->SetTitle("Number of events");
    //dijet_template[5]->GetXaxis()->SetTitleSize(0.032);
    //dijet_template[5]->GetYaxis()->SetTitleSize(0.032);
    dijet_template[5]->GetXaxis()->SetLabelSize(0.025);
    dijet_template[5]->GetYaxis()->SetLabelSize(0.025);
    dijet_template[5]->GetXaxis()->SetMoreLogLabels(); 

    dijetSS_template[5]->GetXaxis()->SetTitle("Mass[GeV]");
    dijetSS_template[5]->GetYaxis()->SetTitleOffset(1.5);
    dijetSS_template[5]->GetYaxis()->SetTitle("Number of events");
    //dijetSS_template[5]->GetXaxis()->SetTitleSize(0.032);
    //dijetSS_template[5]->GetYaxis()->SetTitleSize(0.032);
    dijetSS_template[5]->GetXaxis()->SetLabelSize(0.025);
    dijetSS_template[5]->GetYaxis()->SetLabelSize(0.025);
    dijetSS_template[5]->GetXaxis()->SetMoreLogLabels();

    dijet_template[0]->GetXaxis()->SetTitle("Mass[GeV]");
    dijet_template[0]->GetYaxis()->SetTitleOffset(1.5);
    dijet_template[0]->GetYaxis()->SetTitle("Number of events");
    //dijet_template[0]->GetXaxis()->SetTitleSize(0.032);
    //dijet_template[0]->GetYaxis()->SetTitleSize(0.032);
    dijet_template[0]->GetXaxis()->SetLabelSize(0.025);
    dijet_template[0]->GetYaxis()->SetLabelSize(0.025);
    dijet_template[0]->GetXaxis()->SetMoreLogLabels(); 

    dijet_template[1]->GetXaxis()->SetTitle("Mass[GeV]");
    dijet_template[1]->GetYaxis()->SetTitleOffset(1.5);
    dijet_template[1]->GetYaxis()->SetTitle("Number of events");
    //dijet_template[1]->GetXaxis()->SetTitleSize(0.032);
    //dijet_template[1]->GetYaxis()->SetTitleSize(0.032);
    dijet_template[1]->GetXaxis()->SetLabelSize(0.025);
    dijet_template[1]->GetYaxis()->SetLabelSize(0.025);
    dijet_template[1]->GetXaxis()->SetMoreLogLabels(); 

    dijet_ratio[5]->GetXaxis()->SetTitle("Mass[GeV]");
    dijet_ratio[5]->GetYaxis()->SetTitleOffset(1.5);
    dijet_ratio[5]->GetYaxis()->SetTitle("Number of events");
    //dijet_ratio[5]->GetXaxis()->SetTitleSize(0.032);
    //dijet_ratio[5]->GetYaxis()->SetTitleSize(0.032);
    dijet_ratio[5]->GetXaxis()->SetLabelSize(0.025);
    dijet_ratio[5]->GetYaxis()->SetLabelSize(0.025);
    dijet_ratio[5]->GetXaxis()->SetMoreLogLabels(); 

    dijetSS_ratio[5]->GetXaxis()->SetTitle("Mass[GeV]");
    dijetSS_ratio[5]->GetYaxis()->SetTitleOffset(1.5);
    dijetSS_ratio[5]->GetYaxis()->SetTitle("Number of events");
    //dijetSS_ratio[5]->GetXaxis()->SetTitleSize(0.032);
    //dijetSS_ratio[5]->GetYaxis()->SetTitleSize(0.032);
    dijetSS_ratio[5]->GetXaxis()->SetLabelSize(0.025);
    dijetSS_ratio[5]->GetYaxis()->SetLabelSize(0.025);
    dijetSS_ratio[5]->GetXaxis()->SetMoreLogLabels();

    dijet_ratio[0]->GetXaxis()->SetTitle("Mass[GeV]");
    dijet_ratio[0]->GetYaxis()->SetTitleOffset(1.5);
    dijet_ratio[0]->GetYaxis()->SetTitle("Number of events");
    //dijet_ratio[0]->GetXaxis()->SetTitleSize(0.032);
    //dijet_ratio[0]->GetYaxis()->SetTitleSize(0.032);
    dijet_ratio[0]->GetXaxis()->SetLabelSize(0.025);
    dijet_ratio[0]->GetYaxis()->SetLabelSize(0.025);
    dijet_ratio[0]->GetXaxis()->SetMoreLogLabels(); 

    dijet_ratio[1]->GetXaxis()->SetTitle("Mass[GeV]");
    dijet_ratio[1]->GetYaxis()->SetTitleOffset(1.5);
    dijet_ratio[1]->GetYaxis()->SetTitle("Number of events");
    dijet_ratio[1]->GetXaxis()->SetLabelSize(0.025);
    dijet_ratio[1]->GetYaxis()->SetLabelSize(0.025);
    dijet_ratio[1]->GetXaxis()->SetMoreLogLabels(); 

    setTDRStyle();
    tdrGrid(true);
    lumiTextSize = 0.5;
    cmsTextSize = 0.75;

    TCanvas* canv = new TCanvas("canv","",1200,1200);
    canv->SetFillColor(0);
    canv->SetLeftMargin( L/W );
    canv->SetRightMargin( R/W );
    canv->SetTopMargin( T/H );
    canv->SetBottomMargin( B/H );

    /////
    for(int i=1; i<46; i++) {
    if(wjets_template[0]->GetBinContent(i) < 0) {
      wjets_template[0]->SetBinContent(i,0.0);
      wjets_template[0]->SetBinError(i,0.0);
    }
    }
    for(int i=1; i<46; i++) {
        if(dijet_template[0]->GetBinContent(i) < 0) {
            dijet_template[0]->SetBinContent(i,0.0);
            dijet_template[0]->SetBinError(i,0.0);
        }
    }
    dijet_template[5]->Add(dijet_template[0],-1.0);
    dijet_template[5]->Add(dijet_template[1],-1.0);

    dijet_ratio[5]->Add(dijet_ratio[0],-1.0);
    dijet_ratio[5]->Add(dijet_ratio[1],-1.0);

    for(int i=1; i<46; i++) {
        if(dijet_template[5]->GetBinContent(i) < 0) {
            dijet_template[5]->SetBinContent(i,0.0);
            dijet_template[5]->SetBinError(i,0.0);
        }
        if(dijet_ratio[5]->GetBinContent(i) < 0) {
            dijet_ratio[5]->SetBinContent(i,0.0);
            dijet_ratio[5]->SetBinError(i,0.0);
        }
    }

    //dijet_template[5]->Smooth();
    TLegend* legg = new TLegend(.6,.65,.95,.89);
    legg->AddEntry(dijet_template[5],"Opposite sign", "F");
    legg->AddEntry(dijetSS_template[5],"Same sign", "P");

    dijet_template[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    canv->Print("print/dijet_template.pdf");
    dijetSS_template[5]->Draw("EPSAME");
    legg->Draw("SAME");
    canv->Print("print/dijetBoth_template.pdf");
    canv->Clear();
    dijetSS_template[5]->SetFillColor(7);
    dijetSS_template[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    canv->Print("print/dijetSS_template.pdf");
    canv->Clear();

    dijet_ratio[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    canv->Print("print/dijet_ratio.pdf");
    dijetSS_ratio[5]->Draw("EPSAME");
    legg->Draw("SAME");
    canv->Print("print/dijetBoth_ratio.pdf");
    canv->Clear();
    dijetSS_ratio[5]->SetFillColor(7);
    dijetSS_ratio[5]->Draw("HIST");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    canv->Print("print/dijetSS_ratio.pdf");
    canv->Clear();

    double error = 0;
    dijet_template[5]->IntegralAndError(1,45,error);
    cout<<"QCD(template) = "<<dijet_template[5]->Integral(1,45)<<"+-"<<error<<endl;
    error = 0;
    dijetSS_template[5]->IntegralAndError(1,45,error);
    cout<<"QCD(template) SS = "<<dijetSS_template[5]->Integral(1,45)<<"+-"<<error<<endl;
    error = 0;
    dijet_ratio[5]->IntegralAndError(1,45,error);
    cout<<"QCD(ratio) = "<<dijet_ratio[5]->Integral(1,45)<<"+-"<<error<<endl;
    error = 0;
    dijetSS_ratio[5]->IntegralAndError(1,45,error);
    cout<<"QCD(ratio) SS = "<<dijetSS_ratio[5]->Integral(1,45)<<"+-"<<error<<endl;

    TH1D* dijet = (TH1D*)dijet_template[5]->Clone();
    dijet->Sumw2();
    dijet->SetName("dijet");

    for(int i=1; i<46; i++) {
        if(dijet_template[5]->GetBinContent(i) < 0) {
          dijet_template[5]->SetBinContent(i,0.0);
          dijet_template[5]->SetBinError(i,0.0);
        }
    }

    TH1D* dijet_total      = new TH1D("dijet_total","",binnum,bins);
    TH1D* dijet_systematic = new TH1D("dijet_systematic","",binnum,bins);
    TH1D* dijet_stat       = new TH1D("dijet_stat","",binnum,bins);

    for(int i=1; i<binnum+1; i++) {

        double systematic = fabs(dijet->GetBinContent(i) - dijet_ratio[5]->GetBinContent(i));
        double stat = dijet->GetBinError(i);
        double total = sqrt( systematic*systematic + stat*stat );
        if(dijet->GetBinContent(i)==0) {
          systematic = 0;
          stat = 0;
          total = 0;
        }

        dijet_systematic->SetBinContent(i,systematic);
        dijet_stat->SetBinContent(i,stat);
        dijet_total->SetBinContent(i,total);

        dijet->SetBinError(i,total);
    }

    TFile* gg = new TFile("result/dijet.root","RECREATE");
    dijet->Write();
    dijet_systematic->Write();
    dijet_stat->Write();
    gg->Close();

    //wjets_systematic->Divide(wjets);
    dijet_systematic->Divide(dijet);
    dijet_stat->Divide(dijet);
    dijet_total->Divide(dijet);

    dijet_total->SetMarkerStyle(20);
    dijet_total->SetMarkerSize(3);
    dijet_total->SetMarkerColor(1);

    dijet_systematic->SetMarkerStyle(22);
    dijet_systematic->SetMarkerSize(3);
    dijet_systematic->SetMarkerColor(2);

    dijet_stat->SetMarkerStyle(21);
    dijet_stat->SetMarkerSize(3);
    dijet_stat->SetMarkerColor(4);

    TLegend* leg = new TLegend(.8,.15,.95,.3);
    leg->AddEntry(dijet_total,"Total","P");
    leg->AddEntry(dijet_systematic,"Sys.","P");
    leg->AddEntry(dijet_stat,"Stat.","P");

    massFrame->GetYaxis()->SetTitle("Unceratinty");
    massFrame->SetMinimum(0);
    massFrame->SetMaximum(1);
    massFrame->GetXaxis()->SetLabelSize(0.025);
    massFrame->GetXaxis()->SetTitleSize(0.05);

    massFrame->Draw();
    dijet_total->Draw("HISTPSAME");
    dijet_systematic->Draw("HISTPSAME");
    dijet_stat->Draw("HISTPSAME");
    massFrame->Draw("AXISSAME");
    leg->Draw("SAME");
    CMS_lumi(canv,4,11);
    canv->Update();
    canv->RedrawAxis();
    canv->GetFrame()->Draw();
    canv->SetLogx();
    canv->Print("print/dijet_uncertainty.pdf");
}
