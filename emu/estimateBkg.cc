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

#include "../interface/tdrstyle.C"
#include "../interface/CMS_lumi.C"

using namespace std;

void fillSystematics( TH1D* data_driven, TH1D* stat, TH1D* systematic, TH1D* total );
void removeNegativeBins( TH1D* hist );

void estimateBkg() {

    const int binnum = 43;
    double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1500, 3000};

    TFile* file_data = new TFile("histograms/hist6.root","READ");
    TFile* file[16];
    file[0] = new TFile("histograms/hist0.root","READ");
    file[3] = new TFile("histograms/hist100.root","READ");
    file[1] = new TFile("histograms/hist999.root","READ");
    file[11] = new TFile("histograms/hist11.root","READ");
    file[12] = new TFile("histograms/hist12.root","READ");
    file[13] = new TFile("histograms/hist13.root","READ");
    file[14] = new TFile("histograms/hist14.root","READ");
    file[15] = new TFile("histograms/hist15.root","READ");

    //DYMuMu, ttbar, WJets, WW, tautau, QCD
    double nEvts[16] = {8.1236e+07, 187625980, 1.65412e+07, /*5.57133e+07*/4.52115e+07, 3268361,2.23076e+06,0,0,0,0,0,988416,999996,985598,999999,999399}; 
    double xsec[16] = {2008.4*3, 831.76, 18610, 54.8, 1915, 2.23076e+06,0,0,0,0,0,118.7,47.13,16.523,38.09,38.09};
    double norm_Silver[16];
    double norm_MuPhys[16];

    // UPDATED IN 2017
    double lumi_Silver = 2620.675;
    double lumi_MuPhys = 2759.017;

    TH1D* emu[16];
    TH1D* emuSS[16];
    TH1D* dimu[16];

    for(int i=0;i<16;i++) {
        if( i==2 || (i>3 && i<11) ) continue;
        norm_Silver[i] = (xsec[i]*lumi_Silver)/nEvts[i];
        norm_MuPhys[i] = (xsec[i]*lumi_MuPhys)/nEvts[i];

        emu[i] = (TH1D*)file[i]->Get("emu_mass")->Clone("emu"+TString::Itoa(i,10));
        emuSS[i] = (TH1D*)file[i]->Get("emuSS_mass")->Clone("emuSS"+TString::Itoa(i,10));
        dimu[i] = (TH1D*)file[i]->Get("dimu_mass")->Clone("dimu"+TString::Itoa(i,10));

        emu[i]->Scale(norm_Silver[i]);
        emuSS[i]->Scale(norm_Silver[i]);
        dimu[i]->Scale(norm_MuPhys[i]);

        emu[i]->SetFillColor(i+2);
        emuSS[i]->SetFillColor(i+2);
        dimu[i]->SetFillColor(i+2);

        emu[i]->SetStats(kFALSE);
        emuSS[i]->SetStats(kFALSE);
        dimu[i]->SetStats(kFALSE);
    } 

    TH1D* emu_data = (TH1D*)file_data->Get("emu_mass")->Clone("emu_data");
    TH1D* emu_ttbar = emu[1];
    TH1D* emu_DYtautau = emu[0];
    TH1D* emu_DYtautau_M10to50 = emu[3];
    TH1D* emu_WW = emu[11];
    TH1D* emu_WZ = emu[12];
    TH1D* emu_ZZ = emu[13];
    TH1D* emu_tW = emu[14];
    TH1D* emu_antitW = emu[15];

    TH1D* emuSS_data = (TH1D*)file_data->Get("emuSS_mass")->Clone("emuSS_data");
    TH1D* emuSS_ttbar = emuSS[1];
    TH1D* emuSS_DYtautau = emuSS[0];
    TH1D* emuSS_DYtautau_M10to50 = emuSS[3];
    TH1D* emuSS_WW = emuSS[11];
    TH1D* emuSS_WZ = emuSS[12];
    TH1D* emuSS_ZZ = emuSS[13];
    TH1D* emuSS_tW = emuSS[14];
    TH1D* emuSS_antitW = emuSS[15];

    TH1D* dimu_ttbar = dimu[1];
    TH1D* dimu_DYtautau = dimu[0];
    TH1D* dimu_DYtautau_M10to50 = dimu[3];
    TH1D* dimu_WW = dimu[11];
    TH1D* dimu_WZ = dimu[12];
    TH1D* dimu_ZZ = dimu[13];
    TH1D* dimu_tW = dimu[14];
    TH1D* dimu_antitW = dimu[15];

    emu_data->SetMarkerStyle(33);
    emu_data->SetMarkerSize(3);
    emu_data->SetStats(kFALSE);

    emuSS_data->SetMarkerStyle(33);
    emuSS_data->SetMarkerSize(3);
    emuSS_data->SetStats(kFALSE);

    // DY M50 + M10to50
    emu_DYtautau->Add(emu_DYtautau_M10to50);
    emuSS_DYtautau->Add(emuSS_DYtautau_M10to50);
    dimu_DYtautau->Add(dimu_DYtautau_M10to50);

    // WW + WZ + ZZ
    TH1D* emu_diboson = emu_WW;
    TH1D* emuSS_diboson = emuSS_WW;
    emu_diboson->Add(emu_WZ);
    emu_diboson->Add(emu_ZZ);
    emuSS_diboson->Add(emuSS_WZ);
    emuSS_diboson->Add(emuSS_ZZ);
    // dimu_WW->Add(dimu_WZ);
    // dimu_WW->Add(dimu_ZZ);

    // tW + antitW
    emu_tW->Add(emu_antitW);
    emuSS_tW->Add(emuSS_antitW);
    dimu_tW->Add(dimu_antitW);

    cout<<"data in emu: "<<emu_data->Integral(1,43)<<endl;
    cout<<"ttbar in emu: "<<emu_ttbar->Integral(1,43)<<endl;
    cout<<"DY in emu: "<<emu_DYtautau->Integral(1,43)<<endl;
    cout<<"diboson in emu: "<<emu_diboson->Integral(1,43)<<endl;
    cout<<"tW in emu: "<<emu_tW->Integral(1,43)<<endl;
    cout<<"total MC in emu: "<<emu_DYtautau->Integral(1,43)+emu_ttbar->Integral(1,43)+emu_diboson->Integral(1,43)+emu_tW->Integral(1,43)<<endl;
    cout<<endl;
    cout<<"data in emuSS: "<<emuSS_data->Integral(1,43)<<endl;
    cout<<"ttbar in emuSS: "<<emuSS_ttbar->Integral(1,43)<<endl;
    cout<<"DY in emuSS: "<<emuSS_DYtautau->Integral(1,43)<<endl;
    cout<<"diboson in emuSS: "<<emuSS_diboson->Integral(1,43)<<endl;
    cout<<"tW in emuSS: "<<emuSS_tW->Integral(1,43)<<endl;
    cout<<endl;
    cout<<"ttbar in dimu: "<<dimu_ttbar->Integral(1,43)<<endl;
    cout<<"DY in dimu: "<<dimu_DYtautau->Integral(1,43)<<endl;
    cout<<"diboson in dimu: "<<dimu_WW->Integral(1,43)<<endl;
    cout<<"tW in dimu: "<<dimu_tW->Integral(1,43)<<endl;
    cout<<endl;

    // remove negative bins
    removeNegativeBins( emu_DYtautau );
    removeNegativeBins( emuSS_DYtautau );
    removeNegativeBins( dimu_DYtautau );

    TH1D* emu_sumBkg = new TH1D("emu_sumBkg","",binnum,bins);
    THStack* emu_stackBkg = new THStack("emu_stackBkg","");

    TH1D* emuSS_sumMC = new TH1D("emuSS_sumMC","",binnum,bins);
    THStack* emuSS_stackMC = new THStack("emuSS_stackMC","");

    emu_sumBkg->Add(emu_DYtautau);	
    emu_sumBkg->Add(emu_ttbar);	
    emu_sumBkg->Add(emu_diboson);	
    emu_sumBkg->Add(emu_tW);	

    emu_stackBkg->Add(emu_diboson);
    emu_stackBkg->Add(emu_DYtautau);
    emu_stackBkg->Add(emu_tW);
    emu_stackBkg->Add(emu_ttbar);

    emuSS_sumMC->Add(emuSS_DYtautau);  
    emuSS_sumMC->Add(emuSS_ttbar); 
    emuSS_sumMC->Add(emuSS_diboson); 
    emuSS_sumMC->Add(emuSS_tW);  

    emuSS_stackMC->Add(emuSS_diboson);
    emuSS_stackMC->Add(emuSS_DYtautau);
    emuSS_stackMC->Add(emuSS_tW);
    emuSS_stackMC->Add(emuSS_ttbar);

    // Legend
    TLegend* legend = new TLegend(.65,.55,.95,.89);
    legend->AddEntry(emu_data,"Data");
    legend->AddEntry(emu_ttbar,"ttbar","F");
    legend->AddEntry(emu_tW,"tW+#bar{t}W","F");
    legend->AddEntry(emu_DYtautau,"DY","F");
    legend->AddEntry(emu_diboson,"VV","F");
    legend->SetBorderSize(0);  

    TH1D* emu_QCD = (TH1D*)emuSS_data->Clone();
    emu_QCD->Add(emuSS_DYtautau,-1.0);
    emu_QCD->Add(emuSS_ttbar,-1.0);
    emu_QCD->Add(emuSS_diboson,-1.0);
    emu_QCD->Add(emuSS_tW,-1.0);
    emu_QCD->SetFillColor(7);

    const double RR = 0.57147108645;
    emu_QCD->Scale(1/RR);

    removeNegativeBins(emu_QCD);

    emu_sumBkg->Add(emu_QCD);
    emu_stackBkg->Add(emu_QCD);
    legend->AddEntry(emu_QCD,"QCD","F");

    TH1D* emu_ratio = (TH1D*)emu_data->Clone("emu_ratio");
    emu_ratio->Divide(emu_data,emu_sumBkg,1.0,1.0,"B");

    TH1D* data_driven_ttbar = (TH1D*)emu_ratio->Clone("data_driven_ttbar");
    TH1D* data_driven_tW = (TH1D*)emu_ratio->Clone("data_driven_tW");
    TH1D* data_driven_WW = (TH1D*)emu_ratio->Clone("data_driven_WW");
    TH1D* data_driven_DYtautau = (TH1D*)emu_ratio->Clone("data_driven_DYtautau");

    removeNegativeBins(data_driven_ttbar);
    removeNegativeBins(data_driven_tW);
    removeNegativeBins(data_driven_WW);
    removeNegativeBins(data_driven_DYtautau);

    data_driven_ttbar->Multiply(dimu_ttbar);
    data_driven_DYtautau->Multiply(dimu_DYtautau);
    data_driven_WW->Multiply(dimu_WW);
    data_driven_tW->Multiply(dimu_tW);

    // replace the highest mass bins with MC (becasue data-driven entries are zero)
    data_driven_ttbar->SetBinContent(43,dimu_ttbar->GetBinContent(43));
    data_driven_ttbar->SetBinError(43,dimu_ttbar->GetBinError(43)); 

    TH1D* ttbar_total = new TH1D("ttbar_total","",binnum,bins);
    TH1D* ttbar_systematic = new TH1D("ttbar_systematic","",binnum,bins);
    TH1D* ttbar_stat = new TH1D("ttbar_stat","",binnum,bins);

    TH1D* tW_total = new TH1D("tW_total","",binnum,bins);
    TH1D* tW_systematic = new TH1D("tW_systematic","",binnum,bins);
    TH1D* tW_stat = new TH1D("tW_stat","",binnum,bins);

    TH1D* WW_total = new TH1D("WW_total","",binnum,bins);
    TH1D* WW_systematic = new TH1D("WW_systematic","",binnum,bins);
    TH1D* WW_stat = new TH1D("WW_stat","",binnum,bins);

    TH1D* DYtautau_total = new TH1D("DYtautau_total","",binnum,bins);
    TH1D* DYtautau_systematic = new TH1D("DYtautau_systematic","",binnum,bins);
    TH1D* DYtautau_stat = new TH1D("DYtautau_stat","",binnum,bins);

    ttbar_systematic->Add(dimu_ttbar);
    ttbar_systematic->Add(data_driven_ttbar,-1.0);

    DYtautau_systematic->Add(dimu_DYtautau);
    DYtautau_systematic->Add(data_driven_DYtautau,-1.0);

    tW_systematic->Add(dimu_tW);
    tW_systematic->Add(data_driven_tW,-1.0);

    WW_systematic->Add(dimu_WW);
    WW_systematic->Add(data_driven_WW,-1.0);

    fillSystematics( data_driven_ttbar, ttbar_stat, ttbar_systematic, ttbar_total );
    fillSystematics( data_driven_DYtautau, DYtautau_stat, DYtautau_systematic, DYtautau_total );
    fillSystematics( data_driven_tW, tW_stat, tW_systematic, tW_total );
    fillSystematics( data_driven_WW, WW_stat, WW_systematic, WW_total );

    data_driven_ttbar->SetName("ttbar");  
    data_driven_DYtautau->SetName("DYtautau");
    data_driven_tW->SetName("tW");
    data_driven_WW->SetName("WW");

    dimu_ttbar->SetName("ttbar_MC");  
    dimu_DYtautau->SetName("DYtautau_MC");
    dimu_tW->SetName("tW_MC");
    dimu_WW->SetName("WW_MC");

    TFile* g = new TFile("result/emu.root","RECREATE");

    data_driven_ttbar->Write();
    data_driven_DYtautau->Write();
    data_driven_tW->Write();
    data_driven_WW->Write();

    dimu_ttbar->Write();
    dimu_DYtautau->Write();
    dimu_tW->Write();
    dimu_WW->Write();

    ttbar_systematic->Write();
    DYtautau_systematic->Write();
    tW_systematic->Write();
    WW_systematic->Write();

    ttbar_stat->Write();
    DYtautau_stat->Write();
    tW_stat->Write();
    WW_stat->Write();

    g->Close();

}

void fillSystematics( TH1D* data_driven, TH1D* stat, TH1D* systematic, TH1D* total ) {

    double binSystematic = 0;
    double binStat = 0;
    double binTotal = 0;

    for(int i=0; i<data_driven->GetNbinsX(); i++) {

        if(data_driven->GetBinContent(i+1)!=0) {  
            binStat = data_driven->GetBinError(i+1);
            binSystematic = fabs(systematic->GetBinContent(i+1));
            binTotal = sqrt(binSystematic*binSystematic + binStat*binStat);
        }
        else{
            binSystematic = 0;
            binStat = 0;
            binTotal = 0;
        }

        systematic->SetBinContent(i+1,binSystematic);
        stat->SetBinContent(i+1,binStat);    
        total->SetBinContent(i+1,binTotal);

        data_driven->SetBinError(i+1,binTotal);

    }

}

void removeNegativeBins( TH1D* hist ) {

    for(int i=0; i<hist->GetNbinsX(); i++) {
        if(hist->GetBinContent(i+1)<0) {
            hist->SetBinContent(i+1,0);
            hist->SetBinError(i+1,0);
        }
    }   

}
