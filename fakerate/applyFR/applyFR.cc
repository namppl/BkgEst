#include <TStyle.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
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
#include <TTimeStamp.h>
#include <TColor.h>
#include <TLatex.h>
#include <TEfficiency.h>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "../interface/analysis.h"
using namespace std;

void applyFR(int index) {

    cout<<"Chain"<<endl;
    bool mc = false;
    PhysicsEvent* event = new PhysicsEvent();
    TChain* chain = new TChain("tree/physicsTree");
    if(index==-1) {
        chain->Add("/data4/Users/knam/13TeV/76X_38T/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns/MuonPhys/Run2015C_25ns-16Dec2015-v1/*.root");
        chain->Add("/data4/Users/knam/13TeV/76X_38T/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns/MuonPhys/Run2015D-16Dec2015-v1/*.root");
    }
    else {
        mc = true;
        if(index==0) chain->Add("/data4/Users/knam/13TeV/76X_v1/DYJetsToLL_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/M-50_ext4/*.root");
        else if(index==100) chain->Add("/data4/Users/knam/13TeV/76X_v1/DYJetsToLL_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/M-10to50/*.root");
        else if(index==1) chain->Add("/data4/Users/knam/13TeV/76X_v1/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/*.root");
        else if(index==2) chain->Add("/data4/Users/knam/13TeV/76X_v1/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/*.root");
        else if(index==11) chain->Add("/data4/Users/knam/13TeV/76X_v1/WW_TuneCUETP8M1_13TeV-pythia8/*.root");
        else if(index==12) chain->Add("/data4/Users/knam/13TeV/76X_v1/WZ_TuneCUETP8M1_13TeV-pythia8/*.root");
        else if(index==13) chain->Add("/data4/Users/knam/13TeV/76X_v1/ZZ_TuneCUETP8M1_13TeV-pythia8/*.root");
        else if(index==14) chain->Add("/data4/Users/knam/13TeV/76X_v1/ST_tW_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/top/*.root");
        else if(index==15) chain->Add("/data4/Users/knam/13TeV/76X_v1/ST_tW_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/antitop/*.root");
        else {
            cout<<"Wrong input"<<endl;
            return;
        }
    }
    chain->SetBranchAddress("event",&event);  
    cout<<"# of events = "<<chain->GetEntries()<<endl;

    vector<pair<PhysicsMuon,int>>* passingMuons = new vector<pair<PhysicsMuon,int>>;
    vector<pair<PhysicsMuon,int>>* failingMuons = new vector<pair<PhysicsMuon,int>>;
    pair<PhysicsMuon,PhysicsMuon>* tempMuons = new pair<PhysicsMuon,PhysicsMuon>;

    int binnum = 43;
    double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1500, 3000};

    TH1D* histDijet1 = new TH1D("histDijet1","",binnum,bins);
    TH1D* histDijet2 = new TH1D("histDijet2","",binnum,bins);
    TH1D* histSameDijet1 = new TH1D("histSameDijet1","",binnum,bins);
    TH1D* histSameDijet2 = new TH1D("histSameDijet2","",binnum,bins);

    TH1D* fitDijet1 = new TH1D("fitDijet1","",37,15,200);
    TH1D* fitDijet2 = new TH1D("fitDijet2","",37,15,200);
    TH1D* fitSameDijet1 = new TH1D("fitSameDijet1","",37,15,200);
    TH1D* fitSameDijet2 = new TH1D("fitSameDijet2","",37,15,200);

    TH1D* rapDijet1 = new TH1D("rapDijet1","",48,-2.4,2.4);
    TH1D* rapDijet2 = new TH1D("rapDijet2","",48,-2.4,2.4);
    TH1D* rapSameDijet1 = new TH1D("rapSameDijet1","",48,-2.4,2.4);
    TH1D* rapSameDijet2 = new TH1D("rapSameDijet2","",48,-2.4,2.4);

    histDijet1->Sumw2();
    histDijet2->Sumw2();
    histSameDijet1->Sumw2();
    histSameDijet2->Sumw2();

    fitDijet1->Sumw2();
    fitDijet2->Sumw2();
    fitSameDijet1->Sumw2();
    fitSameDijet2->Sumw2();

    rapDijet1->Sumw2();
    rapDijet2->Sumw2();
    rapSameDijet1->Sumw2();
    rapSameDijet2->Sumw2();

    TH1D* histWJets1 = new TH1D("histWJets1","",binnum,bins);
    TH1D* histWJets2 = new TH1D("histWJets2","",binnum,bins);
    TH1D* histSameWJets1 = new TH1D("histSameWJets1","",binnum,bins);
    TH1D* histSameWJets2 = new TH1D("histSameWJets2","",binnum,bins);

    TH1D* fitWJets1 = new TH1D("fitWJets1","",37,15,200);
    TH1D* fitWJets2 = new TH1D("fitWJets2","",37,15,200);
    TH1D* fitSameWJets1 = new TH1D("fitSameWJets1","",37,15,200);
    TH1D* fitSameWJets2 = new TH1D("fitSameWJets2","",37,15,200);

    TH1D* rapWJets1 = new TH1D("rapWJets1","",48,-2.4,2.4);
    TH1D* rapWJets2 = new TH1D("rapWJets2","",48,-2.4,2.4);
    TH1D* rapSameWJets1 = new TH1D("rapSameWJets1","",48,-2.4,2.4);
    TH1D* rapSameWJets2 = new TH1D("rapSameWJets2","",48,-2.4,2.4);

    TH1D* histPass = new TH1D("histPass","",10,0,10);
    TH1D* histFail = new TH1D("histFail","",10,0,10);

    histWJets1->Sumw2();
    histWJets2->Sumw2();
    histSameWJets1->Sumw2();
    histSameWJets2->Sumw2();

    fitWJets1->Sumw2();
    fitWJets2->Sumw2();
    fitSameWJets1->Sumw2();
    fitSameWJets2->Sumw2();

    rapWJets1->Sumw2();
    rapWJets2->Sumw2();
    rapSameWJets1->Sumw2();
    rapSameWJets2->Sumw2();

    double pt = 0;
    double eta = 0;
    double wt = 1.0;
    double wtsum = 0;
    bool leading = false;

    double FR1_template;
    double FR2_template;
    double FR1_ratio;
    double FR2_ratio;
    double weight_template;
    double weight_ratio;
    double mass;
    double sign;
    double rap;

    int nPass = 0;
    int nFail = 0;

    int nEntries = chain->GetEntries();
    for(int i=0; i!=nEntries; i++) {
        chain->GetEntry(i);
        if(mc) wt = event->weight;
        else wt = 1.0;
        wtsum += wt;

        leading = false; 
        passingMuons->clear();
        failingMuons->clear();

        for(unsigned j=0; j!=event->muons.size(); j++) {
            PhysicsMuon* mu_ = (PhysicsMuon*)&event->muons.at(j);
            if( mu_->highPtMuonID() && mu_->acceptance(10,2.4) ) {
                if( mu_->pt > 22 ) leading = true;

                if( mu_->isolation(0.1) ) passingMuons->push_back({*mu_,j});
                else failingMuons->push_back({*mu_,j});        
            }
        }

        if( !leading ) continue;

        histPass->Fill(passingMuons->size(),wt);
        histFail->Fill(failingMuons->size(),wt);

        nPass += wt*passingMuons->size();
        nFail += wt*failingMuons->size();

        if( failingMuons->size()>1 ) {
            if( dijetDY(event->dimuons,failingMuons,tempMuons) ) {
                FR1_template = FR_template(tempMuons->first);
                FR2_template = FR_template(tempMuons->second);
                FR1_ratio = FR_ratio(tempMuons->first);
                FR2_ratio = FR_ratio(tempMuons->second);

                weight_template = wt*FR1_template*FR2_template/((1-FR1_template)*(1-FR2_template));
                weight_ratio = wt*FR1_ratio*FR2_ratio/((1-FR1_ratio)*(1-FR2_ratio));
				/*
                if(weight_template > 1 || weight_ratio > 1 ) {
                    cout<<"wt_temp = "<<weight_template<<", wt_ratio= "<<weight_ratio<<endl;
                    cout<<"fr_temp1 = "<<FR1_template<<", fr_ratio1= "<<FR1_ratio<<endl;
                    cout<<"fr_temp2 = "<<FR2_template<<", fr_ratio2= "<<FR2_ratio<<endl;
                    cout<<"pt = "<<tempMuons->first.pt<<", "<<tempMuons->second.pt<<endl;
                    cout<<"eta = "<<tempMuons->first.eta<<", "<<tempMuons->second.eta<<endl;
                }*/

                mass = (tempMuons->first.momentum() + tempMuons->second.momentum()).M();
                sign = tempMuons->first.charge * tempMuons->second.charge;
                rap  = (tempMuons->first.momentum() + tempMuons->second.momentum()).Rapidity();

                if( sign < 0 ) {
                    if( mass > 15 && mass < 3000) {
                    	histDijet1->Fill(mass, weight_template);
                    	histDijet2->Fill(mass, weight_ratio);
                    	fitDijet1->Fill(mass, weight_template);
                    	fitDijet2->Fill(mass, weight_ratio);
                        rapDijet1->Fill(rap, weight_template);
                        rapDijet2->Fill(rap, weight_ratio);
                    }
                }
                else {
                    if( mass > 15 && mass < 3000) {
                    	histSameDijet1->Fill(mass, weight_template);
                    	histSameDijet2->Fill(mass, weight_ratio);
                    	fitSameDijet1->Fill(mass, weight_template);
                    	fitSameDijet2->Fill(mass, weight_ratio);
                        rapSameDijet1->Fill(rap, weight_template);
                        rapSameDijet2->Fill(rap, weight_ratio);
                    }
                }
            }
        }
        else if( failingMuons->size()==1 && passingMuons->size()==1 ) {
            if( wjetsDY(event->dimuons, failingMuons->at(0), passingMuons->at(0), tempMuons) ) {
                FR1_template = FR_template(tempMuons->first);
                FR1_ratio = FR_ratio(tempMuons->first);

                weight_template = wt*FR1_template/(1-FR1_template);
                weight_ratio = wt*FR1_ratio/(1-FR1_ratio);
                mass = (tempMuons->first.momentum() + tempMuons->second.momentum()).M();
                sign = tempMuons->first.charge * tempMuons->second.charge;
                rap  = (tempMuons->first.momentum() + tempMuons->second.momentum()).Rapidity();

                if( sign < 0 ) {
                    if( mass > 15 && mass < 3000) {
                    	histWJets1->Fill(mass, weight_template);
                    	histWJets2->Fill(mass, weight_ratio);
                    	fitWJets1->Fill(mass, weight_template);
                    	fitWJets2->Fill(mass, weight_ratio);
                        rapWJets1->Fill(rap, weight_template);
                        rapWJets2->Fill(rap, weight_ratio);
                    }
                }
                else {
                    if( mass > 15 && mass < 3000) {
                    	histSameWJets1->Fill(mass, weight_template);
                    	histSameWJets2->Fill(mass, weight_ratio);
                    	fitSameWJets1->Fill(mass, weight_template);
                    	fitSameWJets2->Fill(mass, weight_ratio);
                        rapSameWJets1->Fill(rap, weight_template);
                        rapSameWJets2->Fill(rap, weight_ratio);
                    }
                }
            }
        } 
    }

    if(index==-1) index=6;
    TFile* f = new TFile("histograms/fake"+TString::Itoa(index,10)+".root","RECREATE");

    histDijet1->Write();
    histDijet2->Write();
    histWJets1->Write();
    histWJets2->Write();
    fitDijet1->Write();
    fitDijet2->Write();
    fitWJets1->Write();
    fitWJets2->Write();
    rapDijet1->Write();
    rapDijet2->Write();
    rapWJets1->Write();
    rapWJets2->Write();

    histSameDijet1->Write();
    histSameDijet2->Write();
    histSameWJets1->Write();
    histSameWJets2->Write();
    fitSameDijet1->Write();
    fitSameDijet2->Write();
    fitSameWJets1->Write();
    fitSameWJets2->Write();
    rapSameDijet1->Write();
    rapSameDijet2->Write();
    rapSameWJets1->Write();
    rapSameWJets2->Write();

    histPass->Write();
    histFail->Write();

    f->Close();

    cout<<"# of passing muons = "<<nPass<<endl;
    cout<<"# of failing muons = "<<nFail<<endl;
    cout<<endl;
    cout<<"# of passing muons per event = "<<nPass/wtsum<<endl;
    cout<<"# of failing muons per event = "<<nFail/wtsum<<endl;
    cout<<endl;

    cout<<"Success"<<endl;
    cout<<wtsum<<endl;
}
