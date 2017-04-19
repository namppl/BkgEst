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
#include "../../interface/analysis.h"

using namespace std;

void selectDenAndNumForFR(int index)
{
    // Event & muon
    cout<<"Event"<<endl;
    PhysicsEvent* event = new PhysicsEvent();
    vector<PhysicsMuon>* passingMuons = new vector<PhysicsMuon>;

    // Histograms  
    const int ptbinnum_endcap = 9;
    double ptbin_endcap[ptbinnum_endcap+1] = {47,52,60,70,80,90,100,150,200,500};
    const int ptbinnum = 17;
    double ptbin[ptbinnum+1] = {47,52,60,70,80,90,100,120,140,160,180,200,250,300,350,400,450,500};

    TH1D* denominator_pt = new TH1D("denominator_pt","",ptbinnum,ptbin); 
    TH1D* denominator_pt_barrel = new TH1D("denominator_pt_barrel","",ptbinnum,ptbin); 
    TH1D* denominator_pt_endcap = new TH1D("denominator_pt_endcap","",ptbinnum_endcap,ptbin_endcap);  

    TH1D* numerator_pt = new TH1D("numerator_pt","",ptbinnum,ptbin); 
    TH1D* numerator_pt_barrel = new TH1D("numerator_pt_barrel","",ptbinnum,ptbin); 
    TH1D* numerator_pt_endcap = new TH1D("numerator_pt_endcap","",ptbinnum_endcap,ptbin_endcap); 

    TH1D* denominator_eta = new TH1D("denominator_eta","",48,-2.4,2.4); 
    TH1D* numerator_eta = new TH1D("numerator_eta","",48,-2.4,2.4); 

    TH1D* denominator = new TH1D("denominator","",100,0,5);
    TH1D* denominator_barrel = new TH1D("denominator_barrel","",100,0,5);
    TH1D* denominator_endcap = new TH1D("denominator_endcap","",100,0,5);

    TH1D* numerator = new TH1D("numerator","",100,0,5);
    TH1D* numerator_barrel = new TH1D("numerator_barrel","",100,0,5);
    TH1D* numerator_endcap = new TH1D("numerator_endcap","",100,0,5);

    denominator_pt->Sumw2();
    denominator_pt_barrel->Sumw2();
    denominator_pt_endcap->Sumw2();

    numerator_pt->Sumw2();
    numerator_pt_barrel->Sumw2();
    numerator_pt_endcap->Sumw2();

    denominator_eta->Sumw2();
    numerator_eta->Sumw2();

    denominator->Sumw2();
    denominator_barrel->Sumw2();
    denominator_endcap->Sumw2();

    numerator->Sumw2();
    numerator_barrel->Sumw2();
    numerator_endcap->Sumw2();

    // Chain
    cout<<"Chain"<<endl;
    bool mc = false;
    TChain* chain = new TChain("tree/physicsTree");    
    if(index==-1) {
        chain->Add("/data4/Users/knam/13TeV/76X/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns/MuonPhys/Run2015C_25ns-16Dec2015-v1/*.root");
        chain->Add("/data4/Users/knam/13TeV/76X/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns/MuonPhys/Run2015D-16Dec2015-v1/*.root");
    }
    else {
        mc = true;
        if(index==0) chain->Add("/data4/Users/knam/13TeV/76X/DYJetsToLL_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/M-50/*.root");
        else if(index==1) chain->Add("/data4/Users/knam/13TeV/76X/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/*.root");
        else if(index==2) chain->Add("/data4/Users/knam/13TeV/76X_v1/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/*.root");
        else if(index==11) chain->Add("/data4/Users/knam/13TeV/76X/WW_TuneCUETP8M1_13TeV-pythia8/*.root");
        else if(index==12) chain->Add("/data4/Users/knam/13TeV/76X/WZ_TuneCUETP8M1_13TeV-pythia8/*.root");
        else if(index==13) chain->Add("/data4/Users/knam/13TeV/76X/ZZ_TuneCUETP8M1_13TeV-pythia8/*.root");
        else if(index==14) chain->Add("/data4/Users/knam/13TeV/76X/ST_tW_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/top/*.root");
        else if(index==15) chain->Add("/data4/Users/knam/13TeV/76X/ST_tW_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/antitop/*.root");
        else {
            cout<<"Wrong input"<<endl;
            return;
        }
    }

    chain->SetBranchAddress("event",&event);  
    cout<<"# of events = "<<chain->GetEntries()<<endl;

    double pt = 0;
    double eta = 0;
    double iso = 0;
    double wt = 1.0;
    double wtsum = 0;
	//int entries = 10000;
    int entries = chain->GetEntries();
    
    for(int i=0; i!=entries; i++) {

        chain->GetEntry(i);
        if(mc) {
            wt = event->weight;
            wtsum += wt;
        }

        if( !event->TriggerSelection("HLT_Mu50_v*") && !event->TriggerSelection("HLT_Mu45_eta2p1_v*") ) continue;

        passingMuons->clear();

        for(unsigned j=0; j!=event->muons.size(); j++) {

            PhysicsMuon* mu_ = (PhysicsMuon*)&event->muons.at(j);

            if( mu_->acceptance(47,2.4) ) {
                passingMuons->push_back(*mu_);
            }
        }

        for(unsigned j=0; j!=passingMuons->size(); j++) {

            PhysicsMuon mu = passingMuons->at(j);
            pt = mu.pt;
            eta = mu.eta;
            iso = mu.isolationR03_sumpt/mu.pt;

            if( !mu.highPtMuonID() ) continue;

            denominator_pt->Fill(pt,wt);
            denominator_eta->Fill(eta,wt);
            denominator->Fill(iso,wt); 

            if( fabs(eta)<1.2 ) {
                denominator_barrel->Fill(iso,wt); 
                denominator_pt_barrel->Fill(pt,wt);
            } 
            else {
                denominator_endcap->Fill(iso,wt); 
                denominator_pt_endcap->Fill(pt,wt);
            } 

            if( !mu.isolation(0.1) ) continue; 

            numerator_pt->Fill(pt,wt);
            numerator_eta->Fill(eta,wt);
            numerator->Fill(iso,wt); 

            if( fabs(eta)<1.2 ) {
                numerator_barrel->Fill(iso,wt); 
                numerator_pt_barrel->Fill(pt,wt);
            } 
            else {
                numerator_endcap->Fill(iso,wt); 
                numerator_pt_endcap->Fill(pt,wt); 
            }
        }
    }

    if(index==-1) index=6;
    TFile* f = new TFile("histograms/hist"+TString::Itoa(index,10)+".root","RECREATE");

    denominator_pt->Write();
    denominator_pt_barrel->Write();
    denominator_pt_endcap->Write();

    numerator_pt->Write();
    numerator_pt_barrel->Write();
    numerator_pt_endcap->Write();

    denominator_eta->Write();
    numerator_eta->Write();

    denominator->Write();
    denominator_barrel->Write();
    denominator_endcap->Write();

    numerator->Write();
    numerator_barrel->Write();
    numerator_endcap->Write();

    f->Close();

    cout<<"Success"<<endl;
    cout<<wtsum<<endl;
}
