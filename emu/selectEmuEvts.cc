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

void selectEmuEvts(int index)
{
    //Event
    TChain*       chain = new TChain("tree/physicsTree");
    PhysicsEvent* event = new PhysicsEvent();
    bool isData = false;
    if(index==6) {
        isData = true;
        chain->Add("/data4/Users/knam/13TeV/76X_v1/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns/Silver/Run2015C_25ns-16Dec2015-v1/*.root");
        chain->Add("/data4/Users/knam/13TeV/76X_v1/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns/Silver/Run2015D-16Dec2015-v1/*.root");
    }
    else if(index==100) chain->Add("/data4/Users/knam/13TeV/76X_zip/DYJetsToLL_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/M-10to50/*.root");
    else if(index==0) chain->Add("/data4/Users/knam/13TeV/76X_v1/DYJetsToLL_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/M-50_ext4/*.root");
    else if(index==1) chain->Add("/data4/Users/knam/13TeV/76X_v1/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8/*.root");
    else if(index==111111) chain->Add("/data4/Users/knam/13TeV/76X_v1/TT_TuneCUETP8M1_13TeV-amcatnlo-pythia8_another/*.root");
    else if(index==11) chain->Add("/data4/Users/knam/13TeV/76X_v1/WW_TuneCUETP8M1_13TeV-pythia8/*.root");
    else if(index==12) chain->Add("/data4/Users/knam/13TeV/76X_v1/WZ_TuneCUETP8M1_13TeV-pythia8/*.root");
    else if(index==13) chain->Add("/data4/Users/knam/13TeV/76X_v1/ZZ_TuneCUETP8M1_13TeV-pythia8/*.root");
    else if(index==14) chain->Add("/data4/Users/knam/13TeV/76X_v1/ST_tW_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/top/*.root");
    else if(index==15) chain->Add("/data4/Users/knam/13TeV/76X_v1/ST_tW_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/antitop/*.root");
    else {
        cout<<"wrong input"<<endl;
        return;
    }
    chain->SetBranchAddress("event",&event);  

    vector<pair<PhysicsMuon,int>>*     passingMuons     = new vector<pair<PhysicsMuon,int>>;
    vector<pair<PhysicsElectron,int>>* passingElectrons = new vector<pair<PhysicsElectron,int>>;
    pair<PhysicsMuon,PhysicsMuon>*     dimuon           = new pair<PhysicsMuon,PhysicsMuon>;
    pair<PhysicsElectron,PhysicsMuon>* emu              = new pair<PhysicsElectron,PhysicsMuon>;

    //Histogram 
    const int binsize = 43;
    const double bins[44] = {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 150, 160, 171, 185,  200, 220, 243, 273, 320, 380, 440, 510, 600, 700,  830, 1000, 1500, 3000};

    TH1D* el_etSC     = new TH1D("el_etSC","",97,30,1000); 
    TH1D* el_etaSC    = new TH1D("el_etaSC","",42,-2.1,2.1); 
    TH1D* el_phi      = new TH1D("el_phi","",64,-3.2,3.2); 
    TH1D* emu_mass    = new TH1D("emu_mass","",binsize,bins);
    TH1D* emuSS_mass  = new TH1D("emuSS_mass","",binsize,bins);
    TH1D* dimu_mass   = new TH1D("dimu_mass","",binsize,bins);
    TH1D* dimuSS_mass = new TH1D("dimuSS_mass","",binsize,bins);

    el_etSC->Sumw2();
    el_etaSC->Sumw2();
    el_phi->Sumw2();
    emu_mass->Sumw2();
    emuSS_mass->Sumw2();
    dimu_mass->Sumw2();
    dimuSS_mass->Sumw2();

    cout<<"# of events = "<<chain->GetEntries()<<endl;

    bool leadingMu = false;
    bool leadingEl = false;
    double weight = 1.0;
    double weightedSum = 0;

    int tryEmu = 0;
    int totalMu = 0;
    int totalEl = 0;

    //Iteration
    //int entries = 100000;
    int entries = chain->GetEntries();
    
    for(int i=0; i!=entries; i++) { 
        chain->GetEntry(i);

        if( !isData ) {
            if( event->weight>0 ) weight = 1.0;
            else weight = -1.0;
        }
        weightedSum += weight;

        //if( ( index==0 || index==100 ) && !event->TauSelection() ) continue;

        if( !event->TriggerSelection("HLT_IsoMu20_v") && !event->TriggerSelection("HLT_IsoTkMu20_v") ) continue;

        leadingMu = false; 
        leadingEl = false; 

        passingMuons->clear();
        passingElectrons->clear();

        for(unsigned j=0; j!=event->muons.size(); j++) {

            PhysicsMuon* mu = (PhysicsMuon*)&event->muons.at(j);

            if( mu->acceptance(10,2.4) && mu->highPtMuonID() && mu->isolation(0.1) ) {
                totalMu += weight;
				//cout<<mu->pt<<endl;
                pair<PhysicsMuon,int> taggedmu = {*mu,j};
                passingMuons->push_back(taggedmu);
                if(mu->pt > 22) leadingMu = true;
            } 
        } 

        for(unsigned j=0; j!=event->electrons.size(); j++) {

            PhysicsElectron* el = (PhysicsElectron*)&event->electrons.at(j);

            if( el->acceptance(20,2.4) && el->WPMedium() ) { 
                totalEl += weight;
                pair<PhysicsElectron,int> taggedel = {*el,j};
                passingElectrons->push_back(taggedel);
                if(el->pt > 25) leadingEl = true;
            } 
        } 

        if( !leadingMu && !leadingEl ) continue;  // try to change acceptance

        if( (event->TriggerSelection("HLT_Ele22_eta2p1_WP75_Gsf_v") || event->TriggerSelection("HLT_Ele22_eta2p1_WPLoose_Gsf_v")) && (event->TriggerSelection("HLT_IsoMu20_v") || event->TriggerSelection("HLT_IsoTkMu20_v")) ) { //Single electron spectra
            for(unsigned j=0; j!=passingElectrons->size(); j++) { 
                PhysicsElectron el = passingElectrons->at(j).first;
                if( el.acceptance(30,2.1) ) { 
                    el_etSC->Fill(el.etSC,weight);
                    el_etaSC->Fill(el.etaSC,weight);
                    el_phi->Fill(el.phi,weight);
                } 
            } 
        } 

        if( event->TriggerSelection("HLT_IsoMu20_v") || event->TriggerSelection("HLT_IsoTkMu20_v") ) { 
            if( passingMuons->size() > 0 && passingElectrons > 0 ) ++tryEmu;
            if( emuDY(event->triggerobjects, event->emus, passingElectrons, passingMuons, emu) ) {
                double mass = ( emu->first.momentum() + emu->second.momentum() ).M();
                if( emu->first.charge * emu->second.charge < 0 ) emu_mass->Fill(mass,weight); 
                else emuSS_mass->Fill(mass,weight);
            } 
            if( dimuonDY(event->triggerobjects, event->dimuons, passingMuons, dimuon) ) { 
                double mass = ( dimuon->first.momentum() + dimuon->second.momentum() ).M();
                if( ( index==0 || index==100 ) && !event->TauSelection() ) continue;
                if( dimuon->first.charge * dimuon->second.charge <0 ) dimu_mass->Fill(mass,weight);
                else dimuSS_mass->Fill(mass,weight);
            } 
        } 
    }

    //Save 
    if(index==1) index=999;
    TFile* f = new TFile("histograms/hist"+TString::Itoa(index,10)+".root","RECREATE");

    el_etSC->Write();
    el_etaSC->Write();
    el_phi->Write();
    emu_mass->Write();
    emuSS_mass->Write();
    dimu_mass->Write();
    dimuSS_mass->Write();

    f->Close();

    cout<<"# of emu tries = "<<tryEmu<<endl;
    cout<<"# of passing muons = "<<totalMu<<endl;
    cout<<"# of passing electrons = "<<totalEl<<endl;
    cout<<"Weighted sum = "<<weightedSum<<endl;
    cout<<"Success"<<endl;
}
