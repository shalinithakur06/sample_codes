#define DBMuon_cxx
#include "DBMuon.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void DBMuon::Loop()
{
   if (fChain == 0) return;
TFile *f = new TFile("data_sept082016.root","RECREATE");
TH1D *h_nMuons = new TH1D("h_nMuons","Number of selected Muons", 20, 0, 20);
TH1D *h_LeadMuonPt = new TH1D("h_LeadMuonPt","Lead Muons Pt", 500, 0, 500);
TH1D *h_2ndMuonPt = new TH1D("h_2ndMuonPt","Lead Muons Pt", 500, 0, 500);
TH1D *h_LeadMuonEta = new TH1D("h_LeadMuonEta","Lead Muons Eta", 100, -5, 5);
TH1D *h_2ndMuonEta = new TH1D("h_2ndMuonEta","Lead Muons Eta", 100, -5, 5);
   Long64_t nentries = fChain->GetEntries();
cout << "number of entries " << nentries << endl;
//nentries=1000;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
int nmuons(0);
for(int imu=0; imu !=nMu; ++imu){
if( (*muPt)[imu] < 24) continue;
if( abs((*muEta)[imu]) > 2.4) continue;
h_LeadMuonPt->Fill( (*muPt)[imu]);
h_LeadMuonEta->Fill( (*muEta)[imu]);
++nmuons;
}
if(nmuons>0) cout << "total muons and selected muons " << nMu << "\t" << nmuons << endl;
if(nmuons) h_nMuons->Fill(nmuons);
   }
f->Write();
f->Close();
}
