#define testTree_cxx
#include "testTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>


bool arrange(TLorentzVector t1, TLorentzVector t2){
return (t1.Pt() > t2.Pt());
}


void testTree::Loop(bool debug)
{
   if (fChain == 0) return;
TFile *f = new TFile("dataWC_2March.root","RECREATE");
TH1D *h_nMuons = new TH1D("h_nMuons","Number of selected Muons", 20, 0, 20);
TH1D *h_nJets = new TH1D("h_nJets","Number of selected Jets", 20, 0, 20);
TH1D *h_LeadMuonPt = new TH1D("h_LeadMuonPt","Muons Pt", 500, 0, 500);
TH1D *h_LeadMuonEta = new TH1D("h_LeadMuonEta","Muons Eta", 100, -5, 5);
TH1D *h_JetPt = new TH1D("h_JetPt","Jets Pt", 500, 0, 500);
TH1D *h_JetMass = new TH1D("h_JetMass","Jets Mass", 500, 0, 500);
TH1D *h_JetEta = new TH1D("h_JetEta","Jets Eta", 100, -5, 5);
TH1D *h_2ndMuonPt = new TH1D("h_2ndMuonPt","Muons Pt", 500, 0, 500);
TH1D *h_2ndMuonEta = new TH1D("h_2ndMuonEta","Muons Eta", 100, -5, 5);

TH1D *h_ZPt = new TH1D("h_ZPt","DiMuon Pt", 500, 0, 500);
TH1D *h_ZMass = new TH1D("h_ZMass","DiMuon Pt", 500, 0, 500);
TH1D *h_ZEta = new TH1D("h_ZEta","DiMuon Eta", 100, -5, 5);
//
TH1F *h_scalerHt = new TH1F("h_scalerHt", "Ht", 980, 40, 2000);
TH1F *h_nPV  = new TH1F("h_nPV", "Number of Primary Vertices", 50, 0, 50);
   Long64_t nentries = fChain->GetEntries();
if(debug) cout << "number of entries " << nentries << endl;


//put some counters to check the efficiency of each cut.
//number of muon hits > 0
//number of matched muon stations > 1
//transversal impact parameter dB < 0.2 cm
//longitudinal impact parameter dZ < 0.5 cm
//number of pixel hits > 0
//number of tracker layers with hits > 5
//relative error on track momentum spT/pT < 0.3
//track-based isolation relative < 0.1
int c_nPV(0), c_muHits(0), c_muStations(0), c_xyIP(0), c_zIP(0), c_muPixelHits(0), c_muTrkLayers(0), c_muTrkRelErr(0), c_muTrkRelIso(0), 
c_muChi2NDF(0); 
//nentries=100000;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
if(debug) cout << "Processing entry " << jentry << endl;
if(jentry % 1000000 ==0) cout << "Processing Event Number " << jentry << endl;
// normalisation
//double wnorm=0.195;
double wnorm=1.;

//Primary Vertex Selection
if(!hasGoodVtx) continue; 

//Muon Selection Cuts
//if(debug) cout << "number of muon intial " << nMu << endl;
vector<TLorentzVector> vm; //vector containers to store the muon related variables.
vm.clear();
for(int imu=0; imu !=nMu; ++imu){
/*if(debug)
cout 
<< (*muPt)[imu]
<< "\t" << abs((*muEta)[imu])
<< "\t" << (*muMuonHits)[imu]
<< "\t" << (*muStations)[imu]
<< "\t" << (*muD0)[imu]
<< "\t" << (*muDz)[imu]
<< "\t" << (*muPixelHits)[imu]
<< "\t" << (*muTrkLayers)[imu]
<< "\t" << (*muBestTrkPtError)[imu]/(*muPt)[imu]
<< "\t" << (*muIsoTrk)[imu]/(*muPt)[imu]
<< endl; */
if( (*muPt)[imu] < 35) continue;
if( fabs((*muEta)[imu]) > 2.4) continue;
++c_nPV; //counting number of events passing primary vertex cut;
//if( (*muChi2NDF)[imu] > 10) continue; ++c_muChi2NDF;
if( (*muMuonHits)[imu] < 1) continue; ++c_muHits;
if( (*muStations)[imu] < 2) continue; ++c_muStations;
if( (*muD0)[imu] >=0.2) continue;     ++c_xyIP;
if( (*muDz)[imu] >=0.5) continue;     ++c_zIP;
if( (*muPixelHits)[imu] < 1) continue;++c_muPixelHits;
if( (*muTrkLayers)[imu] < 6) continue;++c_muTrkLayers;
if( ((*muBestTrkPtError)[imu]/(*muPt)[imu]) >= 0.3) continue;++c_muTrkRelErr;
if( ((*muIsoTrk)[imu]/(*muPt)[imu]) >= 0.1) continue;++c_muTrkRelIso;
TLorentzVector m; 
m.SetPtEtaPhiE((*muPt)[imu], (*muEta)[imu], (*muPhi)[imu], (*muEn)[imu]);
vm.push_back(m);
}
//if(debug && nMu) cout << "Number of Selected Muon " << vm.size() << " selection efficiency " << (vm.size()/nMu) << endl;
if(vm.size() !=2) continue;
 h_nMuons->Fill(vm.size(), wnorm);
//if(debug) cout << "total muons and selected muons " << nMu << "\t" << vm.size() << endl;
sort(vm.begin(), vm.end(), arrange);

if(debug){
for(unsigned int imu=0; imu !=vm.size(); ++imu){
//cout << imu << "\t" << "Pt: " << vm[imu].Pt() << "\t" << vm[imu].Eta() << endl; 
}
}
if(vm.size() > 1){
h_nPV->Fill(nVtx);
h_LeadMuonPt-> Fill(vm[0].Pt(), wnorm);
h_LeadMuonEta->Fill(vm[0].Eta(), wnorm);
h_2ndMuonPt-> Fill(vm[1].Pt(), wnorm);
h_2ndMuonEta->Fill(vm[1].Eta(), wnorm);

//constructing the Z system
TLorentzVector muon_system = vm[0] + vm[1]; 
h_ZPt->Fill(muon_system.Pt(), wnorm);
h_ZEta->Fill(muon_system.Eta(), wnorm);
h_ZMass->Fill(muon_system.M(), wnorm);
}
//jet selection
vector<TLorentzVector> vj; //vector containers to store the muon related variables.
vector<double> jmass;
int nJetsHt(0);
float scalerHt(0.);
for(int iht=0; iht !=nJet; ++iht){
if( (*jetPt)[iht] > 50){
//if(debug) cout << "jet Pt/Eta " << (*jetPt)[iht] << "\t" << (*jetEta)[iht] << endl;
if(abs( (*jetEta)[iht] ) > 2.5 ) continue;
++nJetsHt;
scalerHt += (*jetPt)[iht];
}
}
//if(debug) cout << "Intial Jets " << nJet << "\t selected Jets\t" << nJetsHt << "\tscalerHt\t" << scalerHt << endl; 
if(scalerHt) h_scalerHt->Fill(scalerHt);
//
//
//
int njets(0);
//if(debug) cout << "\n";
for(int ij=0; ij !=nAK8Jet; ++ij){
if( (*AK8JetPt)[ij] <= 100) continue;
if( abs((*AK8JetEta)[ij]) > 2.4) continue;
//cout << "jet is passing the threshold" << nmuons << endl;
//if(debug) cout << "jet pt " << (*AK8JetPt)[ij] << "\t";
TLorentzVector j; 
j.SetPtEtaPhiE((*AK8JetPt)[ij], (*AK8JetEta)[ij], (*AK8JetPhi)[ij], (*AK8JetEn)[ij]);
//loop over the muons to check if muon is within 0.8 of jet 
bool deltaR = false;
for(unsigned int imu=0; imu !=vm.size(); ++imu){
//cout << j.DeltaR(vm[imu]) << endl;
if(j.DeltaR(vm[imu]) > 0.8) deltaR = true;
}
//cout << "deltaR" <<  deltaR << endl;
if(!deltaR) continue;
vj.push_back(j);
jmass.push_back( (*AK8CHSPrunedJetMass)[ij]);
if(jmass[ij] > 40) h_JetMass-> Fill(jmass[ij], wnorm);
h_JetPt-> Fill((*AK8JetPt)[ij], wnorm);
h_JetEta-> Fill((*AK8JetEta)[ij], wnorm);
++njets;
}
h_nJets->Fill(njets, wnorm);
if(njets){
for(unsigned int ij=0; ij !=vj.size(); ++ij){
//h_JetPt-> Fill(vj[ij].Pt());
//h_JetEta-> Fill(vj[ij].Eta());
//cout << "selected jet pt/eta" << vj[ij].Pt() << "\t" << vj[ij].Eta() <<endl;
 //exit(1);
 }
}
if(debug) cout << "*******************" << endl;

   }
f->Write();
f->Close();
}
