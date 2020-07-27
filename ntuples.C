#define ntuples_cxx
#include "ntuples.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdio.h>
#include <algorithm>
#include <math.h>
#include <iostream>
float max(float& x, float& y){
  if(x > y) return x;
  else return y;
}
void ntuples::Loop()
{
  if (fChain == 0) return;
  TFile *f = new TFile("ntuples_RunB.root", "RECREATE");
  TTree *tree = new TTree("tree","Dimuon Tree");
  int     t_nPU;
  int     t_puBX;
  float   t_puTrue;

  int nMuons, nJets, nAK8Jets;
  int nVtxs, nGVtxs;
  float Vtx_x, Vtx_y, Vtx_z;  

  int JetPFLooseId[100]; 
  int JetID[100];
  float JetVtxMass[100];

  float MuPt[100];
  float MuEn[100];
  float MuEta[100];
  float MuPhi[100];
  float MuSIP[100];
  float MuInnerD0[100];
  float MuInnerDz[100];
  float MuIsoTrk[100];
  float MuPFChIso[100];
  float MuPFPhoIso[100];
  float MuPFNeuIso[100];
  float MuPFPUIso[100];
  float MuPFMiniIso[100];
  float MuInnervalidFraction[100];
  float MusegmentCompatibility[100];
  float Muchi2LocalPosition[100];
  float MutrkKink[100];
  float MuBestTrkPtError[100];
  float MuBestTrkPt[100];

  float   JetPt[100];
  float   JetEn[100];
  float   JetEta[100];
  float   JetPhi[100];
  float   aK8JetPt[100];
  float   aK8JetEn[100];
  float   aK8JetEta[100];
  float   aK8JetPhi[100];
  float   aK8JetMass[100];
  float   aK8Jet_tau1[100];
  float   aK8Jet_tau2[100];
  float   aK8Jet_tau3[100];
  float   aK8JetPrunedMass[100];
  float   aK8JetPrunedMassCorr[100];

  int MuPixelLayers[100];
  int MuMatches     [100];
  int MuCharge      [100];
  int MuType        [100];
  int MuTrkQuality  [100];
  unsigned int MuFiredTrgs   [100];
  unsigned int MuFiredL1Trgs [100];
  unsigned short MuIDbit[1000];
  
  //vector<bool> *b;
  tree->Branch("nPU", &t_nPU);
  tree->Branch("nPU_BX", &t_puBX);
  tree->Branch("nPU_true", &t_puTrue);

  tree->Branch("nVtxs",  &nVtxs, "nVtxs/I");
  tree->Branch("nGVtxs", &nGVtxs,"nGVtxs/I");
  tree->Branch("Vtx_x",  &Vtx_x, "Vtx_x/F"); 
  tree->Branch("Vtx_y",  &Vtx_y, "Vtx_y/F");
  tree->Branch("Vtx_z",  &Vtx_z, "Vtx_z/F");
 
  tree->Branch("nMuons", &nMuons, "nMuons/I");
  tree->Branch("nAK8Jets", &nAK8Jets, "nAK8Jets/I");
  tree->Branch("nJets", &nJets, "nJets/I");
  tree->Branch("JetPt",  &JetPt ,"JetPt[nJets]/F");
  tree->Branch("JetEn",  &JetEn ,"JetEn[nJets]/F");
  tree->Branch("JetEta",  &JetEta ,"JetEta[nJets]/F");
  tree->Branch("JetPhi",  &JetPhi ,"JetPhi[nJets]/F");
  tree->Branch("JetID",  &JetID, "JetID[nJets]/I");
  tree->Branch("JetVtxMass",  &JetVtxMass, "JetVtxMass[nJets]/F");
  tree->Branch("JetPFLooseId", &JetPFLooseId, "JetPFLooseId[nJets]/I");
  tree->Branch("aK8JetPt",  &aK8JetPt ,"aK8JetPt[nAK8Jets]/F");
  tree->Branch("aK8JetEn",  &aK8JetEn ,"aK8JetEn[nAK8Jets]/F");
  tree->Branch("aK8JetEta",  &aK8JetEta ,"aK8JetEta[nAK8Jets]/F");
  tree->Branch("aK8JetPhi",  &aK8JetPhi ,"aK8JetPhi[nAK8Jets]/F");
  tree->Branch("aK8JetMass",  &aK8JetMass ,"aK8JetMass[nAK8Jets]/F");
  tree->Branch("aK8Jet_tau1",  &aK8Jet_tau1 ,"aK8Jet_tau1[nAK8Jets]/F");
  tree->Branch("aK8Jet_tau2",  &aK8Jet_tau2 ,"aK8Jet_tau2[nAK8Jets]/F");
  tree->Branch("aK8Jet_tau3",  &aK8Jet_tau3 ,"aK8Jet_tau3[nAK8Jets]/F");
  tree->Branch("aK8JetPrunedMass",  &aK8JetPrunedMass ,"aK8JetPrunedMass[nAK8Jets]/F");
  tree->Branch("aK8JetPrunedMassCorr",  &aK8JetPrunedMassCorr ,"aK8JetPrunedMassCorr[nAK8Jets]/F");
  
  tree->Branch("MuPt", &MuPt, "MuPt[nMuons]/F");
  tree->Branch("MuEn", &MuEn, "MuEn[nMuons]/F");
  tree->Branch("MuEta", &MuEta, "MuEta[nMuons]/F");
  tree->Branch("MuPhi", &MuPhi, "MuPhi[nMuons]/F");
  tree->Branch("MuSIP", &MuSIP, "MuSIP[nMuons]/F");
  tree->Branch("MuInnerD0", &MuInnerD0, "MuInnerD0[nMuons]/F");
  tree->Branch("MuInnerDz", &MuInnerDz, "MuInnerDz[nMuons]/F");
  tree->Branch("MuIsoTrk", &MuIsoTrk, "MuIsoTrk[nMuons]/F");
  tree->Branch("MuPFChIso", &MuPFChIso, "MuPFChIso[nMuons]/F");
  tree->Branch("MuPFPhoIso", &MuPFPhoIso, "MuPFPhoIso[nMuons]/F");
  tree->Branch("MuPFNeuIso", &MuPFNeuIso, "MuPFNeuIso[nMuons]/F");
  tree->Branch("MuPFPUIso", &MuPFPUIso, "MuPFPUIso[nMuons]/F");
  tree->Branch("MuPFMiniIso", &MuPFMiniIso, "MuPFMiniIso[nMuons]/F");
  tree->Branch("MuInnervalidFraction", &MuInnervalidFraction, "MuInnervalidFraction[nMuons]/F");
  tree->Branch("MusegmentCompatibility", &MusegmentCompatibility, "MusegmentCompatibility[nMuons]/F");
  tree->Branch("Muchi2LocalPosition", &Muchi2LocalPosition, "Muchi2LocalPosition[nMuons]/F");
  tree->Branch("MutrkKink", &MutrkKink, "MutrkKink[nMuons]/F");
  tree->Branch("MuBestTrkPtError", &MuBestTrkPtError, "MuBestTrkPtError[nMuons]/F");
  tree->Branch("MuBestTrkPt", &MuBestTrkPt, "MuBestTrkPt[nMuons]/F");
  tree->Branch("MuPixelLayers", &MuPixelLayers ,"MuPixelLayers[nMuons]/I");
  tree->Branch("MuMatches", &MuMatches     ,"MuMatches[nMuons]/I");
  tree->Branch("MuCharge", &MuCharge      ,"MuCharge[nMuons]/I");
  tree->Branch("MuType", &MuType        ,"MuType[nMuons]/I");
  tree->Branch("MuTrkQuality", &MuTrkQuality  ,"MuTrkQuality[nMuons]/I");
  tree->Branch("MuFiredTrgs", &MuFiredTrgs   ,"MuFiredTrgs[nMuons]/i");
  tree->Branch("MuFiredL1Trgs", &MuFiredL1Trgs ,"MuFiredL1Trgs[nMuons]/i");
  tree->Branch("MuIDbit",       &MuIDbit       ,"MuIDbit[nMuons]/s");
  
  TH1F *nvtx = new TH1F("nvtx" , "number of primary vertices ", 50, 0, 50);
  TH1F *muTrkLayers_= new TH1F("muTrkLayers_","", 50, 0, 50);
  TH1F *muPixelHits_= new TH1F("muPixelHits_","", 50, 0, 50);
  TH1F *muMuonHits_= new TH1F("muMuonHits_","", 50, 0, 50);
  TH1F *muChi2NDF_= new TH1F("muChi2NDF_","", 100, -50, 50);
  TH1F *muStations_= new TH1F("muStations_","", 10, 0, 10);
  TH1F *muD0_= new TH1F("muD0_","", 100, -1, 1);
  TH1F *muDz_= new TH1F("muDz_","", 100, -1, 1);
  TH1F *muPt_= new TH1F("muPt_","", 500, 0, 500);
  TH1F *muEta_= new TH1F("muEta_","", 60, -3, 3);
  TH1F *hmu_iso = new TH1F("mu_iso", "", 100, -1,1);
  //TH1F *htrigdecision = new TH1F("trigdecision", "", 100, -1,1);
  Long64_t nentries = fChain->GetEntries();
  //nentries =100;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry%100000 ==0) cout << "processing " << jentry << endl;
    nvtx->Fill(nVtxs);
    /// if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) nGoodVtx_++;
    if(!nGoodVtx) continue; ///atleast one good vertext;
    nGVtxs = nGoodVtx;
    nVtxs  = nVtx;
    Vtx_x= vtx; 
    Vtx_y= vty; 
    Vtx_z= vtz; 
    nvtx->Fill(nGoodVtx);
    nMuons=nJets=nAK8Jets=0;
    ///trigger
    //b = trigDecision;
    //unsigned int triggered = b
    //cout << "size of trigDecision: " << trigDecision->size() <<endl; 
    //if(b->size()==3){
    //cout << (b->at(0))<<"  "<<(b->at(1))<<"  "<<(b->at(2))<<"  "<<endl;
    //cout << (trigDecision->at(0)) << " " <<endl;
    //}
    if(!trigDecision->at(0)==1) continue;
    //cout << "the value of trigger at at(0)= " << trigDecision->at(0) <<endl;
    //}
    //loop over the muonts
    for(unsigned int i=0; i !=nPU->size(); ++i){
    if( (*puBX)[i] ==0){
    t_nPU   = (*nPU)[i];
    t_puBX  = (*puBX)[i];
    t_puTrue= (*puTrue)[i];
      }
    }
    
    for(int im =0; im !=nMu; ++im){
      if( (*muPt)[im] < 35) continue;
      if( (*muEta)[im] >= 2.4) continue;
      if( (*muTrkLayers)[im] < 6) continue;//recoMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
      if( (*muPixelHits)[im] < 1) continue;//recoMu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
      if( (*muMuonHits)[im] < 1)  continue;//recoMu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
      if( (*muChi2NDF)[im] >= 10) continue;//χ2/ndof of the global-muon track fit < 10
      if( (*muStations)[im] < 2)  continue;//recoMu.numberOfMatchedStations() > 1
      if( fabs( (*muD0)[im] ) >= 0.2 )    continue;//fabs(recoMu.muonBestTrack()->dxy(vertex->position())) < 0.2
      if( fabs( (*muDz)[im] ) >= 0.5 )    continue;//fabs(recoMu.muonBestTrack()->dz(vertex->position())) < 0.5
      if( ((*muBestTrkPtError)[im]/(*muPt)[im]) >= 0.3) continue;
      //     muPFChIso_  .push_back(iMu->pfIsolationR04().sumChargedHadronPt);
      //    muPFPhoIso_ .push_back(iMu->pfIsolationR04().sumPhotonEt);
      //    muPFNeuIso_ .push_back(iMu->pfIsolationR04().sumNeutralHadronEt);
      //    muPFPUIso_  .push_back(iMu->pfIsolationR04().sumPUPt);
      // float muon_iso = float((MuonVector[i]->pfIsolationR04().sumChargedHadronPt + max(0., MuonVector[i]->pfIsolationR04().sumNeutralHadronEt + MuonVector[i]->pfIsolationR04().sumPhotonEt - 0.5*MuonVector[i]->pfIsolationR04().sumPUPt)))/MuonVector[i]->pt();
      //cout << (*muPFChIso)[im] << "\t" << (*muPFNeuIso)[im] << "\t" << (*muPFPhoIso)[im] << "\t" << (*muPFPUIso)[im] << "\t" << (*muPt)[im] << endl;
       float muon_iso = float(( (*muPFChIso)[im] + max(0., (*muPFNeuIso)[im] + (*muPFPhoIso)[im] - 0.5*(*muPFPUIso)[im])))/(*muPt)[im];
      //cout << "muon iso " << muon_iso << endl;
      if(muon_iso >= 0.15) continue;
      hmu_iso->Fill(muon_iso);
    
      if( (*muPt)[im] < 35 || (*muEta)[im] >= 2.4 || (*muPixelHits)[im] < 1 || (*muMuonHits)[im] < 1 || (*muChi2NDF)[im] >= 10
      || (*muStations)[im] < 2 || fabs( (*muD0)[im] ) >= 0.2 || fabs( (*muDz)[im] ) >= 0.5 ||  (*muTrkLayers)[im] < 6){
        cout << "error implement the cuts " << endl;
        exit(1);
      }
      MuPt                     [nMuons]  =  (*muPt)[im];
      MuEn                     [nMuons]  =  (*muEn)[im];
      MuEta                    [nMuons]  =  (*muEta)[im];
      MuPhi                    [nMuons]  =  (*muPhi)[im];
      MuSIP                    [nMuons]  =  (*muSIP)[im];
      MuInnerD0                [nMuons]  =  (*muInnerD0)[im];
      MuInnerDz                [nMuons]  =  (*muInnerDz)[im];
      MuIsoTrk                 [nMuons]  =  (*muIsoTrk)[im];
      //MuPFChIso                [nMuons]  =  (*muPFChIso)[im];
      //MuPFPhoIso               [nMuons]  =  (*muPFPhoIso)[im];
      //MuPFNeuIso               [nMuons]  =  (*muPFNeuIso)[im];
      //MuPFPUIso                [nMuons]  =  (*muPFPUIso)[im];
      //MuPFMiniIso              [nMuons]  =  (*muPFMiniIso)[im];
      //MuInnervalidFraction     [nMuons]  =  (*muInnervalidFraction)[im];
      //MusegmentCompatibility   [nMuons]  =  (*musegmentCompatibility)[im];
      MutrkKink                [nMuons]  =  (*mutrkKink)[im];
      MuBestTrkPtError         [nMuons]  =  (*muBestTrkPtError)[im];
      MuBestTrkPt              [nMuons]  =  (*muBestTrkPt)[im];
      Muchi2LocalPosition      [nMuons]  =  (*muchi2LocalPosition)[im];
      MuPixelLayers [nMuons]  = (*muPixelLayers)[im];
      MuMatches     [nMuons]  = (*muMatches    )[im];
      MuCharge      [nMuons]  = (*muCharge     )[im];
      MuType        [nMuons]  = (*muType       )[im];
      MuTrkQuality  [nMuons]  = (*muTrkQuality )[im];
      MuFiredTrgs   [nMuons]  = (*muFiredTrgs  )[im];
      MuFiredL1Trgs [nMuons]  = (*muFiredL1Trgs)[im];
      MuIDbit       [nMuons]  = (*muIDbit      )[im];
                
      ++nMuons;
    
      muTrkLayers_->Fill( (*muTrkLayers)[im]);//recoMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
      muPixelHits_->Fill( (*muPixelHits)[im]);//recoMu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0
      muMuonHits_ ->Fill( (*muMuonHits)[im]);//recoMu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0
      muChi2NDF_  ->Fill( (*muChi2NDF)[im]);//χ2/ndof of the global-muon track fit < 10
      muStations_ ->Fill( (*muStations)[im]);//recoMu.numberOfMatchedStations() > 1
      muD0_       ->Fill( fabs( (*muD0)[im] ) );//fabs(recoMu.muonBestTrack()->dxy(vertex->position())) < 0.2
      muDz_       ->Fill( fabs( (*muDz)[im] ) );//fabs(recoMu.muonBestTrack()->dz(vertex->position())) < 0.5
      muPt_       ->Fill( (*muPt)[im]);//fabs(recoMu.muonBestTrack()->dxy(vertex->position())) < 0.2
      muEta_       ->Fill( (*muEta)[im]);//fabs(recoMu.muonBestTrack()->dxy(vertex->position())) < 0.2
    
    }//loop over muon ends
    
    if(!nMuons) continue;

    for(int ij =0; ij !=nJet; ++ij){
      if((*jetPt)[ij] <= 20) continue;
      if(fabs ((*jetEta)[ij]) >= 2.5) continue;

      JetPt[nJets]                 = (*jetPt)[ij];
      JetEn[nJets]                 = (*jetEn)[ij];
      JetEta[nJets]                = (*jetEta)[ij];
      JetPhi[nJets]                = (*jetPhi)[ij];
      JetID[nJets]                 = (*jetID)[ij];
      if(jentry < 2){
        //cout << "ID is " << ij << "\t" << (*jetID)[ij] << "\t" << (*jetVtxMass)[ij]
        //<< "\t" << (*jetPFLooseId)[ij] << endl;

        //cout << nGoodVtx << " \t" <<Vtx_x << "\t " << Vtx_y  << "\t" << Vtx_z << endl;
      }
      JetVtxMass[nJets] = (*jetVtxMass)[ij];
      JetPFLooseId[nJets] = (*jetPFLooseId)[ij];
      ++nJets;
    }

    for(int ij =0; ij !=nAK8Jet; ++ij){
      aK8JetPt[nAK8Jets]              = (*AK8JetPt)[ij];
      aK8JetEn[nAK8Jets]              = (*AK8JetEn)[ij];
      aK8JetEta[nAK8Jets]             = (*AK8JetEta)[ij];
      aK8JetPhi[nAK8Jets]             = (*AK8JetPhi)[ij];
      aK8JetMass[nAK8Jets]            = (*AK8JetMass)[ij];
      aK8Jet_tau1[nAK8Jets]           = (*AK8Jet_tau1)[ij];
      aK8Jet_tau2[nAK8Jets]           = (*AK8Jet_tau2)[ij];
      aK8Jet_tau3[nAK8Jets]           = (*AK8Jet_tau3)[ij];
      aK8JetPrunedMass[nAK8Jets]      = (*AK8JetPrunedMass)[ij];
      aK8JetPrunedMassCorr[nAK8Jets]  = (*AK8JetPrunedMassCorr)[ij];
      ++nAK8Jets;
    }
    tree->Fill();
  }

  f->Write();
  f->Close();
}

