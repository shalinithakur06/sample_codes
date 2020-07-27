#include <memory>
#include <fstream>
#include <iostream>
#include <vector>
#include <TLorentzVector.h>
#include <stdlib.h>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TFile.h"
// user include files
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
//Triggers 
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

using namespace std;
using namespace edm;
class MuonTuple_data : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit MuonTuple_data(const edm::ParameterSet&);
      ~MuonTuple_data();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      void ClearTreeVectors();
      int run;
      int event;
      int lumi;
      int Vrtx;
      int gprim;
      int bXnum;
      double PU_npT;
      double PU_npIT;
      bool geninfo;
      double MC_EventWeight;
      double Ndof,X,Y,Z,Rho; 
      TTree* T;
      Service<TFileService>    fs;
      vector<double> *Mu_Pt,*Mu_Eta,*Mu_Phi,*Mu_E;
      vector<int>    *Mu_Charge;
      vector<double> *Mu_Chi2Ndf,*Mu_NumMuonHits,*Mu_NumMatchStat,*Mu_NumPixHits;
      vector<bool>   *Mu_PFMuon,*Mu_GlbMuon,*Mu_TrkMuon;
      vector<double> *Mu_pfisor03chargedhadron;
      vector<double> *Mu_pfisor03chargedparticle;
      vector<double> *Mu_pfisor03neutralhadron;   
      vector<double> *Mu_pfisor03photon;
      vector<double> *Mu_pfisor03neutralhadronht;
      vector<double> *Mu_pfisor03photonht;
      vector<double> *Mu_pfisor03pu;             
      vector<double> *Mu_pfisor04chargedhadron;
      vector<double> *Mu_pfisor04chargedparticle;
      vector<double> *Mu_pfisor04neutralhadron;
      vector<double> *Mu_pfisor04photon;
      vector<double> *Mu_pfisor04neutralhadronht;
      vector<double> *Mu_pfisor04photonht;       
      vector<double> *Mu_pfisor04pu;              
      vector<double> *Mu_trackerIsoSumPT;
      vector<double> *Mu_ecalIso;
      vector<double> *Mu_hcalIso;
      vector<double> *Mu_hoIso;
      vector<double> *Mu_ecalVetoIso;
      vector<double> *Mu_hcalVetoIso;
      vector<bool>   *Mu_isMedium;
      vector<bool>   *Mu_isTight;
      vector<int>    *TrigFired;
      vector<double> *Mu_TrkHits;
      vector<double> *Mu_Dxy;
      vector<double> *Mu_Dz;
      vector<std::string>      trigNames;
      vector<double> *trg1_Pt,*trg1_Eta,*trg1_Phi,*trg1_Energy;
      vector<double> *trg2_Pt,*trg2_Eta,*trg2_Phi,*trg2_Energy;
      EDGetTokenT<std::vector<pat::Muon>> muonsToken;
      EDGetTokenT<std::vector<pat::Electron>> eleToken;
      EDGetTokenT<reco::VertexCollection> vtxToken; 
      EDGetTokenT<edm::TriggerResults>    resToken;
      EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsTok_;
      EDGetTokenT<pat::PackedTriggerPrescales> trigPrescalesTok_;
      EDGetTokenT<double> rhoTag_;
      EDGetTokenT<reco::BeamSpot> bsmToken;
      EDGetTokenT<reco::ConversionCollection> hConversionsToken;
      vector<double> *Ele_Pt,*Ele_Eta,*Ele_Phi,*Ele_E;
      vector<int>    *Ele_Charge;
      vector <double> *Ele_fullsigmaIetaIeta;
      vector <double> *Ele_dPhiIn;
      vector <double> *Ele_HoverE;
      vector <double> *Ele_EoverP;
      vector <double> *Ele_misHit;
      vector <double> *Ele_dz;
      vector <double> *Ele_d0;
      vector <double> *Ele_dEtaInSeed; 
      vector <double> *Ele_PFIso;
      vector <bool>   *Ele_MatchConv; 
      vector <double> *Ele_superEta;
      EffectiveAreas areas_;
      float rho;
      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MuonTuple_data::MuonTuple_data(const edm::ParameterSet& iConfig):
  areas_( (iConfig.getParameter<FileInPath>("effAreasConfigFile")).fullPath())
{
   //now do what ever initialization is needed
   trigNames    = iConfig.getParameter<std::vector<std::string> >("Triggers");
   muonsToken   = consumes<std::vector<pat::Muon>>(iConfig.getParameter<InputTag> ("muonsInputTag"));
   eleToken     = consumes<std::vector<pat::Electron>>(iConfig.getParameter<InputTag> ("electronInputTag"));
   vtxToken     = consumes<reco::VertexCollection>(iConfig.getParameter<InputTag> ("vtxInputTag"));
   resToken     = consumes<edm::TriggerResults>(iConfig.getParameter<InputTag>("trgResInputTag"));
   triggerObjectsTok_ = (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<InputTag>("objecttag")));
   trigPrescalesTok_ = consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<InputTag>("prescaletag"));
   rhoTag_ = consumes<double>( iConfig.getParameter<InputTag>("rhoTag")); 
   bsmToken = consumes<reco::BeamSpot>(iConfig.getParameter<InputTag>("bsmInputTag"));
   hConversionsToken  = consumes<reco::ConversionCollection>(iConfig.getParameter<InputTag>("ConversionsInputTag"));
}


MuonTuple_data::~MuonTuple_data()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MuonTuple_data::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   ClearTreeVectors();
   run=0;  event=0;  lumi=0; bXnum=0;
   Vrtx=0; gprim=0;
   PU_npT=0.;PU_npIT=0.;
    //muon 
   Handle<std::vector<pat::Muon> > muon;
   iEvent.getByToken(muonsToken , muon); 
   //vertex
   Handle<reco::VertexCollection> vertexHandle;
   iEvent.getByToken(vtxToken , vertexHandle );

  const auto& pv = vertexHandle->at(0);
  for (std::vector<reco::Vertex>::const_iterator vtx = vertexHandle->begin(); vtx != vertexHandle->end(); ++vtx){
                Ndof=vtx->ndof();
                X=vtx->x();
                Y=vtx->y();
                Z=vtx->z();
                Rho=sqrt(X*X + Y*Y);
                if(vtx->isValid() && !vtx->isFake() && Ndof>4 && abs(Z)<=24 && Rho< 2){
                  gprim++;
               }
            }
   Vrtx=gprim;
   //trigger results
   Handle<edm::TriggerResults> triggerResults;
   iEvent.getByToken(resToken, triggerResults);

   Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
   iEvent.getByToken(triggerObjectsTok_, triggerObjects);
 
   Handle<pat::PackedTriggerPrescales> trigPrescales;
   iEvent.getByToken(trigPrescalesTok_,trigPrescales);

  
  if (!triggerResults.isValid()) {
    cout << "ProcessedTreeProducer::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
   }
    bool ok(false);
    if (triggerResults.isValid()) {
        const edm::TriggerNames & triggerNames = iEvent.triggerNames(*triggerResults);
        const std::vector<std::string> & triggerNames_ = triggerNames.triggerNames();
        for (unsigned int iHLT=0; iHLT<triggerResults->size(); iHLT++) {
          int hlt    = triggerResults->accept(iHLT);
           if(trigNames.size()==0) continue;
            for (unsigned int i=0; i<trigNames.size(); ++i) {
              if (triggerNames_[iHLT].find(trigNames[i].c_str())!=std::string::npos) {
                if (hlt > 0){
                  ok = true;
                  //cout << trigPrescales->getPrescaleForIndex(iHLT) << endl;
                  }
                TrigFired->push_back(hlt);
             }
            }
           }
           if(ok) {
           bool passesTriggerFilter = false;
           for (pat::TriggerObjectStandAlone obj : *triggerObjects) {
             obj.unpackPathNames(triggerNames); 
             for(const auto &filterLabel : obj.filterLabels()) {
               if(filterLabel == ("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")){
                 trg1_Pt->push_back(obj.pt());
                 trg1_Eta->push_back(obj.eta());
                 trg1_Phi->push_back(obj.phi()); 
                 trg1_Energy->push_back(obj.energy());
               }
               if(filterLabel == ("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09")){
                 trg2_Pt->push_back(obj.pt());
                 trg2_Eta->push_back(obj.eta()); 
                 trg2_Phi->push_back(obj.phi()); 
                 trg2_Energy->push_back(obj.energy());
               }
             }
           }
         }
      }
     for(std::vector<pat::Muon>::const_iterator it=muon->begin(); it!=muon->end(); ++it){
      if(!it->isGlobalMuon() && !it->isTrackerMuon()) continue;
      if(it->pt()<15.) continue;
      Mu_Pt->push_back(it->pt()); 
      Mu_Eta->push_back(it->eta());
      Mu_Phi->push_back(it->phi());
      Mu_E->push_back(it->energy());
      Mu_Charge->push_back(it->charge());
      Mu_GlbMuon->push_back(it->isGlobalMuon());
      Mu_TrkMuon->push_back(it->isTrackerMuon());
      Mu_PFMuon->push_back(it->isPFMuon());
      Mu_trackerIsoSumPT->push_back(it->isolationR03().sumPt);
      Mu_ecalIso->push_back(it->isolationR03().emEt);
      Mu_hcalIso->push_back(it->isolationR03().hadEt);
      Mu_hoIso->push_back(it->isolationR03().hoEt);
      Mu_ecalVetoIso->push_back(it->isolationR03().emVetoEt);
      Mu_hcalVetoIso->push_back(it->isolationR03().hadVetoEt);
      Mu_pfisor03chargedhadron->push_back(it->pfIsolationR03().sumChargedHadronPt);
      Mu_pfisor03chargedparticle->push_back(it->pfIsolationR03().sumChargedParticlePt);
      Mu_pfisor03neutralhadron->push_back(it->pfIsolationR03().sumNeutralHadronEt);
      Mu_pfisor03photon->push_back(it->pfIsolationR03().sumPhotonEt);
      Mu_pfisor03neutralhadronht->push_back(it->pfIsolationR03().sumNeutralHadronEtHighThreshold);
      Mu_pfisor03photonht->push_back(it->pfIsolationR03().sumPhotonEtHighThreshold);
      Mu_pfisor03pu->push_back(it->pfIsolationR03().sumPUPt);
      Mu_pfisor04chargedhadron->push_back(it->pfIsolationR04().sumChargedHadronPt);
      Mu_pfisor04chargedparticle->push_back(it->pfIsolationR04().sumChargedParticlePt);
      Mu_pfisor04neutralhadron->push_back(it->pfIsolationR04().sumNeutralHadronEt);
      Mu_pfisor04photon->push_back(it->pfIsolationR04().sumPhotonEt);
      Mu_pfisor04neutralhadronht->push_back(it->pfIsolationR04().sumNeutralHadronEtHighThreshold);
      Mu_pfisor04photonht->push_back(it->pfIsolationR04().sumPhotonEtHighThreshold);
      Mu_pfisor04pu->push_back(it->pfIsolationR04().sumPUPt);
      bool isMedium = muon::isLooseMuon(*it) && it->innerTrack()->validFraction() > 0.49 && muon::segmentCompatibility(*it) >0.451; 
      Mu_isMedium->push_back(isMedium);
      if(it->isGlobalMuon() )
	    {
	      Mu_NumMuonHits->push_back(it->globalTrack()->hitPattern().numberOfValidMuonHits()  );
	      Mu_NumPixHits->push_back(it->globalTrack()->hitPattern().numberOfValidPixelHits() );
	      Mu_Chi2Ndf->push_back(it->globalTrack()->normalizedChi2());
	    }
      Mu_NumMatchStat->push_back(it->numberOfMatchedStations());
      Mu_TrkHits->push_back(it->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      bool isTight = muon::isTightMuon(*it,pv);      
      Mu_isTight->push_back(isTight);

      Mu_Dxy->push_back(it->muonBestTrack()->dxy(pv.position()));
      Mu_Dz->push_back(it->muonBestTrack()->dz(pv.position()));      
   }
   run=iEvent.id().run();
   event=iEvent.id().event();
   lumi=iEvent.luminosityBlock();   
   bXnum= iEvent.bunchCrossing();
  
  //electrons
   Handle<vector<pat::Electron>> electron;
   iEvent.getByToken(eleToken, electron);
   
  Handle<double> rhoHandle;
  iEvent.getByToken(rhoTag_, rhoHandle);
  rho = *rhoHandle;

  Handle<reco::BeamSpot> bsm;
  iEvent.getByToken(bsmToken,bsm);
  const reco::BeamSpot BS = *bsm;
 
  Handle<reco::ConversionCollection> hConversions;
  iEvent.getByToken(hConversionsToken,hConversions);

  for(vector<pat::Electron>::const_iterator ite=electron->begin(); ite!=electron->end(); ++ite){
      Ele_Pt->push_back(ite->pt());
      Ele_Eta->push_back(ite->eta());
      Ele_Phi->push_back(ite->phi());
      Ele_E->push_back(ite->energy());
      Ele_Charge->push_back(ite->charge());
      Ele_superEta->push_back(ite->superCluster()->eta());
      if(vertexHandle->size() > 0){
      double d0 = (-1) * ite->gsfTrack()->dxy((*vertexHandle)[0].position());
      double dz = ite->gsfTrack()->dz( (*vertexHandle)[0].position());
      Ele_d0->push_back(d0);
      Ele_dz->push_back(dz);
    }
     double dEtaInSeed = 0.;
     //ID Variables 
       Ele_fullsigmaIetaIeta->push_back(ite->full5x5_sigmaIetaIeta());
       if(ite->superCluster().isNonnull() && ite->superCluster()->seed().isNonnull()) {
          dEtaInSeed = ite->deltaEtaSuperClusterTrackAtVtx() - ite->superCluster()->eta() + ite->superCluster()->seed()->eta();}
       else dEtaInSeed = -99;
       Ele_dEtaInSeed->push_back(dEtaInSeed);
       Ele_dPhiIn->push_back(ite->deltaPhiSuperClusterTrackAtVtx());
       Ele_HoverE->push_back(ite->hadronicOverEm()); 
       reco::GsfElectron::PflowIsolationVariables pfIso = ite->pfIsolationVariables();
       float abseta =  abs(ite->superCluster()->eta());
       float eA = areas_.getEffectiveArea(abseta);
       Ele_PFIso->push_back((pfIso.sumChargedHadronPt + max(0.0f,pfIso.sumNeutralHadronEt+pfIso.sumPhotonEt-eA*rho))/ite->pt()); 
       Ele_EoverP->push_back(abs(1.0 - ite->eSuperClusterOverP())/ite->ecalEnergy()); 
       float mHits = ite->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
       Ele_misHit->push_back(mHits);
        Ele_MatchConv->push_back(ConversionTools::hasMatchedConversion(reco::GsfElectron (*&*ite),hConversions,BS.position()));
  }
  T->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
MuonTuple_data::beginJob()
{
  T=fs->make<TTree>("T","MyTuple"); 
  TrigFired=new std::vector<int>();
  T->Branch("TrigFired","vector<int>",&TrigFired);
  trg1_Pt=new std::vector<double>();
  T->Branch("trg1_Pt","vector<double>",&trg1_Pt);
  trg1_Eta=new std::vector<double>();
  T->Branch("trg1_Eta","vector<double>",&trg1_Eta);
  trg1_Phi=new std::vector<double>();
  T->Branch("trg1_Phi","vector<double>",&trg1_Phi);
  trg1_Energy=new std::vector<double>();
  T->Branch("trg1_Energy","vector<double>",&trg1_Energy);
  trg2_Pt=new std::vector<double>();
  T->Branch("trg2_Pt","vector<double>",&trg2_Pt);
  trg2_Eta=new std::vector<double>();
  T->Branch("trg2_Eta","vector<double>",&trg2_Eta);
  trg2_Phi=new std::vector<double>();
  T->Branch("trg2_Phi","vector<double>",&trg2_Phi);
  trg2_Energy=new std::vector<double>();
  T->Branch("trg2_Energy","vector<double>",&trg2_Energy);

 
  Mu_Pt=new std::vector<double>();
  T->Branch("Mu_Pt","vector<double>",&Mu_Pt);
  Mu_Eta=new std::vector<double>();
  T->Branch("Mu_Eta","vector<double>",&Mu_Eta);
  Mu_Phi=new std::vector<double>();
  T->Branch("Mu_Phi","vector<double>",&Mu_Phi);
  Mu_E=new std::vector<double>();
  T->Branch("Mu_E","vector<double>",&Mu_E);
  Mu_Charge=new std::vector<int>();
  T->Branch("Mu_Charge","vector<int>",&Mu_Charge);
  //muon ID cut 
  Mu_isMedium=new std::vector<bool>();
  T->Branch("Mu_isMedium","vector<bool>",&Mu_isMedium);
  Mu_isTight=new std::vector<bool>();
  T->Branch("Mu_isTight","vector<bool>",&Mu_isTight);

  //muon ID
  Mu_Chi2Ndf=new std::vector<double>();
  T->Branch("Mu_Chi2Ndf","vector<double>",&Mu_Chi2Ndf);
  Mu_NumMuonHits=new std::vector<double>();
  T->Branch("Mu_NumMuonHits","vector<double>",&Mu_NumMuonHits);
  Mu_NumMatchStat=new std::vector<double>();
  T->Branch("Mu_NumMatchStat","vector<double>",&Mu_NumMatchStat);
  Mu_NumPixHits=new std::vector<double>();
  T->Branch("Mu_NumPixHits","vector<double>",&Mu_NumPixHits);
  Mu_PFMuon=new std::vector<bool>();
  T->Branch("Mu_PFMuon","vector<bool>",&Mu_PFMuon);
  Mu_GlbMuon=new std::vector<bool>();
  T->Branch("Mu_GlbMuon","vector<bool>",&Mu_GlbMuon);
  Mu_TrkMuon=new std::vector<bool>();
  T->Branch("Mu_TrkMuon","vector<bool>",&Mu_TrkMuon);
  Mu_TrkHits=new std::vector<double>();
  T->Branch("Mu_TrkHits","vector<double>",&Mu_TrkHits);
  Mu_Dxy=new std::vector<double>();
  T->Branch("Mu_Dxy","vector<double>",&Mu_Dxy);
  Mu_Dz=new std::vector<double>();
  T->Branch("Mu_Dz","vector<double>",&Mu_Dz);

 //muon ISO
  Mu_trackerIsoSumPT=new std::vector<double>();
  T->Branch("Mu_trackerIsoSumPT","vector<double>",&Mu_trackerIsoSumPT);
  Mu_ecalIso=new std::vector<double>();
  T->Branch("Mu_ecalIso","vector<double>",&Mu_ecalIso);
  Mu_hcalIso=new std::vector<double>();
  T->Branch("Mu_hcalIso","vector<double>",&Mu_hcalIso);
  Mu_hoIso=new std::vector<double>();
  T->Branch("Mu_hoIso","vector<double>",&Mu_hoIso);
  Mu_ecalVetoIso=new std::vector<double>();
  T->Branch("Mu_ecalVetoIso","vector<double>",&Mu_ecalVetoIso);
  Mu_hcalVetoIso=new std::vector<double>();
  T->Branch("Mu_hcalVetoIso","vector<double>",&Mu_hcalVetoIso);

  Mu_pfisor03chargedhadron=new std::vector<double>();
  T->Branch("Mu_pfisor03chargedhadron","vector<double>",&Mu_pfisor03chargedhadron);
  Mu_pfisor03chargedparticle=new std::vector<double>();
  T->Branch("Mu_pfisor03chargedparticle","vector<double>",&Mu_pfisor03chargedparticle);
  Mu_pfisor03neutralhadron=new std::vector<double>();
  T->Branch("Mu_pfisor03neutralhadron","vector<double>",&Mu_pfisor03neutralhadron);
  Mu_pfisor03photon=new std::vector<double>();
  T->Branch("Mu_pfisor03photon","vector<double>",&Mu_pfisor03photon); 
  Mu_pfisor03neutralhadronht=new std::vector<double>();
  T->Branch("Mu_pfisor03neutralhadronht","vector<double>",&Mu_pfisor03neutralhadronht);
  Mu_pfisor03photonht=new std::vector<double>();
  T->Branch("Mu_pfisor03photonht","vector<double>",&Mu_pfisor03photonht);
  Mu_pfisor03pu=new std::vector<double>();
  T->Branch("Mu_pfisor03pu","vector<double>",&Mu_pfisor03pu);
  Mu_pfisor04chargedhadron=new std::vector<double>();
  T->Branch("Mu_pfisor04chargedhadron","vector<double>",&Mu_pfisor04chargedhadron);
  Mu_pfisor04chargedparticle=new std::vector<double>();
  T->Branch("Mu_pfisor04chargedparticle","vector<double>",&Mu_pfisor04chargedparticle);
  Mu_pfisor04neutralhadron=new std::vector<double>();
  T->Branch("Mu_pfisor04neutralhadron","vector<double>",&Mu_pfisor04neutralhadron);
  Mu_pfisor04photon=new std::vector<double>();
  T->Branch("Mu_pfisor04photon","vector<double>",&Mu_pfisor04photon);
  Mu_pfisor04neutralhadronht =new std::vector<double>();
  T->Branch("Mu_pfisor04neutralhadronht","vector<double>",&Mu_pfisor04neutralhadronht);
  Mu_pfisor04photonht=new std::vector<double>();
  T->Branch("Mu_pfisor04photonht","vector<double>",&Mu_pfisor04photonht);
  Mu_pfisor04pu=new std::vector<double>();
  T->Branch("Mu_pfisor04pu","vector<double>",&Mu_pfisor04pu);

 //electrons
  Ele_Pt=new std::vector<double>();
  T->Branch("Ele_Pt","vector<double>",&Ele_Pt);
  Ele_Eta=new std::vector<double>();
  T->Branch("Ele_Eta","vector<double>",&Ele_Eta);
  Ele_Phi=new std::vector<double>();
  T->Branch("Ele_Phi","vector<double>",&Ele_Phi);
  Ele_E=new std::vector<double>();
  T->Branch("Ele_E","vector<double>",&Ele_E);
  Ele_Charge=new std::vector<int>();
  T->Branch("Ele_Charge","vector<int>",&Ele_Charge); 
  Ele_EoverP=new std::vector<double>();
  T->Branch("Ele_EoverP","vector<double>",&Ele_EoverP);
  Ele_misHit=new std::vector<double>();
  T->Branch("Ele_misHit","vector<double>",&Ele_misHit);
  Ele_superEta=new std::vector<double>();
  T->Branch("Ele_superEta","vector<double>",&Ele_superEta);
  Ele_d0=new std::vector<double>();
  T->Branch("Ele_d0","vector<double>",&Ele_d0);
  Ele_dz=new std::vector<double>();
  T->Branch("Ele_dz","vector<double>",&Ele_dz);
  //electron ID variables 
  Ele_fullsigmaIetaIeta=new std::vector<double>();
  T->Branch("Ele_fullsigmaIetaIeta","vector<double>",&Ele_fullsigmaIetaIeta);
  Ele_dEtaInSeed=new std::vector<double>();
  T->Branch("Ele_dEtaInSeed","vector<double>",&Ele_dEtaInSeed);
  Ele_dPhiIn=new std::vector<double>();
  T->Branch("Ele_dPhiIn","vector<double>",&Ele_dPhiIn);
  Ele_HoverE=new std::vector<double>();
  T->Branch("Ele_HoverE","vector<double>",&Ele_HoverE);
  Ele_PFIso=new std::vector<double>();
  T->Branch("Ele_PFIso","vector<double>",&Ele_PFIso);
  Ele_MatchConv=new std::vector<bool>();
  T->Branch("Ele_MatchConv","vector<bool>",&Ele_MatchConv);

//general information
  T->Branch("run",&run,"run/I");
  T->Branch("event",&event,"event/I");
  T->Branch("lumi",&lumi,"lumi/I");
  T->Branch("bXnum",&bXnum,"bXnum/I");
  T->Branch("Vrtx",&Vrtx,"Vrtx/I");
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MuonTuple_data::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonTuple_data::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
void MuonTuple_data::ClearTreeVectors()
    {
      TrigFired->clear(); 
      trg1_Pt->clear();trg1_Eta->clear();trg1_Phi->clear();trg1_Energy->clear();
      trg2_Pt->clear();trg2_Eta->clear();trg2_Phi->clear();trg2_Energy->clear();
      Mu_Pt->clear();Mu_Eta->clear();Mu_Phi->clear();Mu_E->clear();Mu_Charge->clear();Mu_TrkHits->clear();
      Mu_Chi2Ndf->clear();Mu_NumMuonHits->clear();Mu_NumMatchStat->clear();Mu_NumPixHits->clear();
      Mu_PFMuon->clear();Mu_GlbMuon->clear();Mu_TrkMuon->clear();Mu_isTight->clear();
      Mu_isMedium->clear();Mu_pfisor03chargedhadron->clear();Mu_pfisor03chargedparticle->clear();
      Mu_pfisor03neutralhadron->clear();Mu_pfisor03photon->clear();Mu_pfisor03neutralhadronht->clear();
      Mu_pfisor03photonht->clear();Mu_pfisor03pu->clear();Mu_Dxy->clear();Mu_Dz->clear();
      Mu_pfisor04chargedhadron->clear();Mu_pfisor04chargedparticle->clear();Mu_pfisor04neutralhadron->clear();
      Mu_pfisor04photon->clear();Mu_pfisor04neutralhadronht->clear();Mu_pfisor04photonht->clear();Mu_pfisor04pu->clear();            
      Mu_trackerIsoSumPT->clear();Mu_ecalIso->clear();Mu_hcalIso->clear();Mu_hoIso->clear();Mu_ecalVetoIso->clear();Mu_hcalVetoIso->clear();
      Ele_Pt->clear();Ele_Eta->clear();Ele_Phi->clear();Ele_E->clear();Ele_Charge->clear();Ele_superEta->clear();
      Ele_fullsigmaIetaIeta->clear();Ele_dPhiIn->clear();Ele_HoverE->clear();Ele_EoverP->clear();Ele_misHit->clear();
      Ele_d0->clear();Ele_dz->clear();Ele_dEtaInSeed->clear();Ele_PFIso->clear();Ele_MatchConv->clear();
    }

//define this as a plug-in
DEFINE_FWK_MODULE(MuonTuple_data);
