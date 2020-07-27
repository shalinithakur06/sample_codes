#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"


using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
// (local) variables associated with tree branches
//normal jets (ak4)
Int_t         nJet_;
vector<float> jetPt_;
vector<float> jetEn_;
vector<float> jetEta_;
vector<float> jetPhi_;
vector<float> jetRawPt_;
vector<float> jetRawEn_;
vector<float> jetMt_;
vector<float> jetArea_;
vector<float> jetLeadTrackPt_;
vector<float> jetLeadTrackEta_;
vector<float> jetLeadTrackPhi_;
vector<int>   jetLepTrackPID_;
vector<float> jetLepTrackPt_;
vector<float> jetLepTrackEta_;
vector<float> jetLepTrackPhi_;
vector<float> jetCHF_;
vector<float> jetNHF_;
vector<float> jetCEF_;
vector<float> jetNEF_;
vector<int>   jetNCH_;
vector<float> jetHFHAE_;
vector<float> jetHFEME_;
vector<int>   jetNConstituents_;
vector<float> AK8JetMUF_;
vector<float> jetVtxPt_;
vector<float> jetVtxMass_;
vector<float> jetVtxNtrks_;
vector<float> jetVtx3DVal_;
vector<float> jetVtx3DSig_;
vector<float> jetpfCombinedInclusiveSecondaryVertexV2BJetTags_; // recommended
vector<float> jetJetProbabilityBJetTags_;
vector<float> jetpfCombinedMVABJetTags_;
vector<int>   jetPartonID_;
vector<bool>  jetPFLooseId_;
vector<float> jetPUidFullDiscriminant_;
vector<float> jetJECUnc_;
vector<int>   jetFiredTrgs_;
//gen-info for ak4
vector<int>   jetGenJetIndex_;
vector<float> jetGenJetEn_;
vector<float> jetGenJetPt_;
vector<float> jetGenJetEta_;
vector<float> jetGenJetPhi_;
vector<int>   jetGenPartonID_;
vector<float> jetGenEn_;
vector<float> jetGenPt_;
vector<float> jetGenEta_;
vector<float> jetGenPhi_;
vector<int>   jetGenPartonMomID_;
//fat-jets (ak8)
Int_t         nAK8Jet_;
vector<float> AK8JetPt_;
vector<float> AK8JetEn_;
vector<float> AK8JetRawPt_;
vector<float> AK8JetRawEn_;
vector<float> AK8JetEta_;
vector<float> AK8JetPhi_;
vector<float> AK8JetMass_;
vector<float> AK8Jet_tau1_;
vector<float> AK8Jet_tau2_;
vector<float> AK8Jet_tau3_;
vector<float> AK8JetCHF_;
vector<float> AK8JetNHF_;
vector<float> AK8JetCEF_;
vector<float> AK8JetNEF_;
vector<int>   AK8JetNCH_;
vector<int>   AK8Jetnconstituents_;
vector<bool>  AK8JetPFLooseId_;
vector<bool>  AK8JetPFTightLepVetoId_;
vector<float> AK8CHSSoftDropJetMass_;
vector<float> AK8CHSSoftDropJetMassCorr_;
vector<float> AK8CHSPrunedJetMass_;
vector<float> AK8CHSPrunedJetMassCorr_;
vector<float> AK8JetpfBoostedDSVBTag_;
vector<float> AK8JetJECUnc_;
vector<float> AK8JetL2L3corr_;

//gen-info for ak8
vector<int>   AK8JetPartonID_;
vector<int>   AK8JetGenJetIndex_;
vector<float> AK8JetGenJetEn_;
vector<float> AK8JetGenJetPt_;
vector<float> AK8JetGenJetEta_;
vector<float> AK8JetGenJetPhi_;
vector<int>   AK8JetGenPartonID_;
vector<float> AK8JetGenEn_;
vector<float> AK8JetGenPt_;
vector<float> AK8JetGenEta_;
vector<float> AK8JetGenPhi_;
vector<int>   AK8JetGenPartonMomID_;
//soft drop subjets
vector<int>             nAK8softdropSubjet_ ;
vector< vector<float> > AK8softdropSubjetPt_ ;
vector< vector<float> > AK8softdropSubjetEta_ ;
vector< vector<float> > AK8softdropSubjetMass_ ;
vector< vector<float> > AK8softdropSubjetPhi_ ;
vector< vector<float> > AK8softdropSubjetE_ ;
vector< vector<int > >  AK8softdropSubjetCharge_ ;
vector< vector<int > >  AK8softdropSubjetFlavour_;
vector< vector<float> > AK8softdropSubjetCSV_ ;

void ggNtuplizer::branchesJets(TTree* tree) {
  
  tree->Branch("nJet",            &nJet_);
  tree->Branch("jetPt",           &jetPt_);
  tree->Branch("jetEn",           &jetEn_);
  tree->Branch("jetEta",          &jetEta_);
  tree->Branch("jetPhi",          &jetPhi_);
  tree->Branch("jetRawPt",        &jetRawPt_);
  tree->Branch("jetRawEn",        &jetRawEn_);
  tree->Branch("jetMt",           &jetMt_);
  tree->Branch("jetArea",         &jetArea_);
  tree->Branch("jetLeadTrackPt",  &jetLeadTrackPt_);
  tree->Branch("jetLeadTrackEta", &jetLeadTrackEta_);
  tree->Branch("jetLeadTrackPhi", &jetLeadTrackPhi_);
  tree->Branch("jetLepTrackPID",  &jetLepTrackPID_);
  tree->Branch("jetLepTrackPt",   &jetLepTrackPt_);
  tree->Branch("jetLepTrackEta",  &jetLepTrackEta_);
  tree->Branch("jetLepTrackPhi",  &jetLepTrackPhi_);
  tree->Branch("jetpfCombinedInclusiveSecondaryVertexV2BJetTags", &jetpfCombinedInclusiveSecondaryVertexV2BJetTags_);
  tree->Branch("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
  tree->Branch("jetpfCombinedMVABJetTags", &jetpfCombinedMVABJetTags_);
  if (doGenParticles_){
    tree->Branch("jetPartonID", &jetPartonID_);
    tree->Branch("jetGenJetIndex", &jetGenJetIndex_);
    tree->Branch("jetGenJetEn", &jetGenJetEn_);
    tree->Branch("jetGenJetPt", &jetGenJetPt_);
    tree->Branch("jetGenJetEta", &jetGenJetEta_);
    tree->Branch("jetGenJetPhi", &jetGenJetPhi_);
    tree->Branch("jetGenPartonID", &jetGenPartonID_);
    tree->Branch("jetGenEn", &jetGenEn_);
    tree->Branch("jetGenPt", &jetGenPt_);
    tree->Branch("jetGenEta", &jetGenEta_);
    tree->Branch("jetGenPhi", &jetGenPhi_);
    tree->Branch("jetGenPartonMomID", &jetGenPartonMomID_);
  }  
  tree->Branch("jetPFLooseId", &jetPFLooseId_);
  tree->Branch("jetPUidFullDiscriminant", &jetPUidFullDiscriminant_);
  tree->Branch("jetJECUnc",    &jetJECUnc_);
  tree->Branch("jetFiredTrgs", &jetFiredTrgs_);
  tree->Branch("jetCHF",       &jetCHF_);
  tree->Branch("jetNHF",       &jetNHF_);
  tree->Branch("jetCEF",       &jetCEF_);
  tree->Branch("jetNEF",       &jetNEF_);
  tree->Branch("jetNCH",       &jetNCH_);
  tree->Branch("jetVtxPt",     &jetVtxPt_);
  tree->Branch("jetVtxMass",   &jetVtxMass_);
  tree->Branch("jetVtxNtrks",  &jetVtxNtrks_);
  tree->Branch("jetVtx3DVal",  &jetVtx3DVal_);
  tree->Branch("jetVtx3DSig",  &jetVtx3DSig_);
  if (development_) {
    tree->Branch("jetHFHAE",         &jetHFHAE_);
    tree->Branch("jetHFEME",         &jetHFEME_);
    tree->Branch("jetNConstituents", &jetNConstituents_);
  }

  // SubJet
  if (dumpSubJets_) {
    tree->Branch("nAK8Jet",                  &nAK8Jet_);
    tree->Branch("AK8JetPt",                 &AK8JetPt_);
    tree->Branch("AK8JetEn",                 &AK8JetEn_);
    tree->Branch("AK8JetRawPt",              &AK8JetRawPt_);
    tree->Branch("AK8JetRawEn",              &AK8JetRawEn_);
    tree->Branch("AK8JetEta",                &AK8JetEta_);
    tree->Branch("AK8JetPhi",                &AK8JetPhi_);
    tree->Branch("AK8JetMass",               &AK8JetMass_);
    tree->Branch("AK8Jet_tau1",              &AK8Jet_tau1_);
    tree->Branch("AK8Jet_tau2",              &AK8Jet_tau2_);
    tree->Branch("AK8Jet_tau3",              &AK8Jet_tau3_);
    tree->Branch("AK8JetCHF",                &AK8JetCHF_);
    tree->Branch("AK8JetNHF",                &AK8JetNHF_);
    tree->Branch("AK8JetCEF",                &AK8JetCEF_);
    tree->Branch("AK8JetNEF",                &AK8JetNEF_);
    tree->Branch("AK8JetNCH",                &AK8JetNCH_);
    tree->Branch("AK8Jetnconstituents",      &AK8Jetnconstituents_);
    tree->Branch("AK8JetMUF",                &AK8JetMUF_);
    tree->Branch("AK8JetPFLooseId",          &AK8JetPFLooseId_);
    tree->Branch("AK8JetPFTightLepVetoId",   &AK8JetPFTightLepVetoId_);
    tree->Branch("AK8CHSSoftDropJetMass",    &AK8CHSSoftDropJetMass_);
    tree->Branch("AK8CHSSoftDropJetMassCorr",&AK8CHSSoftDropJetMassCorr_);
    tree->Branch("AK8CHSPrunedJetMass",      &AK8CHSPrunedJetMass_);
    tree->Branch("AK8CHSPrunedJetMassCorr",      &AK8CHSPrunedJetMassCorr_);
    tree->Branch("AK8JetpfBoostedDSVBTag",   &AK8JetpfBoostedDSVBTag_);
    tree->Branch("AK8JetJECUnc",             &AK8JetJECUnc_);
    tree->Branch("AK8JetL2L3corr",             &AK8JetL2L3corr_);

    if (doGenParticles_){
      tree->Branch("AK8JetPartonID", &AK8JetPartonID_);
      tree->Branch("AK8JetGenJetIndex", &AK8JetGenJetIndex_);
      tree->Branch("AK8JetGenJetEn", &AK8JetGenJetEn_);
      tree->Branch("AK8JetGenJetPt", &AK8JetGenJetPt_);
      tree->Branch("AK8JetGenJetEta", &AK8JetGenJetEta_);
      tree->Branch("AK8JetGenJetPhi", &AK8JetGenJetPhi_);
      tree->Branch("AK8JetGenPartonID", &AK8JetGenPartonID_);
      tree->Branch("AK8JetGenEn", &AK8JetGenEn_);
      tree->Branch("AK8JetGenPt", &AK8JetGenPt_);
      tree->Branch("AK8JetGenEta", &AK8JetGenEta_);
      tree->Branch("AK8JetGenPhi", &AK8JetGenPhi_);
      tree->Branch("AK8JetGenPartonMomID", &AK8JetGenPartonMomID_);
    }
    tree->Branch("nAK8softdropSubjet",       &nAK8softdropSubjet_);
    tree->Branch("AK8softdropSubjetPt",      &AK8softdropSubjetPt_);
    tree->Branch("AK8softdropSubjetEta",     &AK8softdropSubjetEta_);
    tree->Branch("AK8softdropSubjetPhi",     &AK8softdropSubjetPhi_);
    tree->Branch("AK8softdropSubjetMass",    &AK8softdropSubjetMass_);
    tree->Branch("AK8softdropSubjetE",       &AK8softdropSubjetE_);
    tree->Branch("AK8softdropSubjetCharge",  &AK8softdropSubjetCharge_);
    tree->Branch("AK8softdropSubjetFlavour", &AK8softdropSubjetFlavour_);
    tree->Branch("AK8softdropSubjetCSV",     &AK8softdropSubjetCSV_);
  }
}

void ggNtuplizer::fillJets(const edm::Event& e, const edm::EventSetup& es) {

  // cleanup from previous execution
  jetPt_                                  .clear();
  jetEn_                                  .clear();
  jetEta_                                 .clear();
  jetPhi_                                 .clear();
  jetRawPt_                               .clear();
  jetRawEn_                               .clear();
  jetMt_                                  .clear();
  jetArea_                                .clear();
  jetLeadTrackPt_                         .clear();
  jetLeadTrackEta_                        .clear();
  jetLeadTrackPhi_                        .clear();
  jetLepTrackPt_                          .clear();
  jetLepTrackPID_                         .clear();
  jetLepTrackEta_                         .clear();
  jetLepTrackPhi_                         .clear();
  jetpfCombinedInclusiveSecondaryVertexV2BJetTags_.clear();
  jetJetProbabilityBJetTags_              .clear();
  jetpfCombinedMVABJetTags_               .clear();
  jetPartonID_                            .clear();
  jetPFLooseId_                           .clear();
  jetPUidFullDiscriminant_                .clear();
  jetJECUnc_                              .clear();
  jetFiredTrgs_                           .clear();
  jetCHF_                                 .clear();
  jetNHF_                                 .clear();
  jetCEF_                                 .clear();
  jetNEF_                                 .clear();
  jetNCH_                                 .clear();
  jetVtxPt_                               .clear();
  jetVtxMass_                             .clear();
  jetVtxNtrks_                            .clear();
  jetVtx3DVal_                            .clear();
  jetVtx3DSig_                            .clear();
  if (development_) {
    jetHFHAE_                               .clear();
    jetHFEME_                               .clear();
    jetNConstituents_                       .clear();
  }
  jetGenJetIndex_.clear();
  jetGenJetEn_.clear();
  jetGenJetPt_.clear();
  jetGenJetEta_.clear();
  jetGenJetPhi_.clear();
  jetGenPartonID_.clear();
  jetGenEn_.clear();
  jetGenPt_.clear();
  jetGenEta_.clear();
  jetGenPhi_.clear();
  jetGenPartonMomID_.clear();
  
  // SubJet
  AK8JetPt_           .clear();
  AK8JetEn_           .clear();
  AK8JetRawPt_        .clear();
  AK8JetRawEn_        .clear();
  AK8JetEta_          .clear();
  AK8JetPhi_          .clear();
  AK8JetMass_         .clear();
  AK8Jet_tau1_        .clear();
  AK8Jet_tau2_        .clear();
  AK8Jet_tau3_        .clear();
  AK8JetCHF_          .clear();
  AK8JetNHF_          .clear();
  AK8JetCEF_          .clear();
  AK8JetNEF_          .clear();
  AK8JetNCH_          	.clear();
  AK8Jetnconstituents_     .clear();
  AK8JetMUF_          .clear();
  AK8JetPFLooseId_         .clear();
  AK8JetPFTightLepVetoId_  .clear();
  AK8CHSSoftDropJetMass_   .clear();
  AK8CHSSoftDropJetMassCorr_   .clear();
  AK8CHSPrunedJetMass_   .clear();
  AK8CHSPrunedJetMassCorr_   .clear();
  AK8JetpfBoostedDSVBTag_  .clear();
  AK8JetJECUnc_            .clear();
  AK8JetL2L3corr_            .clear();

  AK8JetPartonID_ .clear();
  AK8JetGenJetIndex_.clear();
  AK8JetGenJetEn_.clear();
  AK8JetGenJetPt_.clear();
  AK8JetGenJetEta_.clear();
  AK8JetGenJetPhi_.clear();
  AK8JetGenPartonID_.clear();
  AK8JetGenEn_.clear();
  AK8JetGenPt_.clear();
  AK8JetGenEta_.clear();
  AK8JetGenPhi_.clear();
  AK8JetGenPartonMomID_.clear();

  nAK8softdropSubjet_ .clear();
  AK8softdropSubjetPt_ .clear();
  AK8softdropSubjetEta_ .clear();
  AK8softdropSubjetPhi_ .clear();
  AK8softdropSubjetMass_ .clear();
  AK8softdropSubjetCharge_ .clear();
  AK8softdropSubjetE_ .clear();
  AK8softdropSubjetFlavour_ .clear();
  AK8softdropSubjetCSV_ .clear();

  nJet_ = 0;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetsAK4Label_, jetHandle);

  if (!jetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
    return;
  }

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  if(doGenParticles_)e.getByToken(genParticlesCollection_, genParticlesHandle);
  
  // Accessing the JEC uncertainties 
  //ak4  
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc=0;
  jecUnc = new JetCorrectionUncertainty(JetCorPar);
  //ak8
  edm::ESHandle<JetCorrectorParametersCollection> AK8JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK8PFchs",AK8JetCorParColl);
  JetCorrectorParameters const & AK8JetCorPar = (*AK8JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *AK8jecUnc=0;
  AK8jecUnc = new JetCorrectionUncertainty(AK8JetCorPar);
  
  
  //start jets Lvdp
  for (edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {
    
    if (isAOD_ && iJet->pt() < 10) continue;
    jetPt_.push_back(    iJet->pt());
    jetEn_.push_back(    iJet->energy());
    jetEta_.push_back(   iJet->eta());
    jetPhi_.push_back(   iJet->phi());
    jetRawPt_.push_back( (*iJet).correctedJet("Uncorrected").pt());
    jetRawEn_.push_back( (*iJet).correctedJet("Uncorrected").energy());
    jetMt_.push_back(    iJet->mt());
    jetArea_.push_back(  iJet->jetArea());
    jetCEF_.push_back(   iJet->chargedEmEnergyFraction());
    jetNEF_.push_back(   iJet->neutralEmEnergyFraction());
    jetCHF_.push_back(   iJet->chargedHadronEnergyFraction());
    jetNHF_.push_back(   iJet->neutralHadronEnergyFraction());
    jetNCH_.push_back(   iJet->chargedMultiplicity());
    if (development_) {
      jetHFHAE_.push_back( iJet->HFHadronEnergy());
      jetHFEME_.push_back( iJet->HFEMEnergy());
      jetNConstituents_.push_back(iJet->numberOfDaughters());
    }

    if (fabs(iJet->eta()) < 5.2) {
      jecUnc->setJetEta(iJet->eta());
      jecUnc->setJetPt(iJet->pt()); // here you must use the CORRECTED jet pt
      jetJECUnc_.push_back(jecUnc->getUncertainty(true));
    } else {
      jetJECUnc_.push_back(-1.);
    }
    
    jetFiredTrgs_.push_back(matchJetTriggerFilters(iJet->pt(), iJet->eta(), iJet->phi()));    

    //Searching for leading track and lepton
    float leadTrkPt  = -99;
    float leadTrkEta = -99;
    float leadTrkPhi = -99;
    int   lepTrkPID  = -99;
    float lepTrkPt   = -99;
    float lepTrkEta  = -99;
    float lepTrkPhi  = -99;

    for (unsigned id = 0; id < iJet->getJetConstituents().size(); id++) {

      const edm::Ptr<reco::Candidate> daughter = iJet->getJetConstituents().at(id);

      if (daughter.isNonnull() && daughter.isAvailable()) {
	if (daughter->charge() != 0 && daughter->pt() > leadTrkPt) {
	  leadTrkPt  = daughter->pt();
	  leadTrkEta = daughter->eta();
	  leadTrkPhi = daughter->phi();
	}

	if (abs(daughter->pdgId()) == 11 || abs(daughter->pdgId()) == 13) {
	  if (daughter->pt() > lepTrkPt) {
	    lepTrkPID = daughter->pdgId();
	    lepTrkPt  = daughter->pt();
	    lepTrkEta = daughter->eta();
	    lepTrkPhi = daughter->phi();
	  }
	}
      }
    }

    jetLeadTrackPt_ .push_back(leadTrkPt);
    jetLeadTrackEta_.push_back(leadTrkEta);
    jetLeadTrackPhi_.push_back(leadTrkPhi);
    jetLepTrackPID_ .push_back(lepTrkPID);
    jetLepTrackPt_  .push_back(lepTrkPt);
    jetLepTrackEta_ .push_back(lepTrkEta);
    jetLepTrackPhi_ .push_back(lepTrkPhi);    
    jetVtxPt_       .push_back(sqrt(pow(iJet->userFloat("vtxPx"),2)+pow(iJet->userFloat("vtxPy"),2)));
    jetVtxMass_     .push_back(iJet->userFloat("vtxMass"));
    jetVtxNtrks_    .push_back(iJet->userFloat("vtxNtracks"));
    jetVtx3DVal_    .push_back(iJet->userFloat("vtx3DVal"));
    jetVtx3DSig_    .push_back(iJet->userFloat("vtx3DSig"));
    
    //b-tagging
    jetpfCombinedInclusiveSecondaryVertexV2BJetTags_.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    jetJetProbabilityBJetTags_.push_back(iJet->bDiscriminator("pfJetProbabilityBJetTags"));
    jetpfCombinedMVABJetTags_.push_back(iJet->bDiscriminator("pfCombinedMVABJetTags"));
  
    //parton id
    jetPartonID_.push_back(iJet->partonFlavour());

    //jet PF Loose ID
    //pat::strbitset retjet = pfLooseId_.getBitTemplate();
    //jetPFLooseId_.push_back(pfLooseId_(*iJet, retjet));
    bool jetID = true;
    if (fabs(iJet->eta()) <= 3.0) {
      if (!(iJet->neutralHadronEnergyFraction() < 0.99))                       jetID = false;
      if (!(iJet->neutralEmEnergyFraction() < 0.99))                           jetID = false;
      if (!((iJet->chargedMultiplicity() + iJet->neutralMultiplicity()) > 1))  jetID = false;
      if (fabs(iJet->eta()) <= 2.4) {
        if (!(iJet->chargedHadronEnergyFraction() > 0))  jetID = false;
        if (!(iJet->chargedMultiplicity() > 0))          jetID = false;
        if (!(iJet->chargedEmEnergyFraction() < 0.99))   jetID = false;
      }
    }
    if (fabs(iJet->eta()) > 3.0) {
      if (!(iJet->neutralEmEnergyFraction() < 0.90))  jetID = false;
      if (!(iJet->neutralMultiplicity() > 10))        jetID = false;
    }
    jetPFLooseId_.push_back(jetID);
    
    //PUJet ID from slimmedJets
    jetPUidFullDiscriminant_.push_back( iJet->userFloat("pileupJetId:fullDiscriminant"));
    
    // gen jet and parton
    int jetGenPartonID = -99;
    int jetGenPartonMomID = -99;
    float jetGenEn = -999.;
    float jetGenPt = -999.;
    float jetGenEta = -999.;
    float jetGenPhi = -999.;
    if (doGenParticles_ && genParticlesHandle.isValid() ) {
      if ((*iJet).genParton()) {
	jetGenPartonID = (*iJet).genParton()->pdgId();
	jetGenEn = (*iJet).genParton()->energy();
	jetGenPt = (*iJet).genParton()->pt();
	jetGenEta = (*iJet).genParton()->eta();
	jetGenPhi = (*iJet).genParton()->phi();
	if ((*iJet).genParton()->mother()) {
	  jetGenPartonMomID = (*iJet).genParton()->mother()->pdgId();
	}
      }
    }
    jetGenPartonID_.push_back(jetGenPartonID);
    jetGenPartonMomID_.push_back(jetGenPartonMomID);
    jetGenEn_ .push_back(jetGenEn);
    jetGenPt_ .push_back(jetGenPt);
    jetGenEta_ .push_back(jetGenEta);
    jetGenPhi_ .push_back(jetGenPhi);
    int jetGenJetIndex = -1;
    float jetGenJetEn = -999.;
    float jetGenJetPt = -999.;
    float jetGenJetEta = -999.;
    float jetGenJetPhi = -999.;
    if (doGenParticles_ && genParticlesHandle.isValid() ) {
      if ((*iJet).genJet()) {
	jetGenJetIndex = 1;
	jetGenJetEn = (*iJet).genJet()->energy();
	jetGenJetPt = (*iJet).genJet()->pt();
	jetGenJetEta = (*iJet).genJet()->eta();
	jetGenJetPhi = (*iJet).genJet()->phi();
      }
    }
    jetGenJetIndex_.push_back(jetGenJetIndex);
    jetGenJetEn_.push_back(jetGenJetEn);
    jetGenJetPt_.push_back(jetGenJetPt);
    jetGenJetEta_.push_back(jetGenJetEta);
    jetGenJetPhi_.push_back(jetGenJetPhi);
    
    nJet_++;
  }
  
  if (dumpSubJets_) {
    edm::Handle<edm::View<pat::Jet> > jetsAK8;
    e.getByToken(jetsAK8Label_, jetsAK8);
    
    if (!jetsAK8.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK8AK8) in event";
      return;
    }
    
    //Double B-tagger 
    edm::Handle<reco::JetTagCollection> pfBoostedDoubleSecondaryVertex; 
    e.getByLabel("pfBoostedDoubleSecondaryVertexAK8BJetTags",pfBoostedDoubleSecondaryVertex); 
    
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetWtagging#Recipes_to_apply_JEC_on_the_prun
    
    std::vector<JetCorrectorParameters> vPar;
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
      JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
    }
    
    // Make the FactorizedJetCorrector
    jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
    
    edm::Handle<double> rhoHandle;
    e.getByToken(rhoLabel_, rhoHandle);
    float rho = *(rhoHandle.product());
    
    edm::Handle<reco::VertexCollection> vtxHandle;
    e.getByToken(vtxLabel_, vtxHandle);
    if (!vtxHandle.isValid()) edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";
    
    
    nAK8Jet_ = 0;
    //jet substructure
    int nsubjets = 0;
    std::vector<float> vecSoftdropSubjetcsv ;
    std::vector<float> vecSoftdropSubjetpt ;
    std::vector<float> vecSoftdropSubjeteta ;
    std::vector<float> vecSoftdropSubjetmass ;
    std::vector<float> vecSoftdropSubjetphi ;
    std::vector<float> vecSoftdropSubjete ;
    std::vector<int > vecSoftdropSubjetcharge ;
    std::vector<int > vecSoftdropSubjetflavour;
    
    edm::View<pat::Jet>::const_iterator beginAK8 = jetsAK8->begin();
    edm::View<pat::Jet>::const_iterator endAK8 = jetsAK8->end();
    edm::View<pat::Jet>::const_iterator ijetAK8 = beginAK8;
    int ijetRef = -1;   
    // Loop over the "hard" jets
    for(ijetAK8 = beginAK8; ijetAK8 != endAK8; ++ijetAK8 ) {
      ijetRef++;
      if( ijetAK8->pt() < 30.0 ) continue;
      nAK8Jet_++;
      AK8JetPt_.push_back( ijetAK8->pt() );
      AK8JetEn_.push_back( ijetAK8->energy() );
      AK8JetRawPt_.push_back( (*ijetAK8).correctedJet("Uncorrected").pt() );
      AK8JetRawEn_.push_back( (*ijetAK8).correctedJet("Uncorrected").energy() );
      AK8JetEta_.push_back( ijetAK8->eta() );
      AK8JetPhi_.push_back( ijetAK8->phi() );
      AK8JetMass_.push_back( ijetAK8->mass() );
      AK8Jet_tau1_.push_back( ijetAK8->userFloat("NjettinessAK8:tau1") );
      AK8Jet_tau2_.push_back( ijetAK8->userFloat("NjettinessAK8:tau2") );
      AK8Jet_tau3_.push_back( ijetAK8->userFloat("NjettinessAK8:tau3") );
      AK8JetCHF_.push_back( ijetAK8->chargedHadronEnergyFraction()); // 0.0
      AK8JetNHF_.push_back( ijetAK8->neutralHadronEnergyFraction()); //0.99
      AK8JetCEF_.push_back( ijetAK8->chargedEmEnergyFraction()); //0.99
      AK8JetNEF_.push_back( ijetAK8->neutralEmEnergyFraction()); //0.99
      AK8JetNCH_.push_back( ijetAK8->chargedMultiplicity()); //0
      AK8Jetnconstituents_.push_back( ijetAK8->chargedMultiplicity() + ijetAK8->neutralMultiplicity()); //1  
      AK8JetMUF_.push_back(ijetAK8->muonEnergyFraction());
      //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
      //LooseID
      bool AK8jetID = true;
      if (fabs(ijetAK8->eta()) <= 3.0) {
	if (!(ijetAK8->neutralHadronEnergyFraction() < 0.99))                       AK8jetID = false;
	if (!(ijetAK8->neutralEmEnergyFraction() < 0.99))                           AK8jetID = false;
	if (!((ijetAK8->chargedMultiplicity() + ijetAK8->neutralMultiplicity()) > 1))  AK8jetID = false;
	if (fabs(ijetAK8->eta()) <= 2.4) {
	  if (!(ijetAK8->chargedHadronEnergyFraction() > 0))  AK8jetID = false;
	  if (!(ijetAK8->chargedMultiplicity() > 0))          AK8jetID = false;
	  if (!(ijetAK8->chargedEmEnergyFraction() < 0.99))   AK8jetID = false;
	}
      }
      if (fabs(ijetAK8->eta()) > 3.0) {
	if (!(ijetAK8->neutralEmEnergyFraction() < 0.90))  AK8jetID = false;
	if (!(ijetAK8->neutralMultiplicity() > 10))        AK8jetID = false;
      }
      AK8JetPFLooseId_.push_back(AK8jetID);

      //TightIDMuon Fraction
      bool AK8jetIDTightLepVeto = true;
      if (fabs(ijetAK8->eta()) <= 3.0) {
        if (!(ijetAK8->neutralHadronEnergyFraction() < 0.90))                       AK8jetIDTightLepVeto = false;
        if (!(ijetAK8->neutralEmEnergyFraction() < 0.90))                           AK8jetIDTightLepVeto = false;
        if (!((ijetAK8->chargedMultiplicity() + ijetAK8->neutralMultiplicity()) > 1))  AK8jetIDTightLepVeto = false;
        if (!(ijetAK8->muonEnergyFraction()<0.8))                                   AK8jetIDTightLepVeto = false;
        if (fabs(ijetAK8->eta()) <= 2.4) {
          if (!(ijetAK8->chargedHadronEnergyFraction() > 0))  AK8jetIDTightLepVeto = false;
          if (!(ijetAK8->chargedMultiplicity() > 0))          AK8jetIDTightLepVeto = false;
          if (!(ijetAK8->chargedEmEnergyFraction() < 0.90))   AK8jetIDTightLepVeto = false;
        }
      }
      AK8JetPFTightLepVetoId_.push_back(AK8jetIDTightLepVeto);

      AK8CHSSoftDropJetMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSSoftDropMass")); 
      AK8CHSPrunedJetMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSPrunedMass")); 
      // AK8JetpfBoostedDSVBTag_.push_back(ijetAK8->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
      AK8JetpfBoostedDSVBTag_.push_back(     (*pfBoostedDoubleSecondaryVertex).value(ijetRef));
      
      const LorentzVector uncorrJet = (*ijetAK8).correctedP4(0);
      jecAK8_->setJetEta(uncorrJet.eta());
      jecAK8_->setJetPt ( uncorrJet.pt() );
      jecAK8_->setJetE  ( uncorrJet.energy() );
      jecAK8_->setJetA  ( (*ijetAK8).jetArea() );
      jecAK8_->setRho   ( rho );
      jecAK8_->setNPV   ( vtxHandle->size() );
      
      float corr = jecAK8_->getCorrection();
      AK8JetL2L3corr_.push_back(corr);
      AK8CHSSoftDropJetMassCorr_.push_back(corr*(ijetAK8->userFloat("ak8PFJetsCHSSoftDropMass")));
      AK8CHSPrunedJetMassCorr_.push_back(corr*(ijetAK8->userFloat("ak8PFJetsCHSPrunedMass")));
      
      //JEC uncertainty
      if (fabs(ijetAK8->eta()) < 5.2) {
	AK8jecUnc->setJetEta(ijetAK8->eta());
	AK8jecUnc->setJetPt(ijetAK8->pt()); // here you must use the CORRECTED jet pt
	AK8JetJECUnc_.push_back(AK8jecUnc->getUncertainty(true));
      } else {
	AK8JetJECUnc_.push_back(-1.);
      }
      
      
      //save gen-info for ak8 jets
      //parton id                                                                                                                                                           
      AK8JetPartonID_.push_back(ijetAK8->partonFlavour());
      int AK8JetGenPartonID = -99;
      int AK8JetGenPartonMomID = -99;
      float AK8JetGenEn = -999.;
      float AK8JetGenPt = -999.;
      float AK8JetGenEta = -999.;
      float AK8JetGenPhi = -999.;
      if (doGenParticles_ && genParticlesHandle.isValid() ) {
	if ((*ijetAK8).genParton()) {
	  AK8JetGenPartonID = (*ijetAK8).genParton()->pdgId();
	  AK8JetGenEn = (*ijetAK8).genParton()->energy();
	  AK8JetGenPt = (*ijetAK8).genParton()->pt();
	  AK8JetGenEta = (*ijetAK8).genParton()->eta();
	  AK8JetGenPhi = (*ijetAK8).genParton()->phi();
	  if ((*ijetAK8).genParton()->mother()) {
	    AK8JetGenPartonMomID = (*ijetAK8).genParton()->mother()->pdgId();
	  }
	}
      }
      AK8JetGenPartonID_.push_back(AK8JetGenPartonID);
      AK8JetGenPartonMomID_.push_back(AK8JetGenPartonMomID);
      AK8JetGenEn_ .push_back(AK8JetGenEn);
      AK8JetGenPt_ .push_back(AK8JetGenPt);
      AK8JetGenEta_ .push_back(AK8JetGenEta);
      AK8JetGenPhi_ .push_back(AK8JetGenPhi);
      int AK8JetGenJetIndex = -1;
      float AK8JetGenJetEn = -999.;
      float AK8JetGenJetPt = -999.;
      float AK8JetGenJetEta = -999.;
      float AK8JetGenJetPhi = -999.;
      if (doGenParticles_ && genParticlesHandle.isValid() ) {
	if ((*ijetAK8).genJet()) {
	  AK8JetGenJetIndex = 1;
	  AK8JetGenJetEn = (*ijetAK8).genJet()->energy();
	  AK8JetGenJetPt = (*ijetAK8).genJet()->pt();
	  AK8JetGenJetEta = (*ijetAK8).genJet()->eta();
	  AK8JetGenJetPhi = (*ijetAK8).genJet()->phi();
	}
      }
      AK8JetGenJetIndex_.push_back(AK8JetGenJetIndex);
      AK8JetGenJetEn_.push_back(AK8JetGenJetEn);
      AK8JetGenJetPt_.push_back(AK8JetGenJetPt);
      AK8JetGenJetEta_.push_back(AK8JetGenJetEta);
      AK8JetGenJetPhi_.push_back(AK8JetGenJetPhi);
      
      
      //save Softdrop subjet info Lvdp
      vecSoftdropSubjetcsv.clear();
      vecSoftdropSubjetpt.clear();
      vecSoftdropSubjeteta.clear();
      vecSoftdropSubjetmass.clear();
      vecSoftdropSubjetphi.clear();
      vecSoftdropSubjete.clear();
      vecSoftdropSubjetcharge.clear();
      vecSoftdropSubjetflavour.clear();
      nsubjets = 0;
      if(dumpSoftDrop_) {
	auto const & sdSubjets = ijetAK8->subjets("SoftDrop");
	for ( auto const & softdropsubjet : sdSubjets ) {
	  nsubjets++;
	  vecSoftdropSubjetpt.push_back(softdropsubjet->pt());
	  vecSoftdropSubjeteta.push_back(softdropsubjet->eta());
	  vecSoftdropSubjetmass.push_back(softdropsubjet->mass());
	  vecSoftdropSubjetphi.push_back(softdropsubjet->phi());
	  vecSoftdropSubjete.push_back(softdropsubjet->energy());
	  vecSoftdropSubjetflavour.push_back(abs(softdropsubjet->partonFlavour()));
	  vecSoftdropSubjetcharge.push_back(softdropsubjet->charge());
	  vecSoftdropSubjetcsv.push_back(softdropsubjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
	}
      }
      nAK8softdropSubjet_.push_back(nsubjets);
      AK8softdropSubjetPt_.push_back(vecSoftdropSubjetpt);
      AK8softdropSubjetEta_.push_back(vecSoftdropSubjeteta);
      AK8softdropSubjetPhi_.push_back(vecSoftdropSubjetphi);
      AK8softdropSubjetMass_.push_back(vecSoftdropSubjetmass);
      AK8softdropSubjetE_.push_back(vecSoftdropSubjete);
      AK8softdropSubjetCharge_.push_back(vecSoftdropSubjetcharge);
      AK8softdropSubjetFlavour_.push_back(vecSoftdropSubjetflavour);
      AK8softdropSubjetCSV_.push_back(vecSoftdropSubjetcsv);
    }
  }
  delete jecUnc;
  delete AK8jecUnc;
}
