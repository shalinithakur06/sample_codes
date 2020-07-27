#define dyjets_cxx
#include "dyjets.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;
void dyjets::Loop()
{

   if (fChain == 0) return;
TFile *f = new TFile("dyjets.root","RECREATE");

TH1D *LeadMuonPt = new TH1D("LeadMuonPt","LeadMuon Pt",500,0,500);
TH1D *LeadMuonEta = new TH1D("LeadMuonEta","LeadMuon Eta",100, -5,5);
TH1D *LeadMuonPhi = new TH1D("LeadMuonPhi","LeadMuon Phi",64, -3.2,3.2);  

TH1D *kLeadMuonPt  = new TH1D("kLeadMuonPt","kLeadMuon Pt",500,0,500);
TH1D *kLeadMuonEta = new TH1D("kLeadMuonEta","kLeadMuon Eta",100,-5,5);
TH1D *kLeadMuonPhi = new TH1D("kLeadMuonPhi", "kLeadMuon Phi",64,-3.2,3.2);

 Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      cout << " No. of MC particles: " << nMC << "\t" << mcStatusFlag->size() << endl;
      vector<int> IndexJets;
      for(unsigned int  ivec=0; ivec !=mcStatusFlag->size(); ++ivec){
         if ( ((*mcStatusFlag)[ivec] ==7)|| ((*mcStatusFlag)[ivec] ==3) || ((*mcStatusFlag)[ivec] ==4) || ((*mcStatusFlag)[ivec] ==20) ){
           if( abs((*mcPID)[ivec]) ==13 && abs((*mcMomPID)[ivec]) !=4000013){
        	LeadMuonPt->Fill((*mcPt)[ivec]);
        	LeadMuonEta->Fill((*mcEta)[ivec]); 
        	LeadMuonPhi->Fill((*mcPhi)[ivec]);
        	if( abs((*mcEta)[ivec]) > 2.4) continue;
        	if(((*mcPt)[ivec]) < 24) continue;
	        kLeadMuonPt->Fill((*mcPt)[ivec]);
                kLeadMuonEta->Fill((*mcEta)[ivec]);
                kLeadMuonPhi->Fill((*mcPhi)[ivec]);
      	  }
       } 
     }
  }
f->cd();
LeadMuonPt->Write();
LeadMuonEta->Write();
LeadMuonPhi->Write();
kLeadMuonPt->Write();
kLeadMuonEta->Write();
kLeadMuonPhi->Write();
f->Close();
}
