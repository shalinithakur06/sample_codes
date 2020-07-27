#include<iostream>
#include<fstream>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
void pile_up(){
  TFile *f1 = TFile::Open("down.root");
  TH1D *h_data =  (TH1D*) f1->Get("pileup");

  TFile *f2 = TFile::Open("Out_Moriond2016.root");
  TH1D *h_mc =  (TH1D*) f2->Get("pileup");

  TH1F *pu_wt = new TH1F("pu_wt","pu_wt",80,0.,80.); 

  h_data->Scale(1./h_data->Integral());
  h_mc->Scale(1./h_mc->Integral());

 for(unsigned int i=1;i<=80;i++){
    double num = h_data->GetBinContent(i)/h_mc->GetBinContent(i);
    pu_wt->SetBinContent(i,num);
   } 
 TFile *fout = new TFile("Down.root","RECREATE");
 fout->cd();
 pu_wt->Write();
 fout->Close();   
}
