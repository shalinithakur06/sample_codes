#include <iostream>
#include <iomanip>
#include "TRandom3.h"
#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TTimeStamp.h"

//https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
const double JEREtaMap[14] = {0., 0.5, 0.8, 1.1, 1.3, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.0, 3.2, 5.0}; 
const double JERSF[13] = {1.109, 1.138, 1.114, 1.123, 1.084, 1.082, 1.140, 1.067, 1.177, 1.364, 1.857, 1.328, 1.16};
const double JERSFUp[13] = {1.109+0.008 , 1.138+0.013, 1.114+0.013, 1.123+0.024, 1.084+0.011, 1.082+0.035, 1.140+0.047, 1.067+0.053, 1.177+0.041, 1.364+0.039, 1.857+0.071, 1.328+0.022, 1.16+0.029};
const double JERSFDown[13] = {1.109-0.008 , 1.138-0.013, 1.114-0.013, 1.123-0.024, 1.084-0.011, 1.082-0.035, 1.140-0.047, 1.067-0.053, 1.177-0.041, 1.364-0.039, 1.857-0.071, 1.328-0.022, 1.16-0.029};

using namespace std;

double UncertaintyComputer::getJERSF(double eta, int jer){
  double SF = 1.0;
  for(size_t i = 0; i < 13; i++){
    if(TMath::Abs(eta) >= JEREtaMap[i] && TMath::Abs(eta) < JEREtaMap[i+1]){
      if(jer == 0)SF = JERSF[i];
      else if (jer == 1) SF = JERSFUp[i];
      else if(jer == -1) SF = JERSFDown[i];  
    }
  }
  return SF;
}

double UncertaintyComputer::jetPtWithJESJER(MyJet jet, int jes, int jer){

  double gen_pt = jet.Genp4.pt();  
  double jet_pt = jet.p4.pt();
  double sigmaJER = jet.resolution ;
  //apply JES uncert scaling 
  jet_pt *= (1+(jet.JECUncertainty*double(jes)));
  //apply JER uncert, scaling
  double delR = DeltaR(jet.Genp4, jet.p4);
  double rCone = 0.4;
  if(gen_pt> 0 && delR<rCone/2 && abs(jet_pt -gen_pt)<3*sigmaJER*jet_pt ){
  //if(gen_pt > 0){
    double SF = getJERSF(jet.p4.eta(), jer);
    double ptscale = max(0.0, 1.0 + (SF - 1)*(jet_pt - gen_pt)/ jet_pt);
    jet_pt *= ptscale;
  }
  return jet_pt;
}


