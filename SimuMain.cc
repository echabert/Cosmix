#include <iostream>

#include <TApplication.h>
#include <TMath.h>

using namespace std;

#include "CosmixSimu.h"

//compilation
//g++ `root-config --glibs --cflags `  CosmixSimu.h Solidmain.C 

int main(int argc, char* argv[]){
   //TApplication app ("app",&argc,argv); 

   //Histograms
   TH1F hEffDistance("hEffDistance","",12,0,12);
   TH1F hThetaResoDistance("hThetaResoDistance","",12,0,12);
   TH1F hEffZenith("hEffZenith","",10,0,TMath::Pi()/2);
   TFile ofile("SimResults.root","RECREATE");

   CosmixSimu omega;
   double eff;
   double err;
   double thetaReso;

   //variation en fonction de l'angle
   for(float i = 0; i<10; i++){
    eff = omega.GetSolidAngle(err, thetaReso, 5,TMath::Pi()/2/10*i);
    hEffZenith.SetBinContent(i,eff);
   }


    //double NPseudoExp = 1e5, bool verbose     = true, bool Draw = false, bool WaitPrimitive=true, bool Write = false, TFile* ofile = 0
   eff = omega.GetSolidAngle(err, thetaReso, 5, 0, 1e5, true, true, false, true, &ofile, true);
   cout<<"We are here!"<<endl;
   for(float d = 0; d<12; d+=1){
    eff = omega.GetSolidAngle(err, thetaReso, d,1.5);
    hEffDistance.SetBinContent(d+1,eff);
    hEffDistance.SetBinError(d+1,err);
    hThetaResoDistance.SetBinContent(d+1,thetaReso);
   }
   cout<<"Width of landau: "<<omega.GetLandauWidth()<<endl;
   ofile.cd();
   hEffDistance.Write();
   hEffZenith.Write();
   hThetaResoDistance.Write();
   ofile.Write();
   ofile.Close();

}
