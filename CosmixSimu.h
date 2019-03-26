#ifndef __CosmixSimu__
#define __CosmixSimu__

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <sstream>

#include "TMath.h"
#include "TFitResult.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

// to do
// - add axis titles
// - add legend
// - add momenta spectra
// - code cleaning
// - code configuration
// - deal with edges
// - tree: add timing info
// - treat case where particle enter in the second detector but not the first one

/*---------------------------------------------------------------------------------------------
|
| Measure detection efficiency of cosmix detector.
| - Z axis is defined by the pointing direction of the cosmix telescope
| - X axis is parallel to the detectors length
| - Y axis is constructed thanks to Z and Y with the usual convention
|
| Define the angle between z axis and vertical direction as "zenith".
| And suppose X axis is horizontal.
|
| Only geometry is taken into account. If both detectors are touched, muon is detected.
|
| I use usual spherical coordinates:
| -phi is the angle between muon direction and Z axis.
| -theta is the angle of the projection of the muon direction in the XY plane and the X axis.
|
| Distance [cm]
| Energy   [MeV]
---------------------------------------------------------------------------------------------*/


class CosmixSimu{

  private:
	//Specify detector Geometry
  	double_t ThickZ_; // Thickness of detectors in the pointing direction
 	double_t ThickY_; // Thickness of detectors in the orthogonal direction
  	double_t Length_; // Length of detectors (i.e: in x direction)
  	//The properties below depends on the material
	double_t Z_;
	double_t A_;
	double_t rho_; //density in g.cm^-3
	double_t dEdx_; //MeV/cm
	double_t ElossThickZ_; // Mean energy loss for MIP traversing the detector with a normal angle [MeV]
	TRandom *tr;
  
        //Histograms
        TH2F *hRat;
        TH2F *hAll;
        TH2F *hPass;
  	TH1F *h1ThetaReso;
  	TH1F *h1Rat;
  	TH1F *h1All;
  	TH1F *h1Pass;
  	TH1F *hCosmic;
  	TH1F *hCosmicAccepted;
    	TH1F *hPathLenght;
	TH1F *hEloss;
	TH1F *hRandEloss;
	//Canvas
	TCanvas *chRat;
  	TCanvas *cCosMuZenith;
    	TCanvas *chCosmic;
	TCanvas *cEloss;
	//TTree
	TTree* tree;

  public:

    CosmixSimu();
    ~CosmixSimu();
 
    double GetSolidAngle( double& err, double& thetaReso, double Distance = 1.2, double Zenith = 0., double NPseudoExp = 1e5, bool verbose = true, bool Draw = false, bool WaitPrimitive=true, bool Write = false, TFile* ofile = 0, bool GenerateTree = false, TString Type = "Cosmic", bool FirstDetectorIsDot = true); // Distance is the distance between detectors in cm
    void SetGeometry(double thickZ, double thickY, double lenght);
    //Z & A of the material - 
    //rho [g/cm^3]
    //dEdX MeV per cm
    void SetMaterialProperties(double Z, double A, double rho, double dEdX);
    //path lengh in cm
    //beta of the incoming particle
    double GetLandauWidth(double pathlenght = 2, double beta = 0.998);
    
  private:
    double GetWeight(double CosMuZenith, TString Type = "Cosmic");
    double GetCosinusMuonToZenith(double phi, double theta, double Zenith);
    double GetPathLenght(double phi, double theta, double Zenith);
    // Does not yet treat the edge effects
    double GetEloss(double phi, double theta, double Zenith);
    double GetRandEloss(double phi, double theta, double Zenith, double beta = 0.998);
};

#endif
