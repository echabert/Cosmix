#include "CosmixSimu.h"

using namespace std;


void CosmixSimu::SetMaterialProperties(double Z, double A, double rho, double dEdX){
  Z_ = Z;
  A_ = A;
  rho_ = rho;
  dEdx_ = dEdX;
}

double CosmixSimu::GetLandauWidth(double pathlenght, double beta){
  /*
  //Cesium 55 137
  //Iode 17 131
  double Z = (55+17.)/2;
  double A = (137+131.)/2;
  double rho = 4.51; // g.cm-3
  */
  double K = 0.307;
  double x = rho_*pathlenght;
  //double beta = 0.998; //mean beta for muons
  //Forume is extracted from http://pdg.lbl.gov/2010/reviews/rpp2010-rev-passage-particles-matter.pdf,
  //Fluctuations in energy loss - 27.2.7
  double w = 4*K/2*Z_/A_*x/(beta*beta);
  return w;
}

CosmixSimu::CosmixSimu(){
  ThickZ_ =  2.;
  ThickY_ =  3.;
  Length_ = 15.;
  //Cesium 55 137
  //Iode 17 131
  Z_ = (55+17.)/2;
  A_ = (137+131.)/2;
  rho_ = 4.51; // g.cm-3
  dEdx_ = 7;
  ElossThickZ_ = dEdx_*ThickZ_;
  tr = new TRandom();
  hRat = 0;
  hAll = 0;
  hPass = 0;
  h1ThetaReso = 0;
  h1Rat = 0;
  h1All = 0;
  h1Pass = 0;
  hCosmic = 0;
  hCosmicAccepted = 0;
  hPathLenght = 0;
  hEloss = 0;
  hRandEloss = 0;
  chRat = 0;
  cCosMuZenith = 0;
  chCosmic = 0;
  cEloss = 0;
  tree = 0;
}

CosmixSimu::~CosmixSimu(){
 delete tr;
}

void CosmixSimu::SetGeometry(double thickZ, double thickY, double lenght){
  ThickZ_ = thickZ;
  ThickY_ = thickY;
  Length_ = lenght;
}


double CosmixSimu::GetWeight(double CosMuZenith, TString Type)
{
  // Cosmics mainly come from above, the distribution is not isotrope.
  // To get the correct cosmic distribution from an isotropic one, one has to weight events.
  // In file Cosmix/Publis/Grieder2001.pdf, you will find info on zenithal dependance of the flux
  // See for exemple eq 1.42, or zenith dependance of atmospheric depth (1.56, 1.57...)
  // Or better eq 3.30:
  //    (I(CosMuZenith) = I(1)*CosMuZenith^n
  //       where n = 1.85+-0.1
  //       Valid for Zenith < 75° ie CosMuZenith > 0.2588
  double CosmicWeight = -99.;

  if( Type == "Cosmic" )
  {
    //if( CosMuZenith < 0.2588 ) CosmicWeight = 0;
    if( CosMuZenith < 0. ) CosmicWeight = 0;
    else CosmicWeight = pow(CosMuZenith,1.85);
  } else if( Type == "Iso" )
  {
    CosmicWeight = 1.;
  } else if( Type == "Vert99" )
  {
    if( CosMuZenith > 0.99 ) CosmicWeight = 1.;
    else CosmicWeight = 0.;
  } else if( Type == "Vert999" )
  {
    if( CosMuZenith > 0.999 ) CosmicWeight = 1.;
    else CosmicWeight = 0.;
  } else
  {
    cout << "Try to weight events with a non existing Type. Weight will be -99." << endl;
  }
  
  return CosmicWeight;
}

double CosmixSimu::GetCosinusMuonToZenith(double phi, double theta, double Zenith)
{
  // In a right handed Cartesian coordinate system where
  // -X is horizonthal and
  // -Z axis has been turned by an angle "Zenith" arround the X axis, from the vertical direction,
  // if "theta" and "phi" are the spherical coordinate of a particle direction  
  // then, the cos of the angle of this particle with respect to the vertical direction is given as:
  // CosMuZenith = sin(theta)sin(phi)sin(Zenith)+cos(phi)cos(Zenith);
  // 
  double CosMuZenith = -99.;
  CosMuZenith = sin(theta)*sin(phi)*sin(Zenith)+cos(phi)*cos(Zenith);
  
  return CosMuZenith;
}

double CosmixSimu::GetPathLenght(double phi, double theta, double Zenith)
{
   return fabs(ThickZ_/GetCosinusMuonToZenith(phi,theta,Zenith));
}
 
double CosmixSimu::GetEloss(double phi, double theta, double Zenith){
   return GetPathLenght(phi, theta, Zenith)*ElossThickZ_/ThickZ_;
}
    
double CosmixSimu::GetRandEloss(double phi, double theta, double Zenith, double beta){
    double Eloss = GetEloss(phi,theta,Zenith);
    //cout<<phi<<" "<<theta<<" "<<Zenith<<" "<<GetPathLenght(phi,theta,Zenith)<<" "<<GetLandauWidth(GetPathLenght(phi,theta,Zenith),beta)<<endl;
    return tr->Landau(Eloss, GetLandauWidth(GetPathLenght(phi,theta,Zenith),beta));
}


double CosmixSimu::GetSolidAngle( double& err, double& thetaReso, double Distance, double Zenith, double NPseudoExp, bool verbose, bool Draw, bool WaitPrimitive, bool Write, TFile* ofile, bool GenerateTree, TString Type, bool FirstDetectorIsDot) // Distance is the distance between detectors in cm
{
  double pi = acos(-1);
  Zenith *= (pi/180.);	// Put angles in rad
  double CosMuZenith = -99.;

  int nbin_phi = 100;
  int nbin_theta = 100;
  
  // Produce 2D histos to have an idea of the efficiency vs angle
  // 2D histo, TH2F("h","h", ntheta,thetamin,thetamax, nphi,phimin,phimax)
  hRat = new TH2F("hRat","hRat", nbin_theta, 0, 2*pi, nbin_phi, 0., pi);
  hAll = new TH2F("hAll","hAll", nbin_theta, 0, 2*pi, nbin_phi, 0., pi);
  hPass = new TH2F("hPass","hPass", nbin_theta, 0, 2*pi, nbin_phi, 0., pi);

  // Produce acceptance vs Cos(Zenith angle)
  h1ThetaReso = new TH1F("h1ThetaReso", "", 100, -pi, pi);
  h1Rat = new TH1F("h1Rat","h1Rat", 100, -1., 1.);
  h1All = new TH1F("h1All","h1All", 100, -1., 1.);
  h1Pass = new TH1F("h1Pass","h1Pass", 100, -1., 1.);

  // Produce 1D histo of accepted muons vs zenith angle.
  hCosmic = new TH1F("hCosmic", "hCosmic", 100, -1., 1.);
  hCosmicAccepted = new TH1F("hCosmicAccepted", "hCosmicAccepted", 100, -1., 1.);
  
  //Produce 1D histo for energy loss
  hPathLenght = new TH1F("hPathLenght","Path lengh",100,0,5);
  hEloss = new TH1F("hEloss", "Energy loss", 100, 0, 50);
  hRandEloss = new TH1F("hRandEloss", "Energy loss", 100, 0, 50);

  


  // Also use coounters to compute solid angle
  double Nmuon = 0;
  double NDetectedmuon = 0;

  // I suppose a muon is starting (P1) at the edge of first detector so its coordinate is P1_z = ThickZ
  // I will then choose random angles and X and Y coordinate inside P1 and chek if muon pass in detector 2.
  double P1_x = -1000;	// X Coordinate opf muon in detector 1... Will be scanned
  double P1_y = -1000;	// Y Coordinate opf muon in detector 1... Will be scanned
  double P1_z = ThickZ_;	// Z Coordinate opf muon in detector 1. Fixed on the face of detector 1 nearest to det 2.

  double P2_x = -1000.;	// X Coordinate opf muon in detector 2... Will be computed
  double P2_y = -1000.;	// Y Coordinate opf muon in detector 2... Will be computed
  double P2_z = ThickZ_+Distance;	// Z Coordinate opf muon in detector 2. Fixed on the face of detector 2 nearest to det 1.

  double phi = -99, theta = -99.;
  double Eloss = -99, RandEloss = -99.; 
  double PathLenght = -99.; 
  bool coinc = false;
  if(GenerateTree){
    tree = new TTree("tree","");
    tree->Branch("phi",&phi);
    tree->Branch("theta",&theta);
    tree->Branch("pathLength",&PathLenght);
    tree->Branch("Eloss",&Eloss);
    tree->Branch("RandEloss",&RandEloss);
    tree->Branch("coinc",&coinc);
  }

  //int PseudoExp = 1e5;
  for( int i = 0; i < NPseudoExp; ++i ) { // Loop for P1_x
      //P1_x = tr->Rndm()*Length;
      if( FirstDetectorIsDot ){ 
        P1_x = Length_/2.; 
    	P1_y = ThickY_/2.;
      }
      else{
        P1_x = tr->Rndm()*Length_;
        P1_y = tr->Rndm()*ThickY_;
      }  
      phi = tr->Rndm()* pi;
      theta = tr->Rndm()* pi;

      CosMuZenith = GetCosinusMuonToZenith(phi, theta, Zenith);
      P2_x = P1_x+(P2_z-P1_z)*tan(phi)*cos(theta);
      P2_y = P1_y+(P2_z-P1_z)*tan(phi)*sin(theta);
      hAll->Fill(theta, phi);
      h1All->Fill(CosMuZenith,sin(phi));
      // Weight by sin(phi) to make isotropic flux, then weight according to muonflux
      Nmuon += sin(phi)*GetWeight(CosMuZenith,Type);
      // and by the cosine of the angle of the detector to the vertical i.e: cos(Zenith)
      // Nmuon += sin(phi)*GetWeight(CosMuZenith,Type)*cos(Zenith);
      hCosmic->Fill(CosMuZenith, sin(phi)*GetWeight(CosMuZenith,Type));
      
      coinc = false;
      if( P2_y >= 0. && P2_y <= ThickY_ && P2_x >=0. && P2_x <= Length_ )
      {
        coinc = true;
        hPass->Fill(theta, phi);
        h1Pass->Fill(CosMuZenith,sin(phi));
        // Weight by sin(phi) to make isotropic flux, then weight according to muonflux
        NDetectedmuon += sin(phi)*GetWeight(CosMuZenith,Type);
        // and by the cosine of the angle of the detector to the vertical i.e: cos(Zenith)
        //NDetectedmuon += sin(phi)*GetWeight(CosMuZenith,Type)*cos(Zenith);
        hCosmicAccepted->Fill(CosMuZenith, sin(phi)*GetWeight(CosMuZenith,Type));
        double vtheta = theta;
	if(theta>TMath::Pi()/2) vtheta=theta-TMath::Pi();
	h1ThetaReso->Fill(vtheta,sin(phi));
        PathLenght = GetPathLenght(phi,theta,Zenith);
	hPathLenght->Fill(PathLenght);
	Eloss = GetEloss(phi,theta,Zenith);
	hEloss->Fill(Eloss);
	RandEloss = GetRandEloss(phi,theta,Zenith);
	hRandEloss->Fill(GetRandEloss(phi,theta,Zenith));
      }
      if(GenerateTree) tree->Fill();
  }
  
  //Compute the efficiency and write a report (if verbose)
  double eff = NDetectedmuon/Nmuon;
  err = sqrt(eff*(1-eff)/Nmuon);
  if(verbose){
    cout << "For detector distance = " << Distance << ", Zenith pointing angle " << Zenith << " the cosmix efficiency is : " << eff << " +/- "<<err<<" |  NDetectedmuon = " << NDetectedmuon << " Nmuon = " << Nmuon << endl;
  }
  hRat->Divide(hPass, hAll, 1, 1, "B");
  h1Rat->Divide(h1Pass, h1All, 1, 1, "B");

  //Estimate the resolution on theta
  TF1 f("f","gaus+[3]",pi/2,pi/2);
  f.SetParameter(0,1000);
  f.SetParameter(1,0);
  f.SetParameter(2,0.5);
  f.SetParameter(3,100);
  TFitResultPtr res = h1ThetaReso->Fit(&f,"S");
  thetaReso = res.Get()->GetParams()[2];
  cout<<"Resolution on theta: "<<h1ThetaReso->GetStdDev()<<endl;
  
  if( Draw )
  {
    chRat = new TCanvas("chRat", "plots",200,10,700,750);
    /*
    hRat->Draw("lego2 sph");
    chRat->Update();
    */
    h1ThetaReso->Draw();
    if( WaitPrimitive )chRat->WaitPrimitive();
  
    cCosMuZenith = new TCanvas("cCosMuZenith", "plots",200,10,700,750);
    cCosMuZenith->Divide(1,2);
    cCosMuZenith->cd(1);
    h1All->SetMinimum(0.);
    h1All->Draw();
    h1Pass->SetMarkerColor(2);
    h1Pass->Draw("same");
    cCosMuZenith->cd(2);
    h1Rat->Draw();
    cCosMuZenith->Update();
    if( WaitPrimitive )cCosMuZenith->WaitPrimitive();
  
    chCosmic = new TCanvas("chCosmic", "plots",200,10,700,750);
    hCosmic->SetMinimum(0.);
    hCosmic->Draw();
    hCosmicAccepted->SetMarkerColor(2);
    hCosmicAccepted->Draw("same");
    chCosmic->Update();
    if( WaitPrimitive )chCosmic->WaitPrimitive();
    
    cEloss = new TCanvas("cEloss","",200,10,700,750);
    hEloss->Draw();
    hRandEloss->SetLineColor(kRed);
    hRandEloss->Draw("same");
    if( WaitPrimitive) cEloss->WaitPrimitive();

  }
  

  if( Draw && Write && ofile ){
      ofile->cd();
      stringstream ss;
      ss<<"d_"<<fixed<<setprecision(2)<<Distance<<"_"<<fixed<<setprecision(2)<<Zenith;
      ofile->mkdir(ss.str().c_str());
      ofile->cd(ss.str().c_str());

      if(GenerateTree) tree->Write();

      hRat->Write();
      hAll->Write();
      hPass->Write();
      h1ThetaReso->Write();
      h1Rat->Write();
      h1All->Write();
      h1Pass->Write();
      hCosmic->Write();
      hCosmicAccepted->Write();
      hPathLenght->Write();
      hEloss->Write();
      hRandEloss->Write();
      chRat->Write();
      cCosMuZenith->Write();
      chCosmic->Write();
      cEloss->Write();
      ofile->cd();
      ofile->Write();
  }
  //delete pointers
  /*
  delete hRat;
  delete hAll;
  delete hPass;
  delete h1ThetaReso;
  delete h1Rat;
  delete h1All;
  delete h1Pass;
  delete hCosmic;
  delete hCosmicAccepted;
  delete hPathLenght;
  delete hEloss;
  delete hRandEloss;
  delete chRat;
  delete cCosMuZenith;
  delete chCosmic;
  delete cEloss;
 */
  return eff;
}

