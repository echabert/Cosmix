#include "DataAnalysis.h"


//superpose ? sauvegarde ? historique ? report ? play with parameters ?
void DataAnalysis::Fit(string model, string xtitle, string ytitle, int type, int nbins, float xmin, float xmax){
   TCanvas c(model.c_str());
   if(type==1){
     TH1F* h = new TH1F("","",nbins,xmin,xmax);
     h->GetXaxis()->SetTitle(xtitle.c_str());
     h->GetYaxis()->SetTitle(ytitle.c_str());
     TFitResultPtr res = h->Fit(functions[model],"S");
     FitResultReport(res, functions[model]->GetNpar());
   }
   if(type==2){
     TGraphErrors* g = new TGraphErrors(xval_.size()-1,xval_.data(),yval_.data(),xvalerr_.data(), yvalerr_.data());
     g->GetXaxis()->SetTitle(xtitle.c_str());
     g->GetYaxis()->SetTitle(ytitle.c_str());
     TFitResultPtr res = g->Fit(functions[model],"S");
     FitResultReport(res, functions[model]->GetNpar());
     g->Draw("AP");
   }
   ofile->cd();
   c.Write();

}

void DataAnalysis::FitResultReport(const TFitResultPtr& res, int npar){
   ///*
   cout<<"######################################"<<endl;
   //cout<<res.Get()->GetFormula()->GetExpFormula()<<endl;
   //Afficher les parametres
   //int npar = res.Get()->GetNpar();
   for(int i=0;i<npar;i++){
    cout<<"["<<i<<"] = "<<res.Get()->GetParams()[i]<<" +/- " <<res.Get()->Error(i)<<endl;
   }
   cout<<"chi2 = "<<res.Get()->Chi2()<<endl;
   //Nombre de degree de liberte
   //int ndof = npoints -2 -1;
   //cout<<"prob = "<<TMath::Prob(p.Get()->Chi2(),ndof)<<endl;
   cout<<"prob = "<<res.Get()->Prob()<<endl;
   cout<<"######################################"<<endl;
   //*/
}

void DataAnalysis::LoadModels(){
    //functions["Poisson"] = new TF1("Poisson","TMath::Poisson",0,100);
    functions["pol1"] = new TF1("pol1","[0]*x+[1]",0,100);
    functions["Poisson"] = new  TF1("pois","[0]*TMath::Poisson(x,[1])",0,100);
    functions["Gauss"]= new TF1("Gauss","gaus",0,100);
    functions["cos2"] = new TF1("cos2","[0]*cos(x)*cos(x)+[1]",0,TMath::Pi());
    functions["cosn"] = new TF1("cosn","[0]*TMath::Power(cos(x),[2])+[1]",0,TMath::Pi());
    functions["exp"] = new TF1("expo","[0]*TMath::Exp(-x/[1])+[2]",0,100);
   
    TF1Convolution fCos2Gaus("cos2","gaus",0,TMath::Pi(),true);
    fCos2Gaus.SetRange(0.,TMath::Pi());
    fCos2Gaus.SetNofPointsFFT(1000);
    functions["cos2Gaus"] = new TF1("cos2Gaus",fCos2Gaus, 0.,TMath::Pi(), fCos2Gaus.GetNpar());
    functions["cos2Gaus"]->SetParameters(600.,500,0,0.5);

}

DataAnalysis::DataAnalysis(){
  ofile = new TFile("AnaResults.root","RECREATE");
  LoadModels();
}

DataAnalysis::DataAnalysis(string ifilename, string ofilename){
    ifile.open(ifilename.c_str());
    ofile = new TFile(ofilename.c_str(),"RECREATE");
    LoadModels();
}

DataAnalysis::~DataAnalysis(){
    ifile.close();
    ofile->Close();
    delete ofile;
}

bool DataAnalysis::OpenFile(string ifilename){
  if(ifile.is_open()) ifile.close();
  ifile.open(ifilename.c_str());
  return ifile.is_open();
}

void DataAnalysis::Reset(){
  xval_.clear();
  yval_.clear();
  xvalerr_.clear();
  yvalerr_.clear();
}


bool DataAnalysis::LoadData(string ifilename, int scheme){
    if(ifilename!=""){
      ifile.open(ifilename.c_str());
    }
    float x,y,xerr,yerr;
    while(!ifile.eof()){
       switch(scheme){
           case 2:
             ifile>>x>>y;
             xval_.push_back(x);
             yval_.push_back(y);
             xvalerr_.push_back(0.);
             yvalerr_.push_back(sqrt(y));
             break;
            case 3:
             ifile>>x>>y>>xerr;
             xval_.push_back(x);
             yval_.push_back(y);
             xvalerr_.push_back(xerr);
             yvalerr_.push_back(sqrt(y));
             break;
            case 4:
             ifile>>x>>y>>xerr>>yerr;
             xval_.push_back(x);
             yval_.push_back(y);
             xvalerr_.push_back(xerr);
             yvalerr_.push_back(yerr);
             break;
            default: 
              cerr<<"DataAnalysis::LoadData:: Selected scheme is not allowed"<<endl; 
              return false;
	}
    }
    return true;

}
    
void DataAnalysis::SetAbsXErrors(float err){
   xvalerr_.clear();
   xvalerr_ = vector<float>(xval_.size(),err);
}
    
void DataAnalysis::SetRelXErrors(float relErr){
  xvalerr_.clear();
  xvalerr_ = vector<float>(xval_.size());
  for(unsigned int i=0;i<xvalerr_.size();i++){
    xvalerr_[i] = xval_[i]*relErr;
  }
}

bool DataAnalysis::SetRelXerrors(vector<float> errs){
  xvalerr_ = errs;
  if(errs.size()!=xvalerr_.size()) return false;
  return true;
}

void DataAnalysis::SetAbsYErrors(float err){
   yvalerr_.clear();
   yvalerr_ = vector<float>(yval_.size(),err);
}
    
void DataAnalysis::SetRelYErrors(float relErr){
  yvalerr_.clear();
  yvalerr_ = vector<float>(yval_.size());
  for(unsigned int i=0;i<yvalerr_.size();i++){
    yvalerr_[i] = yval_[i]*relErr;
  }
}

bool DataAnalysis::SetRelYerrors(vector<float> errs){
  yvalerr_ = errs;
  if(errs.size()!=yvalerr_.size()) return false;
  return true;
}
