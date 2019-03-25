#ifndef __DATAANALYZIS_H_
#define __DATAANALYZIS_H_

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>


#include <TFile.h>
#include <TF1.h>
#include <TF1Convolution.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>

using namespace std;


//to do
// - read directly a file from the card
// - superimpose fit
// - change fit color
// - read config from a file
// - convolution gaus x poisson ... 


class DataAnalysis{

  private:
    ifstream ifile;
    TFile* ofile;
    vector<float> xval_;
    vector<float> yval_;
    vector<float> xvalerr_;
    vector<float> yvalerr_;

    map<string,TF1*> functions;

  public:
    DataAnalysis();
    DataAnalysis(string ifilename, string ofilename = "results.root");
    ~DataAnalysis();


    bool OpenFile(string ifilename);
    bool LoadData(string filename="", int scheme = 2);
    void SetAbsXErrors(float err);
    void SetRelXErrors(float relErr);
    bool SetRelXerrors(vector<float> errs);
    void SetAbsYErrors(float err);
    void SetRelYErrors(float relErr);
    bool SetRelYerrors(vector<float> errs);
    //void Fit(string model, int type); //type 1 (histo) 2 (graph)
    void Fit(string model, string xtitle, string ytitle, int type = 1, int nbins = 20, float xmin = 0, float xmax = 100);
    void Reset();
    map<string,TF1*>& GetModels() {return functions;};
    
  private:
    void LoadModels();
    void FitResultReport(const TFitResultPtr& res, int npar);
};
#endif
