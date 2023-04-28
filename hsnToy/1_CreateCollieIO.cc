/* This script takes in a root file containing histograms for the signal, background and data distributions and creates a CollieIOfile that can be read by the successive scripts.

The input file may contain several histograms for signal, for each of the combination of parameters used in the simulation (e.g. mass and mixing angle).

For each combination of parameters a createMassPoint method is called, which saves a mass point in the output file. Each mass point takes in a signal, data and vector of background histograms.
*/

#include "CollieIOFile.hh"
#include "TRandom.h"

int main()
{
  // Open root input file
  TFile infile("/pc2014-data5/sporzio/collie/hsnTest/Storage/toy.root");
  TH1D* data = (TH1D*)infile.Get("Data");
  TH1D* bkgd1 = (TH1D*)infile.Get("Background1");
  TH1D* bkgd2 = (TH1D*)infile.Get("Background2");

  // Specify outputfile and channel name 
  CollieIOFile* cfile = new CollieIOFile();
  cfile->initFile("/pc2014-data5/sporzio/collie/hsnTest/Storage/toyHSN_CollieIOfile.root", "ToyHsnChannel");  

  // Define your input histograms
  double Xmin = 0.1; 
  double Xmax = 0.5;
  int Nbins = 40;
  cfile->setInputHist(Xmin,Xmax,Nbins);
  cfile->setRebin(1);

  // Define background names
  vector<string> bkgdNames;
  bkgdNames.push_back("Background1");
  bkgdNames.push_back("Background2");
  cfile->createChannel(bkgdNames);
  
  // Vector of backgrounds
  vector<TH1D*> vbkgd;
  vbkgd.push_back(bkgd1);
  vbkgd.push_back(bkgd2);
  // Vector of useless alpha parameters
  // (for smoothing, which we are not using)
  vector<double> valpha;
  valpha.push_back(-1);
  valpha.push_back(-1);

  // Create a parameter point to analyze
  int mass[5] = {150, 200, 250, 300, 350};
  for (int m : mass)
  {
    std::string histName = std::string("Signal_m") + std::to_string(m);
    TString tHistName(histName);
    TH1D* signal = (TH1D*)infile.Get(tHistName);
    cfile->createMassPoint(m, data, signal, -1, vbkgd, valpha);
  }
    
  cfile->storeFile();
}
